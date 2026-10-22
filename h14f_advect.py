"""h14f_advect: does the eddy slowdown slope keep falling from D_mol/D_eff = 1/256 to 1/4096 (log law), or hold near 0.2 (power law)?
(open item 14 / H14, journals/ija_draft.md:252-261; follows h14e_advect.py; 2026-10-21.)

Question. h14e measured S50 at l/w 10 on six ratios, 1/8 ... 1/256: 5.76, 7.71, 10.27, 13.05, 15.28, 18.07 (h14e_advect.log); q(128 -> 256) =
0.24 [0.19..0.29], UNDECIDED. A log law and a power law held near 0.24 fit the ladder but give S50(1e9) = 72 vs 700. h14f runs ONE new cell four
doublings out, V4096_10 (48 runs = 16 jobs x 3), reloads 1/8 ... 1/128 from h14d_advect_raw.pkl and 1/256 from h14e_advect_raw.pkl, and scores
the pre-registered rule: FALLING if the 95% upper bound of q(256 -> 4096) < 0.19, HELD if its lower bound > 0.19, else UNDECIDED; plus BELOW EDGE
SLOPE and LOG LAW INSIDE. Bands P1-P8 and the rules: h14f_predictions.txt.

Nothing existing is edited: h14e_advect.py (spec, d14, reload_check, PKL), h14d, h14c, h14b and advect2d are imported. The one change inside the
model run is the kick: nkick below replaces h14b_advect.kick in each worker (the numpy kick needs ~2 GB of index temporaries per worker at
N = 4021). self_check asserts nkick is bitwise equal to h14b_advect.kick and that a tiny job gives identical output with either kick.

PREDICTIONS: h14f_predictions.txt, sha256 PRED_SHA (asserted on import). Imported files are hash-asserted too (FILE_SHA).
Run: cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv caffeinate -i nice -n 10 uv run --no-project \
       --with numpy --with scipy --with numba --with matplotlib python h14f_advect.py > h14f_advect.log 2> h14f_advect.err
     ("report" re-reads h14f_advect_raw.pkl; "smoke" = self-check and reloads only, no new runs; 4 worker processes, guard 43200 s)
"""
import hashlib
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
PRED_SHA = "d136fe694b8fe890833864590cb910d372424b770275f7519e169280f50c5c99"
FILE_SHA = {"h14e_advect.py": "22e7b592ddfdc859ba87d248725e470a9a96a8646f8548116b7ef2ab421188a2",
            "h14e_advect_raw.pkl": "3338b7347336e71d13e6fee2148f4d3a3582bafbef1fdc23ae39219a86b108f4"}
assert sha("h14f_predictions.txt") == PRED_SHA, "h14f_predictions.txt changed after it was hashed"
for _f, _h in FILE_SHA.items():
    assert sha(_f) == _h, f"{_f} changed since h14f was written"

import advect2d as AD
import h14b_advect as H
import h14c_advect as C
import h14d_advect as D
import h14e_advect as E              # asserts the h14e predictions and the h14d / advect2d / h14d pkl hashes; D asserts h14b / h14c

PKL = os.path.join(HERE, "h14f_advect_raw.pkl")
GUARD_S = 43200
NWORK = 4
NAME = "V4096_10"
RATIO = dict(E.RATIO, **{NAME: 4096})
LADDER = E.LADDER + (NAME,)
Q_THR = 0.19                                                   # 2.5% bound of q(128 -> 256), h14e_advect.log:73
A256 = 18.07                                                   # the h14e 1/256 anchor, taken as a fixed number
LOG6 = -1.98 + 3.59 * np.log(4096)                             # six-point log law (h14e_advect.log:79) at 1/4096
F_SPHERE, F_EDGE = E.F_SPHERE, E.F_EDGE
SCEN = {"A log law over six (-1.98 + 3.59 ln 4096)": LOG6, "B edge path q 0.135 (f 2400)": A256 * 16 ** 0.135,
        "B edge path q 0.141 (f 2600)": A256 * 16 ** 0.141, "C power law q 0.19 from 1/256": A256 * 16 ** 0.19,
        "D power law q 0.24 from 1/256": A256 * 16 ** 0.24, "E power law q 0.29 from 1/256": A256 * 16 ** 0.29,
        "F power law q 0.33 from 1/256": A256 * 16 ** 0.33}


def cases():
    return {f"{NAME}#{j}": E.spec(f"{NAME}#{j}", 4096, 1, 90, 3501 + j) for j in range(16)}


ORDER = list(cases())


# ----------------------------------------------------------------------------------------------- the numba kick
@nb.njit(cache=False)
def _gather(a, m, w, ax2):
    R, N = a.shape[0], a.shape[1]
    out = np.empty_like(a)
    for r in range(R):
        for i in range(N):
            for j in range(N):
                if ax2:
                    k, o = m[r, i], i
                    out[r, i, j] = (w[0, r, o] * a[r, i, (j - k - 2) % N] + w[1, r, o] * a[r, i, (j - k - 1) % N]
                                    + w[2, r, o] * a[r, i, (j - k) % N] + w[3, r, o] * a[r, i, (j - k + 1) % N])
                else:
                    k, o = m[r, j], j
                    out[r, i, j] = (w[0, r, o] * a[r, (i - k - 2) % N, j] + w[1, r, o] * a[r, (i - k - 1) % N, j]
                                    + w[2, r, o] * a[r, (i - k) % N, j] + w[3, r, o] * a[r, (i - k + 1) % N, j])
    return out


def nkick(a, s, h, axis):
    """h14b_advect.kick with the same float32 weights and the same summation order, as one loop (no index temporaries)."""
    sc = s / h
    m = np.floor(sc)
    u = (1 - (sc - m)).astype(np.float32)
    w = np.stack((-u * (u - 1) * (u - 2) / 6, (u + 1) * (u - 1) * (u - 2) / 2, -(u + 1) * u * (u - 2) / 2, (u + 1) * u * (u - 1) / 6))
    return _gather(a, m.astype(np.intp), w, axis == 2)


def job(args):
    H.kick = nkick                                             # run_case calls the module global
    return C.job(args)


# ----------------------------------------------------------------------------------------------- self-check and diagnostics
def kick_diff(N, R, reps, A, seed):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((R, N, N)).astype(np.float32)
    worst = 0.0
    for k in range(reps):
        sh, ax = H.shear(rng, dict(A=A, n=1), R, N, 1.0, 2 if k % 2 == 0 else 1, float(N))
        x, y = H.kick(a, sh, 1.0, ax), nkick(a, sh, 1.0, ax)
        assert x.dtype == y.dtype and x.shape == y.shape and y.flags.c_contiguous
        worst = max(worst, float(np.abs(x - y).max()))
        a = x
    return worst


def self_check():
    t0 = time.time()
    C.self_check()
    cs = cases()
    assert repr(E.spec("V256_10#0", 256, 7, 60, 3401)) == repr(E.cases()["V256_10#0"]), "recipe differs from h14e V256_10#0"
    seeds = [s["seed"] for s in cs.values()]
    taken = set(range(3101, 3133)) | {s["seed"] for s in D.cases().values()} | {s["seed"] for s in E.cases().values()}
    assert len(set(seeds)) == len(seeds) and not set(seeds) & taken, "seeds"
    sp = cs[ORDER[0]]
    fl = sp["flow"]
    assert sp["L"] == 4021 and abs(fl["A"] ** 2 / (8 * fl["Tk"]) - 4095) < 1e-6 and sp["gamma"] == 0.0074 and sp["cad"] == 2.2
    assert sp["D"] * sp["dt"] / sp["h"] ** 2 <= 0.25 + 1e-9 and abs(1 / sp["gamma"] / fl["Tk"] - 5) < 0.7
    for k in ("hold17", "S", "deltas", "h", "dt", "noise", "n", "lw", "D"):
        assert sp[k] == dict(E.cases()["V256_10#0"], S=1, hold17=90)[k], k

    # the kick: bitwise, on random fields, both axes, small and large shifts, and once at the real size and flow
    w = [kick_diff(64, 4, 12, 30.0, 64), kick_diff(257, 3, 12, 234.7, 257), kick_diff(4021, 1, 2, fl["A"], 4021)]
    assert max(w) == 0.0, w
    # end to end: one tiny job (D_eff 16, 3 runs, hold 2 t17) with either kick gives the same output
    tiny = E.spec("T16#0", 16, 1, 2, 999)
    orig, out = H.kick, []
    for kk in (orig, nkick):
        H.kick = kk
        out.append(C.job(("T16#0", tiny, time.time() + 900))[1])
    H.kick = orig
    r0, r1 = out
    same = lambda r: pickle.dumps({k: v for k, v in r.items() if k != "elapsed"})
    assert same(r0) == same(r1), "tiny job differs between the two kicks"
    AD.L = sp["L"]
    dd = float(np.mean([AD.tracer_lagrangian(fl, s) for s in (1, 7)]))
    kap = fl["A"] / fl["l"]
    diag = dict(tracer=(dd, 4095.0), damp=C.kick_damping(fl["Tk"]), lyap=D.lyapunov(kap),
                flow=dict(kappa=kap, gamma=sp["gamma"], kicks_pre=1 / sp["gamma"] / fl["Tk"], L=sp["L"], Tk=fl["Tk"], A=fl["A"], l=fl["l"]))
    print(f"h14f self-check passed ({time.time() - t0:.0f}s): h14c self_check; recipe at D_eff 256 reproduces h14e V256_10#0 exactly; seeds unique and new; "
          f"nkick bitwise equal to h14b kick (max |diff| {max(w)} at N 64, 257, 4021); tiny job identical with either kick "
          f"({int(np.asarray(r0['walled_at_t0']).sum())}/{r0['runs']} walled, t_single {np.round(r0['t_single'], 1).tolist()}); "
          f"{len(cs)} jobs / {sum(s['S'] * len(s['deltas']) for s in cs.values())} runs", flush=True)
    return diag


def ladder6():
    P, st = D.stats(E.d14()["res"])
    r = pickle.load(open(E.PKL, "rb"))["res"]
    P[E.NAME] = C.pool([r[k] for k in sorted(r, key=lambda s: int(s.split("#")[1]))])
    st[E.NAME] = C.Stat(P[E.NAME], seed=sum(map(ord, E.NAME)))
    return P, st


def reload_check():
    E.reload_check()
    _, st = ladder6()
    assert f"{st[E.NAME].pt[0]:.2f}" == "18.07", st[E.NAME].pt[0]
    print("h14e reproduction: S50(1/256) 18.07 reloaded from h14e_advect_raw.pkl with the h14e bootstrap seed", flush=True)


# ----------------------------------------------------------------------------------------------- report
def report(res, elapsed, diag):
    P, st = ladder6()
    parts = [res[k] for k in sorted(res, key=lambda s: int(s.split("#")[1]))]
    P[NAME] = C.pool(parts)
    st[NAME] = C.Stat(P[NAME], seed=sum(map(ord, NAME)))
    fmt, line = H.fmt, H.line
    sv = lambda k, q=0: st[k].pt[q]
    S = lambda k: st[k].bs[:, 0]
    print("hashes: " + "; ".join(f"{f} {h}" for f, h in {"h14f_predictions.txt": PRED_SHA, **FILE_SHA}.items()))

    f, (dd, nom), lam = diag["flow"], diag["tracer"], diag["lyap"]
    print("\nREALISED FLOW (new case; AD.L set to the box, mean of 2 seeds) and regime diagnostic")
    print(f"  D_eff 4096 L {f['L']:4.0f} Tk {f['Tk']:5.1f} A {f['A']:6.1f} eddy l {f['l']:6.1f}: kappa {f['kappa']:.3f} kicks_pre {f['kicks_pre']:.2f} gamma {f['gamma']:g}; "
          f"D_eddy measured {dd:.1f} nominal {nom:.0f} (ratio {dd / nom:.3f}); D_num/D_mol {diag['damp']:.4f}; "
          f"Lyapunov exponent {lam:.3f} per kick ({lam / f['Tk']:.4f} per time unit)")

    print("\nCASES  (L, ratio, l/w, Tk, A, gamma, t17, runs, walled at freeze-out, hold, cens at end, cut, cpu s)")
    for k in (E.NAME, NAME):
        s, fl = P[k]["spec"], P[k]["spec"]["flow"]
        print(f"  {k:8s} L={s['L']:g} Deff={s['Deff']} 1/{s['Deff']} l/w={fl['l'] / np.sqrt(s['Deff']):.2f} Tk={fl['Tk']:g} A={fl['A']:.1f} g={s['gamma']:g} t17={C.t17(s['L'], s['Deff']):.1f} "
              f"runs={P[k]['runs']} walled={st[k].nw} ({st[k].pt[3]:.3f}) hold={s['hold']:g} cens_at_end={int(np.isinf(P[k]['t_single'][P[k]['walled_at_t0']]).sum())} "
              f"cut={P[k]['cut']} cpu={P[k]['elapsed']:.0f}")

    print("\nSLOWDOWN at l/w 10  (time since the end of the sweep / t17; Kaplan-Meier over runs walled at freeze-out; bootstrap over runs, 2000x; 1/8 ... 1/256 reloaded)")
    for k in LADDER:
        s = st[k]
        print(f"  {k:8s} {line(s, 0, 'S50')} | {line(s, 1, 'S80')} | {line(s, 2, 'S_half')}")
        print(f"           walled {s.pt[3]:.3f}  mean domains at t0 {P[k]['ndom0'].mean():.2f}  kicks to clear (S50 t17 / Tk) "
              f"{s.pt[0] * s.t17 / P[k]['spec']['flow']['Tk']:.0f}  cleared by 16 / 24 / 32 / 48 t17: " + " ".join(f"{C.km_cleared(P[k], s, m):.2f}" for m in (16, 24, 32, 48)))

    def ci(b):
        fin = np.isfinite(b)
        return np.percentile(b[fin], [2.5, 97.5]) if (1 - fin.mean()) <= 0.05 and fin.sum() > 100 else np.array([np.nan, np.nan])

    def q2(ka, kb):
        l = np.log(RATIO[kb] / RATIO[ka])
        return np.log(sv(kb) / sv(ka)) / l, np.log(S(kb) / S(ka)) / l

    def sf(pt, b):
        c = ci(b)
        return f"{fmt(pt)} [95% {fmt(c[0])}..{fmt(c[1])}]"

    print("\nLOCAL ratio exponents q(a -> b) = ln(S50_b / S50_a) / ln(r_b / r_a)  (h14e: 0.42, 0.41, 0.34, 0.23, 0.24 [0.19..0.29])")
    for ka, kb in zip(LADDER[:-1], LADDER[1:]):
        print(f"  1/{RATIO[ka]:<4d} -> 1/{RATIO[kb]:<4d}: q = {sf(*q2(ka, kb))}")
    q_pt, q_b = q2(E.NAME, NAME)
    qp_pt, qp_b = q2("V128_10", E.NAME)
    dq_pt, dq_b = q_pt - qp_pt, q_b - qp_b
    c_ = ci(q_b)
    print(f"  q(256 -> 4096) = {sf(q_pt, q_b)}; q(128 -> 256) = {sf(qp_pt, qp_b)}; Delta q = {sf(dq_pt, dq_b)} (paired over the same replicates)")

    def cluster(pts, B=2000, seed=0):
        s = C.Stat(C.pool(pts), B=0)
        jid = np.concatenate([np.full(p["runs"], j) for j, p in enumerate(pts)])
        rng, nj, out = np.random.default_rng(seed), len(pts), np.empty(B)
        for b in range(B):
            idx = np.concatenate([np.flatnonzero(jid == j) for j in rng.integers(0, nj, nj)])
            out[b] = s._est(idx)[0]
        return out
    r256 = pickle.load(open(E.PKL, "rb"))["res"]
    cl = np.log(cluster(parts, seed=4096) / cluster([r256[k] for k in sorted(r256, key=lambda s: int(s.split("#")[1]))], seed=256)) / np.log(16)
    print(f"  job-level bootstrap (printed, NOT the scored interval): q(256 -> 4096) 95% {fmt(ci(cl)[0])}..{fmt(ci(cl)[1])}")

    # descriptive refits over seven ratios (not scored)
    pt = np.array([sv(k) for k in LADDER])
    bs = np.stack([S(k) for k in LADDER], 1)
    x7 = np.log(np.array([RATIO[k] for k in LADDER], float))
    s4, c4 = sv(NAME), ci(S(NAME))
    print(f"\nPER-DOUBLING INCREMENTS of S50: " + ", ".join(f"{x:.2f}" for x in np.diff(pt[:-1])) + f"; 256 -> 4096 per doubling {(s4 - pt[-2]) / 4:.2f} "
          f"[95% {ci((bs[:, 6] - bs[:, 5]) / 4)[0]:.2f}..{ci((bs[:, 6] - bs[:, 5]) / 4)[1]:.2f}] (a log law keeps it near 2.5)")
    print(f"SIX-POINT LOG LAW at 1/4096: {LOG6:.2f}; measured - log law = {s4 - LOG6:+.2f}")
    fin = np.isfinite(bs).all(1)
    bl, al = np.polyfit(x7, pt, 1)
    fits = np.array([np.polyfit(x7, v, 1) for v in bs[fin]])
    res_l = pt - (al + bl * x7)
    bq, aq = np.polyfit(x7, np.log(pt), 1)
    res_q = pt - np.exp(aq + bq * x7)
    print(f"LOG LAW refit over seven ratios: S50 = {al:.2f} + {bl:.2f} ln r  [b 95% {ci(fits[:, 0])[0]:.2f}..{ci(fits[:, 0])[1]:.2f}]; residuals " + ", ".join(f"{x:+.2f}" for x in res_l)
          + f"; rms {np.sqrt((res_l ** 2).mean()):.2f}   (h14e six-point -1.98 + 3.59 ln r, rms 0.20)")
    print(f"POWER LAW refit over seven ratios: q = {bq:.2f}; residuals " + ", ".join(f"{x:+.2f}" for x in res_q) + f"; rms {np.sqrt((res_q ** 2).mean()):.2f}   (h14e six-point: 0.33, rms 0.66)")

    print(f"\nSCENARIOS for S50(1/4096) written in h14f_predictions.txt, against the measured {s4:.2f} [95% {c4[0]:.2f}..{c4[1]:.2f}]")
    for lab, v in SCEN.items():
        print(f"  {lab:42s} {v:6.2f}  implied q {np.log(v / A256) / np.log(16):.3f}  inside the measured 95% interval: {bool(c4[0] <= v <= c4[1])}  (measured - scenario = {s4 - v:+.2f})")

    print("\nEXTRAPOLATION (NOT a measurement; 9 doublings measured (1/8 to 1/4096), 18 more extrapolated to 1e9): S50 at the real ratio 1e9, anchored at the 1/4096 point,"
          f" and the Early Earth threshold S* = f/{F_SPHERE} for f = {F_EDGE[0]}, {F_EDGE[1]}")
    sstar = {fe: fe / F_SPHERE for fe in F_EDGE}
    qstar = {fe: np.log(v / s4) / np.log(1e9 / 4096) for fe, v in sstar.items()}
    print("  S* = " + ", ".join(f"{v:.0f} (f={fe})" for fe, v in sstar.items()) + "; q* (the power-law q that reaches S* at 1e9 from the 1/4096 point): "
          + ", ".join(f"{qstar[fe]:.3f} (f={fe})" for fe in F_EDGE))
    c9 = ci(fits[:, 1] + fits[:, 0] * np.log(1e9))
    rows = [("power law, q(256->4096) 2.5% bound", c_[0]), ("power law, q(256->4096) point", q_pt), ("power law, q(256->4096) 97.5% bound", c_[1]),
            ("log law refit over seven", None), ("constant from 1/4096", 0.0)]
    for lab, q in rows:
        if q is not None and not np.isfinite(q):
            continue
        if q is None:
            v9 = al + bl * np.log(1e9)
            cross = {fe: np.exp((sstar[fe] - al) / bl) for fe in F_EDGE}
        else:
            v9 = s4 * (1e9 / 4096) ** max(q, 0)
            cross = {fe: (4096 * (sstar[fe] / s4) ** (1 / q) if q > 0 else np.inf) for fe in F_EDGE}
        extra = f" [95% {c9[0]:.0f}..{c9[1]:.0f}]" if q is None else ""
        qs = "" if q is None else f" (q {q:5.3f})"
        print(f"  {lab:36s}{qs:>13s}: S50(1e9) = {v9:8.4g}{extra}, f = {F_SPHERE} S = {F_SPHERE * v9:8.4g}; ratio at which S50 reaches S*: "
              + ", ".join(f"f={fe}: {cross[fe]:.2e}" for fe in F_EDGE) + f"; below S* at 1e9: {bool(F_SPHERE * v9 < F_EDGE[1])}")

    # ---------------------------------------------------------------------------------------- predictions
    def band(k, qi, lo, hi):
        v = sv(k, qi)
        return bool(lo <= v <= hi) if np.isfinite(v) else (False if st[k].cens_floor() >= hi else "UNDECIDED")

    def AND(*xs):
        xs = [x if isinstance(x, str) else bool(x) for x in xs]
        if any(x is False for x in xs):
            return False
        return "UNDECIDED" if any(isinstance(x, str) for x in xs) else True
    inband = lambda v, lo, hi: bool(lo <= v <= hi) if np.isfinite(v) else "UNDECIDED"
    out = []

    def verdict(n, ok, txt):
        ok = ok if isinstance(ok, str) else bool(ok)
        out.append((n, ok))
        print(f"P{n} {ok if isinstance(ok, str) else ('HELD' if ok else 'FAILED')}  {txt}")
    print("\nPREDICTIONS (bands from h14f_predictions.txt)")
    verdict(1, band(NAME, 0, 24, 42), f"V4096_10 S50 {fmt(s4)} in [24, 42] (guess 29)")
    verdict(2, band(NAME, 1, 28, 50), f"V4096_10 S80 {fmt(sv(NAME, 1))} in [28, 50] (guess 34)")
    verdict(3, inband(q_pt, 0.10, 0.30), f"q(256 -> 4096) = {fmt(q_pt)} in [0.10, 0.30] (guess 0.17)")
    verdict(4, inband(dq_pt, -0.16, 0.06), f"Delta q = {fmt(dq_pt)} in [-0.16, +0.06] (guess -0.07); q(128 -> 256) = {fmt(qp_pt)}")
    verdict(5, AND(*[sv(kb) > sv(ka) for ka, kb in zip(LADDER[:-1], LADDER[1:])]), "S50 ladder " + " < ".join(fmt(sv(k)) for k in LADDER))
    ok6 = abs(dd / nom - 1) <= 0.05 and diag["damp"] <= 0.02 and 1.45 <= f["kappa"] <= 1.49 and 4.5 <= f["kicks_pre"] <= 5.5 and lam > 0.1
    verdict(6, ok6, f"flow 1/4096: D_eddy ratio {dd / nom:.3f}, D_num/D_mol {diag['damp']:.4f}, kappa {f['kappa']:.3f}, kicks_pre {f['kicks_pre']:.2f}, Lyapunov {lam:.3f} per kick")
    cens = int(np.isinf(P[NAME]["t_single"][P[NAME]["walled_at_t0"]]).sum())
    verdict(7, st[NAME].pt[3] >= 0.99 and cens == 0, f"walled fraction {st[NAME].pt[3]:.3f} (>= 0.99), walled runs still walled at the end of the hold {cens} (0)")
    cut = any(v["cut"] for v in res.values())
    nruns = sum(v["runs"] for v in res.values())
    verdict(8, elapsed < GUARD_S and not cut and nruns == 48, f"total {elapsed:.0f} s ({NWORK} workers, guard {GUARD_S} s), cut by guard {cut}, new runs {nruns} (48)")
    print("\nSUMMARY  held: " + ", ".join(f"P{n}" for n, ok in out if ok is True) + " | failed: " + ", ".join(f"P{n}" for n, ok in out if ok is False)
          + " | undecided: " + ", ".join(f"P{n}" for n, ok in out if isinstance(ok, str)))
    if st[NAME].pt[3] < 1 or cens > 0 or cut:
        print(f"FOUND: walled fraction {st[NAME].pt[3]:.3f}, {cens} walled runs censored at the end of the hold, cut by guard {cut}; no run was dropped")

    ok = np.isfinite(c_).all()
    dec = "UNDECIDED" if not ok else ("FALLING" if c_[1] < Q_THR else ("HELD" if c_[0] > Q_THR else "UNDECIDED"))
    edge = "BELOW EDGE SLOPE" if ok and c_[1] < qstar[F_EDGE[1]] else "NOT BELOW"
    inside = "LOG LAW INSIDE" if np.isfinite(c4).all() and c4[0] <= LOG6 <= c4[1] else ("OUTSIDE (measured above)" if s4 > LOG6 else "OUTSIDE (measured below)")
    print(f"DECISION (pre-registered): {dec}; q(256 -> 4096) = {fmt(q_pt)} [95% {fmt(c_[0])}..{fmt(c_[1])}] (rule: FALLING if upper < {Q_THR}, HELD if lower > {Q_THR})")
    print(f"SECOND LABEL (pre-registered): {edge}; upper bound {fmt(c_[1])} vs q* {qstar[F_EDGE[1]]:.3f} (f {F_EDGE[1]}; {qstar[F_EDGE[0]]:.3f} at f {F_EDGE[0]})")
    print(f"THIRD LABEL (pre-registered): {inside}; six-point log law {LOG6:.2f}, S50(1/4096) {fmt(s4)} [95% {fmt(c4[0])}..{fmt(c4[1])}]")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    t0 = time.time()
    if mode == "report":
        d = pickle.load(open(PKL, "rb"))
        report(d["res"], d["elapsed"], d["diag"])
        return
    print("hashes asserted: " + "; ".join(f"{f} {h}" for f, h in {"h14f_predictions.txt": PRED_SHA, **FILE_SHA}.items()), flush=True)
    diag = self_check()
    reload_check()
    cs = cases()
    print(f"\njobs (load average at start {', '.join(f'{x:.2f}' for x in os.getloadavg())}):", flush=True)
    for k in ORDER:
        s = cs[k]
        print(f"  {k}: L={s['L']:g} gamma={s['gamma']:g} D_eff={s['Deff']} runs={s['S'] * len(s['deltas'])} hold={s['hold']:g} ({s['hold17']} t17) seed={s['seed']} "
              f"Tk={s['flow']['Tk']:g} A={s['flow']['A']:.1f} cad={s['cad']:g}", flush=True)
    if mode == "smoke":
        return
    deadline = t0 + GUARD_S
    res = {}
    with ProcessPoolExecutor(NWORK) as ex:
        futs = [ex.submit(job, (k, cs[k], deadline)) for k in ORDER]
        for fu in as_completed(futs):
            if fu.cancelled():
                continue
            k, r = fu.result()
            res[k] = r
            print(f"  job done: {k} [{r['elapsed']:.0f}s] walled {int(r['walled_at_t0'].sum())}/{r['runs']} cleared {int(np.isfinite(r['t_single'][r['walled_at_t0']]).sum())} "
                  f"cut={r['cut']} (elapsed {time.time() - t0:.0f}s)", flush=True)
            pickle.dump(dict(res=res, elapsed=time.time() - t0, diag=diag, partial=True), open(PKL + ".part", "wb"))
            if time.time() > deadline:
                n = sum(f.cancel() for f in futs)
                print(f"  guard {GUARD_S} s reached: {n} jobs not yet started were cancelled", flush=True)
    elapsed = time.time() - t0
    pickle.dump(dict(res=res, elapsed=elapsed, diag=diag), open(PKL, "wb"))
    report(res, elapsed, diag)


if __name__ == "__main__":
    main()
