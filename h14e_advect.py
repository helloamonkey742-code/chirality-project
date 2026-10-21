"""h14e_advect: does the local ratio exponent q of the eddy slowdown S50 at l/w = 10 keep falling at D_mol/D_eff = 1/256?
(open item 14 / H14, journals/ija_draft.md:252; follows h14d_advect.py; 2026-10-21.)

Question. h14d measured S50 at 1/8, 1/16, 1/32, 1/64, 1/128 = 5.76, 7.71, 10.27, 13.05, 15.28 (h14d_advect.log:96-108); the local exponent fell
0.42, 0.41, 0.34, 0.23 [0.13..0.28] and the pre-registered label for q(64 -> 128) was UNDECIDED. A log law (post hoc: -1.77 + 3.52 ln r) and a
power law that keeps falling give the same S50 at 1/256 (17.75, 17.92); a power law held at q 0.36 / 0.40 gives 19.61 / 20.16. h14e runs ONE new
cell, V256_10 (336 runs = 16 jobs x 21, the V128_10 design), reloads 1/8 ... 1/128 from h14d_advect_raw.pkl (and, through it, h14c_advect_raw.pkl),
and scores the pre-registered rule: SATURATING if the 95% upper bound of q(128 -> 256) < 0.2, GROWING if its lower bound > 0.3, else UNDECIDED, plus
CONSTANT-SLOPE EXCLUDED if the 95% upper bound of S50(1/256) < 15.28 x 2^0.36. Bands P1-P9 and the rules: h14e_predictions.txt.

Nothing existing is edited: h14d_advect.py (cases recipe, stats, cell_sets, lyapunov, repro_check), h14c_advect.py, h14b_advect.py and advect2d.py
are imported. The V256_10 recipe below is h14d_advect.cases().add() with D_eff = 256; self_check asserts that the same recipe at D_eff = 128 and
seed 3301 gives h14d's V128_10#0 spec exactly. Jobs of the new cell are pooled in job order (0 ... 15), not completion order, so that the bootstrap
indices do not depend on which worker finished first (the point estimates do not depend on the order at all).

PREDICTIONS: h14e_predictions.txt, sha256 PRED_SHA (asserted on import). Imported files are hash-asserted too (FILE_SHA).
Run: cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project \
       --with numpy --with scipy --with numba --with matplotlib python h14e_advect.py > h14e_advect.log 2> h14e_advect.err
     ("report" re-reads h14e_advect_raw.pkl; "smoke" = self-check and reloads only, no new runs; 3 worker processes, guard 7200 s)
"""
import hashlib
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
PRED_SHA = "6ea98ce64f21fa6ccc67c9eeb45a3763aefe168a9116a2c8ca7c91a1724cd865"
FILE_SHA = {"h14d_advect.py": "10476e92c21aa052c8134aeaef46dfba4d4c9065cb429ba7e3d7d32b01900721",
            "advect2d.py": "7d84a2e6e8c05c4d807f945a3f1cd3db40115ab87f32ae053e499558e145cf2e",
            "h14d_advect_raw.pkl": "c3c8f8b9b8cafe23027f5ae7232ad59e2c8e7b3aa8c357454bb0043e5ff70288"}
assert sha("h14e_predictions.txt") == PRED_SHA, "h14e_predictions.txt changed after it was hashed"
for _f, _h in FILE_SHA.items():
    assert sha(_f) == _h, f"{_f} changed since h14e was written"

import advect2d as AD
import h14b_advect as H
import h14c_advect as C              # asserts the h14c predictions sha256 on import
import h14d_advect as D              # asserts the h14d predictions sha256 and the h14b / h14c / advect2d file hashes on import

PKL = os.path.join(HERE, "h14e_advect_raw.pkl")
D14 = os.path.join(HERE, "h14d_advect_raw.pkl")
GUARD_S = 7200
NWORK = 3
NAME = "V256_10"
RATIO = dict(D.RATIO, **{NAME: 256})
LADDER = ("R8_10", "R16_10", "V32_10", "V64_10", "V128_10", NAME)
Q_RULE = D.Q_RULE                                              # SATURATING if q_hi < 0.2; GROWING if q_lo > 0.3
A128 = 15.28                                                   # the h14d 1/128 anchor (h14d_advect.log:108), taken as a fixed number
F_SPHERE = 17                                                  # f = 17 S (journals/new_analyses_h13h14_combined.md:26)
F_EDGE = (2400, 2600)                                          # the printed 99.9 share holds below about f = 2600 (h13h14_extreme_locate.log:249-251)
SCEN = {"A log law (-1.77 + 3.52 ln 256)": -1.77 + 3.52 * np.log(256), "B falling power (q 0.23 from 1/128)": A128 * 2 ** 0.23,
        "C power law q 0.36 from 1/128": A128 * 2 ** 0.36, "D power law q 0.40 from 1/128": A128 * 2 ** 0.40}


def spec(name, Deff, S, hold17, seed):
    """The recipe of h14d_advect.cases().add() (flow and gamma exactly as h14c), unchanged."""
    L, fl = C.flow_for(Deff, 10)
    gamma = round(min(0.01, 1 / (5 * fl["Tk"])), 4)           # 5 kicks between lam = 0 and the end of the sweep
    T = C.t17(L, Deff)
    return dict(L=float(L), n=1, h=1.0, dt=0.2, S=S, D=C.D_MOL, flow=fl, noise="x", seed=seed, gamma=gamma, Deff=Deff, lw=10,
                deltas=(0.5, 1.0, 2.0), hold17=hold17, hold=float(np.ceil(hold17 * T / 2) * 2), cad=max(2.0, round(T / 60 / 0.2) * 0.2), name=name)


def cases():
    return {f"{NAME}#{j}": spec(f"{NAME}#{j}", 256, 7, 60, 3401 + j) for j in range(16)}


ORDER = list(cases())


# ----------------------------------------------------------------------------------------------- self-check and diagnostics
def d14():
    return pickle.load(open(D14, "rb"))


def self_check():
    t0 = time.time()
    C.self_check()                                             # h14c's own: nk fix, g from the actual gamma, tiny job()/pool(), 15 worlds
    cs, dc = cases(), D.cases()
    assert repr(spec("V128_10#0", 128, 7, 60, 3301)) == repr(dc["V128_10#0"]), "recipe differs from h14d V128_10#0"
    assert repr(spec("V64_10n#3", 64, 7, 48, 3324)) == repr(dc["V64_10n#3"]), "recipe differs from h14d V64_10n#3"
    seeds = [s["seed"] for s in cs.values()]
    taken = set(range(3101, 3133)) | {s["seed"] for s in dc.values()}
    assert len(set(seeds)) == len(seeds) and not set(seeds) & taken, "seeds"
    sp = next(iter(cs.values()))
    fl = sp["flow"]
    assert abs(fl["A"] ** 2 / (8 * fl["Tk"]) - 255) < 1e-9
    N = sp["L"] / sp["h"]
    assert abs(N - round(N)) < 1e-9 and sp["D"] * sp["dt"] / sp["h"] ** 2 <= 0.25 + 1e-9
    assert abs(1 / sp["gamma"] / fl["Tk"] - 5) < 0.7
    for k in ("flow", "gamma", "hold17", "S", "deltas", "h", "dt", "noise", "n", "lw", "D"):          # same design as V128_10 apart from size
        assert (cs[f"{NAME}#0"][k] == dc["V128_10#0"][k]) or k in ("flow", "gamma"), k
    AD.L = sp["L"]
    dd = float(np.mean([AD.tracer_lagrangian(fl, s) for s in (1, 7)]))
    kap = fl["A"] / fl["l"]
    diag = dict(tracer=(dd, 255.0), damp=C.kick_damping(fl["Tk"]), lyap=D.lyapunov(kap),
                flow=dict(kappa=kap, gamma=sp["gamma"], kicks_pre=1 / sp["gamma"] / fl["Tk"], L=sp["L"], Tk=fl["Tk"], A=fl["A"], l=fl["l"]))
    print(f"h14e self-check passed ({time.time() - t0:.0f}s): h14c self_check; recipe at D_eff 128 / 64 reproduces the h14d V128_10#0 / V64_10n#3 specs exactly; "
          f"seeds unique and outside 3101-3132 and 3301-3344; D_eddy = D_eff - 1 by construction; Euler stability; kicks_pre within 5 +- 0.7; "
          f"{len(cs)} jobs / {sum(s['S'] * len(s['deltas']) for s in cs.values())} runs", flush=True)
    return diag


def reload_check():
    D.repro_check()
    P, ST = D.stats(d14()["res"])
    want = dict(R8_10="5.76", R16_10="7.71", V32_10="10.27", V64_10="13.05", V128_10="15.28")
    for k, w in want.items():
        assert f"{ST[k].pt[0]:.2f}" == w, (k, ST[k].pt[0], w)
    print("h14d reproduction: S50 ladder 5.76, 7.71, 10.27, 13.05, 15.28 reloaded from h14d_advect_raw.pkl with the h14d bootstrap seeds", flush=True)


# ----------------------------------------------------------------------------------------------- report
def report(res, elapsed, diag):
    P, st = D.stats(d14()["res"])
    parts = [res[k] for k in sorted(res, key=lambda s: int(s.split("#")[1]))]
    P[NAME] = C.pool(parts)
    st[NAME] = C.Stat(P[NAME], seed=sum(map(ord, NAME)))
    fmt, line = H.fmt, H.line
    sv = lambda k, q=0: st[k].pt[q]
    S = lambda k: st[k].bs[:, 0]
    print("hashes: " + "; ".join(f"{f} {h}" for f, h in {"h14e_predictions.txt": PRED_SHA, **FILE_SHA}.items()))

    f, (dd, nom), lam = diag["flow"], diag["tracer"], diag["lyap"]
    print("\nREALISED FLOW (new case; AD.L set to the box, mean of 2 seeds) and regime diagnostic")
    print(f"  D_eff 256 L {f['L']:4.0f} Tk {f['Tk']:5.1f} A {f['A']:6.1f} eddy l {f['l']:6.1f}: kappa {f['kappa']:.3f} kicks_pre {f['kicks_pre']:.2f} gamma {f['gamma']:g}; "
          f"D_eddy measured {dd:.2f} nominal {nom:.0f} (ratio {dd / nom:.3f}); D_num/D_mol {diag['damp']:.4f}; "
          f"Lyapunov exponent {lam:.3f} per kick ({lam / f['Tk']:.4f} per time unit)")

    print("\nCASES  (L, ratio, l/w, Tk, A, gamma, t17, runs, walled at freeze-out, hold, cens at end, cut, cpu s)")
    for k in ("V64_10", "V128_10", NAME):
        s, fl = P[k]["spec"], P[k]["spec"]["flow"]
        print(f"  {k:8s} L={s['L']:g} Deff={s['Deff']} 1/{s['Deff']} l/w={fl['l'] / np.sqrt(s['Deff']):.2f} Tk={fl['Tk']:g} A={fl['A']:.1f} g={s['gamma']:g} t17={C.t17(s['L'], s['Deff']):.1f} "
              f"runs={P[k]['runs']} walled={st[k].nw} ({st[k].pt[3]:.3f}) hold={s['hold']:g} cens_at_end={int(np.isinf(P[k]['t_single'][P[k]['walled_at_t0']]).sum())} "
              f"cut={P[k]['cut']} cpu={P[k]['elapsed']:.0f}")

    print("\nSLOWDOWN at l/w 10  (time since the end of the sweep / t17; Kaplan-Meier over runs walled at freeze-out; bootstrap over runs, 2000x; 1/8 ... 1/128 reloaded)")
    for k in LADDER:
        s = st[k]
        print(f"  {k:8s} {line(s, 0, 'S50')} | {line(s, 1, 'S80')} | {line(s, 2, 'S_half')}")
        print(f"           walled {s.pt[3]:.3f}  mean domains at t0 {P[k]['ndom0'].mean():.2f}  kicks to clear (S50 t17 / Tk) "
              f"{s.pt[0] * s.t17 / P[k]['spec']['flow']['Tk']:.0f}  cleared by 12 / 16 / 24 / 32 t17: " + " ".join(f"{C.km_cleared(P[k], s, m):.2f}" for m in (12, 16, 24, 32)))

    def ci(b):
        fin = np.isfinite(b)
        return np.percentile(b[fin], [2.5, 97.5]) if (1 - fin.mean()) <= 0.05 and fin.sum() > 100 else np.array([np.nan, np.nan])

    def q2(ka, kb):
        l = np.log(RATIO[kb] / RATIO[ka])
        return np.log(sv(kb) / sv(ka)) / l, np.log(S(kb) / S(ka)) / l

    def ls(keys):
        x = np.array([RATIO[k] for k in keys], float)
        return H.fit_p(x, np.array([sv(k) for k in keys])), np.array([H.fit_p(x, v) for v in np.stack([S(k) for k in keys], 1)])

    def sf(pt, b):
        c = ci(b)
        return f"{fmt(pt)} [95% {fmt(c[0])}..{fmt(c[1])}]"

    print("\nLOCAL ratio exponents q(a -> b) = ln(S50_b / S50_a) / ln(r_b / r_a)  (h14d: 8 -> 16 0.42; 16 -> 32 0.41; 32 -> 64 0.34; 64 -> 128 0.23 [0.13..0.28])")
    for ka, kb in zip(LADDER[:-1], LADDER[1:]):
        print(f"  1/{RATIO[ka]:<3d} -> 1/{RATIO[kb]:<3d}: q = {sf(*q2(ka, kb))}")
    q_ref, b_ref = ls(("R8_10", "V32_10", "V64_10"))
    q64, b64 = ls(("V64_10", "V128_10", NAME))
    q_all, b_all = ls(LADDER)
    q_hi_pt, b_hi = q2("V128_10", NAME)
    q_prev_pt, b_prev = q2("V64_10", "V128_10")
    dq_pt, dq_b = q_hi_pt - q_prev_pt, b_hi - b_prev
    print(f"  q(8 -> 64), least squares over 8, 32, 64 (h14c definition): {sf(q_ref, b_ref)}   [h14d: 0.40 [0.34..0.48]]")
    print(f"  q(64 -> 256), least squares over 64, 128, 256: {sf(q64, b64)}")
    print(f"  q over all six ratios (8 ... 256), least squares: {sf(q_all, b_all)}   [h14d, five ratios: 0.36 [0.32..0.41]]")
    print(f"  q(128 -> 256) = {sf(q_hi_pt, b_hi)}; q(64 -> 128) = {sf(q_prev_pt, b_prev)}; Delta q = q(128 -> 256) - q(64 -> 128) = {sf(dq_pt, dq_b)} (paired over the same replicates)")
    c_ = ci(b_hi)
    print(f"  q(128 -> 256) 95% interval excludes 0.5 (cellular): {bool(c_[1] < 0.5) if np.isfinite(c_[1]) else 'n/a'}; excludes 0 (chaotic limit): {bool(c_[0] > 0) if np.isfinite(c_[0]) else 'n/a'}; "
          f"upper bound < 0.36 (paired check of the second label): {bool(c_[1] < 0.36) if np.isfinite(c_[1]) else 'n/a'}")
    x6 = np.log(np.array([RATIO[k] for k in LADDER], float))
    quad = lambda v: np.polyfit(x6, np.log(v), 2)[0] if np.isfinite(v).all() else np.nan
    cq_b = np.array([quad(v) for v in np.stack([S(k) for k in LADDER], 1)])
    print(f"  curvature c of ln S50 = a + q ln r + c (ln r)^2 over the six ratios: {sf(quad(np.array([sv(k) for k in LADDER])), cq_b)}   [h14d, five ratios: -0.05 [-0.11..-0.01]]")

    # job-level (cluster) bootstrap of q(128 -> 256): runs in one job share a seed stream
    def cluster(pts, B=2000, seed=0):
        s = C.Stat(C.pool(pts), B=0)
        jid = np.concatenate([np.full(p["runs"], j) for j, p in enumerate(pts)])
        rng, nj, out = np.random.default_rng(seed), len(pts), np.empty(B)
        for b in range(B):
            idx = np.concatenate([np.flatnonzero(jid == j) for j in rng.integers(0, nj, nj)])
            out[b] = s._est(idx)[0]
        return out
    cl = np.log(cluster(parts, seed=256) / cluster(D.cell_sets(d14()["res"])["V128_10"][0], seed=128)) / np.log(2)
    print(f"  job-level bootstrap (printed, NOT the scored interval): q(128 -> 256) 95% {fmt(ci(cl)[0])}..{fmt(ci(cl)[1])}")

    # per-doubling increments and the log law refitted over six points (POST-registration descriptive, not scored)
    pt = np.array([sv(k) for k in LADDER])
    bs = np.stack([S(k) for k in LADDER], 1)
    inc = np.diff(pt)
    ib = bs[:, 5] - bs[:, 4]
    print("\nPER-DOUBLING INCREMENTS of S50: " + ", ".join(f"{x:.2f}" for x in inc) + f"; last (128 -> 256) {inc[-1]:.2f} [95% {ci(ib)[0]:.2f}..{ci(ib)[1]:.2f}]; "
          f"mean of the previous four {inc[:4].mean():.2f} (a log law keeps the increment, a power law with q 0.36 would make it {A128 * (2 ** 0.36 - 1):.2f})")
    bl, al = np.polyfit(x6, pt, 1)
    fits = np.array([np.polyfit(x6, v, 1) for v in bs])
    res_l = pt - (al + bl * x6)
    bq, aq = np.polyfit(x6, np.log(pt), 1)
    res_q = pt - np.exp(aq + bq * x6)
    print(f"LOG LAW refit over six ratios: S50 = {al:.2f} + {bl:.2f} ln r  [b 95% {ci(fits[:, 0])[0]:.2f}..{ci(fits[:, 0])[1]:.2f}]; residuals " + ", ".join(f"{x:+.2f}" for x in res_l)
          + f"; rms {np.sqrt((res_l ** 2).mean()):.2f}   (h14d five-point fit -1.77 + 3.52 ln r, rms 0.19)")
    print(f"POWER LAW refit over six ratios: q = {bq:.2f}; residuals " + ", ".join(f"{x:+.2f}" for x in res_q) + f"; rms {np.sqrt((res_q ** 2).mean()):.2f}   (h14d five-point: 0.36, rms 0.48)")

    s256, c256 = sv(NAME), ci(S(NAME))
    print(f"\nSCENARIOS for S50(1/256) written in h14e_predictions.txt, against the measured {s256:.2f} [95% {c256[0]:.2f}..{c256[1]:.2f}]")
    for lab, v in SCEN.items():
        print(f"  {lab:38s} {v:6.2f}  implied q {np.log(v / A128) / np.log(2):.2f}  inside the measured 95% interval: {bool(c256[0] <= v <= c256[1])}  (measured - scenario = {s256 - v:+.2f})")

    print("\nEXTRAPOLATION (NOT a measurement; 5 doublings measured (1/8 to 1/256), 22 more extrapolated to 1e9): S50 at the real ratio 1e9, anchored at the 1/256 point,"
          f" and the Early Earth threshold S* = f/{F_SPHERE} for f = {F_EDGE[0]}, {F_EDGE[1]}")
    sstar = {fe: fe / F_SPHERE for fe in F_EDGE}
    print("  S* = " + ", ".join(f"{v:.0f} (f={fe})" for fe, v in sstar.items()) + "; the power-law q that reaches S* at 1e9 from the 1/256 point: "
          + ", ".join(f"{np.log(v / s256) / np.log(1e9 / 256):.3f} (f={fe})" for fe, v in sstar.items()))
    rows = [("power law, q(128->256) 2.5% bound", ("p", c_[0])), ("power law, q(128->256) point", ("p", q_hi_pt)), ("power law, q(128->256) 97.5% bound", ("p", c_[1])),
            ("power law, q(64->256) point", ("p", q64)), ("log law refit over six", ("l", None)), ("constant from 1/256", ("p", 0.0))]
    for lab, (kind, q) in rows:
        if kind == "p" and not np.isfinite(q):
            continue
        if kind == "p":
            v9 = s256 * (1e9 / 256) ** max(q, 0)
            cross = {fe: (256 * (sstar[fe] / s256) ** (1 / q) if q > 0 else np.inf) for fe in F_EDGE}
        else:
            v9 = al + bl * np.log(1e9)
            cross = {fe: np.exp((sstar[fe] - al) / bl) for fe in F_EDGE}
        c9 = ci(fits[:, 1] + fits[:, 0] * np.log(1e9))
        extra = f" [95% {c9[0]:.0f}..{c9[1]:.0f}]" if kind == "l" else ""
        qs = "" if kind == "l" else f" (q {q:5.2f})"
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
    print("\nPREDICTIONS (bands from h14e_predictions.txt)")
    verdict(1, band(NAME, 0, 16.2, 20.8), f"V256_10 S50 {fmt(s256)} in [16.2, 20.8] (guess 17.9)")
    verdict(2, band(NAME, 1, 19.5, 25.5), f"V256_10 S80 {fmt(sv(NAME, 1))} in [19.5, 25.5] (guess 21.6)")
    verdict(3, inband(q_hi_pt, 0.05, 0.40), f"q(128 -> 256) = {fmt(q_hi_pt)} in [0.05, 0.40] (guess 0.21)")
    verdict(4, inband(dq_pt, -0.25, 0.15), f"Delta q = {fmt(dq_pt)} in [-0.25, +0.15] (guess -0.02); q(64 -> 128) = {fmt(q_prev_pt)}")
    verdict(5, AND(*[sv(kb) > sv(ka) for ka, kb in zip(LADDER[:-1], LADDER[1:])]), "S50 ladder " + " < ".join(fmt(sv(k)) for k in LADDER))
    verdict(6, inband(q64, 0.12, 0.35), f"q(64 -> 256) = {fmt(q64)} in [0.12, 0.35] (guess 0.23)")
    ok7 = abs(dd / nom - 1) <= 0.05 and diag["damp"] <= 0.02 and 1.45 <= f["kappa"] <= 1.49 and 4.5 <= f["kicks_pre"] <= 5.5 and lam > 0.1
    verdict(7, ok7, f"flow 1/256: D_eddy ratio {dd / nom:.3f}, D_num/D_mol {diag['damp']:.4f}, kappa {f['kappa']:.3f}, kicks_pre {f['kicks_pre']:.2f}, Lyapunov {lam:.3f} per kick")
    cens = int(np.isinf(P[NAME]["t_single"][P[NAME]["walled_at_t0"]]).sum())
    verdict(8, st[NAME].pt[3] >= 0.99 and cens == 0, f"walled fraction {st[NAME].pt[3]:.3f} (>= 0.99), walled runs still walled at the end of the hold {cens} (0)")
    cut = any(v["cut"] for v in res.values())
    nruns = sum(v["runs"] for v in res.values())
    verdict(9, elapsed < GUARD_S and not cut and nruns == 336, f"total {elapsed:.0f} s ({NWORK} workers, guard {GUARD_S} s), cut by guard {cut}, new runs {nruns} (336)")
    print("\nSUMMARY  held: " + ", ".join(f"P{n}" for n, ok in out if ok is True) + " | failed: " + ", ".join(f"P{n}" for n, ok in out if ok is False)
          + " | undecided: " + ", ".join(f"P{n}" for n, ok in out if isinstance(ok, str)))
    if st[NAME].pt[3] < 1 or cens > 0 or cut:
        print(f"FOUND: walled fraction {st[NAME].pt[3]:.3f}, {cens} walled runs censored at the end of the hold, cut by guard {cut}; no run was dropped")

    def rule(lo, hi):
        if not (np.isfinite(lo) and np.isfinite(hi)):
            return "UNDECIDED"
        return "SATURATING" if hi < Q_RULE[0] else ("GROWING" if lo > Q_RULE[1] else "UNDECIDED")
    dec = rule(c_[0], c_[1])
    cd = ci(dq_b)
    slope = "RISING" if np.isfinite(cd[0]) and cd[0] > 0 else ("FALLING" if np.isfinite(cd[1]) and cd[1] < 0 else "UNRESOLVED")
    thr = A128 * 2 ** 0.36
    lab2 = "CONSTANT-SLOPE EXCLUDED" if np.isfinite(c256[1]) and c256[1] < thr else "NOT EXCLUDED"
    print(f"DECISION (pre-registered): {dec}; q(128 -> 256) = {fmt(q_hi_pt)} [95% {fmt(c_[0])}..{fmt(c_[1])}] (rule: SATURATING if upper < {Q_RULE[0]}, GROWING if lower > {Q_RULE[1]})")
    print(f"SECOND LABEL (pre-registered): {lab2}; S50(1/256) {fmt(s256)} [95% {fmt(c256[0])}..{fmt(c256[1])}], upper bound {'<' if c256[1] < thr else '>='} 15.28 x 2^0.36 = {thr:.2f}")
    print(f"LOCAL SLOPE (pre-registered label): {slope}; Delta q = {fmt(dq_pt)} [95% {fmt(cd[0])}..{fmt(cd[1])}]")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    t0 = time.time()
    if mode == "report":
        d = pickle.load(open(PKL, "rb"))
        report(d["res"], d["elapsed"], d["diag"])
        return
    print("hashes asserted: " + "; ".join(f"{f} {h}" for f, h in {"h14e_predictions.txt": PRED_SHA, **FILE_SHA}.items()), flush=True)
    diag = self_check()
    reload_check()
    cs = cases()
    print(f"\njobs (load average at start {', '.join(f'{x:.2f}' for x in os.getloadavg())}):", flush=True)
    for k in ORDER:
        s = cs[k]
        print(f"  {k}: L={s['L']:g} gamma={s['gamma']:g} D_eff={s['Deff']} runs={s['S'] * len(s['deltas'])} hold={s['hold']:g} ({s['hold17']} t17) seed={s['seed']} "
              f"Tk={s['flow']['Tk']:g} A={s['flow']['A']:.1f}", flush=True)
    if mode == "smoke":
        return
    deadline = t0 + GUARD_S
    res = {}
    with ProcessPoolExecutor(NWORK) as ex:
        futs = [ex.submit(C.job, (k, cs[k], deadline)) for k in ORDER]
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
