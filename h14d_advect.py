"""h14d_advect: does the eddy slowdown S50 at l/w = 10 keep growing as (D_eff / D_mol)^q past D_mol/D_eff = 1/64, or saturate?
(open item 14 / H14, journals/ija_draft.md:252; follows h14c_advect.py; 2026-10-21.)

Question. h14c measured S50 at l/w 10 at 1/8, 1/32, 1/64 = 5.76, 10.27, 13.21 (h14c_advect.log:83, :79, :77), q = 0.40 over the three,
local exponents 0.42 (8 -> 32) and 0.36 (32 -> 64). A steady cellular flow gives q -> 0.5, chaotic stirring q -> 0; q = 0.40 carried to the
real ratio 1e9 gives S50 ~ 1e4. h14d adds measured points at 1/16 and 1/128 (same model, same flow family, same per-case gamma rule), adds
runs at 1/64 (pooled with the 42 of h14c: precision of the decision rule), and scores the pre-registered rule for q(64 -> 128):
SATURATING (95% upper bound < 0.2), GROWING (95% lower bound > 0.3), else UNDECIDED. Bands P1-P10 and the rule: h14d_predictions.txt.

Nothing existing is edited: h14b_advect.py, h14c_advect.py (model, Stat, pool, job, kick_damping, self_check) and advect2d.py are imported;
the three cases below are the h14c flow_for / gamma rule at D_eff = 16, 64, 128. The h14c raw pickle is read for the 1/8, 1/32, 1/64 cells.

PREDICTIONS: h14d_predictions.txt, sha256 PRED_SHA (asserted on import). Imported files are hash-asserted too (FILE_SHA).
Run: cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project \
       --with numpy --with scipy --with numba --with matplotlib python h14d_advect.py > h14d_advect.log 2> h14d_advect.err
     ("report" re-reads h14d_advect_raw.pkl; "smoke" = self-check and h14c reproduction only; 3 worker processes)
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

import advect2d as AD
import h14b_advect as H
import h14c_advect as C          # asserts the h14c predictions sha256 on import

HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
PRED_SHA = "5da5f82b7a392c106d9eaee9406d1ebf0d3948a7fc0b4bae1e322a1039f63a0d"
FILE_SHA = {"advect2d.py": "7d84a2e6e8c05c4d807f945a3f1cd3db40115ab87f32ae053e499558e145cf2e",
            "h14b_advect.py": "9308bc499797197a0050933ec7d87a7532ab96ec60aeccf8d5611f6edc5a3d15",
            "h14c_advect.py": "aac6365eb4b8b06a7ac4c2e93e70e36f354888534c3f9257e187f06146c6f35d",
            "h14c_predictions.txt": "522fb6ea19a961bc877d182e6e9f668d96d9f91efc41cc5258d96006b1ce79ea",
            "h14c_advect_raw.pkl": "1df85e4f81fdab7a70b097ce1d40947f0f57bd1ca2666941d012f3b659b2fc36"}
assert sha("h14d_predictions.txt") == PRED_SHA
for _f, _h in FILE_SHA.items():
    assert sha(_f) == _h, f"{_f} changed since h14d was written"

PKL = os.path.join(HERE, "h14d_advect_raw.pkl")
GUARD_S = 150 * 60
NWORK = 3
RATIO = {"R8_10": 8, "R16_10": 16, "V32_10": 32, "V64_10": 64, "V128_10": 128}
Q_RULE = (0.2, 0.3)                                            # SATURATING if q_hi < 0.2; GROWING if q_lo > 0.3


def cases():
    c = {}

    def add(name, Deff, S, hold17, seed):
        L, fl = C.flow_for(Deff, 10)
        gamma = round(min(0.01, 1 / (5 * fl["Tk"])), 4)       # 5 kicks between lam = 0 and the end of the sweep
        T = C.t17(L, Deff)
        c[name] = dict(L=float(L), n=1, h=1.0, dt=0.2, S=S, D=C.D_MOL, flow=fl, noise="x", seed=seed, gamma=gamma, Deff=Deff, lw=10,
                       deltas=(0.5, 1.0, 2.0), hold17=hold17, hold=float(np.ceil(hold17 * T / 2) * 2),
                       cad=max(2.0, round(T / 60 / 0.2) * 0.2), name=name)
    for j in range(16):
        add(f"V128_10#{j}", 128, 7, 60, 3301 + j)
    for j in range(16):
        add(f"V64_10n#{j}", 64, 7, 48, 3321 + j)
    for j in range(4):
        add(f"R16_10#{j}", 16, 16, 36, 3341 + j)
    return c


ORDER = list(cases())                                          # 1/128 first, then 1/64, then 1/16


# ----------------------------------------------------------------------------------------------- self-check and diagnostics
def lyapunov(kappa, nk=20000, nv=200, seed=11):
    """Lyapunov exponent per kick of the tangent map. The shear phase is redrawn every kick, so the tangent matrices are independent:
    [[1, kappa cos th], [0, 1]] and [[1, 0], [kappa cos th, 1]] alternately, th uniform."""
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(nv, 2))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    tot = 0.0
    for j in range(nk):
        c = kappa * np.cos(rng.uniform(0, 2 * np.pi, nv))
        if j % 2 == 0:
            v[:, 0] += c * v[:, 1]
        else:
            v[:, 1] += c * v[:, 0]
        n = np.linalg.norm(v, axis=1)
        tot += np.log(n).sum()
        v /= n[:, None]
    return tot / (nk * nv)


def old_cells():
    return pickle.load(open(C.PKL, "rb"))["res"]


def self_check():
    t0 = time.time()
    C.self_check()                                             # h14c's own: nk fix, g from the actual gamma, tiny job()/pool(), 15 worlds
    cs = cases()
    seeds = [s["seed"] for s in cs.values()]
    assert len(set(seeds)) == len(seeds) and not set(seeds) & set(range(3101, 3133)), "seeds"
    old = old_cells()
    diag = dict(tracer={}, damp={}, lyap={}, flow={})
    for D in (16, 64, 128):
        sp = next(s for s in cs.values() if s["Deff"] == D)
        fl = sp["flow"]
        assert abs(fl["A"] ** 2 / (8 * fl["Tk"]) - (D - 1)) < 1e-9
        N = sp["L"] / sp["h"]
        assert abs(N - round(N)) < 1e-9 and sp["D"] * sp["dt"] / sp["h"] ** 2 <= 0.25 + 1e-9
        assert abs(1 / sp["gamma"] / fl["Tk"] - 5) < 0.7
        AD.L = sp["L"]
        d = float(np.mean([AD.tracer_lagrangian(fl, s) for s in (1, 7)]))
        kap = fl["A"] / fl["l"]
        diag["tracer"][D] = (d, D - 1.0)
        diag["damp"][D] = C.kick_damping(fl["Tk"])
        diag["lyap"][D] = lyapunov(kap)
        diag["flow"][D] = dict(kappa=kap, gamma=sp["gamma"], kicks_pre=1 / sp["gamma"] / fl["Tk"], L=sp["L"], Tk=fl["Tk"], A=fl["A"], l=fl["l"])
    o = old["V64_10#0"]["spec"]                                # the new 1/64 jobs are the h14c design (pooling is legitimate)
    n = cs["V64_10n#0"]
    assert o["flow"] == n["flow"] and o["gamma"] == n["gamma"] and o["L"] == n["L"] and o["hold"] == n["hold"] and o["cad"] == n["cad"] \
        and o["deltas"] == n["deltas"] and o["h"] == n["h"] and o["dt"] == n["dt"] and o["noise"] == n["noise"], (o, n)
    print(f"h14d self-check passed ({time.time() - t0:.0f}s): h14c self_check; seeds unique and outside 3101-3132; D_eddy = D_eff - 1 by construction; "
          f"Euler stability; kicks_pre within 5 +- 0.7; new 1/64 jobs have the h14c V64_10 flow, gamma, L, hold, cadence, Deltas; "
          f"{len(cs)} jobs / {sum(s['S'] * len(s['deltas']) for s in cs.values())} runs", flush=True)
    return diag


def cell_sets(res):
    """name -> (parts, bootstrap-seed name). Old cells are reloaded from the h14c pickle with the h14c seeds."""
    old = old_cells()
    g = lambda d, base: [v for k, v in d.items() if k.split("#")[0] == base]
    return {"R8_10": (g(old, "R8_10"), "R8_10"), "R16_10": (g(res, "R16_10"), "R16_10"), "V32_10": (g(old, "V32_10"), "V32_10"),
            "V64_10o": (g(old, "V64_10"), "V64_10"), "V64_10n": (g(res, "V64_10n"), "V64_10n"),
            "V64_10": (g(old, "V64_10") + g(res, "V64_10n"), "V64_10"), "V128_10": (g(res, "V128_10"), "V128_10")}


def stats(res):
    cs = cell_sets(res)
    P = {k: C.pool(v) for k, (v, _) in cs.items()}
    return P, {k: C.Stat(P[k], seed=sum(map(ord, cs[k][1]))) for k in P}


def repro_check():
    """The reloaded h14c cells give the printed h14c values (5.76, 10.27, 13.21)."""
    old = old_cells()
    g = lambda base: [v for k, v in old.items() if k.split("#")[0] == base]
    for base, want in (("R8_10", "5.76"), ("V32_10", "10.27"), ("V64_10", "13.21")):
        s = C.Stat(C.pool(g(base)), seed=sum(map(ord, base)))
        assert f"{s.pt[0]:.2f}" == want, (base, s.pt[0], want)
    print("h14c reproduction: R8_10 5.76, V32_10 10.27, V64_10 (42 runs) 13.21 reloaded from h14c_advect_raw.pkl with the h14c bootstrap seeds", flush=True)


# ----------------------------------------------------------------------------------------------- report
def report(res, elapsed, diag):
    P, st = stats(res)
    fmt, line = H.fmt, H.line
    sv = lambda k, q=0: st[k].pt[q]
    S = lambda k: st[k].bs[:, 0]
    print("hashes: " + "; ".join(f"{f} {h}" for f, h in {"h14d_predictions.txt": PRED_SHA, **FILE_SHA}.items()))
    print("\nREALISED FLOWS (new cases; AD.L set to the box, mean of 2 seeds) and regime diagnostics")
    for D in (16, 64, 128):
        f, (d, nom) = diag["flow"][D], diag["tracer"][D]
        lam = diag["lyap"][D]
        print(f"  D_eff {D:3d} L {f['L']:4.0f} Tk {f['Tk']:5.1f} A {f['A']:6.1f} eddy l {f['l']:6.1f}: kappa {f['kappa']:.3f} kicks_pre {f['kicks_pre']:.2f} gamma {f['gamma']:g}; "
              f"D_eddy measured {d:.2f} nominal {nom:.0f} (ratio {d / nom:.3f}); D_num/D_mol {diag['damp'][D]:.4f}; "
              f"Lyapunov exponent {lam:.3f} per kick ({lam / f['Tk']:.4f} per time unit); rough mixing time ln(l / 1 cell)/lambda ~ {np.log(f['l']) / lam:.1f} kicks")

    print("\nCASES  (L, ratio, l/w, Tk, A, gamma, t17, runs, walled at freeze-out, hold, cens at end, cut, cpu s)")
    for k in ("R8_10", "R16_10", "V32_10", "V64_10o", "V64_10n", "V64_10", "V128_10"):
        s, fl = P[k]["spec"], P[k]["spec"]["flow"]
        print(f"  {k:8s} L={s['L']:g} Deff={s['Deff']} 1/{s['Deff']} l/w={fl['l'] / np.sqrt(s['Deff']):.2f} Tk={fl['Tk']:g} A={fl['A']:.1f} g={s['gamma']:g} t17={C.t17(s['L'], s['Deff']):.1f} "
              f"runs={P[k]['runs']} walled={st[k].nw} ({st[k].pt[3]:.3f}) hold={s['hold']:g} cens_at_end={int(np.isinf(P[k]['t_single'][P[k]['walled_at_t0']]).sum())} "
              f"cut={P[k]['cut']} cpu={P[k]['elapsed']:.0f}")

    print("\nSLOWDOWN at l/w 10  (time since the end of the sweep / t17; Kaplan-Meier over runs walled at freeze-out; bootstrap over runs, 2000x)")
    for k in ("R8_10", "R16_10", "V32_10", "V64_10o", "V64_10n", "V64_10", "V128_10"):
        s = st[k]
        s50 = s.pt[0]
        print(f"  {k:8s} {line(s, 0, 'S50')} | {line(s, 1, 'S80')} | {line(s, 2, 'S_half')}")
        print(f"           walled {s.pt[3]:.3f}  mean domains at t0 {P[k]['ndom0'].mean():.2f}  kicks to clear (S50 t17 / Tk) "
              f"{s50 * s.t17 / P[k]['spec']['flow']['Tk']:.0f}  cleared by 12 / 16 / 24 t17: " + " ".join(f"{C.km_cleared(P[k], s, m):.2f}" for m in (12, 16, 24)))

    def ci(b):
        fin = np.isfinite(b)
        return np.percentile(b[fin], [2.5, 97.5]) if (1 - fin.mean()) <= 0.05 and fin.sum() > 100 else np.array([np.nan, np.nan])

    def q2(ka, kb):
        l = np.log(RATIO[kb] / RATIO[ka])
        return np.log(sv(kb) / sv(ka)) / l, np.log(S(kb) / S(ka)) / l

    def ls(keys):
        x = np.array([RATIO[k] for k in keys], float)
        pt = H.fit_p(x, np.array([sv(k) for k in keys]))
        return pt, np.array([H.fit_p(x, v) for v in np.stack([S(k) for k in keys], 1)])

    def sf(pt, b):
        c = ci(b)
        return f"{fmt(pt)} [95% {fmt(c[0])}..{fmt(c[1])}]"

    print("\nLOCAL ratio exponents q(a -> b) = ln(S50_b / S50_a) / ln(r_b / r_a)  (h14c: 8 -> 32 0.42 [0.32..0.56]; 32 -> 64 0.36 [0.15..0.53])")
    ladder = ("R8_10", "R16_10", "V32_10", "V64_10", "V128_10")
    for ka, kb in zip(ladder[:-1], ladder[1:]):
        print(f"  1/{RATIO[ka]:<3d} -> 1/{RATIO[kb]:<3d}: q = {sf(*q2(ka, kb))}")
    print(f"  1/8 -> 1/64 (endpoints): q = {sf(*q2('R8_10', 'V64_10'))}")
    q_ref, b_ref = ls(("R8_10", "V32_10", "V64_10"))
    q_all, b_all = ls(ladder)
    print(f"  q(8 -> 64), least squares over 8, 32, 64 (h14c definition, pooled 1/64): {sf(q_ref, b_ref)}   [h14c: 0.40 [0.33..0.50]]")
    print(f"  q over all five ratios (8, 16, 32, 64, 128), least squares: {sf(q_all, b_all)}")
    q_hi_pt, b_hi = q2("V64_10", "V128_10")
    dq_pt, dq_b = q_hi_pt - q_ref, b_hi - b_ref
    print(f"  q(64 -> 128) = {sf(q_hi_pt, b_hi)}; q(8 -> 64) = {sf(q_ref, b_ref)}; Delta q = q(64 -> 128) - q(8 -> 64) = {sf(dq_pt, dq_b)} (paired over the same replicates)")
    c_ = ci(b_hi)
    print(f"  q(64 -> 128) 95% interval excludes 0.5 (cellular): {bool(c_[1] < 0.5) if np.isfinite(c_[1]) else 'n/a'}; excludes 0 (chaotic limit): {bool(c_[0] > 0) if np.isfinite(c_[0]) else 'n/a'}")
    x5 = np.log(np.array([RATIO[k] for k in ladder], float))
    quad = lambda v: np.polyfit(x5, np.log(v), 2)[0] if np.isfinite(v).all() else np.nan
    cq = quad(np.array([sv(k) for k in ladder]))
    cq_b = np.array([quad(v) for v in np.stack([S(k) for k in ladder], 1)])
    print(f"  curvature c of ln S50 = a + q ln r + c (ln r)^2 over the five ratios: {sf(cq, cq_b)}")

    def rule(lo, hi):
        if not (np.isfinite(lo) and np.isfinite(hi)):
            return "UNDECIDED"
        return "SATURATING" if hi < Q_RULE[0] else ("GROWING" if lo > Q_RULE[1] else "UNDECIDED")
    print("  sensitivity of q(64 -> 128) to the 1/64 cell (rule outcome in brackets):")
    sens = {}
    for k64 in ("V64_10", "V64_10o", "V64_10n"):
        l = np.log(2.0)
        qp, qb = np.log(sv("V128_10") / sv(k64)) / l, np.log(S("V128_10") / S(k64)) / l
        cc = ci(qb)
        sens[k64] = rule(cc[0], cc[1])
        print(f"    1/64 cell {k64:8s} ({P[k64]['runs']} runs, S50 {fmt(sv(k64))}): q = {fmt(qp)} [95% {fmt(cc[0])}..{fmt(cc[1])}] -> {sens[k64]}")

    print("\nEXTRAPOLATION (NOT a measurement): S50 at the real ratio 1e9 and the ratio at which S50 reaches f, anchored at the 1/128 point")
    s128 = sv("V128_10")
    cases_q = (("q(64->128) 2.5% bound", c_[0]), ("q(64->128) point", q_hi_pt), ("q(64->128) 97.5% bound", c_[1]), ("q(8->64) point (h14c slope)", q_ref), ("constant from 1/128", 0.0))
    for lab, q in cases_q:
        if not np.isfinite(q):
            continue
        row = [f"S50(1e9) = {s128 * (1e9 / 128) ** max(q, 0):.3g}"]
        for f_ in (100, 200, 500, 1000, 10000):
            row.append(f"f={f_}: " + (f"{128 * (f_ / s128) ** (1 / q):.2e}" if q > 0 and f_ > s128 else ("below S50(128)" if f_ <= s128 else "never")))
        print(f"  {lab:30s} q {q:5.2f}: " + "; ".join(row))

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
    print("\nPREDICTIONS (bands from h14d_predictions.txt)")
    verdict(1, AND(band("R16_10", 0, 6.3, 9.5), st["R16_10"].pt[3] >= 0.9), f"R16_10 S50 {fmt(sv('R16_10'))} in [6.3, 9.5] (guess 7.7), walled {st['R16_10'].pt[3]:.3f}")
    verdict(2, AND(band("V128_10", 0, 12.5, 20.5), band("V128_10", 1, 16, 30), st["V128_10"].pt[3] >= 0.9),
            f"V128_10 S50 {fmt(sv('V128_10'))} in [12.5, 20.5] (guess 16.3), S80 {fmt(sv('V128_10', 1))} in [16, 30] (guess 21), walled {st['V128_10'].pt[3]:.3f}")
    verdict(3, inband(q_hi_pt, 0.10, 0.50), f"q(64 -> 128) = {fmt(q_hi_pt)} in [0.10, 0.50] (guess 0.30)")
    verdict(4, inband(dq_pt, -0.40, 0.05), f"Delta q = {fmt(dq_pt)} in [-0.40, 0.05] (guess -0.10); q(8 -> 64) = {fmt(q_ref)}")
    verdict(5, AND(*[sv(kb) > sv(ka) for ka, kb in zip(ladder[:-1], ladder[1:])]), "S50 ladder " + " < ".join(fmt(sv(k)) for k in ladder))
    r6 = sv("V64_10n") / sv("V64_10o")
    verdict(6, inband(r6, 0.85, 1.18), f"S50(V64_10n {fmt(sv('V64_10n'))}) / S50(V64_10 old {fmt(sv('V64_10o'))}) = {fmt(r6)} in [0.85, 1.18]")
    ok7 = [abs(diag["tracer"][D][0] / diag["tracer"][D][1] - 1) <= 0.05 and diag["damp"][D] <= 0.02 and 1.45 <= diag["flow"][D]["kappa"] <= 1.49
           and 4.5 <= diag["flow"][D]["kicks_pre"] <= 5.5 for D in (16, 64, 128)]
    verdict(7, all(ok7), "flows 1/16, 1/64, 1/128: " + "; ".join(f"D_eddy ratio {diag['tracer'][D][0] / diag['tracer'][D][1]:.3f}, D_num/D_mol {diag['damp'][D]:.4f}, kappa {diag['flow'][D]['kappa']:.3f}, "
                                                                f"kicks_pre {diag['flow'][D]['kicks_pre']:.2f}" for D in (16, 64, 128)))
    verdict(8, all(diag["lyap"][D] > 0.1 for D in diag["lyap"]), "Lyapunov exponent per kick " + ", ".join(f"{diag['lyap'][D]:.3f}" for D in (16, 64, 128)) + " (> 0.1, guess 0.3)")
    verdict(9, inband(cq, -np.inf, 0.0), f"curvature c = {sf(cq, cq_b)} <= 0 (guess -0.02)")
    cut = any(v["cut"] for v in res.values())
    nruns = sum(v["runs"] for v in res.values())
    verdict(10, elapsed < GUARD_S and not cut and nruns == 864, f"total {elapsed:.0f} s ({NWORK} workers, guard {GUARD_S} s), cut by guard {cut}, new runs {nruns} (864)")
    held = [n for n, ok in out if ok is True]
    print("\nSUMMARY  held: " + ", ".join(f"P{n}" for n in held) + " | failed: " + ", ".join(f"P{n}" for n, ok in out if ok is False)
          + " | undecided: " + ", ".join(f"P{n}" for n, ok in out if isinstance(ok, str)))
    dec = rule(c_[0], c_[1])
    cd = ci(dq_b)
    slope = "RISING" if np.isfinite(cd[0]) and cd[0] > 0 else ("FALLING" if np.isfinite(cd[1]) and cd[1] < 0 else "UNRESOLVED")
    print(f"DECISION (pre-registered): {dec}; q(64 -> 128) = {fmt(q_hi_pt)} [95% {fmt(c_[0])}..{fmt(c_[1])}] (rule: SATURATING if upper < {Q_RULE[0]}, GROWING if lower > {Q_RULE[1]}); "
          f"sensitivity pooled / old-only / new-only 1/64: {sens['V64_10']} / {sens['V64_10o']} / {sens['V64_10n']}")
    print(f"LOCAL SLOPE (pre-registered label): {slope}; Delta q = {fmt(dq_pt)} [95% {fmt(cd[0])}..{fmt(cd[1])}]")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    t0 = time.time()
    if mode == "report":
        d = pickle.load(open(PKL, "rb"))
        report(d["res"], d["elapsed"], d["diag"])
        return
    print("hashes asserted: " + "; ".join(f"{f} {h}" for f, h in {"h14d_predictions.txt": PRED_SHA, **FILE_SHA}.items()), flush=True)
    diag = self_check()
    repro_check()
    if mode == "smoke":
        return
    cs = cases()
    print("\njobs (1/128 first):", flush=True)
    for k in ORDER:
        s = cs[k]
        print(f"  {k}: L={s['L']:g} gamma={s['gamma']:g} D_eff={s['Deff']} runs={s['S'] * len(s['deltas'])} hold={s['hold']:g} ({s['hold17']} t17) seed={s['seed']} "
              f"Tk={s['flow']['Tk']:g} A={s['flow']['A']:.1f}", flush=True)
    deadline = t0 + GUARD_S
    res = {}
    with ProcessPoolExecutor(NWORK) as ex:
        futs = [ex.submit(C.job, (k, cs[k], deadline)) for k in ORDER]
        for fu in as_completed(futs):
            k, r = fu.result()
            res[k] = r
            print(f"  job done: {k} [{r['elapsed']:.0f}s] walled {int(r['walled_at_t0'].sum())}/{r['runs']} cleared {int(np.isfinite(r['t_single'][r['walled_at_t0']]).sum())} "
                  f"cut={r['cut']} (elapsed {time.time() - t0:.0f}s)", flush=True)
            pickle.dump(dict(res=res, elapsed=time.time() - t0, diag=diag, partial=True), open(PKL + ".part", "wb"))
    elapsed = time.time() - t0
    pickle.dump(dict(res=res, elapsed=elapsed, diag=diag), open(PKL, "wb"))
    report(res, elapsed, diag)


if __name__ == "__main__":
    main()
