"""Model study for the NaClO3 droplet test (journals/droplet_prereg.md): a better decision rule, simulated drops, not data.

droplet_sim.py showed the prereg's h <= 3 cut calls drying-drop counting FAILED too often.  A three-agent debate
(physicist, statistician, skeptic) agreed on this rule; the full statement and the predictions are in
droplet_sim2_predictions.txt (sha256 PRED_SHA, written before the first run).  Per run of the unchanged protocol:
NS (not scorable) -> FLAT (best fit changes < 0.4 over the levels) -> RISING (trend z > 1.645) -> STEP (L < c*) -> SMOOTH,
L = max ll COUNT - max ll STEP, both blended with the fixed 3% misread.  c* comes from counting cells on a calibration
seed and is frozen before the fresh-seed runs that every prediction is scored on.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python droplet_sim2.py > droplet_sim2.log
"""
import hashlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.special import ndtr

from blind_units import exp_fit
from droplet_power import fit_h, scorable
from droplet_sim import GRID, LADDER, LG, MODELS, calibrate, label, probs

PRED_SHA = "baad47e9bd1e17d97bb9ace8e4e2fa625041cf159584d5032d7a35fee087e750"
EPS_FIT = 0.03
V50 = 3.1
LF = np.linspace(np.log(0.003), np.log(3000), 4001)          # truths live on this ln V grid
Q = np.log(1.6 / 0.2) / 7
XS = np.log(0.2) + Q * np.arange(-2, 10)                     # ln(V / V0_pilot): 2 extension slots each side
FORMS = {"first": lambda m: np.exp(-m / 2), "all": lambda m: 2 / (np.exp(np.minimum(m, 700) / 2) + 1),
         "dry": lambda m: 2 * np.exp(-m / 2) - np.exp(-m)}
MU50 = {"first": 2 * np.log(2), "all": 2 * np.log(3), "dry": -2 * np.log(1 - np.sqrt(0.5))}
AG = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 4 / 3, 1.5, 5 / 3])
SG = np.array([0.02, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6])
LMU = np.arange(-8, 5, 0.025)                                # ln mu at x = 0
CG = np.arange(-6, 4, 0.02)                                  # ln Vc
CONDS = [(cv, eps) for cv in (0.05, 0.15) for eps in (0.03, 0.06)]
MAIN = (0.05, 0.03)
CAL, EVAL = 1000, 500


def curves(x):
    """COUNT and STEP grids at ln V = x: (true P, log P_obs, log(1 - P_obs), parameter) per family."""
    out = {}
    for fam, P, par in (
            ("count", [FORMS[f](np.exp(lm + a * x)) for f in FORMS for a in AG for lm in LMU],
             [a for f in FORMS for a in AG for lm in LMU]),
            ("step", [ndtr(-(x - c) / s) for s in SG for c in CG], [s for s in SG for c in CG])):
        P = np.array(P)
        obs = np.clip(EPS_FIT + (1 - 2 * EPS_FIT) * P, 1e-12, 1 - 1e-12)
        out[fam] = (P, np.log(obs).T, np.log1p(-obs).T, np.array(par))
    return out


CUR = curves(XS)


def stats(K, N, x, cur=CUR):
    """Per run (rows of K, N at ln V = x): range of the best fit, trend z, L, a-hat, s-hat."""
    ll, arg = {}, {}
    for fam, (_, lp, l1, _) in cur.items():
        a, b = [], []
        for i in range(0, len(K), 100):
            m = K[i:i + 100] @ lp + (N - K)[i:i + 100] @ l1
            a.append(m.argmax(1)); b.append(m.max(1))
        arg[fam], ll[fam] = np.concatenate(a), np.concatenate(b)
    L = ll["count"] - ll["step"]
    Pc, Ps = cur["count"][0][arg["count"]], cur["step"][0][arg["step"]]
    Pb = np.where((L >= 0)[:, None], Pc, Ps)
    Pb = np.where(N > 0, Pb, np.nan)
    xb = (N * x).sum(1, keepdims=True) / N.sum(1, keepdims=True)
    pb = K.sum(1) / N.sum(1)
    z = (K * (x - xb)).sum(1) / np.sqrt(pb * (1 - pb) * (N * (x - xb) ** 2).sum(1) + 1e-300)
    return dict(rng=np.nanmax(Pb, 1) - np.nanmin(Pb, 1), z=z, L=L,
                a=cur["count"][3][arg["count"]], s=cur["step"][3][arg["step"]])


def protocol(rng, P, p0, eps, cv):
    """droplet_sim.protocol on truths over LF, returning (V0, k and n at the 12 slots XS) instead of a fit."""
    def draw(vl, m):
        lv = np.log(vl) + cv * rng.standard_normal(m)
        xtal = rng.random(m) >= np.interp(lv, LF, p0)
        single = (rng.random(m) < np.interp(lv, LF, P)) ^ (rng.random(m) < eps)
        return single[xtal].sum(), xtal.sum()

    k, n = np.array([draw(v, 20) for v in LADDER], float).T
    V0 = exp_fit(LADDER[n > 0], k[n > 0] / n[n > 0], n[n > 0])
    kn = {j: draw(V0 * np.exp(XS[j]), 25) for j in range(2, 10)}
    for _ in range(2):                        # extension rule
        k, n = np.array(list(kn.values()), float).T
        Pl = k[n > 0] / n[n > 0]
        j = max(kn) + 1 if Pl.min() > 0.4 else min(kn) - 1 if Pl.max() < 0.6 else None
        if j is None:
            break
        kn[j] = draw(V0 * np.exp(XS[j]), 25)
    out = np.zeros((2, len(XS)))
    for j, (a, b) in kn.items():
        out[:, j] = a, b
    return V0, out


def run_cell(task):
    """Runs of one truth under one condition.  Returns stats, NS flags, and old-rule h if asked (nan = not scored)."""
    P, p0, cv, eps, reps, seed, old = task
    rng = np.random.default_rng(seed)
    V0, KN = zip(*(protocol(rng, P, p0, eps, cv) for _ in range(reps)))
    K, N = np.array(KN)[:, 0], np.array(KN)[:, 1]
    ns = np.array([not scorable(k[n > 0], n[n > 0]) for k, n in zip(K, N)])
    h = np.full(reps, np.nan)
    if old:
        for i in np.flatnonzero(~ns):
            inc = N[i] > 0
            h[i] = fit_h(V0[i] * np.exp(XS[inc]), K[i, inc] / N[i, inc], N[i, inc])[1]
    return dict(stats(K, N, XS), ns=ns, h=h)


def classify(r, cstar):
    lab = np.where(r["L"] < cstar, "STEP", "SMOOTH").astype("<U6")
    lab[r["z"] > 1.645] = "RISING"
    lab[r["rng"] < 0.4] = "FLAT"
    lab[r["ns"]] = "NS"
    return lab


def truths():
    """(name, kind, P, p0) on LF; kind 'null' (counting, sets c*), 'step' (scored), 'report' (shown only)."""
    lv, out = LF - np.log(V50), []
    for f in FORMS:
        for a in (1 / 3, 2 / 3, 1, 4 / 3, 5 / 3):
            mu = MU50[f] * np.exp(a * lv)
            out.append((f"{f} a={a:.2f}", "null", FORMS[f](mu), np.exp(-mu) if f == "all" else 0 * mu))
    with Pool() as pool:
        cs = pool.map(calibrate, MODELS)
        ok = [(m, c) for m, c in zip(MODELS, cs) if c is not None]
        PP = pool.starmap(_model_curve, ok)
    for (m, _), (P, p0) in zip(ok, PP):
        for name, Pc, p0c in (("excl", P, p0), ("dry", p0 + (1 - p0) * P, 0 * p0)):
            out.append((f"run1 {name} {label(m)}", "null", np.interp(LF, LG, Pc), np.interp(LF, LG, p0c)))
    zero = 0 * lv
    for s in (0.1, 0.15, 0.2, 0.3, 0.45, 0.6):
        out.append((f"probit s={s}", "step", ndtr(-lv / s), zero))
    for s in (0.15, 0.3):
        out.append((f"logistic s={s}", "step", 1 / (1 + np.exp(np.clip(1.76 / s * lv, -700, 700))), zero))
    out.append(("plateau 0.3, s=0.15", "report", 0.3 + 0.7 * ndtr(-lv / 0.15), zero))
    for r in (0.3, 0.6):
        V = np.exp(LF)
        mu = MU50["all"] * np.clip((V - r * V50) / ((1 - r) * V50), 0, None)
        out.append((f"Vd = {r} x V50", "report", np.where(V < r * V50, 1.0, FORMS["all"](mu)), zero))
    return out, len(MODELS) - len(ok)


def _model_curve(m, c):
    return np.array([probs(m, c, v, 2000, 1) for v in GRID]).T     # same random numbers as droplet_sim.study


def self_check():
    V0, kn = protocol(np.random.default_rng(0), FORMS["first"](MU50["first"] * np.exp(LF - np.log(V50))),
                      0 * LF, 0.0, 0.0)
    assert 0.8 < V0 / V50 < 2.5 and (kn[1, 2:10] == 25).all()
    x = XS[2:10]
    for P, want in ((FORMS["all"](MU50["all"] * np.exp(x)), "count"), (ndtr(-x / 0.1), "step")):
        r = stats(np.round(P * 1e4)[None], np.full((1, 8), 1e4), x, curves(x))
        assert (r["L"][0] > 0) == (want == "count"), (want, r["L"])
    r = stats(np.array([[5.0] * 8]), np.full((1, 8), 25.0), x, curves(x))
    assert r["rng"][0] < 0.4 and abs(r["z"][0]) < 1e-6


def main():
    assert hashlib.sha256(Path(__file__).with_name("droplet_sim2_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    T, nov = truths()
    null = [t for t in T if t[1] == "null"]
    with Pool() as pool:
        cal = pool.map(run_cell, [(P, p0, cv, eps, CAL, [1, i, j], False)
                                  for i, (_, _, P, p0) in enumerate(null) for j, (cv, eps) in enumerate(CONDS)])
        ev = pool.map(run_cell, [(P, p0, cv, eps, EVAL, [2, i, j], (cv, eps) == MAIN)
                                 for i, (_, _, P, p0) in enumerate(T) for j, (cv, eps) in enumerate(CONDS)])
    pct = {}
    for i, t in enumerate(null):
        for j, c in enumerate(CONDS):
            r = cal[i * len(CONDS) + j]
            Ls = r["L"][(~r["ns"]) & (r["rng"] >= 0.4) & (r["z"] <= 1.645)]
            if Ls.size >= 50:
                pct[(t[0], c)] = np.percentile(Ls, 5)
    cstar = min(pct.values())
    worst = min(pct, key=pct.get)
    res = {(t[0], c): ev[i * len(CONDS) + j] for i, t in enumerate(T) for j, c in enumerate(CONDS)}
    lab = {key: classify(r, cstar) for key, r in res.items()}
    share = lambda key, v: float(np.mean(lab[key] == v))

    print("Model study, not data.  c* from counting cells on the calibration seed (1000 runs per cell and condition);")
    print("everything below is a fresh seed, 500 runs per cell.")
    print(f"c* = {cstar:.2f}, set by {worst[0].strip()} at CV {worst[1][0]}, misread {worst[1][1]}  "
          f"({len(pct)} of {len(null) * len(CONDS)} cell-conditions had >= 50 runs reaching L; {nov} run-1 model(s) no V0)")

    def row(name):
        key = (name, MAIN)
        h = res[key]["h"]
        old = f"  old rule FAILED {np.mean(h > 3):.2f}" if not np.isnan(h).all() else ""
        worst_step = max(share((name, c), "STEP") for c in CONDS)
        return (f"  {name:42} " + " ".join(f"{v} {share(key, v):.2f}" for v in ("STEP", "SMOOTH", "FLAT", "RISING", "NS"))
                + f" | worst STEP over conditions {worst_step:.2f}{old}")

    print(f"\nCounting cells (should not be STEP).  Shares at CV {MAIN[0]}, misread {MAIN[1]}:")
    for t in null:
        print(row(t[0]))
    print(f"\nStep and reported truths.  Shares at CV {MAIN[0]}, misread {MAIN[1]}; then STEP share per condition"
          " and median s-hat of STEP calls:")
    for t in T:
        if t[1] != "null":
            per = " ".join(f"{c[0]}/{c[1]}:{share((t[0], c), 'STEP'):.2f}" for c in CONDS)
            r = res[(t[0], MAIN)]
            sh = r["s"][lab[(t[0], MAIN)] == "STEP"]
            print(row(t[0]) + f"\n      STEP by CV/misread {per}  median s-hat {np.median(sh) if sh.size else np.nan:.3f}")

    probit = [s for s in (0.1, 0.15, 0.2, 0.3, 0.45, 0.6)]
    ok_s = [s for s in probit if min(share((f"probit s={s}", c), "STEP") for c in CONDS) >= 0.95]
    sstar = max(ok_s) if ok_s else float("nan")
    print(f"\ns* = {sstar} (largest probit s called STEP in >= 95% of runs in every condition)")
    rising = [f"run1 excl {label(m)}" for m in (("line", 3, 0.01, 1.0), ("line", 3, 0.1, 10.0), ("bulk", 3, 0.01, 1.0))]
    rising = [n for n in rising if (n, MAIN) in res]
    dry53 = ("dry a=1.67", MAIN)
    checks = [
        ("S1", max(share((t[0], c), "STEP") for t in null for c in CONDS) <= 0.07,
         f"max STEP over counting cells {max(share((t[0], c), 'STEP') for t in null for c in CONDS):.3f}"),
        ("S2", min(share((f"probit s={s}", c), "STEP") for s in (0.1, 0.15) for c in CONDS) >= 0.90,
         f"min {min(share((f'probit s={s}', c), 'STEP') for s in (0.1, 0.15) for c in CONDS):.2f}"),
        ("S3", share(("probit s=0.6", MAIN), "STEP") <= 0.50, f"{share(('probit s=0.6', MAIN), 'STEP'):.2f}"),
        ("S4", share(("logistic s=0.15", MAIN), "STEP") >= 0.90, f"{share(('logistic s=0.15', MAIN), 'STEP'):.2f}"),
        ("S5", len(rising) == 3 and max(share((n, c), "STEP") for n in rising for c in CONDS) <= 0.05,
         f"{len(rising)} curves, max STEP {max(share((n, c), 'STEP') for n in rising for c in CONDS):.2f}"),
        ("S6", np.mean(res[dry53]["h"] > 3) > 0.5 and share(dry53, "STEP") <= 0.07,
         f"old FAILED {np.mean(res[dry53]['h'] > 3):.2f}, new STEP {share(dry53, 'STEP'):.2f}"),
        ("S7", 0.1 <= sstar <= 0.2, f"s* {sstar}"),
        ("S8", -4 <= cstar <= -0.5, f"c* {cstar:.2f}"),
    ]
    print("\nPredictions:")
    for name, held, detail in checks:
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")


if __name__ == "__main__":
    main()
