"""Model study for the NaClO3 droplet test (journals/droplet_prereg.md): the misread rate measured on control drops.
Simulated drops, not data.

droplet_sim2.py failed S2 and S7 because misreads above the fit's fixed 3% hide sharp steps, so the prereg now counts
control drops of known outcome.  This compares two policies over true misread rates 0-10% and 64/128/192 controls:
A, the prereg as written (fixed 3% fit, c* = -1.61, NOT SCORABLE when over 6% of controls are read wrong), and B,
refitting both curves with the measured rate (c*_B calibrated per control count, as droplet_sim2 set c*).
Predictions: misread_sim_predictions.txt (sha256 PRED_SHA, written before the first run).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python misread_sim.py > misread_sim.log
"""
import hashlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from droplet_power import scorable
from droplet_sim2 import CUR, XS, classify, protocol, stats, truths

PRED_SHA = "191ea36486e6577ab56daa77f3661cc7e4b0624737f42d4a854a195706ec89fb"
C_A = -1.61                                   # droplet_sim2.log (calibration seed), frozen
CONDS = [(cv, e) for cv in (0.05, 0.15) for e in (0.0, 0.02, 0.04, 0.06, 0.08, 0.10)]
NC = (64, 128, 192)
STEPS = ("probit s=0.1", "probit s=0.15", "probit s=0.2", "probit s=0.3")
CAL, EVAL = 1000, 500


def cur_eps(e):
    """droplet_sim2.CUR with the misread blend e in place of 3%."""
    out = {}
    for fam, (P, _, _, par) in CUR.items():
        obs = np.clip(e + (1 - 2 * e) * P, 1e-12, 1 - 1e-12)
        out[fam] = (P, np.log(obs).T, np.log1p(-obs).T, par)
    return out


def refit(K, N, eh):
    """droplet_sim2.stats per run, each run fitted with its own measured misread rate eh."""
    r = {}
    for e in np.unique(eh):
        m = eh == e
        for key, v in stats(K[m], N[m], XS, cur_eps(e)).items():
            r.setdefault(key, np.empty(len(K)))[m] = v
    return r


def run_cell(task):
    """Runs of one truth under one condition: policy A stats and control-tier flags, policy B stats, per n_c."""
    P, p0, cv, eps, reps, seed, want_a = task
    rng = np.random.default_rng(seed)
    KN = np.array([protocol(rng, P, p0, eps, cv)[1] for _ in range(reps)])
    K, N = KN[:, 0], KN[:, 1]
    ns = np.array([not scorable(k[n > 0], n[n > 0]) for k, n in zip(K, N)])
    wrong = rng.random((reps, max(NC))) < eps          # controls read wrong; 64 and 128 are the first of the 192
    out = {"A": dict(stats(K, N, XS), ns=ns)} if want_a else {}
    for nc in NC:
        k = wrong[:, :nc].sum(1)
        out[("tier", nc)] = k / nc > 0.06
        out[("B", nc)] = dict(refit(K, N, (k + 0.5) / (nc + 1)), ns=ns)
    return out


def self_check():
    c = cur_eps(0.03)
    assert all(np.allclose(c[f][1], CUR[f][1]) and np.allclose(c[f][2], CUR[f][2]) for f in CUR)
    N = np.zeros((3, len(XS)))
    N[:, 2:10] = 25
    K = np.pad(np.random.default_rng(0).binomial(25, np.linspace(0.95, 0.05, 8), (3, 8)).astype(float), ((0, 0), (2, 2)))
    a, b = stats(K, N, XS), refit(K, N, np.full(3, 0.03))
    assert all(np.allclose(a[key], b[key]) for key in a)
    b = refit(K, N, np.array([0.03, 0.08, 0.03]))
    assert np.allclose(b["L"][[0, 2]], a["L"][[0, 2]]) and not np.isclose(b["L"][1], a["L"][1])


def main():
    assert hashlib.sha256(Path(__file__).with_name("misread_sim_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    T, _ = truths()
    null = [t for t in T if t[1] == "null"]
    ev_t = null + [t for t in T if t[0] in STEPS]
    with Pool() as pool:
        cal = pool.map(run_cell, [(P, p0, cv, e, CAL, [4, i, j], False)
                                  for i, (_, _, P, p0) in enumerate(null) for j, (cv, e) in enumerate(CONDS)])
        ev = pool.map(run_cell, [(P, p0, cv, e, EVAL, [5, i, j], True)
                                 for i, (_, _, P, p0) in enumerate(ev_t) for j, (cv, e) in enumerate(CONDS)])
    cstar, worst = {}, {}
    for nc in NC:
        pct = {}
        for i, t in enumerate(null):
            for j, c in enumerate(CONDS):
                r = cal[i * len(CONDS) + j][("B", nc)]
                Ls = r["L"][(~r["ns"]) & (r["rng"] >= 0.4) & (r["z"] <= 1.645)]
                if Ls.size >= 50:
                    pct[(t[0], c)] = np.percentile(Ls, 5)
        cstar[nc], worst[nc] = min(pct.values()), min(pct, key=pct.get)
    POL = [("A", nc) for nc in NC] + [("B", nc) for nc in NC]
    lab = {}
    for i, t in enumerate(ev_t):
        for j, c in enumerate(CONDS):
            o = ev[i * len(CONDS) + j]
            for nc in NC:
                lab[(t[0], c, ("A", nc))] = classify(dict(o["A"], ns=o["A"]["ns"] | o[("tier", nc)]), C_A)
                lab[(t[0], c, ("B", nc))] = classify(o[("B", nc)], cstar[nc])
    share = lambda name, c, pol, v="STEP": float(np.mean(lab[(name, c, pol)] == v))

    print("Model study, not data.  Policy A: the prereg as written (fixed 3% fit, c* = -1.61, NOT SCORABLE when more")
    print("than 6% of controls are read wrong).  Policy B: refit with the measured rate (k + 0.5) / (n_c + 1).")
    print(f"c*_B from {len(null)} counting cells x {len(CONDS)} conditions, {CAL} runs each on the calibration seed;"
          f" everything below is a fresh seed, {EVAL} runs per cell.")
    for nc in NC:
        print(f"  c*_B({nc}) = {cstar[nc]:.2f}, set by {worst[nc][0]} at CV {worst[nc][1][0]}, misread {worst[nc][1][1]}")
    head = "  CV   misread " + " ".join(f"{p}{nc:>4}" for p, nc in POL)
    print("\nCounting cells: highest STEP share over all cells (should be <= 0.07)")
    print(head)
    for c in CONDS:
        print(f"  {c[0]:.2f} {c[1]:.2f}    " + " ".join(f"{max(share(t[0], c, p) for t in null):5.3f}" for p in POL))
    for p in POL:
        top = max(((share(t[0], c, p), t[0], c) for t in null for c in CONDS))
        print(f"  {p[0]}{p[1]}: max {top[0]:.3f} ({top[1]}, CV {top[2][0]}, misread {top[2][1]})")
    for s in STEPS:
        print(f"\n{s}: STEP share")
        print(head)
        for c in CONDS:
            print(f"  {c[0]:.2f} {c[1]:.2f}    " + " ".join(f"{share(s, c, p):5.2f}" for p in POL))
    print("\nprobit s=0.15: NOT SCORABLE share (scorability plus, for A, the control tier)")
    print(head)
    for c in CONDS:
        print(f"  {c[0]:.2f} {c[1]:.2f}    " + " ".join(f"{share('probit s=0.15', c, p, 'NS'):5.2f}" for p in POL))

    s15 = "probit s=0.15"
    low = [c for c in CONDS if c[1] <= 0.06]
    wB = {nc: min(share(s15, c, ("B", nc)) for c in low) for nc in NC}
    m1 = max(share(t[0], c, ("A", 64)) for t in null for c in CONDS)
    m3 = max(share(t[0], c, ("B", nc)) for t in null for c in CONDS for nc in NC)
    m4 = min(share(s, c, ("B", 64)) for s in STEPS[:2] for c in low)
    m5 = [(c, share(s15, c, ("B", 64)), share(s15, c, ("A", 64))) for c in CONDS if c[1] in (0.04, 0.06)]
    checks = [
        ("M1", m1 <= 0.07, f"max STEP over counting cells, A 64: {m1:.3f}"),
        ("M2", all(share(s15, c, ("A", 64)) < 0.5 for c in CONDS if c[1] == 0.06),
         "A 64 at misread 0.06: " + ", ".join(f"CV {c[0]} {share(s15, c, ('A', 64)):.2f}" for c in CONDS if c[1] == 0.06)),
        ("M3", m3 <= 0.07, f"max STEP over counting cells, B any n_c: {m3:.3f}"),
        ("M4", m4 >= 0.90, f"min STEP for s <= 0.15, B 64, misread <= 0.06: {m4:.2f}"),
        ("M5", all(b > a for _, b, a in m5), " ".join(f"CV {c[0]}/{c[1]}: B {b:.2f} A {a:.2f}" for c, b, a in m5)),
        ("M6", wB[192] - wB[64] <= 0.05, f"worst s=0.15 STEP, B: 64 {wB[64]:.2f}, 128 {wB[128]:.2f}, 192 {wB[192]:.2f}"),
        ("M7", all(-4 <= v <= -0.5 for v in cstar.values()), " ".join(f"{nc}: {v:.2f}" for nc, v in cstar.items())),
    ]
    print("\nPredictions:")
    for name, held, detail in checks:
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")


if __name__ == "__main__":
    main()
