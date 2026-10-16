"""How many amino acids does an abiotic amplifier act on?  Follow-up to life_test.py after a review on 2026-10-16.  Imports
life_test.py (its gen, dataset, fit, pool, load, rates) and edits nothing.

WHY.  life_test.py let the amplifier (scenario A) act on all seven protein amino acids with life's hand and life's share prior.
A review found that two of its readings follow from that choice, not from the test: "life is never confirmed" and "Davila's
rule (life if >= 4 of 7 amino acids show ee >= 5%) has precision 0.52".  Davila's own premise is the opposite: an abiotic excess
would be found in only 2-3 amino acids.  This script makes the amplifier's breadth the variable.  In scenario A, after
life_test.gen draws the ee vector, only `b` randomly chosen protein amino acids keep their ee (the others are set to 0; Iva is
untouched), with b drawn uniformly from a breadth range.  Scenarios B (life) and C (plain abiotic) are unchanged.
 Breadth rows: 2-3 (Davila's premise), 4-5, 6, and 7 (= life_test.py, regression).
 Reported per row (n = 8 with Iva and n = 7 protein only, sigma 2%, equal priors, same classifier as life_test.py):
   acc, acc_V_AB (as defined in life_test.fit), share of true life confidently called life (posterior >= 0.9), and
   Davila's rule on the n = 7 protein panel: precision and recall of its life call.
 Also corrects a misreading in life_test.py's RESULTS block: its log line "Largest drop in acc_V_AB: Iva only (0.136)" is the
 drop from keeping ONLY Iva (dropping the protein ee), not from removing Iva; and acc_V_AB's "visible" set depends on the panel,
 so ablations are compared here on the visible set of the full panel (fixed denominator).

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python life_few.py > life_few.log

PREDICTIONS (written 2026-10-16, before any run of this script).  NOT blind for the 2-3 row and the ablation: the reviewer's
scratch run (same code path, dataset seed 31) gave, for breadth 2-3, acc 0.622, acc_V_AB 0.753, confident life 0.258 (n = 8 with
Iva), Davila precision 0.992, recall 0.217; and on a fixed visible set no-Iva 0.524, Iva-only 0.524.  The 4-5 and 6 rows are new.
 (R1) Regression, breadth 7: Davila precision 0.523 +- 0.01 and confident life 0.000 (life_test.log values).
 (R2) Breadth 2-3: Davila precision >= 0.98, recall 0.20-0.24; confident life (n = 8 with Iva) 0.22-0.30; acc_V_AB 0.72-0.78.
 (R3) Davila precision falls monotonically with breadth (2-3 > 4-5 > 6 > 7), and confident life falls monotonically too.
 (R4) Breadth 4-5: Davila precision between 0.55 and 0.95 (the rule starts to fail once the amplifier can reach 4 amino acids).
 (R5) Fixed visible set (breadth 7, n = 8 with Iva, sigma 2%): acc_V_AB full minus no-Iva 0.05-0.12, and no-Iva and Iva-only
      within 0.03 of each other (so neither the protein ee nor Iva alone is "decisive").
 (R6) Whole script under 10 minutes at 2 threads under nice.
"""
import time

import numpy as np

import life_test as LT

T0 = time.time()
RES = {}
SIG = 0.02
ROWS = [("2-3", (2, 3)), ("4-5", (4, 5)), ("6", (6, 6)), ("7", (7, 7))]


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def narrow(lo, hi):
    """life_test.gen with the amplifier's protein ee kept in only b in [lo, hi] random amino acids."""
    orig = LT.gen

    def gen(rng, m, kind, rates, E, Amask, **kw):
        x = orig(rng, m, kind, rates, E, Amask, **kw)
        if kind == "A" and lo < 7:
            keep = rng.integers(lo, hi + 1, m)
            x[:, :7] *= np.argsort(rng.random((m, 7)), 1) < keep[:, None]
        return x
    return orig, gen


def run_row(lab, lo, hi, EC, AC):
    orig, gen = narrow(lo, hi)
    LT.gen = gen
    try:
        D = LT.dataset(31 if lo < 7 else 11, LT.RATES_C, EC, AC)   # 11 = life_test.py's seed, so breadth 7 is its own dataset
    finally:
        LT.gen = orig
    r8 = LT.fit(D, SIG, LT.cols_iva(8))
    r7 = LT.fit(D, SIG, LT.cols_prot(7))
    cB = lambda r: ((r["P"][:, 1] >= 0.9) & (r["y"] == 1))[r["y"] == 1].mean()
    called = (r7["obs"][:, :7] >= 0.05).sum(1) >= 4
    prec = (r7["y"][called] == 1).mean() if called.any() else float("nan")
    rec = called[r7["y"] == 1].mean()
    print(f"{lab:>7} {r8['acc']:7.3f} {r8['accV']:9.3f} {cB(r8):10.3f} {r7['acc']:7.3f} {r7['accV']:9.3f} {cB(r7):10.3f} {prec:8.3f} {rec:7.3f}")
    return dict(acc=r8["acc"], accV=r8["accV"], cB=cB(r8), prec=prec, rec=rec, D=D, r8=r8)


def ablation(D, r8):
    """acc_V_AB on the full panel's visible set for every ablation (fixed denominator)."""
    sel = (r8["y"] < 2) & r8["vis"]
    out = {"full": (r8["pred"] == r8["y"])[sel].mean()}
    for nm, cols in (("no Iva (n = 7 protein)", LT.cols_prot(7)), ("Iva only", [LT.IVA])):
        rr = LT.fit(D, SIG, cols)
        out[nm] = (rr["pred"] == rr["y"])[sel].mean()
    return out


def main():
    raw = LT.load(0)
    EC, AC = LT.pool(raw, LT.RATES_C, LT.PS.eta_table())
    print(f"life_few.py: amplifier breadth sweep, sigma {SIG:.0%}, M = {LT.M} per class and split, equal priors")
    print(f"{'breadth':>7} {'acc8':>7} {'accV_AB8':>9} {'confLife8':>10} {'acc7':>7} {'accV_AB7':>9} {'confLife7':>10} "
          f"{'Dav prec':>8} {'Dav rec':>7}")
    out = {lab: run_row(lab, lo, hi, EC, AC) for lab, (lo, hi) in ROWS}
    o7, o2 = out["7"], out["2-3"]
    check("R1", abs(o7["prec"] - 0.523) <= 0.01 and o7["cB"] == 0.0, f"precision {o7['prec']:.3f}, confident life {o7['cB']:.3f}")
    check("R2", o2["prec"] >= 0.98 and 0.20 <= o2["rec"] <= 0.24 and 0.22 <= o2["cB"] <= 0.30 and 0.72 <= o2["accV"] <= 0.78,
          f"precision {o2['prec']:.3f}, recall {o2['rec']:.3f}, confident life {o2['cB']:.3f}, acc_V_AB {o2['accV']:.3f}")
    seq = [out[l] for l, _ in ROWS]
    mono = all(a["prec"] > b["prec"] for a, b in zip(seq, seq[1:])) and all(a["cB"] >= b["cB"] for a, b in zip(seq, seq[1:]))
    check("R3", mono, "precision " + " > ".join(f"{s['prec']:.3f}" for s in seq) + "; confident life " +
          " >= ".join(f"{s['cB']:.3f}" for s in seq))
    check("R4", 0.55 <= out["4-5"]["prec"] <= 0.95, f"breadth 4-5 precision {out['4-5']['prec']:.3f}")
    print("\nAblation on the full panel's visible set (breadth 7, n = 8 with Iva, sigma 2%), acc_V_AB:")
    ab = ablation(o7["D"], o7["r8"])
    for nm, v in ab.items():
        print(f"  {nm:24s} {v:.3f}")
    d = ab["full"] - ab["no Iva (n = 7 protein)"]
    check("R5", 0.05 <= d <= 0.12 and abs(ab["no Iva (n = 7 protein)"] - ab["Iva only"]) <= 0.03,
          f"full - no Iva {d:.3f}; no Iva {ab['no Iva (n = 7 protein)']:.3f} vs Iva only {ab['Iva only']:.3f}")
    t = time.time() - T0
    check("R6", t < 600, f"{t:.0f} s")
    print(f"\n{sum(RES.values())} of {len(RES)} checked items HELD: " + ", ".join(f"{k} {'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
