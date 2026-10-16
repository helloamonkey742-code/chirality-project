"""Correction to sugar_measured.py after a review on 2026-10-16. Author's script; imports sugar_measured.py (share_mc, share_grid)
and edits nothing. Two things the review found:
 (1) Wrong sign. sugar_measured.py counted Wagner et al. 2017 (ACS Cent Sci 3:322-328, doi:10.1021/acscentsci.7b00085) as one
     natural sugar -> amino acid result. Its open text (PMC5408343) says D-ribose and D-lyxose favour opposite enantiomers, that
     D-deoxyribose gives the same sense as D-ribose, and that with D-ribose the major aminoamide is R (the unnatural D side); D-lyxose
     and mixtures of D-sugars give the natural L. So W17 is one natural (D-lyxose / D-sugar mixtures) and one against (D-ribose /
     D-deoxyribose), and the fair tally for sugar -> aa is 3 natural / 1 against, not 3 / 0. The 'against' result is the RNA sugar.
 (2) Band choice. sugar_measured's G_sa band (0.11-0.83) takes its whole range from W17's |ee| values, which mix both signs; its
     top (0.83) is probably D-lyxose. The review's sensitivity rows (its own scratch Monte Carlo) moved the share 0.74-0.85.
 This script recomputes those rows inside the repo and adds the mixed-sign case (k_sa < 0 < k_as) by simulation. Mixed signs are
 what D-ribose -> D-aa (W17) together with L-aa -> D-sugar would mean if W17 were the whole sugar -> aa story.
 Band rows (G_as = amino acid -> sugar gain, G_sa = sugar -> amino acid gain; rho = G_as / G_sa; all else as sugar_measured.py):
   OLD PRIMARY  G_as 0.044-0.44, G_sa 0.11-0.83 (sugar_measured.py, for the regression)
   RNA-ONLY     G_as 0.044-0.44, G_sa 0.57-0.75 (T&S04/06 57-60%, Y10 72-75%: the RNA-templating results, sign natural)
   ALA/SER IN   G_as 0.006-0.44, G_sa 0.11-0.83 (B&C's weakest D-sided catalysts, Ala 1.6% and Ser 0.6% ee, kept in)
   BOTH         G_as 0.006-0.44, G_sa 0.57-0.75
   NO LYXOSE    G_as 0.044-0.44, G_sa 0.11-0.60 (W17's top dropped)
 Mixed sign: dna_network.simulate_asymmetric takes sqrt(k_sa k_as), so the Euler loop is copied here with |k_sa k_as| in the
 start time only; k_sa = -0.3 (sugar -> aa unnatural), k_as = +0.3, g_a = 1e-3, g_s = 0.3e-3, 4000 runs.

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python sugar_bands.py > sugar_bands.log

PREDICTIONS (written 2026-10-16, before any run of this script). NOT blind: the reviewer's scratch run (4e6 draws, another seed)
already gave 0.785 (Ala/Ser in), 0.797 (RNA-only), 0.743 (both), 0.848 (no lyxose), and about 25% for each of the four hand
combinations with mixed signs. The ranges below are that run plus Monte Carlo noise; they check that the repo reproduces it.
 (Q1) Corrected sign tally: aa -> sugar 3 natural / 3 against (unchanged); sugar -> aa 3 natural / 1 against.
 (Q2) OLD PRIMARY reproduces sugar_measured's 0.838 within 0.003 (grid and MC).
 (Q3) Shares by grid: RNA-ONLY 0.790-0.805, ALA/SER IN 0.778-0.792, BOTH 0.736-0.750, NO LYXOSE 0.841-0.855; MC within 0.004 of
      grid in every row. Over these five rows the share stays between 0.73 and 0.86, so 'about 74-85%' over band choices.
 (Q4) Rescue (share of r < 0 draws still favouring the natural pair) is lowest for BOTH, below 0.60, and highest for NO LYXOSE.
 (Q5) Mixed sign: same-sign (matched) fraction 0.40-0.60 and each of the four hand combinations 0.18-0.32, i.e. no pairing forms;
      with both couplings positive (control, same loop) matched >= 0.99.
 (Q6) Whole script under 5 minutes at 2 threads under nice.
RESULTS (sugar_bands.log, 2026-10-16): 6 of 6 held (Q1-Q6), 6 s.
 Tally corrected: sugar -> aa 3 natural / 1 against (W17 D-ribose / D-deoxyribose, the RNA sugars). Natural-pair share by band
 choice: old primary 0.839, RNA-only 0.797, Ala/Ser in 0.785, both 0.743, no lyxose 0.847 (MC and grid agree within 0.0006).
 Rescue of an opposing sugar PVED: 0.49 (both) to 0.69 (no lyxose). Mixed signs: the four hand combinations 0.239-0.257,
 matched 0.491; control with both positive matched 1.000.
 Reading: the 83.8% of sugar_measured.py is one choice in a 74-85% span (78-90% over sugar_measured's other variants). The
 sentence "the experiments fix both signs as natural on the RNA path" was wrong for sugar -> aa: RNA templating is natural, but
 D-ribose in Wagner et al.'s direct synthesis is not, and if that route were the whole sugar -> aa coupling no pairing would form.
"""
import time

import numpy as np

from sugar_measured import share_mc, share_grid   # importing it runs nothing (main is guarded)

T0 = time.time()
RES = {}


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def tally():
    print("=== 1. Sign tally, sign stated in an abstract or open text (+ natural, - against) ===")
    aa = [("B&C non-Pro", +1), ("B&C Pro", -1), ("P&W10 ribose", +1), ("P&W10 lyxose", -1), ("P&W10 arabinose", -1), ("P&W06 erythrose", +1)]
    sa = [("T&S04", +1), ("Y10", +1), ("W17 D-lyxose / D-sugar mix", +1), ("W17 D-ribose / D-deoxyribose", -1)]
    for d, rows in (("aa->sugar", aa), ("sugar->aa", sa)):
        print(f"  {d}: " + ", ".join(f"{n} {'+' if s > 0 else '-'}" for n, s in rows))
    n = lambda rows, s: sum(x == s for _, x in rows)
    check("Q1", (n(aa, 1), n(aa, -1), n(sa, 1), n(sa, -1)) == (3, 3, 3, 1),
          f"aa->sugar {n(aa, 1)}/{n(aa, -1)}, sugar->aa {n(sa, 1)}/{n(sa, -1)}")


ROWS = [("OLD PRIMARY", (0.044, 0.44), (0.11, 0.83)), ("RNA-ONLY", (0.044, 0.44), (0.57, 0.75)),
        ("ALA/SER IN", (0.006, 0.44), (0.11, 0.83)), ("BOTH", (0.006, 0.44), (0.57, 0.75)),
        ("NO LYXOSE", (0.044, 0.44), (0.11, 0.60))]
EXPECT = {"RNA-ONLY": (0.790, 0.805), "ALA/SER IN": (0.778, 0.792), "BOTH": (0.736, 0.750), "NO LYXOSE": (0.841, 0.855)}


def bands():
    rng = np.random.default_rng(20261017)
    print("\n=== 2. Natural-pair share (g_a > 0) under each band choice ===")
    print(f"{'row':12} {'G_as':>12} {'G_sa':>12} {'median rho':>10} {'MC':>7} {'grid':>7} {'rescued':>8}")
    out = {}
    for lab, gA, gB in ROWS:
        mc, _ = share_mc(rng, gA, gB)
        gr = share_grid(gA, gB)
        resc = share_grid(gA, gB, p_plus=0.0)
        med = np.sqrt(gA[0] * gA[1]) / np.sqrt(gB[0] * gB[1])   # median of a ratio of independent log-uniforms
        out[lab] = (mc, gr, resc)
        print(f"{lab:12} {gA[0]:.3f}-{gA[1]:.2f} {gB[0]:.2f}-{gB[1]:.2f} {med:10.3f} {mc:7.4f} {gr:7.4f} {resc:8.4f}")
    mc, gr, _ = out["OLD PRIMARY"]
    check("Q2", abs(mc - 0.838) <= 0.003 and abs(gr - 0.838) <= 0.003, f"MC {mc:.4f}, grid {gr:.4f}")
    bad = {k: round(out[k][1], 4) for k, (lo, hi) in EXPECT.items() if not lo <= out[k][1] <= hi}
    gap = max(abs(m - g) for m, g, _ in out.values())
    span = [g for _, g, _ in out.values()]
    check("Q3", not bad and gap <= 0.004 and min(span) >= 0.73 and max(span) <= 0.86,
          f"outside range: {bad or 'none'}; max |MC - grid| {gap:.4f}; span {min(span):.3f}-{max(span):.3f}")
    rs = {k: v[2] for k, v in out.items()}
    check("Q4", rs["BOTH"] < 0.60 and min(rs, key=rs.get) == "BOTH" and max(rs, key=rs.get) == "NO LYXOSE",
          ", ".join(f"{k} {v:.3f}" for k, v in rs.items()))


def sim(ksa, kas, ga=1e-3, gs=0.3e-3, eps=1e-4, gamma=0.01, runs=4000, dt=0.01, seed=0):
    """dna_network.simulate_asymmetric's loop with |k_sa k_as| in the start time, so mixed signs run."""
    rng = np.random.default_rng(seed)
    a, b = np.zeros(runs), np.zeros(runs)
    for t in np.arange(-(1 + np.sqrt(abs(ksa * kas))) / gamma, 1 / gamma, dt):
        lam = gamma * t
        da = (lam * a - a**3 + ksa * b + ga) * dt
        db = (lam * b - b**3 + kas * a + gs) * dt
        a += da + np.sqrt(eps * dt) * rng.standard_normal(runs)
        b += db + np.sqrt(eps * dt) * rng.standard_normal(runs)
    return {"L-aa,D-sugar": np.mean((a > 0) & (b > 0)), "L-aa,L-sugar": np.mean((a > 0) & (b < 0)),
            "D-aa,L-sugar": np.mean((a < 0) & (b < 0)), "D-aa,D-sugar": np.mean((a < 0) & (b > 0))}


def mixed():
    print("\n=== 3. Mixed coupling signs (k_sa = -0.3: sugar -> aa unnatural, as W17's D-ribose; k_as = +0.3) ===")
    res = {}
    for lab, ksa, kas, seed in (("both +0.3 (control)", 0.3, 0.3, 1), ("k_sa -0.3, k_as +0.3", -0.3, 0.3, 2)):
        f = sim(ksa, kas, seed=seed)
        res[lab] = f
        print(f"  {lab:22} " + "  ".join(f"{k} {v:.3f}" for k, v in f.items()) + f"  matched {f['L-aa,D-sugar'] + f['D-aa,L-sugar']:.3f}")
    c, m = res["both +0.3 (control)"], res["k_sa -0.3, k_as +0.3"]
    mm = m["L-aa,D-sugar"] + m["D-aa,L-sugar"]
    check("Q5", c["L-aa,D-sugar"] + c["D-aa,L-sugar"] >= 0.99 and 0.40 <= mm <= 0.60 and all(0.18 <= v <= 0.32 for v in m.values()),
          f"control matched {c['L-aa,D-sugar'] + c['D-aa,L-sugar']:.3f}; mixed matched {mm:.3f}, combos "
          + "/".join(f"{v:.3f}" for v in m.values()))


def main():
    tally()
    bands()
    mixed()
    t = time.time() - T0
    check("Q6", t < 300, f"{t:.0f} s")
    print(f"\n{sum(RES.values())} of {len(RES)} checked items HELD: " + ", ".join(f"{k} {'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
