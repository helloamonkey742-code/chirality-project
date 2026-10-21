"""Blind test, three more amplifier studies: does a bias act against chance the way the model says outside crystals?

WHY.  Every count of deciding units in Section 3.7 comes from crystallizations (crystal_check.py, blind_units.py), so the
paper cannot say whether a reaction in solution reaches the best case of Eq. 13 (each molecule an independent vote).
Three studies of other amplifiers were picked by title on 2026-10-17 from an OpenAlex search, before any was read.
 PG12: stirring against a chiral dopant in chiral domain formation.  The scout found this is NOT a solution: it is a
       Langmuir monolayer of an achiral azobenzene amphiphile (Azo1) spread on water in a stirred trough, with domains
       scored one by one by Brewster-angle microscopy.  It tests how two biases add, not the solution gap.
 MW10: the uncatalysed aldol of acetone and p-nitrobenzaldehyde (3.3% DMSO, 300 rpm, 10 days), claimed to break symmetry
       by product autocatalysis; 16 runs at three volumes (0.15, 0.3, 2.5 mL, one concentration) and 7 runs with extra
       solid aldehyde, none with a chiral additive.  The volume series is the first handle in Section 3.7 on whether a
       solution reaction's chance excess scales with its molecules.
 SH10: the Soai reaction started with an achiral amine (three amines; 75, 20 and 20 runs; no chiral additive), each run
       scored S or R by chiral HPLC.  A chance-zone test like blind_check.py, on runs Section 3.7 has not used.
 PG12's predictions were committed alone (796c2f7).  MW10 and SH10 were added in a second predictions commit, before any
 outcome of the three papers was read; PG12's predictions were not changed.

PROTOCOL.  As blind_units.py.  The predictor (Claude, an AI model) has not read these papers; it saw only setup reports
from separate scout agents told to withhold every outcome.  These predictions are committed before anyone reads an
outcome.  A separate reader then enters the outcomes below, in a later commit, without changing anything above them.
Disclosed exposure, written before any setup report was read: the predictor recalled (i) the Tsogoeva group's claims of
spontaneous symmetry breaking in organocatalytic Mannich and aldol reactions, and that they were disputed; (ii) Soai and
Kawasaki runs without a chiral initiator ending S or R by chance; (iii) Ribo-group work in which vortex sense sets the CD
sign of porphyrin J-aggregates, which it wrongly took to be PG12 (PG12 is from the Sagues group).  No run value of any
of the three was recalled.  The predictor may have seen these papers in training.  PG12 only: (iv) the title says
stirring "competes" with chemical induction; (v) the scout reported that the paper's model has a "compensation" stirring
rate (a signed rate where the hand changes, so a zero crossing is likely in the doped series) and an amplitude lambda < 1
"accounting for domain losses" (a hint that |ee| may stay below 1), and that a statement on behaviour at high rates was
withheld; (vi) the authors' earlier work on stirred undoped monolayers (their ref. 26, Phys. Rev. Lett. 2009) may have
been seen in training; its results are not recalled.  MW10 and SH10 (setup reports read 2026-10-17, after 796c2f7):
(vii) the MW10 title states the claim; its scout withheld which runs had reactant precipitation (so some may have),
"significance-threshold tallies" and a yield range; (viii) the SH10 scout withheld per-amine and pooled S/R counts stated
in the text and whether any run was scored without a hand; (ix) the repo's soai_check.py, written in an earlier session,
lists SH10 under "Not used" and says it reports "stochastic" outcomes; a blind scout confirmed on 2026-10-17 that none of
its runs are among Section 3.7's 91 Soai runs.

WHAT THE MODEL CAN AND CANNOT SAY.
 PG12: each circular domain is one unit whose hand is set when it forms, by a bias b against chance; the run's eeCW =
 (nCW - nCCW)/(nCW + nCCW) is the mean vote of many units.  The dopant gives a fixed bias b_d (+ for one hand, - for its
 mirror, 0 for the racemic 3+3% mix, since biases add); the vortex gives b_s(Omega), odd in the signed rate Omega (CW > 0).
 Per unit P(CW) = Phi(c (b_d + b_s(Omega)) sqrt(n)), so eeCW = 2 Phi(K (Omega - Omega0)) - 1 if b_s is linear in Omega.
 The dopant only moves the zero crossing Omega0; it does not change the steepness K, which is set by chance (n) and the
 stirring coupling.  ASSUMPTION: b_s linear in Omega (the paper's model assumes the same).  If b_s is odd but not linear,
 M1 and M2 still hold and M3 may not.  The model cannot give K or Omega0 (no per-unit bias is known), and so cannot say
 whether |ee| reaches 1 inside the tested range; that is reported, not predicted.  Each Fig. 2 marker is taken as one
 (dopant, Omega) value; replicates per marker are not stated.
 MW10: N is at most the product molecules, about 1e20 (0.15 mL) to 2e21 (2.5 mL), so the chance excess 1/sqrt(N) is about
 1e-10, and a 2% ee from chance alone needs a gain of about 1e8 or more, which the model can neither supply nor rule out.
 If runs end with a hand, the sign is a coin flip unless a stray bias steers (A1), and |ee| does not depend on the sign
 (A3).  Volume: if each molecule were an independent vote (Eq. 13's best case, eps = 1/N) and the gain did not saturate,
 |ee| would scale as N^-1/2, so RMS |ee| would fall by sqrt(2.5/0.15) = 4.1 from the smallest to the largest volume
 (slope -0.5 on log-log axes).  Section 3.7's crystallizations suggest instead that a few early copies decide, a number
 not set by volume (slope 0); A2 predicts that.  Slope 0 is also what a saturating gain, or ee values set by the analysis
 (the authors' 2% threshold), give, so A2 holding does not separate these; A2 failing would be the first sign in Section
 3.7 that a solution reaction reaches eps = 1/N.
 SH10: the stage-1 product, at most 0.2 mmol (1.2e20 molecules), gives a chance excess of about 1e-10; the hand is a coin
 flip only if every stray bias is below that or far fewer units decide.  As blind_check.py, the model predicts the chance
 zone: a 50/50 split (S1), no difference between achiral amines, which cannot carry a hand (S2), no run depending on the
 one before (S3), and |ee| not depending on the sign (S4).  ASSUMPTION for S3: table order is run order; if not, S3 tests
 only the listed order.

SOURCES.
 PG12: Petit-Garrido, N., Claret, J., Ignes-Mullol, J., Sagues, F. 2012, Nature Communications 3, 1001,
       doi:10.1038/ncomms1987 (DOI as printed in the PDF); outcomes from Fig. 2 (eeCW against signed Omega; red circles
       6% w/w S dopant, blue squares 6% R, black triangles racemic 3% R + 3% S; lines are the paper's fits, not data).
 MW10: Mauksch, M., Wei, S., Freund, M., Zamfir, A., Tsogoeva, S.B. 2010, Origins of Life and Evolution of Biospheres 40,
       79-91 (online 2009), doi:10.1007/s11084-009-9177-2; outcomes from Table 1 (16 runs: volume, ee %, configuration),
       Table 2 (7 runs) and Fig. 2 (two seeded reactions, reported only: the seed's hand is unknown to the predictor).
 SH10: Suzuki, K., Hatase, K., Nishiyama, D., Kawasaki, T., Soai, K. 2010, Journal of Systems Chemistry 1, 5,
       doi:10.1186/1759-2208-1-5; outcomes from Tables 1-3 (amines 3, 4, 5: Run, ee, Config.), in table order, not from
       Fig. 2 or the text's counts.  DOIs as printed in the PDFs.
 Values in figures are read by the separate reader from page images, one entry per marker, (Omega rpm, eeCW).

PREDICTIONS (written 2026-10-17, before any outcome was read).  Each series is fitted alone by least squares to
eeCW = 2 Phi(K (Omega - Omega0)) - 1 (K signed, |K| in [0.01, 100] per 1000 rpm).  A series "crosses" if it has markers
with eeCW > 0 and with eeCW < 0.  D = |Omega0(R) - Omega0(S)| / 2, the shift one dopant gives.
 (M1) PG12 mirror: |Omega0(R) + Omega0(S)| <= 0.25 |Omega0(R) - Omega0(S)| (the R and S curves are mirror images; no
      stray lab bias).  NOT SCORABLE unless R and S both cross and D >= 100 rpm.
 (M2) PG12 biases add: |Omega0(rac) - (Omega0(R) + Omega0(S)) / 2| <= 0.25 D (the racemic dopant acts as no dopant).
      NOT SCORABLE unless all three series cross and D >= 100 rpm.
 (M3) PG12 shift, not reshape: max |K| / min |K| over the three series <= 2.  NOT SCORABLE unless every series has >= 4
      markers, >= 2 of them with |eeCW| <= 0.8.
 (M4) PG12 one stirring hand: K has the same sign in all three series.  NOT SCORABLE unless every series has >= 4
      markers.
 MW10 and SH10.  A run "has a hand" if its configuration is given; runs without one are counted, not scored.
 (A1) MW10 coin flip: among Table 1-2 runs with |ee| >= 2% (the authors' reliability threshold) and a hand, the split fits
      50/50 (two-sided binomial p > 0.05).  NOT SCORABLE if fewer than 6 such runs.
 (A2) MW10 volume does not set the excess: the least-squares slope of log RMS ee (Table 1, all runs with a value) on log
      volume over 0.15, 0.3 and 2.5 mL is > -0.25.  NOT SCORABLE unless every volume has >= 3 runs with a value and >= 3
      Table 1 runs have |ee| >= 2%.
 (A3) MW10 mirror: |ee| of S and R runs (Tables 1-2) come from one distribution (Mann-Whitney two-sided p > 0.05).  NOT
      SCORABLE unless each hand has >= 4 runs.
 (S1) SH10 chance zone, pooled: the split fits 50/50 (two-sided p > 0.05) and rules out a lab-wide steering bias of
      Delta >= 1 (chance of the split under it < 1e-3), as blind_lab.py H4.
 (S2) SH10 achiral amines carry no hand: the dispersion D (as blind_check.py B3) over the three amines has
      P(D >= observed | independent coin flips) > 0.05 (200,000 sets, seed 17).
 (S3) SH10 no carry-over: the number of hand changes between consecutive runs with a hand, summed over Tables 1-3, is what
      shuffling each table gives (two-sided P > 0.05, 200,000 shuffles, seed 17).
 (S4) SH10 mirror: |ee| of S and R runs (pooled) come from one distribution (Mann-Whitney two-sided p > 0.05).  NOT
      SCORABLE unless each hand has >= 5 runs.
 Reported, not predicted: K and Omega0 per series; the largest |eeCW| per series.  MW10: RMS |ee| per volume, runs with
 |ee| >= 2%, Table 2 and Fig. 2 values.  SH10: the split per amine, runs with |ee| >= 90%.
 Readings fixed in advance.  M1 fails: the lab or the trough favours one hand on its own.  M2 fails: the R and S dopants
 do not cancel in the racemic mix (they interact), so biases do not add.  M3 fails: the dopant changes how sharply the
 hand switches, not only where, so dopant and stirring do not act as one summed bias against chance.  M4 fails: the same
 stirring sense favours different hands with different dopants, which no summed bias can do.  A1 or S1 fails: a stray bias
 steered the hand across runs (reagents, vessels or analysis); the counting rule allows it (any bias above about 1e-10
 wins here), but the runs are then not spontaneous.  A2 fails (slope <= -0.25): the excess scales with the molecules in
 the vessel, the best case Eq. 13 assumes.  A3 or S4 fails: one hand ends with larger |ee|, a bias acting during
 amplification or in the analysis, not chance at the start.  S2 fails: runs with different amines lean differently; an
 achiral amine cannot, so an amine lot or that set of runs carried a stray bias.  S3 fails low: a bias persisted over
 consecutive runs (a stock, a lot, or product carried over); fails high: runs alternate, which no stray bias gives.
 Every outcome is reported as found.

RESULTS (blind_solution.log, 2026-10-17): 8 of 11 held (M1-M4, A1, A3, S1, S2); 3 FAILED (A2, S3, S4).
Values were read by separate readers after the predictions commits (796c2f7, 2873906) and entered unchanged.
 PG12: Fig. 2 is vector graphics, so the reader took all 39 markers (13 per series) from the drawing commands.  Every
 series crosses.  Omega0 = +222 (R), -211 (S), -26 (rac) rpm; K = +1.55, +1.53, +1.15 per 1000 rpm.  M1 held (|R + S| =
 11 against 108 rpm), M2 held (31 against 54), M3 held (|K| ratio 1.35), M4 held.  The paper's own fit, lambda (Omega -
 Omega_c)/(|Omega| + a), puts the crossings at +120, -110, -11 rpm: a different curve, the same mirror and order.  No
 series passes |eeCW| 0.78.
 MW10: A1 held (8 S / 5 R among 13 runs with |ee| >= 2%, p = 0.58); A3 held (p = 0.37).  A2 FAILED: RMS |ee| 4.80%
 (0.15 mL), 7.54% (0.3 mL), 2.32% (2.5 mL), slope -0.33.  The fall is not monotonic; 0.3 mL is highest because of one
 14.3% run.  Post hoc (not scored): in Table 1 the hand follows the volume (all six 2.5 mL runs R, 9 of 10 smaller runs
 S, Fisher p = 0.0009), and reactant precipitated only in the smaller runs (1-4, 11-16, as the paper says), so volume is
 confounded with hand and with precipitation.  Table 2's largest |ee| is 17.2%; Fig. 2 is reported in the log.
 SH10: Table 1 prints 76 rows for 75 runs (run 26 twice); both rows kept as printed.  S1 held (59 S of 116, p = 0.93;
 a lab-wide steer of Delta >= 1 has chance 7e-17).  S2 held (P = 0.93).  S3 FAILED high: 74 hand changes against about
 58 expected (P = 0.003), mostly amine 3 (48 of 75, P = 0.027 alone).  S4 FAILED: R runs end with larger |ee| (median
 19% against 11%, p = 0.016), from amine 3 (19 against 7); amine 5 goes the other way (S 27, R 15).  Dropping either
 run 26 changes no verdict (S3 P = 0.007 / 0.006, S4 p = 0.017 / 0.009).  No run reaches 90% (largest 57%).
 Reading: PG12 is what Section 3.7 assumes of a summed bias: the dopant shifts where the hand switches and not how
 sharply, the mirror dopants cancel, and the racemic mix acts as none (in a monolayer, not a solution).  A2 failed in
 the direction fixed in advance for "the excess scales with the molecules", but the confound above (one hand per volume,
 precipitation only in small volumes) means it cannot be read as that; the solution gap stays open.  SH10's split is the
 chance zone, as for the 91 runs already in Section 3.7, but its runs are not independent coin flips in table order:
 they alternate too often, and in amine 3 the R runs end with larger |ee|.  The paper does not say the run numbers are
 chronological and lists several experimenters, so S3 may test the listing, not the time series (the stated assumption).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python blind_solution.py > blind_solution.log
"""
import numpy as np
from scipy.special import ndtr
from scipy.stats import mannwhitneyu

from blind_lab import best, dispersion, steer, two_sided

RES = {}

# Outcomes: entered after the predictions commit, as read by a separate reader, unchanged.
PG12 = None   # {"R": [(Omega rpm, eeCW), ...], "S": [...], "rac": [...]}
MW10 = None   # {"T1": [(volume mL, ee % or None, "R"/"S"/None), ...], "T2": [(ee %, config), ...], "F2": "as printed"}
SH10 = None   # {"3": [(ee %, "R"/"S"/None), ...] in table order, "4": [...], "5": [...]}
# PG12: 13 markers per series read from Fig. 2's vector drawing commands (p. 3; x +-0.5 rpm, eeCW +-0.0005); the inset and
# fit lines excluded.  Three markers hidden under another with identical coordinates were recovered (S -750, rac -1600,
# rac -500).
PG12 = {"R": [(-1600.3, -0.7404), (-1000.2, -0.7700), (-750.0, -0.6802), (-499.8, -0.5803), (-249.6, -0.5703),
              (-99.9, -0.6004), (-0.3, -0.5401), (100.2, -0.1501), (250.0, 0.2299), (500.2, 0.3900), (750.4, 0.6802),
              (999.8, 0.7801), (1599.8, 0.7801)],
        "S": [(-1600.3, -0.7801), (-1000.2, -0.7103), (-750.0, -0.6802), (-499.8, -0.4603), (-249.6, -0.2500),
              (-99.9, 0.1300), (-0.3, 0.5401), (100.2, 0.5502), (250.0, 0.5903), (500.2, 0.6205), (750.4, 0.6601),
              (999.8, 0.7404), (1599.8, 0.7404)],
        "rac": [(-1600.3, -0.7404), (-1000.3, -0.6601), (-750.1, -0.6004), (-499.8, -0.4603), (-249.6, -0.4202),
                (-99.8, -0.2801), (-0.3, 0.0402), (100.2, 0.3303), (250.0, 0.4202), (500.2, 0.5301), (750.4, 0.6300),
                (999.8, 0.6802), (1599.8, 0.7605)]}
# MW10: Table 1 (p. 82) and Table 2 (p. 83) as printed; every ee has a configuration.  Fig. 2: the seed's descriptor is not
# printed; the reader derived it from the drawn wedges.
MW10 = {"T1": [(0.3, 14.3, "S"), (0.3, 2.7, "S"), (0.3, 2.3, "S"), (0.3, 3.2, "S"),
               (2.5, 3.2, "R"), (2.5, 4.1, "R"), (2.5, 0.9, "R"), (2.5, 1.7, "R"), (2.5, 0.3, "R"), (2.5, 1.2, "R"),
               (0.15, 3.6, "S"), (0.15, 7.3, "S"), (0.15, 3.6, "S"), (0.15, 7.4, "S"), (0.15, 1.2, "S"), (0.15, 1.7, "R")],
        "T2": [(0.6, "S"), (1.3, "S"), (17.2, "R"), (1.7, "R"), (2.5, "R"), (1.0, "R"), (7.9, "R")],
        "F2": "Seed 4 (30 mol%, 95.4% ee), drawn (2S,1'R) anti (descriptor not printed; read from the wedges). "
              "(a) acetone + aldehyde -> 3: 9.0% ee, configuration not given, only traces. "
              "(b) cyclohexanone + aldehyde in THF -> 4: 4% yield, 96.0% ee, drawn with the seed's configuration."}
# SH10: Tables 1-3 (pp. 3-4) as printed, column by column (= ascending run number); no run lacks a hand.  Table 1 prints
# 76 rows for 75 runs: run 26 appears twice (end of column 1 "26 12 S", start of column 2 "26 3 R"); both are kept here.
# The text and Fig. 2a give 39 S : 36 R for amine 3; the 76 rows give 39 S : 37 R.  The paper does not say whether run
# numbers are chronological.
SH10 = {"3": [(19, "R"), (18, "S"), (8, "R"), (10, "S"), (11, "R"), (8, "S"), (5, "S"), (9, "R"), (6, "S"), (4, "S"),
              (4, "S"), (57, "R"), (2, "R"), (37, "S"), (23, "R"), (4, "S"), (6, "R"), (4, "S"), (7, "S"), (47, "R"),
              (29, "R"), (25, "R"), (10, "S"), (12, "S"), (31, "R"), (12, "S"), (3, "R"), (7, "S"), (4, "S"), (14, "R"),
              (6, "S"), (6, "R"), (4, "S"), (4, "R"), (3, "S"), (15, "S"), (20, "S"), (18, "R"), (9, "R"), (37, "S"),
              (27, "S"), (28, "R"), (3, "R"), (9, "S"), (32, "S"), (29, "R"), (11, "S"), (16, "R"), (52, "R"), (9, "R"),
              (4, "S"), (19, "R"), (36, "S"), (23, "S"), (32, "R"), (32, "S"), (4, "S"), (6, "R"), (12, "S"), (15, "S"),
              (5, "S"), (42, "R"), (52, "R"), (2, "S"), (43, "R"), (35, "R"), (3, "S"), (2, "R"), (4, "S"), (4, "S"),
              (42, "R"), (32, "R"), (18, "R"), (6, "S"), (9, "R"), (19, "R")],
        "4": [(50, "S"), (57, "R"), (47, "R"), (55, "S"), (14, "R"), (27, "R"), (10, "R"), (53, "S"), (39, "R"), (40, "R"),
              (37, "S"), (20, "S"), (47, "R"), (47, "R"), (36, "S"), (24, "S"), (22, "R"), (5, "S"), (35, "R"), (11, "S")],
        "5": [(34, "S"), (37, "R"), (52, "S"), (33, "S"), (20, "R"), (16, "R"), (33, "S"), (15, "R"), (33, "S"), (7, "R"),
              (10, "S"), (6, "S"), (7, "R"), (7, "S"), (6, "S"), (9, "R"), (29, "R"), (12, "S"), (12, "R"), (27, "S")]}


def check(p, ok, detail):
    ok = ok if isinstance(ok, str) else bool(ok)
    RES[p] = ok
    print(f"  {p} {'HELD' if ok is True else 'FAILED' if ok is False else ok}: {detail}")


def fit(pts):
    """Least squares eeCW = 2 Phi(K (W - W0)) - 1, W in 1000 rpm.  Returns (K per 1000 rpm, W0 in rpm)."""
    W, ee = np.asarray(pts, float).T
    W = W / 1000
    k = lambda x: np.sign(x[0]) * np.clip(abs(x[0]), 0.01, 100)
    f = lambda x: float(np.sum((2 * ndtr(k(x) * (W - x[1])) - 1 - ee) ** 2))
    x, _ = best(f, [(s * K, w) for s in (1, -1) for K in (0.5, 2.0, 8.0) for w in (-0.5, 0.0, 0.5)])
    return float(k(x)), float(x[1] * 1000)


def crosses(pts):
    ee = [e for _, e in pts]
    return max(ee) > 0 > min(ee)


def pg12():
    print("\n=== PG12, Fig. 2 ===")
    F = {s: fit(p) for s, p in PG12.items()}
    for s, p in PG12.items():
        print(f"  {s}: {len(p)} markers, K = {F[s][0]:+.3f} per 1000 rpm, Omega0 = {F[s][1]:+.0f} rpm, "
              f"max |eeCW| {max(abs(e) for _, e in p):.2f}, crosses {crosses(p)}")
    (_, wR), (_, wS), (_, wr) = F["R"], F["S"], F["rac"]
    D = abs(wR - wS) / 2
    if crosses(PG12["R"]) and crosses(PG12["S"]) and D >= 100:
        check("M1", abs(wR + wS) <= 0.25 * abs(wR - wS), f"|Omega0 R + S| = {abs(wR + wS):.0f} rpm vs {0.25 * abs(wR - wS):.0f}")
    else:
        check("M1", "NOT SCORABLE", f"crossing or D = {D:.0f} rpm < 100")
    if all(crosses(p) for p in PG12.values()) and D >= 100:
        check("M2", abs(wr - (wR + wS) / 2) <= 0.25 * D, f"|Omega0 rac - mid| = {abs(wr - (wR + wS) / 2):.0f} rpm vs {0.25 * D:.0f}")
    else:
        check("M2", "NOT SCORABLE", f"crossing or D = {D:.0f} rpm < 100")
    if all(len(p) >= 4 and sum(abs(e) <= 0.8 for _, e in p) >= 2 for p in PG12.values()):
        K = [abs(F[s][0]) for s in PG12]
        check("M3", max(K) / min(K) <= 2, f"|K| ratio {max(K) / min(K):.2f}")
    else:
        check("M3", "NOT SCORABLE", "too few markers on the transition")
    if all(len(p) >= 4 for p in PG12.values()):
        check("M4", len({np.sign(F[s][0]) for s in PG12}) == 1, "signs " + ", ".join(f"{s} {F[s][0]:+.2f}" for s in PG12))
    else:
        check("M4", "NOT SCORABLE", "fewer than 4 markers")


def mirror(rows, need):
    """Mann-Whitney two-sided p for |ee| of S runs vs R runs; None if either hand has fewer than `need` runs."""
    a = [abs(e) for e, c in rows if c == "S" and e is not None]
    b = [abs(e) for e, c in rows if c == "R" and e is not None]
    return mannwhitneyu(a, b, alternative="two-sided").pvalue if min(len(a), len(b)) >= need else None


def order_test(seqs, draws=200_000, seed=17):
    """Summed hand changes between consecutive runs, and its two-sided P under shuffles within each sequence."""
    rng = np.random.default_rng(seed)
    obs, sim = 0, np.zeros(draws)
    for q in seqs:
        a = np.array([c == "S" for c in q])
        obs += int((a[1:] != a[:-1]).sum())
        m = rng.permuted(np.tile(a, (draws, 1)), axis=1)
        sim += (m[:, 1:] != m[:, :-1]).sum(1)
    return obs, min(1.0, 2 * min((sim <= obs).mean(), (sim >= obs).mean()))


def mw10():
    print("\n=== MW10, Tables 1-2 ===")
    rows = [(e, c) for _, e, c in MW10["T1"]] + list(MW10["T2"])
    hand = [c for e, c in rows if c in ("R", "S") and e is not None and abs(e) >= 2]
    k, n = hand.count("S"), len(hand)
    print(f"  {len(rows)} runs; {n} with |ee| >= 2% and a hand: {k} S, {n - k} R")
    if n >= 6:
        check("A1", two_sided(k, n) > 0.05, f"p = {two_sided(k, n):.3f} against 50/50")
    else:
        check("A1", "NOT SCORABLE", f"{n} runs with |ee| >= 2% and a hand")
    vols = sorted({v for v, _, _ in MW10["T1"]})
    g = {v: [e for w, e, _ in MW10["T1"] if w == v and e is not None] for v in vols}
    rms = {v: float(np.sqrt(np.mean(np.square(g[v])))) if g[v] else 0.0 for v in vols}
    for v in vols:
        print(f"  {v} mL: {len(g[v])} runs with a value, RMS |ee| {rms[v]:.2f}%")
    big = sum(abs(e) >= 2 for v in vols for e in g[v])
    if len(vols) == 3 and all(len(g[v]) >= 3 for v in vols) and big >= 3 and min(rms.values()) > 0:
        slope = np.polyfit(np.log(vols), np.log([rms[v] for v in vols]), 1)[0]
        check("A2", slope > -0.25, f"slope of log RMS ee on log volume {slope:+.2f} (-0.5 if each molecule votes)")
    else:
        check("A2", "NOT SCORABLE", f"{big} Table 1 runs with |ee| >= 2%, runs per volume {[len(g[v]) for v in vols]}")
    pm = mirror(rows, 4)
    check("A3", pm > 0.05 if pm is not None else "NOT SCORABLE",
          f"Mann-Whitney p = {pm:.3f}" if pm is not None else "fewer than 4 runs of a hand")
    print(f"  Table 2 (reported): {MW10['T2']}\n  Fig. 2 (reported): {MW10['F2']}")


def sh10():
    print("\n=== SH10, Tables 1-3 ===")
    seqs = {a: [c for _, c in rs if c in ("R", "S")] for a, rs in SH10.items()}
    split = {a: (q.count("S"), q.count("R")) for a, q in seqs.items()}
    for a, rs in SH10.items():
        print(f"  amine {a}: {len(rs)} runs, {split[a][0]} S, {split[a][1]} R, {len(rs) - len(seqs[a])} without a hand")
    k, n = sum(x for x, _ in split.values()), sum(map(len, seqs.values()))
    check("S1", two_sided(k, n) > 0.05 and steer(k, n) < 1e-3,
          f"{k} S of {n}, p = {two_sided(k, n):.3f} against 50/50, P(split | Delta = 1) = {steer(k, n):.2g}")
    d, pd = dispersion([split[a] for a in split if sum(split[a]) >= 2])
    check("S2", pd > 0.05, f"D = {d:.2f} over the amines, P = {pd:.3f}")
    obs, po = order_test(list(seqs.values()))
    check("S3", po > 0.05, f"{obs} hand changes between consecutive runs, two-sided P = {po:.3f}")
    rows = [r for rs in SH10.values() for r in rs]
    pm = mirror(rows, 5)
    check("S4", pm > 0.05 if pm is not None else "NOT SCORABLE",
          f"Mann-Whitney p = {pm:.3f}" if pm is not None else "fewer than 5 runs of a hand")
    print(f"  |ee| >= 90%: {sum(e is not None and abs(e) >= 90 for e, _ in rows)} of {len(rows)} runs (reported)")


def main():
    # Self-check: the fit recovers a known curve and its mirror.
    W = np.linspace(-1500, 1500, 13)
    for K, W0 in ((1.5, 300.0), (-0.8, -450.0)):
        k, w = fit(list(zip(W, 2 * ndtr(K * (W - W0) / 1000) - 1)))
        assert abs(k - K) < 1e-3 and abs(w - W0) < 1, (k, w)
    # Self-check: the order test flags blocks (too few changes) and alternation (too many), not a fair shuffle of them.
    assert order_test([list("S" * 10 + "R" * 10)], 20_000)[1] < 0.01
    assert order_test([list("SR" * 10)], 20_000)[1] < 0.01
    if all(x is None for x in (PG12, MW10, SH10)):
        print("Predictions committed; outcomes not yet entered (see the docstring).")
        return
    for data, score in ((PG12, pg12), (MW10, mw10), (SH10, sh10)):
        if data is not None:
            score()
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(
        f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
