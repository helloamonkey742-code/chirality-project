"""Blind test, three more amplifier studies: does a bias act against chance the way the model says outside crystals?

WHY.  Every count of deciding units in Section 3.7 comes from crystallizations (crystal_check.py, blind_units.py), so the
paper cannot say whether a reaction in solution reaches the best case of Eq. 13 (each molecule an independent vote).
Three studies of other amplifiers were picked by title on 2026-10-17 from an OpenAlex search, before any was read.
 PG12: stirring against a chiral dopant in chiral domain formation.  The scout found this is NOT a solution: it is a
       Langmuir monolayer of an achiral azobenzene amphiphile (Azo1) spread on water in a stirred trough, with domains
       scored one by one by Brewster-angle microscopy.  It tests how two biases add, not the solution gap.
 Papers 1-2 (aldol 2009, Soai with achiral amines 2010) are added in a second predictions commit once their setup
 reports are in.  No outcome of PG12 is entered until then, so it cannot steer those predictions.

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
been seen in training; its results are not recalled.

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

SOURCES.
 PG12: Petit-Garrido, N., Claret, J., Ignes-Mullol, J., Sagues, F. 2012, Nature Communications 3, 1001,
       doi:10.1038/ncomms1987 (DOI as printed in the PDF); outcomes from Fig. 2 (eeCW against signed Omega; red circles
       6% w/w S dopant, blue squares 6% R, black triangles racemic 3% R + 3% S; lines are the paper's fits, not data).
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
 Reported, not predicted: K and Omega0 per series; the largest |eeCW| per series.
 Readings fixed in advance.  M1 fails: the lab or the trough favours one hand on its own.  M2 fails: the R and S dopants
 do not cancel in the racemic mix (they interact), so biases do not add.  M3 fails: the dopant changes how sharply the
 hand switches, not only where, so dopant and stirring do not act as one summed bias against chance.  M4 fails: the same
 stirring sense favours different hands with different dopants, which no summed bias can do.  Every outcome is reported
 as found.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python blind_solution.py > blind_solution.log
"""
import numpy as np
from scipy.special import ndtr

from blind_lab import best

RES = {}

# Outcomes: entered after the predictions commit, as read by a separate reader, unchanged.
PG12 = None   # {"R": [(Omega rpm, eeCW), ...], "S": [...], "rac": [...]}


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


def main():
    # Self-check: the fit recovers a known curve and its mirror.
    W = np.linspace(-1500, 1500, 13)
    for K, W0 in ((1.5, 300.0), (-0.8, -450.0)):
        k, w = fit(list(zip(W, 2 * ndtr(K * (W - W0) / 1000) - 1)))
        assert abs(k - K) < 1e-3 and abs(w - W0) < 1, (k, w)
    if PG12 is None:
        print("Predictions committed; outcomes not yet entered (see the docstring).")
        return
    pg12()
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(
        f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
