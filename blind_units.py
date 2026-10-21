"""Blind units test: how many independent units decide a hand, and does volume or stirring change that number as the model says?

WHY.  blind_lab.py's Viedma half failed low (C1): Choi et al. 2024 behaved as if about 34 units decided, not the 100-400
starting crystals.  The model's counting rule says a run's hand is the vote of n_eff independent units, each choosing a
hand by chance, so <ee^2> = 1/n_eff (crystal_check.py).  Three consequences can be checked without any added bias:
n_eff grows with the volume one unit must capture; transport (stirring) carries one unit's copies through the vessel and
lowers n_eff; and n_eff cannot exceed the number of single-handed units in the vessel (unequal sizes only lower it), so
"fewer deciding units than visible ones", the direction C1 failed in, is what the rule expects.  Two papers measure this:
 Bak et al. 2023: NaClO3 crystallizing in evaporating sessile drops under oil (unstirred, no seed); each drop is scored
 all-L, all-D or mixed; P_CSB = share of drops with one hand, plotted against the drop volume at nucleation, Vc.
 Asakura et al. 2002: 1,1'-binaphthyl crystallizing from its racemizing melt (0.50 g, 2.0 cm tube, no seed); ee of the
 whole solid per run, unstirred at 150 and 152 C (10 runs each), stirred at 200 and 500 rpm (150 C), and ee of single
 masses grown from one nucleus (Table 1).

PROTOCOL.  As blind_lab.py.  The predictor (Claude, an AI model) has not read either paper; it saw only two setup reports
from separate scout agents told to withhold every outcome.  These predictions are committed before anyone reads an
outcome.  A separate reader then enters the outcomes below, in a later commit, without changing anything above them.
The author obtained both PDFs on 2026-10-17.  Disclosed exposure: (i) both titles were seen ("Solution volume effects on
spontaneous chiral symmetry breaking of sodium chlorate crystals"; "Probability distributions of enantiomeric excess in
unstirred and stirred crystallization of 1,1'-binaphthyl melt"); (ii) the Bak scout said the paper's own model computes a
"critical volume" (values withheld), which hints at a threshold; the Asakura scout said the paper proposes stirring rate
as the control parameter of a Langevin model whose stationary distribution has one peak below a critical value and two
above it; (iii) frank_citers.csv, read by the predictor in this session, notes "bimodal ee distribution" for two later
Asakura papers on the binaphthyl growth front (2004, 2006), not these runs; (iv) the predictor may have seen these papers
in training; it recalls the general finding of Kondepudi's NaClO3 work that stirring gives large |ee| (already in
Section 3.7 of the paper, n_eff 1.0-1.1 stirred, about 17 unstirred) and does not recall these papers' values.

WHAT THE MODEL CAN AND CANNOT SAY.
 Bak: a drop's first nucleus starts copying its hand; if K more primary nuclei form before the copies fill the drop, each
 a coin flip, the drop ends with one hand with chance E[2^-K].  With primary nuclei arriving at a fixed rate per volume,
 K is Poisson with mean proportional to V, and P_CSB = exp(-V/V0), V0 unknown.  ASSUMPTION: within one series (same Tsat,
 Tinc, oil and substrate) every drop has the same nucleation rate per volume and the same capture time.  The scout reports
 that the supersaturation at onset varies from drop to drop and is linked to Vc; if small drops nucleate at higher
 supersaturation, P_CSB falls less steeply (or rises) at small V, and this prediction ignores it.  The model cannot give
 V0.  It gives the direction (P_CSB falls with V) and the breadth of the fall: 0.8 -> 0.2 over a factor 7.2 in V (0.9 ->
 0.1 over a factor 21.9), not a sharp critical volume.  In the free-width form P = 1/(1 + (V/Vk)^h) that breadth is
 h = 1.40.  P_CSB counts a one-crystal drop as homochiral, as the rule does (one unit).
 Asakura: a run's ee is the mass-weighted vote of its independent single-handed masses, so <ee^2> = sum w_i^2 >= 1/n and
 n_eff = 1/<ee^2> <= n, the number of masses (equal only if all are the same size).  Stirring spreads one mass's secondary
 nuclei through the melt, so fewer units decide.  ASSUMPTION (transfer, not derived): 500 rpm makes the tube about one
 unit, as stirring did for NaClO3 in solution.  More supercooling (150 C vs 152 C) raises primary nucleation faster than
 growth, so more masses form: n_eff(150) >= n_eff(152).  With no bias, ee > 0 and ee < 0 are equally likely.  Sampling:
 with m runs of a Gaussian vote, n_eff_hat * chi2_{0.025,m}/m is the lower end of a 95% range for n_eff (x 0.325 for
 m = 10); for ee near +-1 the Gaussian form is rough.

SOURCES (DOIs checked against each PDF's own text, 2026-10-17).
 BAK23: Bak et al. 2023, J. Phys. Chem. Lett. 14, 785, doi:10.1021/acs.jpclett.2c03873; outcomes from Fig. 2a (P_CSB vs
        Vc; main series = Tsat 50 C, Tinc 25 C, FC-40 oil, PMMA, blue squares), drop counts and all-L / all-D counts from
        the text, captions or Supporting Information if given.
 AS02:  Asakura et al. 2002, Chirality 14, 85, doi:10.1002/chir.10046; outcomes from Fig. 2 (unstirred ee per run, 150 and
        152 C), Fig. 4 (signed ee per run vs stirring rate; ee > 0 = excess S), Table 1 (single masses), and the text's
        count of visible crystal aggregates.
 Values in figures are read by the separate reader from page images (one marker per run; stacked markers counted).

PREDICTIONS (written 2026-10-17, before any outcome of either paper was read).  ee as a fraction; n_eff = 1/mean(ee^2)
over a set's runs.
 (B1) Bak direction: over the main series' markers, the Spearman rank correlation of P_CSB with Vc is < 0 (needs >= 3
      markers; otherwise NOT SCORABLE).
 (B2) Bak breadth: P = 1/(1 + (V/Vk)^h) fitted to the main series (binomial likelihood if drops per marker are given,
      least squares otherwise; h in [0.1, 50]) gives h <= 3, i.e. P_CSB falls from 0.8 to 0.2 over at least a factor 2.5
      in V.  NOT SCORABLE unless the main series has >= 4 markers, one with P_CSB >= 0.6 and one <= 0.4.
 (B3) Bak chance hand: if the paper gives the numbers of all-L and all-D drops, they fit 50/50 (two-sided p > 0.05).
      Otherwise NOT SCORABLE.
 (A1) Asakura chance hand: over every run in Fig. 4 (all rates pooled; ee = 0 counts as no hand), ee > 0 vs ee < 0 fits
      50/50 (two-sided p > 0.05) and rules out a steering bias Delta >= 1 (chance of the split under it < 1e-3).
      NOT SCORABLE if fewer than 20 runs have a hand (below 20, even a perfect split cannot rule out Delta >= 1).
 (A2) Asakura stirring lowers n_eff: n_eff(200 rpm) and n_eff(500 rpm) are both below n_eff(unstirred, 150 C, Fig. 2).
 (A3) Asakura stirring makes about one unit: n_eff(500 rpm) <= 2.
 (A4) Asakura supercooling: n_eff(unstirred, 150 C) >= n_eff(unstirred, 152 C).
 (A5) Asakura one nucleus, one hand: every single-nucleus mass in Table 1 has |ee| >= 50%.  NOT SCORABLE if none.
 (A6) Asakura n_eff no larger than the units seen: where the paper gives the number of visible aggregates (per
      temperature, else pooled over both), the lower end of n_eff's 95% range is <= the mean aggregates per run.
      NOT SCORABLE if no count is given.
 Reported, not predicted: Bak V0 (exponential fit) and the fitted Vk; n_eff for every Asakura set, including Fig. 4's
 0 rpm column.
 Readings fixed in advance.  B1 fails: larger drops end with one hand more often, so the number of deciding units does
 not grow with volume (a supersaturation or capture-time effect dominates).  B2 fails: a sharp critical volume, not
 chance, sets whether a drop ends with one hand.  A2 or A3 fails: stirring does not merge a melt into one unit as it does
 a solution.  A4 fails: more supercooling gave fewer units.  A5 fails: a mass grown from one nucleus is not single-handed,
 so the model's unit is smaller than a mass.  A6 fails: more independent decisions set the ee than there are visible
 aggregates (aggregates are mixed), and Choi's n_eff below the crystal count is not the general case.  B3 or A1 fails: a
 stray bias steered the runs.  Every outcome is reported as found.

RESULTS (blind_units.log, 2026-10-17): 5 of 9 held (B1, A1, A2, A3, A4); 2 FAILED (A5, A6); 2 NOT SCORABLE (B2, B3).
Values were read by separate readers after the predictions commit (c823968) and entered unchanged.
 Bak, main series: only 3 markers, P_CSB 1.0 at 0.164 uL, 8/10 at 0.870 uL, 7/10 at 1.726 uL.  B1 held (rho = -1),
 weakly: a random order of 3 gives rho = -1 one time in 6.  B2 NOT SCORABLE: 3 markers and none at or below 0.4.  The
 exponential fit gives V0 = 4.6 uL (least squares).  B3 NOT SCORABLE: no all-L / all-D counts in the main text, and the
 SI was not reachable (403).  The P = 1 markers of all series share one "(66/66)" label, so their drop counts are unknown.
 Post hoc (printed, not scored): the Tinc 90 C series has 11 markers falling from 1.0 to 5/16; its free-width fit gives
 h = 1.75, a broad fall like the counting rule's 1.40, not a step at the paper's critical volume (about 0.19 uL).
 Asakura, unstirred: n_eff = 272 at 150 C (7 of 10 runs at ee 0, the largest 17%) and 8.6 at 152 C.  Stirred at 150 C:
 n_eff = 3.2 at 200 rpm and 1.96 at 500 rpm.  A1 held: 11 S / 12 R over 23 runs with a hand (P(split | Delta = 1) =
 3.3e-4); without the two 0 rpm markers at |ee| <= 1, 10 S / 11 R and 6.5e-4, still held.  A2 and A3 held.  A3 is close
 to its threshold but does not hang on the readings: the printed means (+70, -72) fix mean |ee| = 0.711, and mean ee^2 >=
 mean |ee|^2 gives n_eff <= 1.98 for any spread.  A4 held, but not by the stated mechanism: the paper counts fewer visible
 aggregates at 150 C (24 in 10 runs) than at 152 C (32).  A5 FAILED: the 13 masses the paper says grew from one nucleus
 include every 150 C mass at ee 0 and 6 of 9 at 152 C at |ee| <= 4%; only 3 (31, 59, 83%) are clearly one-handed.  A6
 FAILED at 150 C: n_eff's 95% low end is 88 against 2.4 aggregates per run (at 152 C, 2.8 against 3.2 holds).  0 rpm
 (reported): 4 separable markers give n_eff = 111; with the uncountable blob as 6 more zeros (10 runs) it is 276, the
 same as Fig. 2's 150 C runs, which those may be.
 Code-only additions after the first run: BAK_TINC90 and the POST HOC line.  No prediction, threshold or value changed.
 Reading: A5 and A6 fail together, in the direction fixed in advance: at 150 C the model's unit is not the visible mass.
 The paper says the same (Discussion): the solid grows as an aggregate of small R and S crystallites, so a mass from one
 nucleus is racemic and the run's ee is the vote of about 90 or more crystallites.  The counting rule then still holds
 (<ee^2> = 1/n_eff, n_eff large, ee near 0); what fails is equating a unit with an aggregate.  Choi's n_eff below the
 crystal count (blind_lab.py C1) is not general: here n_eff is far above the visible count.  Stirring cut n_eff from
 about 270 to about 2 in the same melt at the same temperature, as for NaClO3 in solution.  Bak's main series is too
 short to test the breadth of the fall; the one long series, not pre-registered, falls broadly.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python blind_units.py > blind_units.log
"""
from math import log

import numpy as np
from scipy.stats import chi2, spearmanr

from blind_lab import best, loglik, steer, two_sided

RES = {}

# Outcomes: entered after the predictions commit (c823968), as read by separate readers on 2026-10-17, unchanged.
# Figure values are marker centres from 8x renders against the axis ticks (about +-0.3-1% in Vc, +-0.3-1 point in ee).
# Bak: (k/n) labels on the P < 1 markers; the P = 1 markers of every series share one "(66/66)" label, so drops = None.
# The SI (pubs.acs.org) returned 403; the main text gives no all-L / all-D counts.
BAK_MAIN = [(0.164, 1.0, None), (0.870, 0.80, 10), (1.726, 0.70, 10)]   # Fig. 2a blue squares
BAK_LD = ()
# Fig. 2a red circles (Tsat 50 C, Tinc 90 C, FC-40, PMMA), for the record; used only by the POST HOC line.
BAK_TINC90 = [(0.052, 1.0), (0.066, 1.0), (0.072, 1.0), (0.114, 1.0), (0.136, 1.0), (0.145, 1.0), (0.173, 1.0),
              (0.196, 10 / 12), (0.357, 5 / 6), (0.711, 2 / 6), (2.188, 5 / 16)]
# Asakura Fig. 2 plots |ee|; 150 C has 3 markers plus 7 at 0 ("Seven plots are thus superimposed at ee of 0%").
# Readings reproduce the printed Mean/SD (3/6 and 27/22).
AS_UNSTIRRED = {150: [17.2, 8.4, 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                152: [60.4, 49.1, 44.3, 43.4, 33.3, 19.2, 11.4, 9.4, 1.1, 0.0]}
# Fig. 4, ee > 0 = excess S.  0 rpm: separable markers only; a blob at 0 holds an uncountable number more (its edges
# are the +1.0 and -0.5 markers), and the paper does not say whether these are Fig. 2's 150 C runs.  500 rpm readings
# reproduce the printed means +70 and -72; no mean is printed at 200 rpm.
AS_STIRRED = {0: [17.4, 1.0, -0.5, -7.6],
              200: [75.5, 33.7, 25.4, 19.8, 6.5, -36.8, -62.1, -73.9, -77.7, -82.4],
              500: [77.7, 75.3, 69.1, 58.1, -61.5, -68.5, -73.9, -77.7, -78.5]}
# Table 1: "solid masses that grew from a single nucleus", run = crystallization number in the sample no.
AS_MASSES = [(150, 1, 0.0), (150, 1, 0.0), (150, 2, 0.0), (150, 3, 0.0),
             (152, 1, 0.0), (152, 1, -4.0), (152, 2, -59.0), (152, 3, -3.0), (152, 4, -2.0),
             (152, 5, 0.0), (152, 6, 31.0), (152, 7, 83.0), (152, 7, 0.0)]
AS_AGG = {150: 2.4, 152: 3.2}   # "total number of such aggregates in the 10 crystallizations at 150°C was 24"; 152 C: 32


def check(p, ok, detail):
    ok = ok if isinstance(ok, str) else bool(ok)
    RES[p] = ok
    print(f"  {p} {'HELD' if ok is True else 'FAILED' if ok is False else ok}: {detail}")


def free_width(V, P, n=None):
    """P = 1/(1 + (V/Vk)^h); x = (log10 Vk, log10 h), h clipped to [0.1, 50].  Returns (Vk, h)."""
    lv = np.log10(V)
    curve = lambda x: 1 / (1 + 10 ** np.clip(np.clip(10 ** x[1], 0.1, 50) * (lv - x[0]), -50, 50))
    if n is None:
        f = lambda x: float(np.sum((curve(x) - P) ** 2))
    else:
        f = lambda x: -loglik(P * n, n, curve(x))
    x, _ = best(f, [(np.median(lv), np.log10(h)) for h in (0.5, 1.5, 5.0, 20.0)])
    return 10 ** x[0], float(np.clip(10 ** x[1], 0.1, 50))


def exp_fit(V, P, n=None):
    """P = exp(-V/V0); returns V0."""
    curve = lambda x: np.exp(-V / 10 ** x[0])
    f = (lambda x: float(np.sum((curve(x) - P) ** 2))) if n is None else (lambda x: -loglik(P * n, n, curve(x)))
    x, _ = best(f, [(np.log10(np.median(V)),), (np.log10(V.max()),)])
    return 10 ** x[0]


def neff(ee_pct):
    m2 = float(np.mean((np.asarray(ee_pct, float) / 100) ** 2))
    return np.inf if m2 == 0 else 1 / m2


def neff_low(ee_pct):
    m = len(ee_pct)
    return neff(ee_pct) * chi2.ppf(0.025, m) / m


def setup_report():
    print(f"Bak, counting rule P_CSB = exp(-V/V0): 0.8 -> 0.2 over a factor {log(0.2) / log(0.8):.2f} in V, "
          f"0.9 -> 0.1 over {log(0.1) / log(0.9):.1f}; free-width h equivalent {log(16) / log(log(0.2) / log(0.8)):.2f}; "
          f"h = 3 is a factor {16 ** (1 / 3):.2f}")
    print(f"Asakura, 95% lower factor on n_eff for m = 10 runs: {chi2.ppf(0.025, 10) / 10:.3f}")
    for n in (20, 30, 40):
        print(f"  A1 power: an even split of {n} runs has P(split | Delta = 1) = {steer(n // 2, n):.1e}")


def main():
    # Self-check of the B2 fit: the counting-rule curve gives h near 1.4, a step gives h at the cap.
    V = np.array([0.1, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0])
    h_exp = free_width(V, np.exp(-V))[1]
    h_step = free_width(V, np.where(V < 1.2, 1.0, 0.0))[1]
    assert 1.0 < h_exp < 2.0 and h_step > 10, (h_exp, h_step)
    assert abs(exp_fit(V, np.exp(-V / 0.7)) - 0.7) < 1e-3 and abs(neff([10] * 4) - 100) < 1e-9
    setup_report()
    if any(x is None for x in (BAK_MAIN, BAK_LD, AS_UNSTIRRED, AS_STIRRED, AS_MASSES, AS_AGG)):
        print("Predictions committed; outcomes not yet entered (see the docstring).")
        return

    print("\n=== Bak et al. 2023, Fig. 2a main series ===")
    V = np.array([r[0] for r in BAK_MAIN], float)
    P = np.array([r[1] for r in BAK_MAIN], float)
    n = None if any(r[2] is None for r in BAK_MAIN) else np.array([r[2] for r in BAK_MAIN], float)
    for r in BAK_MAIN:
        print(f"  Vc = {r[0]:.3g} uL: P_CSB = {r[1]:.3f}, drops {r[2]}")
    if len(V) < 3:
        check("B1", "NOT SCORABLE", f"{len(V)} markers")
    else:
        rho = spearmanr(V, P)[0]
        check("B1", rho < 0, f"Spearman rho(Vc, P_CSB) = {rho:+.3f} over {len(V)} markers")
    if len(V) >= 3:
        print(f"  reported: exp fit V0 = {exp_fit(V, P, n):.3g} uL ({'binomial' if n is not None else 'least squares'})")
    if len(V) < 4 or P.max() < 0.6 or P.min() > 0.4:
        check("B2", "NOT SCORABLE", f"{len(V)} markers, P_CSB range {P.min():.2f}-{P.max():.2f}" if len(V) else "none")
    else:
        Vk, h = free_width(V, P, n)
        check("B2", h <= 3, f"free-width fit h = {h:.2f} (Vk = {Vk:.3g} uL); counting rule 1.40, threshold 3")
    if BAK_LD:
        a, b = BAK_LD
        check("B3", two_sided(a, a + b) > 0.05, f"all-L {a}, all-D {b}, p(50/50) = {two_sided(a, a + b):.3f}")
    else:
        check("B3", "NOT SCORABLE", "all-L / all-D counts not given")
    # POST HOC, added after the results were seen; not a prediction and not scored.
    Vr, Pr = np.array(BAK_TINC90).T
    print(f"  POST HOC, Tinc 90 C series ({len(Vr)} markers, least squares): h = {free_width(Vr, Pr)[1]:.2f} "
          f"(Vk = {free_width(Vr, Pr)[0]:.3g} uL), V0 = {exp_fit(Vr, Pr):.3g} uL")

    print("\n=== Asakura et al. 2002 ===")
    ne = {}
    for t, ee in AS_UNSTIRRED.items():
        ne[t] = neff(ee)
        print(f"  unstirred {t} C, {len(ee)} runs: ee % {ee}, n_eff = {ne[t]:.2f} (95% low {neff_low(ee):.2f})")
    for r, ee in AS_STIRRED.items():
        ne[f"{r} rpm"] = neff(ee)
        print(f"  {r} rpm, {len(ee)} runs: ee % {ee}, n_eff = {ne[f'{r} rpm']:.2f}")
    sg = [e for ee in AS_STIRRED.values() for e in ee if e != 0]
    a = sum(e > 0 for e in sg)
    if len(sg) < 20:
        check("A1", "NOT SCORABLE", f"{len(sg)} runs with a hand")
    else:
        check("A1", two_sided(a, len(sg)) > 0.05 and steer(a, len(sg)) < 1e-3,
              f"{a} S / {len(sg) - a} R, p(50/50) = {two_sided(a, len(sg)):.3f}, "
              f"P(split | Delta = 1) = {steer(a, len(sg)):.1e}")
    check("A2", ne["200 rpm"] < ne[150] and ne["500 rpm"] < ne[150],
          f"n_eff 200 rpm {ne['200 rpm']:.2f}, 500 rpm {ne['500 rpm']:.2f}, unstirred 150 C {ne[150]:.2f}")
    check("A3", ne["500 rpm"] <= 2, f"n_eff(500 rpm) = {ne['500 rpm']:.2f}")
    check("A4", ne[150] >= ne[152], f"n_eff 150 C {ne[150]:.2f}, 152 C {ne[152]:.2f}")
    for t, run, ee in AS_MASSES:
        print(f"  single mass {t} C run {run}: ee {ee:+.1f}%")
    if AS_MASSES:
        worst = min(abs(m[2]) for m in AS_MASSES)
        check("A5", worst >= 50, f"{len(AS_MASSES)} masses, smallest |ee| {worst:.1f}%")
    else:
        check("A5", "NOT SCORABLE", "no single-nucleus masses listed")
    if AS_AGG:
        parts = []
        for k, agg in AS_AGG.items():
            ee = AS_UNSTIRRED[k] if k != "pooled" else AS_UNSTIRRED[150] + AS_UNSTIRRED[152]
            parts.append((k, neff_low(ee), agg))
        check("A6", all(lo <= agg for _, lo, agg in parts),
              "; ".join(f"{k}: n_eff low {lo:.2f} vs aggregates {agg:.2f} per run" for k, lo, agg in parts))
    else:
        check("A6", "NOT SCORABLE", "no aggregate count given")
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}"
                                                         for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
