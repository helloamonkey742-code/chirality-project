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

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python blind_units.py > blind_units.log
"""
from math import log

import numpy as np
from scipy.stats import chi2, spearmanr

from blind_lab import best, loglik, steer, two_sided

RES = {}

# Outcomes: entered after the predictions commit by a separate reader.  None = not yet entered.
BAK_MAIN = None   # [(Vc in uL, P_CSB, drops in the marker or None)], Fig. 2a main series
BAK_LD = None     # (all-L drops, all-D drops), or () if not given
AS_UNSTIRRED = None   # {150: [ee %, 10 runs], 152: [...]}, Fig. 2
AS_STIRRED = None     # {0: [signed ee %], 200: [...], 500: [...]}, Fig. 4
AS_MASSES = None      # [(temperature C, run, ee %)], Table 1
AS_AGG = None         # {150: mean aggregates per run, 152: ...} or {"pooled": ...}; {} if not given


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
