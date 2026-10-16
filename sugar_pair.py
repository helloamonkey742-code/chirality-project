"""Does the model also explain D-sugars? A. Davila (NASA), email 2026-10-16: "Does your model also explain life's
preference for D-sugars (DNA/RNA), not just L-amino acids?"  Run through dna_network.py's coupled model (2026-10-16).

THE MODEL USED (imported, never edited): dna_network.py. Amino acids alpha (+ = L) and sugars beta (+ = D) cross the
bifurcation together; k_as = amino acid -> sugar transmission, k_sa = sugar/RNA -> amino acid transmission. The pair
mode is biased by g_eff = (sqrt(k_as) g_a + sqrt(k_sa) g_s) / sqrt(k_as + k_sa)   (dna_network.theoretical_p),
and the natural pair (L-aa + D-sugar) wins with probability Phi(Delta(g_eff)) once the pair is matched. With
r = g_s / g_a and rho = k_as / k_sa:   c = g_eff / g_a = (sqrt(rho) + r) / sqrt(rho + 1).
The natural pair is favoured iff c > 0, i.e. r > 0 (sugar PVED agrees with the amino-acid one) or rho > r^2 (r < 0).
(rho here is the coupling ratio; it is NOT the "rho" in multiwall.run4's extra inputs, which is a different quantity.)
Everything is conditioned on g_a > 0 (the PVED favours L-amino acids; its sign for amino acids in water is itself not
settled, Quack et al. 2022). If g_a's sign were 50/50 as well, the natural pair would be favoured in exactly 50% of
draws for any distribution of r and rho: the coupling explains the PAIRING of the hands, not which pair.

PRIORS (stated, not fitted).
 sign of r: primary 50/50. Evidence both ways: Zanasi et al. (1999, Advances in BioChirality ch. 29, Table 8) find the
 PVED of hydrated D-glyceraldehyde lowers D at all 12 tabulated conformers (same side as L-amino acids); MacDermott (1995)
 finds the natural enantiomer more stable in nearly all cases; but for ribose the sign follows the ring pucker (D for C2'-endo,
 L for C3'-endo; Tranter et al. 1992), and the DNA and RNA double helices are not favoured (Faglioni et al. 2005, 2006,
 abstracts only). Variants p(r > 0) = 0, 0.25, 0.75 are also reported.
 |r|: log-uniform 0.03-3. The only published magnitudes in hand are Zanasi's: D-glyceraldehyde / L-amino acid (PVED = 2 x
 tabulated shift, extended-basis RPA, preferred amino-acid conformers of Gly, Ala, Val, Asp) is computed in Section 1
 below. Tranter et al.'s magnitudes for ribose were not available. The range is Zanasi's ratio widened by 6 times on
 top and on the bottom by 1.3 times; 0.1-1 and 0.01-10 are also reported.
 rho: log-uniform 0.01-100 (symmetric about 1, no direction preferred; k_as / k_sa has never been measured). 0.1-10 and
 0.001-1000 are also reported.
 The coupling strength is not a free input here: any sqrt(k_as k_sa) >> sqrt(gamma) matches the hands (dual.py,
 dna_network.py), so only rho and r matter for which pair.

OCEANS. multiwall.run4 (seed 0, rules main and cascade x10, N = 200,000 per world, seedfix.py's stream correction
re-installed by the same four-line wrapper as in lid_melt.py and g_wide.py) is rerun with the PVED rescaled by |c|
(sweep2.G_NEW multiplied at run time, exactly as g_wide.py moved the top; the draws, rules and seeds are unchanged), for
|c| = 3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001. The sign of c only decides WHICH pair the PVED-set share belongs to:
natural pair share = P_pved x [c > 0] + P_chance / 2; mirror pair share = P_pved x [c < 0] + P_chance / 2. Averaged over
the primary draws of c by interpolating P_pved and P_chance in log10 |c| (flat beyond the grid ends).

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python sugar_pair.py > sugar_pair.log

PREDICTIONS (written 2026-10-16, before any run of this script; the numbers for P3-P6 are closed-form integrals done
by hand, P1/P8-P11 use committed logs only):
 (P1) Regression. (a) For the six spot-check cases of Section 3, c = (sqrt(rho)+r)/sqrt(rho+1) equals
      dna_network's g_eff / g_a to 1e-9. (b) At |c| = 1 the four worlds reproduce lid_melt.py's LOGGED multiwall main
      values (no lid / final) within 0.05 points: early Earth 99.9 / 99.9, Ganymede 44.1 / 44.1, Europa 39.6 / 30.1,
      Enceladus 26.4 / 0.0.
 (P2) Zanasi ratio r_Z = PVED(D-glyceraldehyde conformer) / mean PVED(preferred L-aa conformer): all 12 conformers have
      the same sign as the amino acids (D favoured); r_Z runs from 0.02 to 0.50 (min 0.015-0.035, max 0.40-0.55), mean
      0.22-0.30.
 (P3) Primary prior: the fraction of draws with c > 0 (natural pair favoured) is 0.859-0.869 (closed form 0.864:
      0.5 + 0.5 x 0.727, where 0.727 is the share of r < 0 draws with rho > r^2); MC (10^6 draws) and the closed form
      agree within 0.003.
 (P4) Sensitivity of that fraction (closed forms): rho 0.1-10: 0.875-0.887; rho 0.001-1000: 0.832-0.842; |r| 0.1-1:
      0.870-0.880; |r| 0.01-10: 0.828-0.838; p(r>0) = 0: 0.722-0.732; 0.25: 0.790-0.800; 0.75: 0.928-0.936; 1: 1.000.
      With sign 50/50 every range variant stays inside 0.83-0.89.
 (P5) Rescue condition: among r < 0 draws the median required rho (= r^2) is 0.08-0.10 and 23-25% of them need
      rho > 1 (|r| > 1).
 (P6) Near-cancellation: P(|c| < 0.1) is 2-8% and P(|c| < 0.01) is 0.2-1.2% of the primary draws.
 (P7) Spot checks of dna_network.simulate_asymmetric (6 cases, kappa_bar = 0.3, g_a = 1e-3, 4000 runs): the matched-pair
      rate is >= 0.99 in all six, and the natural-pair rate is within 0.04 of Phi(Delta(g_eff)) (theoretical_p) in all
      six; the case with c = 0 (r = -1, rho = 1) gives 0.50 +/- 0.04.
 (P8) Ocean invariants, all 15 worlds, both lid states, all 8 values of |c|: P_pved is non-decreasing in |c|; P_none is
      the same at every |c| (within 0.05 points); P_pved + P_chance is the same at every |c| (within 0.05).
 (P9) Ocean size, multiwall main P_pved, no lid, loss from |c| = 1: at |c| = 0.3 each of Ganymede, Europa, Enceladus
      loses 1-6 points and early Earth less than 1; at 0.1 each loses 3-13 points and early Earth less than 3; at 0.01
      each loses 8-27 points and early Earth stays at 90% or more. Enceladus with its lid is 0.0 at every |c|.
 (P10) Order. At every |c| in the grid the final (lid) order early Earth > Ganymede > Europa > Enceladus holds
      (multiwall main), and Ganymede is still at 25% or more at |c| = 0.01.
 (P11) Primary-prior ocean outcome, natural-pair share minus mirror-pair share (points; main, interpolated as above):
      early Earth 62-73, Ganymede 22-32, Europa with lid 14-24, Enceladus no lid 12-22, Enceladus with lid 0.0 (+/- 0.05);
      the natural pair is ahead in all four.
 (P12) Cost: the whole script finishes in under 15 minutes at 2 threads under nice.

DISCLOSURE. Before this commit I read dna_network.py, dual.py, sweep2.py, multiwall.py, g_wide.py/log, lid_melt.py, the
Zanasi chapter, the IJA draft and PAPER_DRAFT.md passages named in the brief; I did not run this script or any part of it
(no smoke run, no probe). The closed forms in P3-P6 were computed by hand.

RESULTS: (to be appended after the run, below this line; the predictions above are never edited)
"""
import os
import sys
import time

import numpy as np
from scipy.stats import norm

T0 = time.time()
ORIG_RNG = np.random.default_rng


def fixed_rng(seed=None, *a, **kw):   # seedfix.py's correction, copied (importing seedfix would run it)
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng
import uncertainty as U   # noqa: E402
U.N = int(os.environ.get("SUGAR_N", U.N))   # smoke test only
import sweep2 as S        # noqa: E402
import multiwall as M     # noqa: E402
import dna_network as D   # noqa: E402
from sim import delta     # noqa: E402

FULL = U.N == 200_000
NDRAW = int(os.environ.get("SUGAR_NDRAW", 1_000_000))
RUNS = int(os.environ.get("SUGAR_RUNS", 4000))
RES = {}
PRIMARY = dict(p_plus=0.5, rabs=(0.03, 3.0), rrho=(0.01, 100.0))


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def lap(lab):
    print(f"[{lab}: {time.time() - T0:.0f} s elapsed]")


def cfac(r, rho):
    """c = g_eff / g_a for r = g_s / g_a and rho = k_as / k_sa."""
    return (np.sqrt(rho) + r) / np.sqrt(rho + 1.0)


# ---------------------------------------------------------------- 1. what the literature gives for r
HARTREE_AEV = 27.211386245988e18
ZAN_AA = {   # Zanasi et al. 1999, Tables 3-6, extended-basis RPA, Epv / 1e-20 hartree, L-enantiomer, all angles (as in g_wide.py)
    "glycine":   [-12.733, -14.179, -8.692, 0.000, 14.179, 12.733, 7.240, 0.000],
    "alanine":   [-12.318, -13.170, -10.977, -5.238, 2.949, 14.873, 14.602, 10.245, 3.476, -3.490],
    "valine":    [-16.525, -9.107, -0.761, 9.364, 20.859, 10.619, -9.064],
    "aspartate": [-13.754, -11.028, 4.709, 16.759, 10.383, -4.749],
}
ZAN_GLYC = [-5.024, -5.226, -6.050, -4.686, -2.904, -1.741, -0.742, -0.338, -4.961, -6.804, -3.368, -2.589]  # Table 8, D-glyceraldehyde


def to_aev(shift):   # PVED = 2 x shift (chapter 29, p. 377)
    return 2 * abs(shift) * 1e-20 * HARTREE_AEV


def section1():
    print("=== 1. Sugar / amino-acid PVED ratio in the only calculation in hand (Zanasi et al. 1999) ===")
    pref = [to_aev(min(v)) for v in ZAN_AA.values()]
    ref = float(np.mean(pref))
    print("  preferred L-aa conformers (aeV): " + ", ".join(f"{k} {p:.2f}" for k, p in zip(ZAN_AA, pref)) + f"; mean {ref:.2f}")
    gl = np.array([to_aev(v) for v in ZAN_GLYC])
    same = all(v < 0 for v in ZAN_GLYC)   # negative shift = the tabulated (natural, D) enantiomer is lower
    r_z = gl / ref
    print(f"  D-glyceraldehyde, 12 conformers: {gl.min():.2f} to {gl.max():.2f} aeV, mean {gl.mean():.2f}; all negative (D lower): {same}")
    print(f"  r_Z = glyceraldehyde / amino acid: min {r_z.min():.3f}, mean {r_z.mean():.3f}, max {r_z.max():.3f}  (same sign as the amino acids)")
    check("P2", same and 0.015 <= r_z.min() <= 0.035 and 0.40 <= r_z.max() <= 0.55 and 0.22 <= r_z.mean() <= 0.30,
          f"same sign {same}; min {r_z.min():.3f}, max {r_z.max():.3f}, mean {r_z.mean():.3f}")
    return r_z


# ---------------------------------------------------------------- 2. which pair is favoured
def draw(n, rng, p_plus, rabs, rrho):
    sgn = np.where(rng.random(n) < p_plus, 1.0, -1.0)
    r = sgn * 10 ** rng.uniform(np.log10(rabs[0]), np.log10(rabs[1]), n)
    rho = 10 ** rng.uniform(np.log10(rrho[0]), np.log10(rrho[1]), n)
    return r, rho


def exact(p_plus, rabs, rrho):
    """Closed form: p_plus + (1 - p_plus) * P(rho > r^2 | r < 0), by midpoint integration over log10 |r|."""
    u = np.log10(rabs[0]) + (np.arange(200000) + 0.5) / 200000 * (np.log10(rabs[1]) - np.log10(rabs[0]))
    lo, hi = np.log10(rrho[0]), np.log10(rrho[1])
    return p_plus + (1 - p_plus) * float(np.clip((hi - 2 * u) / (hi - lo), 0, 1).mean())


def section2():
    rng = ORIG_RNG(20261016)
    print("\n=== 2. Fraction of draws in which the coupled system favours the natural pair (c > 0), g_a > 0 ===")
    print(f"{NDRAW} draws per row; sign of r, |r| log-uniform, rho = k_as/k_sa log-uniform")
    print(f"{'variant':44} {'MC':>7} {'closed':>7}")
    rows = [("PRIMARY p(r>0)=0.5, |r| 0.03-3, rho 0.01-100", PRIMARY),
            ("rho 0.1-10", {**PRIMARY, "rrho": (0.1, 10.0)}), ("rho 0.001-1000", {**PRIMARY, "rrho": (1e-3, 1e3)}),
            ("|r| 0.1-1", {**PRIMARY, "rabs": (0.1, 1.0)}), ("|r| 0.01-10", {**PRIMARY, "rabs": (0.01, 10.0)}),
            ("p(r>0)=0 (every sugar PVED opposes)", {**PRIMARY, "p_plus": 0.0}), ("p(r>0)=0.25", {**PRIMARY, "p_plus": 0.25}),
            ("p(r>0)=0.75", {**PRIMARY, "p_plus": 0.75}), ("p(r>0)=1 (sugar agrees, as in Zanasi)", {**PRIMARY, "p_plus": 1.0})]
    out = {}
    for lab, kw in rows:
        r, rho = draw(NDRAW, rng, **kw)
        mc, cf = float((cfac(r, rho) > 0).mean()), exact(**kw)
        out[lab] = (mc, cf)
        print(f"{lab:44} {mc:7.4f} {cf:7.4f}")
    mc, cf = out[rows[0][0]]
    check("P3", 0.859 <= mc <= 0.869 and 0.859 <= cf <= 0.869 and abs(mc - cf) < 0.003, f"MC {mc:.4f}, closed {cf:.4f}")
    g = lambda i: out[rows[i][0]][1]
    ok = (0.875 <= g(1) <= 0.887 and 0.832 <= g(2) <= 0.842 and 0.870 <= g(3) <= 0.880 and 0.828 <= g(4) <= 0.838
          and 0.722 <= g(5) <= 0.732 and 0.790 <= g(6) <= 0.800 and 0.928 <= g(7) <= 0.936 and abs(g(8) - 1) < 1e-9
          and all(0.83 <= g(i) <= 0.89 for i in range(5)))
    check("P4", ok, "closed forms " + ", ".join(f"{g(i):.4f}" for i in range(9)))

    # the primary draws again, for the rescue condition and the small-|c| tail
    r, rho = draw(NDRAW, rng, **PRIMARY)
    c = cfac(r, rho)
    neg = r < 0
    med, big = float(np.median(r[neg] ** 2)), float((np.abs(r[neg]) > 1).mean())
    print("\n  Rescue condition for an opposing sugar PVED (r < 0): the natural pair needs rho = k_as/k_sa > r^2")
    print(f"  {'r':>7} {'rho needed':>11}   c at rho = 1")
    for rv in (-0.03, -0.1, -0.3, -1.0, -3.0):
        print(f"  {rv:7.2f} {rv ** 2:11.4f}   {cfac(rv, 1.0):7.3f}")
    print(f"  among r < 0 draws: median rho needed {med:.4f}; share needing rho > 1: {big:.3f}; "
          f"share rescued: {float((c[neg] > 0).mean()):.3f}")
    check("P5", 0.08 <= med <= 0.10 and 0.23 <= big <= 0.25, f"median {med:.4f}, share |r|>1 {big:.3f}")
    a1, a2 = float((np.abs(c) < 0.1).mean()), float((np.abs(c) < 0.01).mean())
    q = np.percentile(np.abs(c), [10, 25, 50, 75, 90])
    print(f"  |c| percentiles 10/25/50/75/90: " + " ".join(f"{v:.3f}" for v in q) + f";  P(|c|<0.1) = {a1:.4f}, P(|c|<0.01) = {a2:.4f}")
    check("P6", 0.02 <= a1 <= 0.08 and 0.002 <= a2 <= 0.012, f"P(|c|<0.1) {a1:.4f}, P(|c|<0.01) {a2:.4f}")
    return c


# ---------------------------------------------------------------- 3. Monte-Carlo spot checks of dna_network
def section3():
    print(f"\n=== 3. dna_network.simulate_asymmetric vs Phi(Delta(g_eff)), {RUNS} runs/case, g_a = 1e-3, kappa_bar = 0.3 ===")
    print(f"{'r':>6} {'rho':>7} {'c':>7} {'k_as':>6} {'k_sa':>6} {'predicted':>9} {'simulated':>9} {'matched':>8}")
    g_a, kb, eps, gamma = 1e-3, 0.3, 1e-4, 0.01
    cases = [(0.3, 1.0), (-0.3, 1.0), (-1.0, 1.0), (-1.0, 0.25), (-0.3, 16.0), (-0.5, 1 / 16)]
    ok_c, ok_m, ok_p, zero = True, True, True, None
    for i, (r, rho) in enumerate(cases):
        k_as, k_sa = kb * np.sqrt(rho), kb / np.sqrt(rho)
        thy = D.theoretical_p(k_sa, k_as, g_a, r * g_a, eps, gamma)
        c = cfac(r, rho)
        geff = (np.sqrt(k_as) * g_a + np.sqrt(k_sa) * r * g_a) / np.sqrt(k_as + k_sa)
        ok_c &= abs(geff / g_a - c) < 1e-9 and abs(float(norm.cdf(delta(c * g_a, eps, gamma))) - thy) < 1e-9
        m, p = D.simulate_asymmetric(k_sa, k_as, g_a, r * g_a, eps, gamma, runs=RUNS, seed=i + 1)
        print(f"{r:6.2f} {rho:7.3f} {c:7.3f} {k_as:6.3f} {k_sa:6.3f} {thy:9.3f} {p:9.3f} {m:8.3f}")
        ok_m &= m >= 0.99
        ok_p &= abs(p - thy) < 0.04
        if r == -1.0 and rho == 1.0:
            zero = p
    check("P1a", ok_c, "c equals dna_network's g_eff / g_a (and Phi(Delta(c g_a)) equals theoretical_p) to 1e-9 in all six cases")
    check("P7", ok_m and ok_p and zero is not None and abs(zero - 0.5) < 0.04,
          f"matched >= 0.99: {ok_m}; natural within 0.04 of prediction: {ok_p}; c = 0 case {zero:.3f}")


# ---------------------------------------------------------------- 4. oceans with the PVED scaled by |c|
GRID = (3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)
FOUR = ("Early Earth ocean", "Ganymede", "Europa", "Enceladus")
LOGGED = {"Early Earth ocean": (99.9, 99.9), "Ganymede": (44.1, 44.1), "Europa": (39.6, 30.1), "Enceladus": (26.4, 0.0)}


def section4():
    print("\n=== 4. Oceans with the PVED scaled by |c| (multiwall.run4, seed 0, N = %d; Enceladus on stream [s, 0, 1]) ===" % U.N)
    G0 = S.G_NEW
    runs = {}
    for cs in GRID:
        S.G_NEW = (G0[0] * cs, G0[1] * cs)
        runs[cs] = M.run4(0, full=True)
        lap(f"|c| = {cs}")
    S.G_NEW = G0
    names = list(runs[1.0])
    m = lambda cs, n, k="main", lid=0: runs[cs][n][k][lid]   # (P_pved, P_chance, P_none), percent

    print("\nP_pved %, multiwall main, no lid / final")
    print(f"{'world':18}" + "".join(f"{('|c|=%g' % cs):>14}" for cs in GRID))
    for n in sorted(names, key=lambda n: -m(1.0, n, lid=1)[0]):
        print(f"{n:18}" + "".join(f"{m(cs, n)[0]:6.1f} /{m(cs, n, lid=1)[0]:5.1f}" for cs in GRID))
    print("\nThe four worlds, three-way split P_pved / P_chance / P_none, % (no lid), and cascade x10 P_pved (no lid / final)")
    for n in FOUR:
        print(n)
        for cs in GRID:
            a, b, d = m(cs, n)
            x0, x1 = m(cs, n, "x10", 0)[0], m(cs, n, "x10", 1)[0]
            print(f"   |c| = {cs:<6g} {a:5.1f} / {b:5.1f} / {d:5.1f}    x10 {x0:5.1f} / {x1:5.1f}")

    # P1b regression
    if FULL:
        bad = [(n, LOGGED[n], (round(m(1.0, n)[0], 1), round(m(1.0, n, lid=1)[0], 1))) for n in FOUR
               if abs(m(1.0, n)[0] - LOGGED[n][0]) > 0.05 or abs(m(1.0, n, lid=1)[0] - LOGGED[n][1]) > 0.05]
        check("P1b", not bad, "|c| = 1 reproduces lid_melt.py LOGGED for the four worlds" if not bad else f"mismatch {bad}")
    else:
        print("  (smoke run: P1b not evaluated)")

    # P8 invariants
    bad = []
    for n in names:
        for lid in (0, 1):
            pv = [m(cs, n, lid=lid)[0] for cs in sorted(GRID)]
            if any(b < a - 1e-9 for a, b in zip(pv, pv[1:])):
                bad.append((n, lid, "not monotone"))
            for cs in GRID:
                if abs(m(cs, n, lid=lid)[2] - m(1.0, n, lid=lid)[2]) > 0.05:
                    bad.append((n, lid, cs, "P_none"))
                if abs(m(cs, n, lid=lid)[0] + m(cs, n, lid=lid)[1] - m(1.0, n, lid=lid)[0] - m(1.0, n, lid=lid)[1]) > 0.05:
                    bad.append((n, lid, cs, "P_pved+P_chance"))
    check("P8", not bad, "15 worlds x 2 lid states x 8 values: monotone, P_none and P_pved+P_chance constant" if not bad else str(bad[:6]))

    # P9 sizes
    loss = lambda cs, n: m(1.0, n)[0] - m(cs, n)[0]
    three = ("Ganymede", "Europa", "Enceladus")
    ee = "Early Earth ocean"
    print("\nLoss in P_pved (points, no lid) from |c| = 1: " + "; ".join(
        f"|c|={cs:g}: " + ", ".join(f"{n.split()[0]} {loss(cs, n):.1f}" for n in FOUR) for cs in (0.3, 0.1, 0.01)))
    lids = [m(cs, "Enceladus", lid=1)[0] for cs in GRID]
    ok = (all(1 <= loss(0.3, n) <= 6 for n in three) and loss(0.3, ee) < 1
          and all(3 <= loss(0.1, n) <= 13 for n in three) and loss(0.1, ee) < 3
          and all(8 <= loss(0.01, n) <= 27 for n in three) and m(0.01, ee)[0] >= 90 and max(lids) == 0.0)
    check("P9", ok, "see the loss line above; Enceladus with lid max over the grid %.1f" % max(lids))

    # P10 order
    order_ok = all(m(cs, ee, lid=1)[0] > m(cs, "Ganymede", lid=1)[0] > m(cs, "Europa", lid=1)[0] > m(cs, "Enceladus", lid=1)[0]
                   for cs in GRID)
    check("P10", order_ok and m(0.01, "Ganymede", lid=1)[0] >= 25,
          f"order at every |c|: {order_ok}; Ganymede final at 0.01 = {m(0.01, 'Ganymede', lid=1)[0]:.1f}")
    return runs


# ---------------------------------------------------------------- 5. what the primary prior does to the oceans
def section5(runs, c):
    print("\n=== 5. Primary prior (Section 2 draws of c) applied to the oceans: natural-pair vs mirror-pair share, % (main) ===")
    lg = np.log10(np.array(sorted(GRID)))
    lc = np.log10(np.abs(c))
    pos = c > 0
    diffs = {}
    print(f"{'world':24} {'P_pved(|c|=1)':>14} {'natural':>8} {'mirror':>8} {'neither':>8} {'diff':>7}")
    for n in FOUR:
        for lid, tag in ((0, "no lid"), (1, "lid")):
            pv = np.interp(lc, lg, [runs[cs][n]["main"][lid][0] for cs in sorted(GRID)])
            ch = np.interp(lc, lg, [runs[cs][n]["main"][lid][1] for cs in sorted(GRID)])
            nat = float((pv * pos + ch / 2).mean())
            mir = float((pv * ~pos + ch / 2).mean())
            non = float(100 - nat - mir)
            diffs[(n, lid)] = nat - mir
            print(f"{n + ' (' + tag + ')':24} {runs[1.0][n]['main'][lid][0]:14.1f} {nat:8.1f} {mir:8.1f} {non:8.1f} {nat - mir:7.1f}")
    d = lambda n, lid: diffs[(n, lid)]
    ok = (62 <= d("Early Earth ocean", 1) <= 73 and 22 <= d("Ganymede", 1) <= 32 and 14 <= d("Europa", 1) <= 24
          and 12 <= d("Enceladus", 0) <= 22 and abs(d("Enceladus", 1)) <= 0.05
          and all(d(n, 1) > 0 for n in ("Early Earth ocean", "Ganymede", "Europa")) and d("Enceladus", 0) > 0)
    check("P11", ok, "diff: Earth %.1f, Ganymede %.1f, Europa(lid) %.1f, Enceladus(no lid) %.1f, Enceladus(lid) %.2f" % (
        d("Early Earth ocean", 1), d("Ganymede", 1), d("Europa", 1), d("Enceladus", 0), d("Enceladus", 1)))


def main():
    print(f"sugar_pair.py, N = {U.N} per world, {NDRAW} prior draws, {RUNS} runs per spot check")
    section1()
    c = section2()
    lap("sections 1-2")
    section3()
    lap("section 3")
    runs = section4()
    section5(runs, c)
    t = time.time() - T0
    check("P12", t < 900, f"{t:.0f} s")
    held = sum(RES.values())
    print(f"\n{held} of {len(RES)} checked items HELD (P1 is two items: P1a, P1b): "
          + ", ".join(f"{k} {'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items(), key=lambda kv: (len(kv[0]), kv[0]))))


if __name__ == "__main__":
    main()
