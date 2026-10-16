"""Replace the guessed network share f of plume_signal.py (log-uniform 1e-4..1) with one built from the model's own network
concentration and a literature abiotic background (2026-10-16).  Nothing in plume_signal.py is edited; this file imports it.

WHAT THE LITERATURE GIVES (every DOI below resolved on Crossref, 2026-10-16, no email in the request; abstracts read via PubMed)
 - Steel, Davila & McKay 2017, Astrobiology 17:862, doi:10.1089/ast.2017.1673 (abstract).  Serpentinization H2 of 0.6-34 mol/s
   supports abiotic amino-acid synthesis of 1.6-87 g/s and biotic synthesis of 1-44 g/s.  Gly, Ala, AIB and Glu "could all be
   above 0.01 uM just due to abiotic production" (a FLOOR, no ceiling is given).  Biology could reach at most 90 uM (a methanogenic
   ecosystem, not a Frank-type network).  Racemization is fast against production, so no ee in the ambient ocean; plume ee < 10%
   if biology is present.  NOTE: the DOI in the task note, 10.1089/ast.2015.1449, does not resolve on Crossref; 2017.1673 is right.
 - Postberg et al. 2018, Nature 558:564, doi:10.1038/s41586-018-0246-4 (abstract).  Macromolecular organics > 200 u in ice grains;
   a thin organic-rich film on the water table samples the organic inventory "in enhanced concentrations".  So grain organics are
   NOT the bulk ocean mixture, and no total organic carbon or amino-acid abundance is given.
 - Khawaja et al. 2019, MNRAS 489:5231, doi:10.1093/mnras/stz2280 and Khawaja et al. 2025, Nature Astronomy,
   doi:10.1038/s41550-025-02655-y (titles/dates checked on Crossref, abstracts not read in full): N-, O-bearing, aromatic and
   fresh-grain organics; I found no amino-acid abundance or total organic carbon in either.
 - Klenner et al. 2020, Astrobiology 20:1168, doi:10.1089/ast.2019.2188 (abstract): impact-ionization detection limits for
   biosignature molecules are at the uM or nM level.  That is an instrument limit, not an ocean concentration.
 RESULT OF THE SEARCH: no source measures how much of an Enceladus amino acid is made by an amplifying network, and no source
 gives the ocean's total organic carbon.  f is therefore NOT constrained by observation.  What can be done is to put the model's
 network product and the abiotic background in the same unit (molar concentration) and let f vary draw by draw.

THE MODEL SIDE.  sweep2.py (e): the feedstock concentration a = c (uncertainty.py: log-uniform 1e-12..1e-4 M for Enceladus,
ASSUMED; PAPER_DRAFT.md says Steel's >= 10 nM abiotic .. ~90 uM biotic "fall inside the model's range") sets k = k2 c; the
product p = r c, with r = p/a log-uniform 0.01..1 (ASSUMED, ex["r"]), sets the molecule count N = p V f 1000 N_A.  So the
network's product concentration is p = r c per draw, in M.
 f = p / (p + B), B = the abiotic, racemic concentration of the same amino acid.

ASSUMPTIONS (all stated, none measured)
 (a) c is a feedstock (achiral) concentration, separate from the abiotic amino acid B; the network product p adds to B.
     If the draws' c already contained B, f would be smaller.  Not tested here.
 (b) p is the standing network product concentration while it amplifies (it is the quantity that sets N); I treat it as the
     concentration that is also sampled in the plume.  The Postberg film enrichment is taken to act on both alike.
 (c) B is log-uniform between FLOOR 1e-8 M (Steel's 0.01 uM) and a ceiling derived from Steel's own maximum abiotic production
     87 g/s: all of it one amino acid of 75 g/mol, held for the pool's residence tau_res = 1e5 yr (the base tau of
     plume_signal.py) in the model's Enceladus volume 2.7e16 m^3 = 2.7e19 L:
     B_hi = 87/75 mol/s x 3.156e7 s/yr x 1e5 yr / 2.7e19 L = 1.36e-7 M.  The same formula with 1.6 g/s gives 2.5e-9 M (below
     the 0.01 uM floor, which came from Steel's own plume/vent dilution calculation, so I keep the floor).
     B is kept FIXED when tau_res changes below (a longer residence would raise B); a sensitivity scan covers 1e-9..1e-5 M.
 (d) The abiotic pool is racemic, so it only dilutes: ee_plume = f x eta* x s as in plume_signal.py (Sect. 2 there).
 (e) Everything else (gates A, eta*, s, T, Davila's rule) is exactly plume_signal.py's.  P(f >= T/(eta* s)) over the log-uniform B
     is computed in closed form (B <= p (eta* s / T - 1)), so no extra random numbers are drawn.
 (f) The prior c-range, not a measurement, then decides f: with p up to 1e-4 M and B ~ 1e-8..1e-7 M, f is small only for draws
     with small c.  The gates (k2 c t >= 10) favor large c.

Base: same draws as plume_signal.py (multiwall main rule, seed 0, N = 200,000, seedfix stream, ex streams [seed, i]).
Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python plume_share.py > plume_share.log

PREDICTIONS (written 2026-10-16, before any run of this script).  They come from the closed forms above and from a Monte Carlo
over the PRIOR ranges only (c, r, k2, t log-uniform, gates as in plume_signal.py, eta* s set to 0.94), not from the repo's
draws.  Old values are plume_signal.log (tau_res 1e5 yr, plume = lid layer, share of ALL lid draws).
 (P1) Regression: gates pass in 78.7 +- 0.2 % of draws, and the old log-uniform-f rule reproduces 0.360 / 0.231 (+-0.002) for
      T = 1 / 5 %.
 (P2) Derived ceiling B_hi = 1.36e-7 M (+-0.02e-7); the 1.6 g/s version 2.5e-9 M (+-0.2e-9).
 (P3) f over gate-passing draws (primary B): median 0.12-0.20; share with f < 0.01 is 0.25-0.38; share with f > 0.5 is
      0.32-0.45; share with p above Steel's biotic maximum 90 uM is below 0.001.
 (P4) New P(plume ee >= T), all lid draws, tau_res 1e5 yr, primary B:  T = 1%: 0.44-0.54; 2%: 0.40-0.50; 5%: 0.36-0.46;
      10%: 0.32-0.42.  Each is above the old value (0.360 / 0.304 / 0.231 / 0.176); at T = 5% the rise is at least 0.10.
 (P5) Over gate-passing draws only (T = 5%, tau 1e5 yr): 0.50-0.60 (cannot exceed 0.933 = share with eta* > 0).
 (P6) Wide B (log-uniform 1e-9..1e-6 M) moves the all-draw P(ee >= 5%) by less than 0.03 from the primary value.
 (P7) Fixed B scan, all draws, T = 5%, tau 1e5 yr: B = 1e-9: 0.55-0.65; 1e-8: 0.46-0.56; 1e-7: 0.35-0.45; 1e-6: 0.23-0.33;
      1e-5: 0.11-0.21; strictly decreasing in B.
 (P8) tau_res = 1e7 yr (s = 0.148), primary B, all draws: P(ee >= 5%) is 0.12-0.30, above the old 0.083; tau_res = 1e6 yr:
      0.33-0.43.
 (P9) Davila's rule at T = 5%, tau 1e5 yr: the life-read share (n >= 4) equals the abiotic-with-detectable-excess share (n <= 3)
      and both equal the P4 value for T = 5% (0.36-0.46).
"""
import os
import time

import numpy as np

import plume_signal as PS   # importing it patches default_rng (seedfix) and sets U.N exactly as plume_signal.py does

S, M, U = PS.S, PS.M, PS.U
YR = PS.YR
TS, TAUS = PS.TS, PS.TAUS
CHECK = []

B_FLOOR = 1e-8                                             # M, Steel et al. 2017: 0.01 uM
VOL_L = 2.7e16 * 1e3                                       # litres, uncertainty.py WORLDS["Enceladus"]
b_from_g = lambda g_s, tau_yr=1e5, mw=75.0: g_s / mw * YR * tau_yr / VOL_L   # M
B_HI = b_from_g(87.0)
B_PROD_LO = b_from_g(1.6)


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


def pdet_B(p, e, T, lo, hi):
    """P(f >= T/e) over B ~ log-uniform[lo, hi]; f = p/(p+B), e = eta* s.  f >= x <=> B <= p (1/x - 1), x = T/e."""
    with np.errstate(divide="ignore", invalid="ignore"):
        x = T / e
        bmax = p * (1.0 / x - 1.0)
        out = np.clip((np.log(bmax) - np.log(lo)) / (np.log(hi) - np.log(lo)), 0.0, 1.0)
    return np.where((e > 0) & (x < 1), np.nan_to_num(out), 0.0)


def pdet_fixed(p, e, T, B):
    with np.errstate(divide="ignore", invalid="ignore"):
        return ((e > 0) & (p / (p + B) >= T / e)).astype(float)


def main(seed=0):
    t0 = time.time()
    _, W = S.chains(seed)
    i = list(W).index("Enceladus")
    x, win, a, salt = W["Enceladus"]
    rng = np.random.default_rng([seed, i])
    ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
    L = a[1]
    c, k2, Dh, t = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time"))
    k = k2 * c
    krA = S.K_RAC["Ala"]
    A = (k * t >= 10) & (k > 2 * krA) & ((L / U.A) ** 2 / Dh <= t)
    eta = np.where(k > 0, PS.eta_table()(krA / k), 0.0)
    p = ex["r"] * c
    print(f"plume_share.py, N = {U.N}, seed 0, Enceladus; lid in {100 * salt.mean():.1f}% of draws; gates pass in {100 * A.mean():.1f}%.")

    print("\n(0) Regression: old log-uniform f, tau 1e5 yr, share of ALL draws (plume_signal.log: 0.360 / 0.304 / 0.231 / 0.176)")
    e5 = eta * PS.surv(krA, 1e5)
    old = {T: (A * PS.pdet(e5, T)).mean() for T in TS}
    print("  " + " | ".join(f"T {100 * T:g}%: {old[T]:.3f}" for T in TS))

    print("\n(1) Abiotic background B (M)")
    print(f"  floor (Steel 0.01 uM): {B_FLOOR:.2e}; ceiling from 87 g/s, 75 g/mol, tau 1e5 yr, V 2.7e19 L: {B_HI:.3e}; "
          f"same with 1.6 g/s: {B_PROD_LO:.2e}")

    print("\n(2) Network product p = r c (M), and f = p/(p+B), over gate-passing draws, B log-uniform [floor, ceiling]")
    print(f"  p percentiles 10/50/90 (M): " + " ".join(f"{np.percentile(p[A], q):.2e}" for q in (10, 50, 90)) +
          f"; share above 90 uM (Steel biotic max): {(p[A] > 9e-5).mean():.2e}; share below the floor: {(p[A] < B_FLOOR).mean():.3f}")
    Bm = np.sqrt(B_FLOOR * B_HI)
    fm = p / (p + Bm)    # f at the geometric-mid B, enough for a summary of the distribution
    print(f"  f at B = {Bm:.2e} (geometric mid): median {np.median(fm[A]):.3f}; share < 0.01: {(fm[A] < 0.01).mean():.3f}; "
          f"share < 1e-3: {(fm[A] < 1e-3).mean():.3f}; share > 0.5: {(fm[A] > 0.5).mean():.3f}; share > 0.9: {(fm[A] > 0.9).mean():.3f}")

    print("\n(3) NEW P(plume ee >= T), plume = lid layer. Over A draws | over all draws (B log-uniform [floor, ceiling])")
    new, newA = {}, {}
    for tau in TAUS:
        e = eta * PS.surv(krA, tau)
        row = []
        for T in TS:
            q = pdet_B(p, e, T, B_FLOOR, B_HI)
            newA[tau, T], new[tau, T] = q[A].mean(), (A * q).mean()
            row.append(f"T {100 * T:g}%: {newA[tau, T]:.3f} | {new[tau, T]:.3f}")
        print(f"   tau {tau:g} yr: " + " | ".join(row))
    print("  OLD (log-uniform f) for comparison, all draws: " + ", ".join(
        f"tau {tau:g}: " + "/".join(f"{(A * PS.pdet(eta * PS.surv(krA, tau), T)).mean():.3f}" for T in TS) for tau in TAUS))

    print("\n(4) Sensitivity, tau 1e5 yr, all draws, T = 1 / 2 / 5 / 10 %")
    wide = {T: (A * pdet_B(p, e5, T, 1e-9, 1e-6)).mean() for T in TS}
    print("  B log-uniform 1e-9..1e-6 M: " + " / ".join(f"{wide[T]:.3f}" for T in TS))
    fixed = {}
    for B in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        fixed[B] = {T: (A * pdet_fixed(p, e5, T, B)).mean() for T in TS}
        print(f"  B fixed {B:.0e} M: " + " / ".join(f"{fixed[B][T]:.3f}" for T in TS))

    print("\n(5) Davila's rule, tau 1e5 yr, T = 5%: n <= 3 read abiotic with a detectable excess; n >= 4 read life (share of ALL lid draws)")
    w5 = new[1e5, 0.05]
    print(f"  n=1..3: {100 * w5:.1f} / 0.0 / {100 * (1 - w5):.1f};  n=4..6: 0.0 / {100 * w5:.1f} / {100 * (1 - w5):.1f}")

    print("\nPredictions:")
    check("P1", abs(100 * A.mean() - 78.7) <= 0.2 and abs(old[0.01] - 0.360) <= 0.002 and abs(old[0.05] - 0.231) <= 0.002,
          f"gates {100 * A.mean():.1f}%; old {old[0.01]:.3f} / {old[0.05]:.3f}")
    check("P2", abs(B_HI - 1.36e-7) <= 0.02e-7 and abs(B_PROD_LO - 2.5e-9) <= 0.2e-9, f"B_hi {B_HI:.3e}, 1.6 g/s {B_PROD_LO:.2e}")
    ok = 0.12 <= np.median(fm[A]) <= 0.20 and 0.25 <= (fm[A] < 0.01).mean() <= 0.38 and 0.32 <= (fm[A] > 0.5).mean() <= 0.45 \
        and (p[A] > 9e-5).mean() < 1e-3
    check("P3", ok, f"median f {np.median(fm[A]):.3f}; f<0.01 {(fm[A] < 0.01).mean():.3f}; f>0.5 {(fm[A] > 0.5).mean():.3f}; "
                    f"p>90uM {(p[A] > 9e-5).mean():.1e}")
    r = {T: new[1e5, T] for T in TS}
    ok = 0.44 <= r[0.01] <= 0.54 and 0.40 <= r[0.02] <= 0.50 and 0.36 <= r[0.05] <= 0.46 and 0.32 <= r[0.10] <= 0.42 \
        and all(r[T] > old[T] for T in TS) and r[0.05] - old[0.05] >= 0.10
    check("P4", ok, "T 1/2/5/10%: " + " / ".join(f"{r[T]:.3f}" for T in TS) + "; old " + " / ".join(f"{old[T]:.3f}" for T in TS))
    check("P5", 0.50 <= newA[1e5, 0.05] <= 0.60 and newA[1e5, 0.05] <= (eta[A] > 0).mean(), f"over A, T 5%: {newA[1e5, 0.05]:.3f}")
    check("P6", abs(wide[0.05] - r[0.05]) < 0.03, f"wide B {wide[0.05]:.3f} vs primary {r[0.05]:.3f}")
    rng5 = {1e-9: (0.55, 0.65), 1e-8: (0.46, 0.56), 1e-7: (0.35, 0.45), 1e-6: (0.23, 0.33), 1e-5: (0.11, 0.21)}
    v = [fixed[B][0.05] for B in rng5]
    check("P7", all(lo <= fixed[B][0.05] <= hi for B, (lo, hi) in rng5.items()) and all(np.diff(v) < 0),
          "T 5% by B: " + " / ".join(f"{u:.3f}" for u in v))
    check("P8", 0.12 <= new[1e7, 0.05] <= 0.30 and new[1e7, 0.05] > (A * PS.pdet(eta * PS.surv(krA, 1e7), 0.05)).mean()
          and 0.33 <= new[1e6, 0.05] <= 0.43, f"tau 1e7: {new[1e7, 0.05]:.3f}; tau 1e6: {new[1e6, 0.05]:.3f}")
    check("P9", 0.36 <= w5 <= 0.46, f"n >= 4 life-read share {w5:.3f} = n <= 3 abiotic-read share (same quantity by construction)")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held. {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
