"""Bound the network's feedstock concentration c with Enceladus data, then rerun the gates and the plume-ee question (2026-10-16).
Nothing existing is edited; this file imports plume_signal.py and plume_share.py.

WHAT c IS (read from the code, not guessed).  uncertainty.py draws "concentration" c log-uniform 1e-12..1e-4 M for Enceladus
(ASSUMED).  In ocean.py and sweep2.py c is the concentration of the ACHIRAL FEEDSTOCK A of the Frank network (PAPER_DRAFT.md
Sect. 2.4: "A is the achiral feedstock (for example an amino-acid precursor)"): the copying rate is k = k2 c, the
molecule count is N = p V f 1000 N_A with the product concentration p = r c (r = p/a, log-uniform 0.01..1, ASSUMED, sweep2 step e).
c therefore enters (i) the gates k2 c t >= 10 and k2 c > 2 k_rac(Ala), (ii) Delta = F * sqrt(2) pi^(1/4) g (k2 c t)^(1/4) (p V f 1000 N_A)^(1/2)
(the PVED hand), (iii) eta*(k_rac/(k2 c)) (network purity), and (iv) p = r c, which sets the share f = p/(p+B) in plume_share.py.
The lateral-mixing gate (L/A)^2/D_h <= t does not depend on c.  The paper (PAPER_DRAFT.md, Sect. 6) identifies c with the Enceladus
amino-acid-class concentration ("the model's range"), so c and the abiotic racemic background B of plume_share.py are bounded by
the SAME chemistry.  If instead A is an upstream precursor (HCN, aldehydes, NH3), the data bound it only loosely (see below).

SOURCES (every DOI resolved on Crossref 2026-10-16 with a plain GET, no email; abstracts read, nothing downloaded)
 - Steel, Davila & McKay 2017, Astrobiology 17:862, doi:10.1089/ast.2017.1673 (abstract, via Europe PMC): serpentinization H2 of
   0.6-34 mol/s supports abiotic amino-acid synthesis of 1.6-87 g/s; Gly, Ala, AIB and Glu "could all be above 0.01 uM just
   due to abiotic production" (a floor); biology "a maximum of 90 uM"; racemization short against production.
 - Peter et al. 2024, Nature Astronomy 8:164, doi:10.1038/s41550-023-02160-0 (journal DOI checked; numbers read from the open arXiv
   preprint 2301.05259 page, so they may differ from the journal version): plume-gas mixing ratios HCN 0.11 +- 0.02 %, CH4 0.11 %,
   CO2 0.45 %, NH3 1.8 %; earlier HCN limit 0.01-0.2 %; formaldehyde < 0.004 % and acetaldehyde < 0.002 % (not detected).  The
   paper gives NO ocean concentration, and says the HCN may come in part from solid material in the ice shell.
 - Postberg et al. 2018, Nature 558:564, doi:10.1038/s41586-018-0246-4 (abstract): complex organics > 200 u in ice grains; an
   organic-rich film on the water table samples the organic inventory "in enhanced concentrations".  No abundance given.
 - Postberg et al. 2023, Nature 618:489, doi:10.1038/s41586-023-05987-9 (abstract): sodium phosphates; P in plume-forming water
   at least 100-fold Earth's ocean.  Not an amino-acid or organic-carbon constraint.
 - Hao et al. 2022, PNAS, doi:10.1073/pnas.2201388119 (abstract): modelled dissolved inorganic P 1e-7..1e-2 mol/kg; lists organic
   carbon only as a building block, no value.
 - Glein & Waite 2020, GRL, doi:10.1029/2019GL085885 (abstract, Crossref): activity of aqueous CO2 10^-4.6..10^-3.2 (2.5e-5..6.3e-4).
   Inorganic carbon, not an organic feedstock.
 - Waite et al. 2017, Science 356:155, doi:10.1126/science.aai8703 (abstract): H2 in the plume, hydrothermal disequilibrium.
   Khawaja et al. 2019 (MNRAS, doi:10.1093/mnras/stz2280) and 2025 (Nature Astronomy, doi:10.1038/s41550-025-02655-y): organic
   species in grains, abstracts give no abundance.
 RESULT OF THE SEARCH: no source gives the ocean's dissolved organic carbon, an amino-acid abundance, or an HCN or aldehyde
 concentration in the ocean.  The only quantitative handles are Steel's amino-acid floor (0.01 uM) and biotic cap (90 uM), and the
 amino-acid production rate 1.6-87 g/s.  A hard mass-balance ceiling for an HCN-class precursor is x_HCN x 55.5 M (vapor mole fraction
 times the molarity of water, if all plume HCN were dissolved in the water that became vapor) = 0.0011 x 55.5 = 0.061 M, 600x
 ABOVE the model's old ceiling of 1e-4 M, so the HCN data do not tighten c from above.  A Henry's-law conversion would be
 model dependent by two orders of magnitude (Peter et al. suggest non-equilibrium, ice-shell sourcing), so it is not used.

DERIVED c RANGES (log-uniform, same-quantile remap of the existing draws by sweep2.remap, so every other input and the common
random numbers are unchanged; the old draws are the S0 regression)
 S0  old, assumed:                            1e-12 .. 1e-4  M   (plume_share.py's basis)
 S1  PRIMARY, c = the amino-acid-class feedstock held to the abiotic window of Steel et al.:
        floor 1e-8 M (Steel: 0.01 uM), ceiling B_HI = 1.356e-7 M (87 g/s as one 75 g/mol amino acid, residence 1e5 yr, 2.7e19 L;
        plume_share.py).  B uses the same window.  This is the paper's own reading of c, made consistent with B.
 S2  ceiling only (the feedstock has no floor, e.g. an upstream precursor):  1e-12 .. 1.356e-7 M
 S3  Steel's full amino-acid window, floor to biotic cap:                     1e-8 .. 9e-5 M   (B stays at its S1 window)
 S4  S1 with a 100x longer residence of the feedstock pool (1e7 yr; the ceiling scales with residence):  1e-8 .. 1.356e-5 M, and
        B's ceiling is scaled the same way (1e-8 .. 1.356e-5 M)
ASSUMPTIONS (stated, none measured): (a) c and the product p = r c are separate from the racemic abiotic pool B, as in
 plume_share.py; if the draws' c already contained B, f would be smaller.  (b) The floor of S1 is justified only if the feedstock
 is an abiotically made amino-acid-class pool; for a precursor S2 is the honest case.  (c) The ceiling treats all 87 g/s as one
 amino acid and uses residence 1e5 yr, both of which only raise it.  (d) The remaining rules (gates, eta*, s, T) are plume_signal.py's.
 (e) B is kept at its S1 window in S2 and S3.

Base: multiwall main rule, seed 0, N = 200,000, seedfix stream, ex streams [seed, i], exactly as plume_signal.py / plume_share.py.
Reported per scenario: gate-pass share (A); share of A draws with eta* > 0; mean P(PVED hand | A) for the lid layer (f_lid 0.1) and the
bulk; whole-ocean PVED-set share, no lid, multiwall main rule (and with the lid, which is 0); P(plume ee >= 1/2/5/10 %) at tau_res
1e5 yr over all draws and over A, plume = lid layer; and f at B's geometric mid.
Run: cd ~/Documents/chirality-project && OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python plume_feedstock.py > plume_feedstock.log

PREDICTIONS (written 2026-10-16, before any run of this script).  They come from a Monte Carlo over the PRIOR ranges only (400,000
draws of c, k2, t, g, F, V, r with the model's formulas for the gates, eta* >= 0.97 where k2 c > 11.5 k_rac, and the closed form for
P(f >= T/(eta* s))), not from the repo's draws.  That Monte Carlo reproduced S0 (gates 0.785, P(ee >= 5%) 0.451, PVED|A layer
0.928) to within 0.003 of the logged values.  Old values: plume_signal.log and plume_share.log (tau 1e5 yr, share of ALL draws).
 (P1) Regression, S0: gates 78.7 +- 0.2 %; P(ee >= 1 / 5 %) = 0.522 / 0.449 (+-0.003); PVED|A layer 0.928 (+-0.003); no-lid
      PVED-set share 26.36 % (+-0.05).
 (P2) S1 (primary): the gates pass MORE often than in S0, in 0.95-0.99 of draws (the floor removes the low-c draws that fail them).
 (P3) S1: P(ee >= 1 / 2 / 5 / 10 %), all draws, tau 1e5 yr: 0.78-0.88 / 0.70-0.79 / 0.54-0.63 / 0.40-0.48, each above the S0 value
      (0.522 / 0.492 / 0.449 / 0.413).
 (P4) S1: mean P(PVED hand | A): lid layer 0.96-0.99, bulk 0.98-1.00, both above S0 (0.928 / 0.956).
 (P5) S1: whole-ocean no-lid PVED-set share (multiwall main) is 10-28 %; with the lid it is 0.  (No Monte Carlo for the layering,
      so this range is wide.)
 (P6) S2 (ceiling only): gates 0.62-0.71; P(ee >= 1 / 5 %) = 0.21-0.31 / 0.11-0.20, below S0; PVED|A layer 0.83-0.90.
 (P7) S3 (Steel's full window): gates 0.98-1.00; P(ee >= 1 / 5 %) = 0.92-0.98 / 0.83-0.91.
 (P8) S4 (1e7-yr residence): P(ee >= 5 %) = 0.52-0.61; gates 0.97-1.00.
 (P9) Ordering of P(ee >= 5 %): S2 < S0 < min(S1, S4) and max(S1, S4) < S3.  Ordering of the no-lid PVED-set share: S2 < S1.
 (P10) S1, T = 5 %, tau_res 1e6 / 1e7 yr: 0.45-0.54 / 0.13-0.22.
 (P11) S1, T = 5 %, fixed B: 1e-9 M: 0.88-0.94; 1e-7 M: 0.34-0.42; 1e-5 M: below 0.02; strictly decreasing in B.
 (P12) S1, f at B = 3.68e-8 M over A draws: median 0.07-0.12; share f < 0.01: 0.05-0.09; share f > 0.5: 0.05-0.10.
 (P13) HCN mass-balance ceiling 0.0011 x 55.5 = 0.061 M (+-0.001), more than 100x the old c ceiling 1e-4 M.
"""
import time

import numpy as np
from scipy.special import ndtr

import plume_signal as PS    # patches default_rng (seedfix) and sets U.N exactly as plume_signal.py does
import plume_share as SH

S, M, U = PS.S, PS.M, PS.U
TS, TAUS = PS.TS, PS.TAUS
CHECK = []
C_OLD = (1e-12, 1e-4)
B1 = (SH.B_FLOOR, SH.B_HI)
SCEN = {  # name: (c range, B range)
    "S0 old": (C_OLD, B1),
    "S1 primary": ((1e-8, SH.B_HI), B1),
    "S2 ceiling only": ((1e-12, SH.B_HI), B1),
    "S3 Steel full": ((1e-8, 9e-5), B1),
    "S4 residence 1e7": ((1e-8, SH.B_HI * 100), (1e-8, SH.B_HI * 100)),
}
HCN_X, WATER_M = 0.0011, 55.5


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


def main(seed=0):
    t0 = time.time()
    _, W = S.chains(seed)
    i = list(W).index("Enceladus")
    x, win, a, salt = W["Enceladus"]
    L, dz = a[1], a[4]
    rng = np.random.default_rng([seed, i])
    ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
    none = np.zeros(U.N, bool)
    krA = S.K_RAC["Ala"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    print(f"plume_feedstock.py, N = {U.N}, seed 0, Enceladus (index {i}); c enters the gates, Delta, eta* and p = r c.")
    print(f"HCN mass-balance ceiling: {HCN_X} x {WATER_M} M = {HCN_X * WATER_M:.3f} M = {HCN_X * WATER_M / C_OLD[1]:.0f}x the old c ceiling.")

    R = {}
    print("\nPer scenario (tau_res 1e5 yr; plume = lid layer, f_lid 0.1; P(ee >= T) as all draws | over A draws)")
    for name, (cr, br) in SCEN.items():
        x2 = dict(x)
        x2["concentration"] = S.remap(x["concentration"], C_OLD, cr)
        c, k2, Dh, t, V, F = (x2[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
        k = k2 * c
        A = (k * t >= 10) & (k > 2 * krA) & ((L / U.A) ** 2 / Dh <= t)
        eta = np.where(k > 0, PS.eta_table()(krA / k), 0.0)
        p = ex["r"] * c
        phi = lambda f: ndtr(F * M.delta_at(c, g, k2, V * ex["r"], t, f))
        pt, pb = phi(0.1)[A].mean(), phi(0.9)[A].mean()
        w0, _, _ = M.multi(x2, L, dz, ex, none)
        wl, _, _ = M.multi(x2, L, dz, ex, salt)
        e5 = eta * PS.surv(krA, 1e5)
        q = {T: SH.pdet_B(p, e5, T, *br) for T in TS}
        allp = {T: (A * q[T]).mean() for T in TS}
        onA = {T: q[T][A].mean() for T in TS}
        Bm = np.sqrt(br[0] * br[1])
        fm = (p / (p + Bm))[A]
        R[name] = dict(A=A.mean(), etapos=(eta[A] > 0).mean(), pt=pt, pb=pb, w0=100 * w0.mean(), wl=100 * wl.mean(), all=allp, onA=onA,
                       fmed=np.median(fm), flo=(fm < 0.01).mean(), fhi=(fm > 0.5).mean(), x=x2, c=c, p=p, A_=A, eta=eta, br=br)
        r = R[name]
        print(f"\n{name}: c {cr[0]:.0e}..{cr[1]:.3e} M, B {br[0]:.0e}..{br[1]:.3e} M")
        print(f"  gates pass {100 * r['A']:.1f}%; eta*>0 in {100 * r['etapos']:.1f}% of A; P(PVED hand | A): lid layer {pt:.3f}, bulk {pb:.3f}; "
              f"whole-ocean PVED-set, no lid {r['w0']:.2f}%, lid {r['wl']:.2f}%")
        print("  P(ee >= 1/2/5/10%): all " + " / ".join(f"{allp[T]:.3f}" for T in TS) + " | over A " + " / ".join(f"{onA[T]:.3f}" for T in TS))
        print(f"  f at B = {Bm:.2e}, over A: median {r['fmed']:.3f}; < 0.01: {r['flo']:.3f}; > 0.5: {r['fhi']:.3f}")

    r1 = R["S1 primary"]
    print("\nS1 sensitivity, all draws, T = 1 / 2 / 5 / 10 %")
    for tau in TAUS:
        e = r1["eta"] * PS.surv(krA, tau)
        print(f"  tau {tau:g} yr: " + " / ".join(f"{(r1['A_'] * SH.pdet_B(r1['p'], e, T, *B1)).mean():.3f}" for T in TS))
    e5 = r1["eta"] * PS.surv(krA, 1e5)
    fixed = {}
    for B in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        fixed[B] = {T: (r1["A_"] * SH.pdet_fixed(r1["p"], e5, T, B)).mean() for T in TS}
        print(f"  B fixed {B:.0e} M: " + " / ".join(f"{fixed[B][T]:.3f}" for T in TS))
    wide = {T: (r1["A_"] * SH.pdet_B(r1["p"], e5, T, 1e-9, 1e-6)).mean() for T in TS}
    print("  B log-uniform 1e-9..1e-6 M: " + " / ".join(f"{wide[T]:.3f}" for T in TS))
    tau_ee = {tau: (r1["A_"] * SH.pdet_B(r1["p"], r1["eta"] * PS.surv(krA, tau), 0.05, *B1)).mean() for tau in TAUS}

    print("\nSummary vs old (S0), P(plume ee >= 5%), all draws, tau 1e5 yr, plume = lid layer:")
    for name in SCEN:
        print(f"  {name:18} {R[name]['all'][0.05]:.3f}  (1%: {R[name]['all'][0.01]:.3f})")

    print("\nPredictions:")
    s0, s1, s2, s3, s4 = (R[n] for n in SCEN)
    check("P1", abs(100 * s0["A"] - 78.7) <= 0.2 and abs(s0["all"][0.01] - 0.522) <= 0.003 and abs(s0["all"][0.05] - 0.449) <= 0.003
          and abs(s0["pt"] - 0.928) <= 0.003 and abs(s0["w0"] - 26.36) <= 0.05,
          f"gates {100 * s0['A']:.1f}%; ee 1/5% {s0['all'][0.01]:.3f}/{s0['all'][0.05]:.3f}; PVED|A {s0['pt']:.3f}; no-lid {s0['w0']:.2f}%")
    check("P2", 0.95 <= s1["A"] <= 0.99 and s1["A"] > s0["A"], f"S1 gates {100 * s1['A']:.1f}%")
    rng3 = {0.01: (0.78, 0.88), 0.02: (0.70, 0.79), 0.05: (0.54, 0.63), 0.10: (0.40, 0.48)}
    check("P3", all(lo <= s1["all"][T] <= hi and s1["all"][T] > s0["all"][T] for T, (lo, hi) in rng3.items()),
          "S1 T 1/2/5/10%: " + " / ".join(f"{s1['all'][T]:.3f}" for T in TS))
    check("P4", 0.96 <= s1["pt"] <= 0.99 and 0.98 <= s1["pb"] <= 1.0 and s1["pt"] > s0["pt"] and s1["pb"] > s0["pb"],
          f"S1 PVED|A layer {s1['pt']:.3f}, bulk {s1['pb']:.3f}")
    check("P5", 10 <= s1["w0"] <= 28 and s1["wl"] == 0, f"S1 no-lid PVED-set {s1['w0']:.2f}%, lid {s1['wl']:.2f}%")
    check("P6", 0.62 <= s2["A"] <= 0.71 and 0.21 <= s2["all"][0.01] <= 0.31 and 0.11 <= s2["all"][0.05] <= 0.20
          and s2["all"][0.05] < s0["all"][0.05] and 0.83 <= s2["pt"] <= 0.90,
          f"S2 gates {100 * s2['A']:.1f}%; ee 1/5% {s2['all'][0.01]:.3f}/{s2['all'][0.05]:.3f}; PVED|A {s2['pt']:.3f}")
    check("P7", 0.98 <= s3["A"] <= 1.0 and 0.92 <= s3["all"][0.01] <= 0.98 and 0.83 <= s3["all"][0.05] <= 0.91,
          f"S3 gates {100 * s3['A']:.1f}%; ee 1/5% {s3['all'][0.01]:.3f}/{s3['all'][0.05]:.3f}")
    check("P8", 0.52 <= s4["all"][0.05] <= 0.61 and 0.97 <= s4["A"] <= 1.0, f"S4 ee 5% {s4['all'][0.05]:.3f}; gates {100 * s4['A']:.1f}%")
    e5v = [R[n]["all"][0.05] for n in SCEN]
    check("P9", s2["all"][0.05] < s0["all"][0.05] < min(s1["all"][0.05], s4["all"][0.05]) and max(s1["all"][0.05], s4["all"][0.05]) < s3["all"][0.05]
          and s2["w0"] < s1["w0"], "ee 5% by S0..S4: " + " / ".join(f"{v:.3f}" for v in e5v) + f"; no-lid PVED S2 {s2['w0']:.2f} < S1 {s1['w0']:.2f}")
    check("P10", 0.45 <= tau_ee[1e6] <= 0.54 and 0.13 <= tau_ee[1e7] <= 0.22, f"S1 5% at tau 1e6 / 1e7: {tau_ee[1e6]:.3f} / {tau_ee[1e7]:.3f}")
    v = [fixed[B][0.05] for B in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5)]
    check("P11", 0.88 <= fixed[1e-9][0.05] <= 0.94 and 0.34 <= fixed[1e-7][0.05] <= 0.42 and fixed[1e-5][0.05] < 0.02 and all(np.diff(v) < 0),
          "S1 5% by fixed B 1e-9..1e-5: " + " / ".join(f"{u:.3f}" for u in v))
    check("P12", 0.07 <= s1["fmed"] <= 0.12 and 0.05 <= s1["flo"] <= 0.09 and 0.05 <= s1["fhi"] <= 0.10,
          f"S1 f median {s1['fmed']:.3f}; < 0.01: {s1['flo']:.3f}; > 0.5: {s1['fhi']:.3f}")
    check("P13", abs(HCN_X * WATER_M - 0.061) <= 0.001 and HCN_X * WATER_M / C_OLD[1] > 100, f"{HCN_X * WATER_M:.4f} M")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held. {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
