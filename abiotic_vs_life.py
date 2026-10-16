"""If an abiotic amplifier pushed amino-acid excesses near 100%, how would a plume measurement tell it from life?
Run on the model's own Enceladus draws, reusing plume_signal.py (imported, never edited).  Written 2026-10-16.

QUESTION.  A. Davila (NASA Ames, email 2026-10-16): an abiotic L-excess in a plume would be read as abiotic unless it is large
and in more than 2-3 amino acids.  Suppose a Frank-type amplifier does make large excesses.  Which measurable pattern separates
it from life, and which of my first guesses does the model NOT support?

WHAT THE MODEL SAYS (formulas imported or restated from plume_signal.py / persist.py).  An amplifier network sets one hand in
the amino acids it makes with steady purity eta*(r), r = k_rac/(k2 c) (persist.py analytic pure branch; eta* = 0 for r >= 0.087).
A plume amino acid i shows  ee_i = f_i x eta*_i x s_i,  s_i = 1/(1 + 2 k_rac,i tau_res)  (f_i = network share of amino acid i;
tau_res = residence time before sampling).  k_rac at 273 K (sweep2.K_RAC): Ala 9.13e-15 /s, Asp 2.73e-14 /s.  Gates and draws
exactly as plume_signal.py (Enceladus, seed 0, N = 200,000, lid in every draw, plume = lid layer, tau_res = 1e5 yr unless stated,
f log-uniform on 1e-4..1).  The model has NO notion of which amino acids life uses, of abundance patterns, or of which amino
acids the network makes; those come only from the literature below or are stated as assumptions.

SOURCES (DOIs resolved through api.crossref.org on 2026-10-16, no email in the URL; abstracts read via Crossref/Europe PMC, the
McKay page on journals.plos.org; no files downloaded).
 - McKay 2004, PLoS Biol 2:e302, doi 10.1371/journal.pbio.0020302 (web page): "Lego principle": life uses a limited menu of
   molecules, so related molecules appear in very different amounts, while abiotic processes give a smooth spread; Earth proteins
   use only L amino acids; an excess of one hand is a practical first test.  The page I read does not discuss non-protein amino
   acids being racemic; so that claim is NOT sourced to McKay here.
 - Dorn, Nealson, Adami 2011, J Mol Evol 72:283, doi 10.1007/s00239-011-9429-4 (abstract): abiotic amino-acid and carboxylic-acid
   patterns are dominated by small, easily formed, low-formation-energy molecules; life's relative concentrations reflect adaptive
   utility and "belie formation kinetics or equilibrium thermodynamics".  (Their companion: Dorn & Adami 2011, Astrobiology 11:959,
   doi 10.1089/ast.2010.0556, not used further.)
 - Creamer, Mora, Willis 2017, Anal Chem 89:1329, doi 10.1021/acs.analchem.6b04338 (abstract): amino-acid distributions can show
   patterns unique to life; two capillary-electrophoresis methods resolve 17 amino acids (D/L Ala, Asp, Glu, His, Leu, Ser, Val;
   Gly, beta-Ala, GABA), detection limits 5 nM (neutral) and 500 nM (acidic).  No ee limit is given in the abstract.
 - Glavin & Dworkin 2009, PNAS 106:5487, doi 10.1073/pnas.0811618106 (abstract): L-isovaline ee 18.5 +/- 2.6% in Murchison and
   15.2 +/- 4.0% in Orgueil, none in the least altered CR meteorites; isovaline and other alpha-dialkyl amino acids are NOT
   protein amino acids, so an abiotic process can give a large L-excess in a non-protein amino acid, tied to aqueous alteration.
 - Glavin et al. 2020, Chem Rev 120:4660, doi 10.1021/acs.chemrev.9b00474 (abstract): meteorites carry amino acids and aldonic
   acids with large ee of the terrestrial hand from non-biological processes, so ee alone is not a definitive biosignature; they
   propose criteria for life-detection missions (the criteria themselves are not in the abstract and are NOT used here).
 - Pizzarello & Cronin 2000, Geochim Cosmochim Acta 64:329, doi 10.1016/s0016-7037(99)00280-x (title only checked): non-racemic
   amino acids in Murray and Murchison.  Not relied on for any number.
 - Steel, Davila, McKay 2017, Astrobiology 17:862, doi 10.1089/ast.2017.1673 (abstract): abiotic production of Gly, Ala,
   alpha-aminoisobutyric acid (AIB) and Glu above 0.01 uM; racemization is fast against production, so no ee in the ambient ocean,
   and for their methanogen biology the plume ee is "less than 10%".  So Steel's LIFE expectation is a small plume ee.
 - Davila & McKay 2014, Astrobiology 14:534, doi 10.1089/ast.2014.1150: abstract is about chance and necessity in biochemistry;
   it says nothing on ee or abundance patterns, so it is not cited for any discriminator here.
 - Mora et al. 2022, Astrobiology 22:914, doi 10.1089/ast.2021.0091 (abstract): CE-MS limits 0.05-1 uM, no ee limit.

ASSUMPTIONS (stated, not measured).  (a) Everything in plume_signal.py.  (b) For a question of the amplifier's own purity the
draw's gates are those for Ala; for Asp the gate k2 c > 2 k_rac(Asp) replaces k2 c > 2 k_rac(Ala) (the other two gates are unchanged).
(c) Both amino acids come from one network with one eta*(r) curve, evaluated at their own r.  (d) "Independent shares" means each
amino acid has its own f_i, drawn log-uniform on 1e-4..1 and independent; "shared share" means one f for all.  (e) The same
residence-time racemization acts on life's amino acids once they are dissolved (Steel et al.'s premise), so s_i does not
separate life from an amplifier.  (f) k_rac is known in the model only for Ala and Asp, so any ordering test with more than two
amino acids needs k_rac values I did not verify (Bada & Schroeder 1975, not read).

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python abiotic_vs_life.py > abiotic_vs_life.log

PREDICTIONS (written 2026-10-16, before any run of this script).  Disclosure: before writing these I read plume_signal.py and
plume_signal.log (A = 78.7%, eta* > 0 in 93.3% of A, median eta* 0.994, detection chances 0.457/0.294/0.223 for T 1/5/10% at
tau 1e5 yr, survival values); I did not run this script or any probe of it.  Closed forms below were done by hand from those
logged numbers and persist.py's logged shares.  "A" = draws that pass the Ala gates; all shares are over A unless stated.
 (P1) Regression with plume_signal.log: A share 0.786-0.788; eta*_Ala > 0 in 0.931-0.935 of A; median eta*_Ala over A 0.992-0.996.
 (P2) The amplifier's own source purity is life-like: eta*_Ala >= 0.99 in 0.60-0.78 of A draws (so "near 100%" at the source does
      NOT separate an amplifier from life).
 (P3) Rate-gated purity (an amplifier-only feature): among A draws where Ala amplifies, the share where Asp does NOT (eta*_Asp = 0
      because k2 c < 11.5 k_rac(Asp), or the Asp gate fails) is 0.04-0.12.
 (P4) Among A draws where Ala amplifies, the share with eta*_Asp / eta*_Ala < 0.95 (zero counts) is 0.10-0.25; among draws where
      both amplify the median eta*_Asp / eta*_Ala is 0.997-1.000.
 (P5) Survival ratio R = s_Asp / s_Ala = (1 + 2 kA tau)/(1 + 2 kD tau): 0.902 +/- 0.005 at tau 1e5 yr, 0.579 +/- 0.005 at 1e6,
      0.370 +/- 0.005 at 1e7; always between kA/kD = 0.334 and 1; inverting R back to tau reproduces tau to 1e-6 (relative).
 (P6) Magnitude (tau 1e5 yr, f log-uniform, over A): P(plume ee_Ala >= 50%) is 0.055-0.072; P(>= 90%) is 0.002-0.006; of the
      draws with ee >= 5%, the share that also reach 50% is 0.18-0.26.  Near-100% plume ee needs the network to make nearly
      all of that amino acid, so it is rare even if the source is pure.
 (P7) Steel anchor: P(plume ee_Ala >= 10%) over A is 0.19-0.24 (plume_signal's 0.223); that is the share of draws where an
      amplifier gives an excess ABOVE Steel et al.'s <10% expectation for biology (model-dependent; abstract only).
 (P8) Many amino acids (tau 1e5 yr, T = 5%, all with ee_Ala-type survival and the draw's eta*): P(all of n amino acids >= 5%):
      shared share, n = 4: 0.28-0.31; independent shares: n = 3 0.02-0.045, n = 4 0.006-0.013, n = 6 0.0003-0.0020.
 (P9) Ordering by racemization: if the n amino acids' ee were in random order, the chance that they follow the k_rac order is
      1/n!: 0.167 (n = 3), 0.0417 (n = 4), 0.0083 (n = 5) exactly; the smallest n with chance below 0.05 is 4.
 (P10) Chemistry check (RDKit): AIB and glycine have no stereocentre (zero), alanine, isovaline and norvaline have exactly one.
       So AIB cannot carry an ee at all and cannot serve as the "racemic non-protein amino acid" control; isovaline and
       norvaline can.

RESULTS (abiotic_vs_life.log, 2026-10-16): 10 of 10 held.  1 s.  The predictions were committed in c842151 and are unchanged; no
smoke run or probe was made before the run.  P1, P5, P7 are near-regressions of numbers already in plume_signal.log or closed forms,
so they test the reuse, not new physics; P2-P4, P6, P8 are the new model results.
 (P1) HELD  A 0.7872; eta*_Ala > 0 in 0.9335 of A; median 0.9940.
 (P2) HELD  eta*_Ala >= 0.99 in 0.659 of A (>= 0.9: 0.857).  The amplifier's source purity is life-like.
 (P3) HELD  Asp does not amplify in 0.066 of the A draws where Ala does (k2 c between 11.5 and about 34 k_rac(Ala)).
 (P4) HELD  eta*_Asp/eta*_Ala < 0.95 in 0.179 of those draws (10th percentile of the ratio 0.765); median 0.9997 where both amplify.
 (P5) HELD  R = s_Asp/s_Ala 0.9022 (1e5 yr), 0.5788 (1e6), 0.3709 (1e7); limits 1 and 0.3344; inversion exact to 1e-13.
 (P6) HELD  tau 1e5 yr: P(ee >= 50%) 0.061, P(>= 90%) 0.0035; of the draws with ee >= 5% (0.294), 20.7% reach 50%.
 (P7) HELD  P(ee >= 10%) 0.223 (above Steel et al.'s <10% expectation for biology, in their model).
 (P8) HELD  All of n amino acids >= 5%: shared share 0.294 for any n; independent shares 0.0293 (n = 3), 0.00925 (n = 4), 0.00093 (n = 6),
            so a network share that differs by amino acid makes Davila's ">3 amino acids" 32 times rarer at n = 4.
 (P9) HELD  1/n! = 0.167, 0.0417, 0.0083 (n = 3, 4, 5); n = 4 is the smallest n below 0.05.
 (P10) HELD AIB and glycine: 0 stereocentres; alanine, isovaline, norvaline: 1.  AIB cannot show ee.
 Reading.  In this model magnitude does not separate an amplifier from life (source purity is >= 0.99 in two thirds of gate-passing
 draws), and the racemization ordering is shared by both once the amino acids are dissolved.  What the model adds is (i) a
 rate-gated purity, a step: in about 7% of the draws where Ala is pure, Asp is racemic at the source, and in 18% Asp is at least
 5% less pure; (ii) a count: an excess in more than 3 amino acids is rare unless the network's share is common to them
 (0.29 shared vs 0.009 independent, n = 4); (iii) a size: only 6.1% of gate-passing draws give a plume ee of 50% or more, because f must be near 1.
 The Lego/abundance pattern and the racemic-non-protein expectation for life are NOT model results; they rest on McKay 2004,
 Dorn et al. 2011 and Creamer et al. 2017 (abstracts/page only), and meteorites show abiotic ee in a non-protein amino acid
 (isovaline, Glavin and Dworkin 2009), so a large ee in non-protein amino acids points to an abiotic amplifier, not to life.
"""
import os
import sys
import time

import numpy as np

import plume_signal as PS   # installs seedfix's rng correction on import; main() is guarded
import uncertainty as U
import sweep2 as S

YR = 3.156e7
CHECK = []


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


def inv_tau(R, kA, kD):
    """Residence time (s) that gives s_Asp/s_Ala = R: R(1 + 2 kD tau) = 1 + 2 kA tau."""
    return (1 - R) / (2 * (R * kD - kA))


def main(seed=0):
    t0 = time.time()
    _, W = S.chains(seed)
    x, win, a, salt = W["Enceladus"]
    L = a[1]
    c, k2, Dh, t = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time"))
    k = k2 * c
    krA, krD = S.K_RAC["Ala"], S.K_RAC["Asp"]
    mix = (L / U.A) ** 2 / Dh <= t
    A = (k * t >= 10) & (k > 2 * krA) & mix
    AD = (k * t >= 10) & (k > 2 * krD) & mix
    tab = PS.eta_table()
    eA = np.where(k > 0, tab(krA / k), 0.0)
    eD = np.where(k > 0, tab(krD / k), 0.0) * AD
    print(f"abiotic_vs_life.py, N = {U.N}, seed 0, Enceladus; imports plume_signal.py.")

    print("\n(1) Source purity of the amplifier (the part a 'near 100%' excess would come from)")
    mA = A & (eA > 0)
    print(f"  A = {A.mean():.4f}; eta*_Ala > 0 in {(eA[A] > 0).mean():.4f} of A; median eta*_Ala over A {np.median(eA[A]):.4f}; "
          f"eta* >= 0.9: {(eA[A] >= 0.9).mean():.3f}, >= 0.99: {(eA[A] >= 0.99).mean():.3f}, >= 0.5: {(eA[A] >= 0.5).mean():.3f}")

    print("\n(2) Rate-gated purity: the same network applied to Ala (slow racemizer) and Asp (3x faster)")
    noAsp = (eD[mA] == 0).mean()
    ratio = np.where(eA > 0, eD / np.where(eA > 0, eA, 1), np.nan)
    both = mA & (eD > 0)
    low = (ratio[mA] < 0.95).mean()
    print(f"  Among A draws where Ala amplifies: Asp does not {noAsp:.4f}; eta*_Asp/eta*_Ala < 0.95 in {low:.4f}; "
          f"median ratio where both amplify {np.median(ratio[both]):.4f}; 10th percentile of ratio over Ala-amplifying draws "
          f"{np.percentile(ratio[mA], 10):.3f}")

    print("\n(3) Survival ratio R = s_Asp/s_Ala (same for life and amplifier once dissolved) and its inversion to tau_res")
    Rs = {}
    for tau in PS.TAUS:
        Rs[tau] = PS.surv(krD, tau) / PS.surv(krA, tau)
        back = inv_tau(Rs[tau], krA, krD) / YR / tau - 1
        print(f"  tau {tau:g} yr: s_Ala {PS.surv(krA, tau):.4f} s_Asp {PS.surv(krD, tau):.4f} R {Rs[tau]:.4f}; inversion error {back:.1e}")
    print(f"  Limits: R -> 1 as tau -> 0, R -> kA/kD = {krA / krD:.4f} as tau -> infinity.")
    inv_err = max(abs(inv_tau(Rs[tau], krA, krD) / YR / tau - 1) for tau in PS.TAUS)

    print("\n(4) Magnitude of the plume excess (tau_res 1e5 yr, f log-uniform 1e-4..1, over A)")
    e = eA * PS.surv(krA, 1e5)
    pd = {T: PS.pdet(e, T)[A].mean() for T in (0.05, 0.10, 0.50, 0.90)}
    for T, v in pd.items():
        print(f"  P(ee_Ala >= {100 * T:g}%) = {v:.4f}")
    cond50 = pd[0.50] / pd[0.05]
    print(f"  Of draws with ee >= 5%, the share that reach 50%: {cond50:.3f}")
    print("  Steel et al. 2017 (abstract): biology in their model gives a plume ee below 10%; an amplifier gives >= 10% in "
          f"{100 * pd[0.10]:.1f}% of A draws here.")

    print("\n(5) How many amino acids carry a detectable excess (T = 5%, tau 1e5 yr)")
    p5 = PS.pdet(e, 0.05)
    shared = {n: p5[A].mean() for n in (1, 2, 3, 4, 5, 6)}
    indep = {n: (p5[A] ** n).mean() for n in (1, 2, 3, 4, 5, 6)}
    for n in (1, 2, 3, 4, 5, 6):
        print(f"  n = {n}: all n >= 5%: shared share {shared[n]:.4f}, independent shares {indep[n]:.5f}")
    print("  Davila's 'more than 2-3 amino acids' therefore costs a factor " +
          f"{shared[4] / indep[4]:.0f} at n = 4 if the network's share differs by amino acid, none if it is common.")

    print("\n(6) Ordering by racemization rate: chance the ee of n amino acids follow the k_rac order at random")
    import math
    po = {n: 1 / math.factorial(n) for n in (2, 3, 4, 5)}
    print("  " + "; ".join(f"n = {n}: {v:.4f}" for n, v in po.items()))
    print("  Only Ala and Asp have k_rac in the model; with two, the ordering is a coin flip (0.5).")

    print("\n(7) Chirality of candidate control amino acids (RDKit stereocentres)")
    from rdkit import Chem
    smi = {"Gly": "NCC(=O)O", "AIB": "CC(C)(N)C(=O)O", "Ala": "CC(N)C(=O)O", "isovaline": "CCC(C)(N)C(=O)O",
           "norvaline": "CCCC(N)C(=O)O"}
    nc = {n: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True)) for n, s in smi.items()}
    print("  " + "; ".join(f"{n}: {v}" for n, v in nc.items()))

    print("\n(8) Ledger: what separates an amplifier from life, and what supports it")
    for row in (
        "magnitude of source ee near 100%       | NOT a discriminator: model eta* >= 0.99 in the majority of A draws (P2)       | model",
        "plume ee above 10%                     | outside Steel's life expectation (<10%); amplifier gives it in ~22% of A      | model + Steel abstract",
        "ee ordered by k_rac, common share f    | consistency check only: life's dissolved amino acids racemize identically     | model (assumption e)",
        "Asp racemic while Ala is pure         | amplifier-only (threshold k2c > 11.5 k_rac); life keeps proteins ~pure        | model (P3, P4)",
        "ee in non-protein chiral amino acids   | abiotic ee exists (isovaline 18.5%, Murchison); life: expected racemic        | literature (Glavin 2009); life side = logic",
        "AIB as control                         | achiral, cannot show ee; use isovaline or norvaline                           | chemistry (P10)",
        "abundance pattern (Lego, smooth vs peaked) | life peaked on a selected set, abiotic smooth and Gly/Ala heavy          | literature only (McKay 2004, Dorn 2011, Creamer 2017); model silent",
    ):
        print("  " + row)

    print("\nPredictions:")
    check("P1", 0.786 <= A.mean() <= 0.788 and 0.931 <= (eA[A] > 0).mean() <= 0.935 and 0.992 <= np.median(eA[A]) <= 0.996,
          f"A {A.mean():.4f}; eta>0 {(eA[A] > 0).mean():.4f}; median {np.median(eA[A]):.4f}")
    q99 = (eA[A] >= 0.99).mean()
    check("P2", 0.60 <= q99 <= 0.78, f"eta*_Ala >= 0.99 in {q99:.3f} of A")
    check("P3", 0.04 <= noAsp <= 0.12, f"Asp does not amplify in {noAsp:.4f} of Ala-amplifying A draws")
    check("P4", 0.10 <= low <= 0.25 and 0.997 <= np.median(ratio[both]) <= 1.0,
          f"ratio<0.95 in {low:.4f}; median ratio (both) {np.median(ratio[both]):.4f}")
    check("P5", abs(Rs[1e5] - 0.902) <= 0.005 and abs(Rs[1e6] - 0.579) <= 0.005 and abs(Rs[1e7] - 0.370) <= 0.005
          and all(krA / krD <= v <= 1 for v in Rs.values()) and inv_err <= 1e-6,
          f"R {Rs[1e5]:.4f}/{Rs[1e6]:.4f}/{Rs[1e7]:.4f}; inversion error {inv_err:.1e}")
    check("P6", 0.055 <= pd[0.5] <= 0.072 and 0.002 <= pd[0.9] <= 0.006 and 0.18 <= cond50 <= 0.26,
          f"ee>=50%: {pd[0.5]:.4f}; >=90%: {pd[0.9]:.4f}; share of 5% draws that reach 50%: {cond50:.3f}")
    check("P7", 0.19 <= pd[0.10] <= 0.24, f"P(ee >= 10%) over A = {pd[0.10]:.4f}")
    check("P8", 0.28 <= shared[4] <= 0.31 and 0.02 <= indep[3] <= 0.045 and 0.006 <= indep[4] <= 0.013 and
          0.0003 <= indep[6] <= 0.0020,
          f"shared n=4 {shared[4]:.4f}; independent n=3 {indep[3]:.4f}, n=4 {indep[4]:.5f}, n=6 {indep[6]:.6f}")
    n_min = min(n for n, v in po.items() if v < 0.05)
    check("P9", abs(po[3] - 0.167) < 1e-3 and abs(po[4] - 0.0417) < 1e-4 and abs(po[5] - 0.0083) < 1e-4 and n_min == 4,
          f"1/n!: {po[3]:.4f}/{po[4]:.4f}/{po[5]:.4f}; smallest n below 0.05: {n_min}")
    check("P10", nc["AIB"] == 0 and nc["Gly"] == 0 and nc["Ala"] == 1 and nc["isovaline"] == 1 and nc["norvaline"] == 1,
          f"stereocentres {nc}")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held. {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
