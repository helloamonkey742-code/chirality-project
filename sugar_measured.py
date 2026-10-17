"""What do EXPERIMENTS say about amino acids steering sugar chirality (and the reverse)?  Replaces sugar_pair.py's guessed
coupling direction rho = k_as / k_sa (log-uniform 0.01-100, "never measured") by a band built from published experiments.
Author's script, 2026-10-16. Imports sugar_pair.py (cfac) and dna_network.py; edits neither. Question from A. Davila (NASA):
does the model also explain life's D-sugars, not just L-amino acids?

SOURCES (every DOI checked against api.crossref.org on 2026-10-16; abstracts read on PubMed / Europe PMC; numbers marked
"text" were read from the open full text of the paper itself or from a later open paper that reports them, by a web-fetch
summary, and not re-derived from the original tables, so treat them as secondary until checked against the PDF).
 B&C   Breslow R, Cheng Z-L (2010) PNAS 107:5723-5725, doi:10.1073/pnas.1001639107. L-amino acids (not proline) catalyse
       glycolaldehyde + formaldehyde -> glyceraldehyde with D/L up to 60/40; proline is the exception; solution amplified to
       92/8 by solubility. Abstract: yes. Per-amino-acid D% (Ser 50.3, Ala 50.8, Phe 52.2, Val 52.2, Leu 54.4, Glu 60.7,
       Pro 28.9): text. Input amino-acid ee not stated in the abstract; enantiopure L assumed.
 P&W04 Pizzarello S, Weber AL (2004) Science 303:1151, doi:10.1126/science.1093057. L-alanine / L-isovaline catalyse tetrose
       formation. Abstract not on PubMed; "up to about 12% D-threose ee with L-isovaline": text, as restated by Weber and
       Pizzarello 2006 (PMC1568914). Input ee unknown, so 0.12 is a lower bound on the transfer.
 P&W06 Weber AL, Pizzarello S (2006) PNAS 103:12713-12717, doi:10.1073/pnas.0602320103. L-dipeptides (L-Val-L-Val best)
       -> D-erythrose ee > 80% (abstract).
 P&W10 Pizzarello S, Weber AL (2010) Orig Life Evol Biosph 40:3-10, doi:10.1007/s11084-009-9178-1. LL-dipeptide catalysts:
       ribose D-ee up to 44%, lyxose L-ee up to 66%, arabinose L-ee up to 8%, xylose about racemic (abstract).
 H11   Hein JE, Tse E, Blackmond DG (2011) Nat Chem 3:704-706, doi:10.1038/nchem.1108. A small initial imbalance of amino-acid
       enantiomers yields highly enantioenriched RNA precursors from racemic starting materials (abstract). Details (1% ee
       L-proline enriched by a CHCl3 cosolvate, then L-aa binding L-glyceraldehyde and leaving D) are from the open review
       of Blackmond 2019, doi:10.1101/cshperspect.a032540: text, secondary. This route has physical amplification, so it
       is NOT a linear gain and is kept out of the primary rho band (only a sensitivity row).
 H&B12 Hein JE, Blackmond DG (2012) Acc Chem Res 45:2045-2054, doi:10.1021/ar200316n. Amplification of a small imbalance of
       sugars by amino acids or of amino acids by sugars: both directions (abstract).
 T&S04 Tamura K, Schimmel P (2004) Science 305:1253, doi:10.1126/science.1099141. A D-RNA minihelix is aminoacylated with a
       clear preference for L-amino acids; the mirror-image RNA reverses it (abstract). About 4-fold L over D (Ala, Leu, Phe):
       text, from Tamura and Schimmel 2006 PNAS 103:13750-13752, doi:10.1073/pnas.0606070103 (PMC1564265); there Ac-L-Ala /
       Ac-D-Ala 3.7:1 on D-ribose and 1:3.6 on L-ribose.
 Y10   Illangasekare M, Turk R, Peterson GC, Lladser M, Yarus M (2010) RNA 16:2370-2383, doi:10.1261/rna.2385310. D-ribose RNA
       sites for L-histidine are 4-6 times more abundant than for D-histidine; at equilibrium with racemic His at least
       7 of 8 complexes would hold L-His (abstract). Best estimate 86% L sites: text.
 W17   Wagner AJ, Zubarev DY, Aspuru-Guzik A, Blackmond DG (2017) ACS Cent Sci 3:322-328, doi:10.1021/acscentsci.7b00085. D-pentose
       sugars drive the formation of natural L-amino-acid precursors (aminoamides of Ala, Phe, Trp), ee 11-83% across the
       D-sugar reactions (text); D-ribose and D-lyxose push opposite ways, but a mixture of natural D-sugars still favours L
       (abstract). This is the sugar -> amino acid direction by pure chemistry, with no RNA template.
 Cordova et al. (2005) Chem Eur J 11:4772-4784, doi:10.1002/chem.200500139: amino acids catalyse asymmetric sugar formation
       (hexoses, ee > 99% for proline-catalysed allose in organic solvent). Qualitative support only; not used numerically.

HOW THE EXPERIMENTS MAP ONTO THE MODEL (dna_network.py: d alpha = (lam alpha - alpha^3 + k_sa beta + g_a) dt + noise,
d beta = (lam beta - beta^3 + k_as alpha + g_s) dt + noise; alpha > 0 is L-amino acid, beta > 0 is D-sugar).
 SIGN. The coupling terms are +k beta and +k alpha with k > 0: a positive alpha (L) pushes beta positive (D). So the model as
 written ASSUMES the natural direction, L-aa <-> D-sugar. sugar_pair.py never states this; it only varies the size ratio rho. The
 experiments are what justify it: aa -> sugar (k_as): L-aa gives D-glyceraldehyde (B&C), D-tetrose (P&W04, P&W06), D-ribose
 (P&W10, with LL-dipeptides); sugar -> aa (k_sa): D-RNA picks L-aa (T&S04, Y10), D-sugars give L-aa precursors (W17). BOTH measured
 directions are therefore positive for the sugars on the RNA path, which is the case the model's sqrt(k_as k_sa) assumes. If both
 k were negative the growing mode would be (alpha, -beta): the pair that forms is L-aa with L-sugar (or D with D) and the natural
 pair cannot form (checked below with dna_network.simulate_asymmetric). If their signs differed, the linear part has complex
 eigenvalues lam +- i sqrt(|k_as k_sa|): no stable matched pair (not simulated; the function takes sqrt of a negative number).
 Honest caveat: the sign is not uniformly natural. Proline gives L-glyceraldehyde (B&C), and with LL-dipeptides lyxose and
 arabinose come out L-rich (P&W10) and D-ribose / D-lyxose push W17's amino-acid precursor opposite ways. The model has one k
 per direction, so it describes the net sign along the glyceraldehyde -> ribose -> RNA path, not every sugar.
 MAGNITUDE. An experiment gives an output ee for an input ee (gain G = ee_out / ee_in), not a rate constant. In the linear
 regime below the bifurcation (lam < 0), a clamped alpha gives beta = k_as alpha / |lam|, so G_as = k_as / |lam_s| and G_sa = k_sa /
 |lam_a|; assuming the two sectors have the same damping |lam|, rho = k_as / k_sa = G_as / G_sa. This is crude: the gains come
 from different chemistry, at saturation, under different conditions; the only claim is order of magnitude. To keep it honest the
 primary band uses the whole spread of the data and the sensitivity rows move it by a factor 3 either way.
 PRIMARY BAND (no free parameter, no fit): G_as log-uniform 0.044-0.44 (Phe/Val -> D-glyceraldehyde 4.4% ee; LL-dipeptide ->
 D-ribose 44% ee; B&C's best, L-Glu, 21.4%, sits inside) and G_sa log-uniform 0.11-0.83 (W17's 11-83%; T&S04 57% and Y10 72% sit
 inside); rho = G_as / G_sa drawn independently. D-erythrose (80%) is left out: it is not an RNA sugar. Everything else of
 sugar_pair.py is kept: p(r > 0) = 0.5, |r| log-uniform 0.03-3, conditioned on g_a > 0, natural pair favoured iff
 c = (sqrt(rho) + r) / sqrt(rho + 1) > 0.
 VARIANTS: (b) RNA-templating only (G_as 0.044-0.21 monomeric amino acids, G_sa 0.57-0.72); (c) Hein's physical amplification
 allowed (G_as up to 1.0); (d) rho / 3; (e) rho x 3; p(r > 0) = 0, 0.25, 0.75; |r| 0.1-1 and 0.01-10.
 NOT DONE: the oceans (sugar_pair.py section 4); only the share of draws that favour the natural pair is recomputed.

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python sugar_measured.py > sugar_measured.log

PREDICTIONS (written 2026-10-16, before any run of this script). The numbers for P3-P6 come from a one-off scratch Monte Carlo
(2,000,000 draws, no code from this file) done while choosing the band, so they are expectations from that scratch run, not blind;
P1, P2, P7 and P8 are arithmetic or committed-log facts, and the simulation outputs of P7 were not looked at.
 (P1) Regression. sugar_pair's primary closed form exact(0.5, (0.03, 3), (0.01, 100)) = 0.8636 +- 0.0005 (its logged value), and
      cfac equals dna_network's g_eff / g_a to 1e-9 in every spot check below.
 (P2) Arithmetic on the facts table. B&C: 6 of 6 non-proline catalysts D-sided (D% > 50); best L-Glu 21.4% ee (+- 0.1), L-Pro
      -42.2% ee (+- 0.1; the L side); T&S 4:1 -> 60.0% ee, 3.7:1 -> 57.4% ee (+- 0.1); Y10 86% L sites -> 72% ee, 7 of 8 -> 75% ee.
 (P3) The measured rho band: median 0.44-0.48, 5th percentile 0.09-0.12, 95th percentile 1.9-2.1, P(rho < 1) 0.78-0.84. So
      sugar -> amino acid transmission is typically a bit stronger than the reverse (rho < 1), the opposite of "no preference".
 (P4) Primary natural-pair share (g_a > 0, c > 0): 0.828-0.848 by Monte Carlo, the grid integral agrees within 0.004, and it lies
      below sugar_pair's 0.8636 by 0.016-0.036. The measured band does not raise the share; it lowers it slightly.
 (P5) Variants of that share: (b) RNA-templating only 0.765-0.790; (c) Hein included 0.850-0.870; (d) rho / 3 0.770-0.790;
      (e) rho x 3 0.890-0.910; p(r > 0) = 0: 0.660-0.690, 0.25: 0.750-0.770, 0.75: 0.910-0.930; |r| 0.1-1: 0.895-0.915, 0.01-10:
      0.795-0.815; p(r > 0) = 1: 1.000. No variant at p(r > 0) = 0.5 falls below 0.75 or above 0.92.
 (P6) Among r < 0 draws (every sugar PVED opposing), the share of draws still favouring the natural pair (rescued) is 0.665-0.690.
 (P7) Coupling sign, dna_network.simulate_asymmetric (kappa_bar 0.3, g_a = 1e-3, 4000 runs): five positive-coupling cases at the
      measured rho (r, rho) = (0.257, 0.46), (-0.3, 0.46), (-0.678, 0.46), (-1, 0.46), (-0.3, 0.15) have matched >= 0.99 and natural
      rate within 0.04 of Phi(Delta(g_eff)); the c ~ 0 case (-0.678, 0.46) gives 0.50 +- 0.06 (its c is -0.002, not exactly 0). With both
      couplings negative (-0.3, -0.3; r = 0.3) matched <= 0.01: the natural pair never forms.
 (P8) Sign tally of the abstract-verified results (facts table, rows flagged A): aa -> sugar 3 natural / 3 against (B&C non-Pro;
      B&C Pro against; P&W10 ribose natural, lyxose and arabinose against; P&W06 erythrose natural); sugar -> aa 3 natural / 0 against
      (T&S04, Y10, W17).
 (P9) Cost: the whole script finishes in under 5 minutes at 2 threads under nice.

DISCLOSURE. Before this commit I read sugar_pair.py, sugar_pair.log and dna_network.py in full, and the abstracts / open texts named
above. Apart from the scratch Monte Carlo mentioned at P3-P6, I did not run this script or any part of it.
RESULTS (sugar_measured.log, 2026-10-16): 9 of 9 held (P1-P9), whole script 16 s.
 Measured rho = k_as / k_sa band: median 0.46, 5-95% 0.105-2.02, P(rho < 1) = 0.79 (sugar -> amino acid is typically the stronger
 direction). Natural-pair share (g_a > 0): primary 0.838 (grid 0.839) against sugar_pair's 0.864, a change of -0.025. Variants:
 RNA-templating only 0.778; Hein amplification allowed 0.861; rho / 3 0.779; rho x 3 0.898; p(r>0) = 0 / 0.25 / 0.75 / 1:
 0.677 / 0.757 / 0.919 / 1.000; |r| 0.1-1 0.905; |r| 0.01-10 0.805. Rescue of an opposing sugar PVED: 0.677 of r < 0 draws.
 Sign: positive couplings at the measured rho reproduce Phi(Delta(g_eff)) (natural rate within 0.012; c = 0 case 0.501); with both
 couplings negative matched = 0.000, so the natural pair cannot form. Sign tally (abstract-verified): aa -> sugar 3 natural /
 3 against (proline, lyxose, arabinose), sugar -> aa 3 natural / 0 against.
 Reading: the experiments fix the SIGN of both couplings as natural for the sugars on the RNA path (glyceraldehyde, ribose), which
 sugar_pair.py had assumed without saying so, and put rho mostly below 1, not above it. That lowers the share of draws favouring
 the natural pair only slightly (86.4% -> 83.8%; 78-90% over every variant with a random sugar-PVED sign), and it does not touch the
 limit that the pairing is explained, not which pair: g_a > 0 is still taken from the amino-acid PVED. The gain-ratio mapping
 (rho = G_as / G_sa, equal damping) is crude, the sign is not uniformly natural (proline, lyxose, arabinose), and the numbers marked
 "text" in the facts table are from open full texts read through a web-fetch summary and not checked against the PDFs.
 CORRECTION (2026-10-16, after review; the lines above are left as written): the sugar -> aa tally is 3 natural / 1 against,
 not 3 / 0.  Wagner et al. 2017 (W17) gives the natural L precursor only with D-lyxose and mixtures of D-sugars; D-ribose and
 D-deoxyribose, the RNA sugars, give the unnatural D side.  So "fix the SIGN of both couplings as natural for ... ribose" does
 not hold for sugar -> aa.  The P8 line and its count are left as logged.  sugar_bands.py recomputes the tally and the band
 choices (natural-pair share 0.743-0.847, rescue 0.485-0.695); this Reading is superseded by sugar_bands.py.
"""
import time

import numpy as np
from scipy.stats import norm

from sugar_pair import cfac, exact   # importing sugar_pair re-installs seedfix's rng wrapper; harmless here
import dna_network as D              # noqa: E402
from sim import delta                # noqa: E402

T0 = time.time()
RES = {}
NDRAW = 2_000_000


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


# ---------------------------------------------------------------- 1. facts table
# (direction, catalyst/template, product, ee_in %, signed ee_out %, + = natural, flag A = sign in an abstract, T = from text, DOI)
BC = {"L-Ser": 50.3, "L-Ala": 50.8, "L-Phe": 52.2, "L-Val": 52.2, "L-Leu": 54.4, "L-Glu": 60.7, "L-Pro": 28.9}   # D% of glyceraldehyde
ee_bc = {k: 2 * v - 100 for k, v in BC.items()}   # (D - L)/(D + L) in percent; + means D
DOI = dict(BC="10.1073/pnas.1001639107", PW04="10.1126/science.1093057", PW06="10.1073/pnas.0602320103",
           PW10="10.1007/s11084-009-9178-1", H11="10.1038/nchem.1108", HB12="10.1021/ar200316n",
           TS04="10.1126/science.1099141", TS06="10.1073/pnas.0606070103", Y10="10.1261/rna.2385310", W17="10.1021/acscentsci.7b00085")
FACTS = [  # direction, system, product, ee_in, ee_out (+ natural), flag, doi key
    ("aa->sugar", "L-Glu (best)", "glyceraldehyde", "100*", ee_bc["L-Glu"], "A", "BC"),
    ("aa->sugar", "L-Leu", "glyceraldehyde", "100*", ee_bc["L-Leu"], "T", "BC"),
    ("aa->sugar", "L-Phe / L-Val", "glyceraldehyde", "100*", ee_bc["L-Phe"], "T", "BC"),
    ("aa->sugar", "L-Ala", "glyceraldehyde", "100*", ee_bc["L-Ala"], "T", "BC"),
    ("aa->sugar", "L-Ser", "glyceraldehyde", "100*", ee_bc["L-Ser"], "T", "BC"),
    ("aa->sugar", "L-Pro (exception)", "glyceraldehyde", "100*", ee_bc["L-Pro"], "A", "BC"),
    ("aa->sugar", "L-isovaline / L-Ala", "D-threose", "?", 12.0, "T", "PW04"),
    ("aa->sugar", "L-Val-L-Val", "D-erythrose", "100*", 80.0, "A", "PW06"),
    ("aa->sugar", "LL-dipeptides", "D-ribose", "100*", 44.0, "A", "PW10"),
    ("aa->sugar", "LL-dipeptides", "lyxose (L-rich)", "100*", -66.0, "A", "PW10"),
    ("aa->sugar", "LL-dipeptides", "arabinose (L-rich)", "100*", -8.0, "A", "PW10"),
    ("aa->sugar", "1% ee L-Pro (+ amplification)", "ribo-aminooxazoline", "1", 99.0, "T", "H11"),
    ("sugar->aa", "D-RNA minihelix", "Ala/Leu/Phe-RNA", "100", 60.0, "A", "TS04"),   # 4:1 L:D in T&S06 text; sign in abstract
    ("sugar->aa", "D-RNA minihelix, Ac-Ala", "Ac-L-Ala (3.7:1)", "100", 100 * (3.7 - 1) / (3.7 + 1), "T", "TS06"),
    ("sugar->aa", "D-ribose RNA sites", "L-His complexes (86% L)", "100", 2 * 86 - 100, "A", "Y10"),
    ("sugar->aa", "D-pentoses (text range 11-83)", "L-aminoamides", "100*", 83.0, "A", "W17"),
]


def section1():
    print("=== 1. Experimental facts (ee in %, ee out %, + = natural: D-sugar or L-amino acid); * = assumed enantiopure ===")
    print(f"{'direction':10} {'catalyst / template':32} {'product':26} {'ee_in':>6} {'ee_out':>7} {'flag':>4}  DOI")
    for d, s, p, ei, eo, fl, k in FACTS:
        print(f"{d:10} {s:32} {p:26} {ei:>6} {eo:7.1f} {fl:>4}  {DOI[k]}")
    print("  B&C also: D/L 60.7/39.3 at the best, amplified to 92/8 (84% ee) by solubility; flag A = sign stated in an abstract.")
    non_pro = [v for k, v in BC.items() if k != "L-Pro"]
    t1 = (3.7 - 1) / (3.7 + 1) * 100
    ok = (all(v > 50 for v in non_pro) and abs(ee_bc["L-Glu"] - 21.4) <= 0.1 and abs(ee_bc["L-Pro"] + 42.2) <= 0.1
          and abs(100 * (4 - 1) / 5 - 60.0) < 0.1 and abs(t1 - 57.4) <= 0.1 and 2 * 86 - 100 == 72 and 100 * (7 - 1) / 8 == 75.0)
    check("P2", ok, f"6/6 non-Pro D-sided: {all(v > 50 for v in non_pro)}; Glu {ee_bc['L-Glu']:.1f}, Pro {ee_bc['L-Pro']:.1f}; "
          f"4:1 -> 60.0, 3.7:1 -> {t1:.1f}; 86% L -> 72; 7/8 -> 75")
    # P8: sign tally over rows with the sign stated in an abstract
    aa = [("BC non-Pro", +1), ("BC Pro", -1), ("PW10 ribose", +1), ("PW10 lyxose", -1), ("PW10 arabinose", -1), ("PW06 erythrose", +1)]
    sa = [("TS04", +1), ("Y10", +1), ("W17", +1)]
    na, ba = sum(s > 0 for _, s in aa), sum(s < 0 for _, s in aa)
    ns, bs = sum(s > 0 for _, s in sa), sum(s < 0 for _, s in sa)
    print(f"  Sign tally (abstract-verified): aa->sugar {na} natural / {ba} against; sugar->aa {ns} natural / {bs} against")
    check("P8", (na, ba, ns, bs) == (3, 3, 3, 0), f"aa->sugar {na}/{ba}, sugar->aa {ns}/{bs}")


# ---------------------------------------------------------------- 2. rho band and natural-pair share
LU = lambda rng, a, b, n: 10 ** rng.uniform(np.log10(a), np.log10(b), n)


def band(rng, gA, gB, scale=1.0, n=NDRAW):
    return LU(rng, *gA, n) / LU(rng, *gB, n) * scale


def share_mc(rng, gA, gB, p_plus=0.5, rabs=(0.03, 3.0), scale=1.0, n=NDRAW):
    rho = band(rng, gA, gB, scale, n)
    r = np.where(rng.random(n) < p_plus, 1.0, -1.0) * LU(rng, *rabs, n)
    c = cfac(r, rho)
    neg = r < 0
    return float((c > 0).mean()), float((c[neg] > 0).mean()) if neg.any() else float("nan")


def share_grid(gA, gB, p_plus=0.5, rabs=(0.03, 3.0), scale=1.0, m=2000):
    """p_plus + (1 - p_plus) P(rho > r^2 | r < 0), rho = G_as / G_sa x scale; independent numerical integral."""
    la = np.log10(gA[0]) + (np.arange(m) + .5) / m * np.log10(gA[1] / gA[0])
    lb = np.log10(gB[0]) + (np.arange(m) + .5) / m * np.log10(gB[1] / gB[0])
    diff = np.sort((la[:, None] - lb[None, :]).ravel()) + np.log10(scale)   # log10 rho
    u = np.log10(rabs[0]) + (np.arange(m) + .5) / m * np.log10(rabs[1] / rabs[0])
    pr = 1 - np.searchsorted(diff, 2 * u) / diff.size    # P(log10 rho > 2 log10 |r|)
    return p_plus + (1 - p_plus) * float(pr.mean())


G_AS, G_SA = (0.044, 0.44), (0.11, 0.83)
VARIANTS = [("PRIMARY", dict(gA=G_AS, gB=G_SA)),
            ("(b) RNA templating only", dict(gA=(0.044, 0.21), gB=(0.57, 0.72))),
            ("(c) Hein amplification (G_as to 1.0)", dict(gA=(0.044, 1.0), gB=G_SA)),
            ("(d) rho / 3", dict(gA=G_AS, gB=G_SA, scale=1 / 3)), ("(e) rho x 3", dict(gA=G_AS, gB=G_SA, scale=3.0)),
            ("p(r>0) = 0", dict(gA=G_AS, gB=G_SA, p_plus=0.0)), ("p(r>0) = 0.25", dict(gA=G_AS, gB=G_SA, p_plus=0.25)),
            ("p(r>0) = 0.75", dict(gA=G_AS, gB=G_SA, p_plus=0.75)), ("p(r>0) = 1", dict(gA=G_AS, gB=G_SA, p_plus=1.0)),
            ("|r| 0.1-1", dict(gA=G_AS, gB=G_SA, rabs=(0.1, 1.0))), ("|r| 0.01-10", dict(gA=G_AS, gB=G_SA, rabs=(0.01, 10.0)))]


def section2():
    rng = np.random.default_rng(20261016)
    print("\n=== 2. The measured rho = k_as / k_sa band and the natural-pair share (g_a > 0, c > 0) ===")
    rho = band(rng, G_AS, G_SA)
    q = np.percentile(rho, [5, 25, 50, 75, 95])
    p1 = float((rho < 1).mean())
    print(f"  rho percentiles 5/25/50/75/95: " + " ".join(f"{v:.3f}" for v in q) + f";  P(rho < 1) = {p1:.3f}  (sugar_pair assumed log-uniform 0.01-100, median 1)")
    check("P3", 0.44 <= q[2] <= 0.48 and 0.09 <= q[0] <= 0.12 and 1.9 <= q[4] <= 2.1 and 0.78 <= p1 <= 0.84,
          f"median {q[2]:.3f}, 5% {q[0]:.3f}, 95% {q[4]:.3f}, P(rho<1) {p1:.3f}")
    ref = exact(0.5, (0.03, 3.0), (0.01, 100.0))
    check("P1", abs(ref - 0.8636) <= 0.0005, f"sugar_pair primary closed form {ref:.4f}")
    print(f"\n{'variant':40} {'MC':>7} {'grid':>7} {'rescued(MC)':>12}")
    out = {}
    for lab, kw in VARIANTS:
        mc, resc = share_mc(rng, **kw)
        gr = share_grid(**kw)
        out[lab] = (mc, gr, resc)
        print(f"{lab:40} {mc:7.4f} {gr:7.4f} {resc:12.4f}")
    mc, gr, resc = out["PRIMARY"]
    check("P4", 0.828 <= mc <= 0.848 and abs(mc - gr) <= 0.004 and 0.016 <= ref - mc <= 0.036,
          f"MC {mc:.4f}, grid {gr:.4f}, sugar_pair {ref:.4f}, change {mc - ref:+.4f}")
    g = lambda lab: out[lab][1]
    rng_ok = {"(b) RNA templating only": (0.765, 0.790), "(c) Hein amplification (G_as to 1.0)": (0.850, 0.870),
              "(d) rho / 3": (0.770, 0.790), "(e) rho x 3": (0.890, 0.910), "p(r>0) = 0": (0.660, 0.690),
              "p(r>0) = 0.25": (0.750, 0.770), "p(r>0) = 0.75": (0.910, 0.930), "|r| 0.1-1": (0.895, 0.915),
              "|r| 0.01-10": (0.795, 0.815)}
    bad = {k: round(g(k), 4) for k, (lo, hi) in rng_ok.items() if not lo <= g(k) <= hi}
    half = [g(l) for l, _ in VARIANTS if not l.startswith("p(r>0)")]
    ok = (not bad) and abs(g("p(r>0) = 1") - 1) < 1e-9 and min(half) >= 0.75 and max(half) <= 0.92
    check("P5", ok, "outside range: " + (str(bad) if bad else "none") + f"; p(r>0)=0.5 variants span {min(half):.3f}-{max(half):.3f}")
    check("P6", 0.665 <= out["p(r>0) = 0"][1] <= 0.690, f"rescued share (r < 0 only) {out['p(r>0) = 0'][1]:.4f} (MC {resc:.4f} at primary)")


# ---------------------------------------------------------------- 3. coupling sign and spot checks of dna_network
def section3(runs=4000):
    print(f"\n=== 3. Sign of the coupling: dna_network.simulate_asymmetric, {runs} runs/case, g_a = 1e-3, kappa_bar = 0.3 ===")
    print(f"{'r':>7} {'rho':>6} {'c':>7} {'k_as':>7} {'k_sa':>7} {'predicted':>9} {'simulated':>9} {'matched':>8}")
    g_a, kb, eps, gamma = 1e-3, 0.3, 1e-4, 0.01
    cases = [(0.257, 0.46), (-0.3, 0.46), (-0.678, 0.46), (-1.0, 0.46), (-0.3, 0.15)]
    ok_c = ok_m = ok_p = True
    zero = None
    for i, (r, rho) in enumerate(cases):
        k_as, k_sa = kb * np.sqrt(rho), kb / np.sqrt(rho)
        thy = D.theoretical_p(k_sa, k_as, g_a, r * g_a, eps, gamma)
        c = float(cfac(r, rho))
        geff = (np.sqrt(k_as) * g_a + np.sqrt(k_sa) * r * g_a) / np.sqrt(k_as + k_sa)
        ok_c &= abs(geff / g_a - c) < 1e-9 and abs(float(norm.cdf(delta(c * g_a, eps, gamma))) - thy) < 1e-9
        m, p = D.simulate_asymmetric(k_sa, k_as, g_a, r * g_a, eps, gamma, runs=runs, seed=i + 11)
        print(f"{r:7.3f} {rho:6.2f} {c:7.3f} {k_as:7.3f} {k_sa:7.3f} {thy:9.3f} {p:9.3f} {m:8.3f}")
        ok_m &= m >= 0.99
        ok_p &= abs(p - thy) < 0.04
        if r == -0.678:
            zero = p
    check("P1", ok_c and RES.get("P1", True), "cfac equals g_eff / g_a (and Phi(Delta(c g_a)) equals theoretical_p) to 1e-9 in all five cases")
    mneg, _ = D.simulate_asymmetric(-0.3, -0.3, g_a, 0.3 * g_a, eps, gamma, runs=runs, seed=99)
    print(f"  both couplings NEGATIVE (k_sa = k_as = -0.3, r = 0.3): matched (alpha and beta same sign) = {mneg:.3f}")
    print("  -> the stable pair is anti-matched (L-aa with L-sugar, or D with D); the natural pair is never reached.")
    check("P7", ok_m and ok_p and zero is not None and abs(zero - 0.5) < 0.06 and mneg <= 0.01,
          f"positive: matched >= 0.99 {ok_m}, natural within 0.04 {ok_p}, c~0 case {zero:.3f}; negative couplings matched {mneg:.3f}")


def main():
    print(f"sugar_measured.py, {NDRAW} prior draws per row")
    section1()
    section2()
    section3()
    t = time.time() - T0
    check("P9", t < 300, f"{t:.0f} s")
    print(f"\n{sum(RES.values())} of {len(RES)} checked items HELD: " + ", ".join(f"{k} {'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
