"""POST HOC follow-up to blind_units.py: why A5 and A6 failed (Asakura 2002), and whether B2 and B3 (Bak 2023) can be scored.

NOT BLIND.  The project had already read Asakura 2002 (blind_units.py's RESULTS paragraph reports its Fig. 2, Fig. 4, Table 1,
aggregate counts and Discussion).  The writer of the predictions below had read that paragraph and blind_units.log, but not
the PDF text, which was converted to papers/blind_units_posthoc/asakura_2002.txt (pdftotext -layout) only after the
predictions were hashed.  E2 (secondary nucleation inside a mass) was therefore told to the writer, not predicted blind.
blind_units.py and blind_units.log are untouched; this file imports blind_units.py's data and functions unchanged.
One error in the hashed predictions file, left as is: its lines 5-7 put "the vote of about 90 or more crystallites" in the
paper's Discussion.  That sentence is blind_units.py's, not the paper's (X3 below).

PREDICTIONS.  blind_units_posthoc_predictions.txt
  shasum -a 256: 2dea583c0a6942f03325493938caee5afd814f5268112e6ad6f9e1dbdf77d4ba
  date:          Mon Oct 19 23:36:33 PDT 2026 (2026-10-20T06:36:33Z), before asakura_2002.pdf was converted to text.
self_check() asserts the file still has this hash.  Scoring rule fixed there: an explanation SURVIVES if none of its
predictions FAILED and at least one HELD on an item the paper gives that blind_units.py did not use.

RESULTS (blind_units_posthoc.log, run 2026-10-19).  "txt:N" = line N of papers/blind_units_posthoc/asakura_2002.txt.
 Asakura A6 (n_eff at 150 C far above 2.4 aggregates per run):
  E1 merged separate nuclei: does not survive.  P1 not met: the paper describes new nuclei forming near a crystal and joining
     the growing solid (txt:183-185), which is E2, not merging; A5 already showed one nucleus gives a racemic mass at 150 C.
  E2 secondary nucleation inside a mass: SURVIVES, but only as the reason a visible mass is not a unit.  P2a HELD (txt:155,
     275); P2c HELD in direction (first nucleus 979 s vs 3,150 s; Rnuc 3.25e-4 vs 1.01e-4 cm^-2 s^-1, 3.2x), while the
     effective units per aggregate would need to differ ~42x and the visible aggregates go the other way (24 vs 32).  P2b NOT
     TESTABLE: the paper counts no crystallites and says the crystals cannot really be counted (txt:183-185).
  E3 ee measurement: does not survive.  The solid was dissolved in 20 ml benzene and read in a polarimeter (txt:120-123); no
     subsampling is mentioned; precision is not stated (P3 NOT TESTABLE as worded); a detection floor would need to be 43.4 points.
  E4 estimator or tail: not supported.  P4a HELD (derived, from data blind_units.py used): no resample of the 150 C runs gives
     n_eff below 33.8, bootstrap 2.5th percentile 104.  P4b NOT TESTABLE (the paper gives the Langevin form, no fitted values).
 Added after reading, not scored and not pre-registered (X1-X3 below):
  X1 a Gaussian vote with n_eff 272 (sd 6.06 points) puts 7 of 10 runs within +-0.5 point of 0 with probability 5.3e-7, within
     +-1 point 5.5e-5, within +-2 points 4e-3: the 150 C runs look like a spike at 0 plus a tail, so the chi-square low end (88)
     is not a trustworthy bound.  Table 1 reports whole percents, so +-0.5 is the resolution.
  X2 Table 1 weights (unused): at 152 C Spearman(|ee|, weight) = -0.66 over 9 masses, permutation p = 0.07 (not significant at
     0.05); pooling the 150 C masses gives -0.71, p = 0.01, but that mixes temperature with size.  ee^2 x weight ranges from
     0 to 11 mg across the nine 152 C masses, so no single crystallite size fits.
  X3 "about 90 or more crystallites" is not in the paper (the only "90" is "ee greater than 90%" about Pincock et al.); it is
     this project's n_eff low end of 88.  The paper gives two readings of the 150 C solid: aggregates of many small single-handed
     crystals, or crystals with many lattice defects that take up both hands (txt:238-240); the Discussion also hedges that the
     solid "could" be small R and S crystallites with many lattice defects (txt:307).  It chooses neither and gives no count.
 Bak B2 and B3: both stay NOT SCORABLE.  The main text gives only the three main-series points and no N_L, N_D, N_mix; the ACS SI
  returned HTTP 403 to a plain request and Europe PMC holds no open copy or supplement.  Items needed are printed in the log.

Run: /opt/anaconda3/bin/python3 blind_units_posthoc.py > blind_units_posthoc.log     (from the repo directory)
"""
import hashlib
import re
from math import comb, erf, sqrt
from pathlib import Path

import numpy as np
from scipy.stats import chi2, spearmanr

import blind_units as bu

HERE = Path(__file__).resolve().parent
PRED_SHA = "2dea583c0a6942f03325493938caee5afd814f5268112e6ad6f9e1dbdf77d4ba"
ASA = HERE / "papers/blind_units_posthoc/asakura_2002.txt"   # pdftotext -layout of papers/asakura_2002.pdf
BAK = HERE / "papers/blind_units_posthoc/bak_2023.txt"       # pdftotext -layout of papers/bak_2023.pdf

# Unused items from asakura_2002.pdf; the txt line (split on "\n", as grep -n counts) is named in each comment.
T_FIRST = {150: 979, 152: 3150}          # s, average time to first nucleus (txt:168-170)
RNUC_PRINTED = {150: 3.25e-4, 152: 1.01e-4}   # cm^-2 s^-1 (txt:170-172), Rnuc = 1/(T A), A = 3.14 cm^2 (txt:159-166)
AREA = 3.14
# Table 1 (txt:132-151): (T, run, weight mg, ee %).  blind_units.AS_MASSES has the same ee values; weights were unused.
TABLE1 = [(150, 1, 93, 0), (150, 1, 31, 0), (150, 2, 118, 0), (150, 3, 165, 0),
          (152, 1, 37, 0), (152, 1, 16, 4), (152, 2, 13, 59), (152, 3, 35, 3), (152, 4, 27, 2),
          (152, 5, 113, 0), (152, 6, 58, 31), (152, 7, 16, 83), (152, 7, 36, 0)]

F95 = chi2.ppf(0.025, 10) / 10   # 0.325, as blind_units.py


def ln(path, needle):
    """1-based line of the first line of the pdftotext file containing needle (counted on "\\n", like grep -n); asserts it exists."""
    hit = [i + 1 for i, s in enumerate(path.read_text().split("\n")) if needle in s]
    assert hit, f"{needle!r} not in {path.name}"
    return hit[0]


def self_check():
    """Reproduce blind_units.log's n_eff values from blind_units.py's own data and functions; pin the predictions file."""
    assert hashlib.sha256((HERE / "blind_units_posthoc_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    log = (HERE / "blind_units.log").read_text()
    got = {}
    for t, ne, lo in re.findall(r"unstirred (\d+) C, \d+ runs:.*n_eff = ([\d.]+) \(95% low ([\d.]+)\)", log):
        assert round(bu.neff(bu.AS_UNSTIRRED[int(t)]), 2) == float(ne), (t, ne)
        assert round(bu.neff_low(bu.AS_UNSTIRRED[int(t)]), 2) == float(lo), (t, lo)
        got[t] = ne
    for r, ne in re.findall(r"(\d+) rpm, \d+ runs:.*n_eff = ([\d.]+)", log):
        assert round(bu.neff(bu.AS_STIRRED[int(r)]), 2) == float(ne), (r, ne)
        got[r + " rpm"] = ne
    assert set(got) == {"150", "152", "0 rpm", "200 rpm", "500 rpm"}, got
    assert "A5 FAILED" in log and "A6 FAILED" in log and "B2 NOT SCORABLE" in log and "B3 NOT SCORABLE" in log
    assert [abs(x[3]) for x in TABLE1] == [abs(m[2]) for m in bu.AS_MASSES], "Table 1 ee differs from blind_units.AS_MASSES"
    print("self-check OK: n_eff of", ", ".join(f"{k}={v}" for k, v in got.items()), "match blind_units.log; predictions sha256 matches")


def a6_numbers():
    e = bu.AS_UNSTIRRED[150]
    ne, lo, agg = bu.neff(e), bu.neff_low(e), bu.AS_AGG[150]
    print("=== The failure, from blind_units.py's own data (blind_units.py:133, 145) ===")
    print(f"  150 C: n_eff {ne:.1f}, 95% low {lo:.1f}, aggregates {agg}/run -> {ne / agg:.0f} units per aggregate "
          f"(low end {lo / agg:.1f}); 152 C: {bu.neff(bu.AS_UNSTIRRED[152]) / bu.AS_AGG[152]:.2f} per aggregate")
    need = agg / F95
    s3 = sum((x / 100) ** 2 for x in e)
    floor = sqrt((10 / need - s3) / 7)
    rms = np.sqrt(np.mean((np.array(e) / 100) ** 2))
    print(f"  low end <= {agg} needs n_eff <= {need:.2f}, i.e. rms ee >= {100 / sqrt(need):.1f} points (observed {100 * rms:.1f}); "
          f"the 7 zeros would have to be +-{100 * floor:.1f} points")
    assert abs(100 * floor - 43.4) < 0.1 and abs(need - 7.39) < 0.01
    print(f"  ratio of per-aggregate units 150 C / 152 C: {(ne / agg) / (bu.neff(bu.AS_UNSTIRRED[152]) / bu.AS_AGG[152]):.0f}x")
    rng = np.random.default_rng(0)
    a = np.array(e) / 100
    m2 = (a[rng.integers(0, 10, (20000, 10))] ** 2).mean(1)
    lo_bs = 1 / np.quantile(m2, 0.975)
    print(f"P4a (derived): smallest possible n_eff from any resample = {1 / a.max() ** 2:.1f} (max ee^2 {a.max() ** 2:.4f}); "
          f"bootstrap (20000, seed 0) 2.5th percentile n_eff = {lo_bs:.1f} vs chi-square low end {lo:.1f}; both far above {agg}")
    assert 1 / a.max() ** 2 > 33 and lo_bs > agg
    return ne


def score_asakura():
    print("\n=== Asakura 2002: pre-registered predictions (blind_units_posthoc_predictions.txt), scored by the fixed rules ===")
    for t in (150, 152):
        r = 1 / (T_FIRST[t] * AREA)
        print(f"  Rnuc({t} C) = 1/({T_FIRST[t]} s x {AREA} cm^2) = {r:.3g} cm^-2 s^-1 (printed {RNUC_PRINTED[t]:.3g}; txt:168-172)")
        assert abs(r / RNUC_PRINTED[t] - 1) < 0.01
    print(f"  Rnuc ratio 150/152 = {RNUC_PRINTED[150] / RNUC_PRINTED[152]:.2f}x; units per aggregate need ~"
          f"{(bu.neff(bu.AS_UNSTIRRED[150]) / 2.4) / (bu.neff(bu.AS_UNSTIRRED[152]) / 3.2):.0f}x; visible aggregates 24 vs 32 (txt:191-194)")
    L = {n: ln(ASA, n) for n in ("cannot really count the number of crystals", "recurrence of nucleations on or near the surface",
                                 "20 ml of benzene", "7 out of the 10 crystallizations", "polycrystalline aggregates consisting of both R- and S-",
                                 "very small R- and S-crystals", "large number of lattice", "tallites with many lattice defects", "0.50 g",
                                 "mostly racemic", "ee as high as 17%")}
    print("  quote-check lines (txt):", L)
    print("  P1 (E1: paper describes merging of separately nucleated masses at 150 C): NO.  Nearest text says new nuclei form near a crystal and "
          f"join the growing solid (txt:{L['cannot really count the number of crystals'] - 2}-{L['cannot really count the number of crystals']}), "
          "i.e. secondary nucleation.  E1 does not survive.")
    print(f"  P2a (150 C solid = many small crystallites of both hands within a mass): HELD.  Discussion: at high supercooling the solid is a "
          f"polycrystalline aggregate of very small R and S crystals (txt:{L['very small R- and S-crystals']}); the stirred-crystallization section "
          f"allows either many small single-handed crystals or crystals with many lattice defects taking up both hands (txt:{L['large number of lattice']}-"
          f"{L['large number of lattice'] + 2}); the paper does not choose.")
    print("  P2b (independent crystallite count >= 37 per aggregate or 88 per run): NOT TESTABLE.  No count or size at 150 C; Fig. 1 is 152 C only, "
          f"rings 1.0 mm apart, 'cloudy' (txt:77-81); the paper says the crystals cannot really be counted (txt:{L['cannot really count the number of crystals']}).")
    print("  P2c (faster or finer at 150 C): HELD in direction (first nucleus 979 s vs 3,150 s, Rnuc 3.2x); the size of the effect is not matched.")
    print(f"  P3 (ee method and precision <= +-5 points): method found (JASCO polarimeter, the solid 'thus obtained' dissolved in 20 ml benzene, 1.0 dm "
          f"cell, txt:120-{L['20 ml of benzene'] + 1}; no subsampling mentioned); precision NOT stated: NOT TESTABLE as worded.  The only route for E3 "
          "(aliquot averaging) is absent, so E3 does not survive.")
    print("  P4b (paper's fitted low-stirring distribution SD <= 10.7 points): NOT TESTABLE; the paper gives the Langevin form and no fitted values.  "
          "(Its printed SD, 6% at 150 C, is the data blind_units.py already reproduced: blind_units.py:132.)")
    print("  Survival by the fixed rule: E1 no; E2 YES (P2a, P2c HELD, none FAILED); E3 no; E4 not supported (P4a is derived from used data, "
          "not an unused item).")
    print("  Reading: the visible aggregate is not the unit, so 'n_eff 272' cannot be compared with 2.4 aggregates.  The paper does not support "
          "'n_eff counts independent crystallites': it counts none, and the near-zero ee of 150 C masses is also what lattice-defect incorporation "
          "of both hands would give (X3).  What E2 does not explain: the 42x per-aggregate jump over 2 K, and 7 of 10 runs at exactly 0 (X1).")


def x1_spike(ne):
    print("\n=== X1 (added after reading; not pre-registered, not scored): is the 150 C set a Gaussian vote with n_eff = %.0f? ===" % ne)
    sd = 100 / sqrt(ne)
    for r in (0.5, 1.0, 2.0):
        p = erf(r / (sd * sqrt(2)))
        tail = sum(comb(10, j) * p ** j * (1 - p) ** (10 - j) for j in range(7, 11))
        print(f"  Gaussian sd {sd:.2f} points: P(|ee| < {r}) = {p:.3f}; P(>= 7 of 10 runs that close to 0) = {tail:.1e}")
    print(f"  observed: 7 of 10 runs at 0.0 (Fig. 2 caption, txt:{ln(ASA, '7 out of the 10 crystallizations')}); the other three 17.2, 8.4, 1.1.  "
          "Table 1 is in whole percents, so +-0.5 is the resolution.")


def x2_weights():
    print("\n=== X2 (added after reading; not pre-registered): Table 1 weights (unused) vs |ee| of single-nucleus masses (txt:132-151) ===")
    for T in (152, None):
        rows = [r for r in TABLE1 if T is None or r[0] == T]
        w = np.array([r[2] for r in rows], float)
        a = np.array([abs(r[3]) for r in rows], float)
        rho = spearmanr(w, a)[0]
        rng = np.random.default_rng(0)
        perm = [spearmanr(w, rng.permutation(a))[0] for _ in range(5000)]
        p = float(np.mean(np.abs(perm) >= abs(rho) - 1e-12))
        print(f"  {'152 C only' if T else 'both temperatures (mixes temperature with size)'}: n = {len(rows)}, "
              f"Spearman(|ee|, weight) = {rho:+.2f}, permutation p = {p:.2f} (5000 shuffles, seed 0)")
    r152 = [r for r in TABLE1 if r[0] == 152]
    print("  ee^2 x weight (mg), 152 C, in Table 1 order:", ", ".join(f"{(r[3] / 100) ** 2 * r[2]:.2f}" for r in r152),
          "(under a Gaussian vote with one crystallite size these would be equal in expectation)")
    print(f"  150 C masses weigh 31-165 mg (6-33% of the 0.50 g charge, txt:{ln(ASA, '0.50 g')}) and all read 0")


def x3_ninety():
    lines = ASA.read_text().split("\n")
    print("\n=== X3 provenance of 'about 90 or more crystallites' (blind_units.py RESULTS) ===")
    for i, s in enumerate(lines):
        for m in re.finditer(r"\b90\b", s):
            print(f"  '90' at txt:{i + 1}: ...{s[max(0, m.start() - 60):m.end() + 3].strip()}")
    assert not re.search(r"90 or more|about 90 crystal", ASA.read_text())
    print("  none is a crystallite count; the phrase is the project's n_eff low end (88.33), not the paper's.")


def bak():
    print("\n=== Bak 2023: B2 and B3 ===")
    V = np.array([r[0] for r in bu.BAK_MAIN], float)
    P = np.array([r[1] for r in bu.BAK_MAIN], float)
    gate = len(V) < 4 or P.max() < 0.6 or P.min() > 0.4     # blind_units.py:218, unchanged
    print(f"  B2 gate (blind_units.py:218, unchanged): {len(V)} main-series markers, P_CSB {P.min():.2f}-{P.max():.2f} -> "
          f"{'NOT SCORABLE' if gate else 'scorable'}")
    assert gate
    V0 = bu.exp_fit(V, P)
    print(f"  to reach P_CSB <= 0.4 the series would need Vc about {-V0 * np.log(0.4):.1f} uL if exp(-V/V0) holds with V0 = {V0:.2f} uL "
          "(an assumption, not data)")
    txt = BAK.read_text()
    assert re.search(r"P\s*CSB\s*=", txt) and "Nmix" in txt and not re.search(r"N_?L\s*=\s*\d", txt)
    print(f"  bak_2023.pdf main text: P_CSB defined as (N_L + N_D)/N_tot with N_mix (txt:{ln(BAK, 'CSB probability was 100%')}-{ln(BAK, 'number of sessile drops evaluated')}); "
          f"values given only as 0.8 at 0.87 uL and 0.7 at 1.72 uL, with some drops holding both hands (txt:{ln(BAK, 'PCSB = 0.8 at Vc = 0.87')}-"
          f"{ln(BAK, 'contained both L and D crystals')}); no N_L or N_D value anywhere.")
    print(f"  Figure S1a-d are cited as the supporting images (txt:{ln(BAK, 'Figure S1a')}-{ln(BAK, 'kinetics (Figure S1c)')}); they are in the SI, not in the PDF.")
    print("  Fig. 2a labels give (k/n) only for P < 1 markers; the P = 1 markers share one '(66/66)' label (blind_units.py:124).")
    print(f"  Fig. 3 (caption txt:{ln(BAK, 'Figure 3. Determination')}) plots supersaturation vs volume per drop, not N_L / N_D, and names no series.")
    print("  ACS SI https://pubs.acs.org/doi/suppl/10.1021/acs.jpclett.2c03873/suppl_file/jz2c03873_si_001.pdf: plain HEAD request returned HTTP 403 "
          "on 2026-10-19 (Cloudflare); not retried, no workaround.  Europe PMC REST (DOI 10.1021/acs.jpclett.2c03873): pmid 36652610, "
          "inEPMC N, isOpenAccess N, hasSuppl N.")
    print("  B3 NOT SCORABLE unchanged (blind_units.py:223-227 needs BAK_LD = (all-L, all-D)).  B2 NOT SCORABLE unchanged.")
    print("  ITEMS NEEDED (open the SI in a browser, or ask the corresponding authors Su Yeon Bak and Gerard Coquerel):")
    print("   (1) the SI PDF jz2c03873_si_001.pdf: any table giving N_L, N_D, N_mix per Vc bin or per series.  Figure S1 is described as images,")
    print("       so expect to need (2).")
    print("   (2) the per-drop record behind Fig. 2a for the main series (Tsat 50 C, Tinc 25 C, FC-40, PMMA): for each drop Vc and outcome L / D / mixed.")
    print("       That one table makes B3 (all-L and all-D counts) scorable and, with at least 4 volume bins and one at P_CSB <= 0.4, B2 too.")
    print("   (3) for B2 only: if the series has no drop above about 1.7 uL, drops at larger Vc (about 4 uL or more) in the same conditions.")


def main():
    self_check()
    ne = a6_numbers()
    score_asakura()
    x1_spike(ne)
    x2_weights()
    x3_ninety()
    bak()


if __name__ == "__main__":
    main()
