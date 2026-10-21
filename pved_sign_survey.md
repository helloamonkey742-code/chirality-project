# Sign of the amino-acid PVED: what the published calculations say

Compiled 2026-10-16. Sources were abstracts and web pages only (Crossref, OpenAlex, PubMed, PMC page text, Semantic Scholar, a KAKEN grant record). No PDFs were downloaded. Every DOI below was checked directly with `https://api.crossref.org/works/<doi>` (no email in the URL) and returned the expected title, authors, journal and pages.

Sign convention: "L" means the calculation puts the natural L enantiomer lower in energy. "D" means it puts L higher. Epv negative = L stable in the Kikuchi papers; I confirmed that convention only for Kikuchi and Wang 1990 (abstract) and Watanabe 2004 (secondhand, below). For Quack 2022 I report the sign of the parity-violating potential as the review words it and infer the L/D reading.

Evidence grade: P = read the abstract (primary). S = secondhand (search-engine summary of the page, a citing paper, or the paper's own text in the student draft). T = title and bibliographic data only; no result seen.

## Table 1. Calculations that give an amino-acid sign

| # | Paper (DOI) | Molecules | Method | Phase / conformation | Sign | Magnitude | Grade |
|---|---|---|---|---|---|---|---|
| 1 | Mason & Tranter 1983, Chem Phys Lett 94, 34 (10.1016/0009-2614(83)87205-4); 1984, Mol Phys 53, 1091 (10.1080/00268978400102881); 1985, Proc R Soc A 397, 45 (10.1098/rspa.1985.0003) | amino acids, early set | semi-empirical / minimal-basis SCF, first generation | gas-phase-type single geometries | L | about 0.1 aeV scale (per the draft) | S |
| 2 | Tranter 1985, Chem Phys Lett 120, 93 (10.1016/0009-2614(85)87019-6); Mol Phys 56, 825 (10.1080/00268978500102741) | alpha-amino acids | same family of methods | single geometries; not shown to be aqueous | L in each case reported | not extracted | S |
| 3 | Tranter 1986, J Theor Biol 119, 467 (10.1016/s0022-5193(86)80196-5); MacDermott 1995, Orig Life Evol Biosph 25, 191 (10.1007/BF01581582) | biomolecules incl. amino acids | as above / later SCF | mostly single conformers | L in nearly all (MacDermott 1995 as summarised in the draft and a review) | not extracted | S |
| 4 | Kikuchi & Wang 1990, Bull Chem Soc Jpn 63, 2751 (10.1246/bcsj.63.2751) | Gly, Ala, Ser zwitterions | HFO-5G/6G minimal bases vs split-valence | zwitterion, scanned over C-C rotation | basis-set dependent: L (negative) for L-Ala at the stable conformer with HFO bases, opposite sign with split-valence bases; L-Ser L-favoured over nearly the whole rotation range | not extracted | P |
| 5 | Zanasi & Lazzeretti 1998, Chem Phys Lett 286, 240 (10.1016/s0009-2614(98)00108-0) | a few L-alpha-amino acids | CPHF / RPA-level PVES | not specified | L (title and a search summary of the abstract) | not extracted | S |
| 6 | Zanasi, Lazzeretti, Ligabue, Soncini 1999, Phys Rev E 59, 3382 (10.1103/PhysRevE.59.3382) | L-valine (also hydrated D-glyceraldehyde); draft says five amino acids | RPA, about 10x larger than older methods | draft says sign changes with dihedral angle, most negative at the water-preferred conformer | L at the water-preferred conformer; conformer-dependent elsewhere (the dihedral-angle statement is from the draft, not the abstract) | about one order of magnitude above earlier values; draft quotes up to 11.4 aeV | P (abstract) + S (dihedral claim) |
| 7 | Berger & Quack 2000, ChemPhysChem 1, 57 (10.1002/1439-7641(20000804)1:1<57::aid-cphc57>3.0.co;2-j) | alanine, zwitterion model | MC-LR / RPA, conformational scan, solvent effects | gas and condensed phase (continuum) | no clear sign: the abstract says many earlier L-stabilisation claims are unjustified | not extracted | P |
| 8 | Laerdahl, Wesendrup, Schwerdtfeger 2000, ChemPhysChem 1, 60 (10.1002/1439-7641(20000804)1:1<60::aid-cphc60>3.0.co;2-m) | alanine, 13 gas-phase conformers | relativistic (Dirac-Hartree-Fock type) | gas phase, neutral conformers | conformer-dependent: L in only 7 of 13 | not extracted | P |
| 9 | Wesendrup, Laerdahl, Compton, Schwerdtfeger 2003, J Phys Chem A 107, 6668 (10.1021/jp022568v) | several amino acids | relativistic parity-nonconserving calculations | gas phase | no clear evidence that L is stabilised | not extracted | P |
| 10 | Watanabe, Morihashi, Takahashi, Kitayama, Yagi, Kikuchi 2004, J Mol Struct THEOCHEM 671, 119 (10.1016/j.theochem.2003.10.060) | L-alanine zwitterion, L-Ala-(H2O)n, n=1-10 | MC/MO spin-orbit-coupling SCF; 100 Monte Carlo structures; also a generalized-Born continuum map over NH3 and CO2 torsions | explicit water, 298 K | D: averaged Epv positive, L less stable than D; a small change in how water faces CO2 can flip the sign | about 2 x 10^-20 hartree; the exponent sign was garbled in the page I saw, so treat the magnitude as unverified | S (KAKEN grant record; the journal abstract was not reachable) |
| 11 | MacDermott, Fu, Hyde, Nakatsuka, Coleman 2009, Orig Life Evol Biosph 39, 407 (10.1007/s11084-009-9161-x) | alanine and other amino acids | CPHF and DFT, partial optimisation of torsions | gas phase; solution structures taken from other groups' solvated optimisations | gas phase: L, "decisively" (neutral conformers); solution: sign cannot be fixed because of torsion-angle sensitivity | not extracted | P |
| 12 | MacDermott, Fu, Nakatsuka, Coleman, Hyde 2009, Orig Life Evol Biosph 39, 459 (10.1007/s11084-009-9162-9) | four Murchison alpha-methyl amino acids | same | gas phase, neutral | L for all four | not extracted | P |
| 13 | Quack, Seyfang, Wichmann 2022, Chem Sci 13, 10598 (10.1039/d2sc01323a) | zwitterionic L-alanine (Fig. 21), after refs the review attributes to earlier work | RPA gas phase, SCRF-RPA solution, CASSCF-LR | gas and continuum solution, carboxylate rotation | conformer-dependent: parity-violating potential positive for roughly 50-140 degrees of carboxylate rotation and negative elsewhere; the ordinary-energy minimum sits near 60 degrees, inside the positive window (my inference: L is the higher-energy form at the energy minimum if positive means L destabilised; not stated in the text I saw) | gas and solution differ by about as much as a change of theory level; alanine value 6.2 aeV is the draft's table entry, not seen by me | S (PMC page, truncated at the later sections) |
| 14 | Faglioni, Passalacqua, Lazzeretti 2005, Orig Life Evol Biosph 35, 461 (10.1007/s11084-005-3511-0) | alpha-helical polypeptides | chain-polymer algorithm | helix, not aqueous | abstract gives magnitude and an evolutionary model but no sign | not extracted | P |

## Table 2. Checked, no amino-acid sign obtained (do not cite as sign evidence)

| Paper (DOI) | Why not |
|---|---|
| Bakasov, Ha, Quack 1998, J Chem Phys 109, 7263 (10.1063/1.477360); erratum 1999 (10.1063/1.478513) | Method paper; I did not see an amino-acid result. T only. |
| Berger, Quack, Tschumper 2000, Helv Chim Acta 83, 1919 (10.1002/1522-2675(20000809)83:8<1919::aid-hlca1919>3.0.co;2-d) | Cyanoaziridine and cyanooxirane precursors, not amino acids. P. |
| Kitayama et al. 2002, J Mol Struct THEOCHEM 589-590, 183 (10.1016/s0166-1280(02)00255-5) | SOC-SCF method paper; the alanine-in-water application is the Watanabe 2004 paper. T only. |
| Kikuchi et al. 1994, J Mol Struct THEOCHEM 305, 79 (10.1016/0166-1280(94)80143-6) | Generalized-Born conformations, not PVED. T only. |
| Faglioni, Cuesta 2010, Orig Life Evol Biosph 41, 249 (10.1007/s11084-010-9229-7); Faglioni et al. 2020 (10.1007/s11084-020-09598-6) | Protein secondary structure and metal centres; no abstract and no sign seen. T only. |
| Senami & Shimizu 2020, Phys Lett A 384, 126796 (10.1016/j.physleta.2020.126796) | Electron chirality of Ala, Ser, Val, not PVED sign. S. |
| Aucar, Stroppa, Aucar 2023, J Phys Chem Lett 15, 234 (10.1021/acs.jpclett.3c03038); arXiv 2603.01885 (2026) | PVED versus electronic chirality measure; page text gives no amino-acid sign. S. |
| Sitha 2021, Comput Theor Chem 1205, 113446 (10.1016/j.comptc.2021.113446) | Alanine in a dipole and oriented electric field; no abstract available to me. T only. |
| Quack, Stohner, Willeke 2008 (10.1146/annurev.physchem.58.032806.104511); Berger & Stohner 2019 (10.1002/wcms.1396) | Reviews; bibliographic data only, no sign statement extracted. T only. |

Named in the brief but not found: a solvated or zwitterionic amino-acid PVED sign paper from Sahu/Gopakumar, Thierfelder, or Wang/Hennig. Searches returned nothing I could confirm, so none are cited. I also found no 2015-2026 study that computes an amino-acid PVED sign in water; the 2022 review (row 13) is the most recent source.

## Conclusion

Counting only studies that treat a zwitterion, explicit water, or a continuum solvent (rows 4, 6, 7, 10, 11, 13):

- Clearly L: 1 of 6 (Zanasi 1999, at the conformer preferred in water; the abstract itself speaks only of valine and a sugar, so even this rests on the draft's description).
- Basis-set dependent, leaning L with minimal bases: 1 of 6 (Kikuchi and Wang 1990).
- Clearly D: 1 of 6 (Watanabe 2004, explicit water, secondhand).
- Sign depends on conformer or cannot be fixed: 3 of 6 (Berger and Quack 2000; MacDermott 2009 for the solution part; Quack 2022).

So the fraction favouring L is between 1 in 6 and 2 in 6 depending on how Kikuchi and Wang is read, and no aqueous study is both high-level and decisive.

Gas phase is better but not closed either. Early work (rows 1-3) and the 2009 MacDermott papers favour L for neutral gas-phase structures. Relativistic work disagrees: L in 7 of 13 alanine conformers (Laerdahl 2000), and no clear L preference across amino acids (Wesendrup 2003). Berger and Quack 2000 call earlier L claims unjustified.

Does recent high-level work agree with L? No. The higher-level treatments (relativistic Dirac-type, MC-LR, CASSCF-LR, Quack 2022) report a sign that follows the conformation, and the one explicit-water study gives D. The mechanism is the same everywhere: the carboxylate and NH3 torsions flip the sign, and the solution conformation is not known to the accuracy needed. MacDermott 2009 argues amino acids are a special case because their favoured conformers are nearly achiral; that is an argument for why the effect might still be L, not a calculation that shows it in water.

Is "unsettled" still fair? Yes. A firmer sentence in the L direction is not justified. A firmer sentence against L is not justified either, since the D result comes from one secondhand HF-level study. The defensible firmer wording is: the gas-phase sign is mostly but not uniformly L, the aqueous-zwitterion sign is conformer-dependent and has not been fixed by any calculation, and the one explicit-water study points to D. One thing in the draft's favour: this supports the paper's point that the L-amino-acid plus D-sugar pairing rests on a sign that the literature has not established.

## Caveats

- Most magnitudes were not visible in abstracts; I left them blank rather than guess.
- Rows marked S need the full text before any sentence in the paper leans on them, especially rows 5, 6 (dihedral claim), 10 and 13.
- The Watanabe 2004 exponent in the KAKEN text lost its minus sign in the page I read.
- Tranter 1985 exists as two papers (Chem Phys Lett 120, 93 and Mol Phys 56, 825); the draft cites only the first.
- Author initials from Crossref: Watanabe T., Morihashi K., Takahashi O., Kitayama T., Yagi T., Kikuchi O.
