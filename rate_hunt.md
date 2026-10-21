# rate_hunt.md: measured rates for the two "?" cells in scorecard.py

Literature work only. scorecard.py and aped.py were read, not edited. Written 2026-10-19.
Sources: Europe PMC full-text XML (open access), Crossref, PubMed E-utilities. Nothing paywalled was opened.
Downloaded XML is kept in `papers/rate-hunt/xml_PMC*/` (untrusted data, read only with `python3 -I`).

## Result

| cell | letter implied | confidence | why |
|---|---|---|---|
| APED R5 (NCA activation, `a`) | **stays ?** | high that no OA number exists; low on any upgrade | No measured NCA-formation rate was retrievable open-access. The only lab bracket is the Leman 2004 abstract wording, converted below, and it is at lab COS levels, not at a setting's c. Update 2026-10-21: Leman 2004 and its SOM were read in full; still ?, see 1b. |
| Joyce 1984 R5 (template-directed rate) | **magnitude passes every setting; sign argues N** | medium on magnitude, medium-low on N | 2-methylimidazolide RNA extension k' = 2.8e-5 to 1.4e-4 /M/s clears 3.17e-7 (pond, lake) by 88x to 438x and 3.17e-9 (ocean) by about 9e3x to 4e4x. But Joyce 1984 is cross-inhibition (negative feedback), so no autocatalytic k2 exists for the floor to test. |

The Joyce 1984 numbers themselves were not read (Nature 310:602-604, doi:10.1038/310602a0, paywalled). The number used is from a follow-up on the same leaving-group chemistry (below).

**Update 2026-10-21: the Joyce row is superseded.** Joyce 1984 has now been read (`joyce_r5.py`, predictions hashed before reading) and R5 = Y is in scorecard.py; see 2e. The 2d argument for N is kept below as written, with the reason it was not followed.

## 0. What aped.py needs (receipts)

- aped.log L89 `[P10]`: s* = sigma/a = 0.0063166 at ctot = 2 M. Pass test in aped.py: s*·a >= 10/tau.
- aped.log L99-L101 (floors, YR = 3.156e7 s):

| setting | c | tau | floor 10/tau (/s) | k2 floor (/M/s) | `a` needed (/s) |
|---|---|---|---|---|---|
| pond | 1e-2 M | 1e2 yr | 3.17e-9 | 3.17e-7 | 5.02e-7 |
| lake | 1e-4 M | 1e4 yr | 3.17e-11 | 3.17e-7 | 5.02e-9 |
| ocean | 1e-6 M | 1e8 yr | 3.17e-15 | 3.17e-9 | 5.02e-13 |

- Check: a needed = floor / s* = 3.17e-9 / 0.0063166 = 5.02e-7 (pond). Same arithmetic for the others.
- aped.log L91-L97: the unstable window needs ctot >= 1 M (none at 0.5 M), but the scorecard settings are 1e-2 to 1e-6 M. The pass verdicts compare sigma at 2 M against each setting's own floor. Any measured `a` is only meaningful if it is the rate per amino acid at the amino-acid concentration the model assumes. Caveat on every APED verdict, not just this one.

## 1. APED R5: NCA activation rate `a` (/s)

### 1a. Measured, open access: none for formation

| item | value | source and location | use |
|---|---|---|---|
| Val-NCA hydrolysis (a loss rate, not formation) | 0.018 /s at pH ~8.5; t1/2 = ln2/0.018 = 38.5 s | PMC7547885 (Devaraj-group JACS 2016, acylation with fatty acids) main text, "observed rate constant of hydrolysis of 1 at pH ca. 8.5 as 0.018 s-1 (see Figure S7)"; temperature not stated in the passage I read | Not a value for `a`. It says NCAs live seconds to minutes at pH >= 8, so activation must be continuously refreshed (the paper itself says so). |

### 1b. Bracket from an abstract only (derived, NOT a measurement)

Leman, Orgel & Ghadiri 2004, Science 306:283-286, doi:10.1126/science.1102722 (Crossref-verified). The abstract says NCA/peptide conversion "up to 80% in minutes to hours at room temperature". The body (COS pressure, pH, amino acid concentration, tables) is paywalled and **unread**.

**Update 2026-10-21:** the full text has now been read, from a copy sent by C. McKay on 2026-10-20. Table 1 gives only single-time end points, with no time course. The lab dose is about 25 mM dissolved COS. Note 14 caps volcanic COS at 0.09 mol%, and that is the sourced prebiotic level this section lacked. The bracket below is superseded by `aped_r5.py` and `journals/new_analyses_aped_r5.md`. The supporting material (`leman.som.pdf`) was read later the same day, from a copy the author downloaded (science.org had shown me a Cloudflare check, which I did not bypass). It prints no NCA formation rate or NCA time course either. Its dipeptide time courses (Figs. S3, S4, from preformed thiocarbamate) are digitized in `journals/new_analyses_aped_r5.md` section 2b.

If 80% of the amino acid is converted by first-order activation, a_lab = ln(1/0.2)/t = 1.609/t:

- t = 1 h = 3600 s: a_lab = 1.609/3600 = 4.5e-4 /s
- t = 10 min = 600 s: a_lab = 1.609/600 = 2.7e-3 /s

The 10 min figure is my reading of "minutes", not the paper's. Treat the pair as an order-of-magnitude bracket, "best case", at unknown lab COS.

If `a` scales linearly with COS (an assumption, not measured), a setting fails only if prebiotic COS is below a_lab-to-needed times lower than lab COS:

| setting | a needed (/s) | lab COS may fall by (using 4.5e-4) | (using 2.7e-3) |
|---|---|---|---|
| pond | 5.02e-7 | 4.5e-4/5.02e-7 = 9e2 | 2.7e-3/5.02e-7 = 5e3 |
| lake | 5.02e-9 | 9e4 | 5e5 |
| ocean | 5.02e-13 | 9e8 | 5e9 |

So the question is the prebiotic COS level against the lab level, which no source I could read gives for the Leman conditions.

### 1c. Other routes, no rate retrievable

| route | reference | status |
|---|---|---|
| cyanate / N-carbamoyl amino acids | Danger, Boiteau, Cottet, Pascal, JACS 2006, 128:7412-7413, doi:10.1021/ja061339+ ("The Peptide Formation Mediated by Cyanate Revisited. N-Carboxyanhydrides as Accessible Intermediates"; title and DOI Crossref-verified) | paywalled, **unread** at first. Update 2026-10-21: read from a copy sent by C. McKay; end points only, no rate or half-life (section 3 row and `aped_r5.py`). The SI is unread |
| nitrosation of N-carbamoyl amino acids | Pascal and coworkers; I could not recover the source passage for the "quantitative in under 1 h" claim after compaction | **unverified**, not used; a yield-time statement would give only a lower bound on `a` anyway |
| NCA by other chemistry | Collet, Bied, Mion, Taillades, Commeyras, Tetrahedron Lett 1996, 37:9043-9046, doi:10.1016/S0040-4039(96)01974-0 (found in the reference list of PMC12188423) | synthetic route, paywalled, unread, not prebiotic |
| carbodiimides | no OA kinetics found | none |

Searches that came back with nothing usable: Europe PMC for (NCA or N-carboxyanhydride) with COS, cyanate, nitrosation or carbodiimide plus rate or kinetics, open access only (top hits were polymer chemistry or reviews without rate tables).

### 1d. Verdict for APED R5

No setting can be decided from a retrieved number. If the Leman bracket held at prebiotic COS, pond/lake/ocean would all pass (lab margins 9e2x and up). If prebiotic COS is far below lab, only the ocean (needs 5.02e-13 /s) is likely to pass. **Keep R5 = ?.** A full-text read of Leman 2004 (pressure of COS, amino acid concentration, time courses) plus a defensible prebiotic COS level is the one next step that could move it. Update 2026-10-21: that step is done (1b update). R5 stays ?. No NCA formation rate is printed, in the paper or its SOM, and the assumption that `a` scales with dissolved COS is not measured (`journals/new_analyses_aped_r5.md` sections 5 and 7). The SOM time courses raise the lab-dose lower bound to 1.2e-6 /s (post hoc, section 2b of that note); R5 stays ? for the dose and model-window reasons given there.

## 2. Joyce 1984 R5: template-directed rate, converted to effective k2

### 2a. The convention being tested

ocean.py (L9, L44-L46): rate k = k2·c, so k2_eff = rate / c. Floors from scorecard.py are k2 >= 10/(c·tau): 3.17e-7 /M/s (pond, lake), 3.17e-9 /M/s (ocean). Conversion used throughout: M-1 h-1 / 3600 = M-1 s-1.

### 2b. Measured values, open access, same chemistry (2-methylimidazole leaving group)

**Kervio, Sosson & Richert 2016, NAR (PMC4937335), Tables 4 and 5.** Second-order rate constant k' for enzyme-free primer extension, MALDI-MS, 20 degC, HEPES 200 mM, NaCl 400 mM, MgCl2 80 mM. Table 5 (RNA; 36 uM template:primer; monomer 30-60 mM, pH 7.7 for MeIm):

| monomer | k' (M-1 h-1) | k' (/M/s) = k'/3600 | x pond/lake floor 3.17e-7 | x ocean floor 3.17e-9 |
|---|---|---|---|---|
| MeIm-UMP | 0.1 | 2.78e-5 | 88 | 8.8e3 |
| MeIm-GMP | 0.3 | 8.33e-5 | 263 | 2.6e4 |
| MeIm-CMP | 0.5 | 1.39e-4 | 438 | 4.4e4 |
| MeIm-AMP | 0.1 | 2.78e-5 | 88 | 8.8e3 |

Table 4 (DNA template, 3.6 mM monomer, pH 7.0): MeIm-dTMP 8, dGMP 51, dCMP 22, dAMP 13 M-1 h-1 = 2.2e-3, 1.4e-2, 6.1e-3, 3.6e-3 /M/s (7e3x to 4.5e4x the pond floor). Largest in the paper: OAt-dCMP 370 M-1 h-1 = 0.103 /M/s (abstract gives the 0.1-370 range).

The relevant-to-Joyce row is the RNA one, so I use 2.8e-5 to 1.4e-4 /M/s. Monomer Kd for the template is 11 to >500 mM (Tables 2-3), so at the lake and ocean concentrations (1e-4, 1e-6 M) the reaction is in the linear regime, where a second-order k' is the right quantity. For the pond (1e-2 M) it is borderline for the tightest binders.

**Ding, Zhou, Giurgiu & Szostak 2021 (PMC8754633).** kobs for imidazolium-bridged dinucleotides at 5 mM, range 0.35 /h (U*U) to 84 /h (C*G), with 1.7 /h for A*U (text of Results, Figure 2). Dividing by the intermediate concentration: 0.35/3600/5e-3 = 1.9e-2 /M/s; 84/3600/5e-3 = 4.7 /M/s. kobs/c underestimates the low-c second-order constant when binding is saturating (kobs = kcov·c/(Kd+c), so kobs/c = kcov/(Kd+c) < kcov/Kd), so these are lower bounds, but c here is the bridged intermediate, not total monomer, so per-monomer values are lower. Not used for the verdict.

**Walton & Szostak 2017 (doi:10.1021/acs.biochem.7b00792).** Already in the repo (rna_copying.py): kobs = 5.0 to 9.7 /h at 24 mM Cp*, Fig. 8a. Same convention: 5.0/3600/0.024 = 5.8e-2 and 9.7/3600/0.024 = 1.1e-1 /M/s. Not re-read here.

### 2c. Joyce 1984 itself

Joyce, Visser, van Boeckel, van Boom, Orgel & van Westrenen, Nature 310:602-604, doi:10.1038/310602a0 (Crossref-verified). Paywalled, **unread**. I found no kobs for it in any OA follow-up. Schmidt, Nielsen & Orgel 1997 (papers/schmidt1997/19970024833.txt, Fig. 2/3, conditions: 4 degC, 1.2 M NaCl, 0.2 M MgCl2, 5 mM (dC)10 template, 0.2 M 2,6-lutidine pH 7.8, 1 d, 0.1 M D-2-MeImpG) gives yields after 1 d, no rate constants.

### 2d. Verdict for Joyce R5

- Magnitude: every converted value is at least 88x above the pond/lake floor and at least 8.8e3x above the ocean floor, at 20 degC with the leaving group Joyce used. So "no rate was found" no longer explains the "?"; the rate is not the obstacle.
- Sign: R5 in scorecard.py asks for a growth rate of the enantiomeric excess, with N where "a model at measured rates gives an excess that decays (k2 = k/c <= 0)". Joyce 1984 is enantiomeric cross-inhibition on an achiral poly(C) template: the opposite hand terminates the chain, so the feedback on the excess is negative and there is no autocatalytic k2 to compare to the floor. The row's own R1 = N and "inhibition only" comment say the same thing. By that rule the implied letter is **N**, not Y.
- Caveats: Kervio's conditions (0.4 M NaCl, 80 mM Mg, 20 degC, primer pre-annealed at 36 uM) are not prebiotic; its k' is for a defined primer:template, not Joyce's poly(C)/(G)n system; and the N rests on mechanism, not on a number. If the parent wants the cell decided only by measured magnitude, the honest letter is Y-magnitude/N-sign, i.e. leave ? and add a note.

### 2e. Update 2026-10-21: Joyce 1984 read, R5 = Y

`joyce_r5.py` (predictions `joyce_r5_predictions.txt`, sha256 2610b2a1...ac1f, written before the paper was opened; 5 of 6 held, P2 missed because one run is at 0.01 M) reads Table 1 (p. 602; 0 degC, 7 d, poly(C)): k2 = -ln(1-f)/(c t) with f = fraction in oligomers of length >= 4, so k2 is a lower bound. All-D >= 1.7e-5, racemic >= 1.7e-6 /M/s (8.4e-7, `joyce_r5.log:10`, if the DL concentration is per enantiomer, which the paper does not state). Both clear the scorecard floor of 3e-7.

Why 2d's N was not followed: the N clause of the R5 rule (scorecard.py comment) applies where "a model at measured rates gives an excess that decays", and no model has been run at Joyce's rates. R5 scores whether the step is fast enough; the missing amplification is already scored by R1 = N and R9 = N, so scoring it again in R5 would count it twice. Against the paper's own pond value (1 mM, 3e-6, ija_draft.md Section 5) the racemic bound falls short and all-D clears it.

## 3. Paywalled or unread, with locations

| reference | where the number should be |
|---|---|
| Leman, Orgel, Ghadiri 2004 Science 306:283, doi:10.1126/science.1102722 | read in full 2026-10-21 (copy from C. McKay); Table 1 end points and notes 5, 8, 9 and 14, see `aped_r5.py`; SOM (`leman.som.pdf`, in `papers/leman_2004_som/`, untracked) read 2026-10-21: dipeptide time courses in Figs. S3 and S4, no NCA rate, see `journals/new_analyses_aped_r5.md` 2b |
| Danger et al. 2006 JACS 128:7412, doi:10.1021/ja061339+ | read 2026-10-21 (copy from C. McKay, received 2026-10-20); no NCA formation rate or half-life, only end points and a CAA decomposition curve (Fig. 1). `aped_r5.py` uses Table 1 entries 2 (CAA 20 mM alone, 80 C, 28 d, f 0.017) and 3 (CAA 20 mM + Ala 100 mM, 80 C, 53 d, f 0.0275): a_lo 7.088e-9 and 6.090e-9 /s (`aped_r5.log:11-12`). The SI is unread |
| Joyce et al. 1984 Nature 310:602, doi:10.1038/310602a0 | read 2026-10-21 (copy from C. McKay); Table 1, see 2e |
| Collet et al. 1996 Tetrahedron Lett 37:9043, doi:10.1016/S0040-4039(96)01974-0 | NCA synthesis conditions |

Not retrievable without bypassing a bot check or paywall (not attempted): the HAL thesis (Anubis page), Wiley EurJOC 2025 and RSIF 2024 (HTTP 403), Springer Bomba 2019 (bot page), ETH Research Collection (HTTP 500). PMC11077012 returned a 151-byte error from Europe PMC.

## 4. Follow-ups, not done

1. Done 2026-10-21 for Leman 2004 and its SOM (`aped_r5.py`; SOM in `journals/new_analyses_aped_r5.md` 2b). The COS level is at most 0.09 mol% (note 14). This did not decide APED R5: it stays ? until a measured formation rate exists.
2. Done 2026-10-21 for Joyce 1984 (2e). Inoue & Orgel 1982/1983 not read.
