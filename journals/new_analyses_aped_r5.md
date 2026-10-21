# aped_r5: can the scorecard R5 "?" for APED peptides via NCAs be settled from published data?

Unit `aped_r5`, 2026-10-21. Target: `scorecard.py:62` (R5 "?"). Result: **R5 stays "?"**. No source prints an NCA formation rate constant. End-point data give a bracket, and that bracket straddles the pond need. A model check added during the unit shows a second reason the mapping used so far is incomplete (section 4). The Leman SOM, read after the unit, adds dipeptide time courses but no NCA rate (section 2b, post hoc).

Not blind, and said plainly: the v1 predictions file was written after I had read Leman 2004 Table 1 and Danger 2006 Table 1, so P2-P7 are arithmetic checks of a fixed rule. The real tests were P3, P4/P8, P6, P9, P10 and P11 (header of `aped_r5_predictions.txt`).

## 1. What ran

| step | file | pre-registration (sha256, written 2026-10-21 PDT) | log |
|---|---|---|---|
| rates to letters, model check at x* | `aped_r5.py` | `aped_r5_predictions.txt` 2fb0f2759ac8c188fc5319dd42eaecb515d3d7ad01210b6446886aa7a7b1c715 (00:31) | `aped_r5.log` |
| activation slowed at free drive x | `aped_r5_scan.py` | `aped_r5_predictions_v2.txt` 06e8573cba32dd0e5c500529a5c5cbbdf832adfa4374f77be47489c953b5cb00 (00:36) | `aped_r5_scan.log` |
| sigma(a) at fixed large x | `aped_r5_scan2.py` | `aped_r5_predictions_v3.txt` bef0d7c99e480bbcdfc6aec492c2a809415a84c8f99228ee54518a19de4d503d (00:40) | `aped_r5_scan2.log` |

Each script asserts its predictions hash and reproduces `aped.log` L89 (sigma = 0.0024332 /s at x* = 0.3852 M). No file existed before its hash; no smoke runs. v2 was written because P10 failed; v3 because v2's design was flawed (section 3). No predictions file was edited after a run. Compute: about 1 s, 2 s, 14 s.

## 2. Literature result (rate to letter)

Rule (`scorecard.py` L8, L37-43; `aped.py`): fast enough if sigma = s* a >= 10/tau, s* = 0.0063166, giving a_need = 5.016e-7 (pond), 5.016e-9 (lake), 5.016e-13 /s (ocean) (`aped_r5.log:1`).

Sources actually used for numbers:

- Leman, Orgel & Ghadiri, Science 306, 283 (2004), doi:10.1126/science.1102722, Table 1 (entries 1, 5-9) and notes 5, 8, 9, 14.
- Danger, Boiteau, Cottet & Pascal, JACS 128, 7412 (2006), doi:10.1021/ja061339+, Table 1 entries 2, 3.
- Devaraj et al., JACS 2016, doi:10.1021/jacs.6b08801, Fig. S7 (NCA hydrolysis kH = 0.018 /s, PMC7547885).
- Liu, Beaufils, Rossi & Pascal, Sci. Rep. 4, 7440 (2014), doi:10.1038/srep07440 (half-lives, read, not used for a rate).

No source prints a rate constant or time course for NCA formation from a free amino acid with COS, cyanate or NOx in water (P1 held). The Leman SOM, read after this unit, does not change that: its time courses are of dipeptide, from preformed thiocarbamate (section 2b). Only end points and dipeptide time courses exist, so a is bracketed.

| entry | a_lo (/s) | letters (pond, lake, ocean) | log line |
|---|---|---|---|
| L1 Phe 50 mM + COS, 25 C, final pH 8.7, 56 h | 1.716e-07 (molecule basis 3.493e-07) | ? Y Y (same on both bases) | `aped_r5.log:5`, `:14` |
| L1 upper bound, thiocarbamate loss ln2/10 h | a_hi 1.925e-05 (38.4x pond need) | none can be N | `aped_r5.log:15` |
| L5-L8 Pb2+, Pb2+ in Pacific water, Cd2+ | 7.304e-06, 1.565e-05, 1.102e-05, 3.961e-06 | Y Y Y | `aped_r5.log:6-9` |
| L9 K3Fe(CN)6 100 mM, 20 min | 1.669e-04 | Y Y Y | `aped_r5.log:10` |
| D2, D3 cyanate, 80 C | 7.088e-09, 6.090e-09 | ? Y Y | `aped_r5.log:11-12` |
| L1 scaled to volcanic COS (x 9.2e-4 and x 9.0e-4, unmeasured proportionality) | 1.58e-10, 1.54e-10 | N ? Y | `aped_r5.log:16-17` |

Row letter under the existing rule: **?** (the pond need 5.0e-7 lies inside 1.7e-7 to 1.9e-5), P8 held (`aped_r5.log:34`).

Estimator: a_lo = ln(1/(1-f_min))/t with f_min the fraction of residues that sit in a peptide bond (Leman: y2/2 + 2 y3/3; Danger: (y2 + 2 y3)[CAA]/[Ala]). It is a lower bound because every bond used one NCA. Metal and oxidant entries are reported but cannot change the row letter: those additives are not shown to be prebiotic, and the metal runs formed sulfide precipitates.

## 2b. Leman SOM (read 2026-10-21, after the unit; POST HOC, not pre-registered)

Source: `leman.som.pdf`, the Science supporting material for doi:10.1126/science.1102722; 7 pages, sha256 `613706d44e20b9f057577629a0065ed755b310827b82b9bd859183ffc202c98c`. The author downloaded it; it is kept in `papers/leman_2004_som/` (untracked) and was read with `python3 -I` (pypdf). The figures are raster images (about 300 ppi). Figs. S3 and S4 were digitized from them and each point was checked on an overlay. Marker identities were checked against the rendered captions, because the text layer drops the symbols. Plot readings carry about +-0.05 mM.

What it contains:

- No NCA formation rate constant and no NCA time course.
- Fig. S1 (1H NMR, 25 mM Phe, 500 mM borate, pD 8.9): 10 min of slow COS bubbling "cleanly converted" the amino acid to thiocarbamate. After 10 h, free amino acid is back beside the thiocarbamate (pD 8.6). So at the lab dose, thiocarbamate forms within 10 min, consistent with the main text calling cyclization to the NCA "likely" rate-limiting. Reading "cleanly" as >= 95% (my reading) gives k_f >= ln 20 / 600 s = 5.0e-3 /s at the lab dose.
- Fig. S2: dipeptide yield against pH at 52 h, open to air (one time point).
- Figs. S3 and S4: dipeptide concentration against time, from 25 mM Phe + 25 mM preformed thiocarbamate, with no COS gas. S3: 300 mM CHES pH 9.2, with 25 mM CdCl2, PbCl2, FeCl2 or FeCl3, or under argon. S4: 400 mM CHES pH 9.3 under argon, split at 40 h with one half open to air. Both start from 0 mM at t = 0.
- Methods: runs from free amino acid used 0.1 mmol in 2 ml (50 mM). HPLC quantitation against ABA is stated accurate within 15%.

**Yield basis.** Table 1 yields equal [FF]/25 mM, which is the fraction of all residues in dipeptide: PbCl2 at 27 h, 5.37 mM gives 21.5% (Table 1 entry 12: 21.6%); CdCl2 at 16 h, 10.69 mM gives 42.8% (entry 13: 48.0%, inside the 15%); FeCl3 at 42 h, 3.07 mM gives 12.3% (entry 15: 13.6% at 53 h). [FF]/50 mM would give 10.7%, 21.4% and 6.1%. So the residue basis (f = y2/2) is the right one for the 25+25 entries. Entry 1 (plain COS) is not plotted, so its basis is inferred, not shown.

**Argon time courses.** Same estimator as section 2, a_lo = ln(1/(1-f))/t with f = [FF]/50 mM (bonds per residue):

| run | t (h) | [FF] (mM) | yield [FF]/25 mM (%) | a_lo (/s) | a_lo / pond need |
|---|---|---|---|---|---|
| S4 argon | 1.75 | 0.387 | 1.5 | 1.233e-06 | 2.46 |
| S4 argon | 4.76 | 0.565 | 2.3 | 6.632e-07 | 1.32 |
| S4 argon | 15.97 | 0.775 | 3.1 | 2.717e-07 | 0.54 |
| S4 argon | 24.02 | 0.877 | 3.5 | 2.046e-07 | 0.41 |
| S4 argon | 38.96 | 1.202 | 4.8 | 1.735e-07 | 0.35 |
| S4 argon | 47.95 | 1.167 | 4.7 | 1.368e-07 | 0.27 |
| S4 argon | 63.98 | 1.475 | 5.9 | 1.300e-07 | 0.26 |
| S3 argon | 15.02 | 0.483 | 1.9 | 1.795e-07 | 0.36 |
| S3 argon | 24.02 | 0.861 | 3.4 | 2.009e-07 | 0.40 |
| S3 argon | 38.00 | 1.513 | 6.1 | 2.246e-07 | 0.45 |
| S3 argon | 47.98 | 2.081 | 8.3 | 2.461e-07 | 0.49 |

Only half the residues start as thiocarbamate, so per thiocarbamate the bound is about twice these (2.476e-06 /s at 1.75 h).

What this does and does not change:

1. At the lab dose, the best lower bound rises from 1.716e-07 (L1) to 1.233e-06 /s (S4, 1.75 h), 2.46x the pond need; 2.09x if the HPLC reads 15% high. Only the first two S4 points pass the pond need (4.76 h: 1.12x after the same 15% cut). Every point from 15 h on is 1.3e-07 to 2.7e-07 /s, like L1.
2. It is a different experiment from L1: preformed thiocarbamate and no COS. It bounds the cyclization-and-coupling stage, skipping a step that S1 shows takes at most 10 min at the lab dose.
3. The early rate rests on two points from one run. S3's argon run has no point between 0 and 15 h, and at 15-16 h it reads 0.48 mM against S4's 0.78 mM. S3's argon curve is labeled Table 1 entry 3, but it gives 8.3% at 48 h against entry 3's 5.6% at 60 h, and the buffers differ (S3 caption: 300 mM CHES pH 9.2; entry 3 footnote: 300 mM Me3N, pH 9.4). S4 at 64 h gives 5.9%, close to entry 3.
4. Dose. Scaled by the section 5 proportionality (x 9.0e-4 to 9.2e-4), the S4 bound becomes 1.11e-09 to 1.13e-09 /s, still below the lake need (5.0e-09). An alternative, also unmeasured: if only thiocarbamate formation scales with COS (k_f x dose) and cyclization does not, the steady thiocarbamate fraction against the 10 h loss is >= 0.189 to 0.193, giving a >= 4.7e-07 to 4.8e-07 /s (0.93x to 0.95x the pond need). The two dose models differ by about 400x.
5. Model window. Every SOM lower bound is below a_c1 in all three configurations of section 4 (A 2.83e-05, B 1.88e-04, C 7.46e-06 /s). The upper bound a_hi = 1.925e-05 is unchanged.

R5 stays "?". No letter was changed, no script was rerun, and scorecard.py and aped_r5.py are untouched. The SOM rows are arithmetic done after reading the figures; they need a hashed predictions file before they enter `aped_r5.py`.

## 3. Predictions: held and failed

- v1 (`aped_r5.log:27-38`): **10 of 11 held.** P1-P9 and P11 held. **P10 FAILED**: sigma/a at x* for kR = 0 was 6.6e-3, 1.9e-2, 1.3e-2, -7.1e-2, -0.91, -9.4 over a/kH = 0.385 down to 3.85e-6, against the registered band [2.1e-3, 1.9e-2] (`aped_r5.log:36`). The mapping sigma = s* a does not hold at slow activation at x*; sigma levels off near -3.6e-5 instead of going to zero (`aped_r5.log:20-25`).
- v2 (`aped_r5_scan.log:35-41`): **1 of 6 held** (V1 only). **V2-V6 FAILED** with "no boundary". The failure is a design flaw of mine, not a result about the model. Under the model's reading B, lowering kA at free x just moves the optimum to x ~ 1/kA, so a = kA x stays at about 0.38 kH (`aped_r5_scan.log:4-8`), and the grid ends at x = 1e4 M, where kA = 1e-5 reaches only a/kH = 0.1 (`aped_r5_scan.log:9`). Small a/kH was never tested. I report V2-V6 as failed.
- v3 (`aped_r5_scan2.log:65-71`): **6 of 6 held** (W1-W6), the slice that does test small a (next section).

## 4. Model finding: the instability window needs a/kH above a threshold

Slice: x = 1e4 M, kH = 1, kA = a/x (activation effectively irreversible); reverse activation is 1e-7 a. Leading difference-mode Re(lambda) against a/kH, 71 log points from 1e-7 to 1:

| config | a_c1/kH (sigma turns positive) | peak sigma, at a/kH | log line |
|---|---|---|---|
| A: kR = 0, kP = 1 | 1.572e-03 | +2.5545e-03, 0.398 | `aped_r5_scan2.log:19` |
| B: kR = 1e-4, kP = 1 (the aped.py base) | 1.044e-02 | +2.4357e-03, 0.398 | `aped_r5_scan2.log:37` |
| C: kR = 0, kP = 22 (real coupling is >= 0.4 /M/s, so kP/kH >= ~22 /M; Leman note 9) | 4.146e-04 | +2.5008e-02, 1.0 (grid edge) | `aped_r5_scan2.log:55` |

Below a_c1 the difference mode is stable at every a (A: -3.61e-05 at a/kH = 1e-7, `aped_r5_scan2.log:4`; B: -2.35e-04, `:22`), so the pool is not amplified however long tau is. The upper closing a_c2 was not reached (the grid stops at a/kH = 1, where sigma is still positive).

Against the lab values with kH = 0.018 /s (a_c1 in /s = log value x 0.018; arithmetic, not a log line): A 2.83e-05, B 1.88e-04, C 7.46e-06 /s. Plain COS is bracketed at 1.7e-07 to 1.9e-05 /s. So:

- With A and B, plain COS sits below the window: sigma < 0 at both lab ends (W5 held; `aped_r5_scan2.log:58-59` and `:61-62`).
- The ferricyanide lower bound (a/kH = 9.27e-3) is inside the window for A and outside it for B (W6 held; `aped_r5_scan2.log:60`, `:63`).
- With C the plain-COS bracket straddles a_c1, so it decides nothing.

Position of a_c1 depends on unmeasured ratios (kR alone moves it by 6.6x, A to B). The slice agrees with the x* slice, which brackets the zero between a/kH = 3.85e-4 and 3.85e-3 for A and 3.85e-3 and 3.85e-2 for B (`aped_r5.log:21-23`).

Consequence: a measured a is necessary for R5 but not sufficient. The existing floor (a >= 10/(tau s*)) uses the time rescaling sigma = s* a, which P10 and the slice above show is valid only near a/kH ~ 0.4 and not at the lab values for plain COS. Whether the floor is the binding condition or the window threshold is, at plain-COS rates, cannot be told without measured kH, kP, kD, kE, kR in the same conditions.

## 5. Other caveats

- Danger's cyanate route is at 80 C with no activation energy; its bracket has no upper bound.
- The lab COS dose (~400 mM nominal, ~25 mM dissolved) is not a prebiotic level. Volcanic COS is at most 0.09 mol% (Leman note 14); the proportionality of a to dissolved COS is an assumption, not measured.
- The yield basis in Leman Table 1 is unstated in the paper (factor 2, both printed). The SOM plots show it is the residue basis for the 25+25 entries (section 2b); L1 keeps both values because entry 1 is not plotted.
- The unstable window in aped.py needs ctot >= ~1 M, while the pond, lake and ocean settings are 1e-2 to 1e-6 M. This was already in the scorecard comment; the letters here follow its rule and do not remove that gap.
- Upper bound a_hi holds only if the thiocarbamate is not regenerated from COS.

## 6. Searches made (OpenAlex, Europe PMC, direct OA downloads)

OpenAlex works search (polite-pool email), with and without `open_access.is_oa:true`:
"carbonyl sulfide peptide formation water amino acid N-carboxyanhydride"; "N-carboxyanhydride prebiotic peptide formation aqueous kinetics"; "N-carbamoyl amino acid nitrosation N-carboxyanhydride peptide prebiotic"; "Peptide bond formation in water mediated by carbonyl sulfide Leman Huang Ghadiri"; "Pathways for the formation and evolution of peptides in prebiotic environments Danger Plasson Pascal"; "Amino acid N-carboxyanhydrides activated peptide monomers behaving as phosphate-activated nucleotides Biron Pascal"; "Carbonyl sulfide-mediated prebiotic formation of peptides Leman Orgel Ghadiri"; "Peptide formation mediated by cyanate revisited N-carboxyanhydrides accessible intermediates"; "Hitz Luisi spontaneous onset of homochirality oligopeptide N-carboxyanhydride water"; "Huber Wachtershauser peptides by activation of amino acids with CO on NiFeS surfaces"; "kinetics hydrolysis aminolysis N-carboxyanhydride amino acid water pH rate constant"; "carbonyl sulfide amino acid N-carboxyanhydride kinetics rate constant formation thiocarbamate aqueous"; "prebiotic carbonyl sulfide concentration volcanic atmosphere ocean peptide"; "N-carboxyanhydride formation rate aqueous amino acid carbonyl sulfide kinetics"; "cyanate amino acid peptide N-carboxyanhydride prebiotic aqueous"; "N-carboxyanhydride hydrolysis rate constant aqueous pH"; "amino acid N-carboxyanhydride prebiotic activation review Danger Plasson Pascal"; "carbonyl sulfide amino acid thiocarbamate cyclization N-carboxyanhydride rate"; plus title and author filters for Leman (carbonyl sulfide, 2015). Europe PMC DOI lookups for 10.1002/ejoc.202500786, 10.1007/s11084-019-09586-5, 10.1002/chem.71684, 10.26434/chemrxiv.11358596, 10.1098/rsif.2024.0014, 10.1038/s41467-022-30191-0, 10.1038/s42004-026-02149-6; OpenAlex location lookups for the first, third and fourth.

Downloaded (open access, `papers/aped_r5/oa/`, untracked): Pascal BIO Web Conf. 2014 (no rates), Liu 2014 Sci. Rep. 7440, Nat. Commun. 2023 coacervate paper, Commun. Chem. 2026 (irrelevant).

Blocked, not bypassed, not used: Desfoux et al. 2025 EurJOC (Wiley, and the HAL copy behind an Anubis bot challenge); ChemRxiv 11358596 (403, Cloudflare); Springer OLEB papers (Bomba 2019, Hitz & Luisi 2004); Royal Society RSIF 2024. Leman 2004 and Danger 2006 were sent by C. McKay on request on 2026-10-20 (not from Dworkin, not obtained around a paywall); the Science and JACS SI were not received. The Science SOM was later downloaded by the author on 2026-10-21 (science.org had shown me a Cloudflare check, which I did not bypass) and read (section 2b). The JACS SI is still unread. No Dworkin email data was used.

## 7. Proposed wording (NOT APPLIED; scorecard.py, README.md, journals/ija_draft.md untouched)

- `scorecard.py:62` comment, replace "R5 ? (NCA activation rate unmeasured)" with: "R5 ? (NCA activation rate only bracketed: plain-COS end points give 1.7e-7 to 1.9e-5 /s against a pond need of 5.0e-7; lake and ocean pass at the lab COS dose; in the aped.py model the difference mode is also stable below a/kH ~ 2e-3 to 1e-2, so sigma = s* a is not valid there; aped_r5)". Letter stays "?".
- README / paper sentence: "A measured NCA formation rate would settle R5 only together with the hydrolysis, coupling and racemization rates under the same conditions; the published end points bracket the activation rate at 1.7e-7 to 1.9e-5 /s (plain COS, 25 C), which straddles the pond requirement of 5.0e-7 /s." Add: "In the model, amplification needs the activation rate above roughly 2e-3 to 1e-2 of the NCA hydrolysis rate (0.018 /s), so the floor from time scale alone understates the requirement."
- If the SOM rows (section 2b) are adopted after a hashed run, the lab-dose lower bound in the two items above becomes 1.2e-6 /s (preformed thiocarbamate), and "straddles the pond requirement" no longer holds at the lab dose.
- If a COS-dose-scaled estimate is quoted: "assuming (unmeasured) proportionality to dissolved COS, volcanic COS lowers the plain-COS lower bound to 1.6e-10 /s (lake ?, ocean Y)".

## 8. Found outside the unit (nothing edited)

- `aped.py:120-138` (docstring POST-HOC paragraphs) quotes Leman/Danger numbers from the McKay copies; fine as disclosed, but R5 wording there ("every measured uncatalysed a fails the pond and passes the lake and ocean") now also needs the window caveat of section 4.
- `aped.py:152` and `:356` onward: the script runs its full ~86 s scan at import (module-level `T0 = time.time()` and prints, no `if __name__` guard anywhere in the 510 lines), so it cannot be imported; `aped_r5.py` copies the functions for that reason.
- `scorecard.py:62`: R5 comment wording, see section 7.
