# DNA and RNA Chirality, Parity Violation, and Prebiotic Coupling: Literature Survey

**Survey Date:** 2026-10-10  
**Status:** References checked against Crossref on 2026-10-10. One-line findings are paraphrased from titles/abstracts unless noted, not checked against full texts.  
**Context:** Notes on the handedness of DNA/RNA (D-sugars, right-handed B-DNA helix) in relation to the electroweak parity-violating energy difference (PVED), cross-catalysis with L-amino acids, and prebiotic reaction networks.

---

## 1. Parity Violation in Nucleic Acids & Sugars

| Reference | System / Focus | Findings | Verification Status |
| :--- | :--- | :--- | :--- |
| **Faglioni, D'Agostino, Cadioli & Lazzeretti (2005)**, *Parity violation energy of biomolecules – II: DNA*, *Chem. Phys. Lett.* 407:522–526 | B-DNA double helix | *Ab initio* PVED calculations on DNA models. Per its abstract, the weak interaction does not favor the natural right-handed double helix over its mirror image. | Crossref verified; DOI: 10.1016/j.cplett.2005.04.009 |
| **Faglioni, García Cuesta & Lazzeretti (2006)**, *Parity violation energy of biomolecules – III: RNA*, *Chem. Phys. Lett.* 432:263–268 | RNA | *Ab initio* PVED calculations on RNA models; reported as not favoring the natural form. | Crossref verified; DOI: 10.1016/j.cplett.2006.10.030 |
| **Tranter, MacDermott, Overill & Speers (1992)**, *Computational studies of the electroweak origin of biomolecular handedness in natural sugars*, *Proc. R. Soc. Lond. A* 436:603–615 | Ribose / deoxyribose conformers | PVED sign and size depend on the sugar ring conformation, which is the root of the "conformation dependence" problem for sugars. The specific conformer-by-conformer signs previously listed here have not been checked against the paper. | Crossref verified; DOI: 10.1098/rspa.1992.0037 |
| **Quack, Seyfang & Wichmann (2022)**, *Perspectives on parity violation in chiral molecules…*, *Chem. Sci.* 13:10598–10643 | Molecular parity violation review | Review of PVED theory, spectroscopic experiments toward measuring it, and links to biomolecular homochirality. Specific PVED magnitudes previously quoted here were not checked against the review and are removed. | Crossref verified; DOI: 10.1039/d2sc01323a. (A "2026 *Annu. Rev. Phys. Chem.*" review previously cited could not be verified and is removed.) |

---

## 2. Prebiotic Sugar–Amino Acid Coupling & Network Homochirality

| Reference | Mechanism | Relevance to the Hand of DNA |
| :--- | :--- | :--- |
| **Ozturk & Sasselov (2025)**, *Life's homochirality: Across a prebiotic network*, *PNAS* 122(34):e2505126122, doi:10.1073/pnas.2505126122 | Perspective / analysis (not an experiment) | Argues for a unidirectional chiral flow from nucleic acids to coded peptides across a prebiotic network. |
| **Higgs & Blackmond (2025)**, *Autocatalytic symmetry breaking and chiral amplification in a feedback network combining amino acid synthesis and ligation*, *PNAS* 122(20):e2423683122, doi:10.1073/pnas.2423683122 | Model of amino acid synthesis + peptide ligation feedback | The forward reaction alone gives no ee amplification (homochiral dimers catalyse the opposite-hand monomer); closed systems return to racemic, while open flow can keep broken symmetry. |
| **Breslow & Cheng (2010)**, *PNAS* 107:5723–5725, doi:10.1073/pnas.1001639107 | L-amino acids → D-sugar excess | L-amino acid monomers (not dipeptides) catalysed formation of D-glyceraldehyde in excess (up to about 60:40) from formaldehyde and glycolaldehyde. |
| **Tamura & Schimmel (2004)**, *Chiral-selective aminoacylation of an RNA minihelix*, *Science* 305:1253, doi:10.1126/science.1099141 | D-RNA → L-amino acid aminoacylation | RNA minihelix aminoacylation reported to prefer L- over D-amino acids (about 4:1; ratio carried over from the earlier draft, not re-checked). |

Not included: a "Kim et al. 2025" result that mirror-image RNA prefers D-amino acids could not be found and is treated as unverified.

---

## 3. The Opposing-PVED Question (`dna_network.py`, model only)

A standing question in the PVED literature:
> *If the weak force favors L-amino acids ($g_a > 0$) but opposes the natural D-sugar / right-handed helix ($g_s < 0$), don't the two biases cancel?*

The project's model in `dna_network.py` gives one possible answer. It is a model result, checked only for internal self-consistency, not against experiment:
1. **Asymmetric cross-catalysis:** the growing pair mode $u$ has left eigenvector $w = (\sqrt{\kappa_{as}}, \sqrt{\kappa_{sa}})^T$.
2. **Projected bias:**
   $$g_{\text{eff}} = \frac{\sqrt{\kappa_{as}} g_a + \sqrt{\kappa_{sa}} g_s}{\sqrt{\kappa_{as} + \kappa_{sa}}}$$
3. **Condition in the model:**
   When $\kappa_{as} / \kappa_{sa} > (|g_s| / g_a)^2$, $g_{\text{eff}} > 0$ and the natural pair (L-amino acid + D-sugar) is favored in the model despite the opposing sugar PVED.

---

## 4. Synthetic Mirror-Image Molecular Biology (Ting Zhu lab)

*   **Mirror Pfu DNA polymerase:** chemically synthesized D-amino-acid Pfu polymerase used to amplify and store information in L-DNA. *Fan, Deng & Zhu (2021), Nat. Biotechnol. 39:1548–1555, doi:10.1038/s41587-021-00969-6.*
*   **Mirror transcription:** mirror-image T7 RNA polymerase transcribing L-RNAs, including ribosomal RNAs. *Xu & Zhu (2022), Science 378:405–412, doi:10.1126/science.abm0646.*
*   **Mirror protein sequencing:** *Zhang & Zhu (2024), Nat. Chem. 16:592–598, doi:10.1038/s41557-023-01411-x.*
*   **Commentary:** *Zhu (2025), Nature 645:588–591, doi:10.1038/d41586-025-02912-0.*
*   **Implication:** mirror enzymes show comparable activity to their natural counterparts in these demonstrations. This suggests that the D-DNA / L-protein choice is not required for function, which is consistent with (but does not prove) the idea that the hand was set by an early symmetry-breaking event and then locked in by cross-catalysis.
