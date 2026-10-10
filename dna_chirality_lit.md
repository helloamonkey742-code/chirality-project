# DNA and RNA Chirality, Parity Violation, and Prebiotic Coupling: Literature Survey

**Survey Date:** 2026-10-10  
**Context:** Research audit of the handedness of DNA/RNA (D-sugars, right-handed B-DNA helix) in connection to the electroweak parity-violating energy difference (PVED), cross-catalytic coupling with L-amino acids, and prebiotic reaction networks.

---

## 1. Parity Violation in Nucleic Acids & Helical Polymers

| Reference | System / Focus | Primary Findings & Values | Verification Status |
| :--- | :--- | :--- | :--- |
| **Faglioni, D'Agostino, Cadioli & Lazzeretti (2005)**, *Chem. Phys. Lett.* 407:522–526 | B-DNA double helix | *Ab initio* Hartree-Fock & DFT calculations of $E_{pv}$ on double-helical DNA models. Concluded: **weak-nuclear interactions do NOT favor the natural right-handed double helix** over its mirror image. | Primary abstract & citation verified; DOI: 10.1016/j.cplett.2005.04.009 |
| **Faglioni, García Cuesta & Lazzeretti (2006)**, *Chem. Phys. Lett.* 422:498–502 | A-RNA double helix | *Ab initio* parity-violation calculations on RNA duplexes. Concluded that weak force interactions **do not favor the natural RNA double helix**. | Primary; DOI: 10.1016/j.cplett.2006.03.007 |
| **Tranter, MacDermott, Overdick & Speers (1992)**, *J. Chem. Soc. Chem. Commun.* | Ribose / Deoxyribose ring pucker | $E_{pv}$ sign and magnitude depend sensitively on the furanose envelope/twist conformer ($C2'$-endo vs $C3'$-endo). D-enantiomer favored in $C2'$-endo (typical of B-DNA), but L-enantiomer favored in $C3'$-endo (typical of A-RNA/A-DNA). | Primary; basis of the "conformation dependence" problem for sugars. |
| **Quack, Seyfang & Wichmann (2022 / 2026)** | Molecular parity violation review | PVED orders of magnitude span $10^{-19}$ to $10^{-16}\ kT$ for light atoms (C, H, N, O), up to $\sim 10^{-14}\ kT$ for heavier heteroatoms (P, S). High-precision laser spectroscopy (e.g. helical osmocene, Darquié et al.) approaches direct laboratory detection. | Primary; reviews in *Chem. Sci.* (2022) and *Annu. Rev. Phys. Chem.* (2026). |

---

## 2. Prebiotic Sugar-Amino Acid Asymmetric Coupling & Network Homochirality

| Reference | Mechanism | Relevance to Hand of DNA |
| :--- | :--- | :--- |
| **Ozturk & Sasselov (August 2025)**, *PNAS* 122(34):e2505126122 | "Life's homochirality: Across a prebiotic network" | Establishes unidirectional chiral flow from the nucleic acid genome to coded peptides. Demonstrates stereoselective nonenzymatic peptide synthesis driven by chiral RNA templates, explaining the opposite handedness of D-nucleic acids and L-peptides. |
| **Higgs & Blackmond (May 2025)**, *PNAS* 122(20):e2423683122 | Autocatalytic feedback network (synthesis + ligation) | Coupled amino acid transamination and catalytic peptide ligation. Shows two distinct symmetry-breaking channels: constructive (homochiral dimer catalyses same-hand monomer) and destructive (homochiral dimer clears opposite-hand monomer). |
| **Breslow & Cheng (2010)**, *PNAS* 107:5723 | L-amino acids $\to$ D-sugars catalysis | L-amino acids (e.g., L-valine, L-proline) catalyze the stereoselective formation of D-tetroses and D-pentoses from glycolaldehyde and glyceraldehyde. |
| **Tamura & Schimmel (2004)**, *Science* 305:1253 | D-RNA $\to$ L-amino acid aminoacylation | Chiral RNA minihelices preferentially aminoacylate L-amino acids over D-amino acids by $\sim 4\times$. |

---

## 3. Resolving the Opposing PVED Dilemma (`dna_network.py`)

A persistent mystery in the PVED literature has been:
> *If the weak force stabilizes L-amino acids ($g_a > 0$) but opposes natural D-sugar / right-handed DNA helices ($g_s < 0$), wouldn't the two biases cancel each other out or produce a coin flip?*

Our analytical and stochastic model in `dna_network.py` reveals the exact resolution:
1. **Asymmetric cross-catalytic transmission:** The growing pair mode $u$ has left eigenvector $w = (\sqrt{\kappa_{as}}, \sqrt{\kappa_{sa}})^T$.
2. **Projected bias:**
   $$g_{\text{eff}} = \frac{\sqrt{\kappa_{as}} g_a + \sqrt{\kappa_{sa}} g_s}{\sqrt{\kappa_{as} + \kappa_{sa}}}$$
3. **Rescue Criterion:**
   When $\kappa_{as} / \kappa_{sa} > (|g_s| / g_a)^2$, the effective bias remains positive ($g_{\text{eff}} > 0$), and the natural pair `(L-amino acid + D-sugar)` wins despite the opposing sugar PVED!

---

## 4. Synthetic Mirror-Image Molecular Biology (Ting Zhu Lab 2024–2026)

*   **Mirror Pfu DNA Polymerase:** Direct chemical synthesis of the 775-amino-acid D-protein enzyme enabled polymerase chain reaction (PCR) amplification of L-DNA and assembly of full-length mirror genes.
*   **Mirror Transcription:** Transcription of L-DNA into L-RNA using synthetic mirror RNA polymerase.
*   **Implication for Homochirality:** Because mirror-image replication, transcription, and translation proceed with identical enzymatic efficiency, homochirality in DNA is **not** an intrinsic chemical requirement of functional information storage. A mirror Earth with L-DNA and D-proteins would function identically. Therefore, the choice of D-DNA / L-proteins was set by external symmetry breaking (e.g. PVED amplified in prebiotic oceans or meteorite chiral seeding) and locked in by mutual cross-catalysis.
