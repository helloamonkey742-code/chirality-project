# Full-text check A (H7): Kondepudi & Nelson, Yamagata, Cowan & Furnstahl, Viedma

Checked 2026-10-13. Only legitimate open copies were used. Every DOI below was confirmed on Crossref, and each title matches the draft's reference list.
Line numbers refer to `ija_draft.md` as it stands today. The draft has not been edited.

---

## 1. Kondepudi & Nelson 1983, 1984a, 1984b, 1985 (and Kondepudi 1987)

**References**
- Kondepudi & Nelson 1983, PRL 50, 1023–1026, doi 10.1103/physrevlett.50.1023
- Kondepudi & Nelson 1984a, Physica A 125, 465–496, doi 10.1016/0378-4371(84)90065-7
- Kondepudi & Nelson 1984b, Phys Lett A 106, 203–206, doi 10.1016/0375-9601(84)90320-7
- Kondepudi & Nelson 1985, Nature 314, 438–441, doi 10.1038/314438a0
- Kondepudi 1987, Biosystems 20, 75–83, doi 10.1016/0303-2647(87)90022-0

**Source read: ABSTRACT ONLY.** All five papers are closed.
- OpenAlex `is_oa` is false for each, and none has an arXiv, PMC or author copy.
- The 1985 Nature abstract was read at https://www.nature.com/articles/314438a0.
- No full text of any K&N paper was read, so their noise treatment and volume estimate **could not be checked first-hand**.

**What can be confirmed**
- *From the 1985 abstract:*
  - The mechanism is simple and extremely sensitive.
  - Systems spontaneously evolve to an L- or D-dominated state.
  - A systematic chiral interaction no stronger than the weak-neutral-current energy in amino acids can decide which enantiomer dominates in about 15,000 years.
  - This agrees with Mason & Tranter's finding that L is favoured.
  - The abstract gives **no volume, concentration or probability figure**.
- *From open secondary sources (full text):*
  - **Quack, Seyfang & Wichihase 2022**, Chem Sci, doi 10.1039/d2sc01323a, PMC9491092 (already cited in the draft). It says K&N proposed that selection could happen in a lake with a plausible "chemical soup" over roughly 10⁴–10⁵ years (its refs 359, 360, 362). It adds that the autocatalytic mechanism is formal mathematics, with no attempt to simulate early-Earth conditions.
  - **Cowan & Furnstahl 2022**, arXiv 2110.12232, Introduction. It says K&N gave time frames of up to about 100,000 years, depending on reaction rates, concentrations, supply, reversibility and diffusion.
  - **Gleiser & Thorarinson 2006** (arXiv astro-ph/0601399) and **Gleiser, Thorarinson & Walker 2007** (astro-ph/0606593). Both say K&N wrote the bias as k_L/k_R ≈ 1 + g, with g = ΔE/kT ~ 10⁻¹⁷–10⁻¹⁸. This supports the draft's L313.
  - **Gleiser 2022** (arXiv 2205.01193). It describes K&N as a rate-equation model in which a small excess is amplified, with parity violation acting as the bias and with fluctuations included. It notes that their well-mixed assumption was later questioned by spatiotemporal analyses.
- **Not confirmed:** I half-remember a specific lake size and selection probability from the 1985 paper, but no open source confirms it. **Do not use any such number.**

**Draft sentences and proposed fixes**

**L279 (part 1).** Keep it. The abstract supports it.
> "In systems that evolve to an L- or D-dominated state, a systematic chiral interaction no stronger than the weak-neutral-current energy in amino acids can determine which enantiomer dominates over about 15,000 years (Kondepudi and Nelson, 1985), in line with Mason and Tranter's finding for L-amino acids."

**L279 (last sentence).**
> "The details of their noise treatment and volume estimate are [CHECK FULL TEXT]."

Replace with:
> "Later accounts place their estimate in a lake-sized body of water over about 10⁴–10⁵ years (Quack et al., 2022); I could not obtain the original papers to check their noise treatment or the exact volume."

Alternatively, cut the sentence.

**L289.**
> "Kondepudi and Nelson (1985) already estimated the volume and time such a bias needs, and the "large quantities of reactants and a long reaction time" of Tranter (1985b) make the same point."

Replace with:
> "Kondepudi and Nelson (1985) already argued that such a bias needs a large system and a long time (a lake-sized body over about 10⁴–10⁵ years in later accounts; Quack et al., 2022), and the "large quantities of reactants and a long reaction time" of Tranter (1985b) make the same point."

**L721.**
> "The size rule restates the balance of Kondepudi and Nelson (1985), who already estimated the volume and time a weak-current bias needs, as a count of molecules in one mixed pool ..."

Replace the clause with:
> "...the balance of Kondepudi and Nelson (1985), who already argued that a weak-current bias needs a large system and a long time, as a count of molecules..."

**Note 15 (L214).**
> "...their value is [CHECK FULL TEXT] (H7). Quote it once read."

Replace with:
> "**Kondepudi and Nelson (1985) volume estimate.** The original papers are not openly available. Sections 1.1 and 4 therefore credit them only with the abstract's ~15,000-year time scale and the lake-sized system reported by Quack et al. (2022), and quote no volume."

**L313.** Keep it. The k(1 ± g) parameterisation is confirmed by the Gleiser papers.

**Confidence**
- High that the papers are closed and that the abstract says what L279 part 1 claims.
- Medium on "lake, 10⁴–10⁵ yr". It is secondary, but Quack and Cowan independently agree.
- **Low** on anything about their noise model or volume beyond that. Hence the hedge.

---

## 2. Yamagata 1966

**Reference.** Yamagata, Y. 1966, "A hypothesis for the asymmetric appearance of biomolecules on earth", J Theor Biol 11, 495–498, doi 10.1016/0022-5193(66)90110-x.

**Source read: NONE (title only).**
- The paper is closed: OpenAlex lists no OA location, and PubMed 5967446 has no abstract.
- The main later source, Wesendrup et al. 2003, is also closed.
- The open citing papers checked through Europe PMC full text name Yamagata without describing his mechanism or wording.
- One closed 1999 paper has the title "Chirality amplification – the accumulation principle revisited" (doi 10.1023/a:1006646021670, PMID 10666744). Its title shows that later authors call the idea an "accumulation principle".
- That is **later naming, not confirmed as Yamagata's own wording**.

**What the paper says.** It cannot be confirmed from an open source. Neither the mechanism (steps at which a tiny energy difference compounds) nor the words "accumulation" or "amplification" can be attributed to Yamagata himself.

**Draft sentence (L273).**
> "Yamagata (1966) put forward a hypothesis for the asymmetric appearance of biomolecules on Earth; later work calls the link between electroweak energy differences and homochirality "the Yamagata hypothesis" (Wesendrup et al., 2003) and summarises its core as a tiny bias that accumulates through repeated steps or amplification until one hand dominates [CHECK FULL TEXT for this wording]."

**Proposed replacement.**
> "Yamagata (1966) put forward a hypothesis for the asymmetric appearance of biomolecules on Earth; later work calls the link between electroweak energy differences and homochirality "the Yamagata hypothesis" (Wesendrup et al., 2003), and it is often described as an accumulation principle."

This cuts the "tiny bias that accumulates through repeated steps or amplification" summary. If even the last clause is unwanted, end the sentence after "(Wesendrup et al., 2003)".

**L197 (H7 list).** Mark it resolved: "cut — original closed; attribution limited to title and later naming."

**Confidence**
- High that nothing confirms the wording.
- Medium that "accumulation principle" is the common later label. That rests on one 1999 title, which should be checked if anyone has access to it.

---

## 3. Cowan & Furnstahl 2022 (Section 1.2)

**Reference.** Cowan, J.A., Furnstahl, R.J. 2022, "Origin of Chirality in the Molecules of Life", ACS Earth Space Chem 6, 2575–2581, doi 10.1021/acsearthspacechem.2c00032.

**Source read: FULL TEXT (author preprint),** arXiv 2110.12232 (https://arxiv.org/abs/2110.12232), including the supplementary equation.
- This is the preprint, not the typeset ACS version, so page numbers are not given.
- Table and equation numbers refer to the preprint.

**What the paper says**
- **The reaction.** They model the autocatalytic self-replication of glyceraldehyde in the formose reaction: formaldehyde plus glycolaldehyde gives glyceraldehyde, a precursor of ribose. The setting is an RNA-world framework, not amino acids.
- **Metal ions.** The aldol step is promoted by alkaline-earth ions.
- **Computed PVEDs (Table 1), for one and two catalytic ions:**
  - Ca²⁺: about 0.3–0.6×10⁻²⁰ au
  - Sr²⁺: about 0.6–1.2×10⁻¹⁹ au
  - Ba²⁺: about 0.3–0.6×10⁻¹⁸ au
- **How selectivity builds up.** Selectivity compounds deterministically over turnovers as exp(ΔΔE/RT)ⁿ (Eq. S1), with a rate constant of about 1 min⁻¹. Tables 2–3 cover time spans from 10⁴ years to Gyr.
- **Temperature.** They use T = 315 K, taken as an approximate average early-Earth ocean temperature. Lower and higher temperatures are discussed.
- **Conclusion.** An energy difference of at least about 10⁻¹⁹ au is needed. Ca²⁺ is marginal; heavier ions work over geological time.
- **No noise.** There is no fluctuation or noise term. They only mention a simple stochastic model and K&N's time frames in the Introduction.
- They suggest exoplanet chirality statistics as an eventual test.

**Draft sentence (L293).**
> "Cowan and Furnstahl (2022) combined the PVED with metal-catalysed self-amplifying synthesis in an early-Earth ocean [CHECK FULL TEXT]."

**Proposed replacement.** Keep the sentence, made more exact.
> "Cowan and Furnstahl (2022) computed the PVED for an alkaline-earth-ion-catalysed autocatalytic sugar (formose) step and compounded it over turnovers at an early-ocean temperature, without fluctuations."

This keeps the contrast with the draft's noise-limited model.

**Confidence.** High.

---

## 4. Viedma 2005

**Reference.** Viedma, C. 2005, PRL 94, 065504, doi 10.1103/PhysRevLett.94.065504.

**Source read: FULL TEXT (author preprint),** arXiv cond-mat/0407479v5 (https://arxiv.org/abs/cond-mat/0407479), dated 2005-03-13. Its journal reference is PRL 94 065504.
- The preprint's title differs from the PRL title.
- The UCM institutional PDF of the published version could not be downloaded.

**What the paper says**
- **Setup.**
  - 12 g of NaClO₃ crystals in 10 mL water, at 24 °C.
  - Magnetic stirring at about 600 rpm.
  - 2–8 g of 3 mm glass balls. Their continuous abrasion keeps the crystals small, and the paper says the system works like a mill.
  - Each flask holds more than a million crystals.
- **Results.**
  - Stirring without balls left the initial crystal excess unchanged.
  - With balls, which the paper calls the abrasion-grinding process, the excess rose to total chiral purity.
  - A 5% L or D excess went to 100% of the same hand in 8 h (Fig. 5).
  - Symmetric mixtures reached purity with a random hand (Fig. 6).
  - The rate depends on the mass of balls (Fig. 3) and the rpm (Fig. 4).
  - Results are averages over more than 20 repeats.
- **Mechanism.**
  - Grinding acts as induced secondary nucleation, so fragments keep their parent's hand.
  - Recycling is Ostwald-ripening dissolution of the smallest fragments back to achiral solute.

**Draft sentences**

**L283.**
> "Viedma (2005) showed that grinding-driven crystallization with recycling turns a symmetric population of D and L crystals into one of complete chiral purity."

"Grinding" is correct. A slightly more exact version:
> "Viedma (2005) showed that grinding crystals with glass beads in their saturated solution, with dissolution and regrowth as recycling, turns a symmetric population of D and L sodium chlorate crystals into one of complete chiral purity."

The current sentence is also acceptable.

**L710.**
> "A 5% crystal excess over more than 10⁶ crystals always set the outcome (Viedma, 2005)."

"Always" goes slightly beyond what is reported (averaged results). Replace with:
> "A 5% crystal excess among more than 10⁶ crystals went to 100% of the same hand within 8 h (Viedma, 2005)."

**L759 and EXPERIMENT.md L65.** "Viedma grinding positive control" is fine.

**L206 and history-section flag (`ija_history_section.md` L153).** Resolved: "grinding" is confirmed in full text. The paper's own term is "abrasion-grinding".

**Confidence.** High. The only caveat is that the source is the preprint, not the typeset PRL.
