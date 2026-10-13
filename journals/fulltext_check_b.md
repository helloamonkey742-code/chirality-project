# Full-text check B: Section 1.1 [CHECK FULL TEXT] items

Checked 2026-10-13. Draft: `journals/ija_draft.md`, which was not edited.

**Access summary.** OpenAlex and Semantic Scholar list all seven target papers as closed, with no repository copy. ScienceDirect returned 403 and ACS returned a Cloudflare challenge; I did not bypass either. Requests to link.springer.com and pubs.rsc.org from this machine were redirected to google.com, which looks like a local network or DNS filter. IA Scholar requires a bot check, so I did not use it. Springer's free first-page preview server (page-one.springer.com) did work for Rein (1974) and Bonner (1991). The open secondary sources I read in full:

- **Bonner (1974) preprint.** NASA NTRS 19740027018, https://ntrs.nasa.gov/api/citations/19740027018/downloads/19740027018.pdf (24 pp.). It was presented at the 1973 Jülich symposium and published as Bonner, W.A. 1974, *J. Mol. Evol.* 4, 23–39, https://doi.org/10.1007/bf01732769 (Crossref checked: title, volume and pages match). I read the NASA preprint, not the journal typeset version.
- **Quack, Seyfang and Wichmann (2022), *Chem. Sci.* 13, 10598.** PMC9491092, full text via NCBI efetch. It is already in the reference list.
- **Darquié et al. (2010), *Chirality* 22, 870.** arXiv:1007.3352, full text. It is not in the reference list.

All DOIs below were checked with `api.crossref.org/works/<doi>`: title, volume and pages match the draft's entries.

---

## 1. Ulbricht and Vester (1962)

**Reference (draft line 1038).** Ulbricht, T.L.V., Vester, F. 1962, Attempts to induce optical activity with polarized β-radiation, Tetrahedron, 18, 629–637. https://doi.org/10.1016/s0040-4020(01)92714-0

**Source actually read.** The paper itself is closed. OpenAlex reports `oa_status: closed`, ScienceDirect returned 403, and there is no abstract in Crossref, OpenAlex or Semantic Scholar. What I read in full is Bonner's 1974 account (NASA NTRS preprint, PDF pp. 2–4), which summarizes this paper and its 1959 predecessor directly.

**What it says on the point** (Bonner 1974, paraphrased):

- **Mechanism (PDF pp. 2–3).** Vester first suggested the idea in a 1957 Yale seminar, and Ulbricht (1959, 1962) set it out as a three-step chain:
  1. β-decay gives longitudinally polarized electrons.
  2. As these slow down in matter, they emit circularly polarized γ- and X-ray bremsstrahlung.
  3. The circularly polarized photons act asymmetrically on organic matter.
- **Expected size (p. 3).** Ulbricht and Vester (1962) calculated that even under favourable assumptions the effect might not be measurable in a reasonable time. Their favourable case was a 1 Ci source, full absorption, a product with specific rotation of 100° and a quantum yield of 1; even then, measurable activity would take about a year.
- **Experiments and outcome (pp. 3–4).** Bonner describes the 1959 and 1962 work together:
  - About ten synthetic and degradative reactions were run with several β-emitters (including ³²P and ⁹⁰Sr–⁹⁰Y), at 30–2500 mCi and from about −50 to +50 °C.
  - Some 36 experiments were measured on two precision polarimeters.
  - Every rotation was zero within experimental error. The authors put any induced optical activity below 0.02% and said stronger sources and longer runs would be needed.

**Current draft sentence (line 273).**
> "and Ulbricht and Vester (1962) reported attempts to induce optical activity with polarized beta-radiation. Their mechanism (including any role of circularly polarized bremsstrahlung) and outcome are [CHECK FULL TEXT]."

**Proposed replacement** (this replaces the first sentence's last clause and the second sentence, and covers item 2 as well):
> "Vester, Ulbricht and Krauch (1959) and Ulbricht and Vester (1962) proposed and tested a beta-decay route: longitudinally polarized beta electrons slow down in matter and emit circularly polarized bremsstrahlung, which could act asymmetrically on organic molecules. Across about ten reactions with several beta emitters they found no optical activity above 0.02% (Bonner, 1974)."

Add to the reference list: Bonner, W.A. 1974, Experiments on the origin of molecular chirality by parity non-conservation during β-decay, Journal of Molecular Evolution, 4, 23–39. https://doi.org/10.1007/bf01732769

**Confidence.** Medium-high for the mechanism, the null result and the 0.02% limit. They come from a detailed contemporary account by the main later experimenter, but not from the 1962 paper itself. The emitter list and temperatures were partly garbled in the scan, so the replacement leaves them out.

---

## 2. Vester, Ulbricht and Krauch (1959)

**Reference (line 1044).** Vester, F., Ulbricht, T.L.V., Krauch, H. 1959, Optische Aktivität und die Paritätsverletzung im β-Zerfall, Die Naturwissenschaften, 46, 68. https://doi.org/10.1007/bf00599091

**Source actually read.** The paper itself is closed: OpenAlex reports it as closed, and the Springer preview returned 403. What I read in full is Bonner (1974), the NASA preprint, PDF pp. 2–4.

**What it says on the point.** Bonner credits the mechanism to Vester (1957 seminar) and lists this 1959 note as one of the two papers that tested it experimentally (pp. 3–4). The mechanism is the same chain as in item 1: β-decay, then polarized electrons, then circularly polarized bremsstrahlung, then asymmetric photochemistry. Bonner pools the 1959 and 1962 experiments, so I cannot tell which results appeared in the one-page 1959 note.

**Current draft sentence (line 273).**
> "Vester, Ulbricht and Krauch (1959) connected optical activity with parity violation in beta-decay,"

**Proposed replacement.** Merge it into the replacement in item 1 and drop the separate clause. If you want to keep it separate, use: "Vester, Ulbricht and Krauch (1959) proposed that parity violation in beta-decay could produce optical activity through circularly polarized bremsstrahlung." Do not credit the 1959 note alone with a result.

**Confidence.** Medium for the mechanism, since it comes from a secondary account. The present clause is already safe because it paraphrases the title.

---

## 3. Rein (1974) and Letokhov (1975)

**References.**
- Line 1004: Rein, D.W. 1974, Some remarks on parity violating effects of intramolecular interactions, J. Mol. Evol., 4, 15–22. https://doi.org/10.1007/bf01732768
- Line 954: Letokhov, V.S. 1975, On difference of energy levels of left and right molecules due to weak interactions, Phys. Lett. A, 53, 275–276. https://doi.org/10.1016/0375-9601(75)90064-x

**Sources actually read.**
- **Rein.** The publisher's free preview (https://page-one.springer.com/pdf/preview/10.1007/BF01732768), which covers pp. 15–16 in full: the author's Summary, the Introduction and the start of §2. Pages 17–22 were not available.
- **Letokhov.** The paper is closed and ScienceDirect returned 403. I read two later papers in full that describe it:
  - Darquié et al. (2010), arXiv:1007.3352, p. 2.
  - Quack et al. (2022), PMC9491092, §4 "Towards spectroscopic experiments".

**What they say.**
- **Rein (Summary, p. 15).** Parity-violating forces do distinguish mirror-image molecules, but the corrections to binding energy are tiny. In conventional weak-interaction theory, gauge invariance rules out long-range parity-violating forces at first order in G. The effect is therefore suppressed twice: once by the small coupling and again by the ratio of force range to electron-cloud size. His conclusion is that "no appreciable energy difference between optical isomers is expected" (p. 15).
  - On p. 16 he notes that a non-negligible difference could be accumulated by mechanisms such as Yamagata's (1966), acting over geological time.
- **Letokhov (from Darquié et al. 2010, p. 2, and Quack et al. 2022, §4).** Letokhov proposed detecting the parity-violating energy difference spectroscopically, as a difference between the infrared vibration–rotation line frequencies of the two separated enantiomers. Darquié et al. say he expected splittings of about 10⁻¹⁵–10⁻¹⁶ of the transition frequency. Quack et al. call this "a scheme originally proposed by Letokhov in 1975", which is still the basis of the current experiments.
  - Note: the CHFClBr test molecule came from a later Letokhov-group paper (Kompanets et al., 1976), not from the 1975 one.

**Current draft sentence (line 275).**
> "Rein (1974) discussed parity-violating effects of intramolecular interactions, and Letokhov (1975) discussed the resulting difference between the energy levels of left- and right-handed molecules; what each proposed beyond this is [CHECK FULL TEXT]."

**Proposed replacement.**
> "Rein (1974) showed that the weak interaction does distinguish mirror-image molecules, but concluded that within the weak-interaction theory of the time the energy difference would be tiny. Letokhov (1975) proposed measuring it as a small shift between the infrared vibration–rotation lines of the two enantiomers (Quack et al., 2022)."

**Confidence.**
- Rein: high. This is the author's own summary, read in full.
- Letokhov: medium-high. Two independent later papers, both read in full, agree, but I did not read Letokhov's own two pages.

**Side note.** The sentence before this one says the PVED "comes from weak neutral currents inside each molecule". Rein's summary works in "conventional" (pre-neutral-current) theory. Neutral currents entered this question with Garay and Hraskó (1975, *J. Mol. Evol.* 6, 77; Quack et al. 2022 ref. 189). The framing is fine for the field as a whole, but don't imply that Rein used neutral currents.

---

## 4. Bonner (1991)

**Reference (line 852).** Bonner, W.A. 1991, The origin and amplification of biomolecular chirality, Origins of Life and Evolution of the Biosphere, 21, 59–111. https://doi.org/10.1007/bf01809580

**Source actually read.** The publisher's free preview (https://page-one.springer.com/pdf/preview/10.1007/BF01809580), pp. 59–60: the full table of contents and §§1–2. The weak-force sections (pp. 77–89) were **not** read. The paper is closed (OpenAlex), and PubMed (PMID 1758688) has no abstract. The Semantic Scholar auto-summary is garbled and should not be used.

**What it says on the point.** The table of contents (pp. 59–60) confirms the scope. Under "Universal Processes" the review covers:
- the violation of parity (p. 77);
- the Vester–Ulbricht hypothesis (pp. 77–80);
- radioracemization (p. 80);
- direct effects of chiral particles: electrons, protons, positrons and muons (pp. 81–84);
- weak interactions, neutral currents and PVEDs (pp. 84–89);
- conclusions (p. 89).

I could not see what those conclusions say.

**Current draft sentence (line 281).**
> "Bonner (1991) reviewed the origin and amplification of biomolecular chirality [CHECK FULL TEXT for its conclusions on the weak force]."

**Proposed replacement.**
> "Bonner (1991) reviewed the origin and amplification of biomolecular chirality, including the beta-decay and PVED proposals."

The next sentences on Bonner (1996, 2000) already carry his weak-force conclusions, so nothing is lost.

**Confidence.** High for the hedged wording, which the table of contents supports. Do not attribute any conclusion to the 1991 paper unless someone reads pp. 84–89.

---

## 5. Avalos et al. (2000)

**Reference (line 834).** Avalos, M., Babiano, R., Cintas, P., Jiménez, J.L., Palacios, J.C. 2000, From parity to chirality: chemical implications revisited, Tetrahedron: Asymmetry, 11, 2845–2874. https://doi.org/10.1016/s0957-4166(00)00265-2

**Source actually read.** The paper is closed: OpenAlex lists it as closed, and ScienceDirect returned 403 in both curl and the browser. No abstract is available from Crossref, OpenAlex or Semantic Scholar, and the University of Extremadura record (opendata.unex.es/investiga/publicaciones/2000-143) has metadata only. The OpenAlex list of open papers that cite it shows only generic uses. In the ones I read in full:
- Darquié et al. (2010, arXiv:1007.3352, p. 1) cites it only for the statement that the origin of homochirality is unknown.
- Arxiv:0906.3391 cites it for chiral bias from external fields and polarized radiation.

**What it says on the point.** It could not be confirmed beyond the title. It is a 30-page review on the chemical implications of parity violation.

**Current draft sentence (line 281).**
> "Avalos et al. (2000) revisited the chemical implications of parity violation [CHECK FULL TEXT]."

**Proposed replacement.** **Cut.** The sentence adds only the title to the paragraph. If you want to keep the citation, fold it into a list of reviews instead, for example by changing the existing Quack (2002) sentence to begin "Reviews of this period (Avalos et al., 2000; Quack, 2002) note…". That works only if the claim that follows is Quack's alone, so the safest option is a plain cut. If you cut it, also delete reference line 834. It is cited nowhere else (line 203 is the H7 to-do list). The list then has 118 entries.

**Confidence.** High that cutting is safe. Its content was not checked.

---

## 6. Wesendrup et al. (2003)

**Reference (line 1052).** Wesendrup, R., Laerdahl, J.K., Compton, R.N., Schwerdtfeger, P. 2003, Biomolecular Homochirality and Electroweak Interactions. I. The Yamagata Hypothesis, J. Phys. Chem. A, 107, 6668–6673. https://doi.org/10.1021/jp022568v

**Source actually read.** Abstract only: the publisher abstract via OpenAlex, https://api.openalex.org/works/doi:10.1021/jp022568v. The full text is closed: OpenAlex lists it as closed, and ACS returned a Cloudflare challenge, which I did not pass. I checked the open papers that cite it (arXiv:1007.3352, 1407.3998, 1311.3978, 1807.09425, physics/0411072, quant-ph/0411034). All of them cite it generically and none describes its method. MacDermott et al. (2009, two OA Springer papers) and Chem. Soc. Rev. (2022) may describe it, but both publishers were redirected on this machine, so I could not read them.

**What it says (abstract).** It credits Yamagata (1966) with proposing that weak interactions between electrons and nucleons could set biomolecular homochirality. The authors tested this with relativistic parity-non-conserving calculations for amino acids. Their result is "no clear evidence that the naturally occurring l-amino acids are stabilized". The abstract gives no detail on the molecules, conformers or level of theory.

**Current draft sentence (line 281).**
> "Wesendrup et al. (2003) re-examined the Yamagata hypothesis itself [CHECK FULL TEXT for method and result]."

**Proposed replacement.**
> "Wesendrup et al. (2003) tested the Yamagata hypothesis with relativistic calculations of the parity-violating energy for amino acids and found no clear evidence that parity violation stabilizes the natural L-amino acids."

**Confidence.** Medium-high, because the replacement says no more than the abstract. Do not add the number of amino acids, the conformers or the method name unless someone reads the full text.

**Related flag (line 273, the Yamagata sentence).** The draft cites Wesendrup et al. (2003) as the later work that "summarises [Yamagata's] core as a tiny bias that accumulates through repeated steps or amplification".
- The title supports the name "the Yamagata hypothesis", but the abstract does not support the accumulation summary.
- The source I did read, Rein (1974, p. 16), describes Yamagata's idea as an accumulation mechanism acting over geological periods.
- Safer wording: "…calls the link … 'the Yamagata hypothesis' (Wesendrup et al., 2003); Rein (1974) described it as an accumulation mechanism acting over geological time."
- Whoever has the Yamagata item should confirm this against Yamagata (1966) itself.
