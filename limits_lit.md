# What the literature says against the model, and what is left of it

Checked 2026-10-10. A literature agent found 25 papers the project did not cite; every DOI was checked on Crossref
(the preprint on arXiv). Summaries come from each paper's abstract, not its full text. For each limit: what the paper
says, the best answer the project can give, and the unknown that is left. Arithmetic here uses numbers already in
PAPER_DRAFT.md and the logs; no reported number is changed.

## 1. The weak force's energy difference (PVED) has never been measured in a molecule
- **The papers.** Laser spectra of the two hands of CHFClBr matched within 13 Hz (Δν/ν < 4×10⁻¹³; Daussy et al. 1999,
  *PRL* 83:1554, doi:10.1103/physrevlett.83.1554). A detection in Ru/Os complexes looks feasible but has not been made
  (Fiechter et al. 2022, *J. Phys. Chem. Lett.* 13:10011, doi:10.1021/acs.jpclett.2c02434).
- **Answer.** (a) The model never uses one PVED value: it draws g over 4×10⁻¹⁹ to 4×10⁻¹⁶, a factor of 1000, and g moves
  the win rate by only +8/+4 points (Enceladus/Europa, §3.9). A theory error of 10× either way stays inside that range.
  (b) The force itself is measured: parity violation in cesium atoms was measured precisely and tests the Standard
  Model at low energy (Wood et al. 1997, *Science* 275:1759, doi:10.1126/science.275.5307.1759). What is unmeasured is
  its size in molecules, not whether it exists. (c) The experimental limits sit far above the predicted values, so they
  neither confirm nor contradict the theory.
- **Left.** The true PVED of the molecules that actually amplified. It matters only if theory is off by more than the
  1000× range already sampled.

## 2. Minerals can push one hand far harder than the weak force
- **The papers.** Natural pyrite makes D-amino acids with 14.5–42.4% excess (Li et al. 2024, *Nat. Commun.* 15:10130,
  doi:10.1038/s41467-024-54481-x). UV on magnetite gives spin-polarized electrons that favour one hand (Ozturk &
  Sasselov 2022, *PNAS* 119:e2204765119, doi:10.1073/pnas.2204765119). Crystallizing one-handed RAO magnetizes
  magnetite "akin to an avalanche" (Ozturk et al. 2023, *Nat. Commun.* 14:6351, doi:10.1038/s41467-023-42130-8).
  Prebiotic magnetite's vortex-state grains re-magnetize irreversibly with one-handed compounds, which can store and
  reinforce a weak bias (Devienne et al. 2026, arXiv:2605.19387, preprint, simulations). Magnetized surfaces adsorb
  one hand faster (Banerjee-Ghosh et al. 2018, *Science* 360:1331, doi:10.1126/science.aar4265) and separate some
  amino acids on crystallization but not others (Tassinari et al. 2019, *Chem. Sci.* 10:5246, doi:10.1039/c9sc00663j).
- **Answer.** (a) Li 2024 and Ozturk & Sasselov 2022 need light. The seafloors of Enceladus, Europa and Ganymede lie
  under kilometres of ice and are dark, so these routes do not apply there (they can apply on sunlit early Earth).
  Devienne's magnetite also forms by UV photo-oxidation or nitrite oxidation, surface-Earth routes. (b) The avalanche
  and vortex results are amplifiers: in both abstracts the magnetization follows the hand of molecules that are already
  one-handed. They store a bias; they do not choose it. If the weak force tips the molecules first, they amplify the
  weak force's choice. (c) Where the magnetization is set at random, so is the favoured hand, and random patches
  average out: `nudge_realism.py` Part A finds they are harmless if one hand holds over less than about 0.03–40 µm
  (all chemistry on surfaces) or 3 µm–4 mm (1% on surfaces). (d) The one push with a shared hand, the planet's field
  aligning magnetite, is already quantified (`nudge_realism.py` Part B): below ~10⁻¹⁵ of the deciding molecules on
  field-aligned magnetite on Enceladus and Europa.
- **Left.** Two geology numbers, both unmeasured: how far one hand extends across seafloor minerals (the avalanche
  effect could lengthen it), and how much of the deciding chemistry happens on field-aligned magnetite.

## 3. Titan's amino acids may come from warm melt water, not 94 K chemistry
- **The papers.** Aminonitriles turn into glycine, alanine and GABA in liquid water over 6 months; on Titan that would
  be in short-lived melts after impacts or cryovolcanism (Farnsworth et al. 2024, *ACS Earth Space Chem.* 8:2380,
  doi:10.1021/acsearthspacechem.4c00114). Dragonfly targets Selk crater, where melt water mixed with organics (Barnes et
  al. 2021, *PSJ* 2:130, doi:10.3847/psj/abfdcf).
- **Answer.** (a) A second route existing does not stop the cold one. (b) Dragonfly samples both: it first lands at a
  dune site near Selk (Lorenz et al. 2021, *PSJ*, doi:10.3847/psj/abd08f), then reaches the water-rich crater
  material, so it can compare cold-surface organics with melt products. (c) Its mass spectrometer is being adapted to
  separate the two hands (Boulesteix et al. 2024, *J. Chromatogr. A* 1722:464860, doi:10.1016/j.chroma.2024.464860).
- **Left.** Whether chiral molecules form at 94 K at all, and how fast (next item).

## 4. The Titan tunnelling rates are assumed
- **The papers.** Hydrogen tunnels under barriers of 28–30 kcal/mol (about 117–125 kJ/mol) with half-lives of 1–2 hours
  at 11 K (Schreiner et al. 2008, *Nature* 453:906, doi:10.1038/nature07010; 2011, *Science* 332:1300,
  doi:10.1126/science.1203761). H atoms move across ice by tunnelling at 10 K (Kuwahata et al. 2015, *PRL* 115:133201,
  doi:10.1103/physrevlett.115.133201). Reviews: Meisner & Kästner 2016, *Angew. Chem.* 55:5400,
  doi:10.1002/anie.201511028; Hama & Watanabe 2013, *Chem. Rev.* 113:8783, doi:10.1021/cr4000978 (title only).
- **Answer.** (a) That light-atom chemistry keeps running far below where heat alone could drive it is measured, not
  assumed. (b) The bar is low. The model needs only k₂cτ ≥ 10; at `titan_95k.py`'s 3.2 µM that is k₂ ≥ 1.0×10⁻⁹ /M/s
  over 100 Myr, or 2.5×10⁻¹¹ /M/s over 4 Gyr. The hand-set 10⁻⁸ /M/s clears the 100 Myr bar by 10×.
- **Left.** No measured rate for a two-molecule proton transfer in liquid methane/ethane at 94 K. The anchors above are
  single-molecule shifts in frozen argon or atoms on ice. The same kind of unmeasured rate sets how fast products lose
  their hand (P5 of `titan_95k.py`, half-life 2.5×10⁶ yr, ASSUMED).

## 5. How fast Enceladus' seafloor reaches the plume
- **The papers.** Silica nanoparticles point to ongoing >90 °C water-rock reactions whose products reach the plume
  "quickly" (Hsu et al. 2015, *Nature* 519:207, doi:10.1038/nature14262), against "tens of years or more" (Zhang et al.
  2025, already cited) and hundreds of years in a stratified ocean (Kang et al. 2022b; Ames et al. 2025, already cited).
- **Answer.** Fast transport means strong vertical mixing, which helps the model: the high half of the vertical-mixing
  range gives 43%/57% (Enceladus/Europa) against 0%/0% in the low half (§3.9). If Hsu is right, Enceladus' odds go up.
  All three estimates fall inside the range the 22% already samples.
- **Left.** Vertical mixing, already the paper's main unknown (`vertical_mixing_lit.md`).

## What the limits come down to
1. Vertical mixing in icy-moon oceans (unchanged; the main unknown).
2. How far one hand extends across seafloor minerals, and how much chemistry happens on field-aligned magnetite.
3. Reaction and racemization rates at 94 K on Titan.
4. The reverse rates that decide how much Bennu counts against the model (unchanged, §3.10).
The PVED's exact size matters only if theory is off by more than 1000×.

## Papers that support the model without needing an answer
- Open, driven systems are needed for deracemization; closed ones are not enough: Buhse et al. 2021, *Chem. Rev.*
  121:2147, doi:10.1021/acs.chemrev.0c00819.
- Organics in Enceladus' ocean: Khawaja et al. 2019, *MNRAS* 489:5231, doi:10.1093/mnras/stz2280 (amines, carbonyls; no
  amino acid identified); Khawaja et al. 2025, *Nat. Astron.* 9:1662, doi:10.1038/s41550-025-02655-y.
- Enceladus' ocean pH 7.95–9.05 (sets racemization speed): Fifer et al. 2022, *PSJ* 3:191, doi:10.3847/psj/ac7a9f.
- Europa's surface CO₂ comes from inside: Trumbo & Brown 2023, *Science* 381:1308, doi:10.1126/science.adg4155.
- Ryugu's amino acids are racemic or nearly so: Parker et al. 2023, *GCA* 347:42, doi:10.1016/j.gca.2023.02.017;
  aqueous parent-body origin: Naraoka et al. 2023, *Science* 379:eabn9033, doi:10.1126/science.abn9033.
- A lander that could test the plume prediction: MacKenzie et al. 2021, *PSJ* 2:77, doi:10.3847/psj/abe4da (Orbilander;
  the abstract names no chiral measurement).

## Read 2026-10-11 (full text; PDF from Dr. McKay)
Lunine, Beauchamp, Smith & Nikolaev 1999, "The Abiotic Generation of Homochirality on Saturn's Moon Titan",
*Advances in BioChirality* pp. 257–270, doi:10.1016/B978-008043404-9/50017-6. A perspective chapter; no model, no data.
- **Same question, no answer.** It asks whether the PVED decides the hand during amplification (citing Kondepudi &
  Nelson 1985) and proposes Titan as a natural lab for abiotic organic chemistry, too cold for life. It computes nothing.
- **Scale.** It cites Zanasi (same conference) for the idea that seeing PVED-driven symmetry breaking needs a system much
  larger than a lab but smaller than a planet's surface: a qualitative version of the model's 10³³–10³⁴ molecules.
  The Zanasi chapter is not read; it is in the same book.
- **Low-temperature rates.** Arrhenius extrapolation puts racemization at Titan temperatures beyond the age of the
  Universe, but the chapter warns that extrapolating one Arrhenius law that far is not justified: other pathways,
  including rates that rise as temperature falls in hydrocarbon liquids, may take over. This is the same caveat as
  the 94 K result (titan_95k.py): the rates there are assumed, not measured.
- **Liquid water on Titan.** A 1 km impactor could melt crater ice for thousands of years (longer with ammonia). Run
  on such pools 2026-10-11 (impact_melt.py): weak force wins in 33% of input space for a 1 km-impactor crater,
  51% for Selk (Wakita et al. 2023 melt volumes) and 74% for Menrva; concentration decides most, and any inherited
  excess above ~1e-16 overrides it.
- **Effect on the project.** None on any number. The RNAAS sentence ("discussed abiotic homochirality on Titan; we
  consider liquid-water oceans only") is accurate.
