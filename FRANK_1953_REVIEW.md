# Frank (1953) Citation Review: Astrobiological & Mechanistic Synthesis

> **STATUS: DRAFT, not for sending.** All 905 citing records were classified on 2026-10-10 (method and reliability in §1.2). Every DOI in this file was checked against Crossref (`check_refs.py`). The one-line notes are paraphrases of titles and abstracts only; no full texts were read.

*Prepared for:* Prebiotic Chirality Project (internal)
*Author:* Saptarshi Ghosh | Cupertino High School, CA | Project Repository: `~/Documents/chirality-project`
*Date:* October 2026

---

## 1. What the 905 citing papers are about

Frank's paper, **"On spontaneous asymmetric synthesis"** (*Biochimica et Biophysica Acta* 11: 459–463, 1953), set out the classic recipe for chemical homochirality: an open system combining **autocatalysis** ($A + L \to 2L$, $A + D \to 2D$) with **mutual antagonism** ($L + D \to$ waste or an inactive heterochiral complex).

OpenAlex lists **905 citing records** for Frank (1953) (work W2026655238, 2026-10-10). About 50 of them repeat an earlier title (duplicate records, or generic titles such as "Chirality"), so there are roughly 855 distinct works. The full list, with each paper's group and a one-line note, is in [`frank_citers.csv`](file:///Users/saptarshighosh/Documents/chirality-project/frank_citers.csv).

### 1.1 Counts

| Group | Papers | What it covers |
| :--- | ---: | :--- |
| Models (THEORY) | 210 | kinetic and stochastic models of mirror-symmetry breaking: Frank-type networks, bifurcations, noise |
| Not about homochirality (OTHER) | 139 | cite Frank for pattern formation, materials, supramolecular chemistry, etc. |
| Reviews (REVIEW) | 119 | reviews, perspectives, book chapters |
| Crystals (CRYSTAL) | 108 | crystallisation, Viedma ripening, stirred crystallisation, conglomerates |
| Soai reaction (SOAI) | 92 | Soai asymmetric autocatalysis experiments |
| Prebiotic chemistry (PREBIO) | 82 | homochiral peptides, RNA, sugars, polymerisation |
| Sources of a bias (BIAS) | 81 | parity violation (PVED), circularly polarised light, magnetic fields, spin, chiral minerals, stirring |
| Solution catalysis (CATAL) | 42 | amplification in solution: non-linear effects, proline, peptide catalysts |
| Space and planets (ASTRO) | 32 | meteorites, comets, planetary settings, biosignatures |

By decade the citations run 1950s 3, 1960s 1, 1970s 27, 1980s 60, 1990s 66, 2000s 229, 2010s 270, 2020s (to Oct 2026) 249. Frank was barely cited for 25 years. The 1980s burst is the parity-violation debate: 20 of that decade's 60 citers are in the bias group, mostly PVED calculations (Mason, Tranter) and the Kondepudi–Nelson papers. The post-2000 growth is Soai, Viedma and model papers.

Three cross-cutting counts:
* **64 papers** mention parity violation or the weak force (parity, electroweak, weak neutral current, PVED) in their title or abstract. This is a plain keyword count, and a lower bound, because 310 of the 905 have no abstract in OpenAlex.
* **83 papers** treat a small fixed bias competing with random noise, which is the project's question. Another 298 could not be judged from the title alone.
* **35 papers** were rated most relevant to this project (§3).

### 1.2 How this was done, and how far to trust it

* Each record's title and abstract (if any) went to one of three Claude subagents with the same written instructions. Each agent gave the paper one group, two yes/no/unclear flags (bias-vs-noise; planetary or open-flow setting), a 0–3 relevance score and a ≤20-word note, using the title and abstract only. Merged output: 905 records, each once.
* **Spot check (30 random records, checked by hand against title and abstract):**
  * 29 of 30 groups were reasonable. "Why 43 codons?" was put in prebiotic chemistry but belongs in OTHER.
  * The flags are noisier: 4 of 30 were debatable.
    * Three papers were marked "planetary/open" without clear support: a meteorite-impact plasma paper, a reaction-model paper and a peptide-stability paper.
    * Singleton & Vo 2002 was marked "bias vs noise", but it is about amplifying tiny accidental imbalances.
  * Treat the group counts as good to about ±10%, and the flag counts as rough.
* Every note on a record without an abstract says "title only" (checked for all 310).
* These are machine labels on titles and abstracts. They are a map for reading, not a reading.

### 1.3 Does any citer already put this on an icy moon?

Searching all titles and abstracts for Enceladus, Europa, Ganymede, Titan, "icy", "ocean world" and "subsurface ocean" found **no paper on Enceladus, Europa, Ganymede or subsurface oceans**. There is **one Titan paper:** Lunine, Beauchamp, Smith & Nikolaev (1999), *The abiotic generation of homochirality on Saturn's moon Titan*, in *Advances in BioChirality*, pp. 257–270 ([doi:10.1016/B978-008043404-9/50017-6](https://doi.org/10.1016/B978-008043404-9/50017-6)). OpenAlex has no abstract for it, and it has **not been read**. It must be read before anything is said about Titan, and before claiming that no one has done this. This search covers only Frank's citers. The broader prior-art scan in PAPER_DRAFT §1 covered the Kondepudi & Nelson (1985) and Brandenburg & Multamäki (2004) citers.

---

## 2. Selected Papers That Cite Frank (1953)

Every entry below resolves on Crossref to the paper named, and OpenAlex records it as citing Frank (1953).

### Pillar I: Non-Equilibrium Thermodynamics & Open-vs-Closed Systems
1. **Frank, F. C. (1953).** *On spontaneous asymmetric synthesis.* Biochim. Biophys. Acta, 11, 459–463. [doi:10.1016/0006-3002(53)90082-1](https://doi.org/10.1016/0006-3002(53)90082-1) *(the source paper itself)*
   * Autocatalysis plus mutual antagonism in an open system can amplify a small chiral imbalance.
2. **Decker, P. (1973).** *Possible resolution of racemic mixtures by bistability in "bioids", open systems which can exist in several steady states.* J. Mol. Evol., 2, 137–143. [doi:10.1007/BF01653993](https://doi.org/10.1007/BF01653993)
   * Open, multistable "bioid" systems as a route to resolving racemates.
3. **Plasson, R., Kondepudi, D. K., Bersini, H., Commeyras, A., & Asakura, K. (2007).** *Emergence of homochirality in far-from-equilibrium systems: Mechanisms and role in prebiotic chemistry.* Chirality, 19(8), 589–600. [doi:10.1002/chir.20440](https://doi.org/10.1002/chir.20440)
   * Review of far-from-equilibrium mechanisms for homochirality and their prebiotic relevance.
4. **Blanco, C., Stich, M., & Hochberg, D. (2011).** *Temporary mirror symmetry breaking and chiral excursions in open and closed systems.* Chem. Phys. Lett., 505, 140–147. [doi:10.1016/j.cplett.2011.02.032](https://doi.org/10.1016/j.cplett.2011.02.032)
   * In closed systems an enantiomeric excess can appear temporarily (a "chiral excursion") before the system returns to racemic.
5. **Stich, M., Ribó, J. M., Blackmond, D. G., & Hochberg, D. (2016).** *Necessary conditions for the emergence of homochirality via autocatalytic self-replication.* J. Chem. Phys., 145, 074111. [doi:10.1063/1.4961021](https://doi.org/10.1063/1.4961021)
   * Conditions an autocatalytic self-replication scheme must meet to give homochirality.
6. **Higgs, P. G. & Blackmond, D. G. (2025).** *Autocatalytic symmetry breaking and chiral amplification in a feedback network combining amino acid synthesis and ligation.* Proc. Natl. Acad. Sci. USA, 122(20), e2423683122. [doi:10.1073/pnas.2423683122](https://doi.org/10.1073/pnas.2423683122)
   * In this network the forward reaction alone gives no ee amplification (homochiral dimers catalyse the opposite-hand monomer); closed systems return to racemic, while open flow can keep broken symmetry.

### Pillar II: Blackmond and Coworkers (Heterochiral Dimers, Amplification)
7. **Blackmond, D. G. (2004).** *Asymmetric autocatalysis and its implications for the origin of homochirality.* Proc. Natl. Acad. Sci. USA, 101(16), 5732–5736. [doi:10.1073/pnas.0308363101](https://doi.org/10.1073/pnas.0308363101)
   * Interprets Frank-type mutual antagonism through reversible homochiral/heterochiral dimer formation in the Soai reaction.
8. **Mathew, S. P., Iwamura, H., & Blackmond, D. G. (2004).** *Amplification of enantiomeric excess in a proline-mediated reaction.* Angew. Chem. Int. Ed., 43(25), 3317–3321. [doi:10.1002/anie.200453997](https://doi.org/10.1002/anie.200453997)
   * Experimental ee amplification in a proline-mediated reaction.
9. **Klussmann, M., Iwamura, H., Mathew, S. P., Wells, D. H., Pandya, U., Armstrong, A., & Blackmond, D. G. (2006).** *Thermodynamic control of asymmetric amplification in amino acid catalysis.* Nature, 441, 621–623. [doi:10.1038/nature04780](https://doi.org/10.1038/nature04780)
   * Solid–solution phase behavior (eutectic ee) of amino acids sets the solution ee, giving amplification without an autocatalytic reaction.
10. **Blackmond, D. G. (2010).** *The origin of biological homochirality.* Cold Spring Harb. Perspect. Biol., 2(5), a002147. [doi:10.1101/cshperspect.a002147](https://doi.org/10.1101/cshperspect.a002147)
    * Review of chemical and physical routes to homochirality.

### Pillar III: Stochastic Bifurcations & Parity Violation
11. **Kondepudi, D. K. & Nelson, G. W. (1985).** *Weak neutral currents and the origin of biomolecular chirality.* Nature, 314, 438–441. [doi:10.1038/314438a0](https://doi.org/10.1038/314438a0)
    * Estimates how a tiny parity-violating bias can select the branch at a Frank-type bifurcation when noise is small enough.
12. **Lente, G. (2004).** *Homogeneous chiral autocatalysis: A simple, purely stochastic kinetic model.* J. Phys. Chem. A, 108(44), 9475–9478. [doi:10.1021/jp046413u](https://doi.org/10.1021/jp046413u)
    * Purely stochastic model of chiral autocatalysis in small systems.
13. **Gleiser, M. & Walker, S. I. (2008).** *An extended model for the evolution of prebiotic homochirality: A bottom-up approach to the origin of life.* Orig. Life Evol. Biosph., 38(4), 293–315. [doi:10.1007/s11084-008-9134-5](https://doi.org/10.1007/s11084-008-9134-5)
    * Extended Frank-type polymerization model including environmental effects.
14. **Gleiser, M. & Walker, S. I. (2012).** *Life's chirality from prebiotic environments.* Int. J. Astrobiol., 11(4), 287–296. [doi:10.1017/S1473550412000377](https://doi.org/10.1017/S1473550412000377)
    * Role of the prebiotic environment (noise, disturbances) in chiral selection.

### Pillar IV: Spatial Coupling & Crystal Deracemization
15. **Brandenburg, A. & Multamäki, T. (2004).** *How long can left and right handed life forms coexist?* Int. J. Astrobiol., 3(3), 209–219. [doi:10.1017/S1473550404001983](https://doi.org/10.1017/S1473550404001983)
    * Reaction–diffusion model of competing chiral domains and how long they coexist.
16. **Viedma, C. (2005).** *Chiral symmetry breaking during crystallization: Complete chiral purity induced by nonlinear autocatalysis and recycling.* Phys. Rev. Lett., 94(6), 065504. [doi:10.1103/PhysRevLett.94.065504](https://doi.org/10.1103/PhysRevLett.94.065504)
    * Grinding-driven deracemization of crystals ("Viedma ripening").
17. **Noorduin, W. L. et al. (2008).** *Emergence of a single solid chiral state from a nearly racemic amino acid derivative.* J. Am. Chem. Soc., 130(4), 1158–1159. [doi:10.1021/ja7106349](https://doi.org/10.1021/ja7106349)
    * Viedma-type deracemization applied to a racemizable amino acid derivative.
18. **Ribó, J. M. & Hochberg, D. (2020).** *Spontaneous mirror symmetry breaking: an entropy production survey of the racemate instability and the emergence of stable scalemic stationary states.* Phys. Chem. Chem. Phys., 22, 14013–14025. [doi:10.1039/d0cp02280b](https://doi.org/10.1039/d0cp02280b)
    * Entropy-production analysis of when the racemic state becomes unstable.

### Related papers that do NOT cite Frank (1953)

Used by the project but not part of the Frank citation lineage (OpenAlex shows no citation of W2026655238):

* **Hein, J. E., Tse, E., & Blackmond, D. G. (2011).** *A route to enantiopure RNA precursors from nearly racemic starting materials.* Nat. Chem., 3, 704–706. [doi:10.1038/nchem.1108](https://doi.org/10.1038/nchem.1108)
* **Kondepudi, D. K. & Nelson, G. W. (1983).** *Chiral symmetry breaking in nonequilibrium systems.* Phys. Rev. Lett., 50(14), 1023–1026. [doi:10.1103/PhysRevLett.50.1023](https://doi.org/10.1103/PhysRevLett.50.1023)
* **Goldanskii, V. I. & Kuz'min, V. V. (1988).** *Spontaneous mirror symmetry breaking in nature and the origin of life.* AIP Conf. Proc., 180, 163–228. [doi:10.1063/1.37867](https://doi.org/10.1063/1.37867)
* **McKay, C. P. (2016).** *Titan as the abode of life.* Life, 6(1), 8. [doi:10.3390/life6010008](https://doi.org/10.3390/life6010008)
* **Trixler, F. (2013).** *Quantum tunnelling to the origin and evolution of life.* Curr. Org. Chem., 17(16), 1758–1770. [doi:10.2174/13852728113179990083](https://doi.org/10.2174/13852728113179990083) — a general review; gives no rate constants or crossover temperatures.
* **Yu, Darú, Deng & Blackmond (2024).** *Prebiotic access to enantioenriched amino acids via peptide-mediated transamination reactions.* Proc. Natl. Acad. Sci. USA, 121(7), e2315447121. [doi:10.1073/pnas.2315447121](https://doi.org/10.1073/pnas.2315447121) — the dipeptide/ee values the project used as a "Yu benchmark" are not in this paper.
* **Longo, H. W. & Remsing, R. C. (2024).** *Quantum tunneling could enable proton transfer reactions on Titan.* arXiv:2405.19195. [arxiv.org/abs/2405.19195](https://arxiv.org/abs/2405.19195) — path-integral simulations of proton-transfer thermodynamics in liquid ethane; it gives **no rate constants**. Any tunnelling rate used in the project is an assumption, not from this paper.

*Removed (do not exist):* "Yu, Wang & Sutherland 2024, PNAS e2318721121" and a "Higgs & Blackmond, Orig. Life Evol. Biosph. 2025" paper (its DOI is not registered). The real Higgs & Blackmond paper is entry 6 above.

---

## 3. The papers closest to this project

The 35 citer records rated most relevant (relevance 3), plus three rated 2 that bear on a specific script. Several have no abstract in OpenAlex and are marked *title only*: what they conclude is unknown until they are read. **Bold** means the project does not yet cite the paper.

### 3.1 Small bias against noise in Frank-type models (the project's core question)
* **Morozov, Kuz'min & Goldanskii (1983).** *Comparative analysis of the role of statistical fluctuations and factor of advantage (parity nonconservation) in the origins of optical activity.* Orig. Life 13: 119–138. [doi:10.1007/BF00928890](https://doi.org/10.1007/BF00928890). Title only. This is the earliest citer that poses the project's exact question, noise against a parity-violation advantage. **Read first.**
* **Tennakone (1984)**, Chem. Phys. Lett. 105: 444–446 ([doi:10.1016/0009-2614(84)80060-3](https://doi.org/10.1016/0009-2614(84)80060-3)); **Nelson (1984)**, Springer Proc. Phys., pp. 245–248 ([doi:10.1007/978-3-642-46508-6_29](https://doi.org/10.1007/978-3-642-46508-6_29)); **Kondepudi (1986)**, Orig. Life Evol. Biosph. 16: 424–425 ([doi:10.1007/BF02422111](https://doi.org/10.1007/BF02422111)); **Avetisov, Kuz'min & Anikin (1987)**, *Sensitivity of chemical chiral systems to weak asymmetric factors*, Chem. Phys. 112: 179–187 ([doi:10.1016/0301-0104(87)80160-X](https://doi.org/10.1016/0301-0104(87)80160-X)). All title only. These are 1980s models of weak-neutral-current selection, contemporary with Kondepudi & Nelson (1985), which the project already cites.
* **Kondepudi & Asakura (2001).** *Chiral autocatalysis, spontaneous symmetry breaking, and stochastic behavior.* Acc. Chem. Res. 34: 946–954. [doi:10.1021/ar010089t](https://doi.org/10.1021/ar010089t). Account of chirally autocatalytic systems, stochastic behaviour, and how sensitive symmetry breaking is to small biases.
* **Lente (2005)**, J. Phys. Chem. A 109: 11058–11063 ([doi:10.1021/jp054613f](https://doi.org/10.1021/jp054613f)), and Lente (2007), PCCP 9: 6134 ([doi:10.1039/b711546f](https://doi.org/10.1039/b711546f)). These are stochastic (discrete-molecule) models. The 2007 abstract concludes that the PVED "cannot be amplified under conditions reasonable" for biological chirality. **This contradicts the project.** The project's answer, that Lente's limit is for closed systems and does not carry over to fed ones, is in `objections.py` §1 and is a model argument, not a refutation.
* **Cruz, Parmananda & Buhse (2008).** *Noise-induced enantioselection in chiral autocatalysis.* J. Phys. Chem. A 112: 1673–1676. [doi:10.1021/jp077415s](https://doi.org/10.1021/jp077415s). In their model, adding noise within a certain range lets a bias too small to act alone pick the hand consistently. **This runs the opposite way to the project's assumption that noise only works against the bias.** It has not been tested in `frank.py`.
* **Hochberg (2009).** *Mirror symmetry breaking and restoration: the role of noise and chiral bias.* Phys. Rev. Lett. 102: 248101. [doi:10.1103/PhysRevLett.102.248101](https://doi.org/10.1103/PhysRevLett.102.248101). Uses a non-equilibrium effective potential for the Frank model with external noise. Above a critical noise level, mirror symmetry is restored.
* **Crusats, Hochberg, Moyano & Ribó (2009)**, ChemPhysChem 10: 2123–2131 ([doi:10.1002/cphc.200900181](https://doi.org/10.1002/cphc.200900181)). A closed-system Frank variant gives absolute asymmetric synthesis by amplifying the initial statistical deviations. **Silva-Dias & López-Castillo (2017)**, PCCP 19: 29424–29428 ([doi:10.1039/c7cp04674j](https://doi.org/10.1039/c7cp04674j)). Compares stochastic and deterministic symmetry breaking in the Frank model.
* Hochberg et al. (2022, Phys. Rev. Research 4: 033183) and Hochberg et al. (2023, PCCP 25: 31583). PVED bias against temperature-driven noise in an open-flow Frank model, with a selectivity criterion. These are already cited and are **the closest work in method**; PAPER_DRAFT should state plainly what the project adds beyond them.

### 3.2 Experiments that tried to see the PVED, or propose to
* **Viedma (2007).** *Selective chiral symmetry breaking during crystallization: parity violation or cryptochiral environment in control?* Cryst. Growth Des. 7: 553–556. [doi:10.1021/cg060698d](https://doi.org/10.1021/cg060698d). In NaClO₃/NaBrO₃ crystallisation, a non-random handedness appeared. The abstract rules out PVED and blames trace chiral impurities. This is an experimental case of `moon_nudges.py`'s concern, where a local chiral push outvotes the weak force.
* Singleton & Vo (2002, 2003): already cited. Soai amplification of tiny accidental imbalances.
* **Cowan (2024).** *How to experimentally determine the influence of the weak nuclear force in promoting molecular chirality.* ChemistrySelect 9. [doi:10.1002/slct.202402950](https://doi.org/10.1002/slct.202402950). The project cites Cowan & Furnstahl 2022 and Cowan 2023, but not this paper. Proposes heavy-metal-catalysed autocatalytic Strecker chemistry to make a weak-force preference visible within weeks to months. This is a possible lab test.

### 3.3 How big the bias is
* **Zanasi, Lazzeretti, Ligabue & Soncini (1999).** Phys. Rev. E 59: 3382–3385. [doi:10.1103/PhysRevE.59.3382](https://doi.org/10.1103/PhysRevE.59.3382). Their more accurate (random-phase approximation) calculations gave PVED values for L-valine and D-glyceraldehyde about 10× larger than earlier methods. This bears on the size of g the project assumes.
* **Tranter (1985).** *Parity-violating energy differences of chiral minerals and the origin of biomolecular homochirality.* Nature 318: 172–173. [doi:10.1038/318172a0](https://doi.org/10.1038/318172a0). Title only; rated 2. Relevant to chiral mineral surfaces in `nudge_realism.py`.
* Quack et al. (2022) is already cited.

### 3.4 Large-scale and planetary settings
* **Mason (1991).** *Prebiotic sources of biomolecular handedness.* Chirality 3: 223–226. [doi:10.1002/chir.530030403](https://doi.org/10.1002/chir.530030403). Argues for the electroweak PVED combined with selection kinetics in a flow reactor. This is the same open-flow framing as the project's.
* **Ribó & Hochberg (2019).** *Chemical basis of biological homochirality during the abiotic evolution stages on Earth.* Symmetry 11: 814. [doi:10.3390/sym11060814](https://doi.org/10.3390/sym11060814). In dissipative networks, a residual bias from physical chiral forces can make the random sign deterministic.
* **Globus et al. (2021).** *Polarized muons and the origin of biological homochirality.* arXiv:2110.01975. [doi:10.48550/arXiv.2110.01975](https://doi.org/10.48550/arXiv.2110.01975). Proposes a competing weak-force route: spin-polarised cosmic-ray muons. It discusses several environments and other worlds. The project has not compared it with the PVED route.
* **Hoang (2024).** *Magnetically aligned grains as a ubiquitous source of spin-polarized electrons and chiral-symmetry-breaking agents.* Astrophys. J. 976: 26. [doi:10.3847/1538-4357/ad7f52](https://doi.org/10.3847/1538-4357/ad7f52). Rated 2. Field-aligned grains as a chiral agent, in interstellar icy grain mantles. Compare with `nudge_realism.py` Part B, where the field is aligned on magnetite.
* Gleiser & Walker (2012): already cited.

The other relevance-3 papers are Rivera-Islas et al. (2005, Soai kinetics), Gleiser (2006), Mason (2007), Chandrasekhar (2007, already cited), Fuß (2008), MacDermott (2012), Konstantinov (2013, 2022) and Stovbun (2025). They are listed in `frank_citers.csv` and were not summarised here.

**Reading order suggested by this review:** Morozov et al. 1983 → Lunine et al. 1999 (Titan) → Cruz et al. 2008 → Hochberg 2009 → Viedma 2007 → Lente 2007. The 35 records are 34 distinct papers (Kondepudi & Asakura 2001 appears twice). The project already cites 9 of them; 25 are new to it.

---

## 4. Four Debates in the Frank Tradition

### Debate 1: Is "mutual destruction" physically plausible?
* **The criticism:** Frank assumed $L + D \to W$. Why would molecules specifically remove their mirror twins?
* **Blackmond's proposal (2004):** antagonism can come from reversible heterochiral dimers ($L + D \rightleftharpoons [LD]$). If $[LD]$ is stable and inactive, it acts as a sink for the minority hand, reproducing the mathematical form of Frank's antagonism without a waste product.

### Debate 2: Open vs. closed systems (Bennu/Ryugu)
* **The question:** Can homochirality survive in a sealed asteroid or a closed pond?
* **Literature view (Blanco et al. 2011; Higgs & Blackmond 2025):** in a closed system the ee can rise temporarily but returns to racemic; open flow can sustain broken symmetry.
* **Project hypothesis (`review.py` R6, model only):** racemic amino acids in Bennu/Ryugu fit a closed parent body whose heat source died out, while tidally heated ocean worlds could stay driven for longer. This is a model interpretation, not a tested result.

### Debate 3: Counting noise vs. the electroweak bias
* **The concern:** the parity-violating energy difference is tiny (order $10^{-17}\,k_BT$), so in small volumes noise should swamp it.
* **Kondepudi & Nelson (1985)** estimated how bias and noise compete at a slowly crossed bifurcation.
* **Project model (`ocean.py`):** uses a selection parameter that grows with $\sqrt{N}$, so in the model a whole ocean's worth of molecules can make the bias decisive. This is a model result; the scaling and its constants have not been compared with experiment.

### Debate 4: Spatial coarsening vs. ocean mixing
* **The question:** does turbulence destroy chiral domains faster than they merge?
* **Project model (`aniso3d.py`, `sphere_layers.py`):** in these runs, vertical mixing ($D_z$) rather than horizontal mixing limited the final ee. The per-moon percentages these scripts print are model outputs and are not repeated here, because the attribution of some of them (e.g. the Enceladus range) is under review.

---

## 5. Mapping to the Project Codebase

| Script | What it models (model only, not validated against experiment) |
| :--- | :--- |
| [`frank.py`](file:///Users/saptarshighosh/Documents/chirality-project/frank.py) | Open-flow Frank model with Langevin counting noise; computes a noise penalty factor. |
| [`frank2.py`](file:///Users/saptarshighosh/Documents/chirality-project/frank2.py) | Frank network with reversible back-reactions obeying detailed balance. |
| [`review.py`](file:///Users/saptarshighosh/Documents/chirality-project/review.py) | R1: noise-vs-bias coherence volume. R6: closed-system decay to racemic. |
| [`dna_network.py`](file:///Users/saptarshighosh/Documents/chirality-project/dna_network.py) | Coupled cross-catalytic Frank network (internal self-consistency checks only). |
| [`titan_95k.py`](file:///Users/saptarshighosh/Documents/chirality-project/titan_95k.py) | Tunnelling vs. Arrhenius kinetics at 94 K. Tunnelling rates are **assumed**; Longo & Remsing (2024) and Trixler (2013) give none. |
| [`stress_test.py`](file:///Users/saptarshighosh/Documents/chirality-project/stress_test.py) | Domain-wall motion vs. vent racemization. |

---

## 6. Briefing Notes (if contacting Dr. Chris McKay, NASA Ames)

Everything below is a **model result** from this project's scripts. None of it has been validated against experimental or mission data.

1. **Titan at 94 K (model):**
   * In the model, classical chemistry with $E_a \ge 40$ kJ/mol is effectively frozen at 94 K.
   * With **assumed** tunnelling rates, the model's crossover temperature came out at 78–95 K. The pre-stated prediction was 120–240 K, so that prediction **failed**.
   * If the bias scales as $1/T$, it is about 2.9× larger at 94 K than at 273 K (simple ratio, not a measurement).
   * Lunine et al. (1999) wrote on abiotic homochirality on Titan (§1.3). It has not been read yet; read it before raising Titan with McKay.
   * The model estimates very slow racemization of heavy stereocenters at Titan surface temperature. This is a model estimate, not "strictly zero"; the same script prints a much shorter half-life (~2.5 Myr) for the proton-transfer channel.

2. **Open vs. closed (model + literature):**
   * The literature (Blanco et al. 2011; Higgs & Blackmond 2025) says closed systems return to racemic and open ones can stay non-racemic. The project's model uses this to suggest why Bennu/Ryugu amino acids are racemic and why driven ocean worlds might differ. This is a hypothesis.

3. **What has NOT been validated:**
   * **Yu et al. (2024) benchmark: withdrawn.** The dipeptide ee values used were not in the paper, so the R² fit is void.
   * **Hein et al. (2011) time series: withdrawn.** The 0–72 h ee series could not be found in the paper, so that R² fit is void.
   * No rover or other in situ mission data were used anywhere in the project.
