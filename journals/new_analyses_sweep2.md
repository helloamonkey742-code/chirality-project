# sweep2: literature PVED range, rate bias, D_z floor, racemization, feedstock/product, depth pinning (2026-10-12)

Script `sweep2.py`, log `sweep2.log` (200,000 draws per world, seed 0; seeds 1-3 checked). Predictions were committed
before the run (commit d54a387). The all-off step reproduces `uncertainty.sweep` bit for bit and every logged share
(`ocean_worlds.log`, `strat_sweep.log`, `assumed_rules.log`) within 0.3 points. Nothing in the paper drafts was edited.

## 1 PVED values from the literature

g = |ΔE_pv| / (R · 273 K); 1 aeV = 9.6485×10⁻¹⁴ J/mol. "Abstract only" means the value or sign was read from the abstract, not the full text.

| # | Molecule / method | \|ΔE_pv\| | g at 273 K | Sign | Source | How checked |
|---|---|---|---|---|---|---|
| 1 | L-Ala (aqueous conformer), Gly, peptides; first-generation ab initio | ~10⁻¹⁴ J/mol (0.10 aeV) | 4.4×10⁻¹⁸ | L lower | Mason & Tranter 1985, Proc R Soc A, 10.1098/rspa.1985.0003 | abstract only |
| 2 | Ala, Val, Ser, Asp; first generation | ~10⁻¹⁴ J/mol | 4.4×10⁻¹⁸ | L lower in aqueous conformations; conformer-sensitive | Tranter 1985, Mol Phys, 10.1080/00268978500102741 | Crossref + abstract only |
| 3 | Alanine (Bakasov 1998 / Berger & Quack 2000, as tabulated) | 6.2 aeV = 6.0×10⁻¹³ J/mol | 2.64×10⁻¹⁶ | undetermined ("strongly conformer dependent") | Quack, Seyfang & Wichmann 2022, Chem Sci, 10.1039/d2sc01323a | full text (PMC9491092) |
| 4 | Ala gas, 13 conformers | not in abstract | – | L lower in 7 of 13 | Laerdahl, Wesendrup & Schwerdtfeger 2000, ChemPhysChem 1:60 | abstract |
| 5 | Ala | not in abstract | – | L stabilisation claims "unjustified" | Berger & Quack 2000, ChemPhysChem 1:57 | abstract |
| 6 | Amino acids, second generation | ~10× first generation | – | gas L-Ala negative; sign in solution undetermined | MacDermott et al. 2009, OLEB, 10.1007/s11084-009-9161-x | abstract |
| 7 | α-methyl amino acids | not stated | – | L negative | MacDermott 2009, 10.1007/s11084-009-9162-9 | abstract |
| 8 | Amino acids, relativistic | not stated | – | no clear L preference | Wesendrup et al. 2003, JPCA, 10.1021/jp022568v | abstract |
| 9 | Generic light chiral molecule | ~10⁻¹³ J/mol (1.0 aeV) | 4.4×10⁻¹⁷ | – | Lente 2006, JPCA, 10.1021/jp0658344 | abstract only |
| 10 | Hydrated D-glyceraldehyde | not in abstract | – | D lower | Tranter 1986, JCS Chem Commun, 10.1039/c39860000060 | abstract |
| 11 | Deoxyribose, ribose, other pentoses | not in abstract | – | conformer-dependent (D for C2-endo, L for C3-endo) | Tranter et al. 1992, Proc R Soc A, 10.1098/rspa.1992.0037 | abstract |
| 12 | Gly, Ala, Ser basis-set study | – | – | basis-set dependent | Kikuchi & Wang 1990, BCSJ, 10.1246/bcsj.63.2751 | abstract |

New prior: g log-uniform 4.4×10⁻¹⁸ to 2.64×10⁻¹⁶ (0.10–6.2 aeV), replacing 4×10⁻¹⁹ to 4×10⁻¹⁶. The lower end (row 1) is the smallest published amino-acid magnitude, not a physical floor; averaging over aqueous conformers could give less, and the sign is unsettled. The upper end (row 3) is the largest verified amino-acid value, so g_max is the alanine value. No sugar magnitude was accessible. The full text of Quack et al. (2022) gives "aeV to feV" for light-atom molecules. Its abstract's 10⁻¹¹–10⁻¹⁰ J/mol figure is for heavier molecules, so §1.1 of the draft should not apply it to amino acids.

Racemization: free alanine in water at 3 °C racemizes at k = 0.046×10⁻⁵ yr⁻¹ (Bada & Schroeder 1975, 10.1007/bf00592179, as tabulated in Steen et al. 2013 PLoS ONE Table 3, 10.1371/journal.pone.0071648, read in full text). Taken to 273 K with Ea = 93 kJ/mol (a fit to sediment-bound data), this gives k_rac = 9.13×10⁻¹⁵ s⁻¹. The sensitivity case uses Asp: 2.73×10⁻¹⁴ s⁻¹ (Ea 110 kJ/mol). All worlds are taken at 273 K (ASSUMED).

## 2 Ladder (cumulative), %

Steps: 0 = as published; a = PVED range; b = rate bias η (log-uniform 0.1–10); c = D_z floor 10⁻⁹ m²/s; d = racemization (k₂c > 2k_rac); e = product/feedstock p/a (log-uniform 0.01–1); f = depth-dependent D_z pinning (bottom/top ratio ρ log-uniform 1–100, Eq. 20 at the top-layer wall). Each cell is no-lid / final (after the salt lid); where only one number is shown, no lid is possible.

| World | 0 | a | b | c | d | e | f |
|---|---|---|---|---|---|---|---|
| Early Earth | 74.7 | 74.7 | 74.7 | 74.7 | 74.7 | 74.7 | **39.9** |
| Umbriel | 38.7 | 42.2 | 41.6 | 48.5 | 44.9 | 42.3 | **14.4** |
| Ganymede | 28.3 | 29.6 | 29.3 | 34.3 | 30.4 | 29.9 | **14.0** |
| Titania | 36.4 | 39.5 | 38.9 | 45.4 | 41.6 | 39.5 | **13.7** |
| Oberon | 36.5 | 39.5 | 38.9 | 45.5 | 41.7 | 39.6 | **13.5** |
| Callisto | 28.5 | 29.8 | 29.5 | 34.5 | 30.6 | 30.0 | **13.4** |
| Ariel | 37.0 | 40.2 | 39.7 | 46.2 | 42.7 | 40.2 | **11.7** |
| Pluto | 30.2 | 32.3 | 31.8 | 37.2 | 33.5 | 32.3 | **9.5** |
| Titan | 24.3 | 25.4 | 25.1 | 29.3 | 25.9 | 25.5 | **9.4** |
| Triton | 27.7 | 29.5 | 29.1 | 34.0 | 30.4 | 29.5 | **8.4** |
| Europa | 28.6 / 21.7 | 30.4 / 23.0 | 30.0 / 22.7 | 35.0 / 26.6 | 31.3 / 23.8 | 30.4 / 23.1 | **9.3 / 7.1** |
| Dione | 33.6 / 21.0 | 36.6 / 22.9 | 36.0 / 22.5 | 42.1 / 26.3 | 38.7 / 24.2 | 36.5 / 22.8 | **8.9 / 5.5** |
| Miranda | 26.9 / 8.9 | 29.7 / 9.8 | 29.3 / 9.7 | 34.3 / 11.3 | 32.8 / 10.8 | 30.1 / 10.0 | **12.4 / 4.1** |
| Mimas | 14.7 / 4.7 | 16.3 / 5.2 | 16.0 / 5.1 | 18.7 / 6.0 | 18.7 / 6.0 | 16.9 / 5.4 | **5.4 / 1.7** |
| Enceladus | 21.8 / 0.0 | 24.1 / 0.0 | 23.3 / 0.0 | 27.2 / 0.0 | 26.4 / 0.0 | 24.5 / 0.0 | **3.9 / 0.0** |

Sensitivities at step f (no lid / final):
- **Eq. 20 threshold ×5.** Earth 40.0. Icy moons without a lid: 14.1 (Titan) to 17.9 (Umbriel). Europa 16.3 / 12.3, Dione 21.7 / 13.5, Miranda 30.1 / 10.0, Mimas 16.0 / 5.1, Enceladus 16.9 / 0.0.
- **Asp k_rac instead of Ala.** No world moves by more than 0.7 points.

Seeds 1–3: the largest difference from seed 0 is 0.29 points.

## 3 Three-way outcome split and unlimited bias (step f), %

P_pved: the PVED-favoured hand wins. P_chance: one hand forms, but the bias is too weak to choose it (F·Δ < 2). P_none: no single hand forms (unfinished, racemized, layers not merged, wall pinned, or salt lid). Unlimited bias = P_pved + P_chance.

| World | No lid: P_pved / P_chance / P_none | Final: P_pved / P_chance / P_none | Unlimited (final) | Unlimited at step 0 (no lid) |
|---|---|---|---|---|
| Early Earth | 39.9 / 0.0 / 60.1 | 39.9 / 0.0 / 60.1 | 39.9 | 74.7 |
| Umbriel | 14.4 / 1.4 / 84.2 | same | 15.8 | 46.9 |
| Ganymede | 14.0 / 0.3 / 85.7 | same | 14.3 | 29.9 |
| Titania | 13.7 / 1.0 / 85.3 | same | 14.7 | 42.7 |
| Oberon | 13.5 / 1.0 / 85.4 | same | 14.6 | 42.8 |
| Callisto | 13.4 / 0.3 / 86.4 | same | 13.6 | 30.2 |
| Ariel | 11.7 / 1.1 / 87.2 | same | 12.8 | 44.6 |
| Pluto | 9.5 / 0.4 / 90.1 | same | 9.9 | 33.5 |
| Titan | 9.4 / 0.2 / 90.4 | same | 9.6 | 25.6 |
| Triton | 8.4 / 0.3 / 91.3 | same | 8.7 | 30.4 |
| Europa | 9.3 / 0.3 / 90.4 | 7.1 / 0.2 / 92.7 | 7.3 | 31.3 |
| Dione | 8.9 / 0.8 / 90.3 | 5.5 / 0.5 / 93.9 | 6.1 | 40.2 |
| Miranda | 12.4 / 2.9 / 84.7 | 4.1 / 0.9 / 95.0 | 5.0 | 35.5 |
| Mimas | 5.4 / 1.4 / 93.2 | 1.7 / 0.5 / 97.8 | 2.2 | 19.0 |
| Enceladus | 3.9 / 2.0 / 94.1 | 0.0 / 0.0 / 100.0 | 0.0 | 28.6 |

## 4 Scale rule (Eq. 13 at kτ = 10)

N = 2/(√π g² √(kτ)). The 10 km³ lake is taken at 1 mM (ASSUMED). A shortfall below 1 means the volume is enough.

| g | N needed | Cube side at 1 µM / 1 mM / 0.1 M | Shortfall: 1 m³ pond at 1 M / 1 km³ lake at 1 mM / 10 km³ lake at 1 mM |
|---|---|---|---|
| 10⁻¹⁷ (paper) | 3.57×10³³ | 18.1 km / 1.81 km / 390 m | 5.9×10⁶ / 5.9 / 0.59 |
| g_min 4.4×10⁻¹⁸ | 1.84×10³⁴ | 31.3 km / 3.13 km / 674 m | 3.1×10⁷ / 31 / 3.1 |
| g_max 2.64×10⁻¹⁶ | 5.12×10³⁰ | 2.04 km / 204 m / 44 m | 8.5×10³ / 0.0085 / 0.00085 |
| alanine 2.64×10⁻¹⁶ (= g_max) | 5.12×10³⁰ | same | same |

The first-generation alanine value equals g_min.

## 5 Predictions: 7 of 9 held

| # | Prediction | Result |
|---|---|---|
| P1 | (a) Enceladus +1 to +3, Europa +0.3 to +1.5, Earth ±0.3 | **FAILED**: Enceladus +2.3 and Earth 0.0 held, but Europa rose +1.8 |
| P2 | (b) η moves no world by > 1 point | held (largest change −0.8) |
| P3 | (c) icy moons ×1.12–1.21; Earth unchanged | held (×1.16–1.18; 7/6 = 1.167) |
| P4 | (d) Earth unchanged; 4-Gyr worlds −5 to −15 % relative; Enceladus < 6 %; Mimas < 3 % | held (7–12 %, 2.9 %, 0 %) |
| P5 | (e) Enceladus −1 to −4, Europa −0.3 to −2, Earth < 1 | held (−1.9, −0.9, 0.0) |
| P6 | (f) largest step everywhere, with ranges for each world | **FAILED** narrowly: Enceladus no-lid 3.9 (predicted ≥ 4) and Umbriel 14.4 (predicted ≤ 14); all other ranges held |
| P7 | (g) P_chance: Europa < 3, Earth < 1, Enceladus no-lid 1–4 | held (0.2, 0.0, 2.0) |
| P8 | Ganymede > Europa > Enceladus = 0 | held (14.0 > 7.1 > 0.0) |
| P9 | Seeds 1–3 within 0.5 points | held (0.29 max) |

## 6 What changes for the paper

**In plain terms:**
- **The weak-force bias is now taken from published calculations.** That range is a little narrower and higher than the one the paper used. It helps the moons by 1–3 points, and Earth is unchanged.
- **The other single-input changes each move a moon by at most 7 points and mostly cancel:**
  - letting the bias act on reaction rates more or less strongly than the energy suggests;
  - removing implausibly weak vertical mixing;
  - letting molecules lose their handedness slowly in water;
  - letting the product be scarcer than the feedstock.
- **One change decides everything: vertical mixing that varies with depth.**
  - A boundary between the two hands stops moving once mixing changes with depth faster than the boundary's own curvature can push it.
  - With a plausible bottom-to-top mixing contrast of 1–100, layered oceans almost never merge into one hand.
  - In Earth's thin ocean, any contrast above 0.2 % does this.
- **New shares** (fraction of plausible cases in which the PVED-favoured hand wins):

  | | Before | Now |
  |---|---|---|
  | Early Earth | 75 % | 40 % (the paper's old "flat-box" bound) |
  | Europa | 22–29 % | 7–9 % |
  | Enceladus | 22 % | 0–4 % |
  | Best confirmed ocean | – | Ganymede, 14 % |

- **The size of the bias barely matters any more.** Even an unlimited bias adds at most 3 points anywhere. Whether the weak force can pick the hand now depends almost entirely on whether one hand can form at all, not on how strong the push is.
- **This conclusion rests on two assumptions that have not been measured:**
  - the 1–100 range for the depth contrast;
  - the sharp pinning criterion (Eq. 20).

  A five times looser criterion lifts the moons back to 14–30 % without a lid. The headline should therefore be reported with the pinning step as a stated assumption and the ×5 case as its range.
- **The size estimate widens.** A pool of about 5×10³⁰ to 2×10³⁴ molecules is needed, depending on which published alanine value is used. At the larger value, a 1 km³ lake at 1 mM is enough; at the smaller, even 10 km³ falls 3× short.

## 7 Proposed replacement sentences (NOT applied; `journals/ija_draft.md` is being edited by another agent)

**Abstract**
- Replace "About 10³³–10³⁴ molecules must react together as one well-mixed pool. That needs kilometre-scale volumes of well-mixed water. Ponds and rock pores are many orders of magnitude too small" with:
  > "Between about 5×10³⁰ and 2×10³⁴ molecules must react together as one well-mixed pool, depending on which published alanine PVED is used. That needs a well-mixed lake or ocean, from about 1 to more than 10 km³ at millimolar concentration. Ponds and rock pores are at least four orders of magnitude too small."
- Replace "the PVED-favoured hand is selected in 22% (Enceladus), 29% (Europa) and 75% (early Earth) of plausible input sets" with:
  > "the PVED-favoured hand is selected in 0–4% (Enceladus), 7–9% (Europa) and 40% (early Earth) of plausible input sets once vertical mixing is allowed to vary with depth (14–30% for the moons and 40% for Earth if the pinning criterion is five times looser)."
- Replace "Ganymede (28%) and Europa (22–26%) score highest among confirmed oceans, and Enceladus lowest (0–18% when a lid forms)" with:
  > "Ganymede (14%) scores highest among confirmed oceans, then Europa (7–9%) and Enceladus (0–4%)."
- Add:
  > "Even an unlimited bias adds at most 3 points: the share is set by whether a single hand can form, not by the size of the PVED."

**§2.1**
- Replace "g is drawn log-uniformly from 4×10⁻¹⁹ to 4×10⁻¹⁶ (`uncertainty.py: sweep()`), which is about 0.01 to 10 aeV per molecule. This range lies at or below the low end of the values given by Quack et al. (2022) for light-atom molecules, so it does not overstate the bias." with:
  > "In `sweep2.py`, g is drawn log-uniformly from 4.4×10⁻¹⁸ to 2.64×10⁻¹⁶ (0.10–6.2 aeV). The lower end is the first-generation amino-acid value of order 10⁻¹⁴ J/mol (Mason & Tranter, 1985); the upper end is the 6.2 aeV tabulated for alanine by Quack et al. (2022). Neither end is a bound: conformer averaging in water can give less, and the sign for amino acids in solution is undetermined (MacDermott et al., 2009)."
- Add after the transition-state paragraph:
  > "To allow for the transition-state difference, g is multiplied by a factor η drawn log-uniformly from 0.1 to 10 (an assumed range); this moves no world by more than 1 point."

**§3.1**
- Replace "For g = 10⁻¹⁷ (the middle of the modern PVED range)" with:
  > "For g = 10⁻¹⁷ (between the first- and second-generation amino-acid values)"
- Add after the cube-side sentence:
  > "Over the published alanine range the requirement runs from 5.1×10³⁰ molecules (g = 2.64×10⁻¹⁶; a 204 m cube at 1 mM, so a 1 km³ lake suffices) to 1.8×10³⁴ (g = 4.4×10⁻¹⁸; a 3.1 km cube at 1 mM, 3× more than a 10 km³ lake holds)."
- Replace "about 7 orders of magnitude too few" with:
  > "4 to 7 orders of magnitude too few, depending on g."

**§3.3**
- Replace the first sentence with:
  > "Varying every input at once, with the literature PVED range, a rate-bias factor, a D_z floor of 10⁻⁹ m² s⁻¹, racemization, separate feedstock and product, and depth-dependent vertical mixing (`sweep2.py`), the PVED-favoured hand is selected in 3.9% of plausible cases for Enceladus without a lid, 9.3% for Europa and 39.9% for early Earth's open ocean."
- Add:
  > "The depth dependence decides this. A wall between layers is pinned when ln(D_bottom/D_top) ≥ 4H/r (Eq. 20), and with a bottom/top contrast drawn from 1 to 100 that holds for most layered draws: every one on Earth, where the threshold contrast is 1.002, about 85% on Enceladus and 93% on Europa. With a five times looser threshold the moons rise to 14–30%, while Earth stays at 40%."
- Replace "Even an unlimited bias stops at 28.6% and 31.2%" with:
  > "Even an unlimited bias stops at 28.6% and 31.2% in the published model, and at 5.9% and 9.6% (Enceladus and Europa, no lid) and 39.9% (Earth) once the new effects are included. The PVED almost always chooses the hand when a single hand forms; whether one forms sets the share."

**Conclusions**
- Replace bullet 1's first sentence with:
  > "For the PVED, rather than chance, to control the outcome, between about 5×10³⁰ and 2×10³⁴ molecules must react together as one well-mixed pool, depending on the published alanine PVED used."
- Replace bullet 2 with:
  > "Across plausible inputs, including depth-dependent vertical mixing, the PVED-favoured hand is selected in 40% of early-Earth cases, 7–9% for Europa and 0–4% for Enceladus. Whether a single hand forms at all, not the size of the bias, sets these numbers, and the pinning of walls by depth-varying mixing is the main unknown."
- Replace "Ganymede (28%) and Europa (22–26%) score highest and Enceladus lowest (0–18% when a lid forms)" with:
  > "Ganymede (14%) scores highest among confirmed oceans, then Europa (7–9%) and Enceladus (0–4%)."

**Also flag (§1.1):** the draft applies Quack et al.'s 10⁻¹¹–10⁻¹⁰ J/mol to amino acids. That figure is for heavier molecules (full text), so it should not be used for amino acids.
