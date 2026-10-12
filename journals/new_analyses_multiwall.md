# multiwall: which layer's hand wins when several walls are past onset (2026-10-12)

Script `multiwall.py` (predictions committed before any run, 616811c), log `multiwall.log`. Sections 3, 1, 2 and 4 were run in that order, one at a time, with `nice -n 10` (1 s, 82 s, 6 s and 11 s). The log has two extra tables marked "not pre-registered": a domain-thickness check (section 2) and a cascade time ×5 (section 4).

## Method in plain words

sweep3 left one question open. When a reversed wall's deep layer is *also* layered, and its own wall is also past the Eq. 20 onset ("recursion", 53–72 % of draws on the moons), which hand fills the ocean? sweep3 counted these draws as no hand. Its bound (iii) counted them as wins. Neither choice was simulated.

1. **Simulation.** I used pin_onset's radial model (reaction plus depth-dependent vertical mixing, no flux at the floor or the ice). It is solved implicitly on a 0.002 grid with a 0.1 time step. I stacked 3–6 layers of set hands, so 2–5 walls. Mixing rises smoothly with depth, so each layer has its own Eq. 20 strength x.
   - **Section 1.** Every wall is past onset (x > 1). There are 336 runs: 4 layer counts × 6 hand patterns × 7 x patterns × 2 layer spacings.
   - **Section 2.** The onset depth r* falls inside the stack, using sweep3's own geometry (x = r/r*): 108 runs.
   - **Section 3.** Halving the grid and time step changes the end time by −1.8 %.

   For each run I recorded which hand filled the shell and when. I then tested three candidate rules (deepest layer, largest volume, top layer) and three time rules:
   - **T_deep:** the lone transit of the deepest wall (the "max" rule);
   - **T_sum:** the walls taking turns;
   - **T_conv:** sweep3's convention.
2. **Wrapper (section 4).** I imported sweep3 and sweep2 without editing them. I called `sweep3.outcome3` and reclassified only its recursion draws with the measured rule.
   - **Who wins.** Walls are past onset down to r* = 4H/ln ρ. The layer containing r* fills the ocean; that is the bottom layer if every wall is past onset. It holds the fraction f = min(d_w, H − K·d_w)/H, where K is the number of walls past onset.
   - **Which hand.** The deciding layer's hand is PVED-set if F·Δ(f) ≥ 2 on that layer's own N. Otherwise it is chance-set.
   - **Timing.** The deepest rising wall's transit to the ice comes first. It is followed by the longest of three times: the reaction time, the sideways mixing time, and the sink time of the walls below r*. Both transits use a closed form at the drawn (constant) D_z. The draw counts only if all of this fits in τ.
   - **Regression.** With "recursion = no hand", the wrapper reproduces every sweep3.log three-way share exactly. K ≥ 2 matches sweep3's recursion mask with 0 mismatched draws.

## Configurations and outcomes (± one standard error)

**Section 1: every wall past onset (R = 16, H = 1), 336 runs, 0 unfinished.**

| Group | Runs | Deepest hand fills | T_end/T_deep | T_end/T_sum | T_end/T_conv |
|---|---|---|---|---|---|
| 3 layers | 84 | 100.0 ± 0.0 % | 0.846 ± 0.027 | 0.774 ± 0.028 | 1.50 ± 0.39 |
| 4 layers | 84 | 100.0 ± 0.0 % | 0.728 ± 0.034 | 0.624 ± 0.036 | 1.36 ± 0.31 |
| 5 layers | 84 | 100.0 ± 0.0 % | 0.634 ± 0.042 | 0.483 ± 0.038 | 1.21 ± 0.22 |
| 6 layers | 84 | 100.0 ± 0.0 % | 0.756 ± 0.033 | 0.633 ± 0.036 | 1.34 ± 0.24 |
| hands alternating | 56 | 100.0 ± 0.0 % | 0.517 ± 0.050 | 0.275 ± 0.035 | 1.19 ± 0.35 |
| only bottom differs | 56 | 100.0 ± 0.0 % | 0.908 ± 0.019 | 0.908 ± 0.019 | 2.75 ± 0.58 |
| only top differs | 56 | 100.0 ± 0.0 % | 0.771 ± 0.022 | 0.771 ± 0.022 | 0.44 ± 0.08 |
| sandwich | 56 | 100.0 ± 0.0 % | 0.780 ± 0.046 | 0.606 ± 0.041 | 1.15 ± 0.17 |
| random hands (2 sets) | 112 | 100.0 ± 0.0 % | 0.836 ± 0.031 / 0.634 ± 0.055 | 0.691 / 0.520 | 1.16 / 1.42 |
| x uniform 1.1 / 1.5 / 3 / 10 | 48 each | 100.0 ± 0.0 % | 0.79 / 0.83 / 0.85 / 0.72 (± ≤ 0.05) | 0.65 / 0.68 / 0.69 / 0.60 | 0.52 / 0.55 / 0.60 / 0.90 |
| x rising with depth | 48 | 100.0 ± 0.0 % | 0.442 ± 0.043 | 0.368 ± 0.045 | 5.21 ± 0.81 |
| x falling with depth | 48 | 100.0 ± 0.0 % | 0.923 ± 0.035 | 0.841 ± 0.037 | 0.19 ± 0.01 |
| x random | 48 | 100.0 ± 0.0 % | 0.630 ± 0.052 | 0.558 ± 0.053 | 1.51 ± 0.18 |
| **all** | 336 | **100.0 ± 0.0 %** | median 0.927, range 0.007–1.028 | | median 0.645, range 0.007–17.4 |

The other two hand rules do worse: the largest-volume rule picked the winner in 52.1 ± 2.7 % of runs and the top-layer rule in 41.7 ± 2.7 %.

**Section 2: onset inside the stack (R = 2, x = r/r*), 108 runs, 0 unfinished.**

| Quantity | Value |
|---|---|
| Layer containing r* fills | 94/108 = 87.0 ± 3.2 % |
| ... when that domain is ≥ 5 wall widths thick (post hoc) | 63/63 |
| ... when it is ≥ 4 wall widths thick (post hoc) | 69/71 |
| T_end / max(T_up, T_down) | median 0.948, max 1.182, mean 0.758 ± 0.033 |
| T_end / T_c (wrapper closed form) | median 1.196, max 5.05, mean 1.30 ± 0.09 |
| T_conv > T_end | 75.0 % of runs |

## The derived rule and how well it fits

**Who wins.** The hand of the layer that contains the onset depth r* fills the ocean. When every wall is past onset, that is the deepest layer.

**Why.** Walls never cross. Each moves at a speed set only by its depth, so a domain between two walls is squeezed where the speed drops (near the ice, where mixing is weaker). It vanishes when it is about one wall width thick.

**How well it fits.**
- It is exact (336/336) when all walls are past onset.
- With onset inside the stack, it holds whenever the deciding domain is at least about 5 wall widths thick (63/63). All 14 failures were domains 2–4.4 wall widths thick.
- Real layers are kilometres thick and walls are metres wide, so the thin-domain failures do not occur in the sweep.
- Volume and the top layer do not predict the winner.

**How long it takes.** The lone transit of the deepest rising wall (the max rule) is an *upper* bound on the time: never exceeded by more than 3 % (1.028), and tight (0.9–1.0) when one wall must sweep the whole shell. Runs end earlier, down to 0.007 of it, when opposite domains are squeezed out on the way. Letting the walls take turns (the sum rule) is always worse (RMS log error 0.91 vs 0.68).

**Caveat: the timing formula.** The wrapper's constant-D_z closed form runs short of the measured time by a median 1.2× (the T_end/T_c ratio above), and by up to 5×. This happens when r* is deep, where mixing varies most along the path. That is why section 4 also reports the cascade time ×10 (pre-registered) and ×5 (post hoc).

## Ocean worlds: sweep3 main, sweep3 upper bound, multiwall

P_pved %, no lid / final (lid). Lid-free worlds have one number. "Multiwall" is the measured rule; "×10" is the same with a cascade 10× slower. The three-way split is P_pved / P_chance / P_none.

| World | sweep3 main | sweep3 bound (iii) | **multiwall** | multiwall ×10 | multiwall split, no lid | multiwall split, final |
|---|---|---|---|---|---|---|
| Early Earth | 67.0 | 85.6 | **99.9** | 99.7 | 99.9 / 0.0 / 0.1 | same |
| Umbriel | 20.3 | 43.4 | **53.7** | 43.6 | 53.7 / 14.9 / 31.4 | same |
| Oberon | 19.3 | 41.0 | **51.4** | 40.8 | 51.4 / 11.7 / 36.9 | same |
| Titania | 19.5 | 41.0 | **51.3** | 40.9 | 51.3 / 11.6 / 37.1 | same |
| Ariel | 16.9 | 41.0 | **49.6** | 39.7 | 49.6 / 13.5 / 36.9 | same |
| Ganymede | 20.2 | 33.1 | **44.1** | 33.5 | 44.1 / 3.4 / 52.5 | same |
| Callisto | 19.3 | 32.9 | **43.5** | 33.0 | 43.5 / 3.8 / 52.7 | same |
| Pluto | 14.2 | 33.5 | **41.1** | 30.6 | 41.1 / 6.1 / 52.8 | same |
| Triton | 12.4 | 30.3 | **37.0** | 26.5 | 37.0 / 4.7 / 58.3 | same |
| Titan | 14.1 | 27.5 | **35.8** | 25.5 | 35.8 / 2.5 / 61.7 | same |
| Europa | 13.9 / 10.5 | 31.8 / 24.1 | **39.6 / 30.1** | 29.0 / 22.0 | 39.6 / 4.9 / 55.5 | 30.1 / 3.7 / 66.2 |
| Dione | 12.7 / 7.9 | 36.6 / 22.8 | **42.0 / 26.2** | 32.2 / 20.1 | 42.0 / 10.5 / 47.5 | 26.2 / 6.6 / 67.2 |
| Miranda | 13.7 / 4.5 | 27.7 / 9.2 | **27.1 / 9.0** | 21.3 / 7.0 | 27.1 / 8.4 / 64.5 | 9.0 / 2.7 / 88.3 |
| Mimas | 7.2 / 2.3 | 15.3 / 4.9 | **17.3 / 5.6** | 11.6 / 3.7 | 17.3 / 6.1 / 76.6 | 5.6 / 2.0 / 92.5 |
| Enceladus | 6.6 / 0.0 | 23.9 / 0.0 | **28.3 / 0.0** | 19.9 / 0.0 | 28.3 / 8.8 / 63.0 | 0.0 / 0.0 / 100.0 |

**Other variants.**
- The sum rule (walls take turns) gives values close to ×10: Ganymede 36.0, Europa 30.4 / 23.1, Enceladus 20.1 / 0.0.
- ×5 (post hoc) falls in between: Ganymede 36.4, Europa 31.9 / 24.2, Enceladus 22.4 / 0.0.
- Seeds 1–3 reproduce every share within 0.38 points.

**Why multiwall beats the "upper bound".** On the moons, 93–100 % of recursion draws are decided by the bottom layer (Miranda 56 %). In those draws every wall rises, nothing has to sink, and the wait is shorter than the full-depth sink that bound (iii) charged. Of the recursion draws, 25–69 % finish in time on the moons and 99.9 % on Earth.

**Why P_chance grows.** The deciding layer is the thin remnant at the bottom, so f ≤ d_w/H. Its hand is more often set by chance: P_chance reaches 8–15 % on the small moons and the Uranian moons.

## Predictions: 5 of 16 held

| # | Prediction (short) | Outcome |
|---|---|---|
| P1 | Deepest hand fills 336/336; volume and top rules each fail ≥ 25 % | held (volume fails 47.9 %, top 58.3 %) |
| P2 | Deep hand never loses > 0.01 of volume | **failed**: 0.083; the same-hand upper domains are squeezed into the ice |
| P3 | T_end/T_deep in [0.5, 1.15], sum-rule error > 2× max-rule | **failed**: min 0.007 (squeezed domains); error ratio 1.34. Max ≤ 1.03 held |
| P4 | T_end/T_conv in [0.1, 3] | **failed**: 0.007–17.4 (> 3 only for x rising with depth or random) |
| P5 | Layer containing r* fills 108/108 | **failed**: 94/108; every failure is a domain ≤ 4.4 wall widths thick |
| P6 | T_end/max(T_up, T_down) in [0.6, 1.3] | **failed**: min 0.003; max 1.18 held |
| P7 | T_end/T_c in [0.2, 5] | **failed** narrowly: max 5.05, min 0.006 |
| P8 | T_conv > T_end in ≥ 80 % | **failed**: 75.0 % |
| P9 | Halving dt, dz changes T_end < 2 % | held (−1.83 %) |
| P10 | Regression exact; mask mismatch ≤ 10 draws | held (exact; 0) |
| P11 | sweep3 main ≤ multiwall ≤ bound (iii) + 0.5 | **failed**: above the bound in 14 of 15 worlds (Earth +14.3, Ganymede +11.0) |
| P12 | Moons gain ≥ 3; Earth 75–86 % | **failed**: gains ≥ 10.1 held, Earth 99.9 |
| P13 | ×10 costs ≤ 3 points | **failed**: up to 10.6 (Ganymede, Oberon); Earth 0.2 |
| P14 | P_chance ≤ 5 | **failed**: up to 14.9 (Umbriel) in 9 worlds |
| P15 | Final: Ganymede > Europa > Enceladus | held (44.1 > 30.1 > 0.0) |
| P16 | Seeds within 0.5 points | held (≤ 0.38) |

The hand rule itself held (P1). P5 failed only for domains thinner than real layers. The failures are mostly about timing: I predicted the walls would move together at the lone-wall speed, but walls also squeeze domains out. The single-wall transit turned out to be a bound, not an estimate. The ocean-world predictions failed in the optimistic direction (P11, P12): nothing has to sink when every wall rises. P13 and P14 failed because the deciding bottom layer is thin and the constant-D_z timing is loose.

## What the paper should now use

The paper should quote the multiwall numbers, with the cascade-×10 values as the low end of each range:
- Early Earth: about 100 %.
- Ganymede: 34–44 %.
- Europa: 29–40 % without a lid and 22–30 % with one.
- Enceladus: 20–28 % without a lid and 0 % with one.
- The four Uranian moons: 40–54 %.
- Callisto 33–44 %, Titan 26–36 %, Triton 27–37 %, Pluto 31–41 %.
- Dione 20–26 %, Miranda 7–9 % and Mimas 4–6 %, all with a lid.

These replace sweep3 main (Earth 67 %, Ganymede 20 %, Europa 14 / 10.5 %, Enceladus 7 / 0 %). sweep3 counted a stack of reversed walls as no hand. The simulation shows instead that every wall past onset rises and the deepest layer's hand fills the ocean. It did so in all 336 runs, and in every run with onset inside the stack where the deciding layer was at least 5 wall widths thick.

The low end of each range allows for the one loose part, the travel time. The constant-mixing formula ran up to 5× short of the simulated time, so a 10× slower cascade costs up to 11 points. The ranking does not change:
- Among confirmed oceans, Ganymede is highest, Europa next and Enceladus last.
- Lid-free moons now score 26–54 %.
- About one draw in ten on the small and Uranian moons is set by chance, because the deciding layer is a thin bottom remnant.

The final Europa range, 22–30 %, happens to match the draft's current 22–26 %. Ganymede (34–44 %) and Earth (~100 %) move up.

## Proposed replacement sentences (not applied)

These replace the draft's figures (22 / 29 / 75 %, Ganymede 28 %, Europa 22–26 %, Enceladus 0–18 %). They also supersede the text proposed in `journals/new_analyses_sweep3.md`.

**Abstract.** Replace from "Varying all uncertain inputs together, ..." through "... (0–18% when a lid forms)." with:

> Varying all uncertain inputs together, the PVED-favoured hand is selected in nearly all plausible input sets for early Earth, 34–44% for Ganymede, 22–30% for Europa and 0% for Enceladus when a salt-stratified lid forms (20–28% without one). Vertical mixing is the main unknown. Weak mixing freezes the ocean into stacked layers of opposite hand. A stratified top layer pushes each boundary up to the ice, and simulations of stacked boundaries show that the hand of the deepest such layer then fills the ocean. Among confirmed oceans Ganymede scores highest and Enceladus lowest.

**Section 3.1.** Replace its closing sentence, "Only ocean-scale water that is mixed within the selection window can pool enough molecules (Figure 3).", with:

> Only ocean-scale water that is mixed within the selection window can pool enough molecules (Figure 3). Real oceans do not stay mixed. They break into layers, and the layer that decides the outcome holds only part of the ocean. With that included, Sections 3.2–3.5 find that the favoured hand wins in nearly all early-Earth cases and in 0–44% for the confirmed icy moons.

**Section 3.3.** Replace the first paragraph ("Varying every input of Table 1 at once ... (`uncertainty.log`).") with:

> Varying every input of Table 1 at once (Section 2.6), the PVED-favoured hand is selected in 20–28% of plausible cases for Enceladus and 29–40% for Europa before the lid test, and in nearly all cases (≥ 99.7%) for early Earth's open ocean (`multiwall.log`). The lower value of each range assumes the stacked boundaries rise 10× slower than the constant-mixing estimate. If layers never merged (the flat-box bound), the moons would fall to 3% and 7% and early Earth to 40% (`uncertainty.log`).

The later 3.3 paragraphs (sensitivity bars, the bias and swing values, the 87 % at the thermocline) come from `uncertainty.log`, which uses the pre-layering pipeline. Either recompute them or label them "before the layering correction".

**Table 3 (Section 3.5).** Replace the "Success share" column with multiwall, final, as "main (cascade ×10)":

| World | Success share |
|---|---|
| Ganymede | 44% (34%) |
| Europa | 30% (22%) |
| Enceladus | 0% (0%) |
| Callisto | 44% (33%) |
| Triton | 37% (27%) |
| Mimas | 6% (4%) |
| Titan | 36% (26%) |
| Pluto | 41% (31%) |
| Umbriel, Ariel, Oberon, Titania | 50–54% (40–44%) |
| Dione | 26% (20%) |
| Miranda | 9% (7%) |

The "Across ocean depths" column was not recomputed by this script. Drop it or rerun the depth sweep through the wrapper.

Replace "Every lid-free world scores between 24% and 39%." with:

> Every lid-free world scores between 36% and 54% (26–44% with a 10× slower cascade).

Replace "Of the three confirmed oceans, Ganymede (28%) and Europa (22–26%) are the most favourable." with:

> Of the three confirmed oceans, Ganymede (34–44%) is the most favourable, then Europa (22–30%).

The Juno (26 %) and Kang et al. (4–18 %) values in the same section are also pre-layering and need the same rerun.

**Conclusions.** Replace bullets 2 and 3 with:

> - Across plausible inputs, the PVED-favoured hand is selected in nearly all early-Earth cases and in 20–28% (Enceladus) and 29–40% (Europa) before the lid test. Vertical mixing is the main unknown.
> - A stratified top layer pushes each boundary between hands up to the ice, and the deepest layer's hand fills the ocean. Among confirmed oceans, Ganymede (34–44%) scores highest, then Europa (22–30%), and Enceladus lowest (0% when a lid forms).
