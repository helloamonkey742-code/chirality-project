# sweep3: counting reversed walls (2026-10-12)

Script `sweep3.py` (predictions committed before any run, 9c8e133), log `sweep3.log`. N = 200,000 draws per world, seeds 0-3, runtime 9 s.

## Method in plain words

sweep2's last step (f) said: if a layered ocean has a stratified top layer strong enough to pin the wall between hands (Eq. 20, x = ln(ρ)·r_w/(4H) ≥ 1), then no hand wins (P_none). `pin_onset.log` shows that is wrong above the threshold: the wall does not stay put, it *reverses* and rises, and the hand of the layer below fills the ocean instead. sweep3 reruns step f for every world with the same draws and recounts only those draws.

For each pinned draw (no salt lid, and the reaction finishes, k₂cτ ≥ 10):

1. **Deep layer.** I take the "deep layer" to be the whole bulk below the top-layer wall, depth H_b = H − d_w. If that bulk is not itself layered (H_b·√(D_h/D_z) ≤ L), it heals sideways and its share of molecules is H_b/H. If it is layered, its own top sublayer decides (share d_w/H), after its own walls sink in H_b·r_w/(2D_z). The share used is min(H_b, d_w)/H.
2. **Recursion.** If the deep layer is layered and its own top wall (at r_w − d_w) also passes Eq. 20, it reverses too and the question repeats. Counted as P_none, as specified, and tallied.
3. **Time.** The wall must rise out in t_rise = H·r_w/(2·D_z·(x − 1)), and the deep layer must then select its hand (the longest of 10/k, the sideways mixing time and its own sink time). If t_rise + t_sel > τ: P_none.
4. **Which hand.** F·Δ_deep ≥ 2 (same Δ formula as every other step, with the deep share): PVED-set. A hand forms but F·Δ_deep < 2: chance-set.

Regression: with the recount off, sweep3 reproduces sweep2 step f bit for bit, and the sweep2.log step e and f lines within 0.05 points (P1).

## Results: sweep2 step e, step f and sweep3

P_pved %, no lid / final (lid). Lid-free worlds have one number.

| World | sweep2 step e | sweep2 step f | sweep3 | gain over f |
|---|---|---|---|---|
| Early Earth | 74.7 | 39.9 | **67.0** | +27.1 |
| Umbriel | 42.3 | 14.4 | **20.3** | +5.9 |
| Ganymede | 29.9 | 14.0 | **20.2** | +6.2 |
| Titania | 39.5 | 13.7 | **19.5** | +5.8 |
| Callisto | 30.0 | 13.4 | **19.3** | +6.0 |
| Oberon | 39.6 | 13.5 | **19.3** | +5.8 |
| Ariel | 40.2 | 11.7 | **16.9** | +5.2 |
| Pluto | 32.3 | 9.5 | **14.2** | +4.7 |
| Titan | 25.5 | 9.4 | **14.1** | +4.7 |
| Triton | 29.5 | 8.4 | **12.4** | +4.0 |
| Europa | 30.4 / 23.1 | 9.3 / 7.1 | **13.9 / 10.5** | +4.6 / +3.5 |
| Dione | 36.5 / 22.8 | 8.9 / 5.5 | **12.7 / 7.9** | +3.8 / +2.4 |
| Miranda | 30.1 / 10.0 | 12.4 / 4.1 | **13.7 / 4.5** | +1.3 / +0.4 |
| Mimas | 16.9 / 5.4 | 5.4 / 1.7 | **7.2 / 2.3** | +1.7 / +0.6 |
| Enceladus | 24.5 / 0.0 | 3.9 / 0.0 | **6.6 / 0.0** | +2.8 / 0.0 |

## Three-way split (P_pved / P_chance / P_none, %)

| World | sweep3 no lid | sweep3 final | sweep2 step f no lid |
|---|---|---|---|
| Early Earth | 67.0 / 0.0 / 33.0 | same | 39.9 / 0.0 / 60.1 |
| Umbriel | 20.3 / 2.8 / 76.9 | same | 14.4 / 1.4 / 84.2 |
| Ganymede | 20.2 / 0.7 / 79.1 | same | 14.0 / 0.3 / 85.7 |
| Titania | 19.5 / 2.2 / 78.4 | same | 13.7 / 1.0 / 85.3 |
| Callisto | 19.3 / 0.7 / 79.9 | same | 13.4 / 0.3 / 86.4 |
| Oberon | 19.3 / 2.1 / 78.6 | same | 13.5 / 1.0 / 85.4 |
| Ariel | 16.9 / 2.2 / 80.9 | same | 11.7 / 1.1 / 87.2 |
| Pluto | 14.2 / 0.9 / 84.9 | same | 9.5 / 0.4 / 90.1 |
| Titan | 14.1 / 0.4 / 85.4 | same | 9.4 / 0.2 / 90.4 |
| Triton | 12.4 / 0.7 / 86.9 | same | 8.4 / 0.3 / 91.3 |
| Europa | 13.9 / 0.8 / 85.4 | 10.5 / 0.6 / 88.9 | 9.3 / 0.3 / 90.4 |
| Dione | 12.7 / 1.5 / 85.8 | 7.9 / 1.0 / 91.1 | 8.9 / 0.8 / 90.3 |
| Miranda | 13.7 / 3.2 / 83.1 | 4.5 / 1.0 / 94.5 | 12.4 / 2.9 / 84.7 |
| Mimas | 7.2 / 1.9 / 91.0 | 2.3 / 0.6 / 97.1 | 5.4 / 1.4 / 93.2 |
| Enceladus | 6.6 / 2.6 / 90.8 | 0.0 / 0.0 / 100.0 | 3.9 / 2.0 / 94.1 |

Chance-set outcomes stay small everywhere (at most 3.2 %). Most of the "no hand" share is not chance; it is draws where no single hand ever forms in time.

**What happens to reversed draws** (no lid, % of all draws): reversed 55-77 % on the moons and 60 % on Earth. Of those, the deep layer's own wall is also past onset (recursion) in 90.0 % (Umbriel) to 96.5 % (Miranda) on the moons and 54.9 % on Earth. Too late: at most 0.5 %. So the recount is limited almost entirely by recursion, not by the rise time.

## Sensitivities (P_pved %, no lid / final)

| World | main | (i) rise time ×10 | (ii) deep layer k > 2k_rac (Ala) | (iii) recursion counted (bound) | step e |
|---|---|---|---|---|---|
| Early Earth | 67.0 | 66.8 | 67.0 | 85.6 | 74.7 |
| Umbriel | 20.3 | 20.1 | 20.0 | 43.4 | 42.3 |
| Ganymede | 20.2 | 18.7 | 19.6 | 33.1 | 29.9 |
| Callisto | 19.3 | 18.0 | 18.8 | 32.9 | 30.0 |
| Titan | 14.1 | 12.7 | 13.7 | 27.5 | 25.5 |
| Europa | 13.9 / 10.5 | 13.1 / 9.9 | 13.5 / 10.2 | 31.8 / 24.1 | 30.4 / 23.1 |
| Dione | 12.7 / 7.9 | 12.5 / 7.8 | 12.4 / 7.8 | 36.6 / 22.8 | 36.5 / 22.8 |
| Miranda | 13.7 / 4.5 | 13.4 / 4.4 | 13.6 / 4.5 | 27.7 / 9.2 | 30.1 / 10.0 |
| Mimas | 7.2 / 2.3 | 6.4 / 2.1 | 7.2 / 2.3 | 15.3 / 4.9 | 16.9 / 5.4 |
| Enceladus | 6.6 / 0.0 | 6.1 / 0.0 | 6.6 / 0.0 | 23.9 / 0.0 | 24.5 / 0.0 |

All 15 worlds are in `sweep3.log`.

- (i) A rise time 10× slower costs 0.1 (Earth) to 1.5 (Ganymede) points.
- (ii) Requiring the deep layer to beat racemization (alanine) costs at most 0.6 points (Ganymede).
- (iii) Counting the recursive draws as wins (an upper bound, not a result) raises the moons to 15-43 % and Earth to 86 %. This often lands *above* step e, because the reversed path does not wait for the top layer's full-depth sink.
- Seeds 1-3 reproduce every share within 0.23 points of seed 0.

## Predictions: 8 of 11 held

| # | Prediction (short) | Outcome |
|---|---|---|
| P1 | Recount off reproduces step f exactly | held |
| P2 | Earth no lid 60-72 %, P_chance ≤ 2 | held (67.0, 0.0) |
| P3 | Each moon gains +1 to +6 (Enceladus +0.5 to +3); Enceladus final 0 | **failed**: Ganymede +6.21 |
| P4 | Recursion ≥ 80 % of reversed draws on moons; Earth 45-65 % | held (90-97 %; Earth 54.9 %) |
| P5 | Every world stays below step e | held |
| P6 | Rise ×10 costs ≤ 1.5 points; Earth < 0.5 | held (≤ 1.50; Earth 0.12) |
| P7 | Deep-layer racemization moves every world ≤ 0.3 points | **failed**: 8 moons lose 0.32-0.57 |
| P8 | Recursion-counted bound ≥ main; moons ≤ step e; Earth ≥ 70 % | **failed** on the moons (11 of 14 exceed step e, by 0.03-3.2); Earth 85.6 held |
| P9 | Ganymede > Europa > Enceladus (final) | held (20.2 > 10.5 > 0.0) |
| P10 | P_chance rises ≤ 1.5 points | held (≤ 1.40) |
| P11 | Seeds within 0.5 points | held (≤ 0.23) |

The three failures are small in size (P3 by 0.21 points, P7 by up to 0.27 points) or concern a bound, not the headline (P8).

## What the paper should now use

The paper should quote the sweep3 main numbers: early Earth 67 %, Ganymede 20 %, Europa 14 % without a lid and 10.5 % with one, Enceladus 7 % without a lid and 0 % with one, and the four Uranian moons 17-20 %. These replace sweep2 step f (Earth 40 %, Ganymede 14 %, Europa 9.3 / 7.1 %, Enceladus 3.9 / 0 %). Step f counted every reversed wall as "no hand", but the pinning runs show a wall above the threshold reverses and hands the ocean to the layer below. That deep hand is set by the PVED often enough to add 1-6 points on the moons and 27 on Earth. The physics choices that remain open barely move these numbers: a 10× slower rise and racemization in the deep layer each cost at most 1.5 points. The one choice that matters is what happens when the layer below is also pinned. That happens in 90-97 % of reversed moon draws. Counting all of those as P_none (the main numbers) versus as wins gives 15-43 % for the moons and 86 % for Earth, so the paper should state the main number with this as the upper range. A multi-wall run would settle it: pin_onset's radial model with 3-5 stacked walls, all past onset.

## Proposed replacement sentences (not applied)

These replace the draft's figures from before sweep2 (22 / 29 / 75 %, Ganymede 28 %, Europa 22-26 %, Enceladus 0-18 %). They also supersede the text proposed in `journals/new_analyses_sweep2.md`.

**Abstract.** Replace from "Varying all uncertain inputs together, ..." through "... (0–18% when a lid forms)." with:

> Varying all uncertain inputs together, the PVED-favoured hand is selected in 67% of plausible input sets for early Earth, 20% for Ganymede, 10–14% for Europa and 0–7% for Enceladus (the lower values when a salt-stratified lid forms). Vertical mixing is the main unknown. Weak mixing freezes the ocean into stacked layers of opposite hand. A stratified top layer stops the boundary between hands, or pushes it back so that the hand of the layer below wins. Counting all such reversals as wins gives an upper range of 15–43% for the moons and 86% for Earth. Among confirmed oceans Ganymede scores highest and Enceladus lowest.

**Section 3.1.** This section has no sweep numbers. Replace its closing sentence, "Only ocean-scale water that is mixed within the selection window can pool enough molecules (Figure 3).", with:

> Only ocean-scale water that is mixed within the selection window can pool enough molecules (Figure 3). Real oceans do not stay mixed. They break into layers, and the layer that decides the outcome holds only part of the ocean. With that included, Sections 3.2–3.5 find that the favoured hand wins in 67% of early-Earth cases and 0–20% for the confirmed icy moons.

**Conclusions.** Replace bullets 2 and 3 with:

> - Across plausible inputs, the PVED-favoured hand is selected in 67% of early-Earth cases (up to 86% if every reversed boundary hands the ocean to the layer below). Vertical mixing is the main unknown.
> - A stratified top layer stops or reverses the boundary between hands. Among confirmed oceans, Ganymede (20%) scores highest, then Europa (10–14%), and Enceladus lowest (0–7%). These shares could reach 24–33% if every reversal succeeds; a simulation with several stacked boundaries would settle this.

Table 3 (Section 3.5) and the 22 / 29 / 75 % figures in Section 3.3 would need the same update from the table above. That is outside this note's scope.
