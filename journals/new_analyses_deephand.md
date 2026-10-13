# deephand_scale: the deepest-hand rule at ocean-like scale, and a fitted cascade time (2026-10-12)

Script `deephand_scale.py` (predictions committed before any run, 8fe8ebd). Log `deephand_scale.log`.

The sections were run in the order 1, 2, 3, each with `nice -n 10` and 2 threads:
- Section 1 was split into resumable chunks.
- Section 2 took 11 min.
- Section 3 took 12 s.

Two jobs, at w = 0.005 and w = 0.0025, were killed when the machine load reached 17–18. Both were restarted from their first missing run. Every run is deterministic in (w, i), and each (w, i) appears exactly once in the log, 600 rows in all.

`multiwall.py`, `sweep2.py`, `sweep3.py` and `uncertainty.py` are imported and were not edited.

## Method in plain words

multiwall found a rule for which layer's hand fills the ocean when several stacked walls are past onset. The rule picks the hand of the layer that contains the onset depth r*, or the bottom layer if every wall is past onset. multiwall measured it with walls 0.035–0.07 of the ocean depth wide. A real ocean's walls are far thinner than that relative to its depth. Two questions were open:
- **H24:** does the rule still hold at ocean scale?
- **H23:** the wrapper times the cascade with a closed form that holds D_z constant, and the simulated time was 0.006–5.05 times that closed form. Can a fitted time do better?

1. **Scale runs (section 1).** I used multiwall's own radial solver and the wrapper's own geometry:
   - mixing grows exponentially with depth, so each wall has its own strength x = r/r*;
   - walls above r* rise and walls below it sink.

   I varied the wall width w (at the ice, as a fraction of depth) over 0.035, 0.02, 0.01, 0.005 and 0.0025. The smallest is 1/400 of the depth, 14 times thinner than anything in multiwall. The grid shrinks with w (dz = w/12).

   Each w got 150 runs (120 at 0.005 and 30 at 0.0025, for cost). The runs cycled through:
   - shell radius R/H = 2, 3, 6, 12, 100;
   - 3, 6 or 10 layers (2–9 walls);
   - ln ρ drawn as sweep2 draws it.

   Layers were laid out as sweep2 lays them: equal full layers in the √D_z-stretched depth, plus a bottom remnant u ~ U(0, 1) of a full layer. Hands alternate in 2 runs of 3 and are random in the third. In 57 runs the onset falls inside the stack, so a wall sinks.

   For each run I recorded:
   - whether the rule's hand filled the shell;
   - the end time T_end;
   - the thickness of the deciding same-hand domain in local wall widths (d/w), the quantity multiwall had flagged after its runs.
2. **Convergence (section 2).** One configuration at four widths, run at (dt, dz) and again at (dt/2, dz/2).
3. **Wrapper (section 3).**
   - **Timing fit.** On the 486 runs where the rule held, I fitted log10(T_end/T_c) as linear in ln ρ, log10(K+1) and log10(d/w). I then fed the fitted factor h per draw into `multiwall.multi` as its existing `scale` argument, which multiplies the rise and sink times.
   - **Inputs to h.** For each draw, K and d/w come from multiwall's own lines (copied in `geom`). The deciding layer is the bottom remnant f·H. Its d/w is f·H / √(2 D_loc / k), with D_loc the drawn D_z at the layer's middle. K and d/w are clipped to the simulated range.
   - **Regression.** h ≡ 1 must reproduce multiwall main.
   - **"Ocean-scale" variant.** Each draw the rule decides counts only with the probability that section 1 measured for its d/w bin. Below one wall width that probability is 0. The rest of the draw is counted as no hand.

## Results (± one standard error)

**P(rule hand fills) by width.** 0 runs unfinished at any width.

| w (fraction of depth) | runs | fills |
|---|---|---|
| 0.035 | 150 | 67.3 ± 3.8 % |
| 0.02 | 150 | 76.7 ± 3.5 % |
| 0.01 | 150 | 86.0 ± 2.8 % |
| 0.005 | 120 | 93.3 ± 2.3 % |
| 0.0025 | 30 | 96.7 ± 3.3 % |

**By thickness of the deciding domain (d/w) × width.** This table is the whole result. The rule depends on d/w, not on w.

| d/w | all | 0.035 | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|---|---|
| [0, 1) | 8.0 ± 2.9 % (88) | 16.3 ± 5.7 (43) | 0 (27) | 0 (14) | 0 (3) | 0 (1) |
| [1, 2) | 60.8 ± 5.5 % (79) | 68.3 ± 7.4 (41) | 72.0 ± 9.2 (25) | 25.0 ± 16.4 (8) | 0 (5) | – |
| [2, 3) | 97.1 ± 2.1 % (68) | 100 (20) | 96.3 ± 3.7 (27) | 92.3 ± 7.7 (13) | 100 (7) | 100 (1) |
| [3, 5) | 100 % (81) | 100 (26) | 100 (19) | 100 (22) | 100 (12) | 100 (2) |
| [5, 10) | 100 % (112) | 100 (17) | 100 (30) | 100 (39) | 100 (23) | 100 (3) |
| [10, ∞) | 100 % (172) | 100 (3) | 100 (22) | 100 (54) | 100 (70) | 100 (23) |

- **Overall by width.** The overall fill fraction rises as w shrinks only because more domains are many wall widths thick at small w.
- **Thick domains never fail.** At ≥ 3 wall widths the rule held in 365/365 runs. The thickest failing domain was 2.85 wall widths.
- **What happens on failure.** In every failure the deciding domain was a thin remnant (or a thin layer), and the hand of the next layer up filled the shell.

**By the other factors.** Each row gives all d/w, then d/w ≥ 5.

| group | all d/w | d/w ≥ 5 |
|---|---|---|
| K = 1 / 2 / 3 / 4 | 100.0 / 85.4 ± 2.6 / 94.4 ± 3.9 / 66.7 ± 8.0 % | 49/49, 115/115, 24/24, 12/12 |
| K = 5 / 6 / 7 / 8 / 9 | 77.4 ± 3.6 / 83.3 ± 9.0 / 71.4 ± 12.5 / 100 / 66.7 ± 4.7 % | 46/46, 7/7, 3/3, 4/4, 21/21 |
| n = 3 / 6 / 10 layers | 88.0 ± 2.3 / 79.5 ± 2.9 / 75.5 ± 3.0 % | 129/129, 87/87, 68/68 |
| R/H = 2 / 3 / 6 / 12 / 100 | 83.3 / 77.5 / 82.5 / 81.7 / 80.0 % (± 3.4–3.8) | all 100 % |
| sub-onset wall in the stack | 54/57 = 94.7 ± 3.0 % | 44/44 |

The lower all-d/w fractions for K and n come from thin remnants alone. Every row is 100 % at d/w ≥ 5.

**Convergence (section 2: R/H = 6, n = 6, ln ρ = 2, u = 0.6).** Halving dt and dz changes T_end by −1.72, −1.93, −1.93 and −1.93 % at w = 0.035, 0.01, 0.005 and 0.0025. The final hand is the same in every case, and the rule hand fills.

## The fitted timing model

Over the 486 finished runs where the rule held, T_end/T_c had:
- median 1.731;
- maximum 11.1;
- minimum 0.004. That minimum comes from a run whose minority hand sat in thin layers, which died at t = 16 before any wall moved.

The RMS of log10(T_end/T_c) at h = 1 is 0.494. The fit is:

    log10(T_end / T_c) = −0.570 + 0.183 ln ρ − 0.048 log10(K+1) + 0.343 log10(d/w)

It leaves a residual RMS of 0.392 (0.79 of the h = 1 scatter) and a worst residual of ×293 (the t = 16 run). The scatter is concentrated at large w:

| w | median T_end/T_c | residual RMS (log10) |
|---|---|---|
| 0.035 | 1.25 | 0.60 |
| 0.02 | 1.50 | 0.52 |
| 0.01 | 1.80 | 0.18 |
| 0.005 | 2.28 | 0.13 |
| 0.0025 | 2.21 | 0.12 |

At ocean-like widths the closed form is about 2.2 times too fast, and the fit matches individual runs to within a factor of about 1.35 (RMS log10 0.12–0.13 at w ≤ 0.005). The fit has no w term, and the median ratio stops changing at w ≤ 0.005. That plateau is the value to use at ocean scale.

On the wrapper's draws the median h is 0.67–1.19, because d/w and ln ρ partly cancel.

## World table

Shares are P_pved / P_chance / P_none in %, no lid | final, N = 200 000, seed 0.

- **Regression:** h ≡ 1 reproduces multiwall run4(0) main exactly, array for array, and all 15 worlds match multiwall.log to 0.1 ("EXACT").
- **Ocean-scale** is a floor: it counts a failed thin remnant as no hand.

| world | multiwall main | multiwall ×10 | fitted h | ocean-scale (d/w-weighted) |
|---|---|---|---|---|
| Early Earth | 99.9 / 0.0 / 0.1 \| same | 99.7 / 0.0 / 0.3 \| same | 99.9 / 0.0 / 0.0 \| same | 89.0 / 0.0 / 11.0 \| same |
| Umbriel | 53.7 / 14.9 / 31.4 \| same | 43.6 / 10.5 / 45.8 \| same | 52.9 / 17.4 / 29.6 \| same | 38.4 / 3.6 / 58.0 \| same |
| Oberon | 51.4 / 11.7 / 36.9 \| same | 40.8 / 7.9 / 51.3 \| same | 50.5 / 14.0 / 35.5 \| same | 36.7 / 2.8 / 60.5 \| same |
| Titania | 51.3 / 11.6 / 37.1 \| same | 40.9 / 8.0 / 51.1 \| same | 50.4 / 13.9 / 35.8 \| same | 36.8 / 2.8 / 60.4 \| same |
| Ariel | 49.6 / 13.5 / 36.9 \| same | 39.7 / 9.4 / 50.9 \| same | 48.9 / 16.0 / 35.1 \| same | 34.5 / 2.9 / 62.5 \| same |
| Ganymede | 44.1 / 3.4 / 52.5 \| same | 33.5 / 2.0 / 64.5 \| same | 42.7 / 4.2 / 53.1 \| same | 33.2 / 0.9 / 66.0 \| same |
| Callisto | 43.5 / 3.8 / 52.7 \| same | 33.0 / 2.3 / 64.7 \| same | 42.2 / 4.7 / 53.1 \| same | 32.4 / 1.0 / 66.6 \| same |
| Pluto | 41.1 / 6.1 / 52.8 \| same | 30.6 / 3.7 / 65.7 \| same | 40.0 / 7.5 / 52.5 \| same | 28.5 / 1.3 / 70.2 \| same |
| Triton | 37.0 / 4.7 / 58.3 \| same | 26.5 / 2.7 / 70.8 \| same | 35.8 / 5.8 / 58.4 \| same | 25.2 / 0.9 / 73.8 \| same |
| Titan | 35.8 / 2.5 / 61.7 \| same | 25.5 / 1.4 / 73.1 \| same | 34.4 / 3.3 / 62.3 \| same | 25.6 / 0.6 / 73.8 \| same |
| Europa | 39.6 / 4.9 / 55.5 \| 30.1 / 3.7 / 66.2 | 29.0 / 2.9 / 68.1 \| 22.0 / 2.2 / 75.8 | 38.4 / 6.1 / 55.5 \| 29.2 / 4.7 / 66.2 | 27.4 / 1.0 / 71.5 \| 20.8 / 0.8 / 78.4 |
| Dione | 42.0 / 10.5 / 47.5 \| 26.2 / 6.6 / 67.2 | 32.2 / 6.9 / 60.9 \| 20.1 / 4.3 / 75.6 | 41.3 / 12.8 / 46.0 \| 25.8 / 8.0 / 66.3 | 27.9 / 2.0 / 70.0 \| 17.5 / 1.3 / 81.2 |
| Miranda | 27.1 / 8.4 / 64.5 \| 9.0 / 2.7 / 88.3 | 21.3 / 5.9 / 72.8 \| 7.0 / 1.9 / 91.0 | 26.2 / 9.6 / 64.2 \| 8.7 / 3.2 / 88.2 | 20.1 / 3.5 / 76.5 \| 6.7 / 1.1 / 92.2 |
| Mimas | 17.3 / 6.1 / 76.6 \| 5.6 / 2.0 / 92.5 | 11.6 / 3.6 / 84.8 \| 3.7 / 1.1 / 95.1 | 16.9 / 7.4 / 75.7 \| 5.4 / 2.4 / 92.2 | 11.5 / 2.1 / 86.3 \| 3.7 / 0.7 / 95.6 |
| Enceladus | 28.3 / 8.8 / 63.0 \| 0.0 / 0.0 / 100.0 | 19.9 / 5.7 / 74.3 \| 0.0 / 0.0 / 100.0 | 27.1 / 10.9 / 62.0 \| 0.0 / 0.0 / 100.0 | 16.9 / 3.1 / 80.0 \| 0.0 / 0.0 / 100.0 |

**Thin deciding layers.** Of the P_pved draws the rule decides, these fractions have the deciding layer below 5 and below 1 wall widths:

| world | below 5 wall widths | below 1 wall width |
|---|---|---|
| Early Earth | 51.7 % | 29.4 % |
| icy moons | 59–74 % | 39–51 % |

The deciding layer is sweep2's bottom remnant H − K·d_w, which is drawn anywhere between 0 and one full layer. That is why about half of these layers are thin on every world.

## Predictions

| # | Prediction (short) | Outcome |
|---|---|---|
| P1 | d/w ≥ 5: ≥ 97 % pooled, ≥ 93 % each w | **held**: 284/284; every w 100 % |
| P2 | d/w < 2: ≤ 50 % | **held**: 55/167 = 32.9 ± 3.6 % |
| P3 | d/w ≥ 5: small w vs large w within 2 SE | **held**: 100 vs 100, 0.00 SE |
| P4 | w ≤ 0.005: ≥ 85 % and ≥ w = 0.035 | **held**: 141/150 = 94.0 ± 1.9 % vs 67.3 ± 3.8 % |
| P5 | K ≥ 5 and d/w ≥ 5: ≥ 95 % | **held**: 81/81 |
| P6 | sub-onset and d/w ≥ 5: ≤ 1 failure | **held**: 0/44 |
| P7 | ≤ 2 % unfinished per w | **held**: 0 |
| P8 | T_end/T_c median in [1, 3], max ≤ 30 | **held**: median 1.73, max 11.1 |
| P9 | c1 > 0 and residual RMS ≤ 0.6 × h = 1 RMS | **failed** (half): c1 = +0.183 held; ratio 0.79 |
| P10 | dt/2, dz/2: < 3 %, same hand | **held**: −1.72 to −1.93 % |
| P11 | h = 1 reproduces multiwall main exactly | **held**: EXACT, 15 worlds |
| P12 | fitted h within [×10 − 0.2, main + 0.2]; Earth ≥ 99.0 | **held**: shifts 0.0 to −1.4 points; Earth 99.9 |
| P13 | thin deciding layer ≥ 20 % on moons, ≤ 5 % on Earth | **failed** (half): moons 59–74 % held; Earth 51.7 % |
| P14 | ocean-scale Earth ≥ 95 %; every moon loses ≥ 5 points | **failed** (half): Earth 89.0; moons −5.8 to −15.3 held |

11 held and 3 failed. All three failures are halves:
- **P9:** the fit narrows the scatter less than predicted.
- **P13 and P14:** I expected Earth's deciding layers to be thick. They are thin as often as the moons', because the remnant is drawn the same way on every world.

## Do Earth's 99.9 % and the main multiwall values survive at ocean scale?

**The rule.** It survives for any layer at least three wall widths thick: 365 of 365 runs, at wall widths down to 1/400 of the depth. The width itself makes no difference once the layer is thick (P3).

**Timing.** Timing does not move anything: with the fitted cascade time, Earth stays at 99.9 % and every moon is within 1.4 points of multiwall main.

**The weak point is the deciding layer itself.** In the wrapper it is sweep2's bottom remnant, and about half the time that remnant is under five wall widths (29 % of Earth's draws are under one). In the simulations, remnants that thin usually lose, and the hand of the layer above fills the shell instead.

**Ocean-scale floor.** If every such lost draw is counted as no hand, Earth falls from 99.9 % to 89.0 % and the moons lose 6–15 points (Europa 30.1 → 20.8 with a lid, Enceladus no lid 28.3 → 16.9). That is a floor, not an estimate. The layer above is a full layer, so it is more likely than the remnant to be PVED-set itself, and many of those draws would come back.

So Earth's value is safely "at least 89 % and probably near 99.9 %", and the main multiwall values are an upper end rather than a central value until the next-layer-up case is put into the wrapper.

## Proposed draft sentences (not applied)

- **Section 2.4, replacing the last sentence of the multiwall paragraph:**
  > "At ocean-like scale (walls down to 1/400 of the depth, 600 radial runs; `deephand_scale.py`), the layer containing r* filled the shell in all 365 runs in which it was at least three wall widths thick, at every wall width (pre-registered); thinner layers usually lost to the layer above (8.0 ± 2.9 % filled below one wall width, 60.8 ± 5.5 % at one to two)."
- **Section 2.4 or 2.6, timing:**
  > "The constant-D_z closed form underestimates the cascade time by a factor of about 2.2 at ocean-like widths; a fitted correction (log10 T/T_c = −0.570 + 0.183 ln ρ − 0.048 log10(K+1) + 0.343 log10(d/w), residual scatter at most a factor of 1.5 at w ≤ 0.01) changes no share by more than 1.4 points (`deephand_scale.log`)."
- **Early Earth headline (Sections 3.3 and 6):**
  > "... rises to 99.9% when draws with several walls past the onset are given the hand of the layer containing r* (`multiwall.log`). At ocean scale this rule holds whenever that layer is at least three wall widths thick; counting the thinner deciding layers as failures gives a floor of 89% (`deephand_scale.log`)."
- **H23, close:**
  > "Fitted (`deephand_scale.py`); shares move ≤ 1.4 points. The ×10 column can be dropped or kept as a stress test."
- **H24, close or narrow:**
  > "Tested to walls 1/400 of the depth: holds at ≥ 3 wall widths (365/365). Open: the thin bottom remnant (half of all decided draws). Put the next-layer-up hand into the wrapper to replace the 89% floor."
