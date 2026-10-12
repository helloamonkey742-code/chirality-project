# pin_onset: where a wall gets pinned in a smooth depth profile of vertical mixing (2026-10-12)

Script `pin_onset.py`, log `pin_onset.log` (runtime 106 s). I committed the predictions before the run (commit 6ac662f). No existing script, log or draft was edited. `sweep2.py` was imported, not changed.

## 1 Method in plain words

- **Equation and grid.** I used `two_zone.py`'s radial equation and grid: one spherical wall, the + hand under the ice and the − hand below, R = 16 and depth H = 1.
- **New mixing profile.** Instead of two zones, vertical mixing now rises smoothly with depth: D(r) = D_top·exp(2G(R − r)). The tension gradient |d ln σ/dr| = G is then the same at every depth. The mixing rises with depth, which is the direction `sweep2.py`'s bottom/top ratio ρ ≥ 1 assumes.
- **Gradient strength.** I measure it as x = G/(2/r). Eq. 20, the sharp criterion, predicts pinning at x = 1. So κ, the measured onset, is the multiplier that `sweep2.py` calls `eq20`.

Three measurements:

1. **What happens to the wall.** I scanned x from 0.1 to 20, with 6 runs per point. Each run had a random start depth between 0.3 and 0.7 and small random noise, with fixed seeds. Each run ended in one of three ways:
   - **over:** the wall reached the floor and the top hand took the whole depth;
   - **reversed:** the wall rose to the ice and the deep hand took the whole depth;
   - **stalled:** neither happened within the time limit.
2. **κ from the wall's speed.** After the profile had relaxed, I fitted the wall's position against time. Theory predicts a normalised inward speed s = 1 − x. κ is the x where the fitted s crosses zero. I used 8 random configurations and report the mean ± SE.
   - Convergence: halving both dt and dz moved κ by less than 0.0001.
   - Also tested: radius (R = 2, 4, 16, 64), wall width (D_top = 6.25×10⁻⁴, 2.5×10⁻³, 10⁻²) and bias h (±3×10⁻⁴, ±10⁻³).
3. **Effect on `sweep2.py`.** I reran step f with `eq20` = κ and κ ± SE, using the same draws and the same seed, through `sweep2.py`'s own `chains()` and `outcome()`. A regression check reproduces `sweep2.log` exactly at multipliers ×1 and ×5.

## 2 κ table

| Case | κ ± SE (n) |
|---|---|
| **Base: R = 16, D_top = 2.5×10⁻³, h = 0** | **0.9999 ± 0.0001 (8)** |
| Base at dz/2, dt/2 | 1.0000 (shift +0.0000) |
| R = 2, D_top = 6.25×10⁻⁴ / 2.5×10⁻³ | 1.0043 ± 0.0027 (4) / 1.0124 ± 0.0070 (3) |
| R = 4 | 1.0001 ± 0.0000 / 1.0020 ± 0.0014 |
| R = 16 | 0.9998 ± 0.0000 / 1.0000 ± 0.0001 |
| R = 64 | 0.9998 ± 0.0000 / 0.9998 ± 0.0003 |
| D_top = 10⁻² (any R) | not measurable: the wall is too wide for a depth of 1 |
| h = −10⁻³ / −3×10⁻⁴ / +3×10⁻⁴ / +10⁻³ | 0.690 / 0.909 / 1.091 / 1.296 (Eq. 18 predicts 0.690 / 0.908 / 1.090 / 1.296) |

Fate scan (6 runs per point):

- every run at x ≤ 0.9 ended **over**;
- every run at x ≥ 1.1 ended **reversed**;
- **no run stalled**;
- with mixing that falls with depth (the control), x = 5 and 20 ended over.

What this means:

- **κ = 1.000.** The sharp criterion is exact for a smooth profile.
- **Radius.** κ does not depend on radius. It drifts by at most 1% at R = 2, where the wall is wide compared with the radius.
- **Bias.** Bias shifts κ exactly as Eq. 18 predicts. For the real PVED, the shift is below 3×10⁻⁶ even at 10 × g_max.
- **Why `two_zone.log` found 2–5×.** That came from its sharp joint, which is thinner (w = 0.02) than the wall itself (about 0.07). It does not show that Eq. 20 is too strict.

## 3 Per-world shares at κ

The shares are the same at κ − SE, κ and κ + SE. Across the full range of κ in the radius table (0.9998–1.0124), they move by at most 0.2 points.

| World | Step f at κ, no lid / final (%) | `sweep2.log` ×1 | `sweep2.log` ×5 | Layered draws (%) | Layered draws pinned at κ (%) | Median x of layered draws | P_pved / P_chance / P_none, no lid (%) | Same, final (%) |
|---|---|---|---|---|---|---|---|---|
| Early Earth | 39.9 | 39.9 | 40.0 | 60.1 | 100.0 | 1370 | 39.9 / 0.0 / 60.1 | same |
| Umbriel | 14.4 | 14.4 | 17.9 | 81.8 | 97.0 | 16.7 | 14.4 / 1.4 / 84.2 | same |
| Ganymede | 14.0 | 14.0 | 16.2 | 83.1 | 96.6 | 15.0 | 14.0 / 0.3 / 85.7 | same |
| Titania | 13.7 | 13.7 | 17.3 | 83.0 | 96.7 | 14.9 | 13.7 / 1.0 / 85.3 | same |
| Oberon | 13.5 | 13.5 | 17.3 | 83.3 | 96.6 | 14.5 | 13.5 / 1.0 / 85.4 | same |
| Callisto | 13.4 | 13.4 | 15.9 | 83.9 | 96.3 | 13.7 | 13.4 / 0.3 / 86.4 | same |
| Ariel | 11.7 | 11.7 | 17.0 | 86.1 | 95.5 | 11.0 | 11.7 / 1.1 / 87.2 | same |
| Pluto | 9.5 | 9.5 | 16.9 | 90.2 | 92.5 | 6.72 | 9.5 / 0.4 / 90.1 | same |
| Titan | 9.4 | 9.4 | 14.1 | 89.7 | 93.1 | 7.31 | 9.4 / 0.2 / 90.4 | same |
| Triton | 8.4 | 8.4 | 17.7 | 92.4 | 90.1 | 5.11 | 8.4 / 0.3 / 91.3 | same |
| Europa | 9.3 / 7.1 | 9.3 / 7.1 | 16.3 / 12.3 | 90.3 | 92.7 | 6.75 | 9.3 / 0.3 / 90.4 | 7.1 / 0.2 / 92.7 |
| Dione | 8.9 / 5.5 | 8.9 / 5.5 | 21.7 / 13.5 | 92.6 | 89.7 | 4.89 | 8.9 / 0.8 / 90.3 | 5.5 / 0.5 / 93.9 |
| Miranda | 12.4 / 4.1 | 12.4 / 4.1 | 30.1 / 10.0 | 98.8 | 61.9 | 1.32 | 12.4 / 2.9 / 84.7 | 4.1 / 0.9 / 95.0 |
| Mimas | 5.4 / 1.7 | 5.4 / 1.7 | 16.0 / 5.1 | 96.5 | 79.8 | 2.48 | 5.4 / 1.4 / 93.2 | 1.7 / 0.5 / 97.8 |
| Enceladus | 3.9 / 0.0 | 3.9 / 0.0 | 16.9 / 0.0 | 95.0 | 84.7 | 3.29 | 3.9 / 2.0 / 94.1 | 0.0 / 0.0 / 100.0 |

The median x of layered draws is the κ that would free half of them. On Earth κ would have to be about 1,400. On the moons it would have to be 1.3 (Miranda) to 17 (Umbriel).

## 4 Predictions: 8 of 8 held

The first run crashed in Section 3. I fixed two problems in the test code, left the predictions unchanged, and reran. Both are recorded in the `pin_onset.py` docstring:

- the speed window now stops before the wall comes near the floor or the ice;
- the time limit now uses the smallest D in the shell. The first run had used D_top, so the x = 20 control runs (mixing falling with depth) were cut off and logged as "stalled".

| # | Prediction | Result |
|---|---|---|
| P1 | κ = 1.00 ± 0.05, SE < 0.02 | held: 0.9999 ± 0.0001 |
| P2 | speed law slope −1 ± 0.1, intercept 1 ± 0.1 | held: −0.962, 0.962 |
| P3 | x ≤ 0.9 always over; x ≥ 1.1 always reversed; no stalls | held: 42/42 over, 54/54 reversed, 0 stalled |
| P4 | Control (mixing falls with depth): x = 5 and 20 over | held after the time-limit fix (the first run's x = 20 runs were cut off) |
| P5 | κ 0.9–1.1 for R ≥ 4; 0.8–1.2 at R = 2 | held in the 8 of 12 cells that could be measured; D_top = 10⁻² could not be measured |
| P6 | Bias shift within 25% of Eq. 18; the real g shifts κ by < 10⁻³ | held: ratio 0.999–1.002; at most 2.3×10⁻⁶ |
| P7 | dz/2 and dt/2 move κ by < 0.01 | held: +0.0000 |
| P8 | `sweep2.py` shares at κ within 1 point of step f; ×5 not supported; Earth 100% of layered draws pinned | held: within 0.2 points; 100.0% |

## 5 Does the step-f collapse survive?

**The threshold:**

- The threshold part of the collapse survives in full. The measured onset is κ = 1.000 ± 0.0001. Using it changes no world's share by more than 0.2 points.
- The "×5 looser" sensitivity (moons 14–30%) has no support. It came from a sharp-joint geometry that the smooth profile in `sweep2.py` does not have. The paper should drop it as the range on the headline.
- Earth is far past the threshold. Its layered draws sit about 1,400× above onset, so no plausible correction to κ would rescue them.

**What "pinned" means.** This run changes how pinning should be read:

- In a smooth profile where mixing rises with depth, a wall above onset does not stall. It runs upward, the top hand is squeezed out against the ice, and the deep hand fills the ocean. That happened in all 54 runs above onset; none stalled.
- `sweep2.py` counts these draws as "no single hand" (P_none). The simulation instead says one hand forms: the deep one.
- Whether the weak force picks that deep hand was not modelled here.

**Headline numbers:**

- **Lower end (`sweep2.py` step f, now with a measured threshold):** Earth 40%, Ganymede 14%, Europa 9.3% without a lid and 7.1% with one, Enceladus 3.9% and 0.0%.
- **Upper end, if a reversed wall counts like a sunk one:** roughly the step-e shares in `sweep2.log`, which I did not rerun: Earth 75%, Ganymede 30%, Europa 30% and 23%, Enceladus 25% and 0%.

So the paper should quote step f as its headline. It should state plainly that "pinned" means the top hand is lost, not that two hands survive, unless the base of the mixing profile is sharp. For a sharp base, `two_zone.log` showed that the wall does stall.

The one new run that would close this gap is a `sweep2.py` variant that counts reversed draws by the deep layer's molecule count, with a rise time of H·r/(2D_z(x − 1)).
