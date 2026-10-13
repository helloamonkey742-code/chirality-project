# Fronts in 2D and 3D, and whether real oceans cross the transition as a front (H21, H22; 2026-10-12)

Script `front2d.py`, log `front2d.log`. Predictions were committed before any result-producing run (commit 6d5bd4b, "front2d: predictions before run"). The RESULTS block appended to the docstring marks each one held or failed: 11 of 15 held and 4 failed. Nine sections ran in parallel under `nice -n 10` with two threads each and fixed seeds. Each ends with an assert: a 2×2-cell uniform sweep must reproduce `sim.py`'s Φ(Δ). It passed (−0.79 SE).

Nothing was edited in the existing scripts, logs, `EXPERIMENT.md` or `journals/ija_draft.md`. The per-world part imports `sweep2`, `sweep3` and `multiwall` and wraps `ocean.delta_at` in memory.

## 1 Method in plain words

**Part A: the front in 2D and 3D (H21).** `gradient_ramp.py` found that a transition sweeping along a line as a front sets the hand in a thin layer near where it starts. That layer has length ℓ_A = 2.24 (D v τ/k)^(1/3), the Airy length, and it does not grow with the size of the system. The open question was whether this survives off the line. A wide front could break into side-by-side patches, and the deciding volume might then not be "layer × whole width".

The test uses the same field equation (Eq. 14) on a strip that is 8ξ long, L = 51 cells:
- In 2D the strip is W cells wide (W = 1 cell, 2ξ, 8ξ).
- In 3D it is W × W across (W = 2ξ).
- The transition enters at one end and moves along the strip at speed v.
- The sides are periodic, so the front is a flat plane.

Front speed is given as u = v/V0: how many coherence lengths the front advances during one selection window. The settings are D = 4, giving ξ = 6.3 cells and V0 = 0.632 cells per unit time.

For each width, the bias is set so that the whole volume would give Δ_total = 2. Each point uses 250–400 runs. The effective deciding length per unit width, ℓ_eff, is read back from the measured P. It is then compared with ℓ_A ("front rule") and with L ("whole-volume rule").

**Part B: does each world cross as a front? (H22).** What makes λ cross zero in this model is the build-up of the reacting concentration (k = k₂c). Assumptions:
- **F1:** the supply is local (vents) and is spread by horizontal mixing D_h. The crossing time then differs across the ocean by about the mixing time t_mix = (L/A)²/D_h, A = 5.5, which is sweep2's healing-across time. The same spread arises in both regimes:
  - when mixing is slow, as a travelling front;
  - when mixing is fast, as leftover unevenness of order t_mix/τ.
  The front speed is therefore v = L/t_mix.
- **F2:** the hand field mixes with the same D_h. Then ξ = √(D_h/k)(kτ)^(1/4), V0 = √(D_h k)(kτ)^(−1/4), and **u = A² ξ/L**.
- **F3:** a front spans the depth of the deciding layer and a length L across the ocean. The deciding volume shrinks from f·V to f·V·πℓ_A/(4L), which equals f·V·0.058 u^(4/3).
- **F4:** draws with u < u* = 3 get the reduced N; the others keep multiwall's whole-volume N. Sensitivities: u* = 1, and v × A or v / A.
- **F5:** arrival, healing, wall and recursion rules are unchanged.
- **F6:** a globally set drive, such as cooling through the shell, is uniform, so the shares equal multiwall main.

With the front rule switched off, the code must reproduce multiwall Section 4 exactly. It does: equal floats against `multiwall.run4(0)` and the same printed 0.1 values as `multiwall.log` for all 15 worlds.

## 2 Part A: tables

P ± SE is the fraction of runs that end with the favoured hand as the majority. Φ(Δ_tot) is the whole-volume prediction and Φ(Δ_A) is the front-layer prediction. ℓ_eff/ℓ_A should be about 1 if the front rule holds. ℓ_eff/L should be about 1 if the whole-volume rule holds.

**L = 8ξ; Φ(Δ_tot) = 0.977 in every row**

| u | Φ(Δ_A) | W = 1 cell: P / ℓ_eff/ℓ_A / one hand | 2D, W = 2ξ | 2D, W = 8ξ | 3D, 2ξ × 2ξ |
|---|---|---|---|---|---|
| 0.1 | 0.764 | 0.767±0.021 / 1.04±0.20 / 1.000 | 0.757±0.021 / 0.95±0.19 / 1.000 | 0.780±0.021 / 1.16±0.21 / 0.930 | – |
| 0.3 | 0.806 | 0.780±0.021 / 0.80±0.15 / 1.000 | 0.765±0.021 / 0.70±0.13 / 1.000 | 0.838±0.018 / 1.30±0.20 / 0.955 | 0.788±0.026 / 0.86±0.19 / 1.000 |
| 1 | 0.854 | 0.835±0.019 / 0.85±0.13 / 0.998 | 0.868±0.017 / 1.12±0.16 / 1.000 | 0.885±0.016 / 1.30±0.18 / 0.970 | 0.872±0.021 / 1.16±0.21 / 1.000 |
| 3 | 0.897 | 0.920±0.014 / 1.23±0.16 / 0.930 | 0.917±0.014 / 1.20±0.16 / 0.927 | 0.910±0.014 / 1.12±0.15 / 0.953 | 0.892±0.020 / 0.96±0.16 / 0.940 |
| 10 | 0.939 | 0.975±0.008 / 1.61±0.22 / 0.892 | 0.958±0.010 / 1.24±0.16 / 0.905 | 0.965±0.009 / 1.37±0.18 / 0.945 | – |
| ∞ | 0.977 | 0.980±0.007 / ℓ_eff/L 1.05 / 0.890 | 0.983±0.007 / ℓ_eff/L 1.11 / 0.902 | 0.970±0.009 / ℓ_eff/L 0.88 / 0.917 | 0.980±0.009 / ℓ_eff/L 1.05 / 0.912 |

**Size scan (2D, W = 2ξ, same bias at every size)**

| u | L | Δ_tot | P | ℓ_eff (cells) | ℓ_eff/ℓ_A |
|---|---|---|---|---|---|
| 0.3 | 4ξ | 1.41 | 0.797±0.020 | 8.8±1.5 | 0.93 |
| 0.3 | 8ξ | 2.00 | 0.765±0.021 | 6.7±1.3 | 0.70 |
| 0.3 | 16ξ | 2.83 | 0.825±0.019 | 11.1±1.8 | 1.17 |
| 3 | 4ξ | 1.41 | 0.905±0.015 | 21.9±2.9 | 1.07 |
| 3 | 8ξ | 2.00 | 0.917±0.014 | 24.6±3.2 | 1.20 |
| 3 | 16ξ | 2.83 | 0.975±0.008 | 49.0±6.7 | 2.40 |

**Convergence (2D, W = 2ξ, L = 8ξ, u = 0.3)**
- Base: P = 0.765 ± 0.021.
- dt/2: 0.807 ± 0.020 (+1.45 SE).
- dx/2: 0.825 ± 0.019 (+2.12 SE).
- ℓ_eff/ℓ_A goes 0.70 → 1.01 → 1.17.

**u\* in 2D vs 1D.**
- At u = 0.3 the deciding length stays at the Airy length as the system grows fourfold: 8.8 → 11.1 cells, against ℓ_A = 9.5.
- At u = 3 it grows with the system: 21.9 → 49.0 cells, +4.1 SE in P.
- So in 2D u* lies in (0.3, 3]. u = 1 was not size-scanned, so the 1D bracket (1, 3] cannot be narrowed further. The two are consistent.
- At a fixed L = 8ξ the ratio ℓ_eff/ℓ_A cannot show the crossover: ℓ_A(u = 3) is already 0.4 L. This is why P5 and P9 failed as written.
- 3D (2ξ × 2ξ) follows 2D: 0.86, 1.16, 0.96 at u = 0.3, 1, 3, and the whole volume at u = ∞ (1.05).

## 3 Part B: per-world front fraction and shares (seed 0, N = 200 000 draws per world)

| world | % draws u < 3 (all / reaching one hand) | median u | % u < 1 | % t_mix > τ | P_pved/P_chance/P_none no lid, front rule off = on | final, off = on | ΔP_pved u*=3 / u*=1 / v×A / v÷A (points) |
|---|---|---|---|---|---|---|---|
| Early Earth ocean | 2.8 / 2.8 | 38 | 0.1 | 0.0 | 99.9 / 0.0 / 0.1 | same | 0.00 / 0.00 / 0.00 / 0.00 |
| Umbriel | 4.9 / 5.3 | 108 | 1.1 | 0.0 | 53.7 / 14.9 / 31.4 | same | −0.00 / 0.00 / 0.00 / −0.03 |
| Oberon | 6.2 / 6.9 | 84 | 1.6 | 0.0 | 51.4 / 11.7 / 36.9 | same | −0.00 / 0.00 / 0.00 / −0.03 |
| Titania | 6.5 / 7.1 | 81 | 1.7 | 0.0 | 51.3 / 11.6 / 37.1 | same | −0.00 / −0.00 / 0.00 / −0.03 |
| Ariel | 4.7 / 5.1 | 110 | 1.1 | 0.0 | 49.6 / 13.5 / 36.9 | same | 0.00 / 0.00 / 0.00 / −0.03 |
| Ganymede | 17.2 / 19.2 | 25 | 7.2 | 0.0 | 44.1 / 3.4 / 52.5 | same | −0.01 / 0.00 / 0.00 / −0.08 |
| Callisto | 16.2 / 18.1 | 27 | 6.6 | 0.0 | 43.5 / 3.8 / 52.7 | same | −0.01 / −0.00 / 0.00 / −0.06 |
| Pluto | 9.5 / 10.4 | 55 | 3.0 | 0.0 | 41.1 / 6.1 / 52.8 | same | −0.00 / 0.00 / 0.00 / −0.03 |
| Triton | 10.6 / 11.4 | 47 | 3.6 | 0.0 | 37.0 / 4.7 / 58.3 | same | −0.00 / 0.00 / 0.00 / −0.03 |
| Titan | 16.9 / 18.6 | 25 | 6.9 | 0.0 | 35.8 / 2.5 / 61.7 (on: P_chance 2.6) | same | −0.00 / 0.00 / 0.00 / −0.04 |
| Europa | 11.2 / 12.0 | 44 | 3.9 | 0.0 | 39.6 / 4.9 / 55.5 | 30.1 / 3.7 / 66.2 | −0.00 / 0.00 / 0.00 / −0.03 |
| Dione | 4.5 / 5.0 | 114 | 1.0 | 0.0 | 42.0 / 10.5 / 47.5 | 26.2 / 6.6 / 67.2 | 0.00 / 0.00 / 0.00 / −0.01 |
| Miranda | 2.0 / 2.5 | 210 | 0.2 | 0.0 | 27.1 / 8.4 / 64.5 | 9.0 / 2.7 / 88.3 | 0.00 / 0.00 / 0.00 / −0.00 |
| Mimas | 4.9 / 6.5 | 106 | 1.1 | 0.0 | 17.3 / 6.1 / 76.6 | 5.6 / 2.0 / 92.5 | −0.00 / 0.00 / 0.00 / −0.04 |
| Enceladus | 4.2 / 5.0 | 131 | 1.0 | 0.0 | 28.3 / 8.8 / 63.0 | 0.0 / 0.0 / 100.0 | −0.00 / 0.00 / 0.00 / −0.04 |

Seed 1 with the front rule on agrees with seed 0 within 0.36 points.

Why nothing moves (post hoc, scratch diagnostic, not pre-registered): u = A²ξ/L is small only when ξ is small, which needs a fast reaction (large k). The front-flagged draws have a median k of 2·10⁻⁸ to 1.6·10⁻⁵ s⁻¹, against 3·10⁻¹¹ to 6·10⁻⁹ s⁻¹ for all draws that reach one hand. The same large k and c put their F·Δ at 10^3.7 to 10^5.6 times the threshold of 2. The front rule cuts Δ by a median factor of about 0.3, and 99.9–100 % of these draws still win.

## 4 Which predictions held

| # | prediction | outcome |
|---|---|---|
| P1 | uniform 2×2 limit = Φ(Δ) | held (−0.79 SE) |
| P2 | 1D limit at D = 4, ℓ_eff/ℓ_A in [0.67, 1.5] for u ≤ 1 | held (1.04, 0.80, 0.85) |
| P3 | 2D slow front, ℓ_eff per width / ℓ_A in [0.5, 1.5] | held (0.95, 0.70, 1.16, 1.30) |
| P4 | 2D size independence at u = 0.3; far below whole volume at 16ξ | held (max 2.12 SE; 9.1 SE below) |
| P5 | u*_2D in (1, 3]: ratio > 1.5 at u = 3 (L = 8ξ); P grows with L at u = 3 | **failed**: ratio 1.20 / 1.12 (badly sized test); size part held (+4.1 SE) |
| P6 | 2D u = ∞, ℓ_eff/L in [0.64, 1.3] | held (1.11, 0.88) |
| P7 | one-hand fraction drops ≥ 0.1 at W = 8ξ, slow fronts | **failed** (drops 0.07, 0.045) |
| P8 | dt/2 and dx/2 within 2 SE | **failed**, narrowly: dx/2 +2.12 SE |
| P9 | 3D: front rule at u ≤ 1, ratio > 1.5 at u = 3, whole volume at ∞ | **failed** (u = 3 part only: 0.96) |
| P10 | front rule off = multiwall main exactly | held (EXACT) |
| P11 | front fraction ≤ 30 % (Earth ≤ 15 %), Enceladus < Europa < Ganymede | held (2.0–17.2 %; 4.2 < 11.2 < 17.2) |
| P12 | front fraction among one-hand draws ≥ among all draws | held (Earth equal) |
| P13 | shares fall 0–5 points, P_none fixed, Earth ≥ 99 %, ranking kept | held (≤ 0.01 points) |
| P14 | v / A loses ≤ 10 points | held (≤ 0.08) |
| P15 | seed 1 within 0.5 points | held (0.36) |

## 5 Do fronts change the headlines?

No. The front result from the line holds for a flat front in 2D and 3D:
- When the transition sweeps slower than about one to three coherence lengths per selection window, only a thin layer near where it starts decides the hand.
- That layer has the Airy thickness and spans the whole width of the front, and it does not grow with the ocean.
- Wide fronts sometimes split into side-by-side patches (5–7 % of runs at 8ξ), but the majority still follows the layer's total bias.

In the ocean worlds, under the stated assumptions (vent-supplied drive, crossing spread equal to the mixing time), between 2 % (Miranda) and 17 % (Ganymede) of draws cross as a front. The task's simpler test, "is the drive spread across the ocean faster than it changes?", says yes in 100 % of draws for every world: t_mix < τ always.

The front draws are the fast-chemistry draws, and those win by factors of thousands. Shrinking their deciding volume by about tenfold flips almost none of them. No world's share moves by more than 0.01 points (0.08 points in the most front-friendly sensitivity). Early Earth stays at 99.9 %, and Ganymede 44.1 > Europa 30.1 > Enceladus 0.0 at the final share.

The fronts objection is answered for this model. The answer leans on F1–F3. In particular it assumes the crossing spread is set by horizontal mixing and not by something slower, such as a single long-lived vent field in a stagnant ocean, which this model does not include.

### Proposed draft sentences (NOT applied; `journals/ija_draft.md` is read-only here)

- Section 2.3, after the "transition that sweeps as a front" paragraph: "A planar front in 2D and 3D behaves as on a line. Below u* ≈ 1–3 the deciding volume is the Airy layer ℓ_A = 2.24 (D v τ/k)^(1/3) times the full front width, independent of system size (ℓ_eff/ℓ_A = 0.7–1.3 for u ≤ 1; `front2d.py`)."
- Section 2.6 or 4: "Taking the crossing-time spread across an ocean to be the horizontal mixing time (L/A)²/D_h gives u = A²ξ/L. Only 2–17 % of draws per world have u < 3, and they are the fast-reaction draws, whose Δ exceeds the threshold by 10³–10⁵. Replacing their N by the front layer changes no world's P_pved by more than 0.01 points (0.08 with a 5.5× slower front)."
- Notes H21/H22: mark them as tested by `front2d.py`. H21: the 2D/3D front gives N = ℓ_A × width. H22: the criterion is u = A²ξ/L per world, and the front rule leaves every share unchanged.

## 6 References (DOIs checked against api.crossref.org, 2026-10-12)

- Kibble, T. W. B. & Volovik, G. E. (1997). On phase ordering behind the propagating front of a second-order transition. *JETP Letters* 65. doi:10.1134/1.567332
- Dziarmaga, J., Laguna, P. & Zurek, W. H. (1999). Symmetry breaking with a slant: topological defects after an inhomogeneous quench. *Phys. Rev. Lett.* 82. doi:10.1103/PhysRevLett.82.4749
- Zurek, W. H. (1985). Cosmological experiments in superfluid helium? *Nature* 317. doi:10.1038/317505a0
- van Saarloos, W. (2003). Front propagation into unstable states. *Physics Reports* 386. doi:10.1016/j.physrep.2003.08.001

## 7 Caveats

- **Grid.** D = 4 gives ξ = 6.3 cells, and the Airy width at u = 0.1 is under 3 cells. Halving dx raised P by 2.1 SE (ℓ_eff/ℓ_A 0.70 → 1.17), so ℓ_eff at D = 4 may be up to ~30 % low. Every value stays inside the predicted band.
- **u\* in 2D.** It is bracketed only as (0.3, 3]. Narrowing it needs a size scan at u = 1.
- **3D.** Only one cross-section (2ξ × 2ξ) and 250 runs per point were tested.
- **Part B assumptions.** F1–F3 are order-one physical assumptions. In particular, the front length L and the use of D_h for both the drive and the hand field are assumptions. The sensitivities (u* = 1, v × A, v / A) all leave the shares unchanged.
- **Drive scope.** A drive slower than mixing, or one confined to a stagnant region, is outside this model.
- **Run times.** Four sections ran 20–25 min, past the ~20 min budget. The cause was machine load from other jobs (load average ~38), not the code.
