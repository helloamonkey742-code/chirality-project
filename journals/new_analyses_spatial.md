# Spatial tests of "the whole connected volume counts" (f = 1): patches and fronts (2026-10-12)

Scripts `kz_patches.py` and `gradient_ramp.py`; logs `kz_patches.log` and `gradient_ramp.log`. Predictions were committed before any result-producing run (commits 139df02 and 62bc148). The RESULTS blocks appended to each docstring state which predictions held and which failed. Every section ran under `nice -n 10` with a fixed seed. Every section ends with an assert: a 3-cell uniform sweep must reproduce `sim.py`'s Φ(Δ). It passed everywhere (kz +1.13 SE, ramp +0.33 SE). Nothing in the paper drafts was edited. The text in Section 6 is a proposal only.

Two referee objections are tested:

1. **Patches (Kibble–Zurek).** A large ocean does not choose its hand as one unit. It breaks into many patches that each choose separately. Does the Δ of the whole volume still predict the outcome?
2. **Fronts.** If the transition reaches one region first (near a vent or under thin ice) and that region then invades the rest, maybe only the first region counts.

## 1 Method in plain words

The model is the paper's field equation (Eq. 14): the local rule of `sim.py` plus diffusion. It runs on a line, a square or a ramp of cells. Each cell has its own counting noise ε = 10⁻⁴. A system of n cells therefore has total noise ε/n, and the whole-volume prediction is P = Φ(Δ_total) with Δ_total = `sim.delta(g, ε/n, γ)`, γ = 0.01. The coherence length is ξ = √D·γ^(−¼). In 1D the patches at freeze-out are 12.8ξ long (`domains.py`).

- **Patches, 1D (A1).** Line lengths are 0.5ξ to 512ξ (up to 40 patches). The bias g is set so that Δ_total = 0.5, 1, 1.5 and 2.5, with 600 runs per point. The fit is P_sim = Φ(c·Δ_total). c = 1 means the whole volume counts exactly.
- **Patches, 2D (A2).** Squares are 8² to 64² cells (up to 20ξ on a side), with Δ_total = 1 and 2 and 400 runs per point. The majority hand is recorded just after selection (λ = 0.4), at the end of the sweep, and after t = 500 of coarsening at λ = 1.
- **Imbalance (B).** A 64² random patchwork with zero bias and no noise is started with a set excess m of one hand, then left to coarsen. m is measured in units of the patchwork's own natural spread σ_m (exact, from the lattice correlation). The test asks how often the initially larger hand wins. The patch widths are w = 2 and 4 (M_eff = 1/(4σ_m²) = 117 and 29 patches).
- **Fronts (ramp).** A line with closed ends sees λ(x, t) = γ(t − x/v). Every cell gets the same sweep as `sim.py`, but the start time moves along the line at speed v. Here u = v/V₀ with V₀ = √D·γ^¼. u is the number of coherence lengths the front advances in one selection window (t_w = γ^(−½)). One bias is used throughout, giving Δ_total = 2, 4 and 8 at L = 4, 16 and 64ξ. Each simulated P is converted to the length ℓ_eff that `sim.py` would need to give the same P. That length is compared with:
  - the whole line, L;
  - the task's first-crossing length ℓ_first = vγ^(−½);
  - the pre-registered Airy length ℓ_A = 2.24(Dv/γ)^(1/3) = 2.24·u^(1/3)·ξ. This is the width of the lowest unstable mode next to the wall where the front starts.
- **Convergence.** Halving dt and halving dx (with dt/4) moved P by at most 1.70 combined SE (kz) and 1.68 combined SE (ramp).

## 2 Patches: P_sim against the whole-volume prediction

**1D (A1), end of sweep, n = 600 per point (±SE):**

| L | patches M | Δ=0.5 | Δ=1 | Δ=1.5 | Δ=2.5 | fitted c |
|---|---|---|---|---|---|---|
| Φ(Δ) (c = 1) | | 0.691 | 0.841 | 0.933 | 0.994 | 1 |
| Φ(0.914Δ) (Gaussian-field theory) | | 0.676 | 0.820 | 0.915 | 0.989 | 0.914 |
| Φ(0.798Δ) (box patch sum) | | 0.655 | 0.788 | 0.884 | 0.977 | 0.798 |
| 0.5ξ | 0.04 | 0.698±0.019 | 0.828±0.015 | 0.937±0.010 | 0.992±0.004 | 0.986±0.033 |
| 2ξ | 0.16 | 0.690±0.019 | 0.847±0.015 | 0.917±0.011 | 0.985±0.005 | 0.936±0.031 |
| 8ξ | 0.62 | 0.708±0.019 | 0.835±0.015 | 0.937±0.010 | 0.995±0.003 | 1.015±0.034 |
| 32ξ | 2.5 | 0.691±0.019 | 0.817±0.016 | 0.926±0.011 | 0.990±0.004 | 0.941±0.031 |
| 128ξ | 10 | 0.662±0.019 | 0.823±0.016 | 0.907±0.012 | 0.983±0.005 | 0.876±0.029 |
| 512ξ | 40 | 0.667±0.019 | 0.826±0.015 | 0.909±0.012 | 0.982±0.005 | 0.882±0.029 |

**Fitted c (plateau, 128ξ and 512ξ pooled): 0.879 ± 0.021.**
- The theory written before the runs predicted 0.914 (−1.7 SE).
- The box patch-sum value √(2/π) = 0.798 is excluded by 2.9 SE at 512ξ alone and by 4.0 SE pooled.
- c is flat from 10 to 40 patches, so P still grows as √N; only the prefactor changes.
- The majority set at selection was still the majority at the end in ≥ 99.0% of runs at every L.

**2D (A2), n = 400 per point:**

| side | stage | Δ=1 | Δ=2 | c | single hand (Δ=1) |
|---|---|---|---|---|---|
| Φ(Δ) / Φ(0.958Δ) | | 0.841 / 0.831 | 0.977 / 0.972 | 1 / 0.958 | |
| 8 (2.5ξ) | final | 0.838±0.018 | 0.978±0.007 | 0.994±0.051 | 1.00 |
| 16 (5.1ξ) | final | 0.835±0.019 | 0.968±0.009 | 0.944±0.048 | 1.00 |
| 32 (10ξ) | final | 0.875±0.017 | 0.979±0.007 | 1.077±0.055 | 0.89 |
| 64 (20ξ) | early (λ=0.4) | 0.850±0.018 | 0.980±0.007 | 1.031±0.053 | 0.01 |
| 64 (20ξ) | end of sweep | 0.843±0.018 | 0.973±0.008 | 0.979±0.050 | 0.15 |
| 64 (20ξ) | after coarsening | 0.844±0.018 | 0.973±0.008 | 0.981±0.050 | 0.81 |

In 2D, c = 0.98 ± 0.05 after coarsening (theory 0.958). Coarsening from 15% to 81% single-hand runs did not change P at all. The 2D runs are 2.5× less precise than the 1D runs and cannot separate 0.96 from 1.

## 3 Imbalance curve (B): does coarsening keep the majority?

P(the initially larger hand wins) at t = 600, n = 300 per point. The values at t = 150 and 300 are the same within 0.01.

| x = m/σ_m | m·√M_eff | P, w = 2 (M_eff 117) | P, w = 4 (M_eff 29) |
|---|---|---|---|
| 0 | 0 | 0.483±0.029 | 0.533±0.029 |
| 0.5 | 0.25 | 0.790±0.024 | 0.915±0.016 |
| 1 | 0.5 | 0.943±0.013 | 0.967±0.010 |
| 1.5 | 0.75 | 0.985±0.007 | 0.997±0.003 |
| 2 | 1.0 | 0.990±0.006 | 1.000±0.002 |
| 3 | 1.5 | 1.000±0.002 | 1.000±0.002 |
| probit κ in P = Φ(κx) | | 1.45±0.07 | 2.12±0.12 |

- Coarsening keeps the majority more often than not.
- An excess of one natural spread (x = 1) already wins 94–97% of the time.
- The curves do not collapse onto one in units of 1/√M (4.3 SE apart at x = 0.5). Larger patches keep their majority better.
- All the loss happens in the first relaxation of the sharp, thresholded starting pattern. Nothing changes after t = 150.
- In the physical sweep (A2) the field grows smoothly from noise, and the same loss did not appear: P_end = P_final. So B is a pessimistic bound. The factor it implies on Δ, κ/√(1+κ²), is 0.82 (w = 2) to 0.90 (w = 4).

## 4 Fronts: crossover between "front-set" and "whole volume"

P (±SE), with the three predictions, for Δ_total = 2 (4ξ), 4 (16ξ) and 8 (64ξ). n = 600 at 4ξ and 16ξ; n = 400 at 64ξ.

| u | L | P_sim | Φ(Δ_total) | Φ(Δ_first) | Φ(Δ_Airy) | ℓ_eff (cells) | ℓ_A | ℓ_first | single hand |
|---|---|---|---|---|---|---|---|---|---|
| 0.03 | 4ξ | 0.837±0.015 | 0.978 | 0.569 | 0.798 | 12.1±1.5 | 8.8 | 0.4 | 1.000 |
| 0.03 | 16ξ | 0.785±0.017 | 1.000 | 0.569 | 0.798 | 7.9±1.1 | 8.8 | 0.4 | 1.000 |
| 0.1 | 4ξ | 0.833±0.015 | 0.978 | 0.624 | 0.846 | 11.8±1.5 | 13.1 | 1.3 | 1.000 |
| 0.1 | 16ξ | 0.838±0.015 | 1.000 | 0.624 | 0.846 | 12.3±1.5 | 13.1 | 1.3 | 1.000 |
| 0.3 | 4ξ | 0.897±0.012 | 0.978 | 0.708 | 0.890 | 20.1±2.2 | 19.0 | 3.8 | 1.000 |
| 0.3 | 16ξ | 0.887±0.013 | 1.000 | 0.708 | 0.890 | 18.5±2.1 | 19.0 | 3.8 | 1.000 |
| 0.3 | 64ξ | 0.907±0.014 | 1.000 | 0.708 | 0.890 | 22.2±2.9 | 19.0 | 3.8 | 1.000 |
| 1 | 4ξ | 0.928±0.011 | 0.978 | 0.842 | 0.933 | 27.0±2.8 | 28.3 | 12.6 | 1.000 |
| 1 | 16ξ | 0.940±0.010 | 1.000 | 0.842 | 0.933 | 30.5±3.2 | 28.3 | 12.6 | 1.000 |
| 1 | 64ξ | 0.920±0.014 | 1.000 | 0.842 | 0.933 | 24.9±3.2 | 28.3 | 12.6 | 1.000 |
| 3 | 4ξ | 0.958±0.008 | 0.978 | 0.959 | 0.964 | 37.9±4.0 | 40.9 | 37.9 | 0.990 |
| 3 | 16ξ | 0.997±0.002 | 1.000 | 0.959 | 0.964 | 92.9±16.0 | 40.9 | 37.9 | 0.958 |
| 3 | 64ξ | 1.000±0.001 | 1.000 | 0.959 | 0.964 | > 100 | 40.9 | 37.9 | 0.965 |
| 10 | 16ξ / 64ξ | 1.000 / 1.000 | 1.000 | 0.999 | 0.986 | saturated | 61 | 127 | 0.980 / 0.958 |
| 30 | 4ξ | 0.985±0.005 | 0.978 | 0.978 | 0.978 | 59.5±7.2 | 51 (=L) | 51 | 0.983 |
| ∞ | 4ξ | 0.975±0.006 | 0.978 | 0.978 | 0.978 | 48.5±5.4 | 51 (=L) | 51 | 0.990 |

**For a slow front (u ≤ 1):**
- P does not depend on the size of the line (4ξ to 64ξ).
- P falls far below the whole-volume value (8.7–12.6 SE at 16ξ).
- The Airy length predicts the effective length within 0.90–1.37× from u = 0.03 to 1.
- The first-crossing length vγ^(−½) is 9–30× too small at u ≤ 0.1, and 2× too small at u = 1.
- Every slow-front run ends as a single hand: the region where the front starts sets the hand and invades the rest.

**Crossover:**
- Between u = 1 and u = 3, ℓ_eff leaves ℓ_A and rises towards L (16ξ: 30.5 → 92.9), and the first independent second-hand domains appear (single hand 1.000 → 0.958–0.965, 3.8–5.1 SE).
- The crossover is therefore at **u\* ≈ 1–3**, just below the pre-registered pulled-front estimate u\* ≈ 3.5. It was read post hoc, because the pre-registered diagnostic (single-hand fraction crossing 50% at 64ξ) saturated: with Δ_total = 8, even independent patches all chose the favoured hand.
- For u ≥ 10 the line behaves like a uniform sweep (u = 30 matches v = ∞, 1.3 SE).

## 5 Which predictions held

| Script | Held | Failed | Ambiguous |
|---|---|---|---|
| kz_patches | P1 uniform limit; P2 one c at 10–40 patches; P5 walls frozen in 1D; P10 convergence | P3 (0.80 excluded at 2.9 SE, not > 3); P4 (c(2ξ) = 0.936 < 0.97); P6 (c_early = 1.03, outside 0.96 ± 0.05); P7 (coarsening did not lower c into 0.6–0.9); P8 (no collapse between w = 2 and 4) | P9 (holds with the w = 4 κ, fails with the w = 2 κ; its mechanism is contradicted by P7) |
| gradient_ramp | P1; P2 fast front = uniform; P4 slow-front P independent of L, far below whole-volume; P5 Airy length within [0.67, 1.5]; P7 convergence | P3 (inheritance measured at λ(0) = 0.4 is too early for u ≤ 0.1, by the script's own Airy shift); P6 (diagnostic saturated; post hoc u\* ∈ (1, 3]) | — |

The failures cut both ways:
- P6 and P7 of kz_patches failed because the hand was kept *better* than predicted. Coarsening did not cost anything.
- P3 and P6 of gradient_ramp failed because of how the prediction was measured, not because of the physics claim it tested.

## 6 Does f = 1 survive?

**Plain language.**
- **Patches.** Breaking an ocean into many patches does not break the scale rule. The chance of the favoured hand still grows with the square root of the number of molecules in the whole connected volume. It only needs a small discount, because a patch's vote is all-or-nothing while the bias is spread smoothly. In the simulations, the effective Δ was 0.88× the whole-volume Δ on a line and 0.98× (±0.05) on a plane. Theory gives 0.98× in 3D. In molecule terms, the whole volume counts at 77% (1D) to about 96% (2D/3D). The pessimistic "independent box patches" value (64%) is ruled out. Later coarsening kept the majority hand: on the plane it changed nothing.
- **Fronts.** f = 1 *fails* when the transition does not arrive everywhere at once but sweeps in as a slow front. If the front advances less than about one to three coherence lengths per selection window (u ≲ 1–3), the hand is set by a thin layer next to where the front starts. That layer is about 2.2·u^(1/3) coherence lengths thick, and the rest of the ocean simply inherits its hand. Then N is the molecules in that layer, not in the ocean, and it does not grow with the ocean's size. The task's guess for this layer (front speed × window) is far too thin. For fronts faster than u ≈ 3, the whole volume counts again, with the patch discount above.

**Correction factors (to replace f = 1):**

$$f = c_d^2 \times \begin{cases} 1 & u \gtrsim 3 \\ \min\!\left(1,\ 2.24\,u^{1/3}\,\xi / L_\parallel\right) & u \lesssim 1 \end{cases}$$

- c_d² is 0.77 ± 0.04 in 1D (simulated), 0.96 ± 0.10 in 2D (simulated) and 0.96 in 3D (theory, not simulated).
- L_∥ is the extent of the mixed volume along the front's direction.
- Physical units:
  - ξ = √(D/k)·(kτ)^(1/4) (the same scale as `domains.py` check A)
  - selection window t_w = √(τ/k)
  - V₀ = ξ/t_w = √(Dk)·(kτ)^(−1/4)
  - u = v/V₀
  - ℓ_A = 2.24·(Dvτ/k)^(1/3)
- Here D is the diffusivity along the front's direction (after the stretching of Eq. 15, if the front is vertical).
- 1 < u < 3 is a transition region with no simple formula.
- The front result is 1D. In 3D a planar front of cross-section A would make N ≈ c_{2}²·A·ℓ_A·(molecules per volume). That extrapolation is untested.

**What the paper needs to decide.** The paper's sweep is uniform in space (λ = t/τ everywhere). f = 1 is therefore safe, up to the c² factor, exactly when the transition's arrival time varies across the mixed volume by less than about (L/ξ)/3 selection windows. Whether icy-moon drives (vent output, ice-shell thickness, feedstock delivery) are that uniform is not known, and it should be stated as an assumption.

### Proposed replacement text (NOT applied; `journals/ija_draft.md` is being edited by another agent)

**Section 2.3, first bullet (replace "The whole volume counts: f = 1 in Eq. 12. This is an assumption; …").**

> - H' ≤ L (the stretched ocean is wider than it is deep). After the transition the ocean forms a patchwork of regions of each hand that coarsens to a single domain (Section 2.4). The whole volume then counts, with a small discount: each patch votes with its sign, not its size, so the effective Δ is c·Δ with c = [∫ρ / ∫arcsin ρ]^½ for a field correlation ρ. This gives 0.914 in 1D, 0.958 in 2D and 0.979 in 3D. Simulations of Eq. 14 gave c = 0.88 ± 0.02 on a line of 10–40 patches and 0.98 ± 0.05 on a plane, with the √N scaling intact (`kz_patches.py`). I therefore use f = c² ≈ 0.96 in Eq. 12; using the 1D value (0.77) changes Δ by 10%. Section 3.2 reports a shallow-slab test that does not fully support coarsening to a single domain.
>
> This requires the transition to arrive across the mixed volume nearly at once. If it sweeps in as a front of speed v, the hand is set by a layer next to where the front starts, of thickness ℓ_A = 2.24 (D v τ / k)^(1/3), and the rest of the ocean inherits it. This happens when u = v/V₀ ≲ 1–3, with V₀ = √(Dk)(kτ)^(−1/4), i.e. when the front advances less than one to three coherence lengths ξ = √(D/k)(kτ)^(1/4) per selection window √(τ/k) (`gradient_ramp.py`). Then f is the share of the molecules in that layer, which does not grow with ocean size. I assume a uniform transition (u ≳ 3), i.e. arrival times that differ across the ocean by less than about (L/ξ)/3 windows.

**Section 2.4, "Why walls move", after "…coarsens to a single domain in a time (L/A)²/D_h." (insert).**

> Coarsening keeps the majority set at selection. In 1D the majority at the end of the sweep matched the majority just after selection in ≥ 99% of runs. In 2D, holding the field at λ = 1 until 81% of runs were a single hand left the share of the favoured hand unchanged (0.844 vs 0.843 at Δ = 1; `kz_patches.py`). Even from an artificial sharp patchwork with no bias, an initial excess of one natural spread won 94–97% of the time (`kz_patches.py`, part B). The selection set during the window is therefore the outcome; coarsening does not reset it.

**Section 4, "Mixing is the key unknown" bullet (append).**

> The whole-volume count (f = 1 for a shallow ocean) holds for a transition that arrives everywhere at once, up to a factor c² ≈ 0.8–0.96 for selection in patches (Section 2.3). It fails for a transition that sweeps in as a slow front (less than about one coherence length per selection window). Then only a layer about 2–3 coherence lengths thick next to where the front starts sets the hand (Kibble and Volovik, 1997; Dziarmaga et al., 1999). Whether the drive in icy-moon oceans is spatially uniform enough is not known. The front test is one-dimensional.

## 7 References (DOIs checked against api.crossref.org, 2026-10-12)

- Allen, S.M., Cahn, J.W. 1979, A microscopic theory for antiphase boundary motion and its application to antiphase domain coarsening, Acta Metallurgica 27, 1085–1095. https://doi.org/10.1016/0001-6160(79)90196-2
- Dziarmaga, J., Laguna, P., Zurek, W.H. 1999, Symmetry breaking with a slant: topological defects after an inhomogeneous quench, Physical Review Letters 82, 4749–4752. https://doi.org/10.1103/PhysRevLett.82.4749
- Kibble, T.W.B. 1976, Topology of cosmic domains and strings, Journal of Physics A 9, 1387–1398. https://doi.org/10.1088/0305-4470/9/8/029
- Kibble, T.W.B., Volovik, G.E. 1997, On phase ordering behind the propagating front of a second-order transition, JETP Letters 65, 102–107. https://doi.org/10.1134/1.567332
- Spirin, V., Krapivsky, P.L., Redner, S. 2001, Freezing in Ising ferromagnets, Physical Review E 65, 016119. https://doi.org/10.1103/PhysRevE.65.016119 (2D coarsening can stop in stripes; matches the 81–91% single-hand share at 64²)
- Zurek, W.H. 1985, Cosmological experiments in superfluid helium?, Nature 317, 505–508. https://doi.org/10.1038/317505a0

## 8 Caveats

- The 1d and 2d kz sections took 23.5 and 20.8 min under 9-way parallel load, over the 20-minute budget. The gradient_ramp L16 section took 19.6 min.
- 2D precision (±0.05 on c) cannot distinguish 0.96 from 1.0. No 3D patch run was done, so the 3D value is theory only.
- The front runs are 1D, with closed ends and a linear ramp clipped to [−1, 1]. Real fronts may be curved or patchy.
- The crossover u\* ∈ (1, 3] was identified after the runs. The pre-registered estimate was 3.5.
