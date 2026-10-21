# h14d_advect: does the eddy slowdown keep growing as molecular diffusion shrinks past 1/64? (2026-10-21)

Script `h14d_advect.py` (imports `h14c_advect.py`, `h14b_advect.py` and `advect2d.py`; none of them was edited, their sha256 are asserted at the top of the script and printed at `h14d_advect.log:1`), predictions `h14d_predictions.txt` (v1, written and hashed before the results run, sha256 `5da5f82b7a392c106d9eaee9406d1ebf0d3948a7fc0b4bae1e322a1039f63a0d`, asserted at the top of the script). Run line, exactly as given:

```
cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14d_advect.py > h14d_advect.log 2> h14d_advect.err
```

Log `h14d_advect.log` (148 lines, exit 0, `h14d_advect.err` empty), raw results `h14d_advect_raw.pkl` (`h14d_advect_raw.pkl.part` is a leftover partial-save file that differs from the final pickle; the final pickle is the one used). One post hoc script, `h14d_advect_posthoc.py`, log `h14d_advect_posthoc.log` (18 lines, `.err` empty); everything from it is marked POST HOC below. Nothing was committed, pushed or sent; no network. `ija_draft.md`, `README.md`, `scorecard.py`, `advect2d.py` and the `h14b_*`, `h14c_*` files were not touched by this work (git shows the `M` on `advect2d.py`, `README.md`, `scorecard.py` and `ija_draft.md` from earlier edits, and the `advect2d.py` hash asserted here is the one with the h14c `nk` fix).

File hashes (sha256): `advect2d.py` 7d84a2e6..145cf2e, `h14b_advect.py` 9308bc49..3d5b, `h14c_advect.py` aac6365e..f6c6, `h14c_predictions.txt` 522fb6ea..79ea, `h14c_advect_raw.pkl` 1df85e4f..fc36 (full strings at `h14d_advect.log:1`). Outputs: `h14d_advect.py` 10476e92c21aa052c8134aeaef46dfba4d4c9065cb429ba7e3d7d32b01900721, `h14d_advect.log` 7d1c5f9229ecc4e42e888ccccd4e0d17f3469725b7a9b73edd955e17c2859a97, `h14d_advect_raw.pkl` c3c8f8b9b8cafe23027f5ae7232ad59e2c8e7b3aa8c357454bb0043e5ff70288, `h14d_advect_posthoc.py` d365b3fc47ca4729a5254f6dbd3f107b44ba00ec0118c0495bf833ae6f6ea720, `h14d_advect_posthoc.log` 1e0d06c2779f2d21b61212a997aa4eb36629c662a86e122e31f95cf783ec3bd8.

## What I asked, in plain words

h14c measured the eddy slowdown S50 of wall clearing at 10 wall widths (l/w = 10) at molecular shares D_mol/D_eff = 1/8, 1/32, 1/64: 5.76, 10.27, 13.21 times Eq. 17's time t17, a power law with exponent q = 0.40. Carried to the real ratio ~1e9 that exponent gives S50 ~ 1e4, far above the smallest break-even factor of the sampled worlds (100). That extrapolation fails if the growth saturates. A steady cellular flow gives q -> 0.5; chaotic stirring gives q -> 0 (at most a logarithmic growth). h14d adds measured points at 1/16 and 1/128 and a larger sample at 1/64 and asks whether S50 keeps growing as ratio^0.4 or bends over. Same model, same flow family and the same per-case gamma rule as h14c (5 kicks between lam = 0 and the end of the sweep).

Pre-registered rule: q = q(64 -> 128) = ln(S50(1/128) / S50(1/64)) / ln 2 with its 95% bootstrap interval; SATURATING if the upper bound < 0.2, GROWING if the lower bound > 0.3, otherwise UNDECIDED. Local slope: FALLING if the upper bound of Delta q = q(64 -> 128) - q(8 -> 64) is < 0, RISING if the lower bound is > 0, otherwise UNRESOLVED.

## What ran

| case | D_eff | D_mol/D_eff | l/w | box L / cells N / h | gamma | Tk | A | kappa | kicks to the end of the sweep | runs, walled at freeze-out | hold (t17) | seeds | cpu s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R8_10 (h14c, reloaded) | 8 | 1/8 | 10 | 178 / 178 / 1 | 0.0065 | 31 | 41.7 | 1.47 | 5.0 | 96, 96 | 24 | h14c | 52 |
| R16_10 (new) | 16 | 1/16 | 10 | 251 / 251 / 1 | 0.0069 | 28.8 | 58.8 | 1.472 | 5.03 | 192, 192 | 36 | 3341-3344 | 224 |
| V32_10 (h14c, reloaded) | 32 | 1/32 | 10 | 355 / 355 / 1 | 0.0072 | 27.8 | 83.0 | 1.47 | 5.0 | 96, 96 | 30 | h14c | 326 |
| V64_10 old (h14c, reloaded) | 64 | 1/64 | 10 | 503 / 503 / 1 | 0.0073 | 27.4 | 117.5 | 1.468 | 5.00 | 42, 42 | 48 | h14c | 335 |
| V64_10n (new) | 64 | 1/64 | 10 | 503 / 503 / 1 | 0.0073 | 27.4 | 117.5 | 1.468 | 5.00 | 336, 336 | 48 | 3321-3336 | 2476 |
| V64_10 = old + new (scored 1/64) | 64 | 1/64 | 10 | 503 / 503 / 1 | 0.0073 | 27.4 | 117.5 | 1.468 | 5.00 | 378, 378 | 48 | | 2811 |
| V128_10 (new) | 128 | 1/128 | 10 | 711 / 711 / 1 | 0.0074 | 27.2 | 166.2 | 1.469 | 4.97 | 336, 336 | 60 | 3301-3316 | 5354 |

(`h14d_advect.log:87-93`; flows `:82-84`.) Wall clock 2,695 s = 45 min on 3 workers under `nice -n 10`, guard 150 min (9,000 s), no job cut (`:43-78`, `:144`). The cpu total of the new jobs is 8,054 s against the 5,600 s pre-registered estimate: the machine was shared with two other agents and the 1/128 jobs ran at 245-385 s each against the 231 s implied by the estimate of 11 s per run, 15.9 s per run measured. t17 = 130.2 to 130.9 time units for all cells. The 1/128 and 1/64 cells were run with 16 jobs of 21 runs (7 seeds x 3 sweep-end biases Delta = 0.5, 1, 2), the 1/16 cell with 4 jobs of 48.

Disclosed (also in `h14d_predictions.txt`, DISCLOSURE): (i) the h14c S50 values at 1/8, 1/32, 1/64 and the local exponents were known when the bands were written; (ii) the scratch timing scripts (`calib.py`, `timing.py`, `timing2.py` in the session scratchpad, not in the repo) ran `h14b_advect.run_case` at D_eff 128 and 16 with gamma 0.1 and a 400-time-unit hold; their output was discarded except elapsed seconds, steps, ns per cell-step and the number of runs still walled at the end of that hold, and no S50, S80, q or walled fraction at the design gamma and hold was computed or seen for an h14d cell; (iii) the pooling of extra 1/64 runs with the old 42, the run counts, the holds and the guard were chosen before hashing. After hashing, before the results run: the script was written, the cosmetic parts of it fixed (a dead loop expression, a variable name and the wording of one printed diagnostic), and a smoke mode run that executed only the self-checks and the reproduction of the three h14c values (`h14d_advect.log:2-4`); no S50 of an h14d cell was computed before the results run. The Lyapunov exponent (0.13 per kick) was first printed by that self-check, after the hash; the P8 threshold (> 0.1) was set before it, the guess (0.3) was wrong.

Self-checks (`:2-4`): h14c's `self_check` (nk fix, bias from the actual gamma, tiny `job()`/`pool()`), seeds unique and outside 3101-3132, D_eddy = D_eff - D_mol by construction, Euler stability, kicks to the end of the sweep within 5 +- 0.7, the new 1/64 jobs carry exactly the h14c V64_10 flow, gamma, L, hold, cadence and Deltas (so pooling is legitimate), and the reloaded h14c cells reproduce 5.76, 10.27 and 13.21 with the h14c bootstrap seeds.

## Results

S50, S80 and S_half at l/w 10 (time since the end of the sweep / t17; Kaplan-Meier over runs walled at freeze-out, every cell 100% walled and fully cleared; bootstrap over runs, 2000x; `h14d_advect.log:96-108`):

| D_mol/D_eff | runs | S50 [95%] | S80 [95%] | kicks to clear (S50 t17 / Tk) | mean domains at freeze-out | log line |
|---|---|---|---|---|---|---|
| 1/8 | 96 | 5.76 [4.94-6.44] | 8.86 [7.63-9.78] | 24 | 3.4 | `:96-97` |
| 1/16 (new) | 192 | 7.71 [7.34-8.40] | 12.02 [10.75-13.00] | 35 | 4.9 | `:98-99` |
| 1/32 | 96 | 10.27 [9.48-11.32] | 13.62 [12.71-16.04] | 48 | 7.9 | `:100-101` |
| 1/64 old (h14c) | 42 | 13.21 [11.80-14.28] | 17.22 [14.85-20.10] | 63 | 12.2 | `:102-103` |
| 1/64 new | 336 | 13.05 [12.71-13.85] | 16.43 [16.06-17.26] | 62 | 12.5 | `:104-105` |
| 1/64 pooled (scored) | 378 | 13.05 [12.74-13.80] | 16.48 [16.11-17.26] | 62 | 12.5 | `:106-107` |
| 1/128 (new) | 336 | 15.28 [14.93-15.79] | 18.79 [17.86-20.00] | 73 | 22.8 | `:108-109` |

Local ratio exponents q(a -> b) = ln(S50_b / S50_a) / ln(r_b / r_a) (`:111-121`):

| step | q [95%] | h14c, same step | log line |
|---|---|---|---|
| 1/8 -> 1/16 | 0.42 [0.24-0.70] | | `:112` |
| 1/16 -> 1/32 | 0.41 [0.24-0.56] | | `:113` |
| 1/32 -> 1/64 | 0.34 [0.21-0.50] | 0.36 [0.15-0.53] (1/32 -> 1/64 with 42 runs) | `:114` |
| **1/64 -> 1/128** | **0.23 [0.13-0.28]** | | `:115`, `:119` |
| 1/8 -> 1/64 endpoints | 0.39 [0.34-0.48] | | `:116` |
| q(8 -> 64), least squares over 8, 32, 64 (h14c definition, pooled 1/64) | 0.40 [0.34-0.48] | 0.40 [0.33-0.50] | `:117` |
| q over all five ratios, least squares | 0.36 [0.32-0.41] | | `:118` |

Delta q = q(64 -> 128) - q(8 -> 64) = -0.17 [95% -0.32 to -0.09], paired over the same replicates (`:119`). Curvature of ln S50 against ln r over the five ratios, c = -0.05 [-0.11, -0.01] (`:121`). The 95% interval of q(64 -> 128) excludes 0.5 (the cellular limit) and excludes 0 (the saturating limit) (`:120`).

Regime (`:82-84`, `:97-109`): kappa = A/l = 1.47 for all three new flows, so the flow is in the chaotic regime (kappa >~ 1). The exact Lyapunov exponent of the tangent map of the random-phase alternating shear is 0.13 per kick (0.0046 to 0.0048 per time unit) at all three D_eff. The measured Lagrangian D_eddy is 1.1% above nominal (D_eff - 1) at all three, the numerical kick diffusion is 0.0005 D_mol. S50 clearing takes 24 to 73 kicks, 2 to 3 times the rough mixing time ln(l / 1 cell)/lambda (28 to 36 kicks). Mean domains at freeze-out grow from 3.4 to 22.8 (box L grows as sqrt(D_eff) at fixed l/w).

Predictions (`:135-146`): P1 to P10 all held. P1 R16_10 S50 7.71 in [6.3, 9.5] (guess 7.7); P2 V128_10 S50 15.28 in [12.5, 20.5] (guess 16.3), S80 18.79 in [16, 30] (guess 21); P3 q(64 -> 128) 0.23 in [0.10, 0.50] (guess 0.30); P4 Delta q -0.17 in [-0.40, 0.05] (guess -0.10); P5 ladder 5.76 < 7.71 < 10.27 < 13.05 < 15.28; P6 S50(V64_10n) / S50(V64_10 old) 0.99 in [0.85, 1.18]; P7 flow realisation (D_eddy ratio 1.011, D_num/D_mol 0.0005, kappa 1.468-1.472, kicks 4.97-5.03); P8 Lyapunov 0.13 > 0.1 (guess 0.3, off by a factor 2.3); P9 curvature -0.05 <= 0 (guess -0.02); P10 compute 2,695 s, no job cut, 864 of 864 runs. The guesses (not thresholds) of P2 (16.3 against 15.28) and P3 (0.30 against 0.23) were high, P8 low, in the direction of "growth slows".

**Decision (pre-registered, `:147`): UNDECIDED.** q(64 -> 128) = 0.23 [0.13-0.28]: the upper bound is above 0.2 (SATURATING not reached) and the lower bound is below 0.3 (GROWING not reached). Pooled, old-only (0.21 [0.08-0.38]) and new-only (0.23 [0.13-0.29]) 1/64 cells all give UNDECIDED (`:123-125`). **Local slope label (pre-registered, `:148`): FALLING**, Delta q = -0.17 [-0.32, -0.09].

Extrapolation, NOT a measurement (`:128-132`), anchored at the 1/128 point: S50 at ratio 1e9 under a power law with q at its 2.5% bound (0.13) 130, at the point (0.23) 573, at the 97.5% bound (0.28) 1.31e3; with the h14c slope (0.40) 8.2e3 (h14c got 1.0e4 from the 1/64 anchor); a constant from 1/128, 15.3. The ratio at which S50 reaches f = 100 is 1.4e8 (q 0.13), 4.8e5 (0.23), 1.0e5 (0.28) and never for a constant. By the pre-registered meaning of UNDECIDED: item 14 stays open, the five-point ladder and both envelopes are reported without choosing between them; the envelopes at 1e9 are 15 (constant) and 1.3e3 (power law at the upper bound), which straddle the smallest break-even factor of the sampled worlds (100) and the next ones (200, 500, 1000) but not the largest (10,000).

## POST HOC (`h14d_advect_posthoc.log`; not pre-registered)

- **The increments per doubling are about constant, as for a logarithmic law** (`:2-7`). S50 rises by 1.94, 2.57, 2.77 and 2.24 at the four doublings (95% intervals 1.2-3.2, 1.5-3.6, 1.8-3.9, 1.3-2.8); the last divided by the mean of the first three is 0.92. A power law with q = 0.3 to 0.4 would make each increment larger than the one before; the difference between the last two increments is -0.54 [-2.32, 0.61]. Kicks to clear rise by 10, 13, 14 and 11 per doubling (point values).
- A two-parameter log law S50 = a + b ln(ratio) fits the five points (rms residual 0.19, chi-square 1.9 on 3 degrees of freedom) with b = 3.52 per e-fold [3.27-3.84]; the power-law fit over the same five points (q = 0.36) has rms 0.48 and chi-square 18.4, with residuals of sign pattern -, +, +, +, - (`:10-12`). The chi-squares use the bootstrap width of each point as its SE and unweighted fits, so they are indicative only. The log law is a different description of the same ladder, chosen after seeing it.
- Extrapolation of the log law to ratio 1e9 gives S50 = 71 [67-77] (the power-law fit over all five points gives 4.7e3) (`:13`). It crosses f = 100 at ratio 3.7e12, f = 200 at 8e24. This extrapolates 4 measured doublings over 23 more, and it is a fit, not a mechanism; a power law with a slowly falling exponent matches the data about as well within the spread of the points, only the constant-q power law does not.
- The run-level bootstrap is valid for these data (the flow phase is drawn per run, `h14b_advect.py:162-165`), and a job-level bootstrap gives the same interval: q(64 -> 128) 0.23 [0.13-0.28] at run level and at job level (`:17`); the standard deviation of S50 between the 16 jobs at 1/128 (0.93) and the 18 jobs at 1/64 (1.12) implies SEs of 0.23 and 0.26 against run-bootstrap SEs of 0.21 and 0.30 (`:18`).

## What this means

- The pre-registered answer is UNDECIDED on the exponent and FALLING on the slope: the growth of the l/w-10 slowdown with D_eff/D_mol is no longer well described by q = 0.40 once 1/128 is included (0.23 [0.13-0.28] for the last doubling), the curvature is negative, and neither 0.5 (cellular) nor 0 is consistent with the data. The measured regime is chaotic (kappa 1.47, Lyapunov 0.13 per kick), which fits the family of "slowing growth" and not the steady-cellular picture, but the interval still allows a power law with q up to 0.28.
- The h14c extrapolation S50 ~ 1e4 at ratio 1e9 (`h14c_advect_posthoc.log:35`) is not supported by the 1/128 point (the slope that gave it is 0.40, the last measured slope is 0.23 and falling). It is replaced by a range, 15 (constant) to 1.3e3 (q at its upper bound), with POST HOC 71 for a log law. Every measured value (S80 up to 20.0 at the 95% upper bound, S50 up to 15.8) is still below 100 by a factor 5 or more; the real ratio, 1e9, was not reached.
- E86's factor-10 case (`coarsen_slow.py`) covers the median at 1/8 (5.76) and does not cover it at 1/32, 1/64 and 1/128 (10.3, 13.1, 15.3); it covers the 1/8 and 1/16 medians only if one counts 7.7 as "about 10". Against the break-even factors of the sampled worlds (rows of `h14c_advect_posthoc.log:38`: 6, 7, 3, 4, 10 of 30 at f = 100, 200, 500, 1000, 10000), the upper envelope (1.3e3) is above f = 100, 200, 500 and 1000 (20 of 30 rows open) and the constant and log-law values are below all of them. This is arithmetic on the two extrapolations, not a result about the worlds.
- What would settle the exponent: a measurement at a much larger ratio, or a method that does not resolve the molecular wall (a larger cell with a sharp-interface wall, or a 1D front tracked in the eddy flow). Cost scaling from this run (summed cpu per run 7.4 s at 1/64, 15.9 s at 1/128, so about 2.2 times per doubling on a shared machine): a 1/256 cell of 336 runs would be about 12,000 cpu-s, roughly 70 min on 3 workers if the same scaling held (an extrapolation of cost, not a timing). With the same run count its q(128 -> 256) interval would probably be about as wide as 0.13 to 0.28 (not computed).

## Suggested wording (NOT applied to `journals/ija_draft.md`)

These replace the `new_analyses_h14c_advect.md` suggestions for the same places (apply only one of the two). All line references are to the current `ija_draft.md`; log references are `h14d_advect.log`.

Item 14 (`:252`), replace the h14c suggestion with:

> 14. **[Partly resolved by E84; measured to 30 wall widths and to a molecular share of 1/128 by `h14c_advect.py` and `h14d_advect.py`; not closed.]** **Eddy scale separation not shown.** Treating eddies as diffusion needs wall width √(D/k) > eddy size; in 2D, at 10 wall widths, the median wall-removal time was 5.8, 7.7, 10.3, 13.1 and 15.3 times Eq. 17's time at molecular shares of 1/8, 1/16, 1/32, 1/64 and 1/128 (`h14d_advect.log:96-108`). The exponent of the growth with the ratio of total to molecular diffusivity was 0.40 over 1/8 to 1/64 and 0.23 (95% range 0.13 to 0.28) for the last doubling, 1/64 to 1/128 (`:117`, `:119`); the slope fell (`:119`) and neither the cellular value 0.5 nor 0 is consistent with it. Open: the slowdown at the real ratio, about 10^9, which is bracketed only by extrapolation: 15 if growth stops at 1/128 and 1.3 × 10^3 if the power law continues at its upper bound (`:128-132`); the smallest break-even factor of the sampled worlds is 100. Not tested: 1/64 and 1/128 at 30 wall widths, 3D, a sphere, more than one eddy per box.

Section 2.3, "Eddy diffusion" (`:479`), replace the h14c suggestion (the sentences from "For eddies of 3.6 wall widths, however, ..." through "how far it fails at the real ratio is not shown.") with:

> For eddies of 3.6 wall widths, however, walls survived the transition in 29% of runs (1.6% with diffusion; `advect2d.log:51, 54`) and took longer to clear than Eq. 17 gives. In pre-registered follow-ups (`h14c_advect.py`, `h14d_advect.py`, building on `h14b_advect.py`) the median clearing time was 3 times Eq. 17's at 3.6 wall widths and 6 times at 10 wall widths when molecular diffusion is 1/8 of the total (`h14c_advect.log:97`, `h14d_advect.log:96`); at 10 wall widths it was 8, 10, 13 and 15 times at molecular shares of 1/16, 1/32, 1/64 and 1/128 (`h14d_advect.log:98-108`), and 12 times at 30 wall widths with a share of 1/32 (`h14c_advect.log:75`). The factor grew with the ratio of total to molecular diffusivity as a power of 0.40 up to 1/64 (95% range 0.34 to 0.48, `h14d_advect.log:117`) and 0.23 (0.13 to 0.28) from 1/64 to 1/128 (`:119`); the real ratio is about 10^9 and was not reached. The diffusion description of the outcome is therefore supported; the diffusion description of wall motion is not, once eddies exceed the wall width, and how far it fails at the real ratio is not shown.

Section 4, eddies bullet (`:790`), replace the h14c suggestion with:

> The slowdown was measured in 2D at eddy sizes of 3.6, 10 and 30 wall widths, with one eddy per box and molecular diffusion at 1/8 to 1/128 of the total (at 10 wall widths: 1/8, 1/16, 1/32, 1/64, 1/128; `h14d_advect.log:96-108`). At 10 wall widths the median clearing time was 5.8 to 15.3 times Eq. 17's, growing as the ratio to a power of 0.40 up to 1/64 and 0.23 (0.13 to 0.28) from 1/64 to 1/128. The real ratio is about 10^9 and was not reached; the growth is bracketed at 15 (no further growth) to 1.3 × 10^3 (power law at its upper bound), against break-even factors of 100 to 10,000 for the sampled worlds (`h14_eddy.log:81-171`). The test was not done in 3D or on a sphere.

E86, end (`:187`), replace "H14: the slowdown is now bounded for the shares; the eddy-size question itself stays open." with:

> H14: the factor-10 case covers the measured median slowdown only at a molecular share of 1/8; at 1/32, 1/64 and 1/128 the median is 10 to 15 times Eq. 17's (`h14d_advect.log:100-108`), and its growth with the ratio is slowing (power 0.40, then 0.23 for the last doubling). The slowdown at the real ratio is not bounded below 100 (the upper extrapolation is 1.3 × 10^3).

## Limits

- 2D, one eddy per box, no stratification, no sphere. The l/w 10 ladder runs over 4 doublings (8 to 128) against about 27 doublings from 8 to the real 1e9; every statement about 1e9 is an extrapolation. The interval for the last doubling (0.13 to 0.28) straddles the rule thresholds 0.2 and 0.3; whether another doubling would resolve it was not computed.
- Only l/w 10 was extended. l/w 3.6 and 30 are h14c cells at 1/8-1/64 (l/w 30: 8 runs at 1/32 only). The l/w extrapolation of earlier notes is untouched.
- The 1/8, 1/32 and the old 1/64 cells are the h14c runs (96, 96, 42 runs), reloaded and re-scored with the same bootstrap seeds; their intervals are wider than those of the new cells (336 and 192 runs), and the old-only 1/64 cell gives q(64 -> 128) 0.21 [0.08-0.38].
- The molecular wall is 1 cell (tanh half-width 1.41) at h = 1 for all D_eff; the h = 0.5 partner of h14c agrees only at 1/32 and l/w 3.6. At 1/128 the box is 711 cells per side; the grid cost limits the next doubling.
- Walled fraction at freeze-out is 1.00 and every walled run cleared inside the hold in all new cells (`h14d_advect.log:87-93`), so no quantile was censored.
- The log-law reading (POST HOC) was chosen after seeing the ladder; it is a two-parameter fit to five points, not a mechanism, and a power law with falling q is not excluded.
- The pre-registered guesses for P2, P3 and P8 (not thresholds) were off in the direction of slower growth than assumed for the first two and a smaller Lyapunov exponent than assumed for the third.
