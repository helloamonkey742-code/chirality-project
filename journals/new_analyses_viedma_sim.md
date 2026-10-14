# Viedma simulation: what it predicts for the bench experiment

2026-10-13. Script `viedma_sim.py` (predictions committed first as 91b053d), numbers in `viedma_sim.log`.
2440 runs, 16 minutes on 2 cores. Predictions P1-P10 are in the script docstring and were not edited.

## In plain words

I built a small mechanistic model of NaClO3 grinding: crystals dissolve or grow in one shared achiral
solution (smaller ones dissolve first), grinding chips pieces off crystals, and a chipped piece is
sometimes glued back onto a crystal of the same hand. I started each slurry with a small excess of one
hand and counted how often that hand won.

- **The sqrt(N) scaling held.** Win rates at the same Delta = eps*sqrt(N) agree across a 256-fold range of N
  (250 to 64,000 crystals): chi2 = 1.27, 4 df, p = 0.87.
- **Neutral drift is rejected.** The fitted curve beats drift by a factor of 10^67 in likelihood. Drift would
  give about 53% even at Delta = 2; the model gives 86%.
- **But the paper's coefficient is wrong by about a factor of 2.** P = Phi(a*Delta) with
  **a = 0.57, 95% CI 0.51 to 0.63**. P = Phi(Delta) (a = 1) is rejected (chi2 = 157, 1 df).

Pooled win rates, 410 runs per Delta > 0:

| Delta | Phi(Delta) (paper) | Phi(0.565 Delta) (fit) | simulated | drift |
|---|---|---|---|---|
| 0 | 50% | 50% | 56% (chance; see below) | 50% |
| 0.5 | 69% | 61% | 66% | ~51% |
| 1 | 84% | 71% | 71% | ~52% |
| 2 | 98% | 87% | 86% | ~53% |

Per-N a (95% CI): 0.57 (0.45-0.70), 0.61 (0.48-0.74), 0.54 (0.42-0.67), 0.51 (0.38-0.65), 0.61 (0.38-0.86)
for N = 250, 1000, 4000, 16000, 64000. No trend.

## What it means for the bench

- **Plan's Delta values**: expect about 61%, 71% and 87% at Delta = 0.5, 1, 2, not 69%, 84%, 98%. The model
  needs Delta of about 2.3 for 90%: at the plan's N = 6x10^3, 2.5x10^4, 1x10^5 that is eps = 2.9%, 1.4%, 0.7%.
- **Distinguishing the hypotheses**: at Delta = 2 the three candidates sit at 53% (drift), 87% (Phi(a Delta),
  a ~ 0.57), 98% (paper). Eight runs per cell can tell drift (about 4 wins in 8) from the other two (about 7 in 8), but
  not the last two from each other (11 points apart; about 40 runs at Delta = 2 for a 2-sigma split). Delta = 3 would help: 96% vs 100%.
- **Time**: every run finished (0 undecided). Time to a single hand grows only slowly with N (about +1.5 time
  units per 4x crystals; 1.8x between N = 250 and 64,000), consistent with Saito & Hyuga 2009 (time to
  homochirality does not grow with system size under grinding).

## Sensitivity (most uncertain parameters)

The same-hand re-attachment ("chiral feedback") has no source for its form or strength, so it was varied.
It is needed: Ricci, Stillinger and Debenedetti 2013 (J. Chem. Phys., doi:10.1063/1.4827478) say size-
dependent solubility alone cannot reproduce Viedma ripening, and this model agrees.

| Variant (N = 1000) | decided | a (95% CI) |
|---|---|---|
| feedback off (S = 0) | 0 of 400 by t = 60 (mass fractions drift 0.10 to 0.88 only) | not defined |
| S = 0.5 (post-hoc) | 400 / 400 | 0.41 (0.30-0.53) |
| **S = 2 (base)** | 400 / 400 | 0.61 (0.48-0.74) |
| S = 8 (post-hoc) | 400 / 400 | 0.55 (0.43-0.67) |
| grinding 4x slower (LAMBDA = 1) | 400 / 400 | 0.32 (0.21-0.43) |

Reading: a stays between 0.3 and 0.6 in every variant that amplifies at all, so the qualitative result
(sqrt(N) form yes, coefficient below 1, drift no) is stable; the number a is not. Slower grinding lowers a
by about half (that prediction P9 failed: I expected a to be insensitive to the grinding rate).

## Predictions scorecard (details in docstring)

Held: P1 (Delta = 0 is 50/50 at every N), P2 (no undecided), P4 (collapse), P6 (drift rejected), P7, P8.
Failed: P5 (a = 0.57, not 0.9; Phi(Delta) rejected), P9 (a is sensitive to grinding rate), P3 partly (Delta = 1
at 71% just under the stated floor of 72%, Delta = 2 at 86% under 90%). My prior was too close to the paper's
Phi(Delta); the simulation says a smaller coefficient.

## Honest limits

- **Model choices drive the answer.** Assumptions with no measurement behind them: reaction-limited growth,
  area-weighted chip re-attachment probability S*xA/(1+S*xA), breakage rate proportional to r^2, chip size
  uniform 5-50%, lognormal spread 0.25, LAMBDA = 4 (calibrated only to "ripening alone is at least 10x slower
  than grinding"). The size of a shifted by 2x within the variants I tried; other model forms could move it
  more. The sqrt(N) structure is close to built in for any model where a large, fast gain meets counting
  noise, so its survival is weak evidence for the model and only mild evidence about NaClO3.
- **a < 1 has a clear source here**: initial size-sampling noise alone gives a = 1/0.87 = 1.15, and breakage
  noise during amplification lowers it. Bench size matching (same sieve fraction, measured distributions)
  matters because of this: a narrower size spread would raise a.
- **Delta = 0 majority-label rate was 56% (231 of 410)**, CI just excluding 50%. A re-draw of 1800 Delta = 0
  runs with unused seeds gave 50.8%. Treated as chance, but it is a post-hoc check that is not in the
  registered predictions.
- **Not modelled**: crystal agglomeration of whole crystals, chiral impurities, size differences between
  hands (van Dongen 2026, Steendam 2013), nucleation, temperature gradients, bead dynamics, readout error,
  sampling a subset of the slurry. N is the starting count.
- This is a prediction for the bench with the model as written, **not evidence about oceans or the origin
  of homochirality.**

## Reproduce

`OMP_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python viedma_sim.py` (about 16 minutes, rewrites the log;
`--extra` appends the post-hoc S = 0.5 and S = 8 runs, about 1 minute). Seeds are `SeedSequence([20261013, N,
round(10*Delta), rep, variant])`; the per-run table is at the end of the log.
