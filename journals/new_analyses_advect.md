# advect2d: does Eq. 14 survive explicit eddies? (2026-10-13)

Script `advect2d.py` (predictions committed before the real run, 580c19d; the probe-time bug fix and `ADVECT_PKL` option came after, predictions unchanged). Log `advect2d.log`, 775 s, `nice -n 10`, 2 threads. No existing script, log or draft file was edited. The moon pipeline was not rerun.

## What I tested, in plain words

Section 2.3 replaces ocean eddies with a diffusivity and admits, in its "Eddy diffusion" paragraph, that this needs the wall width sqrt(D/k) to exceed the eddy size, and that larger eddies could fold and stretch walls. I put real, moving eddies into the same equation in 2D (a chaotic stirring flow that shears the field in alternating directions) and compared it with the same box and no flow but a diffusivity equal to the stirring strength measured from the flow (D_eff = 8). Three eddy sizes at the same D_eff = 8, relative to the wall width 2.83: small (0.3 wall widths), matched (0.9), big (3.6). Each arm has 192 runs (64 seeds at each of Delta = 0.5, 1, 2). The box is 7.2 coherence lengths, so these are small-box results.

## Results

Calibration held (P1): the stirring strength measured from tracers is 6.95 to 7.10 against a design 7.0 (plus the molecular 1), and 8.05 to 8.35 against 8.0 from the decay of a long wave.

**Which hand wins.** Fit P = Phi(a x Delta), where a = 1 is the draft's whole-volume rule:

| arm | a (end of sweep) |
|---|---|
| small eddies | 1.02 +- 0.12 |
| matched | 0.86 +- 0.11 |
| big eddies | 1.10 +- 0.13 |
| no flow, D = 8 | 0.93 +- 0.11 |
| no flow, D = 1 | 0.86 +- 0.10 |

All three flow arms are within 1 combined SE of the no-flow D_eff arm, all lie in [0.8, 1.2], and the majority hand at the end of the sweep is the same 80 time units later in 99.5% of the big and matched runs. One pre-registered point failed (P4): matched eddies at Delta = 2 won 64 of 64, which sits 5.6 SE above Phi(0.855 x 2) only because the SE of a probability of 1 is floored at 0.008; I read that as a floor artifact, not a misfit, but it is reported as FAILED.

**Walls.** This is where the eddy-diffusion picture fails.

- Fraction of runs with any wall left at freeze-out (more than one domain): small 1% (2/192), matched 2% (4/192), big 29% (56/192), no flow D = 8 1.6% (3/192), no flow D = 1 79% (151/192).
- Wall length per box length at freeze-out: small 0.04, matched 0.15, big 0.92; no flow D = 8 0.03; no flow D = 1 1.63. So big eddies leave a pattern about 30 times longer in wall than diffusion at D_eff predicts, though still 0.56 of the molecular-only pattern (my P10 guess that it would match the molecular arm failed).
- Removal (post hoc, exploratory, from the saved raw results; not part of the file): of the big-eddy runs with a wall at freeze-out, 14%, 25%, 50%, 80% were single-domain by 10, 20, 40, 80 time units later. The draft's Eq. 17 with D_eff = 8 puts the coarsening time near 17 units, so removal is at least 4 times slower than predicted at the 80-unit mark. It is still much faster than molecular diffusion alone (8%, 11%, 17%, 30%).
- Small and matched eddies behave like diffusion within what 192 runs can show (2 and 4 runs with walls are too few to time).

Pre-registered scorecard: held P1, P2, P3, P5, P7, P9, P11. Failed P4 (one floor-artifact point), P6 (the small-eddy half, 0.04 vs 0.03 walls, is noise; the big-eddy half held), P8 (big eddies: 71% were already single at freeze-out, so the fraction single by 80 differs from the D_eff arm by only 0.04, not 0.15), P10.

## Verdict on the assumption

For the outcome (which hand ends up winning, and so the Phi(a x Delta) rule with N set by the whole mixed volume), the eddy-diffusion assumption holds in all three regimes: a is 0.86 to 1.10 against 0.93 for the matching no-flow arm, to the precision available (differences smaller than about 0.2 in a are invisible with these 192 runs). It does not hold for walls once the eddies are bigger than the wall width (here 3.6 times): the pattern at freeze-out is much finer than diffusion at D_eff gives, with walls in 29% of runs, and the walls that are left take at least 4 times longer to clear than Eq. 16 and Eq. 17 with D_eff say. When the eddies are smaller than or about the wall width, the diffusion description of walls held in this test.

## What this means for the moon shares (qualitative; no rerun)

- The success probability that comes from Delta and the whole-volume N (Eq. 11 to 13, Table 1) is supported by this test. I see no sign it is biased.
- The part that could change is any share that depends on the ocean coarsening to one domain in the available time (case (i), the merge-in-time condition using Eq. 17 with D_h). If the real eddies are bigger than the wall width, those draws would clear more slowly than the draft assumes, so such shares could be lower, not higher. How many sampled inputs are in that regime I did not compute; the draft already says it has not shown it either. The size of the effect is not estimated: this is one point (Da = 2.7) in 2D in a small box, and the 4-times figure is a lower bound for that point, not a scaling law.
- Where the draft sets shares by the layer-at-onset rule or by vertical mixing (D_z), this 2D horizontal test says nothing.

## Suggested wording (NOT applied to journals/ija_draft.md)

Section 2.3, replace the last two sentences of the "Eddy diffusion" paragraph ("This needs the wall width ... Eq. 14 does not include."):

> This needs the wall width sqrt(D/k) to exceed the eddy size. I tested this in 2D (`advect2d.py`, pre-registered), replacing the diffusivity by a chaotic stirring flow with the same measured D_eff = 8 and eddies of 0.3, 0.9 and 3.6 wall widths in a box of 7.2 coherence lengths. The probability that the favoured hand wins followed Phi(aΔ) with a = 1.02 ± 0.12, 0.86 ± 0.11 and 1.10 ± 0.13, against 0.93 ± 0.11 with no flow at the same D_eff. For eddies of 3.6 wall widths, however, walls survived the transition in 29% of runs (1.6% with diffusion), and those that did were removed at least four times more slowly than Eqs. 16 and 17 give. The diffusion description of the outcome is therefore supported; the diffusion description of wall motion is not, once eddies exceed the wall width.

Section 4, add to the limitations about mixing (after "Treating eddies as diffusion assumes the wall width exceeds the eddy size (Section 2.3), which is not shown."):

> A 2D test with explicit eddies supports the outcome rule but shows that eddies larger than the wall width leave more walls at the transition and clear them at least four times more slowly than diffusion at the same D_eff would; any share that relies on coarsening in a finite time could be lower if real eddies exceed the wall width. This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere.

## Limits

- Precision: SE(a) is about 0.11, so differences below about 0.2 are not detectable; the matched arm's a = 0.86 is the same as the no-flow D = 1 value, and cannot be separated from 1.
- Only one eddy size above the wall width (Da = 2.7), in 2D, in a box 7.2 coherence lengths across; no 3D, no sphere, no stratification.
- The stirring flow is a kicked shear with operator splitting; flow arms use cell size 0.5, no-flow arms 1; I did not test convergence in cell size, time step or kick period.
- The Eulerian D_eff of the matched and big arms is 4% above 8 (it includes numerical diffusion of the interpolation); P1 allowed this.
- Disclosed: the first full run crashed in the wall report (probe times 5 and 10 were not on the 2-unit cadence); I had seen its outcome table before the fix. The rerun with the same seeds reproduced those numbers exactly; the log is the rerun.
