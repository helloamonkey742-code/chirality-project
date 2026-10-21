# h14b_advect: how much slower is wall removal when eddies are much larger than the wall? (2026-10-21)

Script `h14b_advect.py` (copy of the physics of `advect2d.py`, which was not edited), predictions `h14b_advect_predictions.txt` (v1, written and hashed before any results run, sha256 `d077b49b367a6a7160696d767a2a19dc97497974dc3e61611264144d83aaa869`, asserted at the top of the script). Log `h14b_advect.log` (full pre-registered matrix, 691 s, 4 workers, exit 0), raw results `h14b_advect_raw.pkl`. One post hoc script, `h14b_advect_posthoc.py`, log `h14b_advect_posthoc.log`, reads the raw file; everything from it is marked POST HOC below. Nothing was committed or sent. `scorecard.py`, `README.md`, `journals/ija_draft.md`, `EXPERIMENT.md` and `droplet_prereg.md` were not touched. No paper was downloaded for this unit.

## What I asked, in plain words

Open item 14 (`ija_draft.md:252`) and E86 (`:187`) say eddies larger than the wall width slow wall removal, but the only measurement is `advect2d.py` at one eddy size (l/w = 3.6, one box of one eddy). The sibling analysis `h14_eddy` found that l/w reaches 4 to 9 (median) and 35 to 178 (90th percentile; 10 for Early Earth) at U = 0.02 m/s (`h14_eddy.log:60-74`), and gave, per world, the uniform slowdown f of the sideways coarsening time at which 5% of the {l > w} draws would be lost (the break-even factor, 100 to 10000 across the 15 worlds, lowest for Ganymede, Titan and Callisto at 100 and 200 for Enceladus, Europa and Mimas, `h14_eddy.log:81-171`). So the question was: how much slower is wall removal at l/w = 10 and 30, at the same D_eff = 8, and does that stay under each world's break-even factor?

I measured the slowdown S as the Kaplan-Meier median (S50) or 80th percentile (S80) of the time from freeze-out to a single domain, among runs that still have a wall at freeze-out, divided by Eq. 17's time with D_eff = 8, t17(L) = (L/5.5)^2/8 (`ija_draft.md:493`). That is the same object as the sibling's f, so the two can be compared directly.

## What ran

Same equation and noise as `advect2d.py`, same D_eff = 8 (D_mol = 1 plus eddy 7), same flow strength (kappa = 1.47, the B arm). Eddy size l = L/(2 pi n); l/w = 3.6, 10, 30 at n = 1 means box L = 64, 178, 533 and kick period Tk = 4, 31, 277 (Tk grows as l^2 at fixed D_eff). Runs spread over Delta = 0.5, 1, 2; hold = time after freeze-out (`h14b_advect.log:43-55`):

| case | l/w | n | box L / cells N / cell h | runs, walled at freeze-out | t17 | hold | note |
|---|---|---|---|---|---|---|---|
| A | 3.6 | 1 | 64 / 128 / 0.5 | 192, 56 | 16.9 | 160 t.u. | self-check against `advect2d.log:49-51` |
| B | 3.6 | 1 | 64 / 64 / 1 | 192, 62 | 16.9 | 204 | resolution partner of A |
| B2 | 3.6 | 1 | 64 / 64 / 1 | 192, 59 | 16.9 | 204 | noise control, not pre-registered |
| C | 3.6 | 2 | 128 / 128 / 1 | 96, 91 | 67.7 | 814 | |
| D | 3.6 | 4 | 256 / 256 / 1 | 24, 24 | 270.8 | 2168 | |
| E | 3.6 | 6 | 384 / 384 / 1 | 9, 9 | 609.3 | 3656 | |
| F | 10 | 1 | 178 / 356 / 0.5 | 36, 36 | 130.9 | 1048 | |
| G | 10 | 1 | 178 / 178 / 1 | 96, 96 | 130.9 | 1048 | resolution partner of F |
| H | 10 | 2 | 356 / 356 / 1 | 12, 12 | 523.7 | 3144 | |
| I | 30 | 1 | 533 / 533 / 1 | 18, 18 | 1173.9 | 9392 | |
| J | none | 1 | 64 / 64 / 1, D = 1 | 96, 76 | 16.9 | 272 | molecular reference |
| K | none | 1 | 178 / 178 / 1, D = 1 | 24, 24 | 130.9 | 2096 | molecular reference |

Disclosed, all made after the predictions were hashed and before any results run (docstring of `h14b_advect.py`): (1) a fused numba kernel replaced the numpy update (noise-off agreement 2.8e-9, asserted in the self-check, `h14b_advect.log:1`), and all cases except A and B draw in-kernel xorshift noise (unit variance), with B2 as the control; (2) case I has 18 runs, not 16; (3) finished runs are dropped from the batch, which changes the random stream after freeze-out, so A is a new realisation, not a replay; (4) the sibling's break-even f is the first ladder value with at least 5% lost (ladder at `h14_eddy.py:53`), so the largest safe ladder value is the one before it, and the table shows both rules; (5) rules for undefined (censored) quantiles were not spelled out in the predictions and were fixed in the script before the run. Not re-measured: the realised D_eff of the l/w = 10 and 30 flows (only the formula A^2/(8 Tk) is asserted; `advect2d.py` measured it only at l/w up to 3.6).

## Results

**Self-check (P1 HELD).** Case A reproduces the advect2d B point: walled 0.292 (192 runs, 56 walled; `h14b_advect.log:44`) against 0.292 in `advect2d.log:51` (single-domain fraction 0.708 at t_after = 0), wall length per box length at freeze-out 0.92 against 0.92 (`advect2d.log:49`), single by 80 time units 0.932 against 0.943 (`h14b_advect.log:59`, `:100`).

**Slowdown, in units of Eq. 17's time** (`h14b_advect.log:57-81`; brackets are 68% of the bootstrap over runs):

| l/w | case | S50 | S80 | S_half (all runs) |
|---|---|---|---|---|
| 3.6 | A (h 0.5) | 2.13 [1.89-2.72] | 5.32 [4.73-5.55] | 2.55 |
| 3.6 | B (h 1) | 2.01 [1.77-2.60] | 4.37 [4.02-6.38] | 2.07 |
| 10 | F (h 0.5) | 5.76 [5.16-6.12] | > 8 (not reached) | 3.38 |
| 10 | G (h 1) | 6.27 [6.05-6.57] | > 8 (not reached) | 3.67 |
| 30 | I (h 1) | > 8 (not reached) | > 8 | 2.15 |
| no flow | J (L 64) | 9.33 [8.51-10.28] | > 16 | 7.46 |
| no flow | K (L 178) | 13.36 [11.26-15.22] | > 16 | 2.22 |

- The 3.6 point reproduces: S50 about 2.0 to 2.1 and S80 about 4.4 to 5.3. The "four times" of E84/E86 is the 80th percentile at that size, not the median.
- Resolution does not matter within the spread: S50(B)/S50(A) = 0.94 [95% 0.55-1.47] (P2), S50(G)/S50(F) = 1.09 [0.86-1.34] (P3, `h14b_advect.log:101, 103`). The xorshift noise gives the same answer as numpy noise: S50(B2)/S50(B) = 1.06 (`:102`).
- At l/w = 10 the median wall takes about 6 times longer than Eq. 17 gives, and the 80th percentile is more than 8 (only 67% of walled runs were clear by 8 times Eq. 17's time; POST HOC, `h14b_advect_posthoc.log:3`).
- At l/w = 30 the quantiles are not reached: 14 of 18 runs still had a wall at 8 times Eq. 17's time (`h14b_advect.log:53`; 4 cleared, at 5.0, 7.4, 7.5 and 7.8 times, `h14b_advect_posthoc.log:9`). The fraction cleared by 8 times is 0.22 [0.06-0.44] against 0.29 [0.12-0.46] for no flow at all (K), so within 18 and 24 runs this case cannot be told from molecular diffusion alone (POST HOC, `h14b_advect_posthoc.log:6, 8`). The reason is built into the setup: Tk = 277 time units against a 200-unit sweep, so no kick acts before freeze-out, and about 34 kicks act in the hold (`h14b_advect_posthoc.log:29`).

**Exponent.** The pre-registered fit over B, G, I is undefined because S50 at 30 is censored (`h14b_advect.log:84-85`). The uncensored S_half gives p = 0.01 [-0.19, 0.29] (`:86`), but S_half is the half-life of the total wall length, which many short walls dominate, not the time for the last wall to go; the pre-registered slowdown is S50 and S80. POST HOC, between 3.6 and 10 only: p = 1.11 [68% 0.85-1.24; 95% 0.73-1.54] from G/B, and 0.97 [0.54-1.28] from F/A (`h14b_advect_posthoc.log:12-13`). From 10 to 30 only a floor exists (p >= 0.22, `:14`). The sibling's stress powers (p = 1.1 and p = 2, `h14b_advect.log:87-88`) are therefore not separated: p = 1.1 sits at the centre of the measured range, and p = 2 is outside the 95% range of the 3.6 to 10 piece (0.73-1.54), which says nothing about 10 to 30.

**Dependence on the number of eddies in the box.** At l/w = 3.6, S50 is 2.01 (n = 1), 1.36 (n = 2), 2.17 (n = 4), 1.51 (n = 6), all within 1.4 to 2.2 (`h14b_advect.log:91-94`). The pre-registered prediction that n = 4 falls to 0.75 of n = 1 failed: ratio 1.08 [0.53-1.67] (P8, `:97`, `:108`). At l/w = 10, n = 2 over n = 1 is 0.54 [0.29-0.72] (H over G, POST HOC, `h14b_advect_posthoc.log:15`), and 0.68 [0.44-1.05] at 3.6 (`:16`). So the slowdown does not look like a rate that holds for any box size; two eddies per box clear faster than one, but beyond n = 2 the runs (24 and 9) do not resolve a trend. Real worlds have R/(2l) of about 2 to 7 at U = 1 and 20 to 50 at U = 0.02 (`h14b_advect_predictions.txt`, P8), so an n = 1 box probably overstates the slowdown, but this is not shown for n beyond 2.

**Molecular ceiling.** With D_mol/D_eff = 1/8, the flow cannot make wall removal slower than the no-flow value (J 9.3 and K 13.4, `h14b_advect.log:78, 80`), and G is already 0.47 of K (P9 HELD, `:109`). The l/w = 30 case sits at this ceiling. So this test cannot bound the slowdown for a world where the molecular part is about 1e-9 of D_eff (the predictions file says the bounds "say nothing yet" about that ratio).

**Pre-registered scorecard** (`h14b_advect.log:100-111, 145, 147`):

| | result |
|---|---|
| HELD | P1 self-check, P2 and P3 resolution, P7 walled at freeze-out (G, F, I all 1.00), P9 molecular reference, P12 runtime (691 s, 4 workers, nothing stopped) |
| FAILED | P8 (n = 4 over n = 1 is 1.08, not <= 0.75); P10 (0 of 15 worlds COVERED at both U) |
| UNDECIDED (censored) | P4 (F: S50 5.76 is inside [1.8, 6.0] but S80 > 8 cannot be tested against [3, 12]); P5 (I: S50 > 8 against [1.5, 8]); P6 (p undefined; and p = 1.1 is not excluded, which the prediction required); P11 (needs S50 at I) |

**Per-world verdicts** under the pre-registered rule (COVERED if the upper bound at the 90th percentile l/w is under the break-even f; NOT COVERED if the lower bound at the median is at or above it): all 30 rows (15 worlds, U 0.02 and U 1) are UNDECIDED, because S80 is not reached (upper bound unbounded) and the lower bound is only the censoring floor 8 (`h14b_advect.log:113-144`). Under the stricter rule (previous ladder f) also all UNDECIDED. The pre-registered decision is "narrow / keep open" (`:148`).

## POST HOC: what the numbers allow if the power law continued (not a measurement)

I extended S(l/w) = 6.27 x (l/w / 10.02)^p beyond l/w = 10, clamped at 6.27 below it (a ceiling there, since S at 3.6 is 2), and compared with the break-even f of the sibling (`h14b_advect_posthoc.log:31-65`). Medians: S stays under f_be for all 15 worlds at both U with p = 1.11 and with p = 1.54. At the 90th percentile l/w (10 to 690 across worlds and both U in the sibling's table, `h14_eddy.log:60-74`, far beyond the measured range), with p = 1.11: 12 of 15 worlds stay under f_be at U = 0.02 and 11 of 15 at U = 1; with p = 1.54 (the 97.5% bound): 11 and 5. The worlds that go over are Ganymede, Titan and Callisto at both U and Europa at U = 1 (f_be 100 or 200; `:36-65`). With an S80 stand-in (S50 x 2.18, the S80/S50 ratio of B), the counts at the 90th percentile are 12 and 7 (p = 1.11) and 10 and 5 (p = 1.54) (`:67-69`). This is an extrapolation of a power law over up to a factor 69 in l/w from a fit over a factor 2.8, at a molecular ratio of 1/8, and I do not offer it as a verdict. What it does show: if p is near 1.1 the median draw is safe by a wide margin, and the worlds at risk are the four with f_be of 100 to 200, in their tails.

## What this means

- Item 14 narrows, it does not close. Measured: the median wall takes about 2 times Eq. 17's time at l/w = 3.6, about 6 times at l/w = 10, and more than 8 times at l/w = 30 (not reached). Every world's break-even factor is 100 or more, so the measured medians are far below them where measured (3.6 and 10); the bound above 10 is missing.
- E86's x10 run (`coarsen_slow.py`) covers the measured median up to l/w = 10 (6.3) and the 80th percentile at 3.6 (4.4 to 5.3); its x4 run does not cover the median at 10. Neither is shown to cover the 80th percentile at 10 (more than 8) or the upper tail of l/w (the sibling's 90th percentile, 10 to 690), which I could not measure.
- The measurement cannot reach the realistic regime: a molecular share of 1/8 caps the slowdown near 9 to 13, and one eddy per box overstates it, so the two limits push in opposite directions. The next step that would decide the question is a run with D_mol/D_eff much smaller than 1/8 (so the flow dominates) and enough kicks before freeze-out at l/w = 30, which needs longer sweeps or a smaller kappa.

## Suggested wording (NOT applied to `journals/ija_draft.md`)

Section 2.3, "Eddy diffusion", replace the sentence "For eddies of 3.6 wall widths, however, ... once eddies exceed the wall width." and keep the last sentence:

> For eddies of 3.6 wall widths, however, walls survived the transition in 29% of runs (1.6% with diffusion; `advect2d.log:51, 54`) and took longer to clear than Eq. 17 gives: the median by a factor of about 2 and the 80th percentile by about 4 to 5 (`h14b_advect.log:58-60`). A follow-up at the same D_eff = 8 (`h14b_advect.py`, pre-registered) found a median factor of about 6 at 10 wall widths (80th percentile more than 8) and, at 30 wall widths, a median not reached within 8 times Eq. 17's time, which cannot be told apart from molecular diffusion alone in that box (`h14b_advect.log:70-76`; `h14b_advect_posthoc.log:3-8`). Between 3.6 and 10 wall widths the factor grew as the eddy size to the power 1.1 (95% range 0.7 to 1.5, `h14b_advect_posthoc.log:12`), so a power of 1.1 is not excluded and the test could not go further, because molecular diffusion is one eighth of the total there and caps the slowdown near 9 to 13. The diffusion description of the outcome is therefore supported; the diffusion description of wall motion is not, once eddies exceed the wall width.

Section 4, eddies bullet, replace "This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere." with:

> This was tested in 2D at eddy sizes of 3.6, 10 and 30 wall widths (box lengths 64 to 533, mostly one eddy per box), with molecular diffusion at one eighth of the total, and not in 3D or on a sphere. The measured median slowdown (about 2, 6 and more than 8) lies within the factor-10 case above up to 10 wall widths; the upper tail of the eddy sizes in the sampled worlds (90th percentile of l/w from 10 to 690; `h14_eddy.log:60-74`) lies outside the measured range, so a slower clearing there is not excluded.

E86 closing sentence, "H14: the slowdown is now bounded for the shares; the eddy-size question itself stays open." becomes:

> H14: the factor-10 case covers the measured median slowdown to 10 wall widths (6.3, `h14b_advect.log:72`); beyond 10 wall widths and at a molecular share much below 1/8 it is not measured.

Item 14 (`ija_draft.md:252`):

> 14. **[Partly resolved by E84 and E86; narrowed by `h14b_advect.py`.]** **Eddy scale separation not shown.** Treating eddies as diffusion needs wall width √(D/k) > eddy size; in 2D the wall-removal slowdown was measured to 30 wall widths (median about 2, 6 and more than 8 times Eq. 17 at 3.6, 10 and 30), with molecular diffusion at one eighth of the total. Open: the slowdown above 10 wall widths and at a molecular share far below one eighth, which is the realistic case, and 3D.

## Limits

- Molecular share 1/8 of D_eff, 2D, a box of one to six eddies, no stratification, no sphere. The ceiling above applies to all cases.
- D_MOL = 1 and D_EFF = 8 are hard-coded at `h14b_advect.py:66`, so this script runs one ratio, D_MOL/D_EFF = 1/8 (the real ratio is ~1e-9, `h14c_advect.py:6`); `W_WALL`, `KAPPA2` and the Eq. 17 time (16.9) all follow from D_EFF = 8. The file cannot take a comment: its sha256 is asserted in `h14d_advect.py:37-44`. h14c computes the per-case D_eff itself (`h14c_advect.py:58-59, 72`) but still takes `KAPPA2` from this module (`h14c_advect.py:43`). `h14c_predictions.txt:23` cites this line as `:76`; it is `:66`.
- l/w = 30 had 18 runs, 14 uncleared; Tk = 277 exceeds the sweep, so the flow barely acted before freeze-out. That case is a lower bound only.
- The realised D_eff of the l/w = 10 and 30 flows was not re-measured; the formula is asserted. Numerical diffusion of the interpolation gave +4% at 3.6 in `advect2d.py`.
- n = 4 and n = 6 at l/w = 3.6 have 24 and 9 runs; the n-dependence beyond 2 is unresolved.
- The extrapolation table is illustrative. p was fitted over one factor of 2.8 in l/w with a 95% range of 0.7 to 1.5.
- The sibling's break-even f is a ladder value, the first at which 5% or more are lost; a stricter reading would use the one before it (`h14b_advect.log:113-144`, `strict`). No row changes verdict.
