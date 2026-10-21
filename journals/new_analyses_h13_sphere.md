# H13 on a sphere: does case (i) coarsen to one hand in the time the sweep assumes? (2026-10-21)

Unit h13_sphere. Nothing in `journals/ija_draft.md`, `scorecard.py`, `README.md` or any existing file was edited; the wording in
"Proposed wording" is **not applied**. All files of this unit start with `h13_sphere`.

## Answer

1. **The share f of case (i) stays 1.** On a sphere every random patchwork reached one hand: 663 of 664 within 30 sweep times
   (`h13_sphere_score.log:1`), and the one that did not (seed 136) got there at 36.3 sweep times
   (`h13_sphere_late136.log:2`). The flat periodic control stays at two hands in 51 to 65% of runs (`h13_sphere_t32_f50.log:102`,
   `h13_sphere_t32_f55.log:102`, `h13_sphere_t57_f50.log:102`: 38, 49, 35 of 100 reach one hand). The stuck states of the
   aniso3d slab (4 of 8) are a property of the periodic box, not of the ocean.
2. **The time is what was wrong.** The sweep's sideways term tau_sw = (pi R / 5.5)^2 / D is too short. A sphere needs
   T1 = f_t * tau_sw with f_t median 7.42, p90 17.00, p99 22.84 at R = 16 (`h13_sphere_s16_f50.log:302`); in units of R^2/D
   (R^2/D = 3.065 tau_sw) that is 2.4, 5.5 and 7.5 R^2/D, against 0.33 R^2/D in the sweep.
3. **The f that makes the sweep's term match the sphere: f_t = 17 (rule D2, p90), equivalent to A = 5.5/sqrt(17) = 1.33**
   (`h13_sphere_score.log:7`); f_t = 7 matches the median draw, 23 the p99.
4. **How much the shares move: at most 0.16 points at f_t = 17, 0.32 at 23, 0.57 at 30** (Ganymede, main rule, no lid;
   `h13_sphere_shares.log:38-40`). If every draw gets its own f_t resampled from the 300 sphere runs the largest loss is 0.08
   points (`h13_sphere_shares.log:41`). The existing `coarsen_slow.log` bound (x4, x10, largest loss 0.0 and 0.1,
   `coarsen_slow.log:20`) does not reach f_t = 17 or 23, so this run is what extends it. Rule D3 sent it to p90 and p99; done.
   At f = 1 the script reproduces `lid_melt.LOGGED` to 0.05 points (`h13_sphere_shares.log:19`).

Why 7 and not 1: a sphere has no stable wall except an exact great circle, and that one is unstable at rate D/R^2, so the last
wall is a loop near a great circle that needs a time (R^2/D) ln(0.98/u0) to shrink (u0 = its offset from the equator). The
simulation reproduces this law within 0.0 to 3.8% (`h13_sphere_special.log:2-5`; u0 = 0.1 takes 584 = 6.99 tau_sw, line 3). A
typical u0 near 0.1 gives about 7 tau_sw, which is the median. The sweep's term is the time for the walls to reach spacing
pi R; the endgame after that is not in it.

## What ran

Model: da/dt = a - a^3 + D lap(a), D = 1, no bias, smoothed-noise patchwork cut at area share f0 (0.5, and 0.55 for the biased
case). Sphere: spherical harmonics on a Gauss-Legendre grid, IMEX SBDF2, dt = 0.25, lmax = 4R. Flat control: periodic square,
FFT. One hand = minority area below 1%. Self-check `python h13_sphere.py check` passes (Laplacian eigenvalues l(l+1)/R^2).

| Job | What | Log | Result |
|---|---|---|---|
| S16a | sphere R=16, f0=0.5, 300 runs | `h13_sphere_s16_f50.log` | 299/300, median 7.42, p90 17.00, p99 22.84 (:302) |
| S16b | sphere R=16, f0=0.55, 300 runs | `h13_sphere_s16_f55.log` | 300/300, median 4.96, p90 8.02, p99 14.03 (:302) |
| S32a | sphere R=32, f0=0.5, 32 runs | `h13_sphere_s32_f50.log` | 32/32, median 6.36, p90 11.39, p99 14.03 (:34) |
| S32b | sphere R=32, f0=0.55, 32 runs | `h13_sphere_s32_f55.log` | 32/32, median 4.04, p90 7.72, p99 16.14 (:34) |
| S48a | sphere R=48 | none | **NOT RUN** (CPU shared with six other jobs; R=32 took 834 and 944 s) |
| T32a/b, T57a | flat periodic square, 100 runs each | `h13_sphere_t32_f50.log`, `_t32_f55.log`, `_t57_f50.log` | 38, 49, 35 of 100 reach one hand (:102) |
| T512 | flat 512 square, 4 runs, coarsening law | `h13_sphere_t512.log` | exponent 0.413, A = 6.26 (:17) |
| SPEC | offset-circle cap law, wavy m=2 circle | `h13_sphere_special.log` | ratios 1.003, 1.000, 1.005, 1.038 (:2-5); m=2 decay rate 3.10 and 2.96 (theory 3, :6-7); minority stays 0.5000 |
| SHELL | thin 3D shell R=16, 8 layers, H'=4, 12 runs | `h13_sphere_shell.log` | 12/12, median 5.80, min 1.93, max 17.95 (:14) |
| shares | world shares at f_t = 1, 8, 17, 23, 30 and mix | `h13_sphere_shares.log` | above |
| late136 | seed 136 rerun to 60 tau_sw | `h13_sphere_late136.log` | reaches one hand at 36.30 tau_sw (:2) |

Scoring: `python3 h13_sphere_score.py > h13_sphere_score.log`. Commands are in each script's docstring (`Run:` line).

Other readings: with f0 = 0.5 the favoured hand wins 0.530 pooled (158/300, 18/32; `s16_f50.log:303`, `s32_f50.log:35`),
a fair coin; with f0 = 0.55 it wins 0.877 pooled (260/300, 31/32). The sphere and the flat torus give the same domain length
at l(20) and l(40) within 0.8% and 2.1% (`h13_sphere_score.log:3`): curvature does not matter until l approaches R.
Per world, the first draw to fail the sideways test appears at f_t = 14.4 for Ganymede, 15.0 Titan, 17.0 Callisto and 18.4
Enceladus (`h13_sphere_shares.log:2-6`); at f_t = 17, 0.056% of Ganymede draws fail it (:4).

## Predictions (`h13_sphere_predictions.txt`, sha256 b316fc86..., 2026-10-21 00:37:01, unchanged; hash asserted in every script)

| P | Claim | Verdict | Numbers |
|---|---|---|---|
| P1 | at most 1% of sphere runs not at one hand by 30 tau_sw | **held** | 1 of 664 (`score.log:1`) |
| P2 | flat control stuck shares in bands 0.25-0.55, 0.10-0.50, 0.20-0.55 | **FAILED** | 0.62, 0.51, 0.65, all above the upper edge (`score.log:2`) |
| P3 | coarsening exponent 0.45-0.55, A 4.5-8.5, sphere = torus within 15% | **FAILED** on the exponent | 0.413; A = 6.26 and sphere-torus gaps held (`score.log:3`) |
| P4 | S16a median 2-6, p90 5-14, p99 9-24; S16b median 1.5-5 | **FAILED** on median and p90 | 7.42, 17.00 too slow; p99 22.84 and S16b 4.96 inside (`score.log:4`) |
| P5 | R-independence, tail e-fold 0.7-1.4 R^2/D | **FAILED** | S32a/S16a median 0.857 inside; S16a tail 1.564 outside (S32a 0.796); S48a not run (`score.log:5`) |
| P6 | favoured hand wins 0.40-0.60 at f0=0.5, 0.80-1.00 at 0.55 | **held** | 0.530, 0.877 (`score.log:6`) |
| P7 | cap-law within 10%, m=2 decays at 2.5-3.5, minority >= 0.45 | **held** | (`score.log:8`) |
| P8 | shell: all 12 reach one hand, median within 0.6-1.7 of S16a | **held** | 12/12, 5.80, ratio 0.78 (`score.log:9`) |

4 held, 4 failed. Failed P4 means the pre-registered expectation was too optimistic about how fast a sphere finishes; the
decision rules D0-D3 were fixed before the runs and were applied as written (D1 false, D2 gives 17, D3 triggered).
Reading of the P3 miss (arithmetic on logged values, not a fit): l/sqrt(Dt) is still falling at t = 320 (9.96, 8.30, 7.31,
6.70, 6.30, 6.03 at t = 10 to 320, `t512.log:11-16`), so a single exponent over [20, 500] gives 0.413; the slope of l^2 between
t = 160 and 320 is 33.2, i.e. A = 5.76, close to the 3D 5.5 but a 2D number.

## Proposed wording (not applied)

**Item 13 (`journals/ija_draft.md:251`).** "13. **Case (i) coarsening.** Settled for the share, not for the time. On a sphere
(`h13_sphere.py`) 663 of 664 random patchworks reached one hand within 30 sweep times and the last at 36, so f = 1 stands; the
flat periodic box keeps two hands in 51 to 65% of runs and does not represent an ocean. The sideways time of the sweep is too
short by a median factor of 7 (90th percentile 17); with the factor set to 17 no share moves by more than 0.16 points
(`h13_sphere_shares.log`)."

**Paragraph at `ija_draft.md:645`, last two sentences.** "A flat wall spanning a sideways-periodic box has zero curvature and
does not move (Eq. 16), so half the slab runs kept two hands. A sphere has no stable wall: a great circle has zero mean
curvature but is unstable, and on a sphere all 664 random starts reached one hand, the last one after 36 times the sweep's
sideways time (`h13_sphere_late136.log`)."

**Limitations bullet (`ija_draft.md:789`).** Replace "That a shallow ocean coarsens to a single domain (case i) is only
partly supported in a periodic box and untested on a sphere (Section 3.2)." with "That a shallow ocean coarsens to a single
domain (case i) holds on a sphere in the same model (664 of 664 runs; a thin 3D shell, 12 of 12), but takes a median 7 times
longer than the sweep's sideways term (90th percentile 17), which I did not put into the sweep (Section 3.2)." Also soften
"Sideways eddy mixing ... never limits horizontal coarsening" in the same bullet: it holds up to a time factor of about 14
(Ganymede, `h13_sphere_shares.log:4`).

**Eddies bullet (`ija_draft.md:790`), append.** "The sphere's own endgame is a further factor of 7 (median) to 17 (90th
percentile); with 17 the largest change is 0.16 points, and 0.32 at 23 (`h13_sphere_shares.log:38-39`)."

## Caveats

- Sphere and flat runs are 2D surfaces: a shallow ocean mixed over its depth. The thin 3D shell (12 runs, H' = 4) agrees with
  the sphere (median 5.80 against 7.42, ratio 0.78) but is small. The 2D torus is only an analogue of the 3D slab of
  `aniso3d.py`; the 51 to 65% stuck fractions are not directly comparable to its 4 of 8.
- S48a was not run, so the R-independence in P5 rests on R = 16 and 32 only (median ratio 0.857, with 32 runs at R = 32). The
  R = 32 percentiles (p90 11.39) are from 32 runs and noisy.
- The `h13_sphere_t512.log` run is 4 runs to 0.1 tau_sw: it fixes the coarsening law, not a time to one hand.
- Smoke runs before the pre-registration are listed in the `h13_sphere.py` docstring. Seed 136 was also run inline with a longer
  horizon before being rerun to `h13_sphere_late136.log` (same T1 = 3032.0); that is a post-hoc look, not a pre-registered
  test, and it is the only run that needed more than 30 sweep times.
- Machine load was 28 to 38 during the runs; wall times in the logs are not benchmarks.
- No bathymetry, continents, ice-shell topography, rotation or advection are in the model; any of these could pin a wall or
  slow clearing (H14b, `advect2d.py`: at least 4x for eddies larger than a wall width). Multiplying 17 by 4 gives 68; that
  combined case was not run. `h13_sphere_shares.log:40` shows 0.57 points at 30, so 68 would be larger than anything tested.
- Open-access literature: none was needed, none cited.

## Reproduce

`UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python h13_sphere.py random kind=sphere size=16 f0=0.5 n=300 seed=0 procs=4 > h13_sphere_s16_f50.log`
(other jobs: see the `Run:` lines of `h13_sphere.py`, `h13_sphere_shell.py`, `h13_sphere_shares.py`, `h13_sphere_score.py`).
