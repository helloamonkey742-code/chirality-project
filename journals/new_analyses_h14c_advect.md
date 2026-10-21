# h14c_advect: does the eddy slowdown of wall clearing grow as molecular diffusion shrinks? (2026-10-21)

Script `h14c_advect.py` (imports the physics of `h14b_advect.py` and `advect2d.py`; the `h14b_*` files were not edited), predictions `h14c_predictions.txt` (v1, written and hashed before the results run, sha256 `522fb6ea19a961bc877d182e6e9f668d96d9f91efc41cc5258d96006b1ce79ea`, asserted at the top of the script). Run line, exactly as given:

```
cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14c_advect.py > h14c_advect.log 2> h14c_advect.err
```

Log `h14c_advect.log` (170 lines, exit 0, `h14c_advect.err` empty), raw results `h14c_advect_raw.pkl`. One post hoc script, `h14c_advect_posthoc.py`, log `h14c_advect_posthoc.log`, reads the raw file; everything from it is marked POST HOC below. Nothing was committed, pushed, or sent; no network. The only existing file edited is `advect2d.py:246`. `ija_draft.md`, `README.md`, `scorecard.py` and the `h14b_*` files were not touched.

## What I asked, in plain words

`h14b_advect.py` measured how much slower a wall is cleared when the eddies are larger than the wall, but at one ratio of molecular to total diffusivity, D_mol/D_eff = 1/8, and with the l/w = 30 flow barely acting before freeze-out (first kick 77 time units after it), so it could only say "S50 > 8". The real ratio is about 1e-9, and the no-flow reference scales as D_eff/D_mol (h14b J 9.33, K 13.36), so if the slowdown grows with the ratio, the h14b verdicts do not carry to real worlds. This run asks: (1) what is S50 at l/w 30 when the flow acts (10 kicks across the sweep), uncapped; (2) how does S50 change with D_mol/D_eff = 1/8, 1/32, 1/64 at l/w 3.6 and 10; (3) is the answer under each world's break-even factor f of `h14_eddy.log:81-171` (100 to 10000).

S is as in h14b: the Kaplan-Meier median (S50) or 80th percentile (S80) of the time from the end of the sweep (freeze-out) to a single domain, among runs still walled at freeze-out, divided by Eq. 17's time t17 = (L/5.5)^2/D_eff. t17 is 16.9, 130 and 1174 time units at l/w 3.6, 10 and 30 for every D_eff. Bootstrap over runs, 2000x.

## The nk fix (`advect2d.py:246`)

`nk = max(2, 2 * int(round(100 / flow["Tk"] / 2)))`. Before: Tk >= 100 gave nk = 0 and a divide by zero in `tracer_lagrangian`. Check (asserted in the self-check, `h14c_advect.log:1`, printed at `:49`), nk old -> new: S (Tk 0.10) 1000 -> 1000, M (0.25) 400 -> 400, B (4.0) 24 -> 24, Tk 31: 4 -> 4, and unchanged for all 499 values of Tk on the 0.2 grid below 100; only Tk >= 100 changes (Tk 100: 0 -> 2, Tk 277.4: 0 -> 2). `git diff advect2d.py` is that one line. The fixed function gave the realised eddy diffusivities below.

## What ran

Same equation, noise and flow family as h14b (kicked random-phase alternating sine shear, kappa = A/l = 1.47, n = 1 eddy per box), D_mol = 1 and cell h = 1 for all D_eff, so the molecular wall is one cell (tanh half-width 1.41 cells) and the effective wall sqrt(D_eff) is 2.83, 5.66 and 8.0 cells at D_eff = 8, 32, 64. One h = 0.5 partner (V32_3.6h) gives the molecular wall 2 cells and the effective wall 11.3 cells. Euler D dt/h^2 = 0.2 (limit 0.25). D_eddy = D_eff - D_mol by construction. Per-case gamma = min(0.01, round(1/(5 Tk), 4)), which puts 5 kicks between lam = 0 and the end of the sweep (10 across the sweep) at l/w 10 and 30 and fixes the pre-freeze-out time at 1.05 to 1.18 t17, removing h14b's confound (5.9, 0.77, 0.085 t17 at l/w 3.6, 10, 30). The bias g uses the actual gamma.

| case | D_eff | D_mol/D_eff | l/w | box L / cells N / h | gamma | Tk | runs, walled at freeze-out | hold (t17) | cpu s (`log:61-72`) |
|---|---|---|---|---|---|---|---|---|---|
| R8_3.6 | 8 | 1/8 | 3.6 | 64 / 64 / 1 | 0.01 | 4.0 | 96, 30 | 12 | 2 |
| R8_10 | 8 | 1/8 | 10 | 178 / 178 / 1 | 0.0065 | 31 | 96, 96 | 24 | 52 |
| V32_3.6 | 32 | 1/32 | 3.6 | 128 / 128 / 1 | 0.01 | 3.6 | 192, 167 | 24 | 19 |
| V32_3.6h | 32 | 1/32 | 3.6 | 128 / 256 / 0.5 | 0.01 | 3.6 | 48, 42 | 24 | 55 |
| V32_3.6m (control) | 32 | 1/32 | 3.6 | 128 / 128 / 1 | 0.0556 | 3.6 | 96, 96 | 24 | 7 |
| V32_10 | 32 | 1/32 | 10 | 355 / 355 / 1 | 0.0072 | 27.8 | 96, 96 | 30 | 326 |
| V32_30 | 32 | 1/32 | 30 | 1066 / 1066 / 1 | 0.0008 | 250.6 | 8, 8 | 26 | 1766 |
| V64_3.6 | 64 | 1/64 | 3.6 | 181 / 181 / 1 | 0.01 | 3.6 | 96, 96 | 40 | 23 |
| V64_10 | 64 | 1/64 | 10 | 503 / 503 / 1 | 0.0073 | 27.4 | 42, 42 | 48 | 335 |
| J32, J64, K32 (no flow) | 32, 64, 32 | 1/32, 1/64, 1/32 | none | 128, 181, 355 / 128, 181, 355 / 1 | 0.01 | - | 96, 96, 12 (all walled) | 80, 150, 90 | 23, 85, 143 |

Compute: 718 s wall on 4 workers (2,836 cpu-s summed over `h14c_advect.log:61-72`; the pre-registered estimate was 7,600 cpu-s, about 32 min on 4 workers, with a 50 min guard), no job cut by the guard (`h14c_advect.log:26-47`, `:167`). 1/64 at l/w 30 was not run (L = 2133, about 19,000 cpu-s).

Disclosed (also in `h14c_predictions.txt`, DISCLOSURE): timing smokes printed only the number of runs still walled, no S, p or q; a Lagrangian D_eddy (within 1.7% of nominal) and kick-damping diagnostic of the 7 flows; per-case gamma, 1/64 at l/w 30 dropped, 42 runs at 1/64 l/w 10 and 8 runs at 1/32 l/w 30 were chosen before hashing. After hashing and before the results run, three changes to the script only (none to a band): the hash constant, a normalisation of numpy booleans in the verdict function (a failing numpy comparison would not have been detected), and the damping diagnostic using the kick that `run_case` uses. Also before the results run, a dry run of `report()` on tiny fake cases (box 12 to 24, 2 runs, hold 20; scratchpad script, not in the repo) to check the report path; its output on fake data was not used. No outcome statistic of any design case was computed or seen before the results run.

## Results

**Realised flows (P2 HELD).** The Lagrangian D_eddy, from the fixed `tracer_lagrangian`, is 1.001 to 1.011 of nominal for all 7 flows, and the cubic-kick numerical diffusion is at most 0.0036 of D_mol (`h14c_advect.log:51-58`). So the flows at l/w 10 and 30 do have the intended D_eff (h14b only asserted the formula).

**Self-checks.** R8_3.6 reproduces the h14b B point: walled 0.312 against 0.323, S50 2.95 against 2.01, S80 4.84 against 4.37 (P1a HELD, `:97-98`, `:120`; R8_3.6 has 30 walled runs, so its S50 interval is wide, 1.18 to 3.54). R8_10 gives S50 5.76 against h14b G 6.27, ratio 0.92 (P1b HELD, `:83`, `:121`; gamma 0.0065 against 0.01). Resolution at 1/32, l/w 3.6: S50(h = 1)/S50(h = 0.5) = 0.88, walled 0.870 against 0.875 (P3 HELD, `:123`).

**Slowdown, in units of Eq. 17's time** (`h14c_advect.log:74-98`; 95% bootstrap in brackets):

| l/w | D_mol/D_eff | case | S50 | S80 | S_half |
|---|---|---|---|---|---|
| 3.6 | 1/8 | R8_3.6 | 2.95 [1.18-3.54] | 4.84 [3.54-7.09] | 2.13 |
| 3.6 | 1/32 | V32_3.6 | 3.54 [2.95-4.14] | 6.62 [5.91-7.44] | 2.56 |
| 3.6 | 1/32 | V32_3.6h (h 0.5) | 4.02 [2.95-4.49] | 6.14 [5.44-8.04] | 2.35 |
| 3.6 | 1/64 | V64_3.6 | 5.55 [4.61-7.45] | 9.34 [8.51-9.81] | 3.64 |
| 3.6 | 1/32 | V32_3.6m (control, 1.06 t17 before freeze-out) | 8.74 [7.92-9.93] | 12.64 [10.75-14.06] | 5.30 |
| 10 | 1/8 | R8_10 | 5.76 [4.94-6.44] | 8.86 [7.63-9.78] | 3.54 |
| 10 | 1/32 | V32_10 | 10.27 [9.48-11.32] | 13.62 [12.71-16.04] | 6.38 |
| 10 | 1/64 | V64_10 | 13.21 [11.80-14.28] | 17.22 [14.85-20.10] | 7.71 |
| 30 | 1/32 | V32_30 (8 runs) | 11.84 [10.52-14.24] | 14.41 [11.84-14.61] | 7.87 |
| none (box of l/w 3.6) | 1/32 | J32 (no flow) | 46.56 [40.18-58.96] | undefined, > 80.1 | 15.67 |
| none (box of l/w 10) | 1/32 | K32 (no flow) | 39.12 [29.39-63.37] | 79.88 [40.74-undefined] | 2.93 |
| none (box of l/w 3.6) | 1/64 | J64 (no flow) | 116.53 [96.56-undefined] | undefined, > 150.1 | 20.53 |

- **l/w 30 is no longer capped.** All 8 runs cleared, between 10.05 and 14.61 t17 (POST HOC, `h14c_advect_posthoc.log:15-23`); S50 = 11.84, against "> 8" in h14b (P7a HELD). It is not "> 25" as my guess (40) had it, so P7 as a whole FAILED. The 8 runs are few; the S80 interval (11.84 to 14.61) is set by the observed clear times.
- **At l/w 10 the median grows with the ratio**: 5.76, 10.27, 13.21 at 1/8, 1/32, 1/64; the 80th percentile 8.86, 13.62, 17.22. All walled runs cleared within the hold (cleared by 24 t17: 1.00, 0.99, 1.00 at 1/8, 1/32, 1/64; `:84`, `:80`, `:78`).
- **The flow is still well below the molecular ceiling**: flow / no flow = 0.26 at l/w 10 (V32_10 over K32), 0.08 at l/w 3.6 (V32_3.6 over J32), 0.05 (V64_3.6 over J64) (`:108-110`; P10 HELD). The ceilings scale with D_eff/D_mol as predicted (J32 46.6, J64 116.5, K32 39.1).
- **Wall removal at ratio 1/64, l/w 10 takes 13 times Eq. 17's time at the median and 17 times at the 80th percentile.** The smallest break-even f in the 15 worlds is 100 (`h14_eddy.log:97`), so every measured S is at least 5.8 times under it.

**Ratio exponent q**, slope of ln S50 against ln(D_eff/D_mol) over D_eff = 8, 32, 64 (`:105-111`): q(l/w 10) = 0.40 [95% 0.33-0.50]; q(l/w 3.6) = 0.28 [0.15-0.75]. Local exponents at l/w 10: 0.42 (8 to 32) and 0.36 (32 to 64), with no sign of saturation; at l/w 3.6: 0.13 [-0.06, 0.79] and 0.65 [0.32, 1.10], so no clean power law there (`:113-116`, descriptive). P9 HELD (q(l/w 10) in [0.1, 1.0], q(l/w 3.6) in [-0.1, 0.9], 2.5% bound 0.33 > 0).

**Length exponent p** (`:101-104`): p(3.6 to 10) = 0.65 [0.38-1.53] at 1/8, 1.04 [0.87-1.22] at 1/32, 0.85 [0.58-1.05] at 1/64; p(10 to 30) at 1/32 = 0.13. The pre-registered P8 needed p(3.6 to 10) in [0.5, 1.6] with p = 1.1 not excluded (HELD: 1.04, 97.5% bound 1.22) and p(10 to 30) >= 0.5 (FAILED: 0.13), so P8 FAILED.

**Pre-registered scorecard** (`h14c_advect.log:119-132`, `:166-170`):

| | result |
|---|---|
| HELD | P1a, P1b (self-check against h14b B and G), P2 (flows), P3 (resolution at 1/32), P4 (V32_3.6: S50 3.54, S80 6.62), P5 (V32_10: S50 10.27, S80 13.62), P6 (V64_10: S50 13.21), P9 (q), P10 (molecular ceilings), P11 (walled fractions 0.870 and 1.000), P12 (S_half/S50 0.72 and 0.62), P14 (718 s, no job cut) |
| FAILED | P7 (V32_30 S50 11.84: > 8 yes, > 25 no), P8 (p(10 to 30) = 0.13 < 0.5), P13 (illustrative worlds: U 0.02 12, inside [5, 12]; U 1 9, band <= 8) |
| UNDECIDED | none |

**Pre-registered decision: A** (the q 2.5% bound 0.33 > 0.15; `:170`). h14b's verdicts at 1/8 do not carry to the real ratio. Item 14 stays OPEN. S50(l/w 10) is 5.8, 10.3, 13.2 at 1/8, 1/32, 1/64. The ratio at which S50(l/w 10) = S50(V64_10) (ratio/64)^q would reach f = 100 is 9.9e3 (3.8e3 with q at its 97.5% bound, 3.0e4 at its 2.5% bound; `:111`). The real ratio (about 1e9) is not reachable, and no world is called covered.

**Illustrative world table (P13, not a measurement)** (`:134-165`): S_ext(lw) = S50(V64_10) (lw/10)^p at ratio 64, with p = 1.04 (point) or 1.22 (97.5%), at each world's 90th-percentile l/w, against the smaller main break-even f. Worlds below f_be: 12 of 15 at U 0.02 and 9 at U 1 (point); 10 and 5 (97.5%). At the point p the worlds above f_be are Ganymede, Titan and Callisto (f_be 100) at U 0.02, and those plus Europa, Triton and Pluto at U 1; with p at its 97.5% bound Europa and Triton join at U 0.02 (5 above), and Enceladus, Titania, Oberon and Mimas join at U 1 (10 above) (`:135-164`). This table uses p = 1.04, which POST HOC (below) turns out not to be an eddy-size exponent.

## POST HOC (`h14c_advect_posthoc.log`; not pre-registered)

**The l/w exponent is mostly a pre-freeze-out-time effect.** V32_3.6 runs 5.9 t17 and 27.8 kicks before freeze-out (gamma 0.01), while V32_10 and V32_30 run 1.07 and 1.06 t17 with 5 kicks. The control V32_3.6m has the same flow as V32_3.6 with 1.06 t17 and 5 kicks, and its S50 is 8.74 against 3.54: ratio 2.47 [2.03-2.97] (`h14c_advect.log:117`; `h14c_advect_posthoc.log:11`). With the pre-time matched, S50 is 8.74, 10.27, 11.84 at l/w 3.6, 10, 30 (`posthoc.log:3-5`), and the exponents are p(3.6 to 10) = 0.16 [95% 0.01-0.30], p(10 to 30) = 0.13 [-0.02, 0.32], and the least-squares slope over all three 0.14 [0.07-0.26] (`:7-10`), against 1.04 for the pre-registered comparison (`:6`). So at 1/32, once the sweep conditions are matched, the median slowdown looks almost independent of eddy size from 3.6 to 30 wall widths, and the "p = 1.1" of h14b (gamma 0.01 at all l/w) and of P8 here would mostly be the pre-time confound. This rests on one control cell (96 runs) at one ratio; at 1/8 and 1/64 there is no matched l/w 3.6 cell. The q estimate at l/w 10 does not have this problem: R8_10, V32_10 and V64_10 all have 5 kicks and 1.05 to 1.18 t17 before freeze-out.

**What a power law in the ratio would give (an extrapolation of 7 decades from 1 decade, not a measurement).** The ratio D_eff/D_mol at which S50(l/w 10) reaches f is, for q = 0.40, 9.9e3 (f 100), 5.6e4 (200), 5.4e5 (500), 3.1e6 (1000) and 9.4e8 (10000); the 20 (world, U) rows with f_be 100 to 1000 are passed at ratios of 1e4 to 3e6 (3.3e7 even at the q 2.5% bound), and the 10 rows with f_be 10000 at 9.4e8, about the real ratio (4.0e7 at the 97.5% bound, 3.7e10 at the 2.5% bound) (`posthoc.log:29-33`, `:38`). S50(l/w 10) at 1e9 would be 1.0e4, between 3.1e3 and 4.9e4 for q at its 2.5% and 97.5% bounds (`:35-37`). The molecular ceiling grows as the ratio (K32 39 at 1/32, about 78 at 1/64), so q <= 1 holds, but the flow is only 0.26 of the ceiling at l/w 10 and nothing here shows whether the exponent stays 0.40, falls, or the slowdown saturates. I do not offer this as a verdict.

**Plateau reading of the world table.** If the slowdown stays at its ratio-64 value 13.2 times (lw/10)^0.13 for l/w above 10 (up to the 90th-percentile l/w, 10 to 690, `h14_eddy.log:60-74`, beyond the measured 30), all 15 worlds stay under f_be at both U, also for p at its 97.5% bound 0.32 (`posthoc.log:41-42`). This is a plateau assumption outside the measured l/w and at ratio 64 only.

## What this means

- The cap is gone: at 30 wall widths, with the flow acting, the median wall takes 11.8 times Eq. 17's time (8 runs; 10.5 to 14.2), not "more than 8". At 1/64 and 10 wall widths it takes 13.2 [11.8-14.3] times, and the 80th percentile 17.2 [14.9-20.1].
- The pre-registered answer is A: the slowdown grows with the ratio of total to molecular diffusivity (q = 0.40 [0.33-0.50] at l/w 10), so the h14b "covered" statements at 1/8 are not statements about real worlds. Every measured S (up to 17.2 at S80, 13.2 at S50) is below the smallest break-even factor, 100, by a factor 5.8 or more, but the ratio reached is 64 against a real 1e9.
- E86's factor-10 case (`coarsen_slow.py`) covers the measured median at 1/8 (2.95 to 5.76) but not at 1/32 and 1/64 (10.27, 11.84, 13.21), and not the 80th percentile at 1/32 and 1/64 (13.6, 14.4, 17.2). Whether the shares of E86 move depends on f_be, which is at least 100 per world (`h14_eddy.log:81-171`); the measured values stay under it, the extrapolated ones may not.
- POST HOC, the l/w dependence beyond 3.6 appears weak at matched sweep conditions (p about 0.14), which would make the l/w extrapolation of h14b (p = 1.1) and of P13 pessimistic; the open quantity is the ratio dependence.
- The next step that would decide item 14 is a measurement of S at a ratio of 1e3 or more at l/w 10 (a cell of about 5,000 per side with a wall of 1 cell would cost about 1,000 times the 1/64 cell, so it needs a different method: a larger cell h with a sharp-interface or level-set wall, or a 1D front tracked in the eddy flow), plus a matched-pre-time l/w 3.6 cell at 1/8 and 1/64.

## Suggested wording (NOT applied to `journals/ija_draft.md`)

These replace the `new_analyses_h14b_advect.md` suggestions for the same places (apply only one of the two). All line references are to the current `ija_draft.md`.

Section 2.3, "Eddy diffusion" (`:479`), replace "For eddies of 3.6 wall widths, however, ... once eddies exceed the wall width." and keep the last sentence:

> For eddies of 3.6 wall widths, however, walls survived the transition in 29% of runs (1.6% with diffusion; `advect2d.log:51, 54`) and took longer to clear than Eq. 17 gives. In a pre-registered follow-up (`h14c_advect.py`, building on `h14b_advect.py`) the median clearing time was 3 times Eq. 17's at 3.6 wall widths and 6 times at 10 wall widths when molecular diffusion is 1/8 of the total, and, with molecular diffusion at 1/32 of the total, 10 times at 10 wall widths and 12 times at 30 wall widths (`h14c_advect.log:75, 79, 83, 97`); at 1/64 and 10 wall widths it was 13 times (80th percentile 17 times, `:77`). The factor grew with the ratio of the total to the molecular diffusivity as a power of 0.40 (95% range 0.33 to 0.50, `:107`) over a factor of 8, a range that cannot be extended to the real ratio of about 10^9. The diffusion description of the outcome is therefore supported; the diffusion description of wall motion is not, once eddies exceed the wall width, and how far it fails at the real ratio is not shown.

Section 4, eddies bullet (`:790`), replace "This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere." with:

> The slowdown was measured in 2D at eddy sizes of 3.6, 10 and 30 wall widths, with one eddy per box and molecular diffusion at 1/8, 1/32 and 1/64 of the total (`h14c_advect.log:74-98`). The factor-10 case covers the median at 1/8 (3 to 6) but not at 1/32 and 1/64 (10 to 13), and the measured factor grew with the ratio of total to molecular diffusivity as a power of 0.40. The real ratio is about 10^9 and was not reached; carried over that range the power law would give a factor near 10^4, above the break-even factor of 100 to 1000 for 20 of 30 world-and-input cases (`h14_eddy.log:81-171`), and this extrapolation is untested. The test was not done in 3D or on a sphere.

E86, end (`:187`), replace "H14: the slowdown is now bounded for the shares; the eddy-size question itself stays open." with:

> H14: the factor-10 case covers the measured median slowdown only at a molecular share of 1/8; at 1/32 and 1/64 the median is 10 to 13 times Eq. 17's (`h14c_advect.log:75-79`), growing with the ratio as a power of 0.40, still below the smallest break-even factor of the sampled worlds (100). The slowdown at the real ratio is not bounded.

Item 14 (`:252`):

> 14. **[Partly resolved by E84; measured to 30 wall widths and to a molecular share of 1/64 by `h14c_advect.py`; not closed.]** **Eddy scale separation not shown.** Treating eddies as diffusion needs wall width √(D/k) > eddy size; in 2D the wall-removal slowdown was measured at 3.6, 10 and 30 wall widths: the median was 3 to 6 times Eq. 17's time at a molecular share of 1/8, 10 to 12 times at 1/32 and 13 times at 1/64 (10 wall widths), growing with the ratio of total to molecular diffusivity as a power of 0.40 (95% range 0.33 to 0.50). Open: the slowdown at the real ratio, about 10^9, where an untested extrapolation of the power law (3 x 10^3 to 5 x 10^4) is above the smaller break-even factors (100 to 1000, 20 of 30 world-and-input cases), and 3D.

## Limits

- 2D, one eddy per box (h14b found two eddies clear faster at l/w 10, ratio 0.54), no stratification, no sphere. The ratio range is 8 to 64 against a real 1e9; q is a three-point fit (8, 32, 64) in logs, and at l/w 3.6 the local exponents differ (0.13 and 0.65).
- l/w 30 has 8 runs (all cleared, quantiles from 8 runs; one cell at 1/32 only; 1/64 at l/w 30 was not run). V64_10 has 42 runs, K32 12 (2 censored), R8_3.6 30 walled runs.
- The l/w 3.6 cells (except the control) have 5.9 t17 and 27.8 kicks before freeze-out, the l/w 10 and 30 cells 1.05 to 1.18 t17 and 5 kicks. The matched l/w exponent (0.14) rests on one control cell at 1/32; the pre-registered p (1.04) is not an eddy-size exponent for that reason.
- The molecular wall is 1 cell (tanh half-width 1.41) at h = 1; the h = 0.5 partner agrees at 1/32, l/w 3.6 (P3), and is not tested at 1/64 or at l/w 10 and 30.
- Walled fractions at freeze-out are 1.00 for every cell at l/w 10 and 30; the fraction walled at l/w 3.6 is 0.31 (1/8), 0.87 (1/32), 1.00 (1/64 and the control), so S at 3.6 is conditional on the runs that still had a wall, which depend on the pre-time.
- The sibling's break-even f is the factor at which 5% of the subset is lost (`h14_eddy.log:81-171`), on a ladder of 100 to 10000; the extrapolated world table is illustrative only.
