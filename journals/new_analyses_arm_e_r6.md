# arm_e_r6: R6 for "Transamination + ligation network, fed (2025)" (scorecard.py:60)

Written 2026-10-21. Nothing in this file has been applied to scorecard.py, README.md, journals/ija_draft.md, EXPERIMENT.md or journals/droplet_prereg.md. Proposed wording is at the end, marked "not applied".

Short answer. The R6 "?" depends on two readings of the existing rule, and the notes give both.

| | measured rates (F = 1) | crossing multiplier (F = 1.5 F_c) |
|---|---|---|
| inherited excess in the feed, literal rule (Delta >= 2) | **Y**, but only because the network passes the excess through (arm_e_r6.log:29-32) | **Y** (arm_e_r6.log:99-102) |
| inherited excess in the feed, strict rule (also >= 1% output) | **N** (arm_e_r6.log:33-36) | **Y** (arm_e_r6.log:103-106) |
| inherited excess as a seed in the pool (inherit.py:3-5), literal | **N** (arm_e_r6.log:67) | **N**, missed by 1.25 to 1.85 x e_need (arm_e_r6.log:107-110) |

I would score R6 = N at measured rates (reason in "Letters"). The crossing multiplier is a hypothetical network, not a scorecard row.

## What ran

| file | role |
|---|---|
| `arm_e_r6_predictions.txt` | P1-P12 and the decision rule, written and hashed before the first results run (sha256 `235e93aa8a0e4e93bd8d68cecd97befa58eae407b22873ddb90d0f616ca3b9ab`, Wed Oct 21 00:58:54 PDT 2026; unchanged; asserted in the script and printed at arm_e_r6.log:1) |
| `arm_e_r6.py` | the main run (114 s, arm_e_r6.log:204): `UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with matplotlib python arm_e_r6.py > arm_e_r6.log` |
| `arm_e_r6.log`, `arm_e_r6_results.json` | output and records (sha256 `b2707fa4...`, `2249a5b8...`) |
| `arm_e_r6_analytic.py`, `arm_e_r6_analytic_predictions.txt` | independent check (P10): full 5-species chemical-Langevin calculation, shares no formula with arm_e_r6.py. The predictions file holds 1465 numbers at 31 optimum points, written by `--predict` before the comparison (sha256 `efd13dc79a9342445fadcc98ccfac6a7cbd8058a0873a696e7952d23e301edad`, 01:21:39 PDT) |
| `arm_e_r6_analytic.log` | the P10 comparison against `arm_e_r6_results.json` (P10 FAILED, line 91) |
| `arm_e_r6_fix.py`, `arm_e_r6_fix.log`, `arm_e_r6_fix_results.json`, `arm_e_r6_analytic_fix.log` | POST HOC fix of the defect P10 found (84 s, arm_e_r6_fix.log:221), and the second comparison (max relative deviation 1.110e-10 over 1968 numbers, arm_e_r6_analytic_fix.log:91) |
| `arm_e_r6_langevin.py`, `arm_e_r6_langevin_predictions.txt`, `arm_e_r6_langevin.log` | P11 chemical-Langevin validation and the addendum V2b' (285 s on 3 cores, arm_e_r6_langevin.log:21). The addendum file is sha256 `4f4b5e3051eaa12d1b64a7959ad2bb083fe337145e705bf99b21442ebdaafaa7`, 01:37:55 PDT |

Every number below is quoted as `arm_e_r6.log:LINE` or from the logs named in the table. Total compute about 12 min wall (114 s + 129 s smoke 3 + 84 s + 285 s, plus the short smokes), at most 3 cores.

Disclosures.
1. Before the predictions file I read arm_e_rates.log:24-36 and re-ran arm_e_rates once at 5 mM (arm_e_r6_predictions.txt:6-10). No gain, noise, e_min, F_c (pond, lake, ocean), lambda(F_b) or Langevin number from the new scripts existed then. Arm G's failures (rna_seed_scan.log:105-113) were used to write the ranges.
2. Three debug runs of `arm_e_r6.py` came before the final one (docstring lines 28-36). Smoke 1 and 2 crashed in the self-check (a shape error in the Lyapunov helper, a tolerance). Smoke 3 was a full run whose output I read; it already scored P4-P8 as FAILED. After it I changed only the integrator in self-check (iii), the nan handling of the retention line, extra reporting (unstable-only rows, monomer-only-seed letters, ee_sat >= 1% in the strict (b) test, 5 digits in the s table) and added the post hoc section (5b). No pre-registered definition or P rule changed, and failed predictions stay failed.
3. Section (5b), the f grid extended down to 1e-19 /s, is POST HOC and scores no prediction (arm_e_r6.log:156-174).
4. The label "existing-rule letter = feed" in the predictions file is my choice. inherit.py:3-5 defines the inherited excess as a seed in the pool, so both are reported in every table.
5. Analytic script: one `--predict` attempt stopped at its input check because I had typed the CODATA Avogadro number (ocean.py uses 6.022e23); no number was computed. The first comparison run crashed while formatting `None` in the extras lines, after the table printed; I fixed only the print and added a `--json` option (post hoc). P10 is scored on that first table.
6. `arm_e_r6_fix.py` was written after seeing the P10 mismatch. The first run's files are kept unchanged as the record.
7. Langevin: two smoke runs (`--smoke`, M = 200 per population) came before the final run, both after P11 was frozen; they are not read as verdicts (arm_e_r6_langevin.py docstring). Smoke 1 showed that at N = 3e6 the pond (b) point is saturated at T = 3/lambda (median |L-D|/(2x) 0.107, sd(w) 0.33 of linear). I did not change the pre-registered V2b. I froze the addendum V2b' (N = 1e12) in its own file, then ran smoke 2 on its code path, then the final run. A comment in smoke 1 that guessed ee_sat ~ 0.03 was replaced by the logged 0.129 (arm_e_r6.log:104).
8. A hand-arithmetic number appears once, marked below (Delta >= 1 for the (b) seed). It is not from a run.
9. No paper was downloaded. Deng, Yu and Blackmond 2024 is closed access, so its ligation rates were not obtained and kh, R, K, f are scanned on arm_e_rates.GRID. Yu 2024 SI rates enter through arm_e_rates.py.

## Definitions (fixed before running)

* Rule: scorecard.py:9 "follows an inherited excess >= e_need (inherit.py)". Comment: scorecard.py:40 "R6 here = shown to work from racemic or <= ~1% ee".
* Settings and e_need, recomputed with `inherit.e_need` (arm_e_r6.log:4-9): pond c = 1e-2 M, tau = 100 yr, N = 1.892e25, e_need 5.770e-12; lake 1e-4 M, 1e4 yr, N = 6.022e29, 3.234e-14; ocean 1e-6 M, 1e8 yr, N = 8.021e38, 2.802e-18. e_need sqrt(N) = 25.09, 25.09, 79.36. P12 held (arm_e_r6.log:9).
* Model: Higgs and Blackmond 2025 destructive mechanism as `arm_e_rates.py`, Yu 2024 measured rates, catalyst makes the opposite hand going forward, s = k2b/k1b = 6 (4 and 9.1 as sensitivity). State (x, a, A, y, b): monomer sum and excess, achiral sink, dimer sum and excess. LD heterodimers and A are sinks. The odd sector (a, b) has Jacobian J; the margin m = J12 J21 /(J11 J22) > 1 means the racemic state is lost.
* Case (a): F = 1 (measured). Case (b): F_b = 1.5 F_c, F_c the smallest rate multiplier with lambda_max > 0 on the grid.
* Inherited excess, two variants. Feed: feed L:D = (1+g):(1-g), response read at tau. Seed: both pools start with ee0 (best case), read at tau; monomer-only seed as sensitivity.
* Followed (literal): Delta = mean/sd of the readout at tau >= 2, with counting noise from the linear-noise (chemical Langevin) equation, N = c Omega. Followed (strict): literal AND readout excess >= 1% (case b: also the linear growth reaches 1% within tau and the nonlinear saturated ee >= 1%).
* e_min = the smallest g or ee0 followed, minimised over the grid. Subsets: U = whole guarded grid; E = margin >= half the maximum. Guard for stable points: f tau >= 5 and -lambda_slow tau >= 5 (sensitivity 1, 2, 3, 10).
* Letter: Y = e_min <= e_need in all three settings, P = some, N = none (scorecard.py:9).

## Results

### Stability and crossing (arm_e_r6.log:11-17)

* At F = 1 lambda_max <= 0 everywhere on the grid in all three settings (-9.99e-14, -1.00e-13, -1.00e-13 /s). The largest margins are 2.93e-2, 2.99e-4, 9.54e-7. The racemic state is stable, so the network cannot break symmetry at measured rates (P1 held).
* F_c (s = 6): 368.7, 3.92e4, 2.998e7 (P2 held). At s = 4 the ocean does not cross (inf, arm_e_r6.log:13-15, :143). With kh allowed up to 1e4 /M/s F_c falls to 368.7, 3.687e4, 3.687e6 (arm_e_r6.log:13-15, last column). F_c x C0 = 3687, 3920, 2.998e4 mM, which is consistent with the 735 x at 5 mM quoted in README.md:476.

### (a) Measured rates, F = 1

Feed excess, literal rule (arm_e_r6.log:19-32). Subset U, best readout: e_min sqrt(N) = 1.69 (pond), 1.89 (lake), 2.0 (ocean), that is 0.067, 0.075 and 0.025 x e_need. Letter **Y**. Subset E (the network as engaged as the grid allows): 2.09, 2.15, 2.75 (0.083, 0.086, 0.035 x e_need), letter Y (arm_e_r6.log:46-49).

Why this Y is vacuous.
* The ocean optimum and the pond monomer optimum sit at f = 1e-3 /s, kh = 1e-9, where the catalysed share is 1.7e-17 to 1.7e-9 and the margin 6.2e-20 to 6.2e-12 (arm_e_r6.log:22-24, :32): the network does nothing. e_min sqrt(N) = 2 is just the counting limit of a pool of N molecules (hand case in the self-check, arm_e_r6.log:2).
* Even the engaged points pass the excess through: the gain, readout excess per feed excess, is at most 2.449, 2.084, 1.188 and at the max-margin point 1.0286, 1.0003, 0.9699 (arm_e_r6.log:55-57). A stable linear system multiplies a feed excess by about one; it does not choose a hand. The same happened for arm G (Y in all 9 cells, rna_seed_scan.log:105-113).
* Delta >= 2 is therefore a signal-to-noise test, not a test of selection. The racemic state stays racemic.

Feed excess, strict rule (readout >= 1%). The smallest feed excess that produces a 1% output is e_min = 2.042e-3, 2.400e-3, 4.207e-3 (arm_e_r6.log:34-36), against e_need 5.8e-12, 3.2e-14, 2.8e-18: 3.5e8, 7.4e10, 1.5e15 times too large. Letter **N**. Subset E: 2.678e-3, 2.746e-3, 4.207e-3, N (arm_e_r6.log:51-53). The s = 4 and s = 9.1 values change this by under 0.03% (arm_e_r6.log:152-154).

Seed excess (arm_e_r6.log:62-96). A seed in a stable pool is washed out at rate f, and the pool carries it only while it holds N molecules. With the guard f tau >= 5 the smallest followed seed is e_min sqrt(N) = 1095 (pond), 1096 (lake), 4.6e137 (ocean); that is 43.6, 43.7 and 5.8e135 x e_need (arm_e_r6.log:67). Letter **N**.
* The ocean number is the f grid floor of 1e-13 /s: the pool needs f tau >= 5, and at tau = 1e8 yr that is f >= 1.6e-15 /s, below the floor (the guard-5 optimum has f tau = 316, arm_e_r6.log:94). With the floor extended to 1e-19 /s (POST HOC, arm_e_r6.log:156-174) the ocean is 1095, 13.8 x e_need (arm_e_r6.log:160-161). Still N.
* The letter depends on the guard. At guard 1 (f tau >= 1: the seed has already lost 63%) e_min = 23.6 and 23.61, 0.94 x e_need (arm_e_r6.log:64), so P on the pre-registered grid and Y with the extended floor (0.9405, 0.9407, 0.2973 x e_need, arm_e_r6.log:158-159). At guard 2 and 3 it is 3.74 x e_need (arm_e_r6.log:65-66), N.
* Strict seed: e_min = 2.32 and 2.74 in the pond and lake (arm_e_r6.log:91-92). A seed above 100% ee does not exist, so N.
* Monomer-only seed, retained fraction of ee0 at tau at engaged guard-passing points: 4.98e-8 in the pond; lake and ocean have no engaged point that passes the guard (arm_e_r6.log:95).

Implied letters for (a): literal feed **Y** (U and E), strict feed **N**, seed **N** (literal and strict, guard >= 2).

### (b) At the crossing multiplier, F = 1.5 F_c (369, 3.92e4 and 3.0e7 x measured)

* Feed, literal: e_min sqrt(N) = 0.19, 0.046, 0.479 (0.0076, 0.0018, 0.0060 x e_need), letter Y (arm_e_r6.log:99-102). The lake optimum is a stable point (lambda = -1.08e-10 /s, arm_e_r6.log:101). Restricted to points with lambda > 0 (arm_e_r6.log:116): 0.19, 0.267, 0.479.
* Feed, strict: 0.566, 0.267, 0.479 (0.023, 0.011, 0.006 x e_need), saturated ee_sat 0.129, 0.076, 0.028 (arm_e_r6.log:104-106). Letter Y.
* Seed, literal: 31.4, 44.2, 147, that is **1.25, 1.76 and 1.85 x e_need** (arm_e_r6.log:108-110), saturated ee 0.213, 0.296, 0.259. Letter **N**. Strict seed: 43.6, 44.2, 147 (1.74, 1.76, 1.85), N (arm_e_r6.log:112-114). Monomer-only seed: 42.2, 33.2, 10.6 x e_need (arm_e_r6.log:120-121).
* Hand arithmetic, not a run: the literal e_min is linear in the Delta threshold, so Delta >= 1 would give 0.63, 0.88, 0.93 x e_need (from 1.25/2, 1.76/2, 1.85/2) and a Y. The existing rule says Delta >= 2.
* Sensitivity (arm_e_r6.log:127-149): multipliers 1.1, 3, 10 x F_c, s = 4 and 9.1, kh to 1e4. The literal feed e_min stays between 0.0009 and 0.0098 x e_need. The strict pond value fails only at 1.1 x F_c, where lambda tau = 3.2 and the 1% gate binds (2.12e5 x e_need, arm_e_r6.log:129). At s = 4 the ocean does not cross (arm_e_r6.log:143), so the feed letter there is P (pond and lake met).
* Saturation: ee_sat is 0.026 to 0.296 at every (b) point (P9 held, arm_e_r6.log:125), so the tipped network reaches an output above 1%.

Implied letters for (b): feed **Y** (literal and strict), seed **N** (narrowly).

### Independent checks

* P10 (analytic, 5-species channels, finite-difference Jacobian, 5x5 Lyapunov, LAPACK eig). The 1465 predicted numbers reproduce to 1e-9 (arm_e_r6_analytic.log:4). The comparison against `arm_e_r6_results.json` found a relative deviation of inf over 1900 numbers (arm_e_r6_analytic.log:57, :91): **P10 FAILED**. Cause: `arm_e_r6.seed_stable` (arm_e_r6.py:99-108) wrote the stable-case seed response with factors (J11 - l2) and (J22 - l1), which cancel to a few digits because the fast eigenvalue l2 sits next to J11. It affects only the monomer-only seed readouts: the old routine deviates by up to 3.92 from an eigendecomposition reference (arm_e_r6_fix.log:8) and in the ocean the monomer-only seed readout was lost (`None` in `arm_e_r6_results.json`; true value about 8e-157, arm_e_r6_fix.py:6), so its e_min was `None` and is now 1.2e154 (arm_e_r6_fix.log:250). The fix takes the small factor of each pair from (J11-l1)(J11-l2) = -J12 J21 (arm_e_r6_fix.py:33-50) and agrees with the reference to 2.7e-12 (arm_e_r6_fix.log:10-15). Effect on reported results: the first-run and fixed-run logs are identical but for the runtime line (arm_e_r6.log:204 and arm_e_r6_fix.log:221); the largest relative change of any record's e_min is 7.9e-8 (the ocean `seed_a`, computed from the two JSON files), no letter changes, the feed and (b) records change by at most 7.6e-12 (arm_e_r6_fix.log:224-247, :260-271), and only the monomer-only seed components move (up to 4.1e-4 in the lake, `None` to 1e154 in the ocean, arm_e_r6_fix.log:248-259). Post hoc, the corrected run against the independent code: 1.110e-10 over 1968 numbers (arm_e_r6_analytic_fix.log:91). The even sector is stable at every point (arm_e_r6_analytic.log:92).
* P11 (chemical Langevin simulation; pond point `a|feed_a|E`, N = 3e6, M = 20000 per population; arm_e_r6_langevin.log:4-8). Simulated sd of the monomer excess 1.0737e-3 against 1.0770e-3 linear (ratio 0.9969). P(correct sign) = 0.8380 (se 0.0026) at Delta = 1 (Phi = 0.8413) and 0.9772 at Delta = 2 (0.9772). V2b (unstable pond point, N = 3e6, T = 3/lambda): P(sign) = 0.8405 and 0.9818 against 0.8413 and 0.9772 (arm_e_r6_langevin.log:12-13). **P11 held** (arm_e_r6_langevin.log:20).
* Caveat on V2b. At N = 3e6 that run is saturated: the g = 0 sd of the unstable-mode projection is 0.336 of the linear value and the median |L-D|/(2x) is 0.106 (arm_e_r6_langevin.log:11). P11 V2b is scored on the sign statistic, which it passes, because the sign is fixed early while the amplitude saturates. The linear accumulated-noise formula behind the (b) e_min is tested by the addendum V2b' at N = 1e12 (arm_e_r6_langevin.log:14-18): sd ratio 1.0032, p95 of |L-D|/(2x) = 1.05e-3, P(sign) = 0.8325 and 0.9750 against 0.8413 and 0.9772. **Q1, Q2a, Q2b, Q3 held** (arm_e_r6_langevin.log:21), reported separately from P11.

## Predictions: held and failed (arm_e_r6.log:192-202, arm_e_r6_analytic.log:91, arm_e_r6_langevin.log:20-21)

Held: P1, P2, P3, P9, P11, P12 (and the addendum Q). Failed: P4, P5, P6, P7, P8, P10. The decision rule says partial counts as failed.

* P4 failed. Gain at the max-margin point was predicted in [1.0, 1.1]; the ocean value is 0.9699 (arm_e_r6.log:57, :59). The pond (1.0286), lake (1.0003) and the strict-N part were right. The network leaks a little to the achiral sink at the ocean's quiet optimum.
* P5 failed. Subset E monomer xN was predicted in [8, 400] (central 50) and came out 3.68, 3.72, 2.75 (arm_e_r6.log:60). The dimer range [1, 25] held (2.09, 2.15, 5.2) and so did the best-readout letter Y. I took the monomer pool from the 5 mM max-margin point (x = 1.6e-6 M, arm_e_r6_predictions.txt:7-8); the E-subset optimum has x = 1.44e-3 M at the pond (arm_e_r6_langevin.log:5), so its monomer counting noise is far smaller than I assumed.
* P6 failed. Pond and lake guard-5 xN (1095, 1096) are in [100, 5000] and the letter is N x3 as predicted. The ocean value (4.6e137) is out of range because of the f floor, and so is its guard-1 value (23.6, 23.6, 4.6e137 against [3, 40]); the monomer-only retention is `nan` for lake and ocean (arm_e_r6.log:96).
* P7 failed. Literal feed xN was predicted in [0.5, 100]; observed 0.19, 0.046, 0.479 (arm_e_r6.log:123). The letter Y held. The range was written without a calculation for the unstable case. I did not diagnose why the (b) feed value sits below 0.5; all sensitivity rows fall the same way (xN 0.05 to 0.5, arm_e_r6.log:129-149), so it is not one grid point.
* P8 failed. I expected the pond optimum to have lambda tau below 23 so that the 1% gate would give N. Observed lambda tau = 32.8 and the strict (b) feed letter is Y (arm_e_r6.log:124). The expectation does hold at 1.1 x F_c (lambda tau = 3.2, arm_e_r6.log:129).
* P10 failed on the pre-registered 1e-5 rule through the defect above. The corrected numbers agree (1.1e-10) and no reported e_min or letter moved.
* Letter predictions: held in P3, P4 (strict N), P5, P6 (seed N) and P7 (feed Y); wrong in P8 (strict (b) feed predicted N or P, observed Y).

## Letters

| reading | rule | pond / lake / ocean | letter | source |
|---|---|---|---|---|
| (a) F = 1, feed excess | literal, U | 0.067 / 0.075 / 0.025 x e_need | Y | arm_e_r6.log:29-32 |
| (a) F = 1, feed excess | literal, E | 0.083 / 0.086 / 0.035 | Y | arm_e_r6.log:46-49 |
| (a) F = 1, feed excess | strict (>= 1% out) | 3.5e8 / 7.4e10 / 1.5e15 | N | arm_e_r6.log:33-36 |
| (a) F = 1, seed | literal, guard 5 | 43.6 / 43.7 / 5.8e135 (13.8 post hoc) | N | arm_e_r6.log:67, :160-161 |
| (a) F = 1, seed | literal, guard 1 | 0.94 / 0.94 / 5.8e135 (0.30 post hoc) | P (Y post hoc) | arm_e_r6.log:64, :158-159 |
| (a) F = 1, seed | strict | e_min > 1 | N | arm_e_r6.log:91-93 |
| (b) 1.5 F_c, feed | literal / strict | 0.0076 / 0.0018 / 0.0060; 0.023 / 0.011 / 0.006 | Y | arm_e_r6.log:99-106 |
| (b) 1.5 F_c, seed | literal | 1.25 / 1.76 / 1.85 | N | arm_e_r6.log:108-110 |
| (b) 1.5 F_c, seed | strict | 1.74 / 1.76 / 1.85 | N | arm_e_r6.log:112-114 |

Reading for the scorecard. Under the existing rule as coded, with a feed excess and Delta >= 2, (a) gives Y and (b) gives Y. I would not put that Y in the cell: the same test passes a network that does nothing (cat share 1.7e-17 to 1.7e-9, arm_e_r6.log:22-32). inherit.py:3-5 defines the inherited excess as a seed, and the scorecard comment at scorecard.py:40 asks for work from racemic or <= ~1% ee. Both give N at measured rates. A Y here would also add 1.0 to the row score (5.0 to 6.0 of 9; `?` scores 0 in scorecard.py:71, so N leaves it at 5.0), hand arithmetic from the `_m` string at scorecard.py:60.

## Proposed wording (not applied)

**scorecard.py:60.** Set R6 to `N`: `_m("P P Y Y N N Y Y N")`. Append to the comment: `R6 = N at measured rates (arm_e_r6.log:33-36, :67): the racemic state is stable, so a feed excess is passed through (gain <= 2.4, Delta >= 2 holds only as a counting limit, arm_e_r6.log:55-57), 1% output needs a 0.2-0.4% feed excess (5.8e-12 needed), and a seed in the pool is followed down to 44 x e_need (pond, lake). At 1.5 x the crossing multiplier (369, 3.9e4, 3.0e7 x) a feed excess is followed (0.002-0.02 x e_need) and a seed is not (1.25-1.85 x e_need), arm_e_r6.log:99-114.` If the scorecard owner prefers the literal reading, the cell is Y and the comment should say it is the counting limit.

**README.md, Part 6 transamination bullet, after "R5 = N and R9 = N at measured rates" (README.md:476).** "R6 is N at measured rates: the stable network passes an inherited excess through (gain at most 2.4), so the smallest feed excess that gives a 1% output is 0.2 to 0.4% against the 5.8e-12 needed in the pond, and a seed in the pool is followed only down to 44 times the pond's e_need (arm_e_r6.log:34-36, :67). At 369, 3.9e4 and 3.0e7 times the measured catalysed rate (1.5 x the crossing multiplier; s = 6) the network tips over and follows a feed excess far below e_need, but a seed excess only down to 1.25 to 1.85 times e_need (arm_e_r6.log:99-114)."

**journals/ija_draft.md, section 5 arm E sentence (line 822, "as is arm E").** "Arm E is also seed-following only, and the following is weak: at the measured rates a feed excess is passed through with gain at most 2.4, and a seed in the pool is followed down to about 44 times e_need (pond, lake), so its R6 is N; at 1.5 times the multiplier that makes the racemic state unstable the feed excess is followed and a seed is missed by 1.25 to 1.85 times (Section 5, `arm_e_r6.py`)." Item 105 in the change log (line 206) would add the new script and the three settings' numbers.

## Found outside the unit (not fixed)

1. `scorecard.py:9` against `scorecard.py:40`: R6 reads "follows an inherited excess >= e_need" in the header and "shown to work from racemic or <= ~1% ee" in the comment. For this row the first gives Y (arm_e_r6.log:177) and the second N (arm_e_r6.log:179).
2. `scorecard.py:9` / `inherit.py`: Delta >= 2 is vacuous for a stable linear pass-through. Evidence: arm_e_r6.log:22-32 (Y with catalysed share 1.7e-9 and margin 6e-12) and rna_seed_scan.log:105-113 (Y in all 9 cells).
3. `arm_e_rates.py:153`: the f grid floor of 1e-13 /s is unusable at ocean scale. The guard-5 seed optimum needs f tau >= 5 at tau = 1e8 yr (f >= 1.6e-15 /s). arm_e_r6.log:67 gives 4.6e137 against 1095 with the floor extended (arm_e_r6.log:160).
4. `README.md:476, :528` and `journals/ija_draft.md:822` do not separate a feed excess from a seed, and do not mention R6 for arm E.
5. `inherit.py:3-5` defines the inherited excess as a seed in the pool, while the rate-sweep sections treat it as a persistent feed bias. The two differ by 1.25 to 1.85 x e_need in (b) alone (arm_e_r6.log:99-110). The scorecard should say which one R6 means.
