# droplet_range: where the D2 droplet rule is blind, and how wide the volume range must be (2026-10-21)

Model study on SIMULATED drops, not data. The D2 rule of `journals/droplet_prereg.md` Section 2 is frozen and was not changed. Every amendment below is a proposal, marked "not applied". No existing file was edited; no commit.

Files (all new, repo root): `droplet_range_common.py` (engine, assert self-check), `droplet_range_a.py` / `.log` (part a), `droplet_range_b.py` / `.log` (part b), `droplet_range_predictions.txt` (v1) and `droplet_range_predictions_v2.txt` (v2), post hoc: `droplet_range_b_pool.py` / `.log`, `droplet_range_b_ns.py` / `.log`, `droplet_range_cap.py` / `.log`. Cites below are `file:line` of those logs.

## How it was run, and what to discount

- Predictions P1 to P24 were written first. v1 (sha256 4ddcacbe..., 00:39:04 PDT) was never edited. v2 (sha256 e81b8eb4e8342e99bc0d8c4253c6b1f233290d6005fe1e2f90f19ec10cbdf9dd, 01:11:44 PDT) has the same predictions and decision rules, with lower run counts only, because the load average was 40 to 69. Both result scripts assert the v2 hash. v2 added three [clarified] notes (P18, P22, P23) that fix undefined details; none changes a threshold.
- Smoke runs before the results run (none used for a claim; scripts are in the scratchpad, not the repo; times are file mtimes):
  - `time1.py` (00:31; timing, no verdict shares; named in both predictions files).
  - `time2.py` (00:47; 400 runs per rule set at placeholder c*; its verdict counts showed B STEP 9 of 200 and AM1 198, AM2 173 at h = 0.3, s = 0.1 with c* = -3.0 for AM). It ran after v1 (00:39) and before v2 (01:11). The `droplet_range_a.py` docstring names it; v2's own disclosure line does not, so say it here.
  - `mkcache.py` (00:52; builds the null and exclusion truth curves; its first launch lacked a `__main__` guard and was killed within a minute, no verdicts).
  - `bench_fit.py` (00:58; speed of the plateau fit loop on random binomial data, no verdicts; not named in either docstring, so named here).
- Ordering therefore: the text of P1 to P24 was fixed at 00:39 (v1), before any of these printed a verdict. v2 (01:11) changed only run counts and the three [clarified] notes (I diffed the two files: the thresholds are identical), but it was written after the `time2.py` counts were seen. A reader who wants predictions strictly blind to those counts should read v1.
- The first launch of `droplet_range_b.py` aborted in its own self-check (tolerance 1e-9, but the volume grid is not symmetric about 0). No simulation ran. Tolerance set to 1e-5, rerun.
- After seeing the part (b) log I added three POST HOC scripts (pooling over the grid, all five verdict shares, a binomial tail). They are not tests and no prediction depends on them.
- Wall time 457 s (part a) and 193 s (part b), 4 processes.

## The rule, in one paragraph

Per run, in order: NOT SCORABLE (fewer than 10 mixed or 10 single drops, no level with P_single at least 0.6 or none at most 0.4, or more than 10% of the 128 control drops read wrong), FLAT (best curve changes by under 0.4), RISING, then STEP if L < c* = -2.21, else SMOOTH. L = best COUNT log-likelihood minus best STEP log-likelihood, both with the misread rate measured on the controls. STEP means the counting picture FAILED, SMOOTH means HELD. "Rule A" below is the earlier fixed-3% rule (c* = -1.61), "rule B" is the rule now in the prereg.

## Part (a): the plateau blind spot

### What the rule reports when a sharp step lands on a plateau

Truth: P_single = h + (1 - h) Phi(-ln(V / 3.1 uL) / s). h = 0 is a pure step to zero. Main condition: pipette CV 0.05, true misread 0.03, 600 runs per cell, rule B. STEP share (`droplet_range_a.log:64-91`):

| plateau h | s = 0.05 | s = 0.1 | s = 0.15 | s = 0.3 | what the other runs say (s = 0.1) |
|---|---|---|---|---|---|
| 0 | 1.00 | 1.00 | 0.99 | 0.77 | nothing, all STEP |
| 0.1 | 0.98 | 0.96 | 0.93 | 0.48 | SMOOTH 0.04 |
| 0.2 | 0.73 | 0.73 | 0.59 | 0.20 | SMOOTH 0.27 |
| 0.25 | 0.36 | 0.30 | 0.23 | 0.04 | SMOOTH 0.69 |
| 0.3 | 0.08 | 0.05 | 0.04 | 0.01 | SMOOTH 0.93 |
| 0.4 | 0.00 | 0.00 | 0.00 | 0.00 | SMOOTH 0.60, FLAT 0.39 |
| 0.5 | 0.00 | 0.00 | 0.00 | 0.00 | FLAT 0.77, NOT SCORABLE 0.13, SMOOTH 0.09 |

- The cliff is between h = 0.1 and h = 0.3, and it is not an edge case of s: at s = 0.1 the STEP share is 0.96, 0.73, 0.30, 0.05 for h = 0.1, 0.2, 0.25, 0.3. The prereg blind-spot line (`droplet_prereg.md:242-243`, "plateau at 0.3 ... 5%") is replicated by rule A (`droplet_range_a.log:53`, 0.05; P1 held, 0.053). It understates the range: the loss already starts at h = 0.2.
- For 0.25 <= h <= 0.4 the rule reports SMOOTH, which is HELD. A sharp step to a plateau is therefore scored as support for counting. This is a false HELD, not a NOT DECIDED. At h >= 0.4 the verdicts turn into FLAT or NOT SCORABLE (not decided), so those are safe in direction.
- Other conditions, rule B, h = 0.2, s = 0.1 (`droplet_range_a.log:153,161,165`): 0.82 with no misreads, 0.50 at 6% misreads, 0.56 at pipette CV 0.15.
- The pre-registered test of "does it matter" (decision rule a): B STEP(h = 0.2, s = 0.1) = 0.728, below 0.80, so yes (`droplet_range_a.log:247`).

### Which physical situations give a plateau, and what the rule says

By construction, not from measured NaClO3 data; I read no paper for this unit, so I have no source for how large h or c really are in NaClO3 drops (not found, not searched).

| situation | shape | simulated as | rule B reports | AM2 reports (proposal) |
|---|---|---|---|---|
| A seeded or dominant-site share f of drops go single whatever the volume (seed crystal, dust, a favoured site) | floor h = f | plateau map | f <= 0.1: STEP 0.96 to 0.98; f = 0.2: 0.73; f = 0.25: 0.30; f >= 0.3: SMOOTH (HELD) 0.93 (`droplet_range_a.log:68-91`) | 0.99, 0.98, 0.96, 0.94 at f = 0.1, 0.2, 0.25, 0.3 (`:127,131,135,139`) |
| Number of nuclei saturates at N in large drops (each of random hand) | floor 2^(1-N): N=4 gives 0.125, N=3 gives 0.25, N=2 gives 0.5 | not run exactly; sits between map rows | N=4: between 0.96 and 0.73 STEP; N=3: about 0.30 STEP, rest SMOOTH; N=2: no STEP, FLAT or NOT SCORABLE | N=3: about 0.96 to 0.90; N=2: 0.27 (`:135-136,147`) |
| A second crystal in a mixed drop is missed with probability m (one-way) | floor m | scenarios | m = 0.06: STEP 1.00; m = 0.12: 0.96; m = 0.20: STEP 0.48, NOT SCORABLE 0.51 (`:196-198`) | same (`:222-224`) |
| Unreadable mixed drops dropped, share u (informative exclusion) | stretches the tail | scenarios | u = 0.5: 0.98; u = 0.8: 0.92; u = 0.9: 0.63 with NOT SCORABLE 0.17 (`:199-201`) | 0.98, 0.92, 0.61 (`:225-227`) |
| Small drops already carry two crystals (polycrystalline, two nuclei) | CEILING c < 1: P = c Phi | scenarios | c = 0.9: STEP 0.51 / 0.38 (s = 0.1 / 0.15); c = 0.8: SMOOTH 0.95 / 0.97; c = 0.7: SMOOTH 0.70, FLAT 0.30 (`:190-195`) | no better: 0.50 / 0.37; c = 0.8 same (`:216-221`) |
| First crystal seeds same-hand daughters (secondary nucleation) | counting curve with a floor, NOT a step | 27 counting-with-floor nulls | false STEP up to 0.100 for B at zero true misread (`:25`); 0.028 at 3% (`:16`) | up to 0.040 (`:31`) |

- The misread term absorbs a one-way detection miss up to m = 0.12; beyond that the control cap turns the run into NOT SCORABLE (a safe outcome). Informative exclusion is safe up to u = 0.8.
- A ceiling is the unfixed part. Neither B nor AM2 sees it. At c = 0.8 the rule says SMOOTH in 95 to 97% of runs.
- I did not simulate the seeding guards already in the plan (`droplet_prereg.md:135` lopsided-slide exclusion, `:51` D3 chance hand). They catch a seeded share only if it is lopsided in hand; a share of seeds of both hands would pass them.

### False STEP on the nulls (the cost of an amendment)

Cell groups: 87 counting and drying curves (15 counting forms from `droplet_sim2.truths()`, 72 drying models from `droplet_sim.py`), 27 exclusion curves from `exclusion_bak.py`, 27 counting-with-floor curves (3 forms x a = 2/3, 1, 5/3 x h = 0.1, 0.2, 0.3). Evaluation seed, 400 runs per cell and condition; limit 0.07. Highest STEP share over each group (`droplet_range_a.log:23-31`):

| rule | counting + drying | exclusion | counting with floor |
|---|---|---|---|
| B (frozen, c* = -2.21) | 0.045 | 0.000 | 0.100 (at 0 misreads; 0.085 at CV 0.15 `:19`) |
| AM1 (floor in STEP only, c* = -9.66) | 0.007 | 0.000 | 0.058 |
| AM2 (floor in both, c* = -2.40) | 0.045 | 0.005 | 0.040 |

Per condition AM2 stays at or under 0.045 on every group (`droplet_range_a.log:15-22`). Calibration: c*_B recomputed = -2.31 (`:6`, frozen value -2.21; P10 held), c*_AM2 = -2.40 (`:10`), c*_AM1 = -9.66 (`:9`).

### The amendment I would propose, and the one I would not

AM1 (a floor in the STEP family only) fails. Its c* has to drop to -9.66 to hold the false-STEP rate, and it then breaks pure steps: h = 0 at s = 0.1 / 0.15 / 0.3 gives STEP 0.94 / 0.79 / 0.04 (`droplet_range_a.log:94-96`), and exclusion u = 0.5 gives 0.48 (`:212`). P14 and P15 failed for it (`:248`).

AM2 (a floor in both families) keeps everything B had and covers the blind spot. STEP share, main condition, s = 0.05 / 0.1 / 0.15 / 0.3 (`droplet_range_a.log:122-149`):

| h | B | AM2 |
|---|---|---|
| 0 | 1.00 / 1.00 / 0.99 / 0.77 | 1.00 / 1.00 / 0.99 / 0.77 |
| 0.1 | 0.98 / 0.96 / 0.93 / 0.48 | 1.00 / 0.99 / 0.98 / 0.57 |
| 0.2 | 0.73 / 0.73 / 0.59 / 0.20 | 0.99 / 0.98 / 0.93 / 0.44 |
| 0.3 | 0.08 / 0.05 / 0.04 / 0.01 | 0.97 / 0.94 / 0.81 / 0.25 |
| 0.4 | 0 / 0 / 0 / 0 | 0.86 / 0.73 / 0.56 / 0.12 |
| 0.5 | 0 / 0 / 0 / 0 | 0.39 / 0.27 / 0.18 / 0.03 |

At 6% misreads, h = 0.2, s = 0.1: B 0.50, AM2 0.90 (`:161,163`). Decision rule (a): the blind spot matters (0.728 < 0.80) and AM2 meets P12, P14, P15 (`:247-250`), so AM2 is the proposal. P16 held: AM2 beats AM1 by 0.768 vs 0.042 at h = 0, s = 0.3, so the smaller change (AM1) is not enough.

Costs to state if this is adopted: AM2 does not help a ceiling (above), it does not make s = 0.3 steps easier (0.77, same as B), and at h = 0.5 it is still mostly undecided (FLAT 0.47, NOT SCORABLE 0.13, `:147`). It adds fit freedom, so c* must be re-set (-2.40) and the result tables in Section 10 re-run.

**Proposed text, NOT APPLIED (Section 2, D2, replacing the COUNT and STEP bullets):**

> - **COUNT:** P = h + (1 - h) F(mu), mu = mu0 (V/V0_pilot)^a with a in [0.05, 5/3], F one of the three forms above; h = 0 or h in {0.05, 0.10, ..., 0.40, 0.50}. h is a floor: a share of droplets that are single whatever their size (a seeded or dominant site, or a cap on the number of nuclei).
> - **STEP:** P = h + (1 - h) Phi(-ln(V/Vc)/s), s in [0.02, 0.6], h as above.
>
> Each family is fitted with h = 0 on the existing grids and with each listed h on a coarser grid (0.1 in ln mu, 0.04 in ln Vc); the best likelihood in each family is used. STEP if L < c* = -2.40. c* is the lowest 5th percentile of L over the 87 counting and drying curves and 27 counting-with-floor curves, 128 controls, misreads 0 to 10%, pipette CV 5 or 15%, set on a calibration seed (`droplet_range_a.log:10`).

**Limit sentence, NOT APPLIED (Section 10, to add whether or not AM2 is adopted):**

> Blind spots of the rule in Section 2 (simulated, `droplet_range_a.log`): a sharp step down to a plateau (a share h of droplets single whatever the size) is called STEP in 96%, 73%, 30% and 5% of runs at h = 0.1, 0.2, 0.25 and 0.3 (s = 0.1, pipette CV 5%, 3% misreads) and SMOOTH in most of the rest. A step that stops at P_single >= 0.25 therefore reads as counting. A ceiling, in which small droplets are already mixed in 10 to 20% of cases, is called SMOOTH in 49% to 62% of runs at c = 0.9 and in 95% to 97% at c = 0.8. With a plateau in both fits the h = 0.1 to 0.3 cases are called STEP in 94% to 99% (s = 0.1), and the ceiling is still not seen.

## Part (b): how wide the volume range must be

Setup: u = V*/V0_pilot, where V* is the critical volume of a sharp step. A step at V* gives a pilot fit V0 of about V*/ln 2, so u near 0.69 is the centre. The plan's own allowance (true V0 within 0.5 to 2 x pilot) is u from 0.35 to 1.39 (span S1). I also tried 0.2 to 4 x (u from 0.14 to 2.8, span S2). Rule B, c* = -2.21, 128 controls, 300 runs per cell. Designs (nominal drops = levels x 25 + 128 controls):

| design | levels | extension rule | nominal drops |
|---|---|---|---|
| base (the plan) | 8, 0.2 to 1.6 x V0 | up to 2 extra levels | 328 |
| ext4 | same | up to 4 extra levels | 328 (up to 428) |
| K10 | 10, 0.2 to 2.9 | up to 2 | 378 |
| K12 | 12, 0.11 to 2.9 | up to 2 | 428 |

### Where a step can sit

Base design, STEP share, s = 0.1, true misread 0.03 (`droplet_range_b.log:6`):

| position of the step | u | STEP |
|---|---|---|
| below the range | 0.05, 0.07, 0.1 | 0.00 (all runs end NOT SCORABLE) |
| low end | 0.14, 0.2, 0.28 | 1.00, 0.98, 1.00 |
| middle | 0.35 to 1.6 | 0.99 to 1.00 |
| top | 1.9, 2.26 | 1.00, 1.00 |
| just past the top | 2.7 | 0.95 (s = 0.15: 0.79, `:9`) |
| above the range | 3.2, 4.5, 6.4 | 0.00 |

- Edges at 95% (`:32,35`): s = 0.1: u from 0.14 to 2.7; s = 0.15: 0.14 to 2.26. This is wider than the nominal 0.2 to 1.6 because the extension rule adds up to two levels, to about 0.11 and 2.9 x V0. Past those the run has no level with P_single on the other side, and the rule returns NOT SCORABLE. P19 and P20 held.
- At misread 0.03, ext4 pushes the upper edge to 4.5 and the lower to 0.07 at s = 0.1 (0.10 at s = 0.15) (`:37-42`); at misread 0.06 none of the three reaches a top edge beyond 0.8 on the 0.95 criterion (`:39,42,45,48,51,54`). K10 reaches 4.5 on top (`:43-48`), K12 reaches 4.5 and 0.07 (`:49-54`). The grid stops at 6.4, which no design reaches.
- The known issue: without the misread term in the fits (E0) a sharp step at u = 1.39 is read as broad (SMOOTH) in 92% of runs at 3% misreads, STEP share 0.069; with the measured term (rule B) 0.994; rule A 1.000 (`droplet_range_b.log:62-64`; 2000 runs; P18 held). Along the grid E0 falls from 1.00 at u = 0.4 to 0.23 at u = 1 (`:65`), so the problem grows toward the top of the range. `droplet_power.log:6-9` (0.84 vs 1.00) is a different statistic (the old h <= 3 fit), so the sizes do not compare; the direction agrees. Under rule B the top of the range is fine at 3% misreads.

### Absolute volume, pilot included

Full pipeline: pilot ladder 0.25 to 16 uL doubling, 20 drops per rung, V0 fitted, then the main run. STEP share, s = 0.1, CV 0.05, 3% misreads, 600 runs (`droplet_range_b.log:68-72`):

| V* (uL) | 0.1 | 0.2 | 0.4 | 0.7 | 1 | 2 | 3.1 | 4 | 7 | 12 | 20 | 32 | 64 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| base | 0.90 | 0.96 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.995 | 0.99 | 0.68 | 0.79 | 0.91 |
| ext4 | 0.94 | 0.96 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 1.00 | 0.90 | 0.92 | 0.98 |
| K10 | 0.86 | 0.97 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.66 | 0.77 | 0.93 |
| K12 | 0.94 | 0.96 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 1.00 | 1.00 | 0.93 | 0.93 | 0.97 |

- The base design reaches 95% for V* from 0.2 to 12 uL. Below 0.2 uL and above 12 uL it loses runs to NOT SCORABLE (V* = 0.1: 0.10; V* = 32: 0.20, `:76`); the top is limited by the pilot ladder (16 uL), not by the main run. Non-monotone at 20 to 64 uL for the same reason (the censored pilot fit lands anywhere). P22 held at the grid points it named (0.1 uL: 0.902, 64 uL: 0.908, both under 0.95).
- The pilot places u well inside S1 for V* from 0.7 to 7 uL: the central 95% of u is 0.36 to 0.82, 0.49 to 0.84, 0.62 to 0.88, 0.63 to 0.95 at V* = 0.7, 2, 3.1, 7 uL (`:91`, P21 held). So the wide span S2 only matters for a step outside the pilot's 0.4 to 12 uL window or a biased pilot.

### Lowest STEP share over the span, and what it costs

Lowest STEP share over the grid points of each span (`droplet_range_b.log:57-60`), by (s, misread):

| design | span | (0.1, 0.03) | (0.1, 0.06) | (0.15, 0.03) | (0.15, 0.06) |
|---|---|---|---|---|---|
| base | S1 | 0.99 | 0.93 | 0.97 | 0.89 |
| base | S2 | 0.95 | 0.88 | 0.79 | 0.71 |
| ext4 | S2 | 0.98 | 0.90 | 0.96 | 0.86 |
| K10 | S2 | 0.98 | 0.92 | 0.97 | 0.85 |
| K12 | S2 | 0.99 | 0.93 | 0.96 | 0.87 |

Cost, measured mean drops per run including the 128 controls, over the S2 points at s = 0.1, 3% misreads (`droplet_range_b.log:84`): base 343, ext4 344, K10 383, K12 428. False STEP at the frozen c* on the 87 nulls: base 0.043, ext4 0.037, K10 0.050, K12 0.060 (`:79-82`; P23's null limit held).

### The pre-registered decision, and what I read from it

Decision rule (b) as registered (lowest-cost design whose lowest share is at least 0.95 in the two CV 0.05 conditions at both s, null at most 0.07): no design qualifies for S1 or for S2 (`droplet_range_b.log:96-99`). P17 failed (the base design at s = 0.15, 6% misreads, lowest share 0.893 against 0.90, `:87`) and P24 failed (no design covers S2 in all four conditions, `:94`).

Two things make that outcome look worse than the range question deserves:

1. At a true misread of 6%, 4 to 5% of runs go NOT SCORABLE because more than 10% of the 128 controls are read wrong: pooled over S1, NOT SCORABLE is 0.049, 0.044, 0.041, 0.048, 0.040, 0.044, 0.042, 0.052 across the designs and s (`droplet_range_b_ns.log:3,6,9,12,15,18,21,24`). That matches the binomial tail, 0.044 (`droplet_range_cap.log:6`). Dividing it out, base S1 s = 0.1 is 0.942 / 0.951 = 0.99. So the 6% shortfall is the control cap, not the volume range, and no extra level can fix it.
2. The registered rule takes the minimum of 7 to 14 noisy shares of 300 runs, which is biased low. Pooled over the span (POST HOC, `droplet_range_b_pool.log`): at 3% misreads every design is at least 0.975 on S1 and S2 (lower 95% bound at least 0.970, `:2,5,8,11,14,17,20,23`). The one weak spot is the single point u = 2.7 for the base design at s = 0.15 (0.79, `droplet_range_b.log:9`).

Reading (not a registered result):
- **Span S1 (step within 0.5 to 2 x the pilot V0), misreads at most 3%:** the plan as written (base, 8 levels, up to 2 extra) is enough. Lowest share 0.99 and 0.97 (`:57`), 328 drops nominal, 343 mean.
- **Span S2 (0.2 to 4 x), misreads at most 3%:** allow up to 4 extra levels on the side that needs them (ext4). Lowest share 0.98 and 0.96 (`:58`); mean cost 344, no more than 100 above nominal. K10 is the fixed-cost alternative (+50 drops, lowest 0.98 and 0.97, `:59`), K12 is +100 (lowest 0.99 and 0.96, `:60`) and also reaches 0.07 on the low side. ext4 has the lowest false STEP of the four (0.037).
- **True misread 6%:** no design gets a STEP share of 0.95 on the full grid; the cap above is why. More levels do not help. More controls would change both the cap and the noise in the measured rate; not simulated here.

**Proposed text, NOT APPLIED (Section 6, extension rule):**

> Extension rule: after the main run, if no level has P_single <= 0.4, add one level at 1.35x the top volume and repeat; if none has P_single >= 0.6, add one at the bottom and repeat. Same droplet count per new level, at most four additional levels in all (at most 100 more droplets). Only P_single per level is looked at to make this decision, not h.

## Predictions: held and failed

19 held, 5 failed. Numbers are from the lines named.

| P | result | number |
|---|---|---|
| P1 | held | rule A at h = 0.3, s = 0.15: 0.053 (`droplet_range_a.log:230`) |
| P2 | held | B at h = 0: 1.000 / 0.998 / 0.993 (`:231`) |
| P3 | held | B at h = 0.3: at most 0.08 in all four conditions (`:232`) |
| P4 | held | STEP never rises with h (`:233`) |
| P5 | held | B at h = 0.1 / 0.2: 0.96 / 0.73 and 0.93 / 0.59 (`:234`) |
| P6 | held | B at h = 0.5: 0.000 (`:235`) |
| P7 | FAILED | ceiling c = 0.9: B STEP 0.513 and 0.383, predicted at least 0.80 (`:236`) |
| P8 | held | one-way miss m = 0.06 / 0.12: 1.000 / 0.960; m = 0.20 NOT SCORABLE 0.508 (`:237`) |
| P9 | FAILED | exclusion u = 0.5 / 0.8 / 0.9: 0.980 / 0.918 / 0.628; the first two pass, u = 0.9 does not (`:238`) |
| P10 | held | c*_B recomputed -2.31 (`:239`) |
| P11 | held | c*_AM2 = -2.40 (`:240`) |
| P12 | held | AM2 highest false STEP 0.045 (`:241`) |
| P13 | FAILED | B on counting-with-floor cells: 0.100, predicted at most 0.07; only at zero true misread (`:242`, `:25`) |
| P14 | held | AM2 at h = 0.1 / 0.2 / 0.3: 0.99 0.98 / 0.98 0.93 / 0.94 0.81 (`:243`) |
| P15 | held | AM2 at h = 0: 1.000 1.000 0.992; s = 0.3: 0.768 vs B 0.773 (`:244`) |
| P16 | held | AM2 minus AM1 at s = 0.3: 0.768 vs 0.042 (`:245`) |
| P17 | FAILED | base on S1: all cells pass except s = 0.15 at 6% misreads, 0.893 vs 0.90 (`droplet_range_b.log:87`) |
| P18 | held | E0 0.069 vs measured 0.994 (`:88`) |
| P19 | held | base u_lo 0.14, u_hi 2.7 (`:89`) |
| P20 | held | STEP 0.00 at u = 0.05, 0.07, 4.5, 6.4 (`:90`) |
| P21 | held | central 95% of u inside 0.35 to 1.39 at V* = 0.7, 2, 3.1, 7 (`:91`); the lowest edge, 0.36, passes by 0.01 |
| P22 | held | base at V* = 1, 2, 4, 7 uL at least 0.995; 0.1 and 64 uL under 0.95 (`:92`) |
| P23 | held | K10 u_hi 4.5 vs base 2.7; null 0.050 (`:93`) |
| P24 | FAILED | base fails S2, but ext4, K10 and K12 also fail in the 6% conditions (`:94`) |

Why P7, P9, P13 failed: I expected a ceiling and heavy informative exclusion to look like a step, and a floor to be safe at the frozen c*. The ceiling is read SMOOTH for c <= 0.8 and half the time for c = 0.9; exclusion at 90% pushes runs to NOT SCORABLE (0.17) and SMOOTH (0.21); a floor on a counting curve is called STEP up to 10% of the time only when true misreads are zero (`:25`), 2.8% at 3% misreads.

## Found, outside this unit (not fixed)

- `droplet_sim2.py:88`: `for _ in range(2)` adds at most 2 extension levels in total. The prereg says "up to 100 more" (`droplet_prereg.md:114`) and "at most two per side" (`:118`). The two agree when only one side needs levels, and differ if both do. Every power number in this unit uses the 2-level cap; ext4 is a change to both.
- `droplet_prereg.md:74-75` (steps with s <= 0.15 called STEP in 81 to 97% of runs at misreads of 8 to 10%) come from `misread_sim.log`, whose policy B has no control cap. D2 item 1 (`droplet_prereg.md:32-34`) makes a run NOT SCORABLE above 10% measured misreads. With the cap, P(NOT SCORABLE) at a true 8%, 9%, 10% is 0.224, 0.367, 0.519 (`droplet_range_cap.log:8-10`). The STEP shares in that line are therefore conditional on being scorable, not shares of runs.
- `droplet_prereg.md:205-207` ("no counting form in Section 1 and none of the drying models below" called STEP in more than 5.4%): true for those families (here 0.045, `droplet_range_a.log:23`), but a counting curve with a floor was called STEP up to 0.100 (`droplet_range_a.log:25`, zero true misread). Floors are outside the claim's list, so this is a scope gap, not a wrong number.
- `droplet_prereg.md:242-243` states the blind spot at a plateau of 0.3 only; the loss starts at 0.2 (`droplet_range_a.log:73`, STEP 0.73).
- `droplet_power.log:8` (0.84) is for the old h <= 3 fit. Rule B's counterpart with the term removed (E0) is 0.069 (`droplet_range_b.log:63`). Same direction, not the same number; do not cite one as the other.

## Not done

- Seeding guards, N-nuclei cap exactly, tilted or non-symmetric misreads, 192 controls: not simulated.
- No NaClO3 data and no paper was used. The plateau heights and ceilings are grid values, not estimates.
- Wall time: 457 s (part a) plus 193 s (part b) for the results runs, well under the 45 min budget; the one-off cache build (`mkcache.py`) was not timed in a log (not found).
