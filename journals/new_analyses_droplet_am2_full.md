# Unit droplet_am2_full: does AM2 still fix the plateau blind spot at its full-grid threshold, and what does it cost

Simulated drops, not data. Date 2026-10-21. `journals/droplet_prereg.md` was not edited (mtime 2026-10-20 23:34:26, older than the predictions file, 2026-10-21 03:34:39). `git status` shows as modified only the five files that were already modified when the session began (`README.md` 00:08, `advect2d.py` 01:47, `deephand_scale.py` 01:17, `journals/ija_draft.md` 00:09, `scorecard.py` 00:07, all older than the predictions file at 03:34; not touched by this unit); everything this unit made is new and untracked. Every wording below is **NOT APPLIED**.

**Bottom line.** AM2 at the full-grid threshold c* = -2.976 (AMf) still fixes the blind spot (STEP 0.97 / 0.94 / 0.91 at h = 0.2 / 0.25 / 0.3, against B's 0.73 / 0.30 / 0.05, main condition, s = 0.1) and its false STEP is no worse than B's (0.028 R-world, 0.042 K-world, against B's 0.045 and 0.054). It costs power on broad steps: up to 7 points (with the cap) and 10 points (without it) of STEP share in the Section 10 ranges, so decision rule D1b fails. Recommendation: keep B, add a limit sentence, offer AM2 only as a labelled post hoc sensitivity analysis. **Read section 2 before trusting any of this:** two post-output code changes were made, and the first printed recommendation ("reject AM2") was wrong because of an indexing bug of mine.

## 1. Files, hashes, run

| file | sha256 / state |
|---|---|
| `droplet_am2_full_predictions.txt` (61 lines, mtime 03:34:39, before any result) | `4a0900aea24cc052301df5d906824b25759df1b3c8b6c4ca66527880431bd9ce` (asserted in the script) |
| `droplet_am2_full.py` (612 lines, final version, last edit 03:57:54) | `7ab04a87147460b4efa1da3f995f3c0cf360fd75c7abda0f2a68a987ef2567db` |
| `droplet_am2_full.log` (317 lines, final clean run) | `78fc5f5a5153aa5721f6d8d419a8b5a52b2abfbd6aa81c3635f7ff80e7d327b7` |
| `droplet_am2_full.err` (161 lines, progress only, no Traceback; K 237 s, R-ev 335 s, R-map 352 s, done 356 s) | `f328728a0e59371ca45d8dc718b6a6cc93432fe2cf124f3fa225603748d795ca` |
| `droplet_am2_full_attempt1.log` (13 lines; stopped at the R-world gate) | `8a9c4828ec1d34f9256ca27d892c2582b32ff30016904f0f51e3572aee4e6a8c` |
| `droplet_am2_full_attempt2.log` (317 lines; the log with the wrong D2) | `71206f114f19f991a486b8da3876798b65756942f049728b07f3047651a65152` |

Run line (3 processes, niced; the self-check was run first with `--check` appended and printed "self-check passed"):

```
UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_am2_full.py > droplet_am2_full.log 2> droplet_am2_full.err
```

The script imports `droplet_ctrl_am2`, `droplet_range_a` and `droplet_range_common` and edits none of them; the extra rules (AMf, AMs, Bsub and their no-cap "n" versions) are added to `RULES` in memory only. Same fit groups and cap as AM2 / B, other c*, so one fit per simulated run feeds every threshold. Two worlds, because the earlier units used different control draws: K = `droplet_ctrl_am2.stage2` (91 truths x 12 conditions x 500 runs, seeds [5, i, j]), R = `droplet_range_a` stages 2-3 (141 cells x 8 conditions x 400 runs, seeds [11, i, j]; plateau map 28 truths x 4 conditions x 600 runs, seeds [12, i, j]). A scratchpad checkpoint directory only saved intermediate arrays.

Attempt history, in order (all attempts used the same simulated counts):
1. **Attempt 1** (script sha `d4c649336b322a5d4f6a6bc3e469d4e4caf5791d85d02c33c46b427ed39e5592`, simulated everything in 292 s, then crashed in `gate_r`: `IndexError: list index out of range`, `droplet_am2_full.py` line 232, `if l.split("|")[col].split() != l.split("|")[1].split():`). Its stdout is `droplet_am2_full_attempt1.log` (c* checks and the K gates passed); the traceback was in the `.err`, which the next run overwrote (a copy is in the scratchpad, `.../scratchpad/am2_full/attempt1.err`).
2. **Attempt 2** (script sha `f78963127c279bb68d71a5713790c7aae358df1eeb7da5c33a3c4bd30882ac65`, reusing attempt 1's simulations from checkpoints): completed, printed D2 FAILED and "reject AM2". Log kept as `droplet_am2_full_attempt2.log`. I reconstructed the attempt-1 and attempt-2 scripts by reversing the two edits below; the attempt-1 reconstruction hashes to the recorded `d4c649...`, which confirms the reversal is exact.
3. **Attempt 3** (final script, checkpoints reused): log identical to the final log.
4. **Final run** (final script, checkpoints emptied, everything simulated again in 356 s): `droplet_am2_full.log`. Its text is byte-identical to attempt 3's (`diff` empty), so every number below comes from one clean run of the final script.

## 2. POST-OUTPUT CHANGES (both made after I had seen output)

Neither changes a simulation, a seed, a threshold, a band or a decision rule. Both were made after seeing output, so neither was blind.

**Change 1, after attempt 1 crashed (gate code, `gate_r`).** My own generated false-STEP line has one `|`, so indexing it with `col = 3` (the AM2 column of the old log) raised IndexError. The first comparison was a leftover; the comparison against the log's row `tab[j]` was already there and is the real gate.

```
before:   n += 1
          if l.split("|")[col].split() != l.split("|")[1].split():
              bad.append(("false-STEP row", r, j))
          if tab[j].split("|")[col].split() != l.split("|")[1].split():
after:    n += 1
          if tab[j].split("|")[col].split() != l.split("|")[1].split():
```

**Change 2, after attempt 2 printed D2 FAILED (scoring code, K-world false-STEP figure).** `cnt[r]` has shape (91, 12, 5): the last axis is the verdict counts STEP, SMOOTH, FLAT, RISING, NS. `kcm` took the maximum over all five columns, which is the SMOOTH count (500 of 500 in some cell), so it printed 1.000. The no-cap figure on the same line (`d2[r][2]`) already used `[:87, :, 0]` and printed 0.042; that inconsistency inside one printed line is how I found the bug (`droplet_am2_full_attempt2.log:312-313`).

```
before:   kcm = {r: cnt[r][:87].max() / REPS for r in ("B", "AM2", "AMs", "AMf", "Bsub", "AMfn")}
after:    kcm = {r: cnt[r][:87, :, 0].max() / REPS for r in ("B", "AM2", "AMs", "AMf", "Bsub", "AMfn")}   # STEP column only
```

Cosmetic edits in the same save, output text only: `rep = [round(...)...]` and `rb = [round(...)...]` now wrapped in `float(...)` (the log no longer prints `np.float64(0.98)`); `per item {f2}` and `largest drops in band {f3d}` now print `{[bool(x) for x in f2]}` and `{[bool(x) for x in f3d]}` (no `np.True_`). The attempt-2 log shows the old forms (`droplet_am2_full_attempt2.log:260, 269, 307`).

**The pre-registered wording that defines the K-world false-STEP figure** (`droplet_am2_full_predictions.txt`, quoted exactly; no sentence in the file says "STEP column" in those words):
- line 31: `PREDICTIONS (percent = share of runs called STEP; "main" = CV 0.05, misread 0.03, s = 0.1)`
- line 36: `F2  Highest false STEP over the 87 counting cells, rule AM2 with the cap.`
- line 38: `K-world (12 conditions, 500 runs): -2.40 0.052 (droplet_ctrl_am2.log:157); AMs in [0.038, 0.052] (point 0.048); AMf in [0.017, 0.040] (point 0.029).` (the 0.052 anchor is a STEP share; it is the "5.2" of `droplet_ctrl_am2.log:154`, the line number in the file being 157 is one of its blemishes, section 8)
- line 56: `D2  False STEP: the AMf counting maximum is at most 0.07 in the K-world (12 conditions) and in the R-world (8 conditions).`

**D2 both ways.**

| | as originally printed (max over all five columns), `attempt2.log:312-315` | as fixed, `droplet_am2_full.log:312-315` |
|---|---|---|
| AMf (-2.9763): R-world / K cap / K no cap | 0.028 / **1.000** / 0.042 -> D2 False | 0.028 / 0.042 / 0.042 -> D2 True |
| AMs (-2.4634): R-world / K cap / K no cap | 0.035 / **1.000** / 0.050 -> D2 False | 0.035 / 0.050 / 0.050 -> D2 True |
| printed recommendation | "reject AM2 (D2 fails: false STEP above 0.07 at the full-grid threshold)" | "keep B; add a limit sentence ...; offer AM2 only as a labelled post hoc sensitivity analysis (D1a and D2 hold, D1b fails)" |
| F2 as printed | 3 of 6 items false (the three K-world items printed 1.000) | 1 of 6 items false (AMf K-world 0.042 above its band, which ends at 0.040) |

My scoring of D2 as written: **HELD**. Reason: "False STEP" in the file means the share of runs called STEP (line 31 sets that unit; line 36 and the 0.052 anchor are STEP shares), and the quantity the same line calls "counting maximum" is the STEP share (`droplet_ctrl_am2.log:154`, B 5.4). A maximum over SMOOTH counts is not a false STEP. The cap/no-cap choice is not stated in D2, but it does not matter: both give 0.042 (AMf) and 0.050 (AMs). The file never says "STEP column" in those words, though; a reader who holds that the text does not pin the column would score D2 AMBIGUOUS, and then the recommendation tree has no branch for it. I do not hold that reading, and I flag it so you can overrule me.

## 3. Reproduction gates and checks (`droplet_am2_full.log:6-17`)

| check | result |
|---|---|
| c* literals against recomputation from `am2_gate/cal_*.npy` and `a_cal.pkl` (log 8-11) | AMf -2.9762901306152343, AMs -2.463367462158203, Bsub -2.0587861532873992, AM2 of `droplet_range_a` -2.403298187255859: all equal |
| GATE K: `droplet_ctrl_am2.part2` text (B, AM2 at -2.40, Bn, AM2n) against `droplet_ctrl_am2.log:58-172` (log 12) | PASS, 114 lines, no differences (including that log's failed P7) |
| GATE K3: `droplet_ctrl_am2.stage3` text against `droplet_ctrl_am2.log:173-189`; B reproduces `bak_policy_b.log` (log 13) | PASS, 17 lines, 36 of 36 |
| GATE R: B and AM2 (-2.4033) false-STEP tables, plateau lines and step/NS rows against `droplet_range_a.log`; counts against `a_ev.pkl`, `a_map.pkl` (log 14) | PASS, 92 comparisons, no differences |
| GATE BAK: paired numbers re-averaged from `droplet_am2_gate.log:95-115` (log 15) | PASS |
| structure (log 16): STEP count falls as c* falls, STEP + SMOOTH / FLAT / RISING / NS equal across thresholds, cap STEP <= no-cap STEP | no violations in the K-world, R false-STEP runs or R plateau map |
| D0 (log 17) | results trusted |

## 4. Results

### 4.1 Plateau map (R-world, 600 runs per cell; `log:19-106`)

STEP share, truth P = h + (1 - h) Phi(-ln(V / 3.1 uL) / s), s = 0.1. Columns are the four noise conditions (CV, true misread): (0.05, 0), (0.05, 0.03) = main, (0.05, 0.06), (0.15, 0.03). Rows for B, Bsub, AM2, AMs, AMf are at `log:81-100`; the verdict-share lines of AMf and AMs for the main condition are `log:21-49` and `log:50-78` (h = 0.2 / 0.25 / 0.3 at lines 31 / 35 / 39 for AMf and 60 / 64 / 68 for AMs).

| h | rule (c*) | (.05, 0) | (.05, .03) | (.05, .06) | (.15, .03) |
|---|---|---|---|---|---|
| 0.2 | B (-2.21) | 0.82 | 0.73 | 0.50 | 0.56 |
| 0.2 | Bsub (-2.059) | 0.83 | 0.75 | 0.56 | 0.61 |
| 0.2 | AM2 (-2.403) | 1.00 | 0.98 | 0.90 | 0.91 |
| 0.2 | AMs (-2.463) | 1.00 | 0.98 | 0.90 | 0.90 |
| 0.2 | **AMf (-2.976)** | 1.00 | 0.97 | 0.88 | 0.86 |
| 0.25 | B | 0.29 | 0.30 | 0.23 | 0.21 |
| 0.25 | Bsub | 0.31 | 0.34 | 0.27 | 0.23 |
| 0.25 | AM2 | 1.00 | 0.96 | 0.85 | 0.84 |
| 0.25 | AMs | 1.00 | 0.96 | 0.85 | 0.83 |
| 0.25 | **AMf** | 1.00 | 0.94 | 0.82 | 0.76 |
| 0.3 | B | 0.05 | 0.05 | 0.05 | 0.03 |
| 0.3 | Bsub | 0.06 | 0.06 | 0.06 | 0.04 |
| 0.3 | AM2 | 0.99 | 0.94 | 0.82 | 0.78 |
| 0.3 | AMs | 0.98 | 0.94 | 0.82 | 0.78 |
| 0.3 | **AMf** | 0.98 | 0.91 | 0.77 | 0.70 |

Points that follow from the table:
- The plateau parameter, not the threshold, fixes the blind spot. B recalibrated to its own subset value (Bsub) moves 0.73 / 0.30 / 0.05 only to 0.75 / 0.34 / 0.06 in the main condition, and raises the K-world false STEP to 0.060 (`log:179`). Bsub is informational only.
- The full-grid threshold costs AM2 a few points on the blind spot (main condition: 0.98 / 0.96 / 0.94 at -2.403 -> 0.97 / 0.94 / 0.91 at -2.976) and 8 points at the noisiest corner (CV 0.15, h = 0.3: 0.78 -> 0.70). It still beats B by 0.24 / 0.64 at h = 0.2 / 0.25 in the main condition.
- It does not rescue large plateaus: at h = 0.4, s = 0.1 AMf calls STEP in 0.66 (main) and 0.39 (CV 0.15); at h = 0.5 the dominant verdict is FLAT (0.47) and NS 0.13 (`log:43, 47`).
- The cost of a plateau term sits on wide smooth steps with no plateau: at h = 0, s = 0.3 the STEP share is B 0.98 / 0.77 / 0.58 / 0.62 and AMf 0.94 / 0.67 / 0.46 / 0.53 (`log:81-100`, last column of the h = 0 block), a loss of 4 to 12 points. AMs loses 0 to 2 points there.

### 4.2 False STEP (`log:108-173`; limit 0.07)

| rule | R-world counting max (8 conditions, 400 runs) | R excl | R floor | K-world counting max (12 conditions, 500 runs, cap) | K no cap |
|---|---|---|---|---|---|
| B | 0.045 (`log:119`) | 0.000 | **0.100** (`log:121`) | 0.054 (`log:179`) | 0.054 |
| Bsub | 0.058 | 0.000 | 0.113 | 0.060 | 0.060 |
| AM2 (-2.403 R / -2.40 K) | 0.045 (`log:145`) | 0.005 | 0.040 | 0.052 | 0.052 |
| AMs (-2.463) | 0.035 (`log:158`) | 0.005 | 0.037 | 0.050 | 0.050 |
| AMf (-2.976) | 0.028 (`log:171`) | 0.005 | 0.028 | 0.042 | 0.042 |

The K-world counting maxima by condition are at `log:242-250` (AMf peaks at 0.042 in CV 0.05, misread 0.02). The floor group is not part of D2, but B's 0.100 there is above the 0.07 limit printed in the same table of the older log (`droplet_range_a.log:13, 25`), and every AM2 variant is below 0.07.

### 4.3 Section 10 ranges, K-world (`log:175-195`)

Percent of runs called STEP; counting max over the 87 cells and 12 conditions, then low-high over the cells of each row.

With the cap (`log:177-185`):

| row | B | Bsub | AM2 | AMs | AMf |
|---|---|---|---|---|---|
| counting max | 5.4 | 6.0 | 5.2 | 5.0 | 4.2 |
| s <= 0.15, misread <= 0.06 | 87-100 | 88-100 | 87-100 | 87-100 | 83-100 |
| s <= 0.15, misread 0.08-0.10 | 33-78 | 34-78 | 34-78 | 33-78 | 30-77 |
| s = 0.2, misread <= 0.06 | 80-100 | 81-100 | 79-100 | 79-100 | 73-100 |
| s = 0.2, misread 0.08-0.10 | 28-69 | 30-71 | 29-70 | 29-70 | 25-67 |
| s = 0.3, misread <= 0.06 | 43-98 | 46-98 | 43-98 | 43-97 | 36-95 |
| s = 0.3, misread 0.08-0.10 | 9-36 | 10-38 | 10-37 | 10-36 | 7-29 |

Without the cap (`log:187-195`; columns Bn, Bsubn, AM2n, AMsn, AMfn):

| row | Bn | Bsubn | AM2n | AMsn | AMfn |
|---|---|---|---|---|---|
| counting max | 5.4 | 6.0 | 5.2 | 5.0 | 4.2 |
| s <= 0.15, misread <= 0.06 | 91-100 | 92-100 | 92-100 | 91-100 | 87-100 |
| s <= 0.15, misread 0.08-0.10 | 81-97 | 83-97 | 82-97 | 80-97 | 76-96 |
| s = 0.2, misread <= 0.06 | 84-100 | 85-100 | 84-100 | 83-100 | 77-100 |
| s = 0.2, misread 0.08-0.10 | 68-91 | 71-92 | 69-92 | 68-92 | 60-88 |
| s = 0.3, misread <= 0.06 | 46-98 | 50-98 | 46-98 | 46-97 | 39-95 |
| s = 0.3, misread 0.08-0.10 | 37-54 | 41-56 | 38-54 | 37-53 | 30-44 |

Cell-by-cell differences are at `log:197-225` (AMf minus B, with and without the cap) and `log:227-240` (AMs minus B, cap). AMf minus B at s = 0.1 and 0.15 is within 5 points everywhere (worst -4.6, CV 0.15, misread 0.10, s = 0.15, no cap); at s = 0.2 the worst is -7.8 (no cap) and -6.8 (cap); at s = 0.3 the losses are 2 to 14 points (worst -13.6, CV 0.15, misread 0.02). AMs minus B stays within +1.2 and -3.6 points in every cell.

### 4.4 Bak et al. (2023) Fig. 2a re-fit, POST HOC (`log:252-257`)

The statistic `sharp` has no c*, so `droplet_am2_gate.log` Part 3 is valid at every threshold and was not rerun; the script re-averages its tables and asserts the paired numbers (`log:253`). Rate at which a sharp step is rejected when the truth is smooth (Bak-like, 500 runs): B 0.866 (misread <= 0.04) / 0.374 (>= 0.06), AM2 0.494 / 0.084; rate at which a sharp step is rejected when the truth is sharp: B 0.016 / 0.022, AM2 0.008 / 0.012 (`droplet_am2_gate.log:113-115`). The real Bak 90 C series (Tinc 90 C) is SMOOTH under B and under AM2 at all three thresholds (AM2 L between -0.57 and +0.06), and the main series is NS everywhere (`log:255-256`). My recomputation equals the stage-3 text for all 36 x 2 rows (`log:257`).

## 5. Scorecard (predictions scored as written)

| id | claim | result | verdict |
|---|---|---|---|
| F1 | reproduction at -2.4034 (0.98 / 0.96 / 0.94, B 0.73 / 0.30 / 0.05); AMf in [0.92, 0.98], [0.87, 0.96], [0.80, 0.93]; AMs in [0.96, 0.98], [0.93, 0.96], [0.90, 0.94]; AMf beats B by >= 0.20 and >= 0.55 | reproduction exact; AMf 0.9717 / 0.9433 / 0.9067; AMs 0.9800 / 0.9583 / 0.9400; gaps 0.243 / 0.638 (`log:260-268`) | HELD (8 of 8) |
| F2 | highest false STEP, AM2 0.045 (R) and 0.052 (K); AMs in [0.030, 0.045] and [0.038, 0.052]; AMf in [0.012, 0.035] and [0.017, 0.040] | 0.045, 0.035, 0.028 (R); 0.052, 0.050, **0.042** (K) (`log:269`) | FAILED (5 of 6; AMf K-world is 0.002 above its band) |
| F3 | Section 10 endpoints of AMf (cap) and AMfn in bands; AMs within 3 points below AM2; largest drop of AMf against B in [4, 11] (cap) and [3, 10] (no cap), AMs in [0, 3] | AMf cap 12 of 13 (counting max 4.2 vs [1.7, 4.0] missed); AMfn 11 of 13 (counting max 4.2, and s = 0.3 at misread 0.08-0.10 high 44 vs [45, 53]); AMs within 3 points: yes (largest gap 0.80 and 1.60); largest drops 7.20 (cap, in band), **10.40** (no cap, band ends at 10), 0.60 and 1.20 (`log:272-307`) | FAILED |
| F4 | Bak paired numbers unchanged; real Tinc 90 C labels SMOOTH and main labels NS at -2.40, -2.4634 and -2.9763 | all equal (`log:253-257`, `log:308`) | HELD |
| F5 | gates and structure hold | all hold (`log:309`) | HELD |

3 of 5 held. Every missed item leans the same way (AMf slightly worse than predicted: false STEP 0.002 and 0.2 point high, one power endpoint 1 point low, one cost 0.4 point high), which suggests the predictions' rough slope of 0.12 STEP share per unit of c* (`droplet_am2_full_predictions.txt:12`) was a little optimistic at the full-grid value. No miss changes a decision.

Interpretation choices disclosed: F3 endpoint bands are scored on the printed rounded percents (counting max with one decimal, rows as integers), with exact values printed beside them (`log:274-300`); the unbanded "high 100" no-cap entries are treated as [99, 100]; D1b is scored on exact shares (the rounded variant gives the same verdict: B 43 -> AMf 36 and Bn 54 -> AMfn 44 are both more than 5 points).

## 6. Decision rules (`log:311-315`)

| rule | AMf (-2.9763), decides | AMs (-2.4634), information |
|---|---|---|
| D0 gates and c* literals | HELD | HELD |
| D1a STEP(main, s = 0.1, h = 0.2) >= 0.90 | 0.9717: HELD | 0.9800: HELD |
| D1b every one of 12 endpoints satisfies AM >= B - 0.05, with cap and without | largest drops 7.20 and 10.40 points: **FAILED** | 0.60 and 1.20: HELD |
| D2 false-STEP counting max <= 0.07 (R, K cap, K no cap) | 0.028 / 0.042 / 0.042: HELD (as fixed; section 2) | 0.035 / 0.050 / 0.050: HELD |

Recommendation (tree in the predictions file, lines 57-59): D1a and D2 hold, D1b fails, so **keep B; add a limit sentence that names the plateau blind spot with B's measured shares; offer AM2 only as a labelled post hoc sensitivity analysis.** Expectation check: D1a holds (95% predicted) yes, D2 holds (95%) yes, D1b fails at -2.9763 (80%) yes, D1b holds at -2.4634 (85%) yes.

Two cautions on the subset value. AMs (-2.4634) looks free in this table, but the subset (true misread 0, .06, .10) was chosen after the full-grid result and is calibrated on 6 of the 12 conditions; the recipe's value is the full-grid one. And AMs is still only a simulation of misread levels the real controls have not measured.

## 7. Proposed prereg wording (NOT APPLIED)

Limit sentence for `droplet_prereg.md`:

> Limit of rule B. A step that sits on a plateau of the count, P = h + (1 - h) Phi(-ln(V / 3.1 uL) / s), is read as SMOOTH by rule B. In simulation at s = 0.1, B calls such a truth STEP in 0.50-0.82 of runs at h = 0.2, 0.21-0.30 at h = 0.25 and 0.03-0.05 at h = 0.3 (four noise conditions, 600 runs each; `droplet_am2_full.log:81-100`). A SMOOTH label from rule B therefore does not exclude a step on a plateau of 0.2 or more.

Sensitivity-analysis sentence, to be labelled post hoc wherever it appears:

> Post hoc sensitivity analysis. The data are refit with rule AM2 (a plateau h in both families) with threshold c* = -2.98, calibrated on the full 12-condition grid (`droplet_am2_gate.log:7`; not -2.40, which came from 8 conditions and 300 runs). In simulation AM2 at -2.98 calls STEP in 0.70-0.98 of runs at h = 0.3, s = 0.1 (four noise conditions) where B calls it in 0.03-0.05, with false STEP at most 0.028 (R-world) and 0.042 (K-world). It costs up to 7 points (with the control cap) to 10 points (without it) of STEP share on broad steps (s >= 0.2). A SMOOTH under AM2 rejects a sharp step in about 8% of runs at misread >= 0.06 on a Bak-like smooth truth (B: 37%), so it is not evidence against a sharp step and must not be quoted as such.

If AM2 is ever adopted in place of B, c* must be recalibrated at the control count actually used (at 480 controls the full-grid value is -2.798, `droplet_am2_gate.log:9`), and every AM2 column of `droplet_ctrl_am2.log` must be rerun first.

## 8. Notes on the predictions file (not edited)

Line-number blemishes in `droplet_am2_full_predictions.txt`: F1 says `droplet_range_a.log:151, 156, 161-ish`, the main-condition AM2 lines are 131, 135, 139 (B: 73, 77, 81); F2 says `droplet_ctrl_am2.log:157`, the counting-max line is 154; F4 says `droplet_am2_gate.log:113-116`, the paired rows are 113 (header), 114 (B), 115 (AM2). None changes a score.

Other things noticed: `droplet_ctrl_am2.py` rounds B's c* to -2.21 and AM2's to -2.40, so its P7 fails (56 of 60, misses at `droplet_ctrl_am2.log:67, 68, 69, 98`) and R1 says NOT TRUSTED; the gate here reproduces that text on purpose. `droplet_range_a.log:8-11` still prints AM2 c* = -2.40 as "amended"; the full-grid value is -2.976.

## 9. In plain words

Rule B (the pre-registered rule) is cautious about calling a STEP: it is right about a clean step (STEP in 100% of simulated runs at low misread) and about a smooth ramp, but it misses a step that sits on top of a flat floor of counts. If the floor is 20% of the maximum, B misses the step in a fifth to a half of runs; at 30% it misses it nearly every time (STEP 0.03-0.05). AM2 adds that floor as a free parameter and catches those steps (STEP in 0.70-0.98 of runs at a 30% floor, the low end being the noisiest condition), with false-STEP rates no higher than B's, and it also fixes B's separate 10% false STEP on flat-floor truths. With its honest threshold of -2.98 it pays for this in two places: it loses up to 7 to 10 points of STEP calls on broad steps (s >= 0.2), and a SMOOTH label from AM2 means much less, because on a smooth truth AM2 rules out a sharp step in only about 8% of the noisy runs (B: 37%). Neither rule handles a plateau of 40% or more well (B never calls it STEP; AM2 calls STEP in 0.39-0.75 of runs at 40% and 0.09-0.25 at 50%, where FLAT and NS take over), so the sentence for the prereg is: keep B, say plainly that it cannot see a step above a floor of about 20%, and run AM2 afterwards as a labelled sensitivity analysis, not as the verdict.
