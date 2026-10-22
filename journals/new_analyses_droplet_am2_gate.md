# Unit droplet_am2_gate: the unrounded c*, c* at 480 controls, and what AM2 does to a sharp-step rejection

Simulated drops, not data. Date 2026-10-21. `journals/droplet_prereg.md` was not edited (mtime 2026-10-20 23:34:26, older than the predictions file); no existing file was edited. Every wording below is **NOT APPLIED**.

## 1. Files, hashes, run

| file | sha256 |
|---|---|
| `droplet_am2_gate_predictions.txt` (65 lines, written 02:26:34, before any result) | `db59504ecfae5b78f71118b6d306aecd21b788beaf3bc993e01eaf0878251b0e` |
| `droplet_am2_gate.py` (330 lines, last edited 02:36:26, the version that ran; asserts the hash above) | `588073e7ef607b642d23781826060de5224232209f103e7d6d8ac8118ba22b2f` |
| `droplet_am2_gate.log` (135 lines) | `0aaad5454461de8245f4263fcce48ab05b764b37d1dab433b3e57a94c3396f5b` |
| `droplet_am2_gate.err` (183 lines, progress only, no Traceback) | `5a7ce299e348ea665de9efa90e9a743def329d6424f27f73ffe56a0e9b8df5e5` |

Run line (exit 0, started about 02:36:30, finished 02:46:25, 593 s, 4 processes). One invocation completed everything; the per-task checkpoints in the scratchpad (`.../scratchpad/am2_gate/`) were only a safeguard and no task was resumed:

```
UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_am2_gate.py > droplet_am2_gate.log 2> droplet_am2_gate.err
```

The script imports `droplet_range_common.py`, `misread_sim.py` and `droplet_ctrl_am2.py` without editing them. Before the main stages it checks itself (`droplet_am2_gate.err:1`): its calibration code equals `misread_sim.run_cell` on two cells (L to 1e-8, same eligibility), and its `sharp` equals the printed values of `bak_policy_b.log` (B) and `droplet_ctrl_am2.log` (AM2) at all 18 (n, misread) pairs. The 87 null curves came from the existing `DR_CACHE` file of the earlier unit (results do not depend on it).

Disclosed deviations and additions:
1. The default time budget in the script was raised from 480 s to 2400 s before the run, because the run went to the background (no 10 min foreground cap). This changes no number.
2. After the run, one POST HOC read of the saved calibration checkpoints (not a new simulation, not in the script) gave the c* values restricted to the misread rates 0, 0.06, 0.10 and the number of (cell, condition) pairs below given values; see 3.2. It is not scored against any prediction.
3. The predictions file (lines 8-9) fitted c*(n) = c_inf - a/sqrt(n) to the three printed points of `misread_sim.log:4-6` (64, 128, 192 controls) and extrapolated it to 480; the extrapolation is what K2 and K3 rest on and it was wrong (see 3.2).

## 2. Scorecard (predictions scored as written)

| id | claim | result | verdict |
|---|---|---|---|
| G1 | c*_B(128) in [-2.2100, -2.2050], cell dry a=1.67, CV 0.05, misread 0.04 | -2.205492, that cell (`log:6`) | HELD (strict form) |
| G2 | with the unrounded c*: >= 58/60 equal digits, the four old misses fixed, tolerance, counting max 0.054+-0.004, NS 0.00; with -2.21: the old 56/60 and the same four cells | Bu 60/60; Bn 56/60 with exactly the four cells (`log:92-93`); counting max 0.054; NS 0.000 | HELD |
| K1 | c*_AM2(128) in [-2.50, -2.15] and within 0.10 of -2.40 | -2.976290 (`log:7`) | FAILED (both clauses) |
| K2 | c*_B(480) in [-1.84, -1.64] | -2.019774 (`log:8`) | FAILED |
| K3 | c*_AM2(480) in [-1.95, -1.65]; AM2 - B in [-0.20, +0.05] | -2.797970 (`log:9`); difference -0.778 | FAILED (both clauses) |
| S1 | B smooth truth: >= 15/18 cells in their bands | 9/18 (`log:97-102`, `log:123`) | FAILED |
| S2 | AM2 smooth truth <= 0.65 everywhere and <= 0.35 for e >= 0.06 | max 0.918 (e = 0, n = 10); max 0.150 for e >= 0.06 (`log:98-102`) | FAILED (first clause) |
| S3 | AM2 <= B + 0.03; B - AM2 >= 0.30 for e <= 0.04 | AM2 never above B (max -0.076); smallest B - AM2 at e <= 0.04 is +0.076 | FAILED (second clause) |
| S4 | sharp-truth shares <= 0.15 for both rules | B max 0.040, AM2 max 0.024 (`log:106-111`) | HELD |

Held 3 of 9 (`log:128`). K1 note: its two clauses overlap only in [-2.50, -2.30]; both fail by a wide margin, so the overlap does not matter.

## 3. Results

### 3.1 Point 1: the failed gate P7 was rounding (`log:12-93`)

c*_B(128) with misread_sim's calibration reproduced exactly (87 cells x 12 conditions, 1000 runs, seeds [4, i, j], 192 control draws after the protocol draws): **-2.205492**, set by dry a=1.67 at CV 0.05, misread 0.04, the cell printed in `misread_sim.log:5`. Rerunning the 60-cell gate (500 runs, seeds [5, i, j], no control cap) on the same simulated data:

| c* | cells equal to the printed digits | cells that differ |
|---|---|---|
| -2.21 ("Bn", rounded) | 56/60 | counting max at CV 0.05, misread 0.06, 0.08, 0.10 (`log:19-21`, 1 run lower each); probit s = 0.15 at CV 0.05, misread 0.08 (`log:50`, 1 run lower) |
| -2.205492 ("Bu") | **60/60** | none |

Every cell is within the old tolerance (0.015 / 0.004) for both; counting maximum 0.054; NOT SCORABLE share 0.000 in all 12 conditions (`log:89-91`). So the four misses of `droplet_ctrl_am2.log:67-69, 98` are the two-decimal rounding of c*, not a different calibration. A prereg that quotes power tables to the digit must quote c* with at least four decimals (-2.2055).

### 3.2 Point 2: c* at 480 controls (`log:4-9`)

Same 1044 (cell, condition) pairs, same simulated runs, a second block of 480 control draws per run; e_hat = (k + 0.5)/481. Unrounded c* by the misread_sim recipe:

| rule | 128 controls | 480 controls | set by (480) |
|---|---|---|---|
| B | -2.2055 | **-2.0198** | dry a=1.67, CV 0.05, misread 0.04 (same cell as at 128) |
| AM2 | -2.9763 | **-2.7980** | run1 dry bulk front w 0.1 B 10, CV 0.15, misread 0.08 (same cell as at 128) |

Readings:
- B: 480 controls move c* by only 0.186. Together with the logged c*_B(192) = -2.03 (`misread_sim.log:6`), the worst cell is flat between 192 and 480 controls (-2.03 to -2.02) while the other conditions keep rising (per-condition minima at CV 0.15, misread 0.10: -1.64 at 128, -1.31 at 480). So c*_B at 480 is set by that cell's own spread, not by control noise; my prediction K2 (-1.74) extrapolated a 1/sqrt(n) law that does not hold there.
- AM2: -2.9763 is not -2.40. POST HOC from the saved calibration checkpoints (`droplet_am2_gate.log` has only the global value): the per-condition minimum of the 5th percentile at 128 controls is -2.98 for CV 0.15, misread 0.08, while CV 0.15, misread 0.06 and 0.10 are -2.05 and -2.29. Only 3 of 851 eligible (cell, condition) pairs are below -2.40 and one below -2.60; the second-lowest is -2.5667 (`log:7`). The -2.98 is one pair. Restricting the same checkpoints to misread 0, 0.06, 0.10 (the rates of `droplet_range_a`'s 8 conditions that this grid has) gives -2.4634 at 128 controls (|difference to -2.40| = 0.06) and -2.3433 at 480. So -2.40 came from a condition grid that skips misread 0.08 and from 300 runs; the 12-condition, 1000-run recipe gives -2.98 for AM2 and has no reason to be different for B (B at 128 on the 8-like subset is -2.0588 against -2.2055 for the full grid, so the grid matters for B too, in the same direction).
- D2 (`log:132`): |c*_B(480) - c*_B(128)| = 0.186, below the 0.25 cut, so a c* tied to the control count is not forced by this calibration; AM2 at 480 differs from B at 480 by 0.778, above 0.05, so AM2 needs its own c*. These are calibration values only; false STEP and power at 480 controls with them were not evaluated.

### 3.3 Point 3: rejecting a sharp step when the truth is smooth (`log:95-115`)

Bak (2023) Tinc 90 C series (11 levels, n = 3, 6 or 10 at the seven unprinted P = 1 markers), misread e known, 500 runs per cell, seeds [31, n idx, e idx, truth idx], share with sharp > 3.84 (SE <= 0.022). Smooth truth = best COUNT curve fitted to the series; sharp truth = best STEP curve with s <= 0.15.

Smooth truth (the number that matters: how often a sharp step is rejected; `log:97-102`):

| n | rule | e=0.00 | 0.02 | 0.04 | 0.06 | 0.08 | 0.10 |
|---|---|---|---|---|---|---|---|
| 3 | B | 1.000 | 0.930 | 0.832 | 0.588 | 0.506 | 0.302 |
| 3 | AM2 | 0.764 | 0.450 | 0.282 | 0.150 | 0.124 | 0.056 |
| 6 | B | 0.998 | 0.858 | 0.756 | 0.484 | 0.378 | 0.216 |
| 6 | AM2 | 0.890 | 0.406 | 0.244 | 0.102 | 0.098 | 0.034 |
| 10 | B | 0.994 | 0.780 | 0.648 | 0.440 | 0.252 | 0.204 |
| 10 | AM2 | 0.918 | 0.284 | 0.204 | 0.110 | 0.048 | 0.034 |

Sharp truth (paired number: a true sharp step wrongly rejected; `log:106-111`): B 0.000-0.040, AM2 0.000-0.024 in all 18 cells; no cell above 0.15 for either rule (D4: the "sharp > 3.84" test is at or below the 5% level in this setting, for both).

Paired summary (mean over n = 3, 6, 10; `log:113-115`):

| rule | smooth truth, e <= 0.04 | smooth truth, e >= 0.06 | sharp truth, e <= 0.04 | sharp truth, e >= 0.06 | sharp truth max |
|---|---|---|---|---|---|
| B | 0.866 | 0.374 | 0.016 | 0.022 | 0.040 |
| AM2 | 0.494 | 0.084 | 0.008 | 0.012 | 0.024 |

Readings:
- AM2 rejects a sharp step 0.76-0.92 of the time at e = 0 but only 0.20-0.45 at e = 0.02-0.04 and <= 0.15 at e >= 0.06. So its loss of power is a misread effect: at zero misread AM2 is nearly as sharp as B (B 0.99-1.00). The sharp-truth side gets smaller for AM2 as well (it also accepts a true step), so AM2 is the more conservative test in both directions, not only the one that accepts smooth series.
- B itself rejects less than predicted at n = 6 and n = 10 (n = 10: 0.780 at e = 0.02, 0.648 at e = 0.04; S1 failed in 9 cells, every one on the low side) and falls to 0.20-0.60 once e >= 0.06. Even rule B does not reject a sharp step reliably in smooth Bak-like data when the misread rate is 6-10%.
- D3 (`log:133`): AM2 is at most 0.50 in 15 of 18 cells (the three exceptions are e = 0) and B is at least 0.75 in 8 of the 9 cells with e <= 0.04, so the stage-3 reading of `new_analyses_droplet_ctrl_am2.md` holds in simulation, not only on the one dataset.
- Not done: the plateau-step truth (a sharp step with a plateau) was not simulated, so the claim that AM2 does not reject it rests on construction, not on a run; control-drop noise was not included in the Bak study (e known).

## 4. Decision rules (as fixed in the predictions file; `log:130-135`)

- D1: G1 and the gate hold with the unrounded c* (60/60), so the rounding of c* explains the misses of P7 and R1 "NOT TRUSTED" is lifted for the B column. The AM2 columns of `droplet_ctrl_am2` used c* = -2.40; this unit's recalibration gives -2.9763 (more than 0.10 away), so by the rule they must be rerun with the recalibrated c* before being quoted. The POST HOC subset above says that the large gap is one (cell, condition) pair at misread 0.08, so which value is "right" is a choice of the calibration grid, not a settled number; I do not choose it here.
- D2: see 3.2. D3: see 3.3. D4: no sharp-truth share above 0.15. D5: nothing edited.

## 5. Proposed prereg wording (NOT APPLIED)

1. Quote c* with four decimals and name its grid: "c*_B = -2.2055 (87 counting and drying cells x 12 conditions, CV 0.05/0.15, true misread 0, .02, .04, .06, .08, .10, 1000 runs per cell, 128 controls)".
2. If the control count is raised to 480 (W1): "c*_B(480) = -2.0198 and c*_AM2(480) = -2.7980 are calibration values only; false STEP and power at 480 controls were not evaluated."
3. If AM2 is adopted, replace any claim that a SMOOTH label rejects a sharp step by: "In simulated smooth series built from the Bak 90 C design (n = 3, 6, 10, misread known), rule B rejects a sharp step (sharp > 3.84) in 0.87 of runs at misread <= 4% and 0.37 at >= 6%; rule AM2 in 0.49 and 0.08; a true sharp step is rejected in at most 0.04 (B) and 0.02 (AM2) of runs. A SMOOTH label under AM2 therefore means 'smooth, or a step with a plateau', not 'not a sharp step'."
4. Before AM2 columns of any table are quoted: the c*_AM2 grid (misread 0.08 included or not) must be fixed in the prereg, because it moves c*_AM2(128) between -2.46 and -2.98.
