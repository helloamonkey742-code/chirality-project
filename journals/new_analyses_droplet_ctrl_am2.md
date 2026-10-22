# Unit droplet_ctrl_am2: how many control drops, and rule AM2 against rule B

Simulated drops and exact binomial arithmetic, not data. Date 2026-10-21. `journals/droplet_prereg.md` was not edited (its mtime, 2026-10-20 23:34, is older than the predictions file). Every wording below is **NOT APPLIED**.

## 1. Files, hashes, run

| file | sha256 |
|---|---|
| `droplet_ctrl_am2_predictions.txt` (written 2026-10-21 01:53:18 PDT, before any result) | `e60efb84dfc449661160e4f6becdba3dd277cc7753c9cbc77e6c3c5fd22f1a32` |
| `droplet_ctrl_am2.py` (311 lines, last edited 02:02:25, the version that ran) | `ea8219d9ae0e021feb792234d2ddd1170ec361e5e4ef8bd9c93ed5b6a2126543` |
| `droplet_ctrl_am2.log` (191 lines) | `1d23ce91feb246a5b5e6b6a95ffb08ba74b3b8fc1142bebae633d882179f1caa` |
| `droplet_ctrl_am2.err` (94 lines, progress only, no Traceback) | `85cde35cc197040e7676d38cd9c5fba896e96b60aa9d6feb5df0d8a2a278702f` |

Run line (as executed, exit 0, started 02:02:29 PDT, stage 2 simulation 170 s, total 172 s, 4 processes):

```
DR_CACHE=<scratchpad>/ctrl_am2/cache UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_ctrl_am2.py > droplet_ctrl_am2.log 2> droplet_ctrl_am2.err
```

Part 2 uses `droplet_range_common.py` (imported, not edited), the seeds `[5, i, j]` and 500 runs per cell of `misread_sim.py`.

## 2. Part 1: how many control drops (exact binomial, `scipy.stats.binom`)

K ~ Binomial(n, p), p the true misread rate. NOT SCORABLE (NS) iff K > floor(n/10). False SCORABLE (FS) = 1 - NS. The prereg states no cost for a control drop (Sections 6-8), so only counts are given: n extra drops on top of the 200 main drops. Sawtooth: n = 10m .. 10m+9 share one threshold, so use a multiple of 10 (here a multiple of 80 keeps 8 levels equal and the threshold whole).

Selected rows of `droplet_ctrl_am2.log:6-29` (all n = 128..1024 were computed):

| n | log line | NS 6% | NS 8% | NS 10% | FS 11% | FS 12% | FS 13% | FS 14% |
|---|---|---|---|---|---|---|---|---|
| 128 | 6 | 0.044 | 0.224 | 0.519 | 0.339 | 0.222 | 0.136 | 0.079 |
| 129 | 7 | 0.047 | 0.232 | 0.531 | 0.328 | 0.213 | 0.129 | 0.074 |
| 160 | 11 | 0.016 | 0.141 | 0.434 | 0.402 | 0.261 | 0.156 | 0.085 |
| 240 | 13 | 0.005 | 0.107 | 0.446 | 0.356 | 0.198 | 0.096 | 0.041 |
| 400 | 16 | 0.001 | 0.062 | 0.458 | 0.293 | 0.123 | 0.040 | 0.010 |
| 480 | 17 | 0.000 | 0.048 | 0.462 | 0.269 | 0.098 | 0.026 | 0.005 |
| 500 | 18 | 0.000 | 0.045 | 0.462 | 0.264 | 0.093 | 0.024 | 0.005 |
| 509 | 19 | 0.000 | 0.059 | 0.516 | 0.221 | 0.072 | 0.017 | 0.003 |
| 640 | 23 | 0.000 | 0.029 | 0.467 | 0.230 | 0.065 | 0.012 | 0.001 |
| 800 | 25 | 0.000 | 0.018 | 0.470 | 0.199 | 0.043 | 0.005 | 0.000 |
| 1000 | 26 | 0.000 | 0.010 | 0.473 | 0.169 | 0.027 | 0.002 | 0.000 |
| 1024 | 29 | 0.000 | 0.011 | 0.490 | 0.155 | 0.023 | 0.002 | 0.000 |

Smallest n meeting the target (`log:33-39`). "first" = smallest n; "stable" = smallest n after which it holds for every larger n up to 1024:

| target | first | stable | n multiple of 10: first / stable | n multiple of 80: first / stable | min..max over n |
|---|---|---|---|---|---|
| NS(6%) <= 0.05 | 128 | 128 | 130 / 130 | 160 / 160 | 0.000 .. 0.047 (n=129) |
| NS(8%) <= 0.05 | 470 | 550 | 470 / 470 | 480 / 480 | 0.010 .. 0.232 |
| NS(10%) <= 0.05 | none | none | none | none | 0.427 (n=130) .. 0.531 (n=129) |
| FS(11%) <= 0.10 | none | none | none | none | 0.144 (n=1019) .. 0.425 |
| FS(12%) <= 0.10 | 389 | 471 | 480 / 480 | 480 / 480 | 0.021 .. 0.294 |
| FS(13%) <= 0.10 | 179 | 231 | 240 / 240 | 240 / 240 | 0.001 .. 0.189 |
| FS(14%) <= 0.10 | 128 | 141 | 150 / 150 | 160 / 160 | 0.000 .. 0.114 |

Answers.
- **True 6%:** 128 controls already give NS <= 5% (4.4%; worst over n is 4.7%). No extra controls needed.
- **True 8%:** 480 controls (60 per level) give NS 4.8%. Any multiple of 10 from 470 up works; 550 is the first n from which every n works. That is 352 more control drops than now.
- **True 10%:** no n up to 1024 gets NS under 5%. The cap sits at the mean, so NS stays 0.43-0.53 for every n. No count fixes this; only a different cap does.
- **False SCORABLE at true 11%:** never under 10% with a 10% cap (0.339 at 128; 0.144 at best, n=1019). At 12% it needs about 480 controls, at 13% about 240, at 14% about 150.
- **Lower cap** (`log:41-46`): the pair NS(6%) <= 0.05 and FS(11%) <= 0.10 is met only by a cap of 8% (stable n = 438) or 9% (stable n = 435). Caps of 6% and 7% never meet it up to 1024; the 10% cap meets NS(6%) but never FS(11%). The price (reported, not scored): cap 8%, n=438: NS 0.036 / 0.458 / 0.910 at a true 6% / 8% / 10%, FS(11%) 0.023 (`log:44`). Cap 9%, n=435: NS 0.005 / 0.201 / 0.735, FS(11%) 0.098 (`log:45`). The 10% cap stays per the task, so these are for information.
- n here is a fixed count. Section 6 counts only readable controls, so the realised n is smaller and random.
- A different n changes the c* that Section 2 needs (B: -2.61 / -2.21 / -2.03 at 64 / 128 / 192, `misread_sim.log:4-6`). Not simulated here for n = 480.

## 3. Part 2: B reproduction and AM2 against B

500 runs per cell, 12 conditions (pipette CV 5 or 15% x true misread 0, 2, 4, 6, 8, 10%), 128 controls. "Bn" and "AM2n" = no 10% cap (what `misread_sim.log` and Sections 3 and 10 print); "B cap" and "AM2 cap" = Section 2 item 1 as written. B uses c* = -2.21, AM2 c* = -2.40. The same simulated data is scored under all four.

**B against the log** (`droplet_ctrl_am2.log:64-135`, column "Bn-log"): 56 of 60 cells equal the printed digits. The four that do not:

| cell | log | Bn | log line |
|---|---|---|---|
| counting max, CV 0.05, misread 0.06 | 0.028 | 0.026 | 67 |
| counting max, CV 0.05, misread 0.08 | 0.044 | 0.042 | 68 |
| counting max, CV 0.05, misread 0.10 | 0.034 | 0.032 | 69 |
| probit s=0.15, CV 0.05, misread 0.08 | 0.96 | 0.95 | 98 |

All four are within the stated tolerance, Bn is the lower value in each, and the counting maximum (0.054) and the NOT SCORABLE share (0.00 everywhere without a cap, `log:139-150`) match. The Section 3 / 10 summary rows are identical for the log and Bn (`log:154-160`).

**AM2n against Bn** (`log:154-160`, percent of runs called STEP):

| | log (B 128) | Bn | AM2n |
|---|---|---|---|
| counting max, all 12 conditions | 5.4 | 5.4 | 5.2 |
| s <= 0.15, misread <= 0.06 | 91-100 | 91-100 | 92-100 |
| s <= 0.15, misread 0.08-0.10 | 81-97 | 81-97 | 82-97 |
| s = 0.2, misread <= 0.06 | 84-100 | 84-100 | 84-100 |
| s = 0.2, misread 0.08-0.10 | 68-91 | 68-91 | 69-92 |
| s = 0.3, misread <= 0.06 | 46-98 | 46-98 | 46-98 |
| s = 0.3, misread 0.08-0.10 | 37-54 | 37-54 | 38-54 |

Cell by cell (`log:79-135`, 48 step cells): AM2n minus Bn lies in [-0.024, +0.014]; for s <= 0.15 the largest |difference| is 0.012 (`log:164`). The largest falls are at s = 0.3 (-0.02 at `log:127`, `129`, `131`). The AM2n counting maximum is 0.052 and per condition AM2n differs from Bn by at most 0.008 (`log:165`, `64-75`).

**The 10% cap as Section 2 writes it** (`log:139-160`). The NOT SCORABLE share is the same for B and AM2, close to the binomial (`log:139-150`):

| true misread | NS share, CV 0.05 / CV 0.15 | binomial (Part 1) |
|---|---|---|
| 0.04 | 0.000 / 0.006 | 0.002 |
| 0.06 | 0.052 / 0.044 | 0.044 |
| 0.08 | 0.224 / 0.232 | 0.224 |
| 0.10 | 0.530 / 0.522 | 0.519 |

Percent of all runs called STEP (NOT SCORABLE runs count as not STEP), B cap against no cap (`log:154-160`): s <= 0.15 at misread 0.08-0.10 falls from 81-97 to 33-78 (AM2 cap 34-78); s = 0.2 from 68-91 to 28-69; s = 0.3 from 37-54 to 9-36; at misread <= 0.06 s <= 0.15 goes from 91-100 to 87-100. Example (`log:99`): s = 0.15, CV 0.05, misread 0.10: 0.92 without the cap, 0.41 with it. The counting maximum is unchanged by the cap (5.4, `log:154`).

## 4. Stage 3 (POST HOC, no predictions): Bak et al. 2023 re-fit

B reproduces `bak_policy_b.log`: 36 of 36 (L, sharp, label) entries equal (`log:189`). AM2 (`log:176-187`):
- 90 C series: SMOOTH under both rules at all 6 misread rates x 3 count guesses. B: L = +8.14 .. -0.02, sharp (2 x delta log-likelihood against a step with s <= 0.15) 282.7 .. 1.2, so a sharp step is rejected (> 3.84) at 0-6% and at 8% for n = 3 and 6. AM2: L between -0.57 and +0.06 and sharp at most 3.2 (`log:177`), so a sharp step with a plateau is not rejected at any misread rate, including 0%. The SMOOTH label survives but the series no longer separates COUNT from STEP.
- Main series: NS at every rate under both rules.

## 5. Predictions scored as written (`droplet_ctrl_am2.log:49-55, 163-168, 191`)

| | outcome | detail |
|---|---|---|
| P1 | HELD | 6%: first 128, stable 128, max 0.047 at n=129 |
| P2 | HELD | 8%: first 470, stable 550, sawtooth present |
| P3 | HELD | 10%: min 0.427, max 0.531, no n reaches 0.05 |
| P4 | HELD | 11%: FS(128) 0.339, min 0.144, no n reaches 0.10 |
| P5 | **FAILED** | first n with FS <= 0.10: 12% 389 (band 420-600), 13% 179 (band 200-330), 14% 128 (held); FS(14%, 128) 0.079 (held) |
| P6 | HELD | pair stable n by cap: 6% none, 7% none, 8% 438, 9% 435, 10% none (band [330, 520] held); price clause reported, not scored |
| P7 | **FAILED** | 56/60 cells equal the printed digits (need 58); all within tolerance; counting max 0.054; Bn NS 0.00 everywhere |
| P8 | HELD | 48/48 step cells within 0.05; range [-0.024, +0.014]; s <= 0.15 max 0.012 |
| P9 | HELD | AM2n counting max 0.052; per-condition difference at most 0.008 |
| P10 | HELD | cap NS: see section 3 table; B equals AM2 exactly; at misread <= 0.04 at most 0.006 |
| P11 | HELD | no violation; example cell 0.41 (band 0.30-0.50) |
| P12 | HELD | cap changes the counting maxima by +0.000 for both rules |

10 of 12 held. Decision rules:
- **R1 applies:** the gate P7 failed, so the AM2 comparison is **NOT TRUSTED** and no prereg wording is proposed from it.
- R2/R3 are therefore not evaluated as verdicts. For information only: the R2 conditions (48 of 48 step cells within 0.05, AM2 counting max <= 0.07) are numerically met (section 3).
- R4 applies to Part 1 (section 2 and W1/W2 below). R5 holds.
- Why P7 failed, as an inference and not a verdict: all four misses are one run or so in 500 and Bn is the lower value in each, which is the direction expected when `misread_sim.py` used its unrounded calibrated c* (somewhere in -2.215 .. -2.205, printed as -2.21) and this run used exactly -2.21. The unrounded value is not printed in `misread_sim.log`, so this was not tested.

## 6. Proposed prereg wording (NOT APPLIED)

**W1 (prereg lines 32-34, 114, 123-128, replacing 128 controls with 480; needs c* refitted at 480 first):** "Control drops: 480, 60 at each of the 8 levels. D2 is NOT SCORABLE if more than 48 of the readable controls are read wrong. A run whose true misread rate is 8% or lower is called NOT SCORABLE at most 5% of the time (4.8% at 8%, under 0.1% at 6%). At a true 10% it is about 46%, and no number of controls changes that, because the cap is the rate. A run whose true rate is 12% passes the cap 10% of the time, 13% 3%, 14% 0.5%, and 11% 27% of the time. Unreadable controls lower the count and shift these figures." (figures: `log:17`)

**W2 (if 128 stays; near line 74 or line 123):** "With 128 controls and the 10% cap, a run whose true misread rate is 6% is called NOT SCORABLE 4% of the time, 8%: 22%, 10%: 52%. A true rate of 11% passes the cap 34% of the time, 12%: 22%, 13%: 14%, 14%: 8%." (figures: `log:6`)

**W3 (prereg line 74 and lines 141-143, independent of AM2):** the printed "81-97% at 8-10%" and "68-91%" are without the Section 2 cap. Replace by "with the cap, counting NOT SCORABLE runs as not STEP: steps with s <= 0.15 were called STEP in 87-100% of runs at misreads up to 6% and 33-78% at 8-10%; s = 0.2: 80-100% and 28-69%" (figures: `log:155-158`, B cap column), or state that the percentages are of scorable runs and recompute them.

**W4 (AM2):** none, per R1. Conditional only (not a verdict): if the author accepts the section 5 rounding inference, the sentence would be "with the plateau term (c* = -2.40) the Section 3 and 10 figures move by at most 0.024; the counting maximum is 5.2%". It must come with a rewrite of Section 10, lines 275-279: under AM2 the Bak 90 C series is SMOOTH but a sharp step is not rejected at any misread rate (section 4).

## 7. Disclosures

- A timing smoke run (`<scratchpad>/ctrl_am2/time1.py`, no verdicts printed) came after the predictions file.
- The null-curve cache in `.../scratchpad/ctrl_am2/cache/null.npz` is a copy of an earlier agent's cache of `droplet_range_common` (curves only, no results).
- An earlier launch, minutes before 02:02:29, aborted in `C.self_check()` (`droplet_range_common.py:309`, AssertionError) because my patch of `C.controls` ignored the `mode` argument. No result was produced and the log was empty. The patch is now mode-aware, the predictions file is unchanged, and the script hash above is the fixed version.
- Part 1 treats n as a fixed number of readable controls; c* was not refitted for any n other than 128; the P6 price clause had no numeric band and is reported, not scored.
- Stage 3 was run after the predictions were written and has no prediction.
- Pre-existing modifications in the working tree before this unit (README.md, scorecard.py, journals/ija_draft.md, mtimes 00:07-00:09) were not touched. `ASDRP_PITCH.md` was not opened. No git add, commit or push; no network.
