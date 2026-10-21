# h14e_advect: does the eddy slowdown keep falling in exponent at a molecular share of 1/256? (2026-10-21)

Script `h14e_advect.py` (imports `h14d_advect.py`, which imports `h14c_advect.py`, `h14b_advect.py` and `advect2d.py`; none of them was edited). Predictions `h14e_predictions.txt` (v1, written and hashed before any 1/256 result existed, sha256 `6ea98ce64f21fa6ccc67c9eeb45a3763aefe168a9116a2c8ca7c91a1724cd865`). The script asserts that hash and the sha256 of `h14d_advect.py` (`10476e92..0721`), `advect2d.py` (`7d84a2e6..145cf2e`) and `h14d_advect_raw.pkl` (`c3c8f8b9..0288`); the full strings are at `h14e_advect.log:1` and again at `:40`. The script file itself was written after the predictions were hashed and is not hash-asserted. Run line, exactly as given:

```
cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14e_advect.py > h14e_advect.log 2> h14e_advect.err
```

Log `h14e_advect.log` (111 lines, exit 0, `h14e_advect.err` empty), raw results `h14e_advect_raw.pkl`. `h14e_advect_raw.pkl.part` is a leftover save file next to it (11 bytes larger than the final pickle, `cmp` shows the first difference at byte 5,037,188; not inspected further; the final pickle is the one used). No post hoc script was run; the few POST HOC items below are arithmetic on numbers printed in the log. Nothing was committed, pushed or sent; no network. `ija_draft.md`, `README.md`, `scorecard.py`, `droplet_prereg.md`, `advect2d.py` and every `h14b_`, `h14c_`, `h14d_` file were not touched; `h14d_advect_raw.pkl.part` was left alone.

Output hashes (sha256): `h14e_predictions.txt` 6ea98ce64f21fa6ccc67c9eeb45a3763aefe168a9116a2c8ca7c91a1724cd865, `h14e_advect.py` 22e7b592ddfdc859ba87d248725e470a9a96a8646f8548116b7ef2ab421188a2, `h14e_advect.log` 576b12d85d9662208f864e415bd8a6b4b65d3aee19c37cc2d60da0f87f8b8eb9, `h14e_advect_raw.pkl` 3338b7347336e71d13e6fee2148f4d3a3582bafbef1fdc23ae39219a86b108f4.

## What I asked, in plain words

h14d measured the eddy slowdown S50 of wall clearing at 10 wall widths (l/w = 10) at molecular shares D_mol/D_eff of 1/8, 1/16, 1/32, 1/64, 1/128: 5.76, 7.71, 10.27, 13.05, 15.28 times Eq. 17's time t17. The local ratio exponent fell, 0.42, 0.41, 0.34, 0.23 [0.13-0.28], and the pre-registered label was UNDECIDED with the local slope FALLING (`h14d_advect.log:112-115`, `:147-148`). A post hoc log law and a power law whose exponent keeps falling agree at 1/128 and give nearly the same value at 1/256; a power law held at 0.36 or 0.40 gives a visibly larger one. h14e adds one cell, V256_10 (1/256, l/w 10), and asks whether the exponent keeps falling. The 1/64 (378 runs) and 1/128 (336 runs) cells are reloaded from `h14d_advect_raw.pkl`, not rerun.

Pre-registered rule: q = q(128 -> 256) = ln(S50(1/256) / S50(1/128)) / ln 2 with its 95% bootstrap interval; SATURATING if the upper bound < 0.2, GROWING if the lower bound > 0.3, otherwise UNDECIDED. Second label: CONSTANT-SLOPE EXCLUDED if the 95% upper bound of S50(1/256) is < 15.28 x 2^0.36 = 19.61, else NOT EXCLUDED. Local slope (as h14d): FALLING / RISING / UNRESOLVED from the 95% interval of Delta q = q(128 -> 256) - q(64 -> 128). Stated in advance: with a half-width of about 0.075 the label was expected to be UNDECIDED again (about 75%), and the second label EXCLUDED under a log law or a falling power law (about 85-90%).

## What ran

| case | D_eff | D_mol/D_eff | l/w | box L / cells / h | gamma | Tk | A | kappa | kicks to the end of the sweep | runs, walled at freeze-out | hold (t17) | seeds | cpu s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V64_10 (h14d pooled, reloaded) | 64 | 1/64 | 10 | 503 / 503 / 1 | 0.0073 | 27.4 | 117.5 | 1.468 | 5.00 | 378, 378 | 48 | h14d | 2811 |
| V128_10 (h14d, reloaded) | 128 | 1/128 | 10 | 711 / 711 / 1 | 0.0074 | 27.2 | 166.2 | 1.469 | 4.97 | 336, 336 | 60 | h14d | 5354 |
| **V256_10 (new)** | 256 | 1/256 | 10 | 1005 / 1005 / 1 | 0.0074 | 27.0 | 234.7 | 1.467 | 5.01 | 336, 336 | 60 (7826 time units) | 3401-3416 | 12345 |

(`h14e_advect.log:46-48`, flow `:43`; jobs `:8-23` and `:24-39`.) Sixteen jobs of 21 runs (7 seeds x 3 sweep-end biases Delta = 0.5, 1, 2), 3 workers under `nice -n 10`, wall clock 4,374 s = 73 min against a 7,200 s guard, no job cut, all 336 runs finished (`:106`). cpu 12,345 s against the 12,000 s of the cost estimate in the predictions; the guess for wall time was 4,100 s. Job times 516 to 876 s (`:24-39`). t17 = 130.4.

The 1/64 and 1/128 cells and the 1/8, 1/16, 1/32 cells are the h14c/h14d runs; `:4-5` show that the h14c values (5.76, 10.27, 13.21) and the h14d ladder (5.76, 7.71, 10.27, 13.05, 15.28) are reproduced from the pickles with their own bootstrap seeds before the new cell is scored. Self-checks at `:2-3`: h14c's `self_check`, the recipe at D_eff 128 and 64 reproduces the h14d V128_10#0 and V64_10n#3 specs exactly, seeds unique and outside 3101-3132 and 3301-3344, D_eddy = D_eff - 1 by construction.

Disclosed (also in `h14e_predictions.txt`, DISCLOSURE, and in the script): (i) the h14d ladder, local exponents, post hoc log law and the 1/128 anchor 15.28 were known when the bands were written; (ii) a throwaway timing script in the session scratchpad, not in the repo (`h14e_timing.py`), called `h14c_advect.job` on the 1/256 and 1/128 flows and printed only elapsed seconds, steps, ns per cell-step and memory; no S50, S80, walled fraction, p or q was computed or seen for any 1/256 run; I also read elapsed and clearing times of the h14d 1/128 jobs (h14d data) to turn the ns figure into a cost; (iii) after hashing, a code-path test of the report section (`h14e_report_test.py`, scratchpad, not in the repo) was run on h14d V128 jobs relabelled as V256_10, so that the printing and scoring code ran once before the results run; that output is not a 1/256 result and is not cited; (iv) the 336 new runs are pooled in job order (#0 to #15), not in completion order (e.g. #8 finished before #7 and #11 before #10, `:31-35`), so the bootstrap does not depend on worker timing; the point estimates do not depend on the order; (v) the anchor 15.28 in the second label is taken as a fixed number (its own 95% interval 14.93-15.79 is not carried), and the paired version (q upper bound < 0.36) is printed next to it.

## Results

S50, S80 and S_half at l/w 10 (time since the end of the sweep / t17; Kaplan-Meier over runs walled at freeze-out; every cell 100% walled and every run cleared inside the hold; bootstrap over runs, 2000x; `h14e_advect.log:51-62`):

| D_mol/D_eff | runs | S50 [95%] | S80 [95%] | S_half [95%] | kicks to clear (S50 t17 / Tk) | mean domains at freeze-out | log line |
|---|---|---|---|---|---|---|---|
| 1/8 | 96 | 5.76 [4.94-6.44] | 8.86 [7.63-9.78] | 3.54 [3.07-4.03] | 24 | 3.43 | `:51-52` |
| 1/16 | 192 | 7.71 [7.34-8.40] | 12.02 [10.75-13.00] | 3.93 [3.49-4.36] | 35 | 4.93 | `:53-54` |
| 1/32 | 96 | 10.27 [9.48-11.32] | 13.62 [12.71-16.04] | 6.38 [5.75-7.01] | 48 | 7.91 | `:55-56` |
| 1/64 | 378 | 13.05 [12.74-13.80] | 16.48 [16.11-17.26] | 7.93 [7.69-8.34] | 62 | 12.50 | `:57-58` |
| 1/128 | 336 | 15.28 [14.93-15.79] | 18.79 [17.86-20.00] | 9.12 [8.72-9.55] | 73 | 22.77 | `:59-60` |
| **1/256 (new)** | 336 | **18.07 [17.66-18.52]** | **21.30 [20.41-21.91]** | 9.31 [9.07-9.71] | 87 | 42.74 | `:61-62` |

Fractions cleared by 12 / 16 / 24 / 32 t17 at 1/256: 0.01, 0.29, 0.90, 0.98 (`:62`); the walled fraction at freeze-out is 1.000 and no walled run is still walled at the end of the 60 t17 hold (`:48`, `:105`), so no quantile is censored.

Local ratio exponents q(a -> b) = ln(S50_b / S50_a) / ln(r_b / r_a) (`:65-75`):

| step | q [95%] | log line |
|---|---|---|
| 1/8 -> 1/16 | 0.42 [0.24-0.70] | `:65` |
| 1/16 -> 1/32 | 0.41 [0.24-0.56] | `:66` |
| 1/32 -> 1/64 | 0.34 [0.21-0.50] | `:67` |
| 1/64 -> 1/128 | 0.23 [0.13-0.28] | `:68` |
| **1/128 -> 1/256** | **0.24 [0.19-0.29]** | `:69`, `:73` |
| q(8 -> 64), least squares over 8, 32, 64 | 0.40 [0.34-0.48] | `:70` |
| **q(64 -> 256), least squares over 64, 128, 256** | **0.23 [0.19-0.26]** | `:71` |
| q over all six ratios, least squares | 0.33 [0.30-0.36] | `:72` |

Delta q = q(128 -> 256) - q(64 -> 128) = 0.01 [95% -0.08 to 0.14], paired over the same replicates (`:73`). The 95% interval of q(128 -> 256) excludes 0.5 (cellular) and excludes 0 (`:74`). Curvature of ln S50 against ln r over the six ratios: -0.04 [-0.08, -0.02] (`:75`). A job-level bootstrap (printed, not the scored interval) gives q(128 -> 256) 0.18-0.29 (`:76`).

Per-doubling increments of S50: 1.94, 2.57, 2.77, 2.24, 2.78; the last (128 -> 256) is 2.78 [2.16-3.39] against a mean of the previous four of 2.38; a power law with q 0.36 would give 4.33 (`:78`). Log law refitted over six ratios: S50 = -1.98 + 3.59 ln r, b 3.39-3.82, rms 0.20 (the h14d five-point fit was -1.77 + 3.52 ln r, rms 0.19); power-law refit over six ratios: q 0.33, rms 0.66, residuals -0.37, -0.00, +0.58, +0.87, -0.03, -1.18 (`:79-80`). Both refits are printed, not scored.

Scenarios written in the predictions against the measured 18.07 [17.66-18.52] (`:83-86`): A log law 17.75 inside; B falling power (q 0.23 from 1/128) 17.92 inside; C power law q 0.36 from 1/128 19.61 outside (measured - scenario = -1.55); D power law q 0.40 from 1/128 20.16 outside (-2.10). A and B cannot be told apart by this cell (they differ by 0.17, the half-width is about 0.4); it separates {A, B} from {C, D}, as designed.

Regime (`:43`): kappa = 1.467, kicks to the end of the sweep 5.01, Lagrangian D_eddy 1.011 times nominal (257.71 against 255), numerical kick diffusion 0.0005 D_mol, Lyapunov exponent 0.130 per kick (0.0048 per time unit), the same values as 1/16, 1/64 and 1/128 in h14d.

Predictions (`:98-106`): P1 to P9 all held. P1 S50 18.07 in [16.2, 20.8] (guess 17.9); P2 S80 21.30 in [19.5, 25.5] (guess 21.6); P3 q(128 -> 256) 0.24 in [0.05, 0.40] (guess 0.21); P4 Delta q +0.01 in [-0.25, +0.15] (guess -0.02); P5 ladder 5.76 < 7.71 < 10.27 < 13.05 < 15.28 < 18.07; P6 q(64 -> 256) 0.23 in [0.12, 0.35] (guess 0.23); P7 flow realisation (D_eddy ratio 1.011, D_num/D_mol 0.0005, kappa 1.467, kicks 5.01, Lyapunov 0.130 > 0.1); P8 walled fraction 1.000 and 0 of 336 still walled at the end of the hold; P9 total 4,374 s, no cut, 336 of 336 runs.

**Decision (pre-registered, `:109`): UNDECIDED.** q(128 -> 256) = 0.24 [0.19-0.29]: the upper bound is above 0.2 (SATURATING not reached) and the lower bound is below 0.3 (GROWING not reached). **Second label (pre-registered, `:110`): CONSTANT-SLOPE EXCLUDED**, S50(1/256) 18.07 [17.66-18.52], upper bound below 19.61; the paired check, q upper bound 0.29 < 0.36, agrees (`:74`). **Local slope label (pre-registered, `:111`): UNRESOLVED**, Delta q = 0.01 [-0.08, 0.14]; h14d's FALLING label was not extended.

## Extrapolation to the real ratio 1e9: NOT a measurement (`:88-95`)

The ladder is 5 doublings of the ratio (8 to 256) against 27 from 8 to 1e9, 22 of them beyond the last measured cell. Anchored at the 1/256 point. The combined factor is f = 17 x S (the sphere factor 17 is taken from the task brief, not re-derived here). The Early Earth ocean share stays at 99.9 only while f is below about 2600, i.e. S below S* = 2600/17 = 153 (2400/17 = 141): the share is 99.925 at f = 2400 (`h13h14_extreme_locate.log:249`), 99.897 at 2600 (`:250`), 99.849 at 2800 (`:251`), 98.760 at 4500 (`:255`); then 94.208 at 1e4 (`h13h14_extreme.log:219`) and 82.133 at 3e4 (`:220`); the four share columns are main no lid / main final / x10 no lid / x10 final (section order `h13h14_extreme.log:54, 71, 88, 105`, table heading `:215`, base share 99.931 at f = 1 on `:216`). I quote the first column, which equals the second in these rows; the x10 columns are 0.2 to 0.4 points lower (99.478 at f = 3000 on `:218`, 93.861 at f = 1e4 on `:219`).

| assumption for the 22 doublings beyond 1/256 | S50 at 1e9 | f = 17 S | against f = 2600 | ratio at which S reaches S* (f = 2400 / 2600) | log line |
|---|---|---|---|---|---|
| constant from 1/256 | 18.07 | 307 | 8x below | never | `:95` |
| log law refitted over six points | 72 [69-76] | 1,230 | 2.1x below | 2.2e17 / 5.8e18 | `:94` |
| power law, q(128 -> 256) at its 2.5% bound (0.19) | 305 | 5,177 | 2.0x above | 1.6e7 / 2.5e7 | `:90` |
| power law, q(64 -> 256) point (0.23) | 638 | 1.08e4 | 4.2x above | 1.6e6 / 2.3e6 | `:93` |
| power law, q(128 -> 256) point (0.24) | 703 | 1.20e4 | 4.6x above | 1.3e6 / 1.8e6 | `:91` |
| power law, q(128 -> 256) at its 97.5% bound (0.29) | 1,567 | 2.66e4 | 10x above | 2.8e5 / 3.7e5 | `:92` |
| POST HOC arithmetic, not in the log: power law q 0.36, 0.40 (both excluded for the last doubling) | 4.3e3, 7.8e3 | 7.3e4, 1.3e5 | 28x, 51x above | not computed | scratch calculation |

The exponent that reaches S* at 1e9 from the 1/256 anchor is 0.141 (f = 2600) or 0.135 (f = 2400) (`:89`), the mean exponent over the 22 doublings. The measured q(128 -> 256) has a 95% lower bound of 0.19, above that, so any power law with a constant local exponent that is consistent with this cell passes the threshold; the log law passes below it only because its local exponent keeps falling (mean 0.091 over the 22 doublings, POST HOC arithmetic: ln(72.34 / 18.07) / ln(1e9 / 256)).

Share at the printed S50(1e9) values, read only from the table rows above and not interpolated: constant and log law (f = 307 and 1,230) lie in the flat part (99.9); the power-law rows f = 5,177, 1.20e4 and 2.66e4 lie between the printed rows for 4500 and 3e4, so between about 98.8 and 82.1 (the f = 5,177 row is just above the 4500 row, and 2.66e4 just below the 3e4 row).

### Which S scenarios the cell supports or excludes (against the Early Earth edge, f = 17 S about 2600)

| scenario for S at 1e9 | what the 1/256 cell says | Early Earth 99.9 share |
|---|---|---|
| power law held at 0.36-0.40 from 1/128 (h14c's S ~ 1e4) | excluded: scenarios C and D outside the 95% interval of S50 (-1.55, -2.10), second label EXCLUDED, paired q upper bound 0.29 < 0.36 | would have lost it (f ~ 7e4 to 1.3e5 at 1e9; share 82.1 at f = 3e4 and 61.2 at f = 1e5, `h13h14_extreme.log:220-221`) |
| power law held at the measured 0.19-0.29 | consistent, not excluded (it is the interval itself); the rows above give f from 5.2e3 to 2.7e4 | does not hold: below 99.9 (98.8 down to about 82) |
| log law (local exponent falling about 0.04 per doubling) | consistent, not excluded; fits all six points (rms 0.20), and its implied local exponents 0.54, 0.39, 0.31, 0.25, 0.22 (POST HOC arithmetic from a = -1.98, b = 3.59, rounding gives +-0.01) lie inside every measured interval | holds (f = 1,230) |
| constant S from 1/256 (growth stops now) | not excluded by this cell alone; the last doubling rose by 2.78 [2.16-3.39], so it requires a break after 1/256 | holds (f = 307) |

So the cell removes the two high power laws and narrows the measured exponent for the last doubling to 0.19-0.29 (was 0.13-0.28), but it does not separate the log law from a power law that has stopped falling, and those two sit on opposite sides of the Early Earth edge. The answer to "does the 99.9 share survive" still depends on whether the local exponent averages below about 0.14 over the next 22 doublings, which the data to 1/256 neither show nor exclude.

## POST HOC (not pre-registered; arithmetic or reading of numbers already printed)

- Between 1/128 and 1/256 the S50 gain per doubling returns to the h14d average (2.78 against the earlier increments 1.94-2.77); the exponent plateau reading (0.23 then 0.24) and the slowly-falling-exponent reading (log law: 0.254 then 0.216) predict the same increments within the interval of the last one (2.64 and 2.49 against 2.78 [2.16-3.39]; arithmetic, 15.28 x (2^0.23 - 1) and b ln 2 with b = 3.59). One more doubling at the same size would not separate them either: the log law (refit, POST HOC arithmetic) puts q(256 -> 512) at 0.19 against 0.23 for the constant-exponent power law, a difference of about 0.04 against a half-width of 0.05 to 0.075.
- S_half (the half-time of the mean wall measure; printed, not predicted) behaves differently from S50: 3.54, 3.93, 6.38, 7.93, 9.12, 9.31 (`:51-61`). The last step is +0.19 against +2.78 for S50, with overlapping 95% intervals (8.72-9.55 and 9.07-9.71). I draw no conclusion; it says that the median of the last-wall-gone time is growing while the bulk decay time has changed little between the last two cells, which would matter if the outcome rule depends on the bulk and not on the last walls. Not examined further.
- The cost estimate in `new_analyses_h14d_advect.md` ("a 1/256 cell of 336 runs about 12,000 cpu-s, roughly 70 min on 3 workers") held: 12,345 cpu-s and 73 min. By the same scaling a 1/512 cell of 336 runs would be roughly 2.3 times the cpu of 1/256, about 28,000 cpu-s and 2.7 hours on 3 workers (extrapolation of cost, not a timing).

## What this means

- The pre-registered answers: UNDECIDED on the exponent (0.24 [0.19-0.29]), CONSTANT-SLOPE EXCLUDED, local slope UNRESOLVED. The slope that fell from 0.4 to 0.23 over 1/8 to 1/128 is 0.24 over the next doubling; a further fall of the size a log law predicts (-0.04) is inside the noise of the change (-0.08 to 0.14). The two readings of the ladder, a log law with a falling exponent and a power law with the exponent stopped near 0.23, are both consistent with all six points' local intervals and with the new cell.
- The h14c extrapolation S50 ~ 1e4 at 1e9 is now excluded at the 95% level for a power law that continues with 0.36 or more; the range at 1e9 is 18 (constant) to 72 (log law) to 1.6e3 (power law at the upper bound of its exponent), and 3.0e2 to 1.6e3 for a power law with the measured interval. The real ratio was not reached.
- For the Early Earth edge (f = 17 S about 2600, S about 150), the power-law range lies above it by a factor 2 to 10 and the log law and the constant lie below it by 2 to 8. Whether the share of 99.9 holds is the same question as whether the exponent keeps falling, and the cell does not decide it.
- What would settle it: a cell at 1/512 costs about 2.7 hours and, with the same run count, would not by itself separate the two readings; a different method (a sharp-interface wall in a larger cell, or a 1D front tracked in the eddy flow) is the way to reach a much larger ratio. Neither was run.

## Suggested wording (NOT applied to `journals/ija_draft.md`)

These replace the `new_analyses_h14c_advect.md` and `new_analyses_h14d_advect.md` suggestions for the same places (apply only one of the three). `ija_draft.md` currently still has its original text at all four places (`:252`, `:479`, `:790`, `:187`). Log references are `h14e_advect.log`.

Item 14 (`:252`), replace with:

> 14. **[Partly resolved by E84; measured to 30 wall widths and to a molecular share of 1/256 by `h14c_advect.py`, `h14d_advect.py` and `h14e_advect.py`; not closed.]** **Eddy scale separation not shown.** Treating eddies as diffusion needs wall width √(D/k) > eddy size; in 2D, at 10 wall widths, the median wall-removal time was 5.8, 7.7, 10.3, 13.1, 15.3 and 18.1 times Eq. 17's time at molecular shares of 1/8, 1/16, 1/32, 1/64, 1/128 and 1/256 (`h14e_advect.log:51-62`). The exponent of the growth with the ratio of total to molecular diffusivity was 0.40 over 1/8 to 1/64, 0.23 (95% range 0.13 to 0.28) for 1/64 to 1/128 and 0.24 (0.19 to 0.29) for 1/128 to 1/256 (`:68-70`); a further fall could not be resolved (change 0.01, range -0.08 to 0.14, `:73`), and an exponent of 0.36 or more is excluded for the last doubling (`:110`). Open: the slowdown at the real ratio, about 10^9, which is bracketed only by extrapolation: 18 if growth stops at 1/256, 72 for a logarithmic law fitted to the six points and 1.6 × 10^3 if the power law continues at the upper bound of its exponent (`:88-95`); the smallest break-even factor of the sampled worlds is 100. Not tested: 1/128 and 1/256 at 30 wall widths, 3D, a sphere, more than one eddy per box.

Section 2.3, "Eddy diffusion" (`:479`), replace the sentences from "For eddies of 3.6 wall widths, however, ..." through "...once eddies exceed the wall width." (keep the last sentence, "I have not shown how many sampled inputs have eddies larger than the wall width.") with:

> For eddies of 3.6 wall widths, however, walls survived the transition in 29% of runs (1.6% with diffusion; `advect2d.log:51, 54`) and took longer to clear than Eq. 17 gives. In pre-registered follow-ups (`h14c_advect.py`, `h14d_advect.py`, `h14e_advect.py`, building on `h14b_advect.py`) the median clearing time was 3 times Eq. 17's at 3.6 wall widths and 6 times at 10 wall widths when molecular diffusion is 1/8 of the total (`h14c_advect.log:97`, `h14e_advect.log:51`); at 10 wall widths it was 8, 10, 13, 15 and 18 times at molecular shares of 1/16, 1/32, 1/64, 1/128 and 1/256 (`h14e_advect.log:53-61`), and 12 times at 30 wall widths with a share of 1/32 (`h14c_advect.log:75`). The factor grew with the ratio of total to molecular diffusivity as a power of 0.40 up to 1/64 (`h14e_advect.log:70`) and 0.23 to 0.24 over the next two doublings (`:68-69`); the real ratio is about 10^9 and was not reached. The diffusion description of the outcome is therefore supported; the diffusion description of wall motion is not, once eddies exceed the wall width, and how far it fails at the real ratio is not shown.

Section 4, eddies bullet (`:790`), replace the last sentence ("This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere.") with the paragraph below, and in the first sentence change "at least four times more slowly than diffusion at the same D_eff would" to "several times more slowly than diffusion at the same D_eff would (5.8 to 18 times in the median at 10 wall widths)":

> The slowdown was measured in 2D at eddy sizes of 3.6, 10 and 30 wall widths, with one eddy per box and molecular diffusion at 1/8 to 1/256 of the total (at 10 wall widths: 1/8, 1/16, 1/32, 1/64, 1/128, 1/256; `h14e_advect.log:51-62`). At 10 wall widths the median clearing time was 5.8 to 18.1 times Eq. 17's, growing as the ratio to a power of 0.40 up to 1/64 and 0.23 to 0.24 over the last two doublings. The real ratio is about 10^9 and was not reached; the growth is bracketed at 18 (no further growth) to 1.6 × 10^3 (power law at the upper bound of its exponent), against break-even factors of 100 to 10,000 for the sampled worlds (`h14_eddy.log:81-171`). The test was not done in 3D or on a sphere.

E86, end (`:187`), replace "H14: the slowdown is now bounded for the shares; the eddy-size question itself stays open." with:

> H14: the factor-10 case covers the measured median slowdown only at a molecular share of 1/8; at 1/32 to 1/256 the median is 10 to 18 times Eq. 17's (`h14e_advect.log:55-61`), and its growth with the ratio slowed from a power of 0.40 to 0.23-0.24 over the last two doublings (`:68-69`) without a further fall that the data could resolve (`:73`). The slowdown at the real ratio is not bounded below 100 (the upper extrapolation is 1.6 × 10^3).

## Limits

- 2D, one eddy per box, no stratification, no sphere. The l/w 10 ladder covers 5 doublings (8 to 256) against about 27 to the real 1e9; every statement about 1e9 is an extrapolation. The interval for the last doubling (0.19 to 0.29) reaches neither rule threshold (upper bound above 0.2, lower bound below 0.3); a further doubling at 336 runs would probably be about as wide (not computed).
- Only l/w 10 was extended. l/w 3.6 and 30 are h14c cells at 1/8-1/64 (l/w 30: 8 runs at 1/32 only).
- The molecular wall is 1 cell (tanh half-width 1.41) at h = 1 for all D_eff; the h = 0.5 partner of h14c agrees only at 1/32 and l/w 3.6. At 1/256 the box is 1005 cells per side.
- The 1/64, 1/128 and earlier cells are reloaded, not rerun, and their intervals are those of h14d; the pooled 1/64 cell (378 runs) is shared by q(64 -> 128) and q(64 -> 256), which is why those two intervals are not independent.
- The 1/256 cell is one flow family at one l/w; the molecular wall is resolved by 1 cell at every ratio, the same as at 1/128.
- Early Earth numbers use the combined factor f = 17 S from the task brief; I did not re-derive the 17.
- The log-law refit and the six-point least-squares fits are descriptions chosen after the ladder was known (the log law in h14d post hoc, the six-point refit in the h14e script before the run); they are two-parameter fits, not mechanisms, and a power law with a falling exponent is not excluded.
- The pre-registered guesses (not thresholds) were close: S50 17.9 against 18.07, S80 21.6 against 21.30, q 0.21 against 0.24, Delta q -0.02 against +0.01, q(64 -> 256) 0.23 against 0.23, wall time 4,100 s against 4,374 s.
