# h14f_advect: does the eddy-slowdown exponent keep falling out to a molecular share of 1/4096? (2026-10-21)

Script `h14f_advect.py` (imports `h14e_advect.py`, `h14c_advect.py`, `h14b_advect.py`, `advect2d.py`; none edited). Predictions `h14f_predictions.txt`, written and hashed before any 1/4096 result existed (sha256 `d136fe694b8fe890833864590cb910d372424b770275f7519e169280f50c5c99`), committed with the script in cc5f0bb before the run. The script asserts that hash and the sha256 of `h14e_advect.py` and `h14e_advect_raw.pkl` (`h14f_advect.log:1`, again at `:41`). Run line:

```
UV_CACHE_DIR=$HOME/Library/Caches/uv caffeinate -i nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14f_advect.py > h14f_advect.log 2> h14f_advect.err
```

Started 2026-10-21 17:54:03 on mains power, exit 0 at 22:21:22. Log `h14f_advect.log` (111 lines; `h14f_advect.err` empty), raw results `h14f_advect_raw.pkl`. As with h14e, a `h14f_advect_raw.pkl.part` save file is left next to it (11 bytes larger; `cmp` first difference at byte 1,022,352; not used). No post hoc script was run; POST HOC items below are arithmetic on printed numbers. Nothing in `ija_draft.md`, `README.md` or `scorecard.py` was changed.

Output hashes (sha256): `h14f_advect.py` 45f2919e775be3b3fad37a058381f43907236c18bb2adf06ec74f8f50c9fe71f, `h14f_advect.log` aaa69c1db0b86ec4a3be270d8780b219bce27190483b3d525b0542dbb00faf94, `h14f_advect_raw.pkl` 702b00499f5c3545ca0b4447d200005d11dfa740173c05ef1f396fb17523883c.

## What I asked, in plain words

h14e left the exponent of the slowdown at 0.24 [0.19-0.29] for 1/128 -> 1/256, with two readings that both fit: a log law (exponent keeps falling) and a power law that stopped falling near 0.2. Those two sit on opposite sides of the early-Earth edge at the real ratio 1e9. One doubling more could not separate them, so h14f jumps four doublings to 1/4096 (cell V4096_10: 10 wall widths, box 4021 cells, 48 runs in 16 jobs of 3, seeds 3501-3516, hold 90 t17).

Pre-registered rule on q = q(256 -> 4096) with its 95% bootstrap interval: FALLING if the upper bound < 0.19 (the lower bound of the doubling before), HELD if the lower bound > 0.19, else UNDECIDED. Second label: BELOW EDGE SLOPE if the upper bound < q* (the exponent that reaches the early-Earth edge S* = 2600/17 at 1e9), else NOT BELOW. Third label: LOG LAW INSIDE if the six-point log law value 27.88 lies in the 95% interval of S50(1/4096). Stated in advance: FALLING about 80% under the log law; NOT BELOW over 99%.

## What ran

Flow at 1/4096 (`:44`): kappa 1.470, kicks to the end of the sweep 5.01, D_eddy 1.011 times nominal, numerical diffusion 0.0005 D_mol, Lyapunov 0.131 per kick; the same regime as every earlier rung. The kick was rewritten in numba so a worker fits in memory (~0.8 GB instead of ~2 GB of numpy temporaries); it is bitwise equal to the `h14b_advect.py` kick at N 64, 257 and 4021, and a tiny job pickles identically with either kick (`:3`). Self-checks and reloads of h14c, h14d and h14e reproduce 5.76, 7.71, 10.27, 13.05, 15.28 and 18.07 (`:4-6`). 4 workers, wall clock 16,019 s (4.4 h) against a 43,200 s guard, no job cut, 48 of 48 runs (`:106`); job times 3,417-4,586 s (`:25-40`).

## Results

S50 at l/w 10 (time since the end of the sweep / t17; Kaplan-Meier; bootstrap over runs, 2000x; `:51-64`):

| D_mol/D_eff | runs | S50 [95%] | S80 [95%] | S_half [95%] | mean domains at freeze-out | log line |
|---|---|---|---|---|---|---|
| 1/8 | 96 | 5.76 [4.94-6.44] | 8.86 [7.63-9.78] | 3.54 [3.07-4.03] | 3.43 | `:51-52` |
| 1/16 | 192 | 7.71 [7.34-8.40] | 12.02 [10.75-13.00] | 3.93 [3.49-4.36] | 4.93 | `:53-54` |
| 1/32 | 96 | 10.27 [9.48-11.32] | 13.62 [12.71-16.04] | 6.38 [5.75-7.01] | 7.91 | `:55-56` |
| 1/64 | 378 | 13.05 [12.74-13.80] | 16.48 [16.11-17.26] | 7.93 [7.69-8.34] | 12.50 | `:57-58` |
| 1/128 | 336 | 15.28 [14.93-15.79] | 18.79 [17.86-20.00] | 9.12 [8.72-9.55] | 22.77 | `:59-60` |
| 1/256 | 336 | 18.07 [17.66-18.52] | 21.30 [20.41-21.91] | 9.31 [9.07-9.71] | 42.74 | `:61-62` |
| **1/4096 (new)** | 48 | **28.85 [25.26-29.69]** | **31.90 [30.46-34.02]** | 14.03 [13.02-15.25] | 645.60 | `:63-64` |

Every run walled at freeze-out and cleared inside the hold; cleared by 16 / 24 / 32 / 48 t17: 0.00, 0.21, 0.81, 0.98 (`:64`, `:105`).

Local exponents (`:67-73`): 0.42, 0.41, 0.34, 0.23, 0.24, then **q(256 -> 4096) = 0.17 [0.12-0.18]**; Delta q against 128 -> 256 = **-0.07 [-0.15, -0.02]**, paired. Job-level bootstrap (not scored): 0.13-0.18 (`:74`). Per-doubling S50 increment over 256 -> 4096: 2.70 [1.84-2.92], against 1.94-2.78 before (`:76`).

Fits (printed, not scored; `:77-79`): the six-point log law from h14e predicts 27.88 at 1/4096, measured minus predicted +0.97. Log law refitted over seven ratios: S50 = -2.51 + 3.74 ln r, rms 0.27, largest residual +0.50 (at 1/8). One power law over seven ratios: q 0.25, rms 2.00, residual -4.42 at 1/4096.

Scenarios against the measured 95% interval (`:82-88`): inside, A log law 27.88 and B edge path 26.27 / 26.71; outside, C power law q 0.19 from 1/256 (30.60), D 0.24 (35.15), E 0.29 (40.38), F 0.33 (45.11).

Predictions (`:99-108`): **P1-P8 all held.** P1 S50 28.85 in [24, 42] (guess 29); P2 S80 31.90 in [28, 50] (guess 34); P3 q 0.17 in [0.10, 0.30] (guess 0.17); P4 Delta q -0.07 in [-0.16, +0.06] (guess -0.07); P5 seven-point ladder monotone; P6 flow checks; P7 walled fraction 1.000, 0 runs still walled at the end of the hold; P8 16,019 s, no cut, 48 runs.

**Decision (pre-registered, `:109`): FALLING.** q(256 -> 4096) = 0.17 [0.12-0.18], upper bound below 0.19. **Second label (`:110`): NOT BELOW**, upper bound 0.18 against q* 0.134 (f 2600; 0.128 at f 2400), as expected. **Third label (`:111`): LOG LAW INSIDE**, 27.88 inside 25.26-29.69.

## Extrapolation to 1e9: NOT a measurement (`:90-96`)

Nine doublings now measured (8 to 4096), 18 extrapolated (was 5 measured, 22 extrapolated). f = 17 S as before. Early-Earth share from `h13h14_extreme_locate.log`: 99.926 at f 2000 (`:248`), 99.925 at 2400 (`:249`), 99.897 at 2600 (`:250`), 99.419 at 3600 (`:253`), 98.922 at 4300 (`:254`), 98.760 at 4500 (`:255`), 94.996 at 9000 (`:256`); read from rows, not interpolated.

| assumption for the 18 doublings beyond 1/4096 | S50 at 1e9 | f = 17 S | early-Earth share | h14e value |
|---|---|---|---|---|
| constant from 1/4096 | 28.85 | 490 | 99.9 | 18 / 307 |
| log law refitted over seven | 74.9 [66-78] | 1,274 | 99.9 | 72 / 1,230 |
| power law at q 2.5% bound (0.122) | 131 | 2,231 | 99.9 (between the 2000 and 2400 rows) | 305 / 5,177 (q 0.19) |
| power law at q point (0.169) | 234 | 3,980 | between 99.4 and 98.9 (3600 and 4300 rows) | 703 / 1.20e4 (q 0.24) |
| power law at q 97.5% bound (0.182) | 276 | 4,697 | just below 98.8 (above the 4500 row) | 1,567 / 2.66e4 (q 0.29) |

The edge exponent q* is 0.134 (f 2600). The measured interval 0.12-0.18 now straddles it, where h14e's 0.19-0.29 lay wholly above. The high end of the range fell from 2.7e4 to 4.7e3.

## POST HOC (arithmetic on printed numbers)

- Under the seven-point log law the local exponent at 1/4096 is b / S50 = 3.74 / 28.85 = 0.13, already at q*, and its mean over the next 18 doublings is ln(74.91 / 28.85) / ln(1e9 / 4096) = 0.077. So the 99.9% holds if the exponent keeps falling as it has; it fails only if the exponent stops at about its present 0.13-0.18 for all 18 doublings.
- S_half, flat between 1/128 and 1/256 (+0.19, flagged in the h14e note), grew again: 9.31 -> 14.03, exponent 0.148 over 256 -> 4096. The h14e worry that the bulk decay had stopped while S50 kept growing does not hold up.
- The 95% interval of S50 is lopsided (25.26-29.69 around 28.85; 68% 26.84-29.28). Not examined; it is why the q interval is lopsided too (0.12-0.18 around 0.17).
- P1-P4 guesses: S50 29 against 28.85, S80 34 against 31.90, q 0.17 against 0.17, Delta q -0.07 against -0.07. Wall time 16,019 s against an estimate of about 23,000 s.

## What this means

- The exponent fell again, from 0.24 to 0.17, and the change excludes zero (-0.15 to -0.02). The log law written down before this run predicted the new point. A power law held at 0.19 or more from 1/256, the path that would have cost early Earth its 99.9%, is excluded at 1/4096.
- What is still open: the 18 doublings from 1/4096 to 1e9. The log law, a stop in growth, and a power law below about 0.13 keep early Earth at 99.9%; a power law that holds the measured 0.17 from here on would put it at about 98.9-99.4. The data favour the falling exponent; they cannot prove it continues.
- A further rung would only shift the same question four doublings on. A different method (sharp-interface or front-tracking at much larger ratio) would be needed to close it. Not planned.

## Suggested wording (NOT applied to `journals/ija_draft.md` or `README.md`)

Find -> replace. Log references are `h14f_advect.log` unless named.

Applied 2026-10-21 as E120, plus one change not listed below: the Section 4 bullet's opening "(5.8 to 18 times in the median" became "(5.8 to 29 times in the median". In Section 2.3 the third replacement keeps a comma: "(`h14e_advect.log:70`), 0.23 to 0.24 ...".

**Item 14 (`ija_draft.md:261`), whole item:**

> 14. **[Partly resolved by E84, E112 and E120; measured to 30 wall widths and to a molecular share of 1/4096 by `h14c_advect.py`, `h14d_advect.py`, `h14e_advect.py` and `h14f_advect.py`; not closed.]** **Eddy scale separation not shown.** Treating eddies as diffusion needs wall width √(D/k) > eddy size; in 2D, at 10 wall widths, the median wall-removal time was 5.8, 7.7, 10.3, 13.1, 15.3, 18.1 and 28.9 times Eq. 17's time at molecular shares of 1/8, 1/16, 1/32, 1/64, 1/128, 1/256 and 1/4096 (`h14f_advect.log:51-64`). The exponent of the growth with the ratio of total to molecular diffusivity fell from 0.40 over 1/8 to 1/64 (`h14e_advect.log:70`) to 0.23 and 0.24 over the next two doublings and 0.17 (95% range 0.12 to 0.18) over 1/256 to 1/4096, a fall of 0.07 (0.02 to 0.15) from the doubling before (`h14f_advect.log:70-73`). A logarithmic law fitted to the first six points before that run predicted 27.9 at 1/4096, inside the measured range of 25.3 to 29.7 (`:77`, `:111`). Open: the slowdown at the real ratio, about 10⁹, which is still bracketed only by extrapolation over 18 more doublings: 29 if growth stops at 1/4096, 75 for the logarithmic law fitted to all seven points and 131 to 276 if a power law continues with the exponent measured over the last four doublings (`:90-96`). Times the sphere factor of 17 (an assumed product) that is about 490 to 4.7 × 10³; early Earth keeps 99.9% up to about 2600 (Section 4). Not tested: shares below 1/32 at 30 wall widths, 3D, a sphere with eddies, more than one eddy per box.

**Section 2.3, Eddy diffusion (`:489`):**

- "(`h14c_advect.py`, `h14d_advect.py`, `h14e_advect.py`, building on" -> "(`h14c_advect.py`, `h14d_advect.py`, `h14e_advect.py`, `h14f_advect.py`, building on"
- "8, 10, 13, 15 and 18 times at molecular shares of 1/16, 1/32, 1/64, 1/128 and 1/256 (`h14e_advect.log:53-61`)" -> "8, 10, 13, 15, 18 and 29 times at molecular shares of 1/16, 1/32, 1/64, 1/128, 1/256 and 1/4096 (`h14f_advect.log:53-63`)"
- "and 0.23 to 0.24 over the next two doublings (`:68-69`); the real ratio" -> "0.23 to 0.24 over the next two doublings (`:68-69`) and 0.17 over the four doublings to 1/4096 (`h14f_advect.log:73`); the real ratio"

**Section 4, eddies bullet (`:800`):**

- "molecular diffusion at 1/8 to 1/256 of the total (`h14e_advect.log:51-62`)" -> "molecular diffusion at 1/8 to 1/4096 of the total (`h14f_advect.log:51-64`)"
- "At 10 wall widths the median clearing time was 5.8 to 18.1 times Eq. 17's, growing with the ratio of total to molecular diffusivity as a power of 0.40 up to 1/64 and 0.23 to 0.24 over the last two doublings." -> "At 10 wall widths the median clearing time was 5.8 to 28.9 times Eq. 17's, growing with the ratio of total to molecular diffusivity as a power of 0.40 up to 1/64, 0.23 to 0.24 over the next two doublings and 0.17 (0.12 to 0.18) over the four doublings to 1/4096; a logarithmic law fitted before the 1/4096 run predicted that point (27.9 against 28.9, 95% range 25.3 to 29.7; `h14f_advect.log:73`, `:111`)."
- "extrapolated, the slowdown there is 18 if growth stops, 72 for a logarithmic law and 300 to 1.6 × 10³ for a power law with the measured exponent (`h14e_advect.log:88-95`)" -> "extrapolated over the remaining 18 doublings, the slowdown there is 29 if growth stops, 75 for a logarithmic law and 131 to 276 for a power law with the exponent measured over the last four doublings (`h14f_advect.log:90-96`)"
- "the extrapolated range is about 300 to 2.7 × 10⁴" -> "the extrapolated range is about 490 to 4.7 × 10³"
- "Its 99.9% therefore holds at the real ratio under the logarithmic law (a factor of 1230) and not if the power law continues with an exponent of 0.19 or more (5.2 × 10³ or more); the data do not decide between them." -> "Its 99.9% therefore holds at the real ratio under the logarithmic law (a factor of 1270), if growth stops (490), and under a power law with an exponent below about 0.13 (`h14f_advect.log:91`); a power law that holds the measured 0.17 from 1/4096 on gives 4.0 × 10³ and a share between 98.9 and 99.4% (`h13h14_extreme_locate.log:253-254`). The exponent has fallen over nine doublings, as the logarithmic law requires; that it keeps falling over the next 18 is not shown."

**E86, end (`:188`):**

- "at 1/32 to 1/256 the median is 10 to 18 times Eq. 17's (`h14e_advect.log:55-61`), and its growth with the ratio slowed from a power of 0.40 to 0.23–0.24 over the last two doublings (`:68-69`) without a further fall that the data could resolve (`:73`). The slowdown at the real ratio, about 10⁹, is bracketed only by extrapolation, 18 to 1.6 × 10³ (`:88-95`); see E112." -> "at 1/32 to 1/4096 the median is 10 to 29 times Eq. 17's (`h14f_advect.log:55-63`), and its growth with the ratio slowed from a power of 0.40 to 0.23–0.24 and then 0.17 over the four doublings to 1/4096 (`:70-73`). The slowdown at the real ratio, about 10⁹, is bracketed only by extrapolation, 29 to 276 (`:90-96`); see E112 and E120."

**New E120 (after E119, `:221`):**

> 120. **Eddy slowdown at a molecular share of 1/4096 (2026-10-21; `h14f_advect.py`, `h14f_advect.log`; predictions `h14f_predictions.txt` committed first in cc5f0bb; `journals/new_analyses_h14f_advect.md`).** One cell of 48 runs at 10 wall widths and a molecular share of 1/4096, four doublings past 1/256. Median clearing time 28.85 times Eq. 17's (95% 25.26–29.69; `h14f_advect.log:63`). The exponent over 1/256 to 1/4096 is 0.17 (0.12–0.18), below the 0.19 threshold, so the pre-registered label is FALLING (`:109`); the change from the doubling before is −0.07 (−0.15 to −0.02; `:73`). The six-point logarithmic law predicted 27.88, inside the measured range (`:111`). Second label NOT BELOW (`:110`): the measured exponent, if held to 10⁹, would still pass the early-Earth edge. All eight predictions held (`:99-108`). The kick was rewritten in numba to fit memory and is bitwise equal to the `h14b_advect.py` kick (`:3`). Item 14, the Section 2.3 eddy paragraph, the Section 4 eddies bullet, the end of E86 and README Part 6 updated: the extrapolated slowdown at 10⁹ is now 29 (no further growth), 75 (logarithmic) and 131–276 (power law with the measured exponent), 490 to 4.7 × 10³ with the sphere factor (was 300 to 2.7 × 10⁴); early Earth's 99.9% holds under the logarithmic law and is not settled for a power law that holds the measured exponent. The real ratio was not reached.

**README Part 6 (`README.md:710-712`):**

- "The measured slowdown times the sphere's 17 gives 34–307×; at the real eddy ratio it is extrapolated to about 300–27,000×, so the early-Earth 99.9% is not settled (`h14e_advect.log:88-95`)." -> "The measured slowdown times the sphere's 17 gives 34–490×; at the real eddy ratio it is extrapolated to about 490–4,700×. The early-Earth 99.9% holds under the logarithmic law, which predicted the 1/4096 point, but not if the exponent stops falling at its measured 0.17 (`h14f_advect.log:90-96`, `:111`)."

## Limits

- 2D, one eddy per box, no stratification, no sphere; only l/w 10 extended. Every statement about 1e9 is an extrapolation over 18 doublings.
- 48 runs at 1/4096 against 336 at 1/256; the interval is set by the bootstrap over runs (job-level bootstrap 0.13-0.18, close).
- Mean domains at freeze-out grew 15x from 1/256 to 1/4096 (42.7 -> 645.6), with the box 4x wider. Part of S50's growth may be the time to clear more walls rather than slower walls; that split was not examined at any rung.
- The log law and the power-law fits are descriptions, not mechanisms; the log law was chosen post hoc in h14d and fixed before this run.
- f = 17 S keeps the assumed sphere factor 17; not re-derived.
