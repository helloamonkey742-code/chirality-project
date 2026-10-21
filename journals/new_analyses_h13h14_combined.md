# H13 x H14: world shares under the combined sphere x eddy slowdown of the sideways term (2026-10-21)

Pre-registration hash (sha256 of `h13h14_combined_predictions.txt`, asserted by the wrapper and printed on line 1 of the log):

`5385115994a2be924acfc1f2525e3b1642063b829b3319626197472474d81167`

Order of events (file times, PDT): predictions file last written 01:45:31, wrapper 01:45:49, results log and err written 01:46:27, exploratory log written 01:47:41.
The header of the predictions file says "~01:50 PDT"; that is a wrong estimate of mine, the file times above are the record. Nothing was run for this question before the hash (a 3-line self-check of the wrapper on made-up numbers aside).
After scoring, one extra exploratory run (not pre-registered, not scored) located the 1-point crossing.

## What was run

Wrapper `h13h14_combined.py` runs `h13_sphere_shares.py` unchanged (runpy, same argv) with `coarsen_slow.run` cached, so the script's own lines are identical to a direct run (f = 17 line reproduces `h13_sphere_shares.log:38`), and adds per-world losses in every column, the worlds at 1 point, rank swaps against f = 1 and the Early Earth lines (log lines 46 to 139).

```
cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project \
  --with numpy --with scipy --with numba --with matplotlib python h13h14_combined.py \
  f=1,17,34,100,136,300,1000 mix=h13_sphere_s16_f50.log > h13h14_combined.log 2> h13h14_combined.err
```

The shares are of all N = 200,000 draws (1 point = 2000 draws), seed 0, loss = share at f = 1 minus share at f. `h13h14_combined.err` holds only `exit=0`.
The results run took about 30 s, so nothing was near the 40 min limit. Exploratory crossing run: same wrapper, `f=1,38,40,42,45,50,60` to `h13h14_combined_cross.log`.

## The combined factor and its two assumptions

f = 17 (sphere, p90 of T1/tau_sw, `h13_sphere_score.log`) x S (eddy slowdown, `h14b_advect.log:58-80`).

| f | meaning |
|---|---|
| 17 | sphere alone |
| 34 | S50 about 2.0, l/w 3.6 (sets A to E: 1.4 to 2.2) |
| 100 | S50 5.8 to 6.3, l/w 10 (sets F, G; set H in a larger box 3.4, which would be 57) |
| 136 | S50 > 8, censored, l/w 30 (a lower bound on the factor) |
| 300, 1000 | stress; beyond the molecular-only ceiling S = 9 to 13 (f = 153 to 221) |

- **A1 (assumed, not tested):** the two factors multiply. The sphere factor was measured without flow; S in a flat periodic box with flow, relative to the no-flow time of the same box.
- **A2:** 17 is a p90 and S50 a median, so the product mixes statistics; S80 is censored at 5 to 16 in most sets.
- **A3:** the script slows every draw, though in nature draws with eddy size l <= w would not be slowed (about a third of Ganymede's P_pved draws at U = 1 m/s, `h14_eddy.log:95-98`), so the losses below are a worst case in that respect.

## Results (`h13h14_combined.log`)

Largest loss against f = 1, lines 37 to 44 (all Ganymede, column main no lid; Ganymede has no lid so main final is equal):

| f | 17 | 34 | 100 | 136 | 300 | 1000 |
|---|---|---|---|---|---|---|
| largest loss, points | 0.16 | 0.71 | 2.73 | 3.55 | 6.31 | 11.90 |

`mix` (line 44, 0.08 points) is the sphere factor alone and is not a combined result. The f = 1 check against the logged values is line 20: 0.05 points.

Per-world loss in points, main columns (lines 48 to 80; Ganymede line 51, Titan 52, Callisto 53, Europa 50, Dione 60):

| world | f = 34 | f = 100 | f = 136 | f = 300 | f = 1000 |
|---|---|---|---|---|---|
| Ganymede | 0.71 | 2.73 | 3.55 | 6.31 | 11.90 |
| Callisto | 0.50 | 2.28 | 3.06 | 5.56 | 10.90 |
| Titan | 0.56 | 2.09 | 2.75 | 4.95 | 9.41 |
| Europa (no lid / final) | 0.11 / 0.08 | 0.59 / 0.44 | 0.95 / 0.72 | 2.32 / 1.76 | 5.76 / 4.35 |
| Enceladus (no lid) | 0.11 | 0.52 | 0.71 | 1.31 | 2.69 |
| Triton | 0.09 | 0.45 | 0.74 | 1.91 | 4.96 |
| Pluto | 0.07 | 0.28 | 0.52 | 1.54 | 4.60 |
| Dione (no lid / final) | 0.01 / 0.01 | 0.04 / 0.02 | 0.06 / 0.03 | 0.15 / 0.10 | 1.11 / 0.70 |

- Worlds at 1.00 point or more (lines 117 to 122): none at f = 17 and 34 in any column; at f = 100 exactly Ganymede 2.73, Titan 2.09, Callisto 2.28 (main), and the same three in the x10 columns (2.41, 1.77, 1.99); at f = 136 the same three (Europa 0.95 just under); at f = 300 eight worlds in the no-lid column (Enceladus, Europa, Ganymede, Titan, Callisto, Triton, Pluto, Mimas); at f = 1000 thirteen.
- Rank swaps against f = 1 (lines 125 to 130): none at f = 17 and 34. At f = 100: main no lid and x10 no lid `Ganymede<Dione, Callisto<Dione` (Dione, which has a lid in the final columns, passes the two big moons when no lid is assumed: Ganymede 44.14 - 2.73 = 41.41, Callisto 43.47 - 2.28 = 41.19, Dione 41.95 - 0.04 = 41.91); both final columns unchanged. At f = 136 also `Callisto<Pluto` in main no lid and main final; more at 300 and 1000.
- Early Earth ocean (lines 133 to 139): printed 99.9 / 99.9 / 99.7 / 99.7 at every f; unrounded 99.931 / 99.931 / 99.670 / 99.670 at f = 1, loss at f = 1000 0.0025 / 0.0025 / 0.0100 / 0.0100.

Exploratory (NOT pre-registered, `h13h14_combined_cross.log`): the largest loss is 0.86 at f = 38, 0.93 at 40, 1.00 at 42 (lines 37 to 39 of that log), 1.11 at 45, 1.31 at 50 (Titan reaches 1.00 at 50), 1.64 at 60 (Callisto 1.27); worlds-at-1-point lines 115 to 120, swaps 123 to 128. So among the f tested the first world to reach 1 point (it prints 1.00, so it is on the line) does so at f of about 42, and none is at 1 point at 40 or below. The first rank swap among those f is Callisto<Dione at f = 50 in the x10 no-lid column.

## Scoring (as written, no band changed)

| | prediction | result | mark |
|---|---|---|---|
| P0 | f = 1 check <= 0.05; f = 17 line = 0.16 +/- 0.01 | 0.05 (line 20); 0.16 (line 38) | HELD |
| P1 | largest loss and world: 0.15 to 0.17 / 0.60 to 0.90 / 2.50 to 3.50 / 2.90 to 4.40 / 4.60 to 7.20 / 10.40 to 13.50, Ganymede, non-decreasing | 0.16 / 0.71 / 2.73 / 3.55 / 6.31 / 11.90, Ganymede at all six, column main no lid | HELD |
| P2 | at f = 100 main columns exactly {Ganymede, Titan, Callisto} at >= 1.00; nobody at f = 34 in any column | line 119 and line 118 | HELD |
| P3 | Early Earth printed values unchanged at all f; loss <= 0.01 | lines 133 to 139: unchanged, maximum loss 0.0100 | HELD, at the edge |
| P4 | no swaps at 17 and 34; at 100 both final columns unchanged, no-lid columns only Dione/Ganymede and Dione/Callisto, both in main no lid | lines 125 to 127 | HELD |
| P5 | first grid f with any world >= 1.00 (any column) is 100 | none at 34 (line 118), three at 100 (line 119) | HELD |

Notes on the marks, so the HELDs are not over-read:
- P1: the aside "Callisto is within 0.5 point" at f = 1000 was wrong; the gap is 1.00 (11.90 against 10.90). It was a caveat, not a band, so the mark stands. The found values were above my stated centrals at f >= 136 (centrals 3.3, 5.6, 11; found 3.55, 6.31, 11.90), all inside their bands (43%, 66% and 47% of the way up).
- P3: the band edge is hit exactly in the x10 columns at f = 1000 (20 draws = 0.0100); the parenthetical "expected 0.000 to 0.002" was wrong for those columns (0.010) and right for the main columns (0.0025).
- The x10 columns at f = 100 were predicted "near 1" for Titan and Callisto and are not scored; they are 1.77 and 1.99, so all three large moons are over 1 point in all four columns.
- P4's 85% on "no other pair flips" held with thin margins at f = 100 (Ganymede 41.41 against Callisto 41.19 in main final; Callisto against Pluto 40.8).
- 6 of 6 held; the floors in P1 and P2 came from `h14_eddy.log` counts, which are lower bounds on the loss, so these were not independent guesses.

## Decision rule (pre-registered, applied)

No world is at 1 point at f = 34 (R2 condition) but three are at f = 100: branch R2. "Never limits" cannot stand as written for eddies as large as 10 wall widths (f = 100); it holds for the one measured eddy size (l/w = 3.6, f = 34, largest loss 0.71) and, exploratorily, up to f about 40. R4 also applies: the effect scales with L^2/D_h, so it is not "equal across worlds", and the no-lid ranking changes at f = 100.

## Proposed wording (NOT APPLIED; `journals/ija_draft.md` untouched)

Line 789, replace "Sideways eddy mixing (Zhang et al., 2024) never limits horizontal coarsening." with:

> Sideways eddy mixing (Zhang et al., 2024) moves a share by a point or more only if it slows horizontal coarsening by more than about 40 times. The slowing factor used is the sphere factor of 17 times the eddy slowdown, assumed to multiply: eddies 3.6 wall widths wide slow wall removal about 2 times (a factor of 34), and no share moves by more than 0.71 point; eddies 10 wall widths wide slow it about 6 times (a factor of 100), and Ganymede, Callisto and Titan lose 2.7, 2.3 and 2.1 points, while every other world loses less than 1 point (`h13h14_combined.log`; the 40 is from an exploratory run, `h13h14_combined_cross.log`).

Line 790, keep everything up to "(`coarsen_slow.py`, `coarsen_slow.log`)", then replace ": in this model coarsening finishes well within the ocean's age, so a slower but finite clearing still finishes in time. This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere." with:

> . Sideways coarsening takes at most 7% of the ocean's age in any draw (`h13_sphere_shares.log:2-16`), so it decides a draw only when it is about 14 times slower. Slowing every draw by the sphere factor 17 times the eddy slowdown (an assumed product) moves no share by more than 0.71 point at a factor of 34, and by 2.7 points at most at 100 (Ganymede), where, with no lid, Dione passes Ganymede and Callisto; at 300, six to eight worlds (depending on the lid column) lose more than 1 point, up to 6.3, and the Early Earth ocean stays at 99.9 (99.7 with the times x10) up to 1000 (`h13h14_combined.log`). The eddy slowdown was measured at 3.6 and 10 wall widths and is only a lower bound at 30, in flat periodic boxes, not in 3D, and the sphere factor was measured without eddies.

Line 791 (R4, optional): after "The order of the worlds holds by construction against any effect that hits every world equally;" add "the sideways slowdown does not hit worlds equally (it scales with L^2/D_h), and it changes the no-lid order only at a factor of 100 or more;".

If the wording is adopted, the claim "never" in any other place (search the draft for "never limits") needs the same bound.

## Limits of this result

- The product f = 17 x S is an assumption; the log shows what happens at each f, not that the product is the right f.
- Applying f to every draw is a worst case (A3); applying it only to the draws with l > w would lower every loss above.
- Seed 0 only; the losses are exact counts on the same 200,000 draws, so they carry no sampling error for these draws but the 1-point threshold and the 1 decimal prints can sit either side of a rank by a few hundredths.
