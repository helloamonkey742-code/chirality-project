# H13 x H14: world shares when the sideways-coarsening slowdown f is far beyond 1000, up to f -> infinity (2026-10-21)

Pre-registration hash (sha256 of `h13h14_extreme_predictions.txt`, asserted by the wrapper and printed on line 1 of the log):

`97a76e595ea94f836c9915d2de0761ed5ea1e673f35cc7fd17c79cff57fa203f`

Order of events (file times, PDT): analytic forecast `h13h14_extreme_forecast.py` 03:02:27 and `.out` 03:03:36 (numpy only, no sweep); predictions file 03:12:01; wrapper last written 03:14:23; results log and err 03:15:41; `h13h14_extreme_repro.out` 03:15:45; exploratory log 03:18:19.
Nothing at f > 1000 was run before the hash. The wrapper was first written after the hash and edited twice before the results run: (1) its self-check on made-up numbers expected counts (3, 2, 1) where the right value is (3, 3, 3), an arithmetic slip in the test, not in the report code; (2) an unused variable was deleted. The self-check never touched a sweep.

## What f -> infinity means (from the code, before any run)

`U.A = A0 / sqrt(f)` multiplies the sideways time (L/A)^2/D_h by f in every draw. That time is a necessary condition for success in all three sweeps (`sweep2.py:127,129,137`; `sweep3.py:122-123`; `multiwall.py:421-422`), so a draw passes only if f x ratio <= 1 (in the timely branch 1 - t_rise/t). Above a finite f per world every draw fails and the share is exactly 0. The limit is therefore degenerate: all 15 worlds, all four columns, share 0.000, a 15-way tie, and the Early Earth ocean goes from 99.93 to 0.00. It is run as f = 1e9 and 1e12 (np.inf would make U.A = 0) and the two are compared cell by cell. The informative result is the f at which each world collapses (about 9e7 to 1.8e8, table below), not the limit value.

## What was run

`h13h14_extreme.py` is a copy of the `h13h14_combined.py` wrapper with the new hash, strict tie-aware ranks, counts of worlds losing >= 1 / 5 / 10 points, a monotone check, a limit check and a `check` mode. It runs `h13_sphere_shares.py` unchanged (runpy, same argv) with `coarsen_slow.run` cached. No existing file was edited.

```
UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project \
  --with numpy --with scipy --with numba --with matplotlib python h13h14_extreme.py \
  f=1,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000,100000000,1000000000,1000000000000 \
  mix=h13_sphere_s16_f50.log > h13h14_extreme.log 2> h13h14_extreme.err
```

61 s, `h13h14_extreme.err` holds `exit=0`. Shares are of all N = 200,000 draws (1 point = 2000 draws), seed 0, loss = share at f = 1 minus share at f. One process at a time.
Regression: `python h13h14_extreme.py check` > `h13h14_extreme_repro.out`: 200 of 200 values at f = 1 and f = 1000 are identical to `h13h14_combined.log` (60 base shares, 60 losses, 45 flip entries, 30 P_pved entries, the f = 1 check line, both largest-loss lines, both Early Earth lines).
Exploratory run (NOT pre-registered, not scored): `f=1,1000,2000,2400,2600,2800,3200,3600,4300,4500,9000,9500,16000,17000,1.2e8,1.5e8,1.8e8,2e8,3e8` to `h13h14_extreme_locate.log` (88 s).

## Results (`h13h14_extreme.log`; line numbers are that file)

Largest loss against f = 1, the script's own lines 38-49 (column as printed):

| f | 1000 | 3000 | 1e4 | 3e4 | 1e5 | 3e5 | 1e6 | 3e6 | 1e7 | 1e8 | 1e9 and 1e12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| world | Ganymede | Ganymede | Ganymede | Ganymede | Early Earth | Early Earth | Early Earth | Early Earth | Early Earth | Early Earth | Early Earth |
| points | 11.90 | 17.14 | 22.91 | 28.15 | 38.97 | 62.47 | 82.94 | 94.61 | 99.84 | 99.93 | 99.93 |
| column | main | main | main | main | **x10 no lid** | **x10 no lid** | main | main | main | main | main |

In the main columns Early Earth is 38.75 at 1e5 (line 127) and 62.39 at 3e5 (line 128); its x10 loss is larger there because the x10 base is 99.67, not 99.93. Second place in the main columns: Callisto 16.08 / 21.88 / 27.06 at 3000 / 1e4 / 3e4, then Ganymede 33.88 / 38.55 / 42.18 at 1e5 / 3e5 / 1e6, Titania 44.56 at 3e6, Umbriel 49.23 at 1e7 and 53.67 to 53.68 from 1e8 (lines 127-134).

Main no lid share, points (lines 55-69 minus the base; main final is equal except Enceladus 0.00, Europa 30.07 -> 22.28 / 18.20 / 14.63 / 10.61 / 3.66, Dione 26.20 -> 24.03 / 21.24 / 17.98 / 14.45 / 7.91, Miranda 8.97 -> 8.80 / 7.85 / 6.67 / 5.46 / 3.26, Mimas 5.56 -> 3.59 / 2.84 / 2.19 / 1.52 / 0.33):

| world | f = 1 | 3000 | 1e4 | 3e4 | 1e5 | 1e6 |
|---|---|---|---|---|---|---|
| Early Earth ocean | 99.93 | 99.76 | 94.21 | 82.13 | 61.18 | 16.99 |
| Umbriel | 53.68 | 48.87 | 43.00 | 36.38 | 29.11 | 15.87 |
| Oberon | 51.38 | 44.59 | 38.11 | 31.82 | 25.08 | 12.45 |
| Titania | 51.28 | 44.16 | 37.66 | 31.43 | 24.66 | 12.09 |
| Ariel | 49.59 | 45.19 | 39.79 | 33.69 | 27.04 | 14.75 |
| Ganymede | 44.14 | 27.00 | 21.23 | 15.99 | 10.26 | 1.96 |
| Callisto | 43.47 | 27.39 | 21.59 | 16.41 | 10.84 | 2.31 |
| Dione | 41.95 | 38.46 | 34.01 | 28.82 | 23.16 | 12.69 |
| Pluto | 41.08 | 31.97 | 26.48 | 21.58 | 16.12 | 6.35 |
| Europa | 39.61 | 29.34 | 23.98 | 19.26 | 14.03 | 4.83 |
| Triton | 37.01 | 27.93 | 23.02 | 18.53 | 13.66 | 4.97 |
| Titan | 35.80 | 22.05 | 17.37 | 13.05 | 8.47 | 1.74 |
| Miranda | 27.09 | 26.58 | 23.69 | 20.12 | 16.41 | 9.73 |
| Enceladus | 26.36 | 21.80 | 19.07 | 15.87 | 12.45 | 6.47 |
| Mimas | 17.30 | 11.17 | 8.86 | 6.81 | 4.70 | 1.06 |

- Early Earth ocean (lines 216-228): main 99.760 / 94.208 / 82.133 / 61.181 / 37.537 / 16.988 / 5.321 / 0.096 at f = 3000 / 1e4 / 3e4 / 1e5 / 3e5 / 1e6 / 3e6 / 1e7 and 0.000 at 1e8, 1e9, 1e12; x10 columns 99.478 / 93.861 / 81.683 / 60.701 / 37.203 / 16.864 / 5.288 / 0.094. Main no lid = main final and x10 no lid = x10 final at every f. Rank 1 in all four columns through f = 1e6; 7 / 5 / 6 / 5 at 3e6 and 11 / 10 / 11 / 10 at 1e7.
- Worlds losing >= 1 / 5 / 10 points (lines 137-148, main no lid | main final): 3000: 13 / 9 / 4 | 12 / 8 / 3; 1e4: 15 / 14 / 9 | 14 / 11 / 9; 3e4: 15 / 15 / 14 | 14 / 12 / 11; 1e5: 15 / 15 / 15 | 14 / 12 / 12; 1e6 to 1e12: 15 / 15 / 15 | 14 / 14 / 12.
- Rank of the three large moons, main no lid (lines 155-157): Ganymede 6 -> 10 / 11 / 12 / 12 / 13 at 1000 / 3000 / 1e4 / 3e4 / 1e5; Callisto 7 -> 9 / 10 / 11 / 11 / 12; Titan 12 -> 13 / 13 / 14 / 14 / 14. Dione, Pluto and Europa pass Ganymede and Callisto by 1e4 (ranks 6, 7, 8 against 12 and 11).
- Five highest worlds (lines 201-213, same set in both main columns): {Early Earth, Umbriel, Ariel, Oberon, Titania} from 3000 through 3e5; at 1e6 Dione replaces Titania in main no lid; at 3e6 Early Earth drops out of the top five in main no lid. Reversed pairs against f = 1 (not scored; lines 187-196): main no lid 12 / 16 / 17 / 20 / 22 / 26 / 32 at 3000 / 1e4 / 3e4 / 1e5 / 3e5 / 1e6 / 3e6.
- Approach to the limit (lines 231-238): at f = 1e7 four of 15 worlds are exactly 0 (the rest up to Umbriel 4.44, Early Earth 0.096); at 1e8 only Enceladus 0.094, Dione 0.021, Ariel 0.0135, Umbriel 0.008 are non-zero; at 1e9 and 1e12 all 60 cells are exactly 0.0, identical in the two runs (lines 240-241: monotone yes; limit identical yes, all exactly 0.0 yes).

Collapse f per world (analytic from the forecast, checked in the exploratory log, `h13h14_extreme_locate.log:267-276` and flip block line 3): Miranda 8.7e7 (zero at 1e8, log line 233 of the main log), Umbriel 1.1e8 and Ariel 1.15e8 (zero at 1.2e8), Dione 1.2e8 (0.001 at 1.2e8, zero at 1.5e8), Enceladus 1.8e8 (0.008 at 1.5e8, zero at 1.8e8, 2e8, 3e8).

Exploratory, Early Earth main share (`h13h14_extreme_locate.log:246-259`): 99.925 at 2400, 99.897 at 2600 (prints 99.9), 99.849 at 2800 (prints 99.8), 99.419 at 3600, 98.922 at 4300 (below 99), 94.996 at 9000 (below 95), 89.814 at 16000 (below 90). The printed 99.9 therefore survives to f about 2600, 11.8 times the molecular-only ceiling of 221, and goes between 2600 and 2800. The analytic forecast (`h13h14_extreme_forecast.out`, lam = 1) matched the main-column Early Earth share to within 0.05 points at every f from 3000 to 1e7 (found minus forecast: -0.023, -0.017, -0.005, -0.022, -0.040 at 3000 to 1e6; +0.045 at 3e6; +0.011 at 1e7).

## Scoring (as written, no band changed)

| | prediction | result | mark |
|---|---|---|---|
| P0 | every f = 1 and f = 1000 value equals `h13h14_combined.log` | 200 of 200 identical (`h13h14_extreme_repro.out`); lines 20, 38, 216-217 | HELD |
| P1 | (a) non-increasing in 60 cells; (b) 1e9 = 1e12, all exactly 0.0; (c) 15-way tie in every column, 0 of 105 ordered; (d) largest-loss line at 1e9 and 1e12 is Early Earth 99.93, a main column | line 240; line 241; lines 197-198 (ordered 0, tied 105 in all four columns) and ranks all 1 (lines 153-184, last two columns); lines 48-49 (99.93, main no lid) | HELD |
| P2 | largest loss and world per f: 3000 Ganymede 15.5-19.0; 1e4 Ganymede 21.5-25.5; 3e4 Ganymede 26.8-30.5; 1e5 Early Earth 38.2-39.2; 1e6 Early Earth 82.3-83.5; column main no lid or main final; takeover between 3e4 and 1e5 | 17.14 (line 39), 22.91 (40), 28.15 (41), 82.94 (44), all main no lid; 1e5 is Early Earth 38.97, in band, but the line names column **x10 no lid** (line 42); takeover between lines 41 and 42 | FAILED as written (4 of 5 rows held; the 1e5 row fails on the column clause) |
| P3 | Early Earth share, main / x10: 3000 99.70-99.86 / 99.40-99.65; 1e4 93.7-94.7 / 93.2-94.6; 3e4 81.6-82.7 / 81.2-82.5; 1e5 60.6-61.8 / 60.1-61.7; 1e6 16.5-17.5 / 16.4-17.4; limit 0.000 exact | 99.760 / 99.478; 94.208 / 93.861; 82.133 / 81.683; 61.181 / 60.701; 16.988 / 16.864; 0.000 in all four at 1e9 and 1e12 (lines 218-228); main no lid = main final | HELD (6 of 6 rows) |
| P4 | counts >= 1 / 5 / 10, main no lid and main final, per f (bands in the predictions file) | all 12 row-halves inside the band (lines 138-148, see the counts above) | HELD |
| P5 | (a) Early Earth rank 1 in both main columns at 3000, 1e4, 3e4, 1e5, 3e5 and rank 1-3 at 1e6; (b) top-five set at 3000-1e5 is {Early Earth, Umbriel, Ariel, Oberon, Titania}; (c) at 1e4 main no lid Ganymede and Callisto strictly below Dione, Pluto and Europa | (a) rank 1 at all six f (lines 167, 184); (b) lines 203-206; (c) ranks 12 and 11 against 6, 7 and 8 (lines 155-166) | HELD |
| P6 | at 1e8 exactly Enceladus, Dione, Ariel, Umbriel non-zero in main no lid; Enceladus largest, 0.03-0.15; the other eleven exactly 0.0 | line 233: those four, Enceladus 0.094, eleven exactly 0.0 | HELD |

6 of 7 items HELD; 32 of 33 scored parts held (P2 at 1e5 is the one miss).

Notes on the marks, so the HELDs are not over-read:
- P2 at 1e5: the value (38.97) and the world are inside the prediction; my parenthetical "equal for the worlds named" was wrong for Early Earth, whose x10 loss exceeds the main loss at 1e5 and 3e5 (38.97 against 38.75, 62.47 against 62.39) and which the script then names. By the rule that every part must hold, it fails. The main-column value 38.750 (line 127) is also in band.
- P3: the x10 shares sit below my centrals at every f (99.478 against 99.52 at 3000; 16.864 against 17.0 at 1e6) but inside the bands; the main bands were met with the analytic model, no fitting.
- P6: the 75% confidence was justified in part. About 550 draws pass the sideways test for Enceladus at 1e8 (flip share 99.725%, line 3), and only 188 of them pass the other tests (0.094%).
- P4 and P5 were a joint forecast with P3 and P2 from one analytic model; they are not independent checks.
- Lambda (the share among the draws that do not flip) was extrapolated beyond f = 1000 and the forecast Ganymede-over-Callisto order at 3000 to 3e4 held with margins of 1.06, 1.03 and 1.09 points (lines 123-126), which is not large.

## Decision rule (pre-registered, applied)

- **X1 (Early Earth headline):** P3 holds at 3000 and 1e4, so the headline 99.9 is kept only with a bound. As pre-registered the bound is the logged 1000 (4.5 times the ceiling 221); the exploratory run moves it to about 2600 (prints 99.9 up to 2600, 99.8 at 2800 and above; below 99 by 4300, 95 by 9000, 90 by 16000). Beyond that Early Earth is 99.8 at 3000, 94.2 at 1e4, 82.1 at 3e4, 61.2 at 1e5, 17.0 at 1e6 and 0 at 1.8e8 and above.
- **X2 (Early Earth first):** P5a holds, so "Early Earth has the highest share" is kept with the bound "for f up to 3e5 at least"; the log shows it through 1e6 and not at 3e6 (rank 7 in main no lid, rank 5 in main final). At the limit it is a 15-way tie, which is not a ranking.
- **X3 (`ija_draft.md:789`, "never limits"):** P2 holds at 3000 (Ganymede 17.14), so the sentence cannot stand beyond f about 40 (`h13h14_combined_cross.log`). Give the measured and ceiling range and, as extrapolation, Ganymede 17.1 / 22.9 / 28.1 at 3000 / 1e4 / 3e4 and Early Earth 38.8 and 82.9 at 1e5 and 1e6 (main).
- **X4 (`ija_draft.md:790`, "finishes in time"):** bounded world by world by the margin 1/worst of the flip block (lines 3-17): Ganymede 14.4, Titan 15.0, Callisto 17.0, Enceladus 18.4, Europa 45.8, Triton 53.7, Mimas 50.2, Pluto 68.7, Titania 156, Oberon 168, Umbriel 286, Ariel 295, Dione 308, Miranda 1755, Early Earth 2400. In the limit no draw of any world finishes: every share is 0 above 1.8e8 (a hard gate set by the declared D_h and age ranges, not a result about real oceans).
- **X5 (`ija_draft.md:791`, "holds by construction against any effect that hits every world equally"):** kept as a conditional; the sideways slowdown scales with L^2/D_h and does not hit worlds equally. P5b holds: the top five do not change through 3e5, while Ganymede, Callisto and Titan fall from ranks 6, 7 and 12 to 12, 11 and 14 by 1e4 (no-lid column). In the limit the order is lost entirely (P1c).
- **X6:** no run failed and nothing was NaN, so every f was scored.

## Proposed wording (NOT APPLIED; `journals/ija_draft.md`, `README.md` and `scorecard.py` untouched)

These extend, beyond f = 1000, the wording proposed in `journals/new_analyses_h13h14_combined.md`; the parts up to f = 1000 are unchanged from there.

Line 789, replace "Sideways eddy mixing (Zhang et al., 2024) never limits horizontal coarsening." with:

> Sideways eddy mixing (Zhang et al., 2024) moves a share by a point or more only if it slows horizontal coarsening by more than about 40 times (`h13h14_combined_cross.log`). The slowing factor used is the sphere factor of 17 times the eddy slowdown, assumed to multiply: the measured eddy sizes give 34 to 100 and the molecular-only ceiling is 221, and no share moves by more than 2.7 points up to 100 (`h13h14_combined.log`). Above that the factor is extrapolated from a power law that was not measured there. Slowing every draw by 3000 takes 17.1 points from Ganymede and 16.1 from Callisto, by 10,000 takes 22.9 and 21.9, and by 100,000 takes 38.8 from the Early Earth ocean (`h13h14_extreme.log:39-42,127`).

Line 790, keep everything up to "(`coarsen_slow.py`, `coarsen_slow.log`)", then replace ": in this model coarsening finishes well within the ocean's age, so a slower but finite clearing still finishes in time. This was tested at one eddy size above the wall width, in a box of 7.2 coherence lengths, and not in 3D or on a sphere." with:

> . Sideways coarsening takes at most 7% of the ocean's age in any draw (`h13_sphere_shares.log:2-16`), so a slowdown decides a draw only when it exceeds that draw's margin, as small as 14 for Ganymede and as large as 2400 for the Early Earth ocean (`h13h14_extreme.log:3-17`). Slowing every draw by a factor up to 1000 moves the Early Earth ocean by 0.0025 points and keeps it at 99.9; it prints 99.8 from about 2800, falls below 95 by 9000, and is 61 at 100,000 and 17 at 1,000,000 (`h13h14_extreme.log:216-223`; `h13h14_extreme_locate.log:246-259`). At 1.8e8 and above no draw of any world finishes in time and every share is 0. The eddy slowdown was measured at 3.6 and 10 wall widths and is only a lower bound at 30, in flat periodic boxes, not in 3D or on a sphere, and the sphere factor was measured without eddies.

Line 791 (X5), after "The order of the worlds holds by construction against any effect that hits every world equally;" add: "the sideways slowdown does not hit worlds equally (it scales with L^2/D_h): at a factor of 10,000 Ganymede, Callisto and Titan fall from ranks 6, 7 and 12 to 12, 11 and 14 while the top five are unchanged as a set, and the Early Earth ocean stays first up to 1,000,000 (`h13h14_extreme.log:155-157,167`);".

Line 666 (`26.4 / 39.6 / 99.9%`) is the f = 1 result; if the headline is bounded, add after it: "for a sideways slowing factor up to about 1000 (99.8 at 3000, 94 at 10,000, 61 at 100,000)".

If adopted, search the draft and `README.md:683` ("and never limits healing here") for the same unbounded claim.

## Limits of this result

- f is applied to every draw, a worst case (in nature draws with eddies no larger than the wall width are not slowed); every loss above is an upper bound in that respect.
- The product 17 x S is an assumption, and f above 221 is extrapolation of a fitted power law, not a measurement. The log shows what happens at each f, not which f is right.
- The f -> infinity result (every share 0) follows from the code's structure and the declared D_h and age ranges. It is a hard gate by construction, not a finding about real oceans; the collapse f per world (9e7 to 1.8e8) is the content.
- Seed 0 only. The shares are exact counts on the same 200,000 draws, so differences between f carry no sampling error for these draws, but ranks between worlds within a few hundredths of a point can sit either side.
- The exploratory run is not pre-registered and is not scored; its onset values (2600 to 2800, 4300, 9000, 16000) agree with the analytic numbers written in the predictions file as "not scored".
