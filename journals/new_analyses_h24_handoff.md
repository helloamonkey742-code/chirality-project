# h24_handoff: can a handoff rule predict the filling hand, and what do the "no layer holds" draws do? (2026-10-21)

**Short answer.** Yes, for the range tested. A four-coefficient logistic handoff rule (V3), fitted on the 600 old radial runs and frozen before any new simulation, predicted the final hand in **96.5 ± 0.8 %** of 570 fresh runs and in **86.5 ± 3.4 %** of the 104 runs where it handed the decision up (`h24_handoff_test.log:32`). The pre-registered bar was 95 % overall and 85 % on fired runs, so the criterion is met on the point estimates. The fired margin is thin (the bar sits 0.44 SE below the estimate). The old "93 % ceiling" was a ceiling of rules that use the deciding domain's thickness alone: even with the threshold tuned on the fresh runs, such a rule reaches 94.7 % (`h24_handoff_test.log:61`). Adding the drift of the domain's rising wall and the thickness of the domain above gains 1.8 points over that oracle. Seed noise is not the remainder: 177 configurations rerun with two fresh noise seeds each flipped the final hand 0 times (`h24_handoff_test.log:18-19`).

For the world shares: Early Earth ocean P_pved is **98.2 %** if the "no layer thick enough" draws are counted as no hand and **99.9 %** if the deepest layer decides them (`h24_handoff_worlds.log:12`). The previous values were 94.3 (three-width cut) and 97.6 (one-width cut). Moons move by +1.6 to +5.0 points between the two treatments (`h24_handoff_worlds.log:102-116`).

**Not closed:** stacks where every domain is under one wall width (68 runs, all from artificial thin stacks). No fixed rule gets above 56 % there (`h24_handoff_test.log:64`).

Scripts `h24_handoff.py`, `h24_handoff_worlds.py`, `h24_handoff_misses.py`. Predictions `h24_handoff_predictions.txt` (sha256 `090efdb7b1472b28fc1dcdda7bc21aef55e8249452600ca60a849348ceb85afb`, written Wed Oct 21 00:39:24 PDT 2026, never edited). Logs `h24_handoff_classify.log`, `h24_handoff_runs_0..3.log`, `h24_handoff_test.log`, `h24_handoff_worlds.log`, `h24_handoff_misses.log`.

## What ran

| step | command | what it produced |
|---|---|---|
| fit / classify | `h24_handoff.py classify` | `h24_handoff_classify.log` (reads `deephand_scale.log`, no simulation) |
| fresh simulations | `h24_handoff.py run p 4` (4 shards) | 1044 ROW lines: 570 fresh (F), 120 thin stacks (T), 354 noise reruns (N); summed simulation time 2544 s (F + N rows; T rows time at 0.0 s each) |
| test | `h24_handoff.py test` | `h24_handoff_test.log`, scorecard P1-P11 |
| worlds | `h24_handoff_worlds.py` | `h24_handoff_worlds.log`, N = 200000, seeds 0-3, 832 s (`h24_handoff_worlds.log:161`) |
| miss tables | `h24_handoff_misses.py` | `h24_handoff_misses.log` (descriptive, written after results) |

Fresh runs use seed 20261021 and `deephand_scale.one` with configurations drawn as in `deephand_scale.config`: w = 0.035, 0.02, 0.01, 0.005 with 150, 150, 150, 120 runs. w = 0.0025 was skipped fresh (about 10 min per run).

## The rules

Domains are maximal runs of equal hand. The rule domain contains r*. If it does not hold, the role passes to the domain above, and so on. Quantities: d/w (thickness in local wall widths at the domain middle), eps (dimensionless drift of its rising top wall), dab (d/w of the domain above).

| rule | domain holds iff |
|---|---|
| THR | d/w >= THR (handoff.py's rule; THR = 1 and 3 published) |
| V1 | d/w >= 1.25 |
| V2 | d/w >= 3, or -2.67 + 5.65 d/w + 1.50 ln eps > 0 (top domain: d/w >= 1.25) |
| V3 | d/w >= 3, or -1.76 + 7.33 d/w + 1.68 ln eps - 1.98 ln dab > 0 (top domain: d/w >= 1.25) |

ln eps and ln dab are clipped to the training range (`h24_handoff_test.log:3-4`). The decision rule (highest old-run CV agreement, ties within 0.2 points go to fewer parameters) picked V3 (`h24_handoff_test.log:12, 31`). If the chain is exhausted the fallback is the rule domain's hand ("deep"), fixed before testing.

## What the misses look like (old runs, descriptive; `h24_handoff_classify.log`)

- The plain deepest-hand rule misses 114 of 600; in all 114 the observed hand is the hand of the domain above the rule domain, and unfinished runs number 0 (`classify.log:2`). The handoff idea is right in structure.
- The 40 misses of THR = 1 (the 93.3 % ceiling) are 33 "kept a domain that lost" (type S) and 7 "handed up a domain that held" (type U). 37 of the 40 have the rule domain at the bottom, 31 have rule d/w in [1, 2), and by n they are 5 / 19 / 16 for n = 3 / 6 / 10 (`classify.log:47`).
- **Thickness.** P(rule domain holds) by its d/w: 1.8 % below 0.5, 19.4 % in [0.5, 1), 45.5 % in [1, 1.5), 80.0 % in [1.5, 2), 97.1 % in [2, 3), 100 % at 3 or more (`classify.log:50-55`). No single threshold can separate a band that holds 45-80 % of the time.
- **Sweep timing.** Domains that held have median Te/Tc 1.35, those that lost 0.83 (`classify.log:59-60`). The largest Te/tmax is 0.347, so no run was cut off (`classify.log:61`). Timing shows a tendency but does not separate the classes.
- **Layer count.** Hold probability for 1 <= d/w < 3 is 81.5 / 68.4 / 84.1 % at n = 3 / 6 / 10 (`classify.log:64-66`); by number of domains there is no monotone trend (`classify.log:67-71`).
- **Noise / seed.** None. All 177 reran configurations reproduce the original hand in every rerun (80 / 80 and 274 / 274; `test.log:18-19`).
- **What predicts holding** (235 training rows, CV accuracy): d/w 84.5 %; plus ln eps (V2) 90.3 %; plus ln dab instead 92.4 %; V3 92.9 % (`classify.log:74-78`).
- **d/w-only ceiling on old runs:** best threshold 1.25 gives 95.2 % in sample, THR = 1 gives 93.3 % (`classify.log:81`). The 95.2 is tuned on the same runs, so it is not a test; the fresh value for the same threshold is 93.9 % (`test.log:28`).

## Fresh-run results (out of sample; `h24_handoff_test.log`)

| rule | overall (570) | fired | exhausted chains |
|---|---|---|---|
| plain | 83.2 ± 1.6 | 0 | 0 |
| THR = 1 | 93.5 ± 1.0 | 86.4 ± 3.8 (n = 81) | 0 |
| THR = 3 | 74.2 ± 1.8 | 33.8 ± 3.2 (n = 222) | 76 |
| V1 (1.25) | 93.9 ± 1.0 | 81.4 ± 4.0 (n = 97) | 4 |
| V2 | 94.7 ± 0.9 | 82.7 ± 3.7 (n = 104) | 0 |
| **V3** | **96.5 ± 0.8** | **86.5 ± 3.4 (n = 104)** | 0 |

(`test.log:25-30`.) V3 against V1: 23 runs only V3 gets, 8 only V1 gets; against V2: 14 vs 4 (`test.log:36-37`). The fallback choice does not matter on these runs (96.5 for deep, top and none; `test.log:35`).

By w: 92.7 (0.035), 94.7 (0.02), 100.0 (0.01), 99.2 (0.005). By n: 100.0, 98.4, 91.1 for n = 3, 6, 10 (`test.log:33-34`). Of V3's 20 fresh misses, 17 are at n = 10, 3 at n = 6, 0 at n = 3 (`misses.log:58-61`); 14 are over-fires and 6 under-fires (`misses.log:30`); 15 sit at rule d/w between 0.5 and 2 (`misses.log:33-35`).

**What this does and does not show.**
- The rule works for w/H between 1/200 and 1/29 at the ice, up to 10 layers, in the 1-D radial model. Worst strata are w = 0.035 and n = 10, at 92.7 % and 91.1 %; the overall 96.5 does not hold in every stratum.
- The V3 gain over the best d/w-only rule is +1.8 points on the fresh runs (96.5 vs 94.7, `test.log:61`); no paired test against that oracle was run (V3 against V1: 23 vs 8 discordant runs, `test.log:36`).
- V3 has four fitted numbers, two of them on features (ln eps, ln dab) clipped to the training range. It is a fitted predictor, not a derived law; the mechanism behind the ln dab term (a thicker neighbour above makes a thin domain more likely to lose) is a reading of the coefficient sign, not a tested cause.
- The old-run V3 numbers (98.2 in sample, 97.7 CV; `test.log:11-12`) are development numbers and are not quoted as validation.

## Thin regime (all domains under one wall width)

Thin stacks (n = 10, w = 0.10 and 0.07, 60 runs each) were the only source of runs where every domain from the top to the rule domain has d/w < 1: n = 68, none from old or fresh runs (`test.log:64`). In that group the final hand equals the rule domain's hand in 55.9 ± 6.1 % and the top domain's hand in 48.5 ± 6.1 %; neither fixed fallback reaches 85 %. V3 gets 65.0 ± 4.4 % on all 120 thin stacks (`test.log:68`), 55.0 % at w = 0.10 and 75.0 % at w = 0.07 (`misses.log:63`). This regime is where the radial rule is not closed. The walls here are 7-10 % of the depth, far from the ocean-like range, so what this says about Earth is limited to: in the wrapper, draws where no layer holds are exactly this kind of draw.

## Worlds: how much the "no layer thick enough" draws move the shares

N = 200000, seed 0, no lid, fitted h (`h24_handoff_worlds.log`). B = V3.

| world | no handoff | THR = 3 | THR = 1 | B none | **B deep** | B top |
|---|---|---|---|---|---|---|
| Early Earth ocean | 99.9 | 94.3 | 97.6 | 98.2 | **99.9** | 100.0 |
| Europa | 38.4 | 32.0 | 34.9 | 35.7 | **39.3** | 41.9 |
| Enceladus | 26.0 | 21.3 | 23.7 | 24.2 | **27.1** | 29.1 |
| Mimas | 16.9 | 14.1 | 15.6 | 16.2 | **17.8** | 19.7 |
| Umbriel | 52.9 | 44.3 | 48.4 | 48.8 | **53.8** | 55.4 |

(P_pved %, `worlds.log:12, 13, 21, 25, 26`; all 15 worlds are in that table.)

- **Movement of the no-layer draws** (B deep minus B none, P_pved points): Early Earth +1.8; moons +1.6 (Mimas) to +5.0 (Umbriel, Ariel); across the 15 worlds +1.6 to +5.0 (`worlds.log:102-116`).
- **How many draws are involved** (B none, chain exhausted, % of recursion draws): Early Earth 5.4; moons 12.9 (Titan) to 30.6 (Miranda) (`worlds.log:84-98`). Under THR = 3 it was 17.1 to 41.1 (same lines).
- **Early Earth bracket:** 98.2 (none) to 99.9 (deep), 100.0 (top). All three are above the 97.6 of the one-width cut and the 94.3 of the three-width cut. The bracket is the honest statement; the thin-regime result says the "deep" assignment is right only about 56 % of the time on those draws, so the middle of the bracket is not more likely than the ends.
- **Moons** under B deep sit 0.0 to 1.1 points above the no-handoff value, not below it (`worlds.log:121-134`, "B deep-off"). That differs from the THR = 3 result of `new_analyses_handoff.md`, where every moon fell 3-9 points below the main value.
- **Rank stability:** Spearman 0.996 between B deep and THR = 1 across the 15 worlds (`worlds.log:158`). Seed spans are at most 0.36 points (`worlds.log:159`).
- **Not pre-registered, exploratory:** with unclipped V3 features Early Earth is 99.5 with 1.4 % of recursion draws exhausted, moons 7.5 to 22.0 % exhausted (`worlds.log:48, 84-98`). Part of the exhausted share under the clipped rule is therefore a clipping effect of wrapper features that lie outside the radial training range. This has not been tested against simulations.

## Predictions

| # | prediction | result |
|---|---|---|
| P1 | V3 on fresh: overall >= 95 and fired >= 85 | **held** (96.5; 86.5; `test.log:71`) |
| P2 | decision rule picks V3 | **held** |
| P3 | V3 >= V1 + 1.0 | **held** (+2.6) |
| P4 | V1 >= THR = 1 + 0.5 | **failed** (+0.4; 93.9 vs 93.5) |
| P5 | THR = 3 < 85 overall | **held** (74.2) |
| P6 | d/w-only oracle <= 96.0 and V3 beats it by >= 0.5 | **held** (94.7; +1.8) |
| P7 | V3 >= 90 in every w and n stratum | **held** (worst 91.1) |
| P8 | reruns reproduce hand: THR = 1 misses >= 90 %, borderline >= 85 % | **held** (100.0; 100.0) |
| P9 | thin group G >= 20 runs | **held** (68) |
| P10 | in G, top-hand share exceeds rule-hand share by >= 10 points | **failed** (48.5 vs 55.9, the other way) |
| P11 | in G, neither fixed fallback reaches 85 % | **held** (best 55.9) |
| W1 | hand2 reproduces handoff.hand exactly | **held** |
| W2 | bisection shortcut mismatches <= 0.5 % | **held** (0 / 22500 each) |
| W3 | Early Earth final P_pved under B/deep in [97.0, 99.9] | **failed** (printed 99.9; the script compares the unrounded value, which is just above 99.9; `h24_handoff_worlds.py:284`, `worlds.log:156`) |
| W4 | M in [0, 10] everywhere, Early Earth <= 2.5 | **held** (+0.0 to +5.0 over no-lid and final; Earth +1.8) |
| W5 | Spearman >= 0.98 | **held** (0.996) |
| W6 | seed span <= 0.5 | **held** (0.36) |

14 held, 3 failed. P4 and W3 are small misses. P10 failed in direction as well: my prior that the chain keeps passing the role upward does not describe the thin regime; the deepest domain does slightly better than the top one, and both are near a coin flip.

## Proposed wording (not applied)

- **H24 (replacing the "Still open" part of item 24):**
  > "Out-of-sample test (`h24_handoff.py`): a handoff rule built from the deciding domain's thickness in local wall widths, the drift of its rising wall and the thickness of the domain above (four coefficients fitted on the 600 earlier runs, frozen before testing) predicted the filling hand in 96.5 ± 0.8% of 570 fresh runs (86.5 ± 3.4% of the 104 where it handed up; `h24_handoff_test.log:32`). Rules using thickness alone reach at most 94.7% even when tuned on the fresh runs (`h24_handoff_test.log:61`). Range tested: walls 1/200 to 1/29 of the depth, up to 10 layers, one-dimensional radial model. Open: stacks where every domain is thinner than one wall width (no fixed rule above 56% in 68 test runs, `h24_handoff_test.log:64`), and walls thinner than 1/200 of the depth were not re-tested."
- **E58 follow-up / Section 2.4:**
  > "A one-wall-width handoff predicts 93.3 ± 1.0% of the 600 earlier runs; adding the drift of the rising wall and the thickness of the layer above raises this to 96.5 ± 0.8% on 570 new runs (`h24_handoff_test.log:32`), while the best thickness-only cut gets 94.7% (`h24_handoff_test.log:61`)."
- **Section 3.3, Early Earth:**
  > "With the validated handoff rule, Early Earth is 98.2% if draws in which no layer is thick enough count as no hand and 99.9% if the deepest layer decides them (`h24_handoff_worlds.log:12`), against 99.9% with no handoff."
- **Moons caption:**
  > "Under the validated handoff rule the moons' P_pved is 0.0-1.1 points above the no-handoff value when draws with no thick layer are decided by the deepest layer, and 1.6-5.0 points lower if they count as no hand (`h24_handoff_worlds.log:102-134`); 13-31% of recursion draws are of this kind (`h24_handoff_worlds.log:84-98`)."

## Disclosures

- **Development before the hash.** The 600 old runs were read in full and fitted in scratchpad exploratory fits and one `h24_handoff.py fit` run before `h24_handoff_predictions.txt` was hashed. The variants were chosen and fitted on those runs only. Old-run agreement numbers (including 98.2 for V3) are development numbers.
- **Smoke runs.** `h24_handoff_worlds.py` was run once with `SMOKE=1` (N = 2000, Europa) before the hash: hand2 matched handoff.hand and the bisection check showed 0 mismatches in 1438 draws, so W2 was not blind. Both docstrings say so.
- **Post-hash debug.** One noise-rerun task was run once in a scratch script to check the numpy-2 patch below (disclosed in the `h24_handoff.py` docstring).
- **numpy 2.** `deephand_scale.one` fails under numpy 2 (`deephand_scale.py:165`, boolean numpy scalar indexing a str). `h24_handoff.py` wraps `DH.config` to return Python-int hands (types only, no numeric change). `deephand_scale.py` itself is untouched.
- **`h24_handoff_misses.py`** was written after the test results were seen; it is exploratory and nothing above is a test result from it (it re-tabulates the same runs).
- **Not tested:** w = 0.0025 fresh; more than 10 layers; non-radial geometry.
- **No predictions file was edited** after the results run; no v2 file was needed. Existing files were not edited, nothing was committed.

## Found outside the unit

- `deephand_scale.py:165`: `"+-"[h < 0]` raises a TypeError under numpy 2 when `h` comes from `rng.choice` (numpy integer). The earlier logs were produced under an older numpy, so `deephand_scale.py` cannot be re-run as is.
