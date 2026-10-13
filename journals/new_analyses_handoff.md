# handoff: replacing the 89 % ocean-scale floor with an estimate (2026-10-12)

**Short answer.** If a deciding layer thinner than three wall widths hands the decision to the layer above it, Early Earth comes out at **94.3 %** (the floor was 89.0, multiwall main 99.9), and every icy moon lands 2.5–5.9 points above its floor and 3.2–9.4 points below multiwall main. The simulations do not support the three-wall-width cut, though. It gets 76.8 % of the 600 simulated stacks right, which is worse than doing no handoff at all (81.0 %). A one-wall-width cut gets 93.3 % right and gives Earth **97.6 %**. Neither cut passes the pre-registered bar (≥ 95 % overall and ≥ 85 % on the runs where the handoff fires).

Script `handoff.py`, log `handoff.log`. Predictions were committed before the run in eea83bd. There was one smoke run beforehand, which only ran self-checks and printed "smoke ok". The result run took 18 s.

## Method in plain words

- **Seed fix.** All world numbers use seedfix's correction: the seed `[s, 0]` becomes `[s, 0, 1]`, so Enceladus gets its own extra draws. Worlds 1–14 are unchanged. Nothing in the project was edited. The wrapper installs the fix and then imports multiwall, sweep2, sweep3 and deephand_scale.
- **The handoff rule (in multiwall's wrapper).**
  - Layers are d_w thick, counted from the ice down, with a bottom remnant. The layer just below the last wall past onset decides.
  - Its thickness in local wall widths (d/w) is measured the way deephand_scale did it.
  - If that layer is under THR = 3 wall widths, the decision passes to the full layer above it. If that one is also too thin, it passes up again, and so on.
  - If no layer up to the ice is thick enough, the draw counts as no hand.
  - Whichever layer ends up deciding gets its own PVED test, F·Δ ≥ 2, with Δ computed from that layer's own molecule count. It also gets its own rise time: the deepest wall above it has to reach the ice.
  - The sinking term is kept as it was. This is an assumption: the thin layers below are taken to die before any wall moves, which is what happened in the simulations.
  - Timing is run two ways: with multiwall's closed form (h = 1), and with deephand_scale's fitted correction h, evaluated for the new deciding layer.
- **Validation against deephand_scale's 600 simulated stacks.**
  - Each stack's layers are rebuilt from its logged R, n, ln ρ and u. Neighbouring layers of the same hand are merged into one domain.
  - Starting from the domain the old rule picks, the decision is passed up the same way.
  - The prediction is the hand of the first domain thick enough, and it is compared with the hand that actually filled.
  - Checks: THR = 0 (no handoff) reproduces deephand's own "rule held" flag in all 600 runs, and the rebuilt d/w matches the logged d/w in all 600.
- **Regressions, all EXACT.**
  - (A) With the handoff off and h = 1, the result equals multiwall.multi array for array, and matches seedfix's multiwall main to the printed 0.1 for all 15 worlds.
  - (B) With the handoff off and fitted h, it matches deephand_scale's fitted-h column for worlds 1–14.
  - (C) The recomputed ocean-scale floor matches deephand_scale's column for worlds 1–14.

## Validation (± one standard error)

| handoff cut | overall agreement | runs where it fires | of those, agree | predicted "no hand" |
|---|---|---|---|---|
| none (plain rule) | 81.0 ± 1.6 % | — | — | 0 |
| THR = 3 (as specified) | **76.8 ± 1.7 %** | 235 | **40.9 ± 3.2 %** | 63 (all wrong) |
| THR = 1 (secondary, pre-registered) | **93.3 ± 1.0 %** | 88 | **92.0 ± 2.9 %** | 0 |

**Why the three-width cut fails.** It hands off layers that went on to win. Here is how often the cut-3 prediction was right among the runs where it fired, split by the thickness of the original deciding domain:

| original domain thickness | runs | cut-3 prediction right |
|---|---|---|
| 0–1 wall widths | 88 | 68.2 % |
| 1–2 wall widths | 79 | 44.3 % |
| 2–3 wall widths | 68 | 1.5 % |

The thin domain itself filled 97 % of the time at 2–3 widths and 61 % at 1–2. So "≥ 3 widths always fills" (365/365) is a safe sufficient condition, but it is the wrong boundary for deciding when to hand off. The real boundary sits near one wall width.

**By width w.** With THR = 3, agreement is 54.0 % at w = 0.035 and 96.7 % at w = 0.0025. With THR = 1, it is 86.7 % and 100.0 %.

**Exploratory scan (post hoc, not pre-registered).** Overall agreement peaks at 94.8 % at THR = 1.25, where 88.3 % of fired runs are right. No cut reaches 95 %.

## World table: P_pved / P_chance / P_none %, seed 0, N = 200 000, seed fix applied

No lid:

| world | multiwall main | ×10 | ocean-scale floor | handoff, h = 1 | **handoff, fitted h** | THR = 1, fitted h |
|---|---|---|---|---|---|---|
| Early Earth ocean | 99.9 / 0.0 / 0.1 | 99.7 / 0.0 / 0.3 | 89.0 / 0.0 / 11.0 | 94.3 / 0.0 / 5.7 | **94.3 / 0.0 / 5.7** | 97.6 / 0.0 / 2.4 |
| Umbriel | 53.7 / 14.9 / 31.4 | 43.6 / 10.5 / 45.8 | 38.4 / 3.6 / 58.0 | 46.6 / 4.8 / 48.6 | **44.3 / 4.7 / 50.9** | 48.4 / 6.8 / 44.9 |
| Oberon | 51.4 / 11.7 / 36.9 | 40.8 / 7.9 / 51.3 | 36.7 / 2.8 / 60.5 | 44.8 / 3.8 / 51.4 | **42.4 / 3.8 / 53.8** | 46.2 / 5.5 / 48.3 |
| Titania | 51.3 / 11.6 / 37.1 | 40.9 / 8.0 / 51.1 | 36.8 / 2.8 / 60.4 | 44.7 / 3.9 / 51.4 | **42.3 / 3.8 / 53.9** | 46.1 / 5.5 / 48.4 |
| Ariel | 49.6 / 13.5 / 36.9 | 39.7 / 9.4 / 50.9 | 34.5 / 2.9 / 62.5 | 42.6 / 4.1 / 53.2 | **40.4 / 4.0 / 55.6** | 44.3 / 5.9 / 49.8 |
| Ganymede | 44.1 / 3.4 / 52.5 | 33.5 / 2.0 / 64.5 | 33.2 / 0.9 / 66.0 | 40.2 / 1.3 / 58.6 | **37.2 / 1.2 / 61.6** | 39.8 / 1.8 / 58.4 |
| Callisto | 43.5 / 3.8 / 52.7 | 33.0 / 2.3 / 64.7 | 32.4 / 1.0 / 66.6 | 39.4 / 1.4 / 59.2 | **36.5 / 1.3 / 62.1** | 39.1 / 2.0 / 58.9 |
| Dione | 42.0 / 10.5 / 47.5 | 32.2 / 6.9 / 60.9 | 27.9 / 2.0 / 70.0 | 35.8 / 3.0 / 61.2 | **33.5 / 2.9 / 63.6** | 37.2 / 4.4 / 58.5 |
| Pluto | 41.1 / 6.1 / 52.8 | 30.6 / 3.7 / 65.7 | 28.5 / 1.3 / 70.2 | 36.1 / 1.9 / 62.0 | **33.3 / 1.9 / 64.8** | 36.5 / 2.8 / 60.7 |
| Europa | 39.6 / 4.9 / 55.5 | 29.0 / 2.9 / 68.1 | 27.4 / 1.0 / 71.5 | 34.9 / 1.6 / 63.5 | **32.0 / 1.5 / 66.4** | 34.9 / 2.4 / 62.7 |
| Triton | 37.0 / 4.7 / 58.3 | 26.5 / 2.7 / 70.8 | 25.2 / 0.9 / 73.8 | 32.6 / 1.5 / 66.0 | **29.7 / 1.4 / 68.9** | 32.5 / 2.1 / 65.3 |
| Titan | 35.8 / 2.5 / 61.7 | 25.5 / 1.4 / 73.1 | 25.6 / 0.6 / 73.8 | 32.3 / 0.9 / 66.7 | **29.3 / 0.9 / 69.8** | 31.7 / 1.3 / 67.0 |
| Miranda | 27.1 / 8.4 / 64.5 | 21.3 / 5.9 / 72.8 | 20.1 / 3.5 / 76.5 | 24.3 / 4.1 / 71.6 | **22.6 / 4.0 / 73.4** | 24.2 / 4.9 / 70.9 |
| Enceladus | 26.4 / 9.8 / 63.9 | 19.4 / 6.4 / 74.1 | 17.4 / 2.6 / 80.0 | 22.9 / 3.6 / 73.5 | **21.3 / 3.5 / 75.2** | 23.7 / 4.9 / 71.4 |
| Mimas | 17.3 / 6.1 / 76.6 | 11.6 / 3.6 / 84.8 | 11.5 / 2.1 / 86.3 | 15.6 / 2.9 / 81.4 | **14.1 / 2.8 / 83.0** | 15.6 / 3.8 / 80.6 |

Final (with lid), P_pved only: multiwall main | floor | **handoff fitted h** | THR = 1

| world | main | floor | handoff | THR = 1 |
|---|---|---|---|---|
| Europa | 30.1 | 20.8 | **24.4** | 26.5 |
| Dione | 26.2 | 17.5 | **21.0** | 23.2 |
| Miranda | 9.0 | 6.7 | **7.5** | 8.0 |
| Mimas | 5.6 | 3.7 | **4.5** | 5.0 |
| Enceladus | 0.0 | 0.0 | **0.0** | 0.0 |

The other worlds have no lid, so their final values equal the no-lid values. Notes on the table:

- multiwall main with the seed fix is exactly seedfix's value. Only Enceladus differs from multiwall.log, which had 28.3 / 8.8 / 63.0 because of the seed collision.
- The ocean-scale floor here is recomputed with the seed fix. Enceladus changes from 16.9 to 17.4.
- **Seeds 0–3:** handoff P_pved never spreads by more than 0.37 points on any world.

**Where the new P_none comes from (THR = 3).** On 17 % (Earth) to 41 % (Miranda) of recursion draws, no layer up to the ice is three wall widths thick. The full layers near the floor are thin too, so these draws count as no hand. All 5.7 points of Earth's P_none come from this group.

In the simulations, every stack still ended with one hand. At THR = 3, the 63 "no hand" predictions were all misses. So the no-hand convention is a conservative choice.

The handoff fires on 44–65 % of recursion draws and usually moves up 1–2 layers. Miranda's median move is 17 layers.

## Predictions

| # | prediction | result |
|---|---|---|
| P1 | THR 3: overall 74–84 % and below 81.0; fired 40–58 % | **held** (76.8; 40.9) |
| P2 | THR 3 not validated | **held** |
| P3 | THR 1: overall 89–96, fired ≥ 80, fails ≥ 95 | **held** (93.3; 92.0) |
| P4 | THR 0 = ROW ok flag; rebuilt d/w = ROW d/w | **held** (600/600) |
| P5 | regression A exact, 15 worlds | **held** |
| P6 | regressions B, C exact, worlds 1–14 | **held** |
| P7 | handoff fires on 25–65 % of recursion draws, every world | **failed** (narrowly): Enceladus 65.3 |
| P8 | no-hand ≤ 5 % of recursion draws | **failed**: 17.1–41.1 % |
| P9 | floor ≤ handoff (fitted h) ≤ main + 3, every world | **held** |
| P10 | Early Earth ≥ 97 % | **failed**: 94.3 |
| P11 | fitted h vs h = 1 within 2 points | **failed**: up to 3.0 (Ganymede, Titan) |
| P12 | seeds 0–3 span ≤ 0.5 | **held** (≤ 0.37) |
| P13 | THR 1 within 3 points of THR 3 | **failed**: up to 4.1 (Umbriel); Earth 3.3 |

8 held and 5 failed. The main surprise (P8): I expected the full layers above a thin remnant to be thick. On a large share of draws they are thin all the way up, so the three-width cut leaves those draws with no hand.

## What this changes

- **Earth's estimate moves from "≥ 89 %" to 94–98 %.** The value depends on where the handoff cut sits: 94.3 % at three widths, 97.6 % at one width.
  - The three-width cut is the one specified for this run, but the simulations reject it (76.8 % agreement).
  - The one-width cut agrees 93.3 % of the time, still short of the pre-registered 95 %.
  - The defensible statement is 94–98 %, nearer the top end. Both cuts still count every "no layer thick enough" draw as no hand, which the simulations never showed.
- **Moons.** Each moon now sits between its floor and multiwall main, roughly in the middle. On every moon the ×10 column sits 0.7–3.8 points below this estimate. It is a pessimistic stand-in, not a central value.
- **Timing.** The fitted correction now moves shares by up to 3 points, against 1.4 before. The new deciding layer is thicker, so its fitted h is larger.

## Proposed draft sentences (not applied)

- **Early Earth headline (Sections 3.3 and 6), replacing the floor sentence proposed in new_analyses_deephand.md:**
  > "At ocean scale, a deciding layer thinner than a few wall widths hands the outcome to the layer above it; putting that handoff into the wrapper gives Early Earth 94–98% (94.3% if the cut is three wall widths, 97.6% at one; `handoff.log`), against 99.9% without it."
- **Section 2.4, after the multiwall paragraph:**
  > "In 600 radial simulations, passing the decision to the next layer up when the deciding layer is under one local wall width predicted the final hand in 93.3 ± 1.0% of runs (92.0 ± 2.9% of the runs where it applied); a three-width cut predicted only 76.8 ± 1.7%, because layers of one to three wall widths usually won (`handoff.py`)."
- **World table caption:**
  > "Main values assume the deciding layer always decides; with the handoff rule (three-width cut, fitted cascade time) each moon's P_pved falls by 3–9 points (e.g. Europa 39.6 → 32.0 no lid, 30.1 → 24.4 final; Enceladus 26.4 → 21.3 no lid)."
- **H24, close:**
  > "Handoff tested (`handoff.py`): thin deciding layers pass the outcome up; the right cut is about one wall width (93% agreement), not three (77%). Open: draws where no layer up to the ice is thick enough (17–41% of recursion draws) are counted as no hand."
