# Enceladus seed collision: fix and regenerated numbers (seedfix.py, 2026-10-12)

Script: `seedfix.py` (predictions committed before the run, f8ef1c6). Log: `seedfix.log` (N = 200000, 137 s).
No existing script, log or md was edited. `journals/ija_draft.md` was only read.

## The bug in plain words

Each world gets three extra random inputs: the rate-bias factor eta, the feedstock share r and the D_z ratio rho.
They are drawn from `np.random.default_rng([seed, i])`, where i is the world's index. numpy drops trailing zeros
from a seed, so for world 0 (Enceladus) `[seed, 0]` is the same seed as `seed`. That is the generator that already
drew Enceladus's own bias g, rate k2 and concentration. So Enceladus's eta followed g exactly, r followed k2, and
rho followed concentration (correlation 1.000000 at every seed; seedfix.log:5, 17, 29, 41). The other 14 worlds
use `[seed, i]` with i ≥ 1, which are independent streams, so only Enceladus was affected.

The bug was in five places: sweep2.py:165, sweep3.py:156, multiwall.py:435, rho2side.py:96 and rerun_base.py:170.

## The fix

`seedfix.py` wraps `np.random.default_rng` before importing the project modules. Only the exact seed `[s, 0]` is
changed, to `[s, 0, 1]`. Every other call is passed through untouched, so worlds 1–14 keep their streams.

I chose this so the regression would be exact. Every non-Enceladus number has to come out the same as in the logs,
to the printed digit, which shows the wrapper changed nothing else. `[s, 0, 1]` is also the stream
`rho2side.posthoc` used, so its post-hoc Enceladus values served as predictions.

After the fix, the Enceladus extras are uncorrelated with its base draws: |corr| ≤ 0.0038 at seeds 0–3
(seedfix.log:5, 17, 29, 41). The remapped call sites are listed at seedfix.log:475.

## Old → new, Enceladus, seed 0 (P_pved / P_chance / P_none %, no lid; with a lid it is 0.0 everywhere)

| quantity | old | new | seeds 0–3, new | seedfix.log |
|---|---|---|---|---|
| sweep2 step b (with η) | 23.3 | 23.8 | | 72 |
| sweep2 step c (D_z floor) | 27.2 | 27.8 | | 72 |
| sweep2 step d (racemization) | 26.4 | 27.1 | | 72 |
| sweep2 step e (feed/product) | 24.5 / 6.8 / 68.6 | 24.7 / 6.6 / 68.6 | 24.5–24.7 | 6, 72 |
| sweep2 step f (D_z rising) | 3.9 / 2.0 / 94.1 | 6.7 / 1.5 / 91.8 | 6.6–6.7 | 7, 75 |
| step f, Eq. 20 ×5 / Asp | 16.9 / 3.8 | 19.9 / 6.5 | | 78 |
| sweep3 main | 6.6 / 2.6 / 90.8 | 9.0 / 2.1 / 88.8 | 9.0–9.0 | 8, 99 |
| sweep3 bound (iii) | 23.9 / 6.9 / 69.2 | 23.6 / 7.3 / 69.1 | 23.4–23.6 | 9, 105 |
| sweep3 reversed table | 72.1 / 68.5 / 95.1 / 0.2 / 0.5 | 66.6 / 63.3 / 95.1 / 0.3 / 0.7 | | 102 |
| multiwall main | 28.3 / 8.8 / 63.0 | 26.4 / 9.8 / 63.9 | 26.2–26.4 | 10, 122 |
| multiwall cascade ×10 | 19.9 / 5.7 / 74.3 | 19.4 / 6.4 / 74.1 | 19.3–19.4 | 11, 125 |
| multiwall cascade ×5 | 22.4 / 6.6 / 71.0 | 21.5 / 7.3 / 71.2 | 21.3–21.5 | 12 |
| multiwall sum rule | 20.1 / 5.8 / 74.0 | 19.9 / 6.3 / 73.8 | 19.7–19.9 | 13, 125 |
| multiwall recursion | 68.5 / 40.6 / 97.9 / 74 | 63.3 / 39.5 / 96.6 / 73 | | 128 |
| rho2side two-sided main | 26.5 / 7.1 / 66.4 | 25.6 / 8.2 / 66.2 | 25.3–25.6 | 14, 148 |
| rho2side two-sided ×10 | 20.8 / 6.6 / 72.6 | 22.1 / 6.5 / 71.4 | 21.8–22.1 | 15, 151 |
| ×10 gap, two-sided vs mixture | 1.45 | 0.00 (all-world max 1.45 → 0.07) | 0.00–0.04 | 16, 160 |
| rho2side real-profile sink | 26.5 → 27.2 (+0.75) | 25.6 → 27.3 (+1.71) | | 157 |
| PVED share, P_pved/(P_pved+P_chance) | 76% | 73% | | 122 |

Seeds 1–3 move every Enceladus no-lid column by at most 0.28 points from seed 0 (seedfix.log:55–64).

### rerun_base (H19/H20/H27); only Enceladus rows changed (seedfix.log:169–453)

- **Section 0** (line 4): main 26.4 / 0.0, ×10 19.4 / 0.0, three-way split 26.4 / 9.8 / 63.9. Section 0 still
  compares against the real `multiwall.log`, so its P1 now prints FAILED for Enceladus (line 174). That is expected.
  The seed-check line now reads `False` with corr ≈ 0 (line 177), because the wrapper is active.
- **Section 1 sensitivities, main**: D_z +44.9 → +44.0 (still the top effect); concentration +32.0 → +29.7;
  time +16.6 → +16.6; k2 +9.5 → +6.4; eta +12.1 → +6.3; rho +32.0 → +6.0; g +12.1 → +5.9; D_h −3.4 → −3.6;
  r +9.5 → +3.0; volume −3.0 → −2.7; F +0.8 → +0.8. The three duplicated pairs no longer give identical effects.
- **Section 2 breakeven**: Enceladus D_z crossing 1.6e-5 → 2.6e-5 (×10: 7.8e-5 → 8.9e-5). Concentration and
  rho no longer reach 50% (they were artifacts of the collision), so rerun_base P5 flips FAILED → HELD (line 339).
- **Section 3 assumed rules** (lines 342–369): Enceladus drops 0.9–1.9 points per row (e.g. whole ocean
  31.8 → 30.4; bias unlimited 37.0 → 36.1). rerun_base P7 flips HELD → FAILED, because its reference 37.1 was a
  collided value.
- **Section 4 wobble** (375–384): 28.3 / 27.6 / 25.6 / 21.7 → 26.4 / 25.7 / 23.6 / 19.9.
- **Section 5 H20** (390): 28.3 / 27.9 / 28.2 || 19.9 / 19.7 / 19.9 → 26.4 / 26.0 / 26.3 || 19.4 / 19.2 / 19.4.
- **Section 6 Kang / strat-alt finals** (393–417): 0, 0, 7.9, 5.5, 22.9, 19.0, 2.2, 9.5, 9.9 →
  0, 0, 7.4, 5.2, 21.3, 17.7, 2.0, 8.8, 9.3. The pinned shares are unchanged. rerun_base P10 flips HELD → FAILED,
  because 22.9 and 19.0 are hard-coded as expected values.
- **Section 9 k2 / Bennu** (420–432): main k2 1e-12..1 wins 18.6 → 16.7; R2 P(win|rac) 1.9 → 1.7;
  R6 step 20.8 / 23.3 → 18.5 / 20.5.
- **Section 10 magnetite** (435, 438): main one src 2.7e-16 → 1.6e-16, f*_G 4.8e-16 → 2.9e-16,
  M=1e3 1.5e-14 → 9.0e-15, field 4.0e-15 → 2.6e-15. The ×10 row moves similarly. rerun_base P14 flips
  FAILED → HELD.
- **Prediction line** (450): P1 held → FAILED (section-0 reference), P5 FAILED → held, P7 held → FAILED,
  P10 held → FAILED, P14 FAILED → held.

Every proposal row in `journals/new_analyses_rerun.md` that quotes these Enceladus values is stale by the same
amounts.

## Draft lines that change (journals/ija_draft.md, read only)

| line | old | new | seedfix.log |
|---|---|---|---|
| 115 (E29) | "0% with a lid; 20–28% without"; "no lid 19.9 cascade ×10 to 28.3 main" | "19–26% without"; "19.4 … 26.4" | 125 |
| 126 (E40) | Enceladus "20–28% to 0%" | "19–26% to 0%" | 125 |
| 127 (E41, Table 3) | "Enceladus 0 → 0 (20–28) (:556)" | "0 → 0 (19–26)"; the order Ganymede > Europa > Enceladus holds | 125 |
| 132 (E46) | Enceladus chance 8.8; PVED share "76% Enceladus and 74% at the minimum (Mimas)" | 9.8; "73% Enceladus, now the minimum (Mimas 74%)" | 122 |
| 133 (E47 ladder) | 21.8 → 24.1 → 23.3 → 27.2 → 26.4 → 24.5 → 3.9 → 6.6 → 28.3 (×10 19.9) | 21.8 → 24.1 → 23.8 → 27.8 → 27.1 → 24.7 → 6.7 → 9.0 → 26.4 (×10 19.4) | 72, 96, 125 |
| 236 (abstract) | "0% for Enceladus … (20–28% without it)" | "(19–26% without it)" | 125 |
| 580 | "20–28% for Enceladus" | "19–26%" | 125 |
| 585 | "with η: 23.3" | 23.8 (published PVED range 24.1 unchanged) | 72 |
| 586 | D_z floor 27.2 | 27.8 | 72 |
| 587 | racemization 26.4 | 27.1 | 72 |
| 588 | separate feedstock and product 24.5 | 24.7 | 72 |
| 589 | D_z rising with depth 3.9 | 6.7 | 72 |
| 590 | reversed walls counted 6.6 | 9.0 | 96 |
| 591 | layer fills the ocean 28.3 (×10 19.9) | 26.4 (×10 19.4) | 125 |
| 595 | "Enceladus 8.8%"; "74% (Mimas) to 100% … 76% for Enceladus" | 9.8%; "73% (Enceladus) to 100% …"; Mimas 74% | 122 |
| 617 | "Enceladus drops from 20–28% to 0%" | "19–26%" | 125 |
| 632 (Table 3) | "0% (20–28%)" | "0% (19–26%)" | 125 |
| 740 | "20–28% for Enceladus" | "19–26%" | 125 |

Lines that do not change:

- 21.8 (step 0) and 24.1 (step a): these steps do not use the extras.
- "chance at most 14.9% (Umbriel)".
- Every Europa and early Earth number.
- The "Asp at most 0.7 points" claim: for Enceladus the Asp move is now 0.2 points.
- The uniform-D_z base-model values: lines 196, 597–607, 644, 662 and 741, including the Kang 4–6%, 15–18% and 0–18%.

There is one new qualitative change. Enceladus replaces Mimas as the world with the lowest PVED share: 72.9% against
Mimas's 73.9% (`multiwall.log:538`).

## Predictions (pre-registered in seedfix.py; 14 held, 2 failed; seedfix.log:457–474)

- **Held:**
  - P1: extras uncorrelated, |corr| ≤ 0.0038.
  - P3: with a lid, 0.0 everywhere.
  - P4: main 26.4.
  - P5: step f 6.7.
  - P6: ×10 gap 0.00; all-world maximum 0.07.
  - P7: post-hoc values reproduced.
  - P8: sweep3 9.0.
  - P9: bound 23.6.
  - P10: sum rule 19.9.
  - P11: chance 9.8 and 1.5.
  - P12: ladder moves ≤ 0.7.
  - P14: Kang rows scale with main, maximum deviation 0.07.
  - P15: f* within a factor 2 (0.6×).
  - P16: seed spread ≤ 0.28.
- **P2 failed as registered.** The classifier flagged one line, rerun_base.log:305: "P14 FAILED: Earth L 46.6 …".
  Its Earth-named first segment changed only in the verdict word, FAILED → HELD, which flipped because the
  Enceladus f* values later on the line moved. Every number naming another world on that line is identical, and so
  is every non-Enceladus field in all five logs.
- **P13 failed.** g and eta now have nearly equal effects (+5.9 vs +6.3, 0.4 apart, short of the > 1 required).
  k2 vs r (+6.4 vs +3.0) and concentration vs rho (+29.7 vs +6.0) do separate. D_z stays on top.

## Risks

- `front2d.py:201`, `deephand_scale.py:299` and `pin_onset.py:254` have the same `[seed, 0]` pattern. Their
  Enceladus numbers are still collided; seedfix does not cover them.
- seedfix is a wrapper: running the project scripts directly still gives the collided numbers. The permanent fix that
  keeps worlds 1–14 exact is `default_rng([seed, i] if i else [seed, 0, 1])` at each of the 8 sites.
