# Enceladus seed collision, part 2: front2d, deephand_scale, pin_onset (seedfix2.py, 2026-10-13)

Script: `seedfix2.py`. The predictions were committed before the run in d9cf389. The results and the log were committed in 3df97ea.
Log: `seedfix2.log` (N = 200000, 70 s). No existing script or log was edited. `journals/ija_draft.md` was only read.

## Bottom line

**No number quoted in `journals/ija_draft.md` changes.** Every manuscript aggregate that includes Enceladus reproduces
exactly (seedfix2.log:53, P13). The Enceladus rows of the three scripts change as expected. None of those rows is
quoted in the draft, but one E-note and one H-note describe the old state, so the text needs updating (see "Draft
lines").

## What was done

The wrapper is the same one seedfix.py uses. Only the exact seed `[s, 0]` becomes `[s, 0, 1]`, so worlds 1–14 keep
their streams. That makes the regression exact: every non-Enceladus field must come out to the printed digit.

- **rerun_base.py: skipped.** seedfix.py Section F already reran all of its result sections (0–11) under this
  wrapper. E55 and E56 already quote those values (seedfix.log Section F).
- **front2d.py:** reran Part B (`worlds()`) only. Part A draws no world inputs. `worlds()` asserts against
  multiwall.log, so it ran on a temp copy of multiwall.log. In that copy the six Enceladus lines were replaced by
  seedfix.log's regenerated lines, and each old line was checked against the real log first (seedfix2.log:2).
- **deephand_scale.py:** reran Section 3 only. Sections 1–2 are world-agnostic simulations, and their ROW lines come
  from the real log.
- **pin_onset.py:** reran sections 2–3 as they are, to get the exact κ, SE and band. Section 4 ran from its own
  source with one change: its hard-coded sweep2.log regression target for Enceladus was itself a seed-collided
  number, so (3.9, 16.9) was replaced by seedfix's (6.7, 19.9).

Regression (seedfix2.log:6, 20, 31, 42):

- 0 violations and 0 aggregate lines changed in all three logs (front2d.log 57/57, deephand_scale.log 199/199,
  pin_onset.log 71/71 lines).
- pin_onset sections 2–3 are identical.
- Both front2d regressions print EXACT, deephand_scale prints EXACT, and the patched pin_onset regression passed.
- Remapped call sites are exactly front2d.py:201 ×6, multiwall.py:435 ×2, deephand_scale.py:299 ×1 and
  pin_onset.py:254 ×1 (seedfix2.log:56).

## Old → new, Enceladus, seed 0 (P_pved / P_chance / P_none %, no lid; with a lid it is 0.0 everywhere)

| quantity | old (log line) | new | seedfix2.log |
|---|---|---|---|
| front2d: % u < 3, all draws | 4.2 (front2d.log:125) | 4.2 | 9 |
| front2d: % u < 3, draws reaching one hand | 5.0 (:125) | 4.4 | 9 |
| front2d: u quartiles / % u < 1 / t_mix > τ | 28.1 131 612 / 1.0 / 0.0 | unchanged | 9 |
| front2d: shares, front rule off = on | 28.3 / 8.8 / 63.0 (:143) | 26.4 / 9.8 / 63.9 | 12 |
| front2d: P_pved change, main / u*=1 / v×A | −0.00 / +0.00 / +0.00 (:160) | unchanged | 15 |
| front2d: P_pved change, v/A | −0.04 (:160) | −0.03 | 15 |
| front2d: seed-1 maximum over worlds | 0.36 (:162) | 0.36 | 46 |
| deephand: multiwall main | 28.3 / 8.8 / 63.0 (deephand_scale.log:883) | 26.4 / 9.8 / 63.9 | 23 |
| deephand: multiwall ×10 | 19.9 / 5.7 / 74.3 (:883) | 19.4 / 6.4 / 74.1 | 23 |
| deephand: fitted h | 27.1 / 10.9 / 62.0 (:883), −1.2 vs main | 26.0 / 11.8 / 62.2, −0.4 | 23, 48 |
| deephand: ocean-scale | 16.9 / 3.1 / 80.0 (:883), −11.4 vs main | 17.4 / 2.6 / 80.0, −9.0 | 23, 49 |
| deephand: deciding layer < 5 / < 1 widths | 68.6 / 45.5 (:900) | 69.2 / 46.6 | 26 |
| deephand: median h | 0.79 (:900) | 0.69 | 26 |
| pin_onset: step f at κ, κ±SE, band, ×1 | 3.9 / 0.0 (pin_onset.log:73) | 6.7 / 0.0 | 34 |
| pin_onset: step f ×5 | 16.9 / 0.0 (:73) | 19.9 / 0.0 | 34 |
| pin_onset: layered / pinned at κ / median x | 95.0 / 84.7 / 3.29 (:73) | 95.0 / 84.9 / 3.29 | 34 |
| pin_onset: three-way split at κ | 3.9 / 2.0 / 94.1 (:90) | 6.7 / 1.5 / 91.8 | 37 |

The pin_onset step-f values equal seedfix's regenerated sweep2 step f (seedfix.log:75, 78), as they must.

## Manuscript aggregates that include Enceladus: all unchanged

| draft line | quoted | what sets it | after the fix | seedfix2.log |
|---|---|---|---|---|
| 213 (H22), 556 (assumption 8), 714 (Section 4), 158 (E59) | P_pved change ≤ 0.01 points | several worlds at ±0.01 | Enceladus stays at 0.00 | 15, 46, 53 |
| 714, 158 | ≤ 0.08 with a front 5.5× slower (v/A) | Ganymede −0.08 | Enceladus −0.04 → −0.03 | 15, 53 |
| 556, 714, 158 | fronts in 2–17% of draws (2.0–17.2), early Earth 2.8% | Miranda, Ganymede, Earth | Enceladus 4.2 unchanged | 9, 53 |
| 214 (H23), 464 (Section 2.4), 157 (E58) | fitted time moves no share by > 1.4 points | Ganymede, Titan −1.4 | Enceladus −1.2 → −0.4 | 48, 53 |
| 215 (H24) | 52–74% of decided draws under five wall widths | Earth 51.7, Mimas 73.6 | Enceladus 68.6 → 69.2 | 26, 53 |
| 598, 723 | "about half" of deciding layers thinner | all worlds | unchanged | 26 |

Ocean-scale losses (not quoted in the body; the range is Mimas −5.8 to Umbriel −15.3): Enceladus −11.4 → −9.0 stays
inside it (seedfix2.log:49, 53). pin_onset's quoted numbers (κ, E48) are world-agnostic, and sections 2–3 are
identical.

## Draft lines that change (journals/ija_draft.md, read only; main session integrates)

No result sentence changes. Three notes describe the old state:

1. **H30 (line 221)**, "Seed reuse in older scripts." It says `rerun_base.py`, `front2d.py`, `deephand_scale.py` and
   `pin_onset.py` "still build Enceladus's extra-input stream with `default_rng([s, 0])`… Their Enceladus rows carry
   the collision." Proposed: prefix the note with **[Resolved by E66.]** and replace the second and third sentences
   with:
   > Their Enceladus rows were regenerated on the fixed stream without editing the scripts: rerun_base in
   > `seedfix.log` Section F (E55, E56), the other three in `seedfix2.log` (E66). No quoted number changed. The scripts
   > themselves still contain `default_rng([s, 0])`; making the fix permanent in their source is an author decision.
2. **E59 (line 158)**, last sentence but one: "Enceladus rows in `front2d.py` carry the seed collision (H30); only
   changes, not Enceladus levels, are quoted from it." Proposed:
   > Enceladus rows in `front2d.log` carry the seed collision (H30); `seedfix2.log` regenerates them on the fixed
   > stream (E66) and none of the quoted changes moves.
3. **E58 (line 157)**, header: "…both use the seed fix for Enceladus". This is true for `handoff.py` but not for
   `deephand_scale.py`: deephand_scale.log:883 shows the collided 28.3 / 8.8 / 63.0. Proposed header ending:
   > "…`handoff.py` applies the seed fix for Enceladus; `deephand_scale.log`'s Enceladus rows are regenerated with
   > it in `seedfix2.log` (E66)".

### Proposed E-note (E65 is taken by the H7 full-text checks, so this is 66)

> 66. **Seed fix carried to front2d, deephand_scale and pin_onset (2026-10-13; H30; `seedfix2.py`, `seedfix2.log`,
> `journals/new_analyses_seedfix2.md`; predictions committed first in d9cf389, 12 of 13 held).** The `seedfix.py`
> wrapper (world 0's extra draws from `[s, 0, 1]`) was rerun over `front2d.py` Part B, `deephand_scale.py` Section 3
> and `pin_onset.py` Section 4, with no script edited. `rerun_base.py` was already covered by `seedfix.log` Section F.
> Every other world reproduces to the printed digit (0 violations in 327 compared lines; `seedfix2.log:6, 20, 31`).
> Enceladus no lid: front rule off/on 28.3 → 26.4; fitted h 27.1 → 26.0 (shift −1.2 → −0.4); ocean-scale 16.9 →
> 17.4; thin deciding layers 68.6 → 69.2%; pin_onset step f at κ 3.9 → 6.7, pinned at κ 84.7 → 84.9%
> (`seedfix2.log:12, 23, 26, 34`). With a lid it stays 0.0. No quoted number changes: ≤ 0.01 and ≤ 0.08 points
> (fronts), 2–17%, ≤ 1.4 points (fitted time) and 52–74% all reproduce (`seedfix2.log:53`). The failed prediction
> (P8) is that the fitted-time shift was smaller than registered (−0.4, not −0.5 to −1.4). H30 resolved; E58's and
> E59's Enceladus caveats are superseded.

## Not in the draft, but changed

- `pin_onset.py` RESULTS block (P8 there) says "Enceladus 84.7" pinned. After the fix it is 84.9. Not quoted.
- `front2d.py` and `deephand_scale.py` RESULTS blocks and `journals/new_analyses_front.md` /
  `new_analyses_deephand.md` quote the collided Enceladus rows (28.3, 27.1, 16.9, 68.6). They are historical records
  and were left as they are.

## Predictions (pre-registered in seedfix2.py; 12 held, 1 failed; seedfix2.log:41–55)

P1–P7 and P9–P13 held. P8 failed: the fitted-h shift for Enceladus is −0.4, not within [−1.4, −0.5]. It is smaller
than predicted, so the manuscript's 1.4-point bound is unaffected.

## Risks

- The patched multiwall.log copy depends on seedfix.log's Section D lines. Each old line was checked against the
  real multiwall.log before replacement, so the copy is exact.
- pin_onset's section 4 ran from its source with one string substituted (the Enceladus regression target; asserted
  to occur exactly once). Everything else is the original code, and it was compiled under its own file name and
  line numbers, which is why the remap shows pin_onset.py:254.
- The old logs (front2d.log, deephand_scale.log, pin_onset.log) still show the collided Enceladus rows. A reader who
  reruns those scripts without the wrapper gets the old numbers. Making the fix permanent would mean editing the
  scripts, which is an author decision (H30).
