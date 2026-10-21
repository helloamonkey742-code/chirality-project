# dh_range: how far do the per-world shares move if the icy-moon D_h upper bound goes from 1e2 to 3e2 and 1e3? (2026-10-21)

Script `dh_range.py`, log `dh_range.log`, predictions `dh_range_predictions.txt` (written before any results run, never edited after; sha256 `5b75060fbfaca2f1efb64b7e5be7132050dfe5e538409178a44d97d8cfe591e2`, hashed Wed Oct 21 00:59:38 PDT 2026, asserted in the script). No existing file was edited, renamed or deleted: `uncertainty.py`, `ocean_worlds.py`, `multiwall.py` and `seedfix.py`'s wrapper idea are used by import and monkeypatch in `dh_range.py` only.

Disclosure. (1) An earlier hash of the predictions file, `b203cdbe...` (00:59:04), was replaced by the one above after I corrected wrong log line numbers in it; no run happened in between. (2) Before the full run I made two smoke runs with `DH_SMOKE=1` (N = 2000): the first stopped at `No module named matplotlib` (`uncertainty.py` imports `ocean.py`, which imports matplotlib; the run command needs `--with matplotlib`), the second ran to the end in 4 s. They were for code paths and layout only; no smoke number is used here. (3) After the smoke runs I changed the script in three places before the full run: run4 outputs are compared with exact array equality instead of `repr` (which truncates the `cnt` arrays), the block-order half of P6 was added to the scoring code (it was in the predictions file already), and the docstring was updated. The predictions file was not touched. (4) One cosmetic flaw in the log: line 102 prints `np.float64(0.0)` instead of `0.0`.

## What the source says (I read it myself, 2026-10-21)

`uncertainty.py:38-39` says Ashkenazy and Tziperman 2021 "use 30" as a horizontal D_h. Their paper (Nat. Commun. 12, 6376, doi:10.1038/s41467-021-26710-0, text in `papers/h14_eddy/ashkenazy2021.txt`) has no 30 as a diffusivity:

- `ashkenazy2021.txt:402-405`: the explicit horizontal eddy viscosity and tracer diffusion are set to 50 and 5 m2/s (2D) and 20 and 2 m2/s (3D). So the tracer diffusion in the model is 5 (2D) and 2 (3D) m2/s.
- `ashkenazy2021.txt:553-559`: from the momentum budget they diagnose diffusion/viscosity coefficients "of about 300 m2/s" and say a rough lower bound for the global mean eddy mixing coefficient is "about 200 m2/s". The same passage says the viscosity is usually larger than the diffusion coefficient, so 300 is a viscosity-type number and may overstate the tracer diffusivity.

Both diagnosed numbers are above the top of the sweep's bracket (1e2). The other two sources in the comment: Zhang, Kang and Marshall 2024 (Science Advances, "Ocean weather systems on icy moons, with application to Enceladus"; `papers/h14_eddy/zhang2024_epmc.xml`) write that K "should be of order ~0.1 m2 s-1" for Enceladus; I did not find "0.03-2" in that file. Zeng and Jansen 2021 (0.1-0.25, "~1000 yr to mix a hemisphere") I did not re-check here. So the bracket 1e-2..1e2 is an assumption, as its own comment says, and one of the three sources in the comment is quoted wrongly.

## Where the quote and the bracket appear (file:line)

States "use 30" or "0.03-30":
- `uncertainty.py:38-39` ("Ashkenazy & Tziperman 2021 use 30").
- `AUTONOMOUS_LOG.md:178` ("ASSUMED bracket around the 0.03-30 literature values").
- Already flagged at `journals/new_analyses_h14_eddy.md:145`.

Relies on the bracket 1e-2..1e2 (no "30"): `uncertainty.py:43-44` (Enceladus and Europa entries), `ocean_worlds.py:82` (the 12 other moons take Europa's D_h), `README.md:605`, `PAPER_DRAFT.md:54`, `journals/psj_draft.md:227`, `journals/astrobiology_draft.md:439`, `journals/ija_draft.md:549` (Table 1 row), `journals/ija_draft.md:571` (the other moons use Europa's ranges), `journals/ija_draft.md:672` (Section 3.3, "-1 to -4 points" for faster sideways mixing).

## What I ran

Same code path as `multiwall.log` Section 4: `multiwall.run4(0, full=True)`, N = 200,000, seed 0, common random numbers. `uncertainty.logu(lo, hi)` takes the same uniform draw for any range, so only D_h differs between runs (a draw's log10 D_h + 2 scales by exactly (4 + log10(hi/1e2)) / 4; `dh_range.log:2`). The D_h range was overridden for Enceladus and Europa in `U.WORLDS` and for the 12 other moons in `ocean_worlds.DH`. Early Earth keeps its own 1e2..1e4 and does not move. Enceladus was run as logged (seed collision `[0, 0]` = `0`) and with the `seedfix.py` wrapper (the paper's Enceladus); all numbers below are the fixed Enceladus unless I say otherwise. Seeds 1-3 were run with the fix.

Self-check (`dh_range.log:4-7`): at the 1e2 bound the 15 worlds reproduce `multiwall.log:525-539` (three-way split) and `multiwall.log:542-556` (P_pved no lid / final, five rules) to the printed 0.1, and fixed Enceladus reproduces `seedfix.log:122` and `:125` (26.4 / 9.8 / 63.9, x10 19.4). Worlds 1-14 are identical between the collided and fixed runs at all bounds.

## Result 1: headline share (multiwall P_pved, main rule), `dh_range.log:13-27`

Change in points, bound 1e2 -> 3e2 and 1e2 -> 1e3, as "no lid / final". Final equals no lid for worlds without a salt lid.

| world | 1e2 level (no lid / final) | change at 3e2 | change at 1e3 |
|---|---|---|---|
| Early Earth | 99.9 / 99.9 | 0.00 | 0.00 |
| Umbriel | 53.7 | -0.65 | -1.31 |
| Oberon | 51.4 | -0.54 | -1.11 |
| Titania | 51.3 | -0.49 | -1.03 |
| Ariel | 49.6 | -0.62 | -1.25 |
| Ganymede | 44.1 | -0.18 | -0.49 |
| Callisto | 43.5 | -0.21 | -0.46 |
| Pluto | 41.1 | -0.33 | -0.73 |
| Triton | 37.0 | -0.33 | -0.69 |
| Titan | 35.8 | -0.22 | -0.46 |
| Europa | 39.6 / 30.1 | -0.29 / -0.23 | -0.65 / -0.49 |
| Dione | 42.0 / 26.2 | -0.47 / -0.30 | -0.96 / -0.60 |
| Miranda | 27.1 / 9.0 | -0.42 / -0.13 | -0.91 / -0.29 |
| Mimas | 17.3 / 5.6 | -0.28 / -0.09 | -0.61 / -0.20 |
| Enceladus (fixed) | 26.4 / 0.0 | -0.39 / 0.00 | -0.88 / 0.00 |

- Every icy world falls, by 0.18 to 0.65 points at 3e2 and 0.46 to 1.31 points at 1e3 (no lid). After the salt lid the lid worlds fall by less (Europa -0.23 / -0.49, Mimas -0.09 / -0.20) and Enceladus stays at 0.0.
- Cascade x10 (`dh_range.log:31-45`) moves by about the same amount: -0.31 to -0.53 at 3e2 and -0.60 to -1.05 at 1e3 (no lid), Enceladus 19.4 -> 19.1 -> 18.7.
- Enceladus as logged (seed collision, `dh_range.log:48-50`): 28.3 -> 27.9 -> 27.5 (change -0.4 / -0.8), almost the same change as the fixed one.
- Order of the 15 worlds by final share (`dh_range.log:71-75`): unchanged at 3e2 and 1e3. The closest pairs at 1e2 are Oberon/Titania (gap 0.09), Ganymede/Callisto (0.67), Triton/Titan (1.21), Titania/Ariel (1.70); none swaps.

## Result 2: where the lost share goes (`dh_range.log:54-68`)

At 1e3 the drop in P_pved is matched by a rise in P_chance (+0.50 to +1.64, every icy world), so the one-hand share P_pved + P_chance changes by at most 0.33 points (Umbriel, Titania; Oberon 0.32) and by 0.03 to 0.19 for Enceladus, Miranda and Mimas, where P_none rises slightly instead. The reading: a larger D_h makes the deciding layer thinner, which makes the draw weaker, not later or more mixed. I checked this only as bookkeeping on the log, not by a separate calculation.

## Result 3: seeds (`dh_range.log:77-94`)

At seeds 1, 2 and 3 the 1e2 -> 1e3 change differs from seed 0 by at most 0.10, 0.12 and 0.08 points in any icy world. The changes are not sampling noise.

## Predictions (v1, `dh_range_predictions.txt`), scored by the script (`dh_range.log:97-107`)

All 10 rules held: P1 self-check, P2 Early Earth unchanged, P3 direction (14/14 fall, 14/14 monotone), P4a size ranges (all within -1.0 at 3e2 and -2.0 at 1e3; largest are Umbriel -0.65 and -1.31), P4b Enceladus -0.39 / -0.88 and Europa -0.29 / -0.65 against central -0.4 / -0.9 and -0.2 / -0.5, P5 lid (Enceladus final 0.0; the four other lid worlds move by less after the lid), P6 order (no flip, block order kept), P7 mass goes to P_chance (14/14 rise, largest one-hand change 0.33 against limit 0.5), P8 dose-response (9/9 qualifying worlds have a 1e3/3e2 ratio between 2.00 and 2.25 against the 1.7-2.6 window), P9 seeds.

Misses I should report even though no rule failed. The per-world central estimates for the 12 other moons (not pass/fail rules) were built by scaling to level and are off for some: Ganymede predicted -0.4 / -1.0, actual -0.18 / -0.49; Callisto -0.4 / -1.0, actual -0.21 / -0.46; Titan -0.4 / -0.8, actual -0.22 / -0.46 (all about half the size); Miranda -0.3 / -0.6, actual -0.42 / -0.91 and Mimas -0.2 / -0.4, actual -0.28 / -0.61 (about 1.5 times the prediction); Umbriel and Ariel were 0.1 to 0.15 too small in size. Mean absolute error of the central estimates 0.10 points at 3e2 and 0.20 at 1e3, largest 0.22 and 0.54 (Callisto at 1e3). The prediction that x10 changes are "about the same size, slightly smaller" held only loosely: x10 changes are within about 0.3 points of the main-rule changes (largest gaps Titan 0.30, Ganymede 0.29 at 1e3), and for Ganymede, Callisto and Titan they are larger than the main-rule changes. The calibration (decile curves for Enceladus and Europa only) cannot predict which of the other moons are less sensitive; I did not investigate why Ganymede, Callisto and Titan are.

## Decision rule, applied

The rule set before the run: if direction (P3), size (P4a) and order (P6) hold, the D_h upper bound is not what the headline shares depend on; the comment is wrong about the source but the numbers do not move by more than about 2 points, and the fix is wording plus one sentence in the paper. The other branch (any change below -3 at 1e3, or a swap in a pair with a 1e2 gap of 1.2 or more) did not trigger: the largest change is -1.31 (Umbriel) and no pair swaps. So the fix is wording only.

## Proposed wording (NOT APPLIED; I edited none of these files)

`uncertainty.py:36-39`, replace the D_h comment with:

```
# horizontal D_h, icy moons: 1e-2..1e2 (ASSUMED bracket, not taken from one source).
#   Zhang, Kang & Marshall 2024 (Sci. Adv.) write K ~ 0.1 m2/s for Enceladus; Zeng & Jansen 2021 use 0.1-0.25 (not re-checked).
#   Ashkenazy & Tziperman 2021 (Nat. Commun. 12, 6376) set explicit horizontal tracer diffusion to 5 (2D) and 2 (3D) m2/s and
#   diagnose an eddy coefficient of about 300 m2/s (rough lower bound about 200 for the global mean; a viscosity-type
#   number), i.e. above this bracket. dh_range.log: moving the top to 3e2 / 1e3 lowers each icy world's share by
#   0.18-0.65 / 0.46-1.31 points (no lid) and does not change the order of the worlds.
```

A paper sentence, after the Table 1 vertical-range paragraph (`journals/ija_draft.md` near :558) or in Section 3.3 near :672:

> The horizontal range for the icy moons (10⁻² to 10² m² s⁻¹) is an assumed bracket. An eddy-resolving Europa ocean simulation diagnoses a mixing coefficient of about 300 m² s⁻¹ (Ashkenazy and Tziperman, 2021), above its top. Raising the top of the bracket to 3 × 10² or 10³ m² s⁻¹ lowers the share of every icy world by 0.18 to 0.65 or 0.46 to 1.31 points (no lid; `dh_range.py`) and leaves the order of the worlds unchanged; the early Earth ocean, with its own range of 10² to 10⁴ m² s⁻¹, does not move.

`AUTONOMOUS_LOG.md:178`: change "the 0.03–30 literature values" to "literature values of 0.1 to 300 m2/s (Zhang et al. 2024: about 0.1; Ashkenazy and Tziperman 2021: about 300, a viscosity-type diagnosis)".

## Caveats

- This tests only the upper bound. The lower bound (1e-2) is untouched, and so are the median-split D_h effect (`rerun_base.log:33,45,58,69`; paper Section 3.3 "-1 to -4 points") and the 22% / 29% win rates in `uncertainty.py`. The mechanism, a larger D_h thins the deciding layer through the assumed rule f = L/H' with H' proportional to the square root of D_h, is the one the paper already calls "not a physical result" (`ija_draft.md:672`); the new numbers are a further run of that rule, not independent evidence.
- A larger D_h also widens the wall, which helps the l > w condition in `journals/new_analyses_h14_eddy.md`; I did not combine the two.
- The 3e2 column is the one that matches Ashkenazy and Tziperman's number if their 300 is read as a tracer diffusivity; their text suggests it is a viscosity, so 3e2 is more likely an upper estimate of the right value than a central one.
- All numbers are seed 0 with 200,000 draws per world; seeds 1-3 agree within 0.12 points.

## Found (outside this unit, not fixed)

- `journals/ija_draft.md:571` says "Eleven more worlds use Europa's ranges", but `ocean_worlds.py:86-98` defines 12 such worlds (Ganymede, Titan, Callisto, Triton, Pluto, Titania, Oberon, Umbriel, Ariel, Dione, Miranda, Mimas) and the next sentences cite ocean-depth sources for all 12; `multiwall.log:525-539` has 15 worlds in all.
- `AUTONOMOUS_LOG.md:178` repeats the "0.03-30" figure (see above).
- `seedfix.log:168` and `:451-453`: `rerun_base.log:305` (P14) reads FAILED in the committed log and HELD with the seed fix (Enceladus M=1e3 1.5e-14 -> 9.0e-15). It is the seed-collision unit's issue, listed here only so it is not lost.
