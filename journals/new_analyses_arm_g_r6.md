# arm_g_r6: R6 for "Nonenzymatic RNA copying, 2-AI (2026)" (scorecard.py:61)

Written 2026-10-21. Nothing in this file has been applied to scorecard.py, README.md, journals/ija_draft.md, EXPERIMENT.md or journals/droplet_prereg.md. Proposed wording is at the end, marked "not applied".

## What ran

| file | role |
|---|---|
| `arm_g_r6_predictions.txt` | predictions P1-P13 and the decision rule, written and hashed before the results run (sha256 `79bea991b1ec5e91a3ab008b4fee5c85bece239e2e153392edf458286d8c64ce`, 2026-10-21 00:35:05 PDT; re-checked after the run, unchanged; asserted in the script) |
| `arm_g_r6.py` | the run (11 s, one core): `UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with matplotlib python arm_g_r6.py > arm_g_r6.log` |
| `arm_g_r6.log` | output; every number below is quoted from it as `arm_g_r6.log:LINE`, or from `rna_seed_scan.log:LINE` |
| `arm_g_r6.err` | stderr: two ResourceWarnings (unclosed files), nothing else |

Disclosures. (1) `python scorecard.py` was tried once before the predictions were written and failed on import (no matplotlib; ocean.py:21). (2) My first launch of the results run died on a BSD `sed -i` typo in the same shell command before Python started (no output produced). (3) Between writing `arm_g_r6.py` and its only results run I changed two things in my own script and nothing else: a variable name that shadowed the positive-control label, and the line filter that echoes `rna_seed_scan.log`. No smoke run of the unit's quantities was made. (4) The script uses `rna_copying.py`, `rna_seed_scan.py` (helpers imported), `scorecard.py` and `ocean.py` unchanged; no existing file was edited.

The model is the one in `rna_copying.py` (n = 3, c = 0.5, w = 0.02, Walton and Szostak 2017 rates, T_in = 3 uM, 1 h cold phase assumed). The primer-template duplex T is the model's unit: it has no free primer pool, so "seeded templates/primers" means an excess in T_D against T_L. Cycle length is not measured; 1 h and 24 h bracket it.

## Definitions (fixed before running)

* R6 rule (scorecard.py:9, inherit.py): follows an inherited excess >= e_need. e_need: pond 5.8e-12, lake 3.2e-14, ocean 2.8e-18 (arm_g_r6.log:2).
* Earlier definition of "followed" (rna_seed_scan.py): Delta = readout excess / s >= 2, with s the stationary linear-noise s of the readout and Omega = N_setting / m_in; e_min = 2 s / gain; letter Y if e_min <= e_need in all three settings, P in some, N in none. Judged here: it is sound as a signal-to-noise test that the readout tracks a continuously supplied feed excess. It is not a test that the system selects a hand, because the racemic state is stable for r < 1 (|lambda| < 1 in every cell, arm_g_r6.log:27-35), and e_need comes from a bifurcation ramp (inherit.py), a different decision problem. I kept it as the "weak" rule so (a) and (b) are comparable.
* Added "strong" rule for (b): weak and the seed is not LOST. Outcome at tau: GROWN if the seed ever exceeds its start, KEPT if at least half of it is left at tau, otherwise LOST. LOST is declared only when an upper bound (the envelope of the decaying complex pair) is below 0.5. The weak rule counts a seed that has decayed by e^-19 as "followed" if the remnant is still above counting noise; the strong rule does not.

## Results

### (a) Excess in the feed: letters taken from the earlier log, not re-derived

* rna_seed_scan.log:26-34 (Table 2) and :115: letter Y in all 9 cells (3 r x 3 feeds). e_min, pond 1.5e-13 to 3.6e-12 (e_need 5.8e-12); lake 8.2e-16 to 2.0e-14 (3.2e-14); ocean 2.2e-20 to 5.5e-19 (2.8e-18); e_min / e_need 0.008 to 0.62 (the log's own "x need" column). Copied into arm_g_r6.log:8-18.
* I reproduced the pieces it rests on. |lambda| at 20 mM = 0.97943, 0.98797, 0.99948 for r = 0.40, 0.56, 0.98 (arm_g_r6.log:28, 31, 34; all 9 cells within 4.6e-6 of rna_seed_scan.log Table 1, P1). Gain at a 1% feed excess, 20 mM: 3.4, 5.25, 79 (arm_g_r6.log:129; ee 0.034, 0.0525, 0.79). The response is bounded and monostable: the product excess approaches 1 only at a feed excess of 0.5 (r = 0.40, ee 0.998; arm_g_r6.log:112), 0.2 (r = 0.56, ee 0.95; arm_g_r6.log:118) or 0.1 (r = 0.98, ee 0.99; arm_g_r6.log:124).
* Setting-consistent supplement (arm_g_r6.log:138-147, not part of the earlier definition). The earlier scan combined each setting's N with a 5-200 mM feed. Using each setting's own concentration (pond 10 mM, lake 0.1 mM, ocean 1e-3 mM):
  * with T_in fixed at 3 uM: e_min / e_need = 0.06-0.12 (pond), 0.153 (lake), 0.0048 (ocean); letter Y at all three r (arm_g_r6.log:142-144);
  * with T_in scaled with the feed (3 uM x m_in / 20 mM): pond 0.05-0.12, lake 2.18, ocean 0.689; letter P at all three r (arm_g_r6.log:145-147);
  * at lake and ocean there is no copying at all: the template level exceeds the no-copying level by 9e-6 (lake) and 9e-12 (ocean) (column Q-1, same lines), and the gain is 3.00 = n, the dilute limit, independent of r. The pure-pool onset feed is 2.857 mM (arm_g_r6.log:95), 29 times the lake feed and 2860 times the ocean feed.
* Implied letter for (a): **Y** under the existing rule and the earlier setting convention (rna_seed_scan.log:115); **P** if T_in is taken to scale with the setting's concentration (arm_g_r6.log:145-147). In both cases the Y is "the readout tracks a feed excess", not "a hand is selected".

### (b) Excess only in the seeded templates, racemic feed

* Per-cycle decay of a template seed is the difference-mode eigenvalue: |lambda| = 0.97943 (r = 0.40), 0.98797 (0.56), 0.99948 (0.98) at 20 mM (arm_g_r6.log:28, 31, 34). All are complex pairs, so the seed oscillates in sign while it decays. No transient growth: the largest value of rho(t >= 1) over 1e5 cycles is 0.999993 (arm_g_r6.log:27-35, P2).
* **Outcome at tau: 54 of 54 cells LOST, 0 KEPT, 0 GROWN** (3 r x 3 feeds x 3 settings x 2 cycle lengths; arm_g_r6.log:46). 20 mM rows, arm_g_r6.log:39-44: the envelope of rho at tau is 10^-8.27 at its largest (r = 0.98, 24 h, pond) and 10^-192 or smaller in every other cell.
* **Letters for (b), weak and strong, per r and cycle length (20 mM): N in all 12 (6 weak, 6 strong)** (arm_g_r6.log:39-44, :189). The only cell with a finite smallest followed seed is r = 0.98, 24 h cycles, pond: e_min^T = 3.7e-3, against e_need = 5.8e-12, nine orders away (arm_g_r6.log:47-48).
* Nonlinear check, 20000 direct iterations of the 4-variable map, seeds 0.1, 0.5, 0.99: |ee_T| ends at 4.7e-15 to 4.8e-15 (r = 0.40), 7.7e-15 (0.56), and 3.3e-06 to 6.0e-05 (0.98) (arm_g_r6.log:65-67), consistent with the linear picture and with rna_seed_scan.log:99-102.
* Positive control: at r = 1.141 (the reading arm_g_rates.log:44 flags) |lambda| = 1.00220 and a 1e-3 seed reaches 0.5 after 2913 cycles; the classifier says GROWN (arm_g_r6.log:70). So the all-LOST result is not a classifier that cannot say GROWN.
* Sensitivity, one parameter at a time at (pond, 24 h cycles, 36,525 cycles, 20 mM): w = 0.005 or 0.05, c = 0.1 or 1.0, T_in = 0.3 or 30 uM, all three r: 21 of 21 LOST, no racemic state became unstable (arm_g_r6.log:151-171).
* How long would a cycle have to be? Break-even cycle length to hold a seed to half over tau, per-cycle parameters unchanged (a scale, not a model result, because w and the cold-phase loss would change with the cycle length; envelope, so these are the shortest lengths that could hold the seed; arm_g_r6.log:53-61, 20 mM):

  | r | pond (100 yr) | lake (1e4 yr) | ocean (1e8 yr) |
  |---|---|---|---|
  | 0.40 | 1.13 yr | 113 yr | 1.1e6 yr |
  | 0.56 | 0.76 yr | 76 yr | 7.6e5 yr |
  | 0.98 | 13.7 d | 3.7 yr | 3.7e4 yr |

  A day-long cycle (the longest tried) is far below all of them except the r = 0.98 pond entry, which is itself 14 d.
* Implied letter for (b): **N** at every measured r (0.40, 0.56, 0.98), both cycle lengths, weak and strong. A template seed grows only where r > 1, as in the 1.141 control, which the gel fit disfavours.

### (c) Does a racemic monomer pool poison copying so that copying only matters after symmetry breaking?

Objection (paraphrased): a racemic monomer pool poisons copying, so copying only matters after symmetry breaking. Tested inside the model as wrong-hand monomer terminating strands (rate ratio r) plus halving the right-hand concentration.

* Per attempt (arm_g_r6.log:76-85). Fresh 20 mM feed (m = 10 mM each hand): full-length probability with racemic pool / pure pool of the same total monomer = 0.31 (r = 0.40), 0.24 (0.56), 0.13 (0.98); against the same right-hand level alone 0.49, 0.38, 0.21. In the running reactor the per-attempt probability is pinned at 0.0406 by washout (c S = w/(1-w) = 0.0408), so the cost appears as a higher monomer level and feed needed to reach it, not as a lower S.
* Onset feed (copying just keeps up with washout), racemic against pure, same total feed (arm_g_r6.log:91-95): 6.29, 6.58, 7.58 mM against 2.857 mM, ratios 2.20, 2.30, 2.65.
* Steady reactor (arm_g_r6.log:99-104). At 5 mM, between the two onsets: racemic template level 1.8-2.3 x the no-copying level (copying essentially off), pure pool 79 x. At 20 mM: racemic 137-192 x against 627 x pure; yield ratio racemic/pure 0.610, 0.551, 0.435.
* Does an inherited feed excess survive the racemic pool? Yes (arm_g_r6.log:108-128). Gain at e = 0.01 is 3.40, 5.25, 79.0; gain >= 2.41 in all 9 (r, feed) cells (rna_seed_scan.log:26-34 gains, reproduced in P11); yield changes by under 0.7% (arm_g_r6.log:184). The wrong hand also raises the gain: with the wrong hand inert (r = 0) the gain is 1.40, at r = 0.06 1.59, at r = 0.40 3.40, at r = 0.98 115 (arm_g_r6.log:132-136). So the cross-inhibition that costs a factor 2-3 in feed is what amplifies a feed excess in this model.
* Pre-registered verdict rule: V1 (racemic onset >= 1.8 x pure) true; V2 (racemic Q >= 50 at 20 mM) true, Q = 192, 174, 137; V3 (gain >= 2 in 9 of 9) true (arm_g_r6.log:188). Result per the rule: **the objection is partly true**. Racemic monomer costs a factor of 2.2-2.7 in feed (onset 6.3-7.6 mM against 2.9 mM) and 39-57% of the yield at 20 mM, and in the feed window between the two onsets (about 2.9 to 6.3-7.6 mM) a pure pool copies and a racemic pool does not. It does not block copying at 20 mM, and it does not erase a feed excess. At the lake and ocean concentrations neither pool copies (above), so there the objection does not matter either way.
* Limit of this test: the model's poisoning is the single parameter r (wrong-hand addition rate ratio) plus chain termination, the same parameter that sets the stability of the racemic state; it is not an independent measurement of how much a racemic pool slows real copying.

## Predictions: held and failed (arm_g_r6.log:174-186)

* Held: P1 (|lambda| matches, max difference 4.6e-6), P2, P3 (54/54 LOST), P5 a, b, c, P6, P7, P8 (onsets 6.29, 6.58, 7.58 and 2.857 mM), P9, P10 (0.610, 0.551, 0.435), P11, P12 a, b, c, P13.
* **Failed: P4** (break-even cycle length, pond, 20 mM). Predicted 3.0 yr (2.2-3.9), 1.75 yr (1.3-2.3), 27 d (20-35 d); observed 1.13 yr, 0.76 yr, 13.7 d (arm_g_r6.log:177). The prediction used ln 2 / |ln lambda| and ignored the amplitude of the leading complex pair (2 |a| = 3.1, 2.5, 2.0), which makes the seed take 88, 132 and 2671 cycles rather than 33, 57 and 1300 to fall to half (arm_g_r6.log:54, 57, 60). It does not change any letter: P3 and P5 held, and the observed lengths are shorter than predicted, so the seed is held at slightly shorter cycle lengths than predicted, still far longer than 24 h.
* Informed by logs and hand estimates, so not blind to rna_seed_scan.log; blind to this script's output (stated in the predictions file).

## Letters

| | letter | source |
|---|---|---|
| (a) excess in the feed | Y (earlier convention); P if T_in scales with the feed | rna_seed_scan.log:115; arm_g_r6.log:142-147 |
| (b) excess in seeded templates, racemic feed | N (all 3 r x 2 cycle lengths, weak and strong) | arm_g_r6.log:39-44, :189 |
| racemic-pool objection | partly true (costs 2.2-2.7 x in feed, 39-57% of yield; does not block copying at 20 mM or erase a feed excess) | arm_g_r6.log:91-104, :188 |

## Proposed wording (not applied)

**scorecard.py:61 comment, and the cell.** The two readings give different letters, so I would put "P" and say why. Suggested comment tail: `R6 = Y for an excess in the feed (rna_seed_scan.log:115: tracks down to 1.5e-13-3.6e-12 in the pond; signal-to-noise only, racemic state stable), N for an excess in templates with a racemic feed (arm_g_r6.log:46: 54 of 54 cells lost, |lambda| 0.979-0.9995 per cycle)`. If the scorecard owner defines the inherited excess as arriving in the feed, the cell is Y; if it must sit in the replicators, N. The comment at scorecard.py:40 ("R6 here = shown to work from racemic or <= ~1% ee") and the header rule at scorecard.py:9 are not the same test; this analysis used the header rule.

**README.md, Part 6 RNA copying bullet, after "(0.67 for a +20% feed at r = 0.40)" (about line 501).** "A feed excess is followed down to a smallest seed of 1.5e-13 to 3.6e-12 in the pond (needed: 5.8e-12), but only as a signal-to-noise test, since the 50/50 state is stable. An excess placed in the templates with a racemic feed is not kept: it decays by a factor 0.979, 0.988 and 0.9995 per cycle at r = 0.40, 0.56 and 0.98, and is lost in 54 of 54 cases over 100 yr to 1e8 yr at 1 h or 24 h cycles; holding it in the pond would take cycles longer than 1.1 yr (r = 0.40) or 14 days (r = 0.98) (`arm_g_r6.py`). A racemic monomer pool does cost: the onset feed rises from 2.9 mM to 6.3-7.6 mM and the yield at 20 mM falls to 44-61% of a pure pool, but copying is not blocked above the onset, and a feed excess is still copied with a gain of 3.4 to 79 at 1% ee."

**journals/ija_draft.md, section 5, arm G paragraph (line 822), after "so arm G is expected to show seed-following only (a +20% feed reaches a plateau of 0.67 at r = 0.40)".** "The seed has to be in the feed: an excess placed in the templates of a racemic feed decays by 0.979 to 0.9995 per cycle and is lost in every case tested (54 of 54, cycles of 1 h or 24 h, 100 yr to 1e8 yr), whereas a feed excess of 1% is copied to 3.4% to 79% at 20 mM. A racemic monomer pool slows copying in the model (onset feed 6.3 to 7.6 mM against 2.9 mM for a pure pool; yield at 20 mM 44% to 61% of a pure pool) but does not stop it or erase a feed excess; at the 0.1 mM and 1 µM concentrations of the lake and ocean settings the model does not copy for either pool."

## Found outside the unit (not fixed)

1. rna_seed_scan.py (docstring lines 12-18, 36; line 197 `s = np.sqrt(c["su2"] / (N[n] / m_in))`): each setting's N is combined with a 5-200 mM feed, but the settings' own concentrations are 10 mM, 0.1 mM and 1e-3 mM. At the lake and ocean values the model has no copying (arm_g_r6.log:143-144, 146-147: Q-1 about 1e-5 and 1e-11, gain 3.0), and the lake letter flips from Y (T_in fixed, 0.153 x e_need) to failing (T_in scaled, 2.18 x e_need; arm_g_r6.log:142-147). rna_seed_scan.log:115 does not say this.
2. rna_seed_scan.log:105-113: four of seven registered predictions failed or partly failed (P1, P3, P4, P5), and the log's letters (Y at all r) are opposite to the registered ones (N, N, P; line 111). Already known; listed so the (a) letters are not quoted without it.
3. scorecard.py:9 against scorecard.py:40: R6 is "follows an inherited excess >= e_need" in the header and "shown to work from racemic or <= ~1% ee" in the comment. For this row the first gives Y (a), the second gives a bounded response (ee 3.4%, 5.25% and 79% from 1% in the feed; a racemic start stays racemic).
4. README.md ~line 497-501 and journals/ija_draft.md:822 do not separate a feed seed from a template seed; the template-seed case is not mentioned.
