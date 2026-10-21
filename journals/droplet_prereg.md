# NaClO3 droplet test: pre-registration plan

Draft, 2026-10-20. Not committed. Not part of EXPERIMENT.md. Power numbers: `droplet_power.py`, `droplet_power.log`. Decision rule
and its power: `droplet_sim2.py`, `droplet_sim2.log`. Misread handling and c*: `misread_sim.py`, `misread_sim.log`.

This is the plan, not the pre-registration itself. Before the main run, Sections 2, 3 and 6 are copied into
`droplet_predictions.txt`, hashed and committed, as was done for `blind_jang.py` (Section 8).

## 1. The question (one sentence)

As sodium chlorate droplets get larger, does the share that ends with one hand **fall gradually**, as the paper's
counting rule says, or does it **drop off at a sharp critical size**?

- **Counting rule:** a droplet ends single-handed only if every independent nucleus shares one hand. The number of
  nuclei grows with volume, so P_single falls smoothly. Every counting form has a free-width exponent h ≤ 1.64:
  - first nucleus plus Poisson extras, P = exp(−V/V0): h = 1.40, so P falls from 0.8 to 0.2 over a factor of 7.2 in V;
  - all nuclei at once: h = 1.64, a factor of 5.4;
  - nucleation through the surface: h = 0.94–1.09.
- **Critical size** (Jang et al., 2024; Bak et al., 2023): below a fixed volume a droplet is single, above it mixed.
  That is a step, h → ∞. h = 3 is the dividing line used before: P then falls from 0.8 to 0.2 over a factor of 2.52.
- **Why the published data could not decide:** both studies were tested blind and both were NOT SCORABLE on this
  point. Bak had 3 volumes. In Jang no 50 µm size bin had more than 47% mixed droplets. The data never reached the
  part of the curve where the two pictures differ. This design puts droplets there on purpose.

## 2. Predictions (fixed before the main run)

V is the pipetted volume. "Single" = L only or D only. "Mixed" = both hands in one droplet.

- **D1, direction.** Mixed droplets are larger than single ones (one-sided Mann-Whitney on V, p < 0.05).
- **D2, no sharp critical size (the main test).** The main run gets one label, from the first test it meets
  (`droplet_sim2.stats` and `classify`, fitted with the measured misread rate as in `misread_sim.refit`; Section 10):
  1. **NOT SCORABLE** unless there are at least 10 mixed and 10 single droplets, one volume level has P_single ≥ 0.6
     and one has P_single ≤ 0.4. Also NOT SCORABLE if more than 10% of the readable control drops are read wrong
     (Section 6).
  2. **FLAT** if the best-fitting curve below changes by less than 0.4 across the levels run.
  3. **RISING** if P_single rises with V (Cochran-Armitage trend z > 1.645).
  4. **STEP** if L < c* = −2.2055, otherwise **SMOOTH**.

  Two curves are fitted to the single and mixed counts at each level (binomial likelihood), both through
  P_obs = ε̂ + (1 − 2ε̂) P. ε̂ = (k + 0.5)/(n + 1) is the misread rate measured on the control drops: k of the n
  readable controls read wrong (Section 6). The misread term is measured, not assumed (Section 3).
  - **COUNT:** P = F(μ), μ = μ0 (V/V0_pilot)^a with a in [0.05, 5/3], and F one of exp(−μ/2) (first nucleus plus
    Poisson extras), 2/(exp(μ/2) + 1) (all nuclei, drops with none excluded) or 2 exp(−μ/2) − exp(−μ) (drops with
    none crystallize once at dryness).
  - **STEP:** P = Φ(−ln(V/Vc)/s), Φ the normal distribution function, s in [0.02, 0.6]. s = 0.6 is about the old
    h = 3 line (h ≈ 1.76/s).

  L = best log-likelihood of COUNT − best of STEP. c* is the lowest 5th percentile of L over simulated counting
  curves (87 curves × true misread rates 0–10% × pipette CV 5 or 15%, each run fitted with its own ε̂ from 128
  controls; 12 conditions in all; 1000 runs per curve and condition; only scorable runs with range ≥ 0.4 and
  z ≤ 1.645 count, and curves with under 50 such runs are skipped). On this grid c* = −2.2055, set on a calibration
  seed and frozen (`droplet_am2_gate.log:4,6`, `droplet_range_a.log:5`); `misread_sim.log` printed −2.21. The
  prediction is **SMOOTH**.
- **D3, chance hand.** Among single droplets, L and D fit 50/50 (two-sided binomial p > 0.05).
- **D4, hand does not depend on size.** The volumes of L-only and D-only droplets do not differ (two-sided
  Mann-Whitney p > 0.05).
- **Reported, not predicted:** V0 (fit P = exp(−V/V0)), L, â and ŝ, the earlier fit's Vk and h with a 95% profile
  interval (`droplet_power.fit_h`), P_single per level, the L/D split per slide, the number of visible crystals per
  droplet, and the control-drop results (Section 6).

## 3. Readings fixed in advance

| Outcome | Reading |
|---|---|
| D2 SMOOTH | HELD. No sharp critical size. This fits the counting rule, **or** a critical size that varies from droplet to droplet by more than about 1.35–1.6× (Section 10). It does not prove counting. |
| D2 STEP | FAILED. A sharp critical size decides. The counting rule's breadth is wrong for these droplets, and the paper says so in Section 3.7. |
| D2 FLAT | NOT DECIDED. P_single changes by less than 0.4 over the levels run: no sharp step inside the range, but too shallow to compare the two pictures. Reported as found. |
| D2 RISING | NOT DECIDED. P_single rises with V, which neither picture predicts (D1 then likely fails too). |
| D2 not scorable (after the extension rule, or over 10% misreads) | NOT SCORABLE. Reported as such; no more levels are added. |
| D1 fails | Larger droplets are not mixed more often, so the number of deciding nuclei does not grow with volume. |
| D3 or D4 fails | Something steered the hand (dust or a seed crystal of one hand), possibly by size. |

**Misreads favour the prediction.** A misread droplet pulls the ends of the curve away from 1 and 0, and that makes
a true step look broader. Both D2 fits therefore use the misread rate measured on 128 control drops (Section 6), and
c* was set with that refit (`misread_sim.log`: 500 runs per case, pipette CV 5 or 15%, true misreads 0–10%):
- Counting curves were called STEP in at most 5.4% of runs at any misread rate.
- Steps with s ≤ 0.15 were called STEP in 91–100% of runs at misreads up to 6%, and 81–97% at 8–10%. Steps with
  s = 0.2: 84–100% and 68–91%. These figures leave out the control cap (Section 2, item 1). With the cap, a NOT
  SCORABLE run is not a STEP, and the figures are 87–100% up to 6% and 33–78% at 8–10% for s ≤ 0.15, and 80–100% and
  28–69% for s = 0.2 (`droplet_ctrl_am2.log:155-158`).
- **Above 10%** measured misreads, the largest rate simulated: D2 is NOT SCORABLE. The cap acts on the measured
  rate, so it does not switch at a true rate of 10%. With 128 controls, D2 is NOT SCORABLE in 4.4% of runs at a true
  rate of 6%, 22.4% at 8% and 51.9% at 10% (`droplet_range_cap.log:6,8,10`). At a true 11%, 33.9% of runs are still
  scored, and at 12%, 22.2% (`droplet_ctrl_am2.log:6`). Rates above 10% were not simulated, so the false-STEP rate in
  those runs is not known.

A SMOOTH at a measured rate above 6% is reported with its weaker power.

## 4. Method

Sessile drops of NaClO3 solution pipetted onto a hydrophobic surface (Parafilm on a glass slide) and left to
evaporate, the method of Bak et al. (2023), the closest published setup that needs no microfluidics. Each drop's
hand is read under crossed polarizers with the analyzer slightly rotated, as in `viedma_plan.md` Section 6.

**V is set by the pipette, not measured from images.** A volume error blurs a step: a step seen through a log-volume
error σ looks like h ≈ 1.7/σ, so σ = 0.3 would already give h ≈ 6. A pipette with a 5–10% CV keeps a true step at
s ≈ 0.05–0.10 (h ≈ 17–34), well inside the steps the rule catches (s ≤ 0.15, Section 3).

## 5. Pilot first (not counted, not pre-registered)

- **Volume ladder:** 0.25, 0.5, 1, 2, 4, 8 and 16 µL (doubling), 20 drops each, 140 drops. Bak's three points
  give V0 ≈ 4–5 µL. The ladder covers P_single 0.8 → 0.2 for any V0 from about 1 to 10 µL.
- **Fit V0** from P = exp(−V/V0) (`blind_units.exp_fit`). In simulation (counting rule, pipette CV 5%, true V0 of
  1, 4 or 10 µL), the fitted V0 fell within 0.73–1.38 × the true V0 in 95% of pilots. The main-run power (Section 7)
  allows the true V0 to be anywhere from 0.5 to 2 × V0_pilot.
- **Fix, then do not change:**
  - solution concentration and evaporation temperature;
  - when the hand is read (for example at full dryness, or a fixed time after crystals appear);
  - drops per slide and the spacing between them;
  - what counts as unreadable.
- **Gates before the main run** (all must pass, or fix and repeat the pilot):
  - pipette CV ≤ 10%, by weighing 20 dispenses of water at each level on the 0.001 g scale;
  - blind re-read of 40 pilot drops disagrees on ≤ 3%;
  - 16 control drops made as in Section 6 (8 of each kind), scored blind among pilot photos, at most 1 read wrong;
  - at least 90% of drops readable;
  - no slide where all single drops share one hand (a sign of seeding).
- Pilot data are never used in D1–D4.

## 6. Main run

- **Levels:** K = 8 volumes, log-spaced from 0.2 to 1.6 × V0_pilot (ratio 1.35 between levels). Under
  P = exp(−V/V0) that runs from P_single 0.82 to 0.20.
- **Droplets:** 25 per level, 200 in all (Section 7), plus 128 control drops. The extension rule can add up to 50 more.
- **Layout:** every slide carries every level, in an order drawn by a script with a fixed seed, so slide, position
  and time are not tied to volume.
- **Extension rule (fixed now):** after the main run, if no level has P_single ≤ 0.4, add one level at 1.35× the
  top volume; if none has P_single ≥ 0.6, add one at the bottom (the top first if both apply). The check runs twice,
  so at most two levels are added in all. Same droplet count per new level. Only P_single per level is looked at to
  make this decision, not h.
- **Blinding:** each drop is photographed under its ID. A script renames the photos to random codes and crops each
  to its drop, rescaled to the same frame, so the reader cannot see the volume. The reader scores L, D, mixed or
  unreadable. The key is joined only after every photo is scored.
- **Control drops (the misread rate):** at every level, 8 extra drops get one seed crystal of known hand and 8 get
  one crystal of each hand, 128 in all. Seeds are pilot crystals whose hand both readers agree on. Control drops
  use saturated solution, so the seed does not dissolve, and they are placed in the drawn order on the same
  slides. They are photographed, renamed and scored blind with the rest; only the key marks them. The misread
  rate is the share of controls read wrong for single or mixed (a one-seed drop read as mixed, a two-seed drop
  read as single). Both D2 fits use it, and above 10% D2 is NOT SCORABLE (Sections 2 and 3). A one-seed drop
  read with the wrong hand is reported separately (it bears on D3 and D4). Unreadable controls are reported, not
  counted. If fewer than 128 controls can be read, c* is −2.61, the value simulated for 64 controls
  (`misread_sim.log:4`), and with fewer than 64 readable controls D2 is NOT SCORABLE. The 10% cap is on the readable
  controls. A control can grow a stray crystal of its own, which looks like a misread, and a seed crystal may be
  easier to read than one that formed by itself, so the estimate can err either way (Section 10). Controls are
  never used in D1–D4.
- **Re-read:** a second reader re-scores a random 20% of photos blind. The disagreement rate is reported.
- **Exclusions, logged before unblinding:** drops with no crystal, unreadable drops, drops that merged with a
  neighbour or were touched. Any slide whose single drops are lopsided at two-sided binomial p < 0.001 is excluded
  as seeded and reported.

## 7. How many droplets (`droplet_power.log`)

**This table was made for the earlier h ≤ 3 rule.** For the rule now in Section 2, `misread_sim.py` ran the same
protocol at 25 per level with 128 control drops, 500 runs per truth, at pipette CV 5 or 15% and true misreads 0–10%.
Counting curves were called STEP in at most 5.4% of runs, and steps with s ≤ 0.15 in 91–100% at misreads up to 6%
(81–97% at 8–10%). With the control cap of Section 2, the step figures are 87–100% and 33–78%
(`droplet_ctrl_am2.log:155-156`). The count of 25 per level stands.

The simulation runs the whole main run 400 times per case: 8 levels, drawing droplets, applying the extension rule
and the scorability rule, then the D2 fit. A run counts as **called right** if it is scorable and gives h ≤ 3 under
a counting form, or h > 3 under a step.

- **Truths:** both counting forms; a sharp step; and a step blurred to h = 6.
- **True V0:** 0.5, 0.7, 1, 1.4 or 2 × V0_pilot.
- **Conditions:** pipette CV 5% with no misreads, CV 5% with 3% misreads, and CV 10% with 6% misreads.

Each cell is the **worst** share called right over every true V0 and every condition.

| Droplets per level (total) | Counting forms | Sharp step | Step blurred to h = 6 | Worst |
|---|---|---|---|---|
| 5 (40) | 0.71 | 0.00 | 0.86 | 0.00 |
| 10 (80) | 0.91 | 0.93 | 0.94 | 0.91 |
| 15 (120) | 0.94 | 0.99 | 0.94 | 0.94 |
| **20 (160)** | **0.97** | **0.99** | **0.97** | **0.97** |
| 30 (240) | 0.99 | 1.00 | 0.99 | 0.99 |

- **The choice.** 20 per level is the smallest count that is called right at least 95% of the time in every case.
  Pipetting **25 per level** leaves 20 even if a fifth of the drops are excluded.
- **Why 5 per level fails.** If the true V0 is 2 × V0_pilot, a step sits above all but the top level. That level
  gives only 5 mixed droplets, short of the 10 that scoring needs, and the extension rule does not fire because
  one level already has P_single ≤ 0.4.
- **Why the fit has a misread term.** Take a step at 2 × V0_pilot, CV 5%, 20 per level, read with the earlier h ≤ 3
  line (`droplet_power.py`). Without a fixed 3% term, the step is called right in 84% of runs at 3% misreads and 30%
  at 6% misreads. With it, the figures are 100% and 99% (`droplet_power.log:8-9`). The fit in Section 2 uses the
  measured rate instead of a fixed 3%.
- **Reported, not predicted.** At 20 per level with 3% misreads, the fitted h under the counting forms falls between
  0.85 and 2.64 in 95% of runs (median 1.44 for first-plus-extras, 1.59 for all-at-once). The 3 cut-off sits above
  that range.
- **Not simulated here.** The pilot precision assumes the counting rule. Evaporation is not in this table; it is
  simulated separately in `droplet_sim.py` and changes the picture (Section 10).

## 8. Pre-registration steps

1. Pilot passes every gate in Section 5. V0_pilot is written down.
2. `droplet_predictions.txt` is written: Sections 2, 3 and 6 with V0_pilot, the eight volumes and the droplet count
   filled in. Its sha256 goes in the docstring of `droplet_test.py`, which asserts it. Both are committed **before
   the first main-run drop is pipetted**.
3. Main run with control drops, blind scoring, re-read. The misread rate from the controls is computed before D2.
4. Outcomes are entered in `droplet_test.py` in a later commit, without changing the predictions file. Every
   outcome is reported as found.

## 9. Kit, safety and forms

- **Kit beyond `viedma_plan.md` Section 9:** an adjustable micropipette covering about 0.2–20 µL (two pipettes may
  be needed) and tips, Parafilm, glass slides, a slide box to keep dust off while drops dry.
- **Safety:** NaClO3 is a strong oxidizer. The rules in `viedma_plan.md` Section 10 apply in full: keep it wet,
  away from paper, sugar, cloth and sulfur; wipe spills wet. Dried drops are solid NaClO3 on Parafilm, so
  discard slides as chemical waste, not with paper.
- **Forms:** the work may not be done at school (ruling of 2026-10-16), and a qualified Direct Supervisor must
  be present for all handling. If this runs as part of a science-fair project, it must be in the Research Plan that
  the SRC approves before any experiment, including the pilot.

## 10. Limits

- **HELD does not prove counting.** A sharp critical size that varies between droplets also gives a broad curve.
  A step blurred by a log-normal spread σ in the critical volume is a STEP curve with s ≈ σ (more with pipette
  error). The earlier fixed-3% rule called STEP in 84% of runs at s = 0.3 and 17% at s = 0.45 (CV 5%, 3% misreads,
  `droplet_sim2.log`). The rule in Section 2 calls s = 0.3 STEP in 69–98% at misreads of 0–4% (CV 5%); s = 0.45 was
  not rerun. So SMOOTH needs σ above about 0.3–0.45, a critical volume varying by about 1.35–1.6× (one standard deviation) between drops made
  the same way. The test rules out a critical size that is sharp *and* the same for every drop.
- **FAILED is the strong outcome.** Under the old h ≤ 3 line that held only in drops that do not evaporate (next
  point). Under the rule in Section 2, no counting form in Section 1 and none of the drying models below was called
  STEP in more than 5.4% of simulated runs at any misread rate from 0 to 10% (`misread_sim.log`). A counting curve
  with a floor (a share of drops single at every volume, from a seed crystal or a favoured site) is outside that
  list. It was called STEP in up to 10.0% of runs at zero true misread and CV 5%, and up to 2.8% at 3%
  (`droplet_range_a.log:15,16`).
- **Evaporating drops can make counting look sharp** (`droplet_sim.py`, `droplet_sim.log`, 2026-10-20). This is a
  model study, not data. Predictions were hashed before the run (`droplet_sim_predictions.txt`).
  - **The model.** A drop dries, and a larger drop takes longer to dry (time ∝ V^(2/3)). Nuclei form in the bulk
    or at the contact line, and each crystal may stop nucleation near it (quench). Thirty-eight models were run,
    each through 200 simulated runs of Sections 5–7.
  - **Nuclei in the bulk.** The count grows like V^(5/3), not V. The true curve has h = 2.66, and the protocol
    returns FAILED in 34% of runs.
  - **Drops that never nucleate.** If they crystallize once at dryness and count as single, rather than being
    excluded, the true h is 3.54 and FAILED comes up in 83% of runs. In the bulk model, 78% of drops at 0.2 × V0
    never nucleate before drying, so the exclusion rule would remove most of the bottom level.
  - **Nuclei at the contact line.** The count grows like V: perimeter V^(1/3) times drying time V^(2/3). Then
    h = 1.59, and every run is HELD. This replaces the earlier guess of h ≈ 0.5, which left out drying time.
  - **A wide quench** behaves like no quench. Bulk models have h 2.1–2.5 and FAILED in up to 28% of runs (up to
    68% with the dryness rule).
  - **A narrow quench with a steep supersaturation ramp** flattens the curve (h 0.1–1.5), and up to 16% of runs
    are NOT SCORABLE. In 3 of 36 quench models P_single rises with V, so D1 would fail. In 2 models no rate gives
    P_single = 1/2 at 3.1 µL.
  - **Predictions.** Four of six held. M4 failed: a steep ramp did not always flatten the curve, and with a wide
    quench h reached 2.52. M6 failed in the 3 models above.
  - **Consequence.** In evaporating drops a FAILED could come from counting itself, so h ≤ 3 is not a fair line
    there. Section 2 now uses a rule that allows nuclei ∝ V^a with a up to 5/3 (next point). Drops that do not
    evaporate (fixed volume, cooled or antisolvent, as in Jang et al.) remain the cleaner option.
- **A rule that allows for drying** (`droplet_sim2.py`, `droplet_sim2.log`, 2026-10-20). Also a model study;
  predictions were hashed before the run (`droplet_sim2_predictions.txt`). Three agents (physicist, statistician,
  skeptic) designed the rule in Section 2 over two rounds. c* (then −1.61, with a fixed 3% misread term) was
  calibrated on every counting curve: the three
  forms at a = 1/3 to 5/3, and every `droplet_sim.py` model, with drops that never nucleate either excluded or
  counted single at dryness. Every result below is from a fresh seed.
  - **Six of eight predictions held.** Counting curves were called STEP in at most 6.4% of runs in any condition,
    where the old line called the worst of them (dry, a = 5/3) FAILED in 85%. Steps with s ≤ 0.2 were called STEP
    in every run at pipette CV 5% and 3% misreads.
  - **S2 and S7 failed, for one reason.** At 6% misreads, twice the fit's 3%, a step with s = 0.15 was caught in
    only 78% of runs at pipette CV 15%. Hence the control drops (Section 6) and the refit with the measured rate
    (below).
  - **Blind spot (reported, not predicted).** A sharp step down to a plateau at P_single = 0.3, rather than to 0,
    was called STEP in only 4% of runs at s = 0.15 (5% under the earlier rule A; `droplet_range_a.log:82,53`). At
    s = 0.1, rule B called a step on a plateau h STEP in 50–82% of runs at h = 0.2, 21–30% at h = 0.25 and 3–5% at
    h = 0.3 (four noise conditions, 600 runs each; `droplet_am2_full.log:19,81-100`). A SMOOTH from rule B does not
    exclude a step on a plateau of 0.2 or more. The rule looks for a step to full mixing only, and rule B stays.
  - **Plateau fit AM2 (POST HOC sensitivity, not predicted).** Fitting a plateau h in both families, with
    c* = −2.9763 set on the full 12-condition grid (`droplet_am2_gate.log:7`), calls the step at h = 0.3, s = 0.1 STEP
    in 70–98% of runs (`droplet_am2_full.log:81-100`). Its false STEP is at most 2.8% and 4.2% in the two simulation
    sets (`droplet_am2_full.log:171,179`). It loses up to 7.2 points (with the cap) and 10.4 points (without it) of
    STEP share on steps at s = 0.3 (`droplet_am2_full.log:303-304`), so rule B stays. On a Bak-like smooth series at
    misreads of 6% or more, AM2 rules out a sharp step in 8% of runs, against 37% for B
    (`droplet_am2_gate.log:114-115`), so a SMOOTH under AM2 says little against a sharp step. If AM2 is ever used in
    place of B, c* must be recalibrated at the control count in use (−2.7980 at 480 controls, `droplet_am2_gate.log:9`).
  - **Exclusion zones cannot fake a step** (`exclusion_bak.py`, `exclusion_bak.log`; POST HOC). If each nucleus
    blocks later ones near it, P_single depends on drop size through one ratio, which in drying drops grows like
    V^e with e ≤ 1/2. Over 27 such curves (nuclei in the volume, on the base or at the contact line; three zone
    growth rates; e = 1/6, 1/3, 1/2) the steepest had h = 1.58, and none was ever called STEP. At e = 1/6 most runs came out FLAT, so in drops like
    that D2 would not be decided.
  - **Bak et al. (2023), re-fit (POST HOC).** Their 90 °C series is SMOOTH for any guess of the unprinted drop
    counts at P = 1 (3, 6 or 10), and a step with s ≤ 0.15 is rejected (2 Δlog-likelihood 14–17). Counting and
    exclusion fit it equally well (within 0.4), so it rules out a sharp step but cannot choose between those two.
    Their main series is NOT SCORABLE (5 mixed drops).
- **Misreads measured, not assumed** (`misread_sim.py`, `misread_sim.log`, 2026-10-20). Also a model study;
  predictions were hashed before the run (`misread_sim_predictions.txt`). Two policies at true misread rates 0–10%,
  pipette CV 5 or 15%, and 64, 128 or 192 control drops: A, the earlier plan (fixed 3% term, c* = −1.61, NOT
  SCORABLE when over 6% of controls are read wrong), and B, the rule now in Section 2 without its control cap (both
  fits use the measured rate; c* set per control count on a calibration seed: −2.61 at 64, −2.2055 at 128, −2.03 at
  192; the 64 and 192 values are to two decimals). 500 runs per case on a fresh seed.
  - **Four of seven predictions held.**
  - **The fixed 3% term failed when misreads are few (M1).** With no misreads, a drying model (nuclei in the bulk,
    no quench, drops that never nucleate counted single) was called STEP in 21% of runs at CV 5% (17% at CV 15%),
    and in 10% at 2% misreads. `droplet_sim2.py` never ran below 3%.
  - **The tiers threw out the hard cases.** At 6% misreads and 64 controls, 51–54% of runs with a step at s = 0.15
    were NOT SCORABLE, and the step was called in only 37–42%.
  - **The refit fixes both.** Counting curves were called STEP in at most 5.6% of runs at any rate and control
    count (5.4% at 128). A step with s = 0.15 at 6% misreads was called in 86–94% with 64 controls, against
    37–42% under A (M5).
  - **M4 and M6 failed.** With 64 controls the s = 0.15 step fell to 86% (CV 15%, 6% misreads), short of the 90%
    predicted; hence 128 controls (91%). Going from 64 to 192 controls gained 0.06, not at most 0.05.
  - **The cost.** Broad steps are caught less often when misreads are few: s = 0.3 at 2% misreads and CV 5% was
    called STEP in 83% of runs (B, 128 controls) against 90% (A, 64).
  - **Not simulated.** Misreads are taken as equally likely both ways, and control drops as misread as often as
    the rest; a seed crystal may be easier to read, which would make the measured rate too low. c* was set for 128
    readable controls; fewer would make ε̂ noisier. Rates above 10% were not run, hence NOT SCORABLE there.
  - **Bak et al. (2023) under this rule** (`bak_policy_b.py`, `bak_policy_b.log`; POST HOC, their data were read
    long before). Bak had no control drops, so the misread rate was tried at 0–10%, and the unprinted drop counts at
    P = 1 at 3, 6 and 10. Their 90 °C series is SMOOTH in every case (L from +8.1 down to −0.02, all above c*). A
    sharp step (s ≤ 0.15) is rejected (2 Δlog-likelihood 5.8–283) at misread rates up to 6%, at 8% for two of the
    three count guesses (3.0–4.8), and not at 10% (1.2–2.6). Their main series is NOT SCORABLE at every rate. No
    unread published dataset was found that could test the rule blind (searched 2026-10-20).
- **One salt, one method.** A result here says nothing yet about the antisolvent droplets of Jang et al.
