# NaClO3 droplet test: pre-registration plan

Draft, 2026-10-20. Not committed. Not part of EXPERIMENT.md. Power numbers: `droplet_power.py`, `droplet_power.log`.

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
- **D2, no sharp critical size (the main test).** Fit P_obs = 0.03 + 0.94/(1 + (V/Vk)^h) to the single and mixed
  counts at each level (binomial likelihood, h in [0.1, 50], `droplet_power.fit_h`). The fixed 3% term stands for
  misread droplets (Section 3). The prediction is **h ≤ 3**. Scorable only if there are at
  least 10 mixed and 10 single droplets, one volume level has P_single ≥ 0.6 and one has P_single ≤ 0.4.
- **D3, chance hand.** Among single droplets, L and D fit 50/50 (two-sided binomial p > 0.05).
- **D4, hand does not depend on size.** The volumes of L-only and D-only droplets do not differ (two-sided
  Mann-Whitney p > 0.05).
- **Reported, not predicted:** V0 (fit P = exp(−V/V0)), Vk and h with a 95% profile interval, P_single per level, the
  L/D split per slide, and the number of visible crystals per droplet.

## 3. Readings fixed in advance

| Outcome | Reading |
|---|---|
| D2 scorable, h ≤ 3 | HELD. No sharp critical size. This fits the counting rule, **or** a critical size that varies from droplet to droplet by more than about 1.8× (Section 10). It does not prove counting. |
| D2 scorable, h > 3 | FAILED. A sharp critical size decides. The counting rule's breadth is wrong for these droplets, and the paper says so in Section 3.7. |
| D2 not scorable after the extension rule | NOT SCORABLE. Reported as such; no more levels are added. |
| D1 fails | Larger droplets are not mixed more often, so the number of deciding nuclei does not grow with volume. |
| D3 or D4 fails | Something steered the hand (dust or a seed crystal of one hand), possibly by size. |

**Misreads favour the prediction.** A misread droplet pulls the ends of the curve away from 1 and 0, and that makes
a true step look broader. The fit therefore allows a fixed 3% misread rate at both ends (D2). So D2 is NOT SCORABLE
if the estimated misread rate, half the blind re-read disagreement (Section 6), is above 6%, the largest rate simulated.

## 4. Method

Sessile drops of NaClO3 solution pipetted onto a hydrophobic surface (Parafilm on a glass slide) and left to
evaporate, the method of Bak et al. (2023), the closest published setup that needs no microfluidics. Each drop's
hand is read under crossed polarizers with the analyzer slightly rotated, as in `viedma_plan.md` Section 6.

**V is set by the pipette, not measured from images.** A volume error blurs a step: a step seen through a log-volume
error σ looks like h ≈ 1.7/σ, so σ = 0.3 would already give h ≈ 6. A pipette with a 5–10% CV keeps a true step at
h ≈ 17–34, far above 3.

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
  - at least 90% of drops readable;
  - no slide where all single drops share one hand (a sign of seeding).
- Pilot data are never used in D1–D4.

## 6. Main run

- **Levels:** K = 8 volumes, log-spaced from 0.2 to 1.6 × V0_pilot (ratio 1.35 between levels). Under
  P = exp(−V/V0) that runs from P_single 0.82 to 0.20.
- **Droplets:** 25 per level, 200 in all (Section 7). The extension rule can add up to 100 more.
- **Layout:** every slide carries every level, in an order drawn by a script with a fixed seed, so slide, position
  and time are not tied to volume.
- **Extension rule (fixed now):** after the main run, if no level has P_single ≤ 0.4, add one level at 1.35× the
  top volume; if none has P_single ≥ 0.6, add one at the bottom. Same droplet count per new level, at most two per
  side. Only P_single per level is looked at to make this decision, not h.
- **Blinding:** each drop is photographed under its ID. A script renames the photos to random codes and crops each
  to its drop, rescaled to the same frame, so the reader cannot see the volume. The reader scores L, D, mixed or
  unreadable. The key is joined only after every photo is scored.
- **Re-read:** a second reader re-scores a random 20% of photos blind. The disagreement rate is reported; half of it
  estimates the misread rate, which sets the 6% rule in Section 3.
- **Exclusions, logged before unblinding:** drops with no crystal, unreadable drops, drops that merged with a
  neighbour or were touched. Any slide whose single drops are lopsided at two-sided binomial p < 0.001 is excluded
  as seeded and reported.

## 7. How many droplets (`droplet_power.log`)

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
- **Why the fit has a misread term.** Take a step at 2 × V0_pilot, CV 5%, 20 per level. Without the 3% term, the
  step is called right in 84% of runs at 3% misreads and 30% at 6% misreads. With it, the figures are 100% and 99%.
- **Reported, not predicted.** At 20 per level with 3% misreads, the fitted h under the counting forms falls between
  0.85 and 2.64 in 95% of runs (median 1.44 for first-plus-extras, 1.59 for all-at-once). The 3 cut-off sits above
  that range.
- **Not simulated.** The pilot precision assumes the counting rule. Contact-line nucleation (h ≈ 0.5; Section 10) is
  not simulated either.

## 8. Pre-registration steps

1. Pilot passes every gate in Section 5. V0_pilot is written down.
2. `droplet_predictions.txt` is written: Sections 2, 3 and 6 with V0_pilot, the eight volumes and the droplet count
   filled in. Its sha256 goes in the docstring of `droplet_test.py`, which asserts it. Both are committed **before
   the first main-run drop is pipetted**.
3. Main run, blind scoring, re-read.
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
  A step blurred by a log-normal spread σ in the critical volume looks like h ≈ 1.7/σ, so h ≤ 3 needs σ ≥ 0.57, a
  critical volume varying by about 1.8× (one standard deviation) between drops made the same way. The test rules
  out a critical size that is sharp *and* the same for every drop.
- **FAILED is the strong outcome.** A fitted h > 3 cannot come from counting under any form in Section 1.
- **Evaporating drops are not Jang's drops.** Supersaturation rises as the drop dries, and nucleation may start at
  the contact line. If the nucleus count grows like the perimeter (V^(1/3)) rather than V, the fall is broader
  still (h about 0.5), which is fine for D2 but may need the extension levels to reach P_single ≤ 0.4.
- **One salt, one method.** A result here says nothing yet about the antisolvent droplets of Jang et al.
