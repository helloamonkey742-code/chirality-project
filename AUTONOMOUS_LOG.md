# Autonomous log — chirality project

## 2026-09-23 — setup (interactive session, not an autonomous wave)
- Project moved to ~/Documents/chirality-project, local git repo created (baseline f6e4ea5).
- Done and self-checked: Parts 1–7 of the README (sim, ocean, domains, patches, coarsen3d, k3, scenarios,
  inherit, dual, scorecard, design, frank). `./run_all.sh` fast set passed 9/9; the 3D runs and frank.py passed separately.
- Fact-check corrections already applied: Enceladus mixing source is Zeng & Jansen 2021 (5e-5 m²/s in
  simulations, plausible 1e-10–1e-3), not "Kang"; the README was updated.
- In flight at hand-off: uncertainty.py (ran: Earth 100%, Enceladus/Europa 18%, mixing dominant),
  frank2.py and shell3d.py (logs may or may not be complete). → wave (a).
- Next recommendation: wave (a).

## 2026-09-23 — interactive session note (not an autonomous wave)
- Wave (a) is being done by the interactive session. **Do not start wave (a).** If README still has no Part 8 when
  you read this, AND `shell3d.log` ends in `exit 0` or `exit 1`, AND the interactive session has not
  committed, then take it over.
- MEASURED: frank2.py passed (F = 0.71, same as scheme 1). The first shell3d run (256² × 8, one seed) FAILED its
  majority assert (0.52 → 0.44). It's being re-run at 512² × 8 with 4 seeds + a 50/50 control, because 256 held only ~8 patches.

## 2026-09-23 00:30 — wave b (literature firm-up)
- Why: wave (a) is still being done interactively (shell3d.log was empty and still running). Load was 10.3 and a
  train_grasp job was running, so this wave did literature work only, with no simulations.
- Prediction: all three would be confirmed from primary abstracts. The "revised downward" clause might have no source.
- Result (sonnet builder + separate sonnet critic, both agreed; all DOCUMENTED):
  - DNA PVED sign: primary source is Faglioni, D'Agostino, Cadioli & Lazzeretti 2005, CPL 407:522,
    doi:10.1016/j.cplett.2005.04.009. The abstract says the weak force does "not favor" natural double helices.
    CONFIRMED, but it is only "not the natural helix", not a strong mirror-helix preference. Full text is paywalled,
    so the solvent/effect size are unchecked. None of the 21 citing papers reverses it.
  - Ozturk 2023: ~60% ee and ~25% → enantiopure are both CONFIRMED in the PMC10246896 full text. The paper does argue
    for "statistically uniform" magnetization on a HEMISPHERE scale. Critic caveat: flipping the magnet flips the
    hand, so the hemispheres may disagree (ASSUMED inference, labelled as such in the README).
  - Stribling & Miller 1987 (PMID 2819806): amino acids 3×10⁻⁴ M and HCN 4×10⁻⁶ M CONFIRMED from the abstract.
    CORRECTION: "later work revised it downward" had no source and none was found. Cleaves 2008 and Aubrey 2009 do
    not give a lower number, and the Chyba & Sagan 1992 text was inaccessible. The claim is withdrawn in the README,
    with a dated correction line. No numbers changed.
- Failed/unfinished: the Faglioni full text and Chyba & Sagan 1992 are behind paywalls.
- Next: wave (a) if the interactive session hasn't finished it (check README for Part 8 + shell3d.log exit line);
  otherwise wave (c) when load < 6.

## 2026-09-23 — wave (a) DONE (interactive session)
- README Part 8 added: uncertainty sweep (Earth 100%, Enceladus/Europa 18%; P(win | heals) 72% / 86%; mixing
  dominant), frank2 (F = 0.71, matches within 1.7 pts), thin shell 512²×8 × 4 seeds (control 49.8 ± 1.4;
  52 → 56.1 ± 1.4; 55 → 65.4 ± 0.5; exponent 0.41). All MEASURED. Summary updated; run_all.sh now includes
  uncertainty.py (fast) and frank2.py + shell3d.py (slow).
- The earlier single-seed shell3d failure was finite-size noise; kept on record above.
- Next recommendation: wave (c), a full re-verification when load < 6; then (d)/(e).

## 2026-09-23 01:05 — wave (c) full re-verification
- Why: waves (a) and (b) are done; no (c) run logged yet. Load 2.2, no train_grasp job.
- Prediction (written before running): all 10 fast + 6 slow scripts pass their self-checks, with numbers matching the README.
  The 3D/random runs may drift slightly if a seed isn't fixed. `--full` is ~30 min, which is over the 20-min per-command limit, so each slow script runs as its own `nice -n 15` command.
- (Added 10:09 by the next wave) This (c) run left no result line: plots were regenerated 01:02–01:04 but no pass/fail
  was logged and nothing was committed. Treat wave (c) as NOT done. Re-run it when load < 6.

## 2026-09-23 10:09 — wave d (Enceladus ocean volume sourced + checked)
- Why: wave (c) is unfinished and 1-min load was 9.0, so no heavy runs. An uncommitted, unlogged README edit (07:46)
  sourced the Enceladus volume to Čadek 2016. Wave (d) verified it rather than leaving it orphaned.
- Prediction: the citation checks out, and the shell volume from the paper's ranges brackets 2.7e16 m³.
- Result (sonnet builder + separate sonnet critic; critic verdict HOLDS):
  - DOCUMENTED: Čadek et al. 2016, GRL 43(11):5653–5660, doi:10.1002/2016GL068634 (Crossref). The abstract gives
    core 180–185 km and ice shell 18–22 km, matching the README.
  - MEASURED: the spherical-shell volume (R = 252.1 km, taken from the literature; ASSUMED precision) over the range
    corners is 2.45–2.93e16 m³, which brackets 2.7e16. ocean.py has a new assert that reads the value from BODIES and
    would fail at 1e16 or 5e16. `ocean.py` exits 0. The patches/uncertainty L = π·232 km is consistent (252 − 20).
  - README wording tightened from "lands close" to the explicit bracket. No reported number changed.
- Failed: nothing. Wikipedia was not reachable, so R = 252.1 km was not re-fetched to the exact digit.
- Next: wave (c) (full re-verification, each slow script separately, when load < 6); then (e).

## 2026-09-23 10:45 — wave (c) full re-verification (interactive session, started at the student's request)
- Load 2.5, no train_grasp job. Each script runs separately under `nice -n 15`, one at a time → verify.log.
- Prediction (written before running): all 10 fast + 6 slow scripts exit 0. Seeded numbers match the README exactly;
  frank*/3D numbers match the README tables to the printed digit (all seeds fixed).
- (12:03, in progress) The fast set passed 10/10 (sim, ocean, domains, patches, scenarios, inherit, dual, scorecard,
  design, uncertainty; all exit 0). The slow set is still running in the interactive session. **Do not start wave (c).**
  Observed: under `nice -n 15`, domains3d.py alone has run 76 min (vs ~10 min un-niced earlier). macOS puts
  low-priority jobs on efficiency cores. So no single 3D script fits the routine's 20-min limit under nice:
  routine waves must skip the 3D scripts (domains3d, coarsen3d, k3, shell3d) and frank*.py.

## 2026-09-23 13:35 — wave d (uncertainty.py headline numbers self-checked)
- Why: wave (c) is still running interactively (coarsen3d under nice, load 4.0), so no heavy runs, and (e) is blocked for the same reason.
  The Summary's headline "100% early Earth / 18% Enceladus and Europa, mixing decides" was printed by uncertainty.py but no assert checked it.
- Prediction: the script reproduces 100/18/18% and a mixing-D split of 0%/36–37% exactly (the RNG seed is fixed).
- Result (sonnet builder + separate sonnet critic; critic verdict HOLDS, and its coverage-gap suggestion was applied):
  - MEASURED: `nice -n 15 uncertainty.py` exits 0 and prints Enceladus 18%, Europa 18%, early Earth 100%; mixing D is the top input for both moons,
    0% (low half) vs 36% (Enceladus) / 37% (Europa) (high half). All match README lines 28–30 and 377–381.
  - New asserts: each world's win rate is within 3 points of the README value, mixing D is ranked first for both moons, and the low/high split matches.
    Perturbation test: changing the README constant to 50% makes the script exit 1; restoring it gives exit 0.
  - No README text or number changed.
- Failed: nothing. Not covered: the "~72% / ~86%" column in the Part 8 table is still printed but not asserted.
- Next: let the interactive wave (c) finish and log its slow-set result. Then (e) needs an un-niced or interactive run,
  because 3D scripts under nice exceed the 20-min limit. Otherwise take another (d) item: the Part 8 ~72%/~86% column.

## 2026-09-23 16:10 — wave c finish + wave d (routine; started with prediction)
- Orient: machine rebooted ~15:37 (uptime 25 min at 16:02). The interactive wave (c) died inside shell3d.py; verify.log
  shows domains3d, coarsen3d, k3, frank, frank2 all exit 0 and shell3d with no result. Load 2.3, no train_grasp.
- Prediction: shell3d.py re-run under nice reproduces shell3d.log exactly (seeded) but may exceed the 19-min kill.
  The five finished slow scripts match their logs digit-for-digit.
- Result (sonnet builder + separate sonnet critic; critic verdict: domains3d check HOLDS, README numbers HOLD except one):
  - MEASURED: shell3d.py under nice exits 0 in 929 s and matches shell3d.log line for line. So wave (c) is now complete:
    fast set 10/10 and slow set 6/6 exit 0 (verify.log, now committed). coarsen3d, k3, frank, frank2 match their logs exactly.
    domains3d prints the same four ratios as before. Its repo log was a stale crash traceback from an older script version; it is now replaced with the verified output.
  - Wave d: domains3d.py had no self-check behind README's "1D patch law fails in 3D (31% spread)". I added law_spread() and
    `assert spread > 0.25`, the mirror image of domains.py's holds-test. Tested on the logged ratios (0.314, passes) and on a
    law-following set (0.014, raises). The full script was not re-run (76 min under nice).
  - CORRECTION: README Part 3 said "A = 5.5 at both mixing strengths"; coarsen3d prints 5.48 and 5.60. The README now has both
    values and a dated correction line. The conclusion is unchanged.
- Failed / caveat: the 31% spread comes from one seed with 2 runs per case, so the 6-point margin over 25% is thin evidence (critic).
- Next: wave (e), a 3D replication with a different seed. k3.py takes 75 s, so it fits the time limit; a k3 seed sweep is the natural choice.

## 2026-09-23 19:05 — wave e (k3.py seed replication; routine; prediction written before running)
- Orient: load 3.8 (1-min), no train_grasp job at start. Waves a–d are done; (e) is the first open wave. k3.py was chosen because it is the only 3D script that fits the 20-min limit under nice.
- Prediction: with 3 new seed pairs, each run still passes its own <30% within-run assert. K3 across seeds stays within ±20% of the logged 3.642. The 1% within-run spread was partly luck; expect 5–15% on new seeds.
- Result (sonnet builder + separate sonnet critic). Prediction partly held: the spread came out even smaller than predicted.
  - MEASURED: k3.py now takes optional seeds (`python k3.py 11 12`). With no arguments it reproduces k3.log digit for digit (critic re-ran it).
    New seeds 11/12, 21/22, 31/32 give K3 = 3.698, 3.481, 3.652, each exit 0 in 159–178 s under nice. Across all four runs:
    mean 3.618, spread (std/mean) 2.3%, range 3.481–3.698. Within-run spreads are 1–2%.
  - Critic HOLDS on reproduction and on the stats. Two real points: (1) the three g values in one run reuse the same noise (domains3d.py:9 re-seeds each call),
    so the within-run 1% is partly built in, and the cross-seed 2.3% is the fair figure; (2) "≈57 cells" came from rounded l.
  - CORRECTION (README Part 3): effective decision volume "≈ 57 cells" → 60.2 cells (MEASURED). The old value multiplied K3 by 2.5³, but the unrounded l = 2.548.
    k3.py now prints K3·l³ from unrounded l. Across seeds, the volume reconstructed from printed fractions is 57–61. The factor-of-1.5 conclusion is unchanged.
  - Critic's claim that 57 carries "±6% rounding uncertainty" was wrong in form: K3·l³ doesn't depend on l at all (K3 ∝ 1/l³). It was simply computed from the rounded l.
- Failed / caveat: no cross-seed assert was added. It would need 4 runs (~11 min under nice) inside one script, too slow for run_all. The evidence is k3_seeds.log.
  NOTE: an unrelated commit 8b5e85b (19:13, generic message, "Fable" co-author) captured this wave's half-finished k3.py/k3_seeds.log/log edits mid-run.
  It was not made by this routine. Left as is (no history rewrite).
- Next: wave (f), the PAPER_DRAFT.md skeleton built only from verified README content.

## 2026-09-23 23:25 — wave f (PAPER_DRAFT.md skeleton; routine; prediction written before drafting)
- Orient: load 2.35, no train_grasp job. Waves a–e are done; (f) is the first open wave. Doc-only wave, no simulations.
- Prediction: a sonnet builder can draft the paper from README alone; a separate sonnet critic plus a mechanical
  number-by-number cross-check will find at least one number or claim that drifted from README (to be fixed before commit).
- Result (sonnet builder + separate sonnet critic; critic verdict NEEDS-FIXES → all 7 fixed). Prediction held: drift was found.
  - DOCUMENTED: PAPER_DRAFT.md (~3,150 words) has an abstract, methods, results (Parts 1–8), limitations, the planned experiment, reproducibility, and references.
    It is labelled DRAFT and not submitted; every number is meant to trace to README.md or EXPERIMENT.md.
  - MEASURED (mechanical check): all 208 numeric strings in the draft occur in README/EXPERIMENT. The only apparent misses were section numbers and DOI-fragment regex artefacts,
    each confirmed by hand. The check cannot see numbers written as words, so every number-word was also checked by hand.
  - Critic found (and I confirmed at the cited README lines) the following. All were fixed in the draft; README and EXPERIMENT were not touched:
    "nine orders of magnitude" for mixing → seven (10⁻¹⁰–10⁻³; README Part 8 says 7); frank2.py mislabelled Part 7 → Part 8;
    Brandenburg & Multamäki called a "founding paper" → closest prior spatial work; the Faglioni "not the natural helix ≠ mirror push" caveat was dropped → restored;
    Europa 79% was missing "at 1 µM"; the 1 M arm was missing "where solubility allows"; the negative controls listed 2 of 3 → all 3.
    I also reworded one paraphrase I found myself ("time limits only if k2 ≲ 1e-6" → README's actual statement at k2 = 1e-6, ~3 nM Enceladus).
  - No README number changed. No scripts touched, so run_all was not needed.
- Failed / caveat: the critic did not check every reference entry exhaustively (spot-checked). Word count is above the 2,500 target.
- Next: waves a–f are all done. The next routine run should log "IDLE: nothing worth doing" unless README changes (e.g. wet-lab results arrive).
- PUSH FAILED (23:50): `git push origin main` → "remote: Repository not found." gh's active account is helloamonkey742-code, but the repo belongs to
  saptarshighosh10-oss (private), so the likely cause is the wrong account (ASSUMED). I did not switch accounts or change credentials. Local main is ahead of origin by 3.
  The student needs to run `gh auth switch` (or log in as the owner) and then `git push origin main`.

## 2026-09-24 08:30 — IDLE: nothing worth doing
- Waves a–f are all done. Wave c ran last at 2026-09-23 16:18, less than 48 h ago. The tree was clean and load was 2.85. No builder or critic was spawned.
- The push retry failed again with "Repository not found" (same gh-account cause as 2026-09-23 23:50, ASSUMED). Local main is 3 commits ahead of origin (4 after this entry). The student needs to run `gh auth switch` and then `git push origin main`.

## 2026-09-24 13:45 — IDLE: nothing worth doing
- Waves a–f are all done. Wave c ran last at 2026-09-23 16:18, less than 48 h ago. The tree was clean and load was 2.52. No builder or critic was spawned.
- The push retry failed again, this time with "Recv failure: Connection reset by peer" (network). `gh auth status` also reports "Failed to log in" for both keyring accounts (helloamonkey742-code active, saptarshighosh10-oss inactive), so the account issue from 2026-09-23 is probably still there (ASSUMED). Local main is 4 commits ahead of origin (5 after this entry). The student needs to re-authenticate as saptarshighosh10-oss (`gh auth login`, then `gh auth switch`) and run `git push origin main`.

## 2026-09-24 18:23 — IDLE: nothing worth doing
- The last two entries were IDLE, so this run ended right away as the rules require. Load was 17.0, so no heavy work would have been allowed anyway. The untracked OUTREACH_EMAILS.md was left alone because it is not the agent's file.

## 2026-09-24 23:22 — IDLE: nothing worth doing
- The last two entries were IDLE, so this run ended right away as the rules require (load 4.05). The untracked OUTREACH_EMAILS.md was left alone.

## 2026-09-25 09:47 — IDLE: nothing worth doing
- The last two entries were IDLE, so this run ended right away as the rules require (load 5.55). The untracked OUTREACH_EMAILS.md was left alone.

## 2026-09-25 13:46 — IDLE: nothing worth doing
- The last two entries were IDLE, so this run ended right away as the rules require (load 1.62). The untracked OUTREACH_EMAILS.md was left alone.

## 2026-09-25 18:05 — wave d (interactive; the student said "continue — is this truly it?")
- **Why:** checked whether "done" was really done. The icy-moon headline (18%, "mixing decides") rests on one input range.
  Sonnet literature agent + direct check of arXiv:2101.10530v2 Sec. II.2 found that Zeng & Jansen's
  "3×10⁻¹⁰ to 3×10⁻³ m²/s" is κ_z, the **vertical** diffusivity (DOCUMENTED). The model used it as the diffusivity for healing
  **across** the ocean (L = πR, horizontal). Same paper: horizontal diffusivity is "much larger than the vertical", and
  horizontal mixing across a hemisphere takes ~1000 yr in their simulation (DOCUMENTED). Zhang, Kang & Marshall 2024
  (Sci. Adv., doi:10.1126/sciadv.adn6857) diagnose lateral eddy diffusivity ~0.03–2 m²/s for Enceladus (DOCUMENTED, via agent).
  Earth's check already separated horizontal and vertical; only the icy moons mixed them up.
- **Fix:** healing needs both horizontal healing over πR with D_h and vertical healing over the ocean depth H with D_z.
  H = V / (4π(L/π)²) from the project's own volumes (≈40 km Enceladus, ≈123 km Europa), so no new inputs. D_h log-uniform
  10⁻²–10² m²/s (ASSUMED bracket around the 0.03–30 literature values); D_z keeps 10⁻¹⁰–10⁻³.
- **Prediction (written before running):** horizontal healing never binds; Enceladus win rate rises from 18% to ~40–50%,
  Europa to ~60–75%; vertical mixing D_z stays the (or a top-two) most decisive input.
- **Result (MEASURED):** Enceladus **43%** (was 18%), Europa **58%** (was 18%), early Earth 100% (unchanged). Sideways
  healing never binds (needs D_h ≥ 5.6×10⁻⁴ / 2.2×10⁻⁴ m²/s, 50–450× below the literature values). Vertical D_z is the most
  decisive input (18%→68% Enceladus, 28%→87% Europa, low vs high half). If the patchwork heals, wins are ~70% / ~88%.
  Prediction check: Enceladus inside the predicted 40–50%. Europa came in **below** the predicted 60–75% (58%), because
  concentration and time still limit it. D_z stayed on top, as predicted.
- **Critic (Sonnet, separate):** HOLDS_WITH_FIXES. All numbers reproduced, win|heals recomputed (70.5%/87.5%). A mutation
  test (old range put back) fails the asserts. Fixed: added the 50–140× conservative margin (0.03 m²/s low end), stated the
  A = 5.5 per-direction extrapolation as untested, and strengthened the layering caveat (300–1400× wider than deep, so
  flat layers stalling may be the normal case → different hands at different depths is possible). Also fixed the stale
  "mixing D" in the Part 8 intro and added a pointer under the Part 2 table. **Rejected:** the claim that Zhang et al. 2024
  never say "of order ~0.1 m²/s" for sideways mixing. The full text (PMC11540039) says "K should be of order ∼0.1 m² s⁻¹",
  where K is their lateral eddy coefficient (checked directly).
- **Checks:** `./run_all.sh` fast 10/10 pass. New asserts in uncertainty.py back every new README number.
  README/PAPER_DRAFT corrections are dated, and the old numbers are kept in the correction notes. EXPERIMENT.md is untouched.
- **Not touched:** OUTREACH_EMAILS.md (the student's untracked file). It still says 18% and should be updated before sending.
- **Next:** the one real open computational question is now whether a layered (strongly anisotropic) ocean coarsens or
  stalls. It can be tested with shell3d.py-style runs using a smaller vertical coupling (e.g. D_z/D_h = 0.01) under nice.
  That test is a worthwhile new wave g; the project was not truly finished.

## 2026-09-25 18:30 — wave g (layered ocean; interactive, the student said "just do that computer step")
- **Idea:** the patch model is a*(1-a²) + D_h∇²_h a + D_z ∂²_z a. The reaction term has no derivatives, so stretching depth
  by √(D_h/D_z) turns it *exactly* into the equal-mixing model in a box of depth H' = H·√(D_h/D_z). So "weak vertical mixing"
  = "equal mixing in a taller box". For Enceladus H/L ≈ 1/18, so H' > L whenever D_z < D_h/333 (Europa: D_h/1386), i.e. for
  most of the plausible vertical range. Test: aniso3d.py, periodic sideways (L = 32), closed top/bottom, depth
  H' = 8 (slab), 32 (cube), 128 (column), 4 seeds, 50% and 55% starts, to t = 1000.
- **Prediction (before running):** slab and cube heal toward one hand at 55% (majority wins). The column collapses into
  stacked horizontal layers in most seeds, and these stop changing between t = 500 and t = 1000 (flat walls don't move),
  so the ocean does NOT reach one hand. If so, Part 3's vertical-healing rule is too optimistic in the column regime.
- **First run (MEASURED, aniso3d.log):** column (H' = 4L) froze into 3–5 stacked layers in 7/8 runs; nothing changed
  t = 500 → 1000. Cube at 55%: 4/4 one hand. Slab: 2/4 one hand; the other 2 are flat side-by-side stripes that wrap around the
  periodic box (exactly 512 side-wall cells = 2 straight walls), a known artifact of small periodic boxes, not layering.
  Prediction for the column confirmed. Prediction for the slab only half right (stripe artifact).
- **Next thought (before running more):** on a real moon a "layer" wall is a sphere (radius ≈ R), not a flat sheet.
  In this model walls move at speed D·curvature (Allen & Cahn 1979), so each wall should creep inward at
  v ≈ 2·√(D_h·D_z)/R (derived through the stretch), and the TOP layer takes over in t ≈ H·R / (2√(D_h·D_z)). Adding (1) a 2D
  shrinking-circle check that this code's walls move at D·curvature and (2) the hand of the top layer.
  Prediction: circle check passes within 10%; with sphericity, win rates drop somewhat from 43%/58% (the deciding volume
  shrinks to the top layer), and if layers never merge (flat bound) they collapse to near 0%.
- **Own error caught before use:** the "√(D_h·D_z)/R" creep speed above is wrong. Stretching depth also steepens the wall's
  curvature; the two effects cancel. The direct result: a nearly flat wall obeys dh/dt = D_h∇²h, so a spherical layer wall
  sinks at 2·D_h/R, independent of D_z. MEASURED (wall_check in aniso3d.py): a wavy wall flattens at 0.00967 vs D_h·k² = 0.00964
  (D_h = 4, D_z = 1) and 0.00239 vs 0.00241 (D_h = 1, D_z = 4). Equal mixing D = 1 came out 22% slow (grid drag on the sharpest
  wall), reported but not asserted. Consequence: layers on a real moon merge in H·R/(2D_h) ≈ 10⁴–10⁵ yr, far shorter than the ages.
- **Result (MEASURED):** column 7/8 frozen layers (3–5 each); cube 55% start 4/4 one hand; slab failures = periodic stripes.
  Circle check +7.1%. Wall check: D_h·k² within 0.3–1% at D_h/D_z = 4 and 1/4. With layering in uncertainty.py (layers merge on
  a round moon; the top layer decides with share L/H'): Enceladus **55%**, Europa **78%**, Earth 100%. If layers never merged
  (flat bound): 3% / 7% / 40%. The most decisive input is now **concentration** (21%→89%, 57%→99%), then bias g (+18 to +23),
  rate k2 (+14), vertical D_z (+10 to +14). Sideways D_h has a slight negative effect (−6 to −8: thinner deciding layer).
  Layers form in 96% / 92% of draws. Merge time ≤ 1.5×10⁴ yr / 2.8×10⁵ yr.
- **Prediction check:** column stall confirmed. My "win rates drop somewhat" prediction was **wrong**: they rose (43→55%, 58→78%),
  because layers merge quickly on a sphere and the old top-to-bottom time rule no longer applies. The "flat bound near 0%"
  prediction was right (3% / 7%).
- **Critic (Sonnet):** HOLDS_WITH_FIXES. All numbers match; the asserts fail under both mutations (wrong creep speed 55→52%,
  no layering 55→69%). Fixed: 7% → 7.1%, a comment on the across∧creep condition, and a great-circle-band argument
  (flagged ASSUMED; no spherical simulation). High-severity note: the student's OUTREACH_EMAILS.md email #5 still quotes 18%
  and "mixing decides"; not edited (the student's file), flagged to the student.
- **Not done:** the full aniso3d.py rerun after adding wall_check was skipped (load ~10, the student was on the machine).
  The 3D code was unchanged; wall_check was run separately and appended to aniso3d.log. aniso3d.py is added to run_all's slow list.
- **Next:** a full 3D spherical-shell layer test would turn the remaining ASSUMED steps (top-layer share, band instability)
  into measurements, but it is heavy. Otherwise the computation really is at a natural stop; the wet lab remains.

## 2026-09-25 18:40 — wave d (routine; light only: load 10–17, then routine-gate check failed on swap 94%)
- **Why:** waves a–g done; the newest claims (Part 3b / Part 8 from wave g) were the least audited. No simulations allowed at this load.
- **Prediction:** most wave-g numbers are already asserted; one or two secondary numbers are not.
- **Result (DOCUMENTED, by reading code and logs; nothing run):** backed already: the aniso3d box table and circle/wall checks
  (aniso3d.log + asserts), 55%/78%/100%, 3%/7%, the concentration splits, D_h/333 and D_h/1386, and the merge times (uncertainty.py asserts).
  **Unbacked:** "layers form in 96%/92%", "bias g +18 to +23", and Earth's "40% if layers never merged". Analytic check (ASSUMED
  log-uniform ranges as in the code): Enceladus layers when log D_h − log D_z > 2.52 → 95.9%; Earth H ≈ 2.6 km, L/H ≈ 7.7e3, layers when
  D_h/D_z > 5.9e7 → ~61% of draws, so the flat bound ≈ 39–40%. All three are consistent with the README; nothing contradicted.
- **Failed:** the builder (opus55-low) wrongly flagged Earth's 40% as contradicted (its stretched-depth estimate was ~3.4e5 m instead of ~2.6e8 m);
  I caught it and reverted its README note. Asserts for the three numbers were written but **reverted uncommitted** because the gate
  blocked running uncertainty.py; untested asserts could break run_all. No critic was run (nothing substantive was committed).
- **Next:** when the machine is free, add asserts to uncertainty.py for README_FLAT_PCT Earth 0.40 (checked for all worlds),
  layered share 0.96/0.92, and the bias-g effect in +18..+23. Then run it and save uncertainty.log. After that the project is at a natural stop.

## 2026-09-25 23:10 — wave d finish (interactive; the student said "do it")
- **Result (MEASURED, uncertainty.log):** uncertainty.py now asserts the three numbers flagged at 18:40. Layered share: Enceladus 96%, Europa 92%.
  Bias-g effect: +23 / +18 points. Early Earth "if layers never merged": 40% (Earth layers in 60% of draws). All match the README. No number changed.
- **Mutation test:** changing each README value (Earth 0.40→0.30, Enceladus 0.96→0.80, g upper 0.23→0.19) makes the script fail. Fast suite: 10/10 pass.
- **Next:** the computational side is at a natural stop. The optional full 3D spherical-shell layer test would measure the two remaining ASSUMED steps
  (top-layer share, band instability), but it is heavy. The wet lab is the remaining real step. OUTREACH_EMAILS.md still says 18%.

## 2026-09-25 23:25 — IDLE: nothing worth doing
- Waves a–g done (g 18:30, d finish 23:10). The only open item is the optional full 3D spherical-shell layer test. It is heavy and not a listed wave, so it was skipped. The load was 5.5 (15-min average 12). The wet lab remains. OUTREACH_EMAILS.md (the student's file) still says 18%.

## 2026-09-26 13:25 — wave c (full re-verification; routine)
- **Why:** last full run was 2026-09-23 16:10 (>48 h); uncertainty.py, aniso3d.py and README Parts 3/3b/8 changed since. Load 3.0, no train_grasp, gate ok (42% free).
- **Prediction (written before running):** fast set 10/10 pass and uncertainty.py matches uncertainty.log; slow scripts coarsen3d, k3, frank, frank2, shell3d, aniso3d exit 0 and match their logs (seeded). domains3d (76 min under nice) exceeds the 20-min limit and is skipped.
- **Result (MEASURED):** fast set 10/10 pass (180 s); uncertainty.py output is identical to uncertainty.log. Slow set under nice: k3 (172 s), frank (277 s), frank2 (216 s), coarsen3d (740 s), aniso3d (316 s) all exit 0 with numbers identical to their logs. aniso3d is the first full run that includes wall_check, and aniso3d.log is replaced with that clean output. shell3d was killed at the 19-min limit (exit 142) after its exponent (0.41) and 50% control (0.498 ± 0.014), both identical; its 52%/55% rows rest on the 2026-09-23 run (929 s then; slower today).
- **Skipped:** domains3d (76 min, over the limit). A first slow pass was skipped because my own fast run pushed the 1-min load to 9. It ran after a back-off.
- **Critic (opus55-low):** HOLDS, with no mismatches in any README number tied to these scripts. Fixed: the README's "within ~1 point" for frank now notes that the no-bias control lands 1.5 points low (48.5% vs 50%). No number changed. Not fixed (minor): Part 3b's "3–5 layers" is read from the wall areas, not printed directly.
- **Prediction check:** right for 9 of 10 scripts tested in full. shell3d did not finish within the limit, so it is only partly verified.
- **Next:** waves a–g done and wave c re-verified. Next routine should log IDLE unless the machine is idle enough for a shell3d rerun (~16–19 min under nice).
- **Push failed:** `git push origin main` returned "remote: Repository not found" for github.com/saptarshighosh10-oss/chirality-project. The commit is local only. It could be an auth or account mismatch, or the repo is missing; not investigated (no credential changes).
- **Push resolved:** the documented gh-account workaround (memory file) worked; origin is in sync.

## 2026-09-27 01:15 — sphere-shell layer test + correction (interactive; the student said "send the subagents out, test everything in parallel")
- **Why:** the only open computational item was the spherical-shell layer test (the ASSUMED top-layer and band steps from wave g). The 22:30 routine pass was blocked by the memory gate (34% free). The student then asked to proceed; the gate passed at 43% on retry. Agents ran in parallel at the student's request, with only one heavy simulation at a time.
- **Prediction (written in sphere_layers.py before running):** wall sink speed 2·D_z/R. (This contradicts the README's 2·D_h/R; the builder derived it before running.)
- **Result (MEASURED, sphere_layers.log):** wall speed 0.991/0.991/0.984 × 2·D_z/r, 0.0% change at D_h ×4. **Wave g's 2·D_h/R was wrong:** on a sphere "down" follows the radius, so a constant-depth wall feels only vertical mixing. Band and hemisphere splits are unstable (rate 1.017× predicted); a centred split stays at 0.50. Top-layer share was 0.52/0.07 against the assumed 0.25 (2 runs; stays ASSUMED).
- **Correction (MEASURED, uncertainty.log):** physics wins Enceladus 55→**22%**, Europa 78→**29%**, early Earth 100→**75%**. The flat bound (3/7/40%) is unchanged. The merge-in-time condition holds in 33%/34% of draws. **Vertical mixing is again the deciding input** (0% vs 43%/57% across its range). Bias g effect is now +8/+4. Worst-case merge is 1.5e12/2.8e13 yr at D_z = 1e-10. Reverting to D_h fails the asserts. Fast suite 10/10. Dated corrections are in README and PAPER_DRAFT (commit 6e7cd5c, pushed).
- **Critics (opus55-low, separate):** physics HOLDS: D_h drops out of the radial operator exactly, and the code has no bugs. Docs HOLD_WITH_FIXES: 2 overclaims in Part 3b ("not a tilt" was an argument, not a run; the creep match was post-hoc from 2 runs) and loose "merge in time" wording; all fixed before commit.
- **OUTREACH_EMAILS.md** (the student's untracked draft; edited at their request, nothing sent): numbers updated to 22/29/75%; A ≈ 5.5 → 5.48/5.60; email 5 now asks about real vertical-mixing strength instead of whether layers merge.
- **shell3d re-verification:** the rerun reproduced the exponent 0.41, the f=0.5 control (0.498 ± 0.014) and the **f=0.52 row (0.561 ± 0.014) identically**. The f=0.55 row was not re-verified: the full run hit the 25-min cap, the separate driver died with its agent (API error), and a retry was blocked by the gate (30% free, load 12). The 55% row still rests on the 2026-09-23 run.
- **Failed/limits:** sphere_layers varied D_z only together with grid spacing (critic note; a fixed-dz D_z sweep would remove that confound). The sphere grids are coarse (R = 6–16 cells), with 2–3 seeds per case.
- **Next:** (1) re-verify shell3d's f=0.55 row alone when the machine is free; (2) optional sphere_layers run with D_z varied at fixed dz; (3) the key open science question is now the real vertical-mixing strength in icy-moon oceans. That is a literature wave (b-style), and it decides the answer.

## 2026-09-27 08:25 — wave c (partial: shell3d f=0.55 row only; routine)
- **Why:** the last log's first "Next" item; it is the only README number from shell3d not re-run since 2026-09-23. Load 1.5, gate ok (64% free), no train_grasp.
- **Prediction (written before running):** seeded, so identical to shell3d.log: f=0.55 goes 0.551 -> 0.654 +- 0.005 over 4 seeds at 512x512x8, t=1 -> 200.
- **Result (MEASURED):** 0.551 -> 0.654 +- 0.005 (change +0.103), per-seed ends 0.658/0.645/0.659/0.654; identical to shell3d.log. Ran shell3d.run() for f=0.55 only (270 s under nice). Every shell3d row (exponent, 50%, 52%, 55%) has now been re-run since the 2026-09-23 original.
- **Critic:** skipped by choice. This is a seeded reproduction compared string-for-string with a committed log, so there is nothing to interpret. No README number changed.
- **Prediction check:** right.
- **Next:** optional sphere_layers D_z sweep at fixed dz (removes the grid-spacing confound); the open science question is the real vertical-mixing strength in icy-moon oceans (literature). Otherwise IDLE.

## 2026-09-27 13:28 — wave b (literature: real icy-moon vertical mixing; routine, light only)
- **Why:** the last two logs name real vertical-mixing strength as the open question that decides the answer. Load 9.4 at start → no simulations. Gate ok (53% free).
- **Prediction (written before results):** published κ_z spans most of the sampled range with no consensus; several model values sit ≥ ~1e-5; the answer stays conditional.
- **Result (DOCUMENTED, vertical_mixing_lit.md):** Zeng & Jansen 2021 Sec. 2.2 (arXiv:2101.10530v2, re-read directly) is still the only primary estimate: κ_z ≈ 3e-10–3e-3, matching the project range. Ames et al. 2025 quote "1e-7–1e-3", but their ref. 31 is the same paper, so it is not a second estimate. Model runs use 5e-5 (Z&J), 5e-3 (Kang et al. 2022, doi:10.1126/sciadv.abm4665) and 1e-5–1e-3 (Ames). These are numerical choices, and all of them clear the merge bar. No Europa-specific value was found; the Enceladus range is used for Europa (ASSUMED).
- **Merge bar (MEASURED arithmetic, H·R/(2t)):** Enceladus D_z ≳ 1.5e-7 (1 Gyr) to 1.5e-4 (1 Myr); Europa ≳ 7e-7 (4 Gyr) to 3e-5 (100 Myr). The published range straddles it, so the literature does not settle the answer.
- **New caveat (DOCUMENTED, same Sec. 2.2):** the molecular tracer floor is ~1e-9, so the bottom 14% of the sampled D_z range is below what molecules do by themselves. No percentage changes, since those draws fail every bar in uncertainty.py either way (checked by reasoning, not rerun). The quoted worst-case merge times (1.5e12/2.8e13 yr at 1e-10) would be ~10× shorter at the floor, which is still far beyond the Solar System's age. This was added as a dated README note; no reported number was edited.
- **Failed/caught:** the builder (opus55-low) compared against the retired heal-across-depth thresholds (2e-9..2e-6), which are 70–180× too low, and wrongly counted Ames as independent. The critic (opus55-low, separate; verdict HOLDS_WITH_FIXES) caught both, and both were fixed before commit. The IOP PDF was bot-walled; the arXiv copy worked.
- **Prediction check:** mostly right (no consensus; model values ≥1e-5). Wrong in one detail: there is only ONE primary estimate, not several.
- **Next:** optional: floor D_z at 1e-9 in uncertainty.py and rerun (the percentages should stay identical; the worst-case merge line changes). Optional sphere_layers D_z sweep at fixed dz. Otherwise IDLE; the wet lab remains.

## 2026-09-27 18:30 — wave d (quality: sphere wall speed vs D_z at fixed grid; routine)
- **Why:** waves a–g are done. The last logs name the open confound: sphere_layers varied D_z only together with grid spacing dz (D_z = 4·dz² hard-wired). The 2·D_z/R result that set the 22/29/75% numbers rests on it. Load 2.7, gate ok (42% free), no train_grasp.
- **Prediction (written in code before running):** at fixed dz = 0.02, speed/(2·D_z/r) in 0.8–1.25 for D_z = 0.0008/0.0016/0.0032, and speed doubles within 15% per doubling.
- **Result 1c (MEASURED, pre-registered, sphere_layers_1c.log):** 0.967 / 0.994 (control; matches section 1's 0.991) / **1.393 → FAILED**, exit 1. The D_z=0.0032 wall moved 0.098 vs 0.065 and ended ~1.6 wall-widths above the no-flux floor.
- **Result 1d (MEASURED, POST-HOC, 1 seed, sphere_layers_1d.log):** twice the depth (nz=48) → 0.967 / 0.986 / 0.995, doubling ×2.04 / ×2.02, walls 10/7/5 widths above the floor, exit 0, 154 s. This is consistent with floor pull as the cause of 1c; it is not proof (depth and start radius changed too).
- **Conclusion:** on a fixed grid the speed follows D_z, so the confound does not undo 2·D_z/R. No reported number changed. README Part 3b has a dated bullet; section 1 code path is unchanged (Dz=None default; not re-run, checked by reading + compile).
- **Critic (opus55-low, separate): HOLDS_WITH_FIXES.** Fixed: "~2.5 widths" → ~1.6; "shown" → "seems"; scope limited to 1 seed/R=16; ratios stated as using mid-depth r; D_z=0.0008 case is under-resolved (~1.4 cells), so its 3% shortfall may be a grid effect.
- **Prediction check:** wrong for the pre-registered 1c (floor artifact); right for the post-hoc 1d.
- **Next:** optional: repeat 1d with a second seed (pre-registered, ~3 min) and a fixed-D_z clearance sweep to pin down the floor-pull distance. Optional D_z 1e-9 floor in uncertainty.py. Otherwise IDLE; the wet lab remains.

## 2026-09-27 23:30 — wave e (replication: sphere 1d with seeds 1 and 2; routine)
- **Why:** waves a–g are done. The last log's first "Next" was to repeat the post-hoc 1d (sphere wall speed follows D_z on a fixed grid) with another seed. Load 3.5, gate ok (69% free), no train_grasp.
- **Prediction (written in code before running):** seeds 1 and 2 give every ratio in 0.95–1.00 (seed 0: 0.967/0.986/0.995) and each doubling in 1.9–2.1; the 1d asserts pass.
- **Result (MEASURED, sphere_layers_1e.log, 179 s under nice, exit 0):** both seeds give 0.967/0.986/0.995 and doublings ×2.038/×2.018, identical to 3 decimals (spread 0.000). The seed did change the run (start radii differ in the 4th decimal).
- **Interpretation:** the seed only sets the bump pattern (std 1.5·dz), which D_h = 1 flattens quickly, so the match was expected. It rules out seed luck for this setup only; it is not an independent test. The open limits stay grid and depth (1c). No reported number changed. README Part 3b has a dated note. Code: section1d takes a seed (default 0, so the 1d path is unchanged apart from the header text), and a new section "1e" was added.
- **Critic (opus55-low, separate): HOLDS_WITH_FIXES.** Fixed: "±1.5 grid steps" → typical size (it is a std, not a bound); "never the weak spot" → limited to R=16, D_h=1, one grid.
- **Prediction check:** right.
- **Next:** remaining optional items are low value (the D_z 1e-9 floor in uncertainty.py leaves every % unchanged). Next routine should log IDLE; the wet lab remains.

## 2026-09-28 09:37 — IDLE: nothing worth doing
- Waves a–g done; last log recommended IDLE; wave c re-verify not due until 2026-09-28 13:25 (48 h rule). Gate ok (62% free), load 4.1.

## 2026-09-28 18:30 — wave c (partial full re-verification; routine)
- **Why:** last full run 2026-09-26 13:25 (>48 h); uncertainty.py and sphere_layers.py changed since. Load 1.9 at start, gate ok (82% free), no train_grasp.
- **Prediction (written 18:24 before the slow runs):** sphere_layers 1/2/2b/3, k3, frank, frank2, aniso3d, coarsen3d all exit 0 with numbers identical to their committed logs (seeded, code unchanged). shell3d and domains3d skipped (over the 20-min limit; every shell3d row was re-run 09-26/27).
- **Result (MEASURED):** fast set 10/10 pass (86 s); uncertainty.py output identical to uncertainty.log (the 22/29/75% numbers). sphere_layers sections 1 and 3 (49 s; the 2·D_z/R wall-speed test behind the 2026-09-27 correction), k3 (89 s) and frank (177 s) all exit 0, and every output line matches its log.
- **Skipped:** frank2, aniso3d, sphere_layers 2 and 2b, coarsen3d. The 1-min load reached 6.35 after frank, mostly the student's Chrome (GPU helper ~260% CPU), so heavy work stopped. No commit of any half-done output.
- **Critic:** skipped. These are seeded reproductions compared line by line with committed logs; no README number changed.
- **Prediction check:** right for every script that ran.
- **Next:** finish the skipped five (≈ 5 + 5 + 14 + 8 + 12 min under nice) when the machine is idle; otherwise IDLE. The wet lab remains.

## 2026-09-29 08:29 — wave c (finish skipped re-verify; routine) — ABORTED, nothing verified
- **Why:** the 09-28 partial re-verify left frank2, aniso3d, sphere_layers 2/2b and coarsen3d un-run. Gate ok (76% free). Load 8.4 at start, 5.7 twenty seconds later, so one slow script was tried. No train_grasp.
- **Prediction:** frank2 exits 0 in about 5 min, and its output matches frank2.log line by line (F = 0.71).
- **Result:** frank2.py ran from 08:30 to about 11:10. That is 9486 s of wall time under nice -n 15, and the load hit 12–14 (the student was using the machine; the Mac may also have slept). It was killed (exit 143) with no output. MEASURED: nothing. No number changed, and no repo file changed except this log.
- **Mistake:** the 20-min limit was only an intention. The command had no hard timeout, so it overran by about 8×. Next time wrap each slow run in `timeout 1200` (or check the clock), and re-check load between scripts.
- **Next:** finish frank2, aniso3d, sphere_layers 2/2b and coarsen3d one at a time with a hard 20-min kill, only when load is < 3. Otherwise IDLE. The wet lab remains.

## 2026-09-29 14:05 — IDLE: nothing worth doing
- Only open item is finishing the wave c re-verify (frank2, aniso3d, sphere_layers 2/2b, coarsen3d), which needs load < 3. Load 6.2 now (student on the machine), so no heavy runs. Gate ok (37% free). No number changed.

## 2026-09-30 08:45 — wave c (finish skipped re-verify; routine) — frank2 only
- **Why:** frank2, aniso3d, sphere_layers 2/2b and coarsen3d were still un-re-run since 09-28. Gate ok (83% free), load 2.5 at start, no train_grasp.
- **Prediction:** frank2 exits 0 and its table matches frank2.log line by line (F = 0.71).
- **Result (MEASURED):** frank2.py exit 0; all output lines are identical to frank2.log (F = 0.71, rows 0.493/0.785/0.952/0.948). No number changed.
- **Mistake (again):** wall time was 12 097 s (3.4 h), not ~5 min. The hard kill was `perl -e 'alarm 1200; exec @ARGV' nice -n 15 python ...`; it did not fire (the alarm is lost across the double exec through `nice`, or real time paused during sleep — not diagnosed). Load rose to 8–10 during the run (the student on Chrome). Next time run python in the background and kill it from a watcher (`( sleep 1200; kill $PID ) &`), which does not depend on alarm inheritance.
- **Skipped:** aniso3d, sphere_layers 2/2b, coarsen3d (load 6–10 at 13:45).
- **Critic:** skipped — seeded reproduction compared line by line with a committed log.
- **Next:** aniso3d, sphere_layers 2/2b, coarsen3d one at a time with a watcher kill, only at load < 3. Otherwise IDLE. The wet lab remains.

## 2026-09-30 23:42 — wave c (finish skipped re-verify; routine) — DONE
- **Why:** aniso3d, sphere_layers 2/2b and coarsen3d were the last scripts not re-run since 09-28. Gate ok (43% free), load 2.5–4.0 during the pass, no train_grasp.
- **Prediction:** each exits 0 and its output matches its committed log line by line (seeded, code unchanged).
- **Result (MEASURED):** aniso3d (155 s), sphere_layers 2 (330 s), sphere_layers 2b (218 s), coarsen3d (390 s, A = 5.48/5.60) all exit 0; every output line matches the log except timing stamps and the log's hand-written trailer. No number changed.
- **Fix that worked:** each run went in the background with a separate `( sleep 1200; pkill ... ) &` watcher; none needed killing.
- **Critic:** skipped (seeded reproductions compared line by line with committed logs).
- **Note:** student's uncommitted EXPERIMENT.md and pdf changes left untouched.
- **Next:** the wave c re-verify cycle is complete (frank2 09-30, rest now). Newer reviewer-revision scripts (review.py, shell_res.py) were not in this cycle; the next routine could re-verify those, otherwise IDLE. The wet lab remains.

## 2026-10-01 04:45 — wave c (re-verify review.py; routine)
- **Why:** review.py (reviewer revision R1–R8, Bennu 2/4/15%, scale rule, flips, beta bias) changed 2026-09-30 and was never re-run by the routine. Gate ok (51% free), load 3.4–3.6, no train_grasp.
- **Prediction (written before the run):** exit 0 and output identical to review.log line by line (deterministic, code unchanged since 28ca8e9).
- **Result (MEASURED):** exit 0 in 1070 s under nice -n 15 (background run + `( sleep 1200; kill ) &` watcher, not needed). All 90 output lines identical to review.log. No number changed.
- **Note:** 1070 s is close to the 20-min limit; on a busier machine review.py may need splitting by section. shell_res.py (>1 h) still un-re-run — over the limit, left for a by-hand run. Student's uncommitted EXPERIMENT.md/pdf changes left untouched.
- **Critic:** skipped (deterministic reproduction compared line by line with a committed log).
- **Next:** everything re-verifiable within limits is now re-verified. IDLE unless the README/paper changes again. The wet lab and the student's paper revision remain.

## 2026-10-01 08:35 — wave b (literature: compartment citation; routine, light only)
- **Why:** the 2026-09-30 reviewer revision cited "Milner-White & Russell 2005" for connected compartments, marked "not yet read in full", with no title/volume/DOI. It was the only unverified reference left from that revision. Gate ok (55% free), load 3.0, no train_grasp. No simulations.
- **Prediction:** the citation exists in OLEB 2005 but may not be the canonical source for the compartment-mound idea.
- **Result (DOCUMENTED, Crossref + OpenAlex):** Milner-White & Russell 2005 = "Sites for phosphates and iron-sulfur thiolates in the first membranes: 3 to 6 residue anion-binding motifs (nests)", OLEB 35(1) 19–27, doi:10.1007/s11084-005-4582-7. Its subject is peptide motifs in early membranes. Abstract not retrievable (OpenAlex has none; Springer/PubMed/Europe PMC failed or timed out). Added Martin & Russell 2003, Phil. Trans. R. Soc. Lond. B 358(1429) 59–85, doi:10.1098/rstb.2002.1183, the standard hydrothermal-mound compartment source, alongside it in PAPER_DRAFT §3.11 and the reference list. No number changed; README untouched (it does not cite either paper).
- **Critic (opus55-low, separate): HOLDS_WITH_FIXES.** Fixed: journal abbreviation "Lond.", issue 1429, and dropped the unread negative claim "not a compartment-network model".
- **Open:** "(Russell, *Scientific American*)" in the same sentence has no year and no reference entry — the student should complete or drop it. pdf/ not rebuilt.
- **Prediction check:** right.
- **Next:** IDLE unless the paper changes again. The wet lab and the student's paper revision remain.

## 2026-10-01 13:50 — wave d (quality: arm E claims added today; routine, light only)
- **Why:** the paper/EXPERIMENT changed after the last routine wave (commits 7678840, 41a2487 added arm E) with an unchecked number: "transamination has a measured rate near 10⁻⁴ /M/s (our estimate from Yu et al. 2024, Table 2), well above the 1-year detection floor", and "E is the only network with a published prediction". Gate ok (54% free), load 1.4, but a train_grasp-related job was running → no simulations.
- **Prediction:** DOIs fine; the rate is a back-calculated estimate, not "measured"; "only" unsupported.
- **Result (DOCUMENTED, Crossref + Europe PMC PMC10873602 full text):** all three DOIs (Yu 2024 PNAS, Deng 2024 Nature, Higgs & Blackmond 2025 PNAS) CONFIRMED. Higgs & Blackmond abstract says the network "may exhibit symmetry breaking" (model result). Yu et al. Table 2 gives % conversion, not rate constants, and is the **reverse** reaction (alanine + pyridoxal → pyruvate + pyridoxamine). Second-order back-calculation (ASSUMED model k = x/(0.1·(1−x)·t)): 2.2e-5 uncatalysed, 1.4e-4 Pro-Pro, down to ~6e-6 (Val-Ile, 9%/48 h). The "well above the detection floor" comparison was invalid: the floor is for the amplification rate k2, not conversion.
- **Fix:** EXPERIMENT.md arm table + priority paragraph and PAPER_DRAFT §limitations + arm list reworded with dated corrections: "a" (not "the only") network with a model prediction; range 6e-6–1e-4 /M/s from the reverse reaction, our estimate; no floor claim. Hypothesis and Decision rule untouched. No model number changed; README untouched (it does not contain this claim).
- **Critic (opus55-low, separate): HOLDS_WITH_FIXES** — caught the reverse-reaction and k2-vs-conversion errors in the first fix; both applied.
- **Prediction check:** right on "estimate" and "only"; missed the reverse-direction and wrong-benchmark problems (critic found them).
- **Open for the student:** pdf/ not rebuilt; the forward-rate (Yu Table 1 / SI Fig. 2) was not read.
- **Next:** IDLE unless the paper changes again.
