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
