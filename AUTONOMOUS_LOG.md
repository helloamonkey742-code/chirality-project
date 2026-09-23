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
