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
