#!/usr/bin/env bash
# Re-run every model and check. Default: fast set (~2 min). --full adds the 3D runs (~50 min).
# --paper adds the scripts behind the paper's headline tables and figures (~70 min); --full --paper runs both.
# A failed assertion stops the run with a non-zero exit (verify_seeds.py instead prints ASSERT FAILED per seed). A failed
# prediction (Section 2.7 of the paper) does not stop a script: read the HELD/FAILED lines and each RESULTS docstring.
set -euo pipefail
cd "$(dirname "$0")"
if [ -z "${PYTHON:-}" ]; then
  if [ -x "/opt/anaconda3/bin/python" ]; then
    PY="/opt/anaconda3/bin/python"
  else
    PY="python3"
  fi
else
  PY="$PYTHON"
fi

fast=(sim.py ocean.py domains.py patches.py scenarios.py inherit.py dual.py dna_network.py stress_test.py chiral_search.py scorecard.py design.py uncertainty.py sketch.py quench.py strat_sweep.py jupiter_moons.py ocean_worlds.py robust_ranking.py crystal_check.py soai_check.py viedma_check.py blind_check.py blind_units.py breakeven.py earth_breakeven.py assumed_rules.py nudge_realism.py prebiotic_search.py test_candidate_chemistry.py extensive_testing.py lethal_stress_tests.py planetary_chirality_audit.py titan_95k.py titan_rates.py impact_melt.py resolve_challenges.py)
# sphere_layers.py alone takes ~22 min (layer-wall, patchwork and band tests on a spherical shell).
slow=(domains3d.py coarsen3d.py k3.py frank.py objections.py moon_nudges.py frank2.py shell3d.py aniso3d.py sphere_layers.py review.py coarsen_fig.py)
# shell_res.py (grid-resolution check, 512x512x16) takes over an hour; run it by hand.
# Paper tier (fresh-clone times, 2 threads): each output line matches the committed <script>.log apart from timings.
# kz_large.py alone takes ~14 min. rerun_base.py rewrites rerun_base.log and journals/new_analyses_rerun.md (identical).
paper=(sweep2.py sweep3.py "multiwall.py 1" "multiwall.py 2" "multiwall.py 3" "multiwall.py 4" seedfix.py seedfix2.py
  rerun_base.py tidal_mix.py tidal_mix_g005.py handoff.py "deephand_scale.py 3" rho2side.py "rho2side.py posthoc"
  "front2d.py worlds" pin_onset.py two_zone.py polar_gap.py currents.py realistic.py earth.py field_walls.py
  field_flips.py wobble.py wobble_sweep.py verify_seeds.py kz_large.py journals/make_fig_schematic2.py
  journals/make_fig_sensitivity3.py make_fig_shares3.py make_fig_threshold2.py
  # Added 2026-10-16 (~17 min; advect2d.py alone ~13 min). life_test.py and life_few.py need scikit-learn,
  # abiotic_vs_life.py needs rdkit (requirements.txt).
  abiotic_vs_life.py advect2d.py coarsen_slow.py g_wide.py kcheck.py lid_melt.py life_few.py life_test.py persist.py
  plume_feedstock.py plume_share.py plume_signal.py sugar_bands.py sugar_measured.py sugar_pair.py
  blind_lab.py)   # ~1 min
# Not in any tier (hours each, or extra packages); run by hand with the section arguments in each docstring:
# front2d.py W1..3Db, gradient_ramp.py, kz_patches.py, deephand_scale.py 1 <w> (section 3 reads its ROW lines),
# realistic.py 5, sphere_layers.py 2b, verify_sstar.py (needs gsw), frank_sbml.py (needs antimony + libroadrunner).

scripts=("${fast[@]}")
for a in "$@"; do
  case "$a" in
    --full) scripts+=("${slow[@]}") ;;
    --paper) scripts+=("${paper[@]}") ;;
    *) echo "unknown option $a (use --full and/or --paper)"; exit 2 ;;
  esac
done

for s in "${scripts[@]}"; do
  echo "== $s"
  "$PY" $s 2>&1 | { grep -v Warning || true; }  # $s unquoted: some entries carry a section argument
done
echo "== finished ${#scripts[@]} scripts with no error exit (failed predictions do not stop a script; see each RESULTS docstring)"
