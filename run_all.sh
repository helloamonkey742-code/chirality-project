#!/usr/bin/env bash
# Re-run every model and check. Default: fast set (~4-7 min). --full adds the 3D runs (~30 min).
# Any failed self-check stops the script with a non-zero exit.
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

fast=(sim.py ocean.py domains.py patches.py scenarios.py inherit.py dual.py dna_network.py stress_test.py chiral_search.py scorecard.py design.py uncertainty.py sketch.py quench.py strat_sweep.py jupiter_moons.py ocean_worlds.py robust_ranking.py crystal_check.py soai_check.py viedma_check.py breakeven.py earth_breakeven.py assumed_rules.py nudge_realism.py prebiotic_search.py test_candidate_chemistry.py extensive_testing.py lethal_stress_tests.py planetary_chirality_audit.py titan_95k.py titan_rates.py resolve_challenges.py)
# sphere_layers.py alone takes ~22 min (layer-wall, patchwork and band tests on a spherical shell).
slow=(domains3d.py coarsen3d.py k3.py frank.py objections.py moon_nudges.py frank2.py shell3d.py aniso3d.py sphere_layers.py review.py coarsen_fig.py)
# shell_res.py (grid-resolution check, 512x512x16) takes over an hour; run it by hand.

scripts=("${fast[@]}")
[[ "${1:-}" == "--full" ]] && scripts+=("${slow[@]}")

for s in "${scripts[@]}"; do
  echo "== $s"
  "$PY" "$s" 2>&1 | grep -v Warning
done
echo "== all checks passed (${#scripts[@]} scripts)"
