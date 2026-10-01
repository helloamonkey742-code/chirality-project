#!/usr/bin/env bash
# Re-run every model and check. Default: fast set (~3 min). --full adds the 3D runs (~30 min).
# Any failed self-check stops the script with a non-zero exit.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

fast=(sim.py ocean.py domains.py patches.py scenarios.py inherit.py dual.py scorecard.py design.py uncertainty.py sketch.py)
# sphere_layers.py alone takes ~22 min (layer-wall, patchwork and band tests on a spherical shell).
slow=(domains3d.py coarsen3d.py k3.py frank.py frank2.py shell3d.py aniso3d.py sphere_layers.py review.py coarsen_fig.py)
# shell_res.py (grid-resolution check, 512x512x16) takes over an hour; run it by hand.

scripts=("${fast[@]}")
[[ "${1:-}" == "--full" ]] && scripts+=("${slow[@]}")

for s in "${scripts[@]}"; do
  echo "== $s"
  "$PY" "$s" 2>&1 | grep -v Warning
done
echo "== all checks passed (${#scripts[@]} scripts)"
