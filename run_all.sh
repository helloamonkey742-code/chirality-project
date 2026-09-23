#!/usr/bin/env bash
# Re-run every model and check. Default: fast set (~3 min). --full adds the 3D runs (~30 min).
# Any failed self-check stops the script with a non-zero exit.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

fast=(sim.py ocean.py domains.py patches.py scenarios.py inherit.py dual.py scorecard.py design.py uncertainty.py)
slow=(domains3d.py coarsen3d.py k3.py frank.py frank2.py shell3d.py)

scripts=("${fast[@]}")
[[ "${1:-}" == "--full" ]] && scripts+=("${slow[@]}")

for s in "${scripts[@]}"; do
  echo "== $s"
  "$PY" "$s" 2>&1 | grep -v Warning
done
echo "== all checks passed (${#scripts[@]} scripts)"
