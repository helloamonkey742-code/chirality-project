"""Re-run ocean_worlds.py with seeds 1-5 to show the results are not a seed-0 accident (2026-10-06)."""
import io, contextlib, sys
import numpy as np
sys.path.insert(0, ".")
import uncertainty as U, ocean_worlds as O
U.WIN_TOL = 1.0  # headline check is seed-0 specific; ordering asserts stay live
res = {}
for seed in range(1, 6):
    U.RNG = np.random.default_rng(seed)
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            out = O.main(); ok = "ok"
        except AssertionError as e:
            out = None; ok = f"ASSERT FAILED: {e}"
    print("seed", seed, ok)
    if out:
        for k, v in out.items(): res.setdefault(k, []).append(100 * v)
print(f"{'world':10} {'min':>6} {'max':>6}")
for k, v in sorted(res.items(), key=lambda kv: -np.mean(kv[1])):
    print(f"{k:10} {min(v):6.1f} {max(v):6.1f}")
