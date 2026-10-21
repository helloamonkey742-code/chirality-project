"""h13_sphere_shell: the thin 3D version of h13_sphere.py - a spherical shell of equal mixing, using sphere_layers.Shell.

Same rule (da/dt = a - a^3 + D lap(a), no bias, patchwork start at volume share 0.5), shell R = 16, nz = 8 layers of
dz = 0.5, D_z = 4 dz^2 = 1 = D_h, so the stretched depth H' = 4 is far below L = pi R = 50 (case (i) of the draft).
T1 = first time the minority volume share is < 1% (checked every 2 time units); tau_sw = (pi R / 5.5)^2 / D.
Pre-registered as P8 in h13_sphere_predictions.txt (hash checked through h13_sphere.need_pred).
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python h13_sphere_shell.py n=12 seed=0 procs=4 > h13_sphere_shell.log
Smoke runs before the results run: none.
"""
import sys
import time
from multiprocessing import Pool

import numpy as np

import h13_sphere as H
import sphere_layers as SL

R, NZ, DZ, F0 = 16, 8, 0.5, 0.5


def run_one(seed):
    sh = SL.Shell(R, NZ, DZ, 1.0)
    SL.check_stability(sh)
    a, t = SL.patchwork(sh, F0, np.random.default_rng(seed)), 0.0
    assert abs(sh.fav(a) - F0) < 0.02, "patchwork volume share is not the requested one"
    tmax = 30 * H.tau("sphere", R)
    while t < tmax:
        a, t = sh.step_to(a, t, t + 2.0)
        f = sh.fav(a)
        if min(f, 1 - f) < 0.01:
            return seed, t, 1 if f > 0.5 else -1
    return seed, np.nan, 0


if __name__ == "__main__":
    H.need_pred()
    kw = dict(a.split("=") for a in sys.argv[1:])
    n, seed, procs = int(kw.get("n", 12)), int(kw.get("seed", 0)), int(kw.get("procs", 4))
    ta = H.tau("sphere", R)
    print(f"# shell R={R} nz={NZ} dz={DZ} D_z=D_h=1 f0={F0} n={n} seed0={seed} tau_sw={ta:.2f} tmax={30 * ta:.0f}")
    t0, res = time.time(), []
    with Pool(procs) as p:
        for s, T1, w in p.imap_unordered(run_one, range(seed, seed + n)):
            res.append(T1)
            print(f"run seed={s} T1={T1:.1f} T1/tau_sw={T1 / ta:.3f} winner={w:+d}", flush=True)
    T = np.array(res) / ta
    print(f"SUMMARY shell: reached one hand {np.isfinite(T).sum()}/{n} by 30 tau_sw; T1/tau_sw median {np.median(T):.2f} "
          f"min {np.nanmin(T):.2f} max {np.nanmax(T):.2f}")
    print(f"# wall time {time.time() - t0:.0f} s")
