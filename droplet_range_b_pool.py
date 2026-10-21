"""POST HOC pooling of droplet_range_b.py's u-scan (read from its cache; no new simulation).  Written after droplet_range_b.log was seen,
because the pre-registered rule (b) takes the minimum of 7-14 noisy shares (300 runs each) and so is biased low.  Not a pre-registered test.
Prints, per design / s / misread: STEP share pooled over the grid points of S1 = (0.35, 1.39) and S2 = (0.14, 2.8), with the pooled
Wilson 95% lower bound, and the pooled share over the S2 points outside S1 (u <= 0.28 and u >= 1.6).

Run: DR_CACHE=<scratch dir> UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy python droplet_range_b_pool.py > droplet_range_b_pool.log
"""
import os
import pickle
from pathlib import Path

import numpy as np

U = (0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.35, 0.4, 0.57, 0.8, 1.0, 1.13, 1.39, 1.6, 1.9, 2.26, 2.7, 3.2, 4.5, 6.4)
S1, S2 = (0.35, 1.39), (0.14, 2.8)
DS, SS, EPS = ("base", "ext4", "K10", "K12"), (0.1, 0.15), (0.0, 0.03, 0.06)


def wilson_lo(k, n, z=1.96):
    p = k / n
    return (p + z * z / (2 * n) - z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


def pooled(outs, sel):
    k = sum(o["counts"]["B"][0] for u, o in zip(U, outs) if sel(u))
    n = sum(o["counts"]["B"].sum() for u, o in zip(U, outs) if sel(u))
    return k / n, wilson_lo(k, n), n


assert abs(wilson_lo(100, 100) - 100 / (100 + 1.96 ** 2)) < 1e-12 and abs(wilson_lo(50, 100) - 0.4038) < 1e-3
R = pickle.loads((Path(os.environ["DR_CACHE"]) / "b_all.pkl").read_bytes())
for d in DS:
    for s in SS:
        for e in EPS:
            o = R[("scan", d, s, e)]
            in1, in2 = (lambda u: S1[0] <= u <= S1[1]), (lambda u: S2[0] <= u <= S2[1])
            out = lambda u: in2(u) and not in1(u)
            a, b, c = pooled(o, in1), pooled(o, in2), pooled(o, out)
            print(f"{d:4s} s={s} misread {e:.2f}: S1 {a[0]:.3f} (lo {a[1]:.3f}, n {a[2]})  S2 {b[0]:.3f} (lo {b[1]:.3f}, n {b[2]})  S2 outside S1 {c[0]:.3f} (n {c[2]})")
