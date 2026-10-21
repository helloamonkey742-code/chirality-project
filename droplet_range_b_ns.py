"""POST HOC (after droplet_range_b.log and droplet_range_b_pool.log were seen; no new simulation): all five verdict shares of rule B pooled over
the S1 = (0.35, 1.39) grid points of droplet_range_b.py's u-scan, per design / s / misread, to see what takes the STEP share below 0.95
at true misread 0.06.  Reads the cache of droplet_range_b.py.  Not a pre-registered test.

Run: DR_CACHE=<scratch dir> UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy python droplet_range_b_ns.py > droplet_range_b_ns.log
"""
import os
import pickle
from pathlib import Path

import numpy as np

U = (0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.35, 0.4, 0.57, 0.8, 1.0, 1.13, 1.39, 1.6, 1.9, 2.26, 2.7, 3.2, 4.5, 6.4)
NAMES = ("STEP", "SMOOTH", "FLAT", "RISING", "NS")


def shares(outs, lo=0.35, hi=1.39):
    c = sum(o["counts"]["B"] for u, o in zip(U, outs) if lo <= u <= hi)
    return c / c.sum(), int(c.sum())


assert abs(shares([dict(counts={"B": np.array([3, 1, 0, 0, 1])})] * 20)[0][0] - 0.6) < 1e-12
R = pickle.loads((Path(os.environ["DR_CACHE"]) / "b_all.pkl").read_bytes())
for d in ("base", "ext4", "K10", "K12"):
    for s in (0.1, 0.15):
        for e in (0.0, 0.03, 0.06):
            p, n = shares(R[("scan", d, s, e)])
            print(f"{d:4s} s={s} misread {e:.2f}: " + " ".join(f"{k} {v:.3f}" for k, v in zip(NAMES, p)) + f"  (n {n})")
