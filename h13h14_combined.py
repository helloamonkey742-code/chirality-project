"""h13h14_combined: world shares under the COMBINED sphere x eddy factor f on the sweep's sideways term. (H13 x H14)

Thin wrapper. It runs h13_sphere_shares.py UNCHANGED (runpy, same argv) with coarsen_slow.run wrapped by a cache, so
that the script's own output is reproduced line for line, and then prints what the script does not: the loss of EVERY
world in EVERY column at each f (2 decimals), the worlds at or over 1 point, the rank-order changes against f = 1, and
the Early Earth lines. Nothing existing is edited; no sweep is rerun (each f is computed once).

f = 17 x S is an ASSUMED product (sphere p90 factor 17 times the eddy slowdown S, h14b_advect.log); see
h13h14_combined_predictions.txt. Pre-registered there; the hash is asserted here and printed at the top of the log.

Run: cd <repo> && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy \
       --with numba --with matplotlib python h13h14_combined.py f=1,17,34,100,136,300,1000 \
       mix=h13_sphere_s16_f50.log > h13h14_combined.log 2> h13h14_combined.err
Smoke runs before the results run: selfcheck() below on made-up numbers only (no sweep run).
"""
import hashlib
import runpy
import sys

import numpy as np

PRED = "h13h14_combined_predictions.txt"
PRED_SHA = "5385115994a2be924acfc1f2525e3b1642063b829b3319626197472474d81167"
COLS = ("main no lid", "main final", "x10 no lid", "x10 final")
SCRIPT = "h13_sphere_shares.py"


def losses(R, base):
    """{world: [loss in points per column]} = base - R[f]."""
    return {n: [b - a for a, b in zip(R[n], base[n])] for n in base}


def order(R, k):
    """Worlds by share in column k, descending (stable on ties)."""
    return sorted(R, key=lambda n: -R[n][k])


def swaps(o0, o1):
    """Pairs (higher at f = 1, lower at f = 1) whose order is reversed in o1."""
    p0, p1 = {n: i for i, n in enumerate(o0)}, {n: i for i, n in enumerate(o1)}
    return [(a, b) for i, a in enumerate(o0) for b in o0[i + 1:] if p1[a] > p1[b]]


def report(cache, fs):
    base = cache[1.0]
    print("\n=== per-world loss against f = 1, points (base - f), by column; N = 200,000 so 1 point = 2000 draws ===")
    for k, c in enumerate(COLS):
        print(f"\n[{c}]  base share, then loss at f = " + ", ".join(f"{f:g}" for f in fs if f != 1.0))
        for n in base:
            print(f"  {n:18} {base[n][k]:6.2f}  " + "  ".join(f"{losses(cache[f], base)[n][k]:7.3f}" for f in fs if f != 1.0))
    print("\n=== worlds at or over 1.00 point, per f and column (main columns are the scored ones) ===")
    for f in fs:
        if f == 1.0:
            continue
        L = losses(cache[f], base)
        print(f"f = {f:g}: " + " | ".join(f"{c}: " + (", ".join(f"{n} {L[n][k]:.2f}" for n in L if L[n][k] >= 1.0) or "none")
                                          for k, c in enumerate(COLS)))
    print("\n=== rank order against f = 1 (pairs that swap: [higher at f = 1] now below [lower at f = 1]) ===")
    for f in fs:
        if f == 1.0:
            continue
        print(f"f = {f:g}: " + " | ".join(f"{c}: " + (", ".join(f"{a}<{b}" for a, b in swaps(order(base, k), order(cache[f], k))) or "none")
                                          for k, c in enumerate(COLS)))
    print("\n=== Early Earth ocean: share and loss in each column ===")
    for f in fs:
        n = "Early Earth ocean"
        print(f"f = {f:g}: share " + " / ".join(f"{cache[f][n][k]:.3f}" for k in range(4)) + "  loss " +
              " / ".join(f"{base[n][k] - cache[f][n][k]:.4f}" for k in range(4)))
    print("\nnot a combined result: the 'mix' line above is the sphere factor alone (every draw its own f_t).")


def selfcheck():
    a = {"X": [50.0, 50.0, 40.0, 40.0], "Y": [49.0, 49.0, 39.0, 39.0], "Z": [10.0, 10.0, 5.0, 5.0]}
    b = {"X": [47.5, 47.5, 38.0, 38.0], "Y": [49.0, 49.0, 39.0, 39.0], "Z": [10.0, 10.0, 5.0, 5.0]}
    assert abs(losses(b, a)["X"][0] - 2.5) < 1e-12 and losses(b, a)["Y"][0] == 0.0
    assert order(a, 0) == ["X", "Y", "Z"] and swaps(order(a, 0), order(b, 0)) == [("X", "Y")]
    assert swaps(order(a, 0), order(a, 0)) == []
    assert hashlib.sha256(b"x").hexdigest().startswith("2d711642")


if __name__ == "__main__":
    selfcheck()
    h = hashlib.sha256(open(PRED, "rb").read()).hexdigest()
    assert h == PRED_SHA, f"{PRED} does not match the recorded hash: {h}"
    print(f"# {PRED} sha256 {h} (asserted)", flush=True)
    import h13_sphere  # noqa: F401  (same import order as h13_sphere_shares.py: seedfix is applied by coarsen_slow -> lid_melt)
    import coarsen_slow as C
    cache, _run = {}, C.run

    def run(f, seed=0):
        out = _run(f, seed)
        if np.isscalar(f):
            cache[float(f)] = out
        return out

    C.run = run
    kw = dict(a.split("=") for a in sys.argv[1:])
    fs = [float(v) for v in kw["f"].split(",")]
    assert 1.0 in fs, "f = 1 is the base"
    sys.argv = [SCRIPT] + sys.argv[1:]
    runpy.run_path(SCRIPT, run_name="__main__")
    report(cache, fs)
