"""APED R5 follow-up 2: sigma against the activation rate a at fixed large drive x (activation effectively irreversible).
Scored against aped_r5_predictions_v3.txt (written after aped_r5_scan.py, whose free-x design could not reach small a).

Predictions hashed before this script existed or ran:
  sha256 bef0d7c99e480bbcdfc6aec492c2a809415a84c8f99228ee54518a19de4d503d, 2026-10-21 00:40 PDT (07:40 UTC)
Disclosure: no smoke run before the hash.  The x* slice (aped_r5.log "Model check") was known; W2/W3 are partly known from it.
Copied aped.py functions come from aped_r5.py (importing it runs nothing).  kH = 1, kA = a / XB, so a/kH = a and the reverse
activation is 1e-3 kA = 1e-7 a.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python aped_r5_scan2.py > aped_r5_scan2.log
"""
import hashlib
import math

import numpy as np
from scipy.optimize import brentq

from aped_r5 import BASE, racemic, lam

PRED_SHA = "bef0d7c99e480bbcdfc6aec492c2a809415a84c8f99228ee54518a19de4d503d"
XB = 1e4
AS = np.geomspace(1e-7, 1.0, 71)
CFG = {"A": dict(kR=0.0, kP=1.0), "B": dict(kR=1e-4, kP=1.0), "C": dict(kR=0.0, kP=22.0)}
KH_REAL = 0.018
LAB = {"L1 a_lo": 1.716e-7 / KH_REAL, "L1 a_hi": 1.925e-5 / KH_REAL, "L9 a_lo": 1.669e-4 / KH_REAL}  # aped_r5.log
Y = {}


def sig(c, a, y0=None):
    """Leading difference-mode Re(lambda) at x = XB, kA = a/XB; nan if no symmetric fixed point."""
    p = dict(BASE, kA=a / XB, **CFG[c])
    try:
        y = racemic(XB, p, y0)
    except RuntimeError:
        return float("nan"), None
    return lam(XB, p, y).real, y


def curve(c):
    s, y = [], None
    for a in AS:
        v, y2 = sig(c, a, y)
        y = y2 if y2 is not None else y
        s.append(v)
    return np.array(s)


def crossings(c, s):
    """Refined (a_c1, a_c2): lowest upward and next downward sign change of sigma over AS (None if absent)."""
    f = lambda la: sig(c, 10 ** la)[0]
    up = dn = None
    for i in range(len(AS) - 1):
        u, v = s[i], s[i + 1]
        if np.isnan(u) or np.isnan(v):
            continue
        r = 10 ** brentq(f, math.log10(AS[i]), math.log10(AS[i + 1]), xtol=1e-6) if u * v < 0 else None
        if r and u < 0 < v and up is None:
            up = r
        elif r and u > 0 > v and up is not None and dn is None:
            dn = r
    return up, dn


def main():
    with open("aped_r5_predictions_v3.txt", "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == PRED_SHA, "predictions file changed"
    print(f"sigma(a) at x = {XB:g} M, kH = 1 (a/kH = a), kA = a/x; reverse activation 1e-3 kA")
    S, C = {}, {}
    for c, cfg in CFG.items():
        S[c] = curve(c)
        C[c] = crossings(c, S[c])
        print(f"\nconfig {c}: kR = {cfg['kR']:g}, kP = {cfg['kP']:g}; skipped {int(np.isnan(S[c]).sum())} of {len(AS)}")
        for a, s in zip(AS[::5], S[c][::5]):
            print(f"  a/kH {a:9.3e}: sigma {s:+.4e}  sigma/a {s / a:+.3e}")
        print(f"  a_c1/kH {C[c][0]}  a_c2/kH {C[c][1]}  max sigma {np.nanmax(S[c]):+.4e} at a/kH {AS[np.nanargmax(S[c])]:.3e}")
    print("\nlab points (kH = 0.018 /s):")
    lab = {c: {k: sig(c, v)[0] for k, v in LAB.items()} for c in "AB"}
    for c in "AB":
        for k, v in LAB.items():
            print(f"  {c} {k} a/kH {v:.3e}: sigma {lab[c][k]:+.4e}")
    sB, _ = sig("B", 0.3852)
    ok = lambda r, lo, hi: r is not None and lo <= r <= hi
    res = [
        ("W1", abs(sB - 0.0024332) / 0.0024332 < 0.01, f"B sigma(0.3852) {sB:.5g} vs 0.0024332 (aped.log L89)"),
        ("W2", ok(C["B"][0], 0.005, 0.02), f"B a_c1/kH {C['B'][0]}"),
        ("W3", ok(C["A"][0], 8e-4, 3e-3), f"A a_c1/kH {C['A'][0]}"),
        ("W4", C["A"][0] is not None and C["C"][0] is not None and C["C"][0] < C["A"][0], f"C a_c1/kH {C['C'][0]} vs A {C['A'][0]}"),
        ("W5", lab["A"]["L1 a_hi"] < 0 and lab["A"]["L1 a_lo"] < 0, f"A sigma at L1 a_hi {lab['A']['L1 a_hi']:+.3e}, a_lo {lab['A']['L1 a_lo']:+.3e}"),
        ("W6", lab["A"]["L9 a_lo"] > 0 > lab["B"]["L9 a_lo"], f"L9 a_lo sigma A {lab['A']['L9 a_lo']:+.3e}, B {lab['B']['L9 a_lo']:+.3e}"),
    ]
    print()
    for n, k, why in res:
        print(f"{n} {'HELD  ' if k else 'FAILED'} {why}")
    print(f"{sum(k for _, k, _ in res)} of {len(res)} held (scorecard.py not changed).")


if __name__ == "__main__":
    main()
