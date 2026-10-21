"""APED R5 follow-up: does the instability window survive slow activation at ANY drive x?
Scored against aped_r5_predictions_v2.txt (written after aped_r5.py's run, because P10 failed there).

Predictions hashed before this script existed or ran:
  sha256 06e8573cba32dd0e5c500529a5c5cbbdf832adfa4374f77be47489c953b5cb00, 2026-10-21 00:36 PDT (07:36 UTC)
Disclosure: no smoke run before the hash.  The x* slice of the answer (aped_r5.log "Model check") was known; V2/V3's open
side and V4's centre are partly known from it (stated in the predictions file).
Copied aped.py functions come from aped_r5.py (importing it runs nothing).  smax is aped.py's, with the grid as a constant
and a count of grid points where no symmetric fixed point was found (they are skipped, not hidden).
Real-rate inputs (read earlier): kH = 0.018 /s (Devaraj 2016, doi:10.1021/jacs.6b08801, Fig. S7); kP >= ~0.4 /M/s
(Leman 2004, note 9) so kP/kH >= ~22 /M.  Lab a from aped_r5.log (L1 a_hi, L9 a_lo).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python aped_r5_scan.py > aped_r5_scan.log
"""
import hashlib
import math

import numpy as np
from scipy.optimize import minimize_scalar

from aped_r5 import BASE, racemic, lam

PRED_SHA = "06e8573cba32dd0e5c500529a5c5cbbdf832adfa4374f77be47489c953b5cb00"
XS = np.geomspace(1e-4, 1e4, 65)
KA = (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
CFG = {"A": dict(kR=0.0, kP=1.0), "B": dict(kR=1e-4, kP=1.0), "C": dict(kR=0.0, kP=22.0)}
KH_REAL = 0.018
A_HI_L1, A_LO_L9 = 1.925e-5, 1.669e-4  # aped_r5.log
CACHE = {}


def smax(p):
    """Largest Re(lambda) over XS (continuation), refined by a bounded search; returns (value, x, n_failed)."""
    y, L, ys, bad = None, [], [], 0
    for x in XS:
        try:
            y = racemic(x, p, y)
            L.append(lam(x, p, y))
        except RuntimeError:
            bad += 1
            L.append(complex(np.nan))
        ys.append(y)
    L = np.array(L)
    i = int(np.nanargmax(L.real))
    if L[i].real <= 0:
        return L[i].real, XS[i], bad
    f = lambda lx: -lam(10 ** lx, p, racemic(10 ** lx, p, ys[i])).real
    r = minimize_scalar(f, bounds=(math.log10(XS[max(i - 1, 0)]), math.log10(XS[min(i + 1, len(XS) - 1)])),
                        method="bounded", options=dict(xatol=1e-5))
    return -r.fun, 10 ** r.x, bad


def sm(c, kA):
    if (c, kA) not in CACHE:
        CACHE[c, kA] = smax(dict(BASE, kA=kA, **CFG[c]))
    return CACHE[c, kA]


def closing(c):
    """a_c/kH: log-bisection (6 steps) of the largest closed kA below the lowest open grid kA."""
    op = [k for k in KA if sm(c, k)[0] > 0]
    if not op or len(op) == len(KA):
        return None
    hi = min(op)
    lo = max(k for k in KA if k < hi)
    for _ in range(6):
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if sm(c, mid)[0] > 0 else (mid, hi)
    return math.sqrt(lo * hi) * sm(c, hi)[1], lo, hi


def main():
    with open("aped_r5_predictions_v2.txt", "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == PRED_SHA, "predictions file changed"
    assert abs(lam(0.3852, BASE, racemic(0.3852, BASE)).real - 0.0024332) < 2e-7  # aped.log L89
    print("sigma_max = max over x in 1e-4..1e4 M of the leading difference-mode Re(lambda); kH = 1, a = kA x")
    for c, cfg in CFG.items():
        print(f"\nconfig {c}: kR = {cfg['kR']:g}, kP = {cfg['kP']:g}")
        for k in KA:
            s, x, bad = sm(c, k)
            a = k * x
            print(f"  kA {k:7.0e}: sigma_max {s:+.4e} /s at x {x:.4g} M  a/kH {a:.3e}  sigma/a {s / a:+.3e}  "
                  f"{'open' if s > 0 else 'closed'}  skipped {bad}")
    ac = {c: closing(c) for c in CFG}
    for c, v in ac.items():
        print(f"\nconfig {c} closing: " + (f"a_c/kH {v[0]:.3e}  (kA between {v[1]:.3e} and {v[2]:.3e})" if v else "no open/closed boundary on the grid"))
    ahi, alo9 = A_HI_L1 / KH_REAL, A_LO_L9 / KH_REAL
    print(f"\nlab a/kH with kH = {KH_REAL}: plain COS a_hi {ahi:.3e} (a_lo {1.716e-7 / KH_REAL:.3e}); ferricyanide L9 a_lo {alo9:.3e}")
    o = lambda c, k: sm(c, k)[0] > 0
    s1, x1, _ = sm("B", 1.0)
    res = [
        ("V1", abs(s1 - 0.0024332) < 1e-6 and abs(x1 - 0.3852) < 7e-4, f"sigma_max {s1:.5g} at x {x1:.4g} (aped.log L89 0.0024332 at 0.3852)"),
        ("V2", all(o("A", k) for k in (1.0, 1e-1, 1e-2)) and not any(o("A", k) for k in (1e-4, 1e-5)),
         "A open: " + ", ".join(f"{k:.0e}={o('A', k)}" for k in KA)),
        ("V3", all(o("B", k) for k in (1.0, 1e-1)) and not any(o("B", k) for k in KA[2:]),
         "B open: " + ", ".join(f"{k:.0e}={o('B', k)}" for k in KA)),
        ("V4", bool(ac["A"] and ac["B"]) and 3e-4 <= ac["A"][0] <= 6e-3 and 3e-3 <= ac["B"][0] <= 1e-1,
         f"a_c/kH A {ac['A'] and ac['A'][0]:.3g}, B {ac['B'] and ac['B'][0]:.3g}" if ac["A"] and ac["B"] else "no boundary"),
        ("V5", bool(ac["A"] and ac["C"]) and ac["C"][0] < ac["A"][0],
         f"a_c/kH kP=22 {ac['C'] and ac['C'][0]:.3g} vs kP=1 {ac['A'] and ac['A'][0]:.3g}" if ac["A"] and ac["C"] else "no boundary"),
        ("V6", bool(ac["A"]) and ahi < ac["A"][0] < alo9, f"a_hi(L1)/kH {ahi:.3g} < a_c/kH {ac['A'] and ac['A'][0]:.3g} < a_lo(L9)/kH {alo9:.3g}" if ac["A"] else "no boundary"),
    ]
    print()
    for name, ok, why in res:
        print(f"{name} {'HELD  ' if ok else 'FAILED'} {why}")
    print(f"{sum(ok for _, ok, _ in res)} of {len(res)} held (scorecard.py not changed).")


if __name__ == "__main__":
    main()
