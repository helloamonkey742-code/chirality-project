"""Bound on the advect2d.py finding: eddies larger than the wall width cleared walls >= 4x more slowly than Eq. 17
gives (exploratory, advect2d.py RESULTS). Here the sideways coarsening time (L/A)^2/D_h is multiplied by f = 4 and 10
for EVERY draw (A = 5.5/sqrt(f) in uncertainty.py, sweep2.py, sweep3.py and multiwall.py, which all read U.A), a worst
case: how many sampled inputs have eddies larger than the wall width is not known. The cascade x10 of Table 3 scales
only the rise/sink times, not this term (multiwall.py:421), so the two are run together here (2026-10-14).

Base: multiwall.py Section 4, seed 0, rules main and x10, no lid and final; Enceladus with the seedfix.py correction
(imported from lid_melt.py, which also holds the logged values).
Run: nice -n 10 /opt/anaconda3/bin/python coarsen_slow.py > coarsen_slow.log

PREDICTIONS (written 2026-10-14, before any run):
 (P1) f = 1 reproduces multiwall.log / seedfix.log (main and x10, no lid and final) within 0.05 points.
 (P2) f = 4: every world loses 0-6 points (main, final); f = 10: 0-10 points.
 (P3) Larger worlds lose more (L is larger): Ganymede loses more points than Europa at f = 4 and 10.
 (P4) Ganymede > Europa > Enceladus (final) at f = 4 and 10, main and x10.
 (P5) Early Earth ocean stays >= 95% at f = 10 (main).

RESULTS (coarsen_slow.log, 2026-10-14): the script prints 5 of 5; counted 4 of 5 (P3 below). Before the predictions
were committed, an N = 2,000 smoke run and a diagnostic (largest (L/A)^2/D_h / age over all draws) were seen; the
predictions were not changed.
 (P1) held, 0.05.  (P2) held: largest loss 0.0 points at f = 4 and 0.1 at f = 10.  (P3) the check passes on
 Ganymede 0.1 vs Europa 0.0 at f = 10, a difference below the printed rounding: counted as NOT held.
 (P4) held, 44.1 / 30.1 / 0.0 (f = 4).  (P5) held, 99.9.
 Why nothing moves: sideways coarsening takes at most 7% of the ocean's age in any draw of any world (Ganymede;
 median 0.01%), so it must be about 14x slower before it decides a single draw. The advect2d slowdown (>= 4x) is far
 inside that margin. The bound is for this model's ages (10^8-10^9 yr, Mimas 2-25 Myr), not for short-lived pools.
"""
import numpy as np

import lid_melt as LM   # applies the seedfix RNG correction before the project modules load
import uncertainty as U
import sweep2 as S
import sweep3 as S3   # noqa: F401  (multiwall imports it)
import multiwall as M

FS = (1.0, 4.0, 10.0)
A0 = U.A


def run(f, seed=0):
    U.A = A0 / np.sqrt(f)
    _, W = S.chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        out[name] = tuple(100 * M.multi(x, a[1], a[4], ex, l, scale=sc)[0].mean() for sc in (1.0, 10.0) for l in (none, salt))
    U.A = A0
    return out   # name: (main no lid, main final, x10 no lid, x10 final)


def main():
    R = {f: run(f) for f in FS}
    names = list(R[1.0])
    print(f"coarsen_slow.py, N = {U.N}, seed 0. P_pved %: main no lid / final | x10 no lid / final, at coarsening time x f")
    for n in names:
        print(f"{n:18} " + "   ".join(f"f {f:>4g}: {r[0]:5.1f}/{r[1]:5.1f} | {r[2]:5.1f}/{r[3]:5.1f}"
                                      for f, r in ((f, R[f][n]) for f in FS)))
    C = []

    def check(p, ok, d):
        C.append(ok)
        print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {d}")

    print("\nPredictions:")
    d = max(abs(a - b) for n in names for a, b in zip(R[1.0][n], LM.LOGGED[n]))
    check("P1", d <= 0.05, f"max |f=1 - logged| = {d:.2f}")
    loss = {f: {n: R[1.0][n][1] - R[f][n][1] for n in names} for f in (4.0, 10.0)}
    check("P2", all(0 <= v <= 6 for v in loss[4.0].values()) and all(0 <= v <= 10 for v in loss[10.0].values()),
          "max loss (main, final): f 4 " + f"{max(loss[4.0].values()):.1f} ({max(loss[4.0], key=loss[4.0].get)}), "
          f"f 10 {max(loss[10.0].values()):.1f} ({max(loss[10.0], key=loss[10.0].get)})")
    check("P3", all(loss[f]["Ganymede"] > loss[f]["Europa"] for f in (4.0, 10.0)),
          "loss Ganymede / Europa: " + ", ".join(f"f {f:g} {loss[f]['Ganymede']:.1f}/{loss[f]['Europa']:.1f}" for f in (4.0, 10.0)))
    check("P4", all(R[f]["Ganymede"][k] > R[f]["Europa"][k] >= R[f]["Enceladus"][k] for f in (4.0, 10.0) for k in (1, 3)),
          "final: " + ", ".join(f"f {f:g} {R[f]['Ganymede'][1]:.1f}/{R[f]['Europa'][1]:.1f}/{R[f]['Enceladus'][1]:.1f}"
                                for f in (4.0, 10.0)))
    check("P5", R[10.0]["Early Earth ocean"][0] >= 95, f"early Earth f 10: {R[10.0]['Early Earth ocean'][0]:.1f}")
    print(f"\n{sum(C)} of {len(C)} held.")


if __name__ == "__main__":
    main()
