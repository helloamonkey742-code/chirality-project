"""h13_sphere_shares: how far would the world shares move if case (i) used a sideways-time factor f_t > 1? (H13, D3)

f_t = T1 / tau_sw from the h13_sphere logs is the factor on the sweep's sideways term (L/A)^2/D_h. coarsen_slow.run(f) sets
U.A = 5.5/sqrt(f), which multiplies that term by f in uncertainty.py, sweep2.py, sweep3.py and multiwall.py, so it is the
same quantity; nothing is edited, coarsen_slow is imported. Also printed, per world, from the same draws (seed 0):
  worst = the largest (L/A)^2/D_h / age over all draws (a draw's sideways test fails only when f * ratio > 1, so 1/worst
  is the f at which the first draw flips), and the share of draws with f * ratio > 1 at each f.
Pre-registered: D3 in h13_sphere_predictions.txt (hash checked through h13_sphere.need_pred).
Also "mix": every draw gets its own f_t, resampled from the T1/tau_sw values of a h13_sphere log (a run that never
reached one hand counts as 30, its horizon), the physical analogue of each ocean starting from a random patchwork.
Run: nice -n 10 /opt/anaconda3/bin/python h13_sphere_shares.py f=1,8,17,23,30 mix=h13_sphere_s16_f50.log > h13_sphere_shares.log
Smoke runs before the results run: none.
"""
import sys

import numpy as np

import h13_sphere as H
import coarsen_slow as C
import lid_melt as LM
import uncertainty as U
import sweep2 as S


def ratio(L, Dh, t):
    """Sideways time over age for the sweep's rule (L/A)^2/D_h <= t; 1.0 means exactly at the threshold."""
    return (L / U.A) ** 2 / Dh / t


if __name__ == "__main__":
    assert abs(ratio(5.5, 1.0, 1.0) - 1) < 1e-12 and abs(ratio(11.0, 1.0, 2.0) - 2) < 1e-12  # hand-computed
    H.need_pred()
    kw = dict(a.split("=") for a in sys.argv[1:])
    fs = [float(v) for v in kw["f"].split(",")]
    _, W = S.chains(0)
    r = {n: ratio(a[1], x["horizontal D_h"], x["time"]) for n, (x, win, a, salt) in W.items()}
    print("sideways time / age over all draws (seed 0): worst, first-flip f = 1/worst, share of draws with f*ratio > 1 at f = "
          + ", ".join(f"{f:g}" for f in fs))
    for n, v in r.items():
        print(f"  {n:18} worst {v.max():.4f}  1/worst {1 / v.max():8.1f}   " + "  ".join(f"{100 * (f * v > 1).mean():6.3f}%" for f in fs))
    R = {f: C.run(f) for f in fs}
    base = R[1.0] if 1.0 in R else C.run(1.0)
    if "mix" in kw:
        T = [float(l.split("T1/tau_sw=")[1].split()[0]) for l in open(kw["mix"]) if l.startswith("run ")]
        T = np.nan_to_num(np.array(T), nan=30.0)
        R["mix"] = C.run(np.random.default_rng(12345).choice(T, U.N))
        print(f"mix: {len(T)} values from {kw['mix']}, median {np.median(T):.2f}, mean {T.mean():.2f}")
    d = max(abs(a - b) for n in base for a, b in zip(base[n], LM.LOGGED[n]))
    print(f"\nf = 1 against the logged values (lid_melt.LOGGED): max |difference| {d:.2f} points")
    print("P_pved %: main no lid / final | x10 no lid / final")
    for n in base:
        print(f"{n:18} " + "   ".join(f"f {f:>4}: {R[f][n][0]:5.1f}/{R[f][n][1]:5.1f} | {R[f][n][2]:5.1f}/{R[f][n][3]:5.1f}" for f in R))
    for f in R:
        loss = {n: [b - a for a, b in zip(R[f][n], base[n])] for n in base}
        k = max(range(4), key=lambda i: max(v[i] for v in loss.values()))
        w = max(loss, key=lambda n: loss[n][k])
        print(f"f = {f}: largest loss against f = 1 is {loss[w][k]:.2f} points ({w}, column {('main no lid', 'main final', 'x10 no lid', 'x10 final')[k]})")
