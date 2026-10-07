"""Two-zone ocean (Soderlund, 2026-10-06 email): fast vertical mixing in a convecting region, slow in a stratified layer.

Same local rule as sphere_layers.py, radial only (sphere_layers.py showed D_h does not move a spherical layer wall):
    da/dt = a - a^3 + (1/r^2) d/dr (r^2 D(r) da/dr),   no flux at the sea floor and the ice.
D(r) = D_hi in the convecting zone, D_lo = D_hi/ratio in the stratified layer, joined by a tanh of width w.
Start: + hand in a thin layer under the ice, - hand below. The wall sinks; the run asks whether + takes the whole depth.

PREDICTIONS (written before any run, 2026-10-06):
 (P1) With one D everywhere (ratio 1) the wall sinks at 2D/r, as in sphere_layers.py, so the takeover time is
      T = integral of r/(2 D(r)) dr. If the zones only change the local speed, T follows the same integral
      with the two-zone D(r): the slow layer adds (its thickness) * r / (2 D_lo), a thickness-weighted harmonic mean.
 (P2) Stratified layer at the BOTTOM: the wall crosses the fast zone, enters the slow layer, finishes. T matches P1.
 (P3) Stratified layer on TOP: the wall crosses the slow layer, then must enter the fast zone below. The equation is a
      gradient flow of E = 4 pi r^2 sigma per wall, with wall tension sigma ~ sqrt(D). Moving inward shrinks the area
      (helps) but raises sigma (hurts); the wall passes only if |d ln sigma / dr| < 2/r. A ratio-10 jump over a few
      wall widths breaks that by ~100x, so the wall should STOP at the interface and two hands stay forever.
      Guess for the onset at w = 0.02 (smoothed by the wall's own width ~0.1-0.2): ratio ~1.03-1.05.
 (P4) Slow layer in the MIDDLE: the wall enters it from above (allowed) and stops at its lower edge (as P3).

RESULTS (two_zone.log, 2026-10-06):
 P1 held (control 1.04x predicted; dt/2 within 3%). P3 and P4 held: with the slow layer on top or in the middle the
 wall stopped ~0.08 above the slow layer's lower edge and never moved again (ratio 10 and 100, run to 3x the
 predicted time). Pinning onset at w=0.02 is between ratio 1.02 (passed, 1.38x slower) and 1.05 (stuck), inside the guess.
 P2 was WRONG in detail: with the slow layer at the bottom the takeover came 0.67-0.71x the predicted time, faster,
 because the wall is pulled INTO low-mixing water (lower tension) on its way down. Post-hoc (2b), not predicted:
 with a smooth joint (w = 0.1, 0.3) the wall slid UP toward the low-mixing water under the ice; at w = 0.3 the top
 hand was squeezed out against the ice and the deep hand took the whole ocean. Walls move toward weak mixing.
"""
import sys
import time

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

R, H, DZ, DT = 16.0, 1.0, 0.002, 0.1
D_HI = 2.5e-3


def run(ratio, zone, w=0.02, layer=0.4, t_max=None, dt=DT, h=0.0):
    """zone: 'top' | 'bottom' | 'middle' | 'none'; h: bias favouring the top (+) hand. Returns (times, wall radii, T_pred, takeover time or None)."""
    n = int(round(H / DZ))
    r_in = R - H
    rc = r_in + (np.arange(n) + 0.5) * DZ
    rf = r_in + np.arange(1, n) * DZ
    s = lambda x: 0.5 * (1 + np.tanh(x / w))                             # 0 below x=0, 1 above
    if zone == "top":
        slow = s(rf - (R - layer))
    elif zone == "bottom":
        slow = s((r_in + layer) - rf)
    elif zone == "middle":
        mid = r_in + H / 2
        slow = s(rf - (mid - layer / 2)) * s((mid + layer / 2) - rf)
    else:
        slow = np.zeros_like(rf)
    lnD = np.log(D_HI) - slow * np.log(ratio)                            # blend in log D
    Df = np.exp(lnD)
    up = Df * rf ** 2 / (DZ ** 2 * rc[:-1] ** 2)                         # flux coupling i <-> i+1, seen from i
    dn = Df * rf ** 2 / (DZ ** 2 * rc[1:] ** 2)                          # seen from i+1
    main = np.zeros(n); main[:-1] -= up; main[1:] -= dn
    A = diags([main, up, dn], [0, 1, -1], format="csc")
    lu = splu((diags(np.ones(n)) - dt * A).tocsc())

    r0 = R - 0.2                                                         # ~3 wall widths below the ice (a first run at R-0.1 let the no-flux ice eat the + layer)
    a = np.tanh((rc - r0) / np.sqrt(2 * D_HI / ratio if zone == "top" else 2 * D_HI))
    end = r_in + 3 * np.sqrt(2 * D_HI)                                   # "done": wall within 3 widths of the floor
    Dc = np.exp(np.interp(rc, rf, lnD))
    sel = (rc < r0) & (rc > end)
    T_pred = float(np.sum(rc[sel] / (2 * Dc[sel])) * DZ)
    t_max = t_max or 3 * T_pred
    t, ts, rw, done = 0.0, [], [], None
    nrec = max(1, int(T_pred / 200 / dt))
    k = 0
    while t < t_max:
        a = lu.solve(a + dt * (a - a ** 3 + h))
        t += dt; k += 1
        if k % nrec == 0:
            i = np.flatnonzero(a > 0)
            wall = rc[i[0]] if i.size else R
            ts.append(t); rw.append(wall)
            if wall < end and done is None:
                done = t
                break
    return np.array(ts), np.array(rw), T_pred, done


def show(label, ratio, zone, **kw):
    t0 = time.time()
    ts, rw, Tp, done = run(ratio, zone, **kw)
    fin = rw[-1]
    q = [rw[min(len(rw) - 1, int(f * len(rw)))] for f in (0.25, 0.5, 0.75)]
    res = f"took over at t={done:.0f}, {done / Tp:.3f} x predicted" if done else f"NOT taken over by t={ts[-1]:.0f} (3x predicted)"
    print(f"  {label:34s} ratio {ratio:<6g} T_pred {Tp:9.0f}  wall depth 1/4,1/2,3/4,end of run: "
          + " ".join(f"{R - x:.3f}" for x in q + [fin]) + f"  -> {res}  [{time.time() - t0:.0f}s]", flush=True)
    return done, Tp, R - fin


def main():
    t0 = time.time()
    sec = sys.argv[1:] or ["1", "2"]
    if "1" in sec:
        print(f"(1) R={R}, depth {H}, slow layer 0.4 thick, D_hi={D_HI}, sharp joint w=0.02, dz={DZ}, dt={DT}")
        c, Tp, _ = show("control: one zone", 1, "none")
        assert c and abs(c / Tp - 1) < 0.15, "control wall does not sink at 2D/r"
        c2, _, _ = show("control: dt/2", 1, "none", dt=DT / 2)
        assert abs(c2 / c - 1) < 0.03, "time step not converged"
        for ratio in (10, 100):
            for zone in ("bottom", "middle", "top"):
                done, _, _ = show(f"slow layer at {zone}", ratio, zone)
                assert (done is not None) == (zone == "bottom"), f"{zone} ratio {ratio}: pinning result changed"
    if "2" in sec:
        print("\n(2) onset of pinning, slow layer on top, w=0.02")
        for ratio in (1.02, 1.05, 1.1, 1.3, 2, 4):
            done, _, _ = show("slow layer at top", ratio, "top")
            assert (done is not None) == (ratio < 1.03), f"ratio {ratio}: pinning onset moved"
        print("\n(2b) wider joint, ratio 10")
        for w in (0.1, 0.3):
            show(f"slow layer at top, w={w}", 10, "top", w=w)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
