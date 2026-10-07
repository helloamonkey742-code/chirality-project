"""Currents: can an overturning circulation carry a pinned wall through a stratified lid? (Soderlund's "organized flows")

polar_gap.py showed that mixing alone pins the wall under a lid unless the polar gap is huge (60-87 deg for Earth), so
Earth's 75% hinges on currents. Same model as polar_gap.py (R = 4, depth 1, lid 0.4 thick everywhere, no gap) plus an
axisymmetric overturning current: sinking at the pole (speed ~2W), equatorward along the floor, rising at low latitudes
(speed ~W) THROUGH the lid, poleward under the ice. Exactly divergence-free on the grid (fluxes from a streamfunction).
Lid Peclet number Pe = W * 0.4 / D_lid: how strongly the current crosses the lid compared with mixing across it.

PREDICTIONS (written before any run, 2026-10-06):
 (P1) W = 0: pinned (polar_gap.py). No lid + current: one hand wins.
 (P2) Advection moves water, not walls' tension: a wall is carried with the water whatever the mixing jump. So the lid
      stops pinning once the current crosses it faster than mixing does: Pe << 1 (W = 1e-4) stays pinned, Pe >> 1
      (W >= 1e-2) ends one-handed. Threshold Pe ~ 1 (W ~ 2.5e-3 at ratio 10), within a factor ~3.
 (P3) Which hand wins once unpinned: low confidence. Guess the top (+) hand, because the sinking branch injects it into
      the deep while the rising branch only brings - water into the lid, where the slow mixing keeps it thin.
 (P4) Ratio 100 (lid 100x slower): same rule in Pe, so W = 1e-3 (Pe 4) already unpins.

RESULTS (currents.log, 2026-10-06):
 P1 held (pinned without current). P2 held for lids 10x and 30x slower: pinned at Pe 0.04 and 0.4, unpinned at Pe 1.2 and
 above (dt/2 agrees: 135 vs 140). P3 WRONG: every unpinned run (and the no-lid control, and the reversed current) ended
 with the DEEP hand, the volume majority (~62%), not the top hand. Stirring hands the ocean to whichever hand holds more
 water; without currents, the sphere's curvature hands it to the top. P4 WRONG: a lid 100x slower stays pinned up to
 Pe 80 (W = 2e-2, the strongest the grid resolves); the + hand shrinks (0.30 -> 0.16 -> 0.06 of the ocean) but keeps a
 pocket in the lid at the pole. Why P2's reasoning failed there: the lid is fixed in space, not in the water, so water can
 flow THROUGH a pinned wall and be converted as it crosses; the pin acts like static friction that a weak flow cannot beat.
 Ceilings: one overturning shape, axisymmetric, R/depth = 4; the wall width / lid thickness ratio (0.35 here) is far larger
 than in a real ocean, so real pinning may be stronger.
"""
import time

from polar_gap import run, DT, D_HI


def show(label, ratio, W, **kw):
    t0 = time.time()
    out, t, frac, _, _ = run(ratio, 0, W=W, t_max=8000.0, **kw)
    pe = W * 0.4 / (D_HI / ratio) if kw.get("lid", True) else float("nan")
    print(f"  {label:26s} ratio {ratio:<4g} W {W:<7g} Pe {pe:7.2f} -> {out:13s} t={t:6.0f}  + share {frac:.3f}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
    return out


def main():
    t0 = time.time()
    print("(1) controls")
    assert show("lid, no current", 10, 0) == "STUCK"
    show("no lid, current", 1, 1e-2, lid=False)
    print("\n(2) lid 10x slower, current strength sweep")
    for W in (1e-4, 1e-3, 3e-3, 1e-2, 3e-2):
        show("lid + current", 10, W)
    show("lid + current, dt/2", 10, 1e-2, dt=DT / 2)
    show("lid + current reversed", 10, -1e-2)
    print("\n(3) stronger lids")
    for ratio, W in ((30, 3e-3), (30, 1e-2), (100, 1e-3), (100, 1e-2), (100, 2e-2)):
        assert show("lid + current", ratio, W) == ("STUCK" if ratio == 100 else "DEEP HAND WON")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
