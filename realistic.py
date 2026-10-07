"""More realistic currents: do thinner walls need stronger currents to break free? Plus an Earth-like three-zone ocean.

currents.py freed a wall under a 10x slower lid once the lid Peclet number Pe = W * 0.4 / D_lid passed ~1, and that is
why Earth's 75% stands (thermocline Pe ~ 5-10). Its biggest ceiling: the hand wall (width ~ sqrt(2 D)) was thick next
to the lid. In a real ocean the wall width is sqrt(D / amplifier rate), metres to tens of metres against a ~1 km
thermocline, so walls are far thinner. Here D_hi drops 2x and 4x (walls 1.4x and 2x thinner relative to the lid),
the grid is refined to keep 7 cells per wall, and the threshold Pe is re-measured. Same model and flow as currents.py
(polar_gap.run). dt = 0.4 (checked: same outcomes as dt = 0.1, takeover times 12-14% later).

PREDICTIONS (written before any run, 2026-10-06):
 (P1) The pin's grip is a speed ~ D_lid * d(ln sigma)/dr, set by how sharply tension jumps across the wall, i.e. ~1/wall
      width; the current must beat it, so the threshold Pe grows as walls thin, about x1.4 and x2 here. Base walls:
      threshold between Pe 0.4 and 1.2 (currents.py). Walls 2x thinner: stuck at 1.2, freed at 4. Low confidence;
      the alternative (threshold stays ~1, set by the lid alone) would also fit currents.py.
 (P2) If P1 holds, the threshold keeps rising for real-ocean walls (10-100x thinner than the lid), so Earth's Pe 5-10
      may NOT be enough and the 75% would need a different rescue (or a slow amplifier with wide walls).
 (P3) Earth-like column: fast mixed layer (top 0.1) over the lid (0.3 thick) over the abyss. Without currents: pinned
      (two_zone.py P4). With Pe 4: freed, deep hand wins, as with no mixed layer.
 (P4, added 2026-10-06 after section 2 showed the mixed layer raises the threshold to Pe 12-40, before section 3 ran)
      Earth's thermocline opens to the surface in the Southern Ocean, so water can pass it without crossing the slow
      layer. Opening a 20 deg polar gap (below the 25-30 deg that frees the wall by mixing alone, polar_gap.py) in the
      three-zone column: Pe 4 now frees the wall in both current directions. Low confidence.
 (P5, added after section 3, before section 4 ran) With the outcrop the current goes around the slow layer instead of
      through it, so thinner walls (2x) should NOT raise the threshold: Pe 1.2 still frees the wall, both directions.
 Not run: an Enceladus-like 1000x lid; its lid wall (width 0.0045) is below the grid spacing, so the run would not be
 resolved (ratio 100 already is only marginally).

RESULTS (realistic.log, 2026-10-06):
 P1 held: freeing threshold Pe between 0.4 and 1.2 for base and 1.4x thinner walls, between 1.2 and 4 for 2x thinner
 walls (the 2x-thinner Pe 1.2 run lost half its + water, 0.37 -> 0.18, then stopped). Thinner walls grip harder.
 P3 WRONG: a fast mixed layer above the lid raises the threshold to between Pe 12 and 40 (Pe 4 and 12 stuck, 40 freed;
 dt = 0.1 gives the same Pe 4 share, 0.250). So P2's worry is real for a closed lid: Earth's Pe 5-10 would not be enough.
 P4 held: with a 20 deg polar outcrop (too small to free the wall by mixing alone; stuck with no current) Pe 1.2 and 4
 free it in both directions. NEW: the direction picks the hand. Sinking at the pole -> the TOP hand takes the ocean;
 rising at the pole -> the deep hand.
 P5 half wrong: with 2x thinner walls, sinking at the pole still frees it at Pe 1.2 (top hand wins), rising at the pole
 stays stuck (share 0.31).
 Meaning for Earth: its deep overturning SINKS at the poles where the thermocline outcrops, the direction that frees
 the wall at the lowest Pe tested even with thinner walls, and hands the ocean to the top hand, which matches the
 "top layer decides" assumption of uncertainty.py. Earth's 75% survives in this model, on the outcrop + sinking branch.
 Ceilings: one flow cell, axisymmetric, R/depth 4, walls still ~5-50x too thick relative to a real thermocline, and a
 diffusive (not adiabatic, along-isopycnal) outcrop.
"""
import time

from polar_gap import run, D_HI

DT = 0.4
LEVELS = ((1.0, "base walls"), (2.0, "walls 1.4x thinner"), (4.0, "walls 2x thinner"))   # D_hi divisor


def show(label, ratio, pe, div, **kw):
    t0 = time.time()
    d_hi = D_HI / div
    lid = 0.4 - kw.get("top", 0.0)
    W = pe * (d_hi / ratio) / lid
    out, t, frac, _, _ = run(ratio, kw.pop("gap", 0), W=W * kw.pop("sign", 1), dt=DT, t_max=8000.0 * div, d_hi=d_hi, grid=div ** 0.5, **kw)
    print(f"  {label:28s} ratio {ratio:<3g} Pe {pe:5.1f} (W {W:.2e}) -> {out:13s} t={t:7.0f}  + share {frac:.3f}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
    return out


def main():
    t0 = time.time()
    print(f"(1) threshold Pe vs wall thickness, lid 10x slower, dt={DT}")
    for div, label in LEVELS:
        for pe in ((0.4, 1.2) if div == 1 else (0.4, 1.2, 4.0)):
            show(label, 10, pe, div)
    print("\n(2) Earth-like column: mixed layer 0.1 / lid 0.3 (10x slower) / abyss, base walls")
    show("three zones, no current", 10, 0.0, 1, top=0.1)
    for pe in (4.0, 12.0, 40.0):
        show("three zones, current", 10, pe, 1, top=0.1)
    print("\n(3) Earth-like column with a 20 deg polar outcrop of the lid, base walls")
    show("three zones + gap, no current", 10, 0.0, 1, top=0.1, gap=20)
    for pe in (1.2, 4.0):
        for sign in (1, -1):
            show(f"three zones + gap, {'sink' if sign > 0 else 'rise'} at pole", 10, pe, 1, top=0.1, gap=20, sign=sign)
    print("\n(4) same with walls 2x thinner")
    for sign in (1, -1):
        show(f"three zones + gap, {'sink' if sign > 0 else 'rise'} at pole", 10, 1.2, 4, top=0.1, gap=20, sign=sign)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
