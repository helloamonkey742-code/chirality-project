"""Ganymede and Callisto in the uncertainty.py sweep: a real percentage instead of a ranking (2026-10-06).

SUPERSEDED for Callisto by ocean_worlds.py V2 (literature check: Callisto's ocean is ~100 km, not 30 km -> ~29%).

Same model, same input ranges as Europa in uncertainty.py; only the size, ocean depth and the stratified-lid check change.
uncertainty.py's own numbers (Enceladus 22% / Europa 29% / Earth 75%) are untouched; this script reports beside them.

Inputs (2026-10-06 lookup):
  Ganymede  radius 2634 km. Ocean between ~150 and ~250 km depth, i.e. ~100 km thick (Saur et al. 2015 JGR Space
            Phys. 10.1002/2014ja020778, HST auroral ovals). Thicker oceans (several hundred km) are allowed by
            interior models (Vance et al. 2014 PSS 10.1016/j.pss.2014.03.011): depth sensitivity 50-300 km.
  Callisto  radius 2410 km. Ocean likely but thickness poorly constrained ("thick conductive ocean" favoured,
            Cochrane et al. 2025 AGU Adv. 10.1029/2024av001237). ASSUMED 30 km headline, sensitivity 10-100 km.
  Both: concentration, D_h, D_z, time ranges = Europa's (ASSUMED; nothing in the lookup sets them apart).
  Stratified lid: ice ~150 km thick -> ocean top at ~200 MPa (Ganymede, 920 kg/m^3 x 1.43 m/s^2 x 150 km) and
            ~170 MPa (Callisto). TEOS-10 (gsw): fresh water at its freezing point already sinks when cooled above
            28 MPa, so there is no thermal stratified lid at any salinity (strat_sweep.py's S* does not exist).
            TEOS-10 is only fitted to 100 MPa, so this is an extrapolation. Not modelled: salt-driven layering
            ("club sandwich", Vance et al. 2014), shown as an explicit "if pinned with probability p" line.
  Volume V = 4 pi R^2 H, matching uncertainty.depth() (area of the moon's surface, as for Europa).

PREDICTIONS (written before running, 2026-10-06):
 (P1) Ganymede at 100 km: about Europa's 29% (H ~ Europa's 123 km; the larger radius slows layer merging by ~1.8x,
      worth a few points) -> 22-29%. With no stratified lid that is also its final number, so Ganymede >= Europa's 22%.
 (P2) Deeper Ganymede ocean (300 km) lowers it (layers take longer to merge); 50 km raises it.
 (P3) Callisto at 30 km: higher than Ganymede (thin ocean, layers merge fast), ~35-45%; but its ocean is the least
      certain to exist, which this number does not include.
 (P4) Vertical mixing D_z is still the deciding input for both, as for Europa.

RESULTS (jupiter_moons.log, 2026-10-06): P1-P4 all held.
 Ganymede 28.5% (100 km ocean); 32.5% at 50 km, 21.8% at 300 km. No stratified lid, so it stays ~29%, above Europa's
 21.7% after its salinity lid. If salt layering pinned it half the time: 14%.
 Callisto 35.9% (30 km, ASSUMED); 42.9% at 10 km, 28.5% at 100 km; times the unknown chance its ocean exists at all.
 Both: D_z decides (low half 0%, high half 57% / 71%), exactly as for Europa.
"""
import numpy as np

import uncertainty as U

RANGES = U.WORLDS["Europa"][2:]           # concentration, D_h, D_z, time: same as Europa
MOONS = {"Ganymede": (2634e3, 100e3, (50e3, 300e3)), "Callisto": (2410e3, 30e3, (10e3, 100e3))}  # R, H, H sensitivity


def args(R, H):
    return (4 * np.pi * R ** 2 * H, np.pi * R, *RANGES)


def main():
    _, win, _ = U.sweep(*U.WORLDS["Europa"])  # first draw in the module RNG, as in uncertainty.py
    assert abs(win.mean() - U.README_WIN_PCT["Europa"]) < U.WIN_TOL, "uncertainty.py headline drifted"
    print(f"Europa (check): {win.mean():.1%}  (21.7% after strat_sweep.py's salinity lid)")
    out = {}
    for name, (R, H, (h_lo, h_hi)) in MOONS.items():
        assert abs(U.depth(*args(R, H)[:2]) - H) < 1, "depth() round trip"
        x, win, parts = U.sweep(*args(R, H))
        out[name] = win.mean()
        dz = x["vertical D_z"]
        lo = dz <= np.median(dz)
        print(f"\n{name}, ocean {H / 1e3:.0f} km: physics wins in {win.mean():.1%}"
              f"  (D_z low/high half {win[lo].mean():.0%} / {win[~lo].mean():.0%})")
        print("  each condition alone: " + ", ".join(f"{k} {v:.0%}" for k, v in parts.items()))
        for h in (h_lo, h_hi):
            print(f"  ocean {h / 1e3:.0f} km: {U.sweep(*args(R, h))[1].mean():.1%}")
        print("  if salt layering pins it with probability p: " +
              ", ".join(f"p={p}: {win.mean() * (1 - p):.1%}" for p in (0.25, 0.5, 0.75)))
    assert 0.22 <= out["Ganymede"] <= 0.30 and out["Callisto"] > out["Ganymede"], "RESULTS above drifted"
    return out


if __name__ == "__main__":
    main()
