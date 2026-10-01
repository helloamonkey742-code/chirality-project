"""Uncertainty sweep: vary every uncertain input at once; which conclusions survive, and which
input decides the answer?

For each world, draw 200,000 random input sets (log-uniform over each plausible range) and
count "physics wins" only if ALL hold:
  finishes      k2*c*tau >= 10                                   (ocean.py time rule)
  one ocean     patchwork heals in time, both ways:              (coarsen3d.py, A = 5.5)
                  across:  (L/A)^2/D_h <= tau   (horizontal mixing, L = half circumference)
                  down:    (H/A)^2/D_z <= tau   (vertical mixing, H = ocean depth = V / area)
                (Correction 2026-09-25: this used one D over L, with Zeng & Jansen's VERTICAL range
                 10^-10..10^-3 m^2/s. Mixing is anisotropic, so each direction gets its own D.)
  strong enough F * Delta_whole >= 2                             (ocean.py x real-chemistry factor F)
Layering (added 2026-09-25, aniso3d.py): stretching depth by sqrt(D_h/D_z) makes the patch model equal-mixing
in a box of depth H' = H*sqrt(D_h/D_z). If H' > L (a "column"), the patchwork freezes into stacked layers of
opposite hands instead of healing (MEASURED in aniso3d.py). Two cases:
  flat (pessimistic)  layers never merge -> physics can't win ocean-wide.
  curved (real moon)  layer walls are spheres of radius R = L/pi, so they sink at 2*D_z/R, whatever D_h is,
                      and the TOP layer takes over within t = H*R/(2*D_z). (corrected 2026-09-26, sphere_layers.py:
                      walls sink at 2*D_z/R, not 2*D_h/R; the old D_h rule came from a flat-box argument in aniso3d.py
                      that does not apply when the vertical direction follows the radius.) The top layer (thickness ~L
                      in stretched units, ASSUMED from the 3-5 layers seen in aniso3d.py; sphere_layers.py top-layer
                      share 0.52/0.07 vs L/H' 0.25 from 2 runs - no support either way) decides alone, so only a
                      fraction L/H' of the molecules count.
Then rank inputs by how much P(win) changes between the low and high half of each range.
"""
import numpy as np
from ocean import delta_at, YR

A = 5.5
RNG = np.random.default_rng(0)
N = 200_000


def logu(lo, hi):
    return 10 ** RNG.uniform(np.log10(lo), np.log10(hi), N)


# horizontal D_h, icy moons: 1e-2..1e2 (ASSUMED bracket; Zhang, Kang & Marshall 2024 diagnose 0.03-2,
#   Zeng & Jansen 2021 use 0.1-0.25 and find ~1000 yr to mix a hemisphere, Ashkenazy & Tziperman 2021 use 30)
# vertical D_z, icy moons: 1e-10..1e-3 (Zeng & Jansen 2021 Sec. II.2 kappa_z, 3e-10..3e-3, rounded outward)
# Earth: horizontal eddies 1e2..1e4 (Abernathey & Marshall 2013); vertical 1e-6..1e-4 (ASSUMED around Ledwell 1.1e-5)
WORLDS = {  # volume m^3, size L m, conc range M, D_h range m^2/s, D_z range m^2/s, time range yr
    "Enceladus": (2.7e16, np.pi * 232e3, (1e-12, 1e-4), (1e-2, 1e2), (1e-10, 1e-3), (1e6, 1e9)),
    "Europa": (3.3e18, np.pi * 1460e3, (1e-12, 1e-4), (1e-2, 1e2), (1e-10, 1e-3), (1e8, 4e9)),
    "Early Earth ocean": (1.332e18, np.pi * 6.371e6, (1e-7, 3e-4), (1e2, 1e4), (1e-6, 1e-4), (1e7, 5e8)),
}


def depth(V, L):
    """Mean ocean depth: volume over the area of a sphere whose half circumference is L."""
    return V / (4 * np.pi * (L / np.pi) ** 2)


def sweep(V, L, c_rng, dh_rng, dz_rng, t_rng, k2_rng=(1e-6, 1.0)):
    x = {"bias g": logu(4e-19, 4e-16), "rate k2": logu(*k2_rng), "concentration": logu(*c_rng),
         "horizontal D_h": logu(*dh_rng), "vertical D_z": logu(*dz_rng), "time": logu(*t_rng) * YR, "volume": V * RNG.uniform(0.7, 1.3, N),
         "chem. factor F": RNG.uniform(0.3, 0.7, N)}
    finishes = x["rate k2"] * x["concentration"] * x["time"] >= 10
    H = depth(x["volume"], L)  # a bigger volume means a deeper ocean, so depth varies with it
    Dh, Dz, t = x["horizontal D_h"], x["vertical D_z"], x["time"]
    across = (L / A) ** 2 / Dh <= t
    column = H * np.sqrt(Dh / Dz) > L                   # stretched depth deeper than wide -> stacked layers
    # layers first need the patchwork to span sideways (across), then the inner layers sink away
    creep = across & (H * (L / np.pi) / (2 * Dz) <= t)  # walls sink at 2*D_z/R (sphere_layers.py)
    top = np.minimum(1.0, L / (H * np.sqrt(Dh / Dz)))  # ASSUMED share of molecules in the deciding (top) layer
    delta = lambda f: x["chem. factor F"] * delta_at(x["concentration"], x["bias g"], x["rate k2"], x["volume"], t, f)
    strong = delta(1.0) >= 2
    win_flat = finishes & across & ~column & strong
    win = finishes & np.where(column, creep & (delta(top) >= 2), across & strong)
    return x, win, {"creep in time": creep.mean(), "finishes": finishes.mean(), "heals across": across.mean(), "layers (column)": column.mean(),
                    "strong enough": strong.mean(), "flat-layer bound": win_flat.mean()}


# README headline numbers (Summary + Part 8 table), checked below against this script's own output.
# (Corrections 2026-09-25: first 18% / 18% with one mixing D; then 43% / 58% before layering was modelled;
#  2026-09-26: 55% / 78% / 100% while layers merged at the sideways rate 2*D_h/R.)
README_WIN_PCT = {"Enceladus": 0.22, "Europa": 0.29, "Early Earth ocean": 0.75}
README_FLAT_PCT = {"Enceladus": 0.03, "Europa": 0.07, "Early Earth ocean": 0.40}  # if layers never merged
WIN_TOL = 0.03  # percentage points; N=200,000 with a fixed seed makes this tight
README_TOP_SPLIT = {"Enceladus": (0.00, 0.43), "Europa": (0.00, 0.57)}  # Part 8: vertical D_z low / high half
README_CREEP_PCT = {"Enceladus": 0.33, "Europa": 0.34}  # share where layers merge in time
README_COLUMN_PCT = {"Enceladus": 0.96, "Europa": 0.92}  # Part 3b: share of input space where the ocean layers


def main():
    for name, args in WORLDS.items():
        x, win, parts = sweep(*args)
        print(f"\n{name}: physics wins in {win.mean():.0%} of plausible input space")
        print("  each condition alone: " + ", ".join(f"{k} {v:.0%}" for k, v in parts.items()))
        effects = []
        for k, v in x.items():
            lo = v <= np.median(v)
            effects.append((win[~lo].mean() - win[lo].mean(), k, win[lo].mean(), win[~lo].mean()))
        effects.sort(key=lambda e: -abs(e[0]))
        for d, k, plo, phi in effects:
            print(f"  {k:15} low half {plo:5.0%}  high half {phi:5.0%}   (effect {d:+.0%})")
        # the README states this world's headline win-rate explicitly (Summary / Part 8 table) -
        # catch it silently drifting if a range or condition changes.
        assert abs(win.mean() - README_WIN_PCT[name]) < WIN_TOL, (
            f"{name}: win rate {win.mean():.0%} no longer matches README's {README_WIN_PCT[name]:.0%}")
        assert abs(parts["flat-layer bound"] - README_FLAT_PCT[name]) < WIN_TOL, f"{name}: flat bound drifted"
        if name in ("Enceladus", "Europa"):
            # README: vertical mixing decides most (layers must merge in time); sideways healing always happens.
            assert effects[0][1] == "vertical D_z", (
                f"{name}: most decisive input is {effects[0][1]!r}, not 'vertical D_z' as the README claims")
            lo, hi = README_TOP_SPLIT[name]
            assert abs(effects[0][2] - lo) < WIN_TOL and abs(effects[0][3] - hi) < WIN_TOL, (
                f"{name}: D_z split {effects[0][2]:.0%}/{effects[0][3]:.0%} no longer matches README")
            assert parts["heals across"] == 1.0, f"{name}: healing now binds"
            assert abs(parts["creep in time"] - README_CREEP_PCT[name]) < WIN_TOL, f"{name}: merge-in-time share drifted"
            assert abs(parts["layers (column)"] - README_COLUMN_PCT[name]) < WIN_TOL, f"{name}: layered share drifted"
            g = [e[0] for e in effects if e[1] == "bias g"][0]  # Part 8: "bias g (+4 to +8 points)"
            assert 0.04 - WIN_TOL < g < 0.08 + WIN_TOL, f"{name}: bias g effect {g:+.0%} outside README's +4..+8"
    # README Part 3 note: when the ocean layers, and how long the layers take to merge
    for name, ratio, t_worst in (("Enceladus", 334, 1.5e12), ("Europa", 1386, 2.8e13)):
        V, L, _, _, (dz_lo, _), _ = WORLDS[name]
        H = depth(V, L)
        t_merge = H * (L / np.pi) / (2 * dz_lo) / YR
        print(f"{name}: layers form if D_z < D_h/{(L / H) ** 2:.0f}; they merge within {t_merge:.1e} yr at D_z = {dz_lo}")
        assert abs((L / H) ** 2 / ratio - 1) < 0.05 and abs(t_merge / t_worst - 1) < 0.1
    # sanity: a best-case input set must win, a hopeless one must lose
    assert delta_at(1e-6, 1e-16, 1e-3, 1.3e18, 1e8 * YR) * 0.7 > 2
    assert delta_at(1e-6, 1e-18, 1e-3, 3.0, 1e2 * YR) * 0.7 < 2


if __name__ == "__main__":
    main()
