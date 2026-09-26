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


def sweep(V, L, c_rng, dh_rng, dz_rng, t_rng):
    x = {"bias g": logu(4e-19, 4e-16), "rate k2": logu(1e-6, 1.0), "concentration": logu(*c_rng),
         "horizontal D_h": logu(*dh_rng), "vertical D_z": logu(*dz_rng), "time": logu(*t_rng) * YR, "volume": V * RNG.uniform(0.7, 1.3, N),
         "chem. factor F": RNG.uniform(0.3, 0.7, N)}
    finishes = x["rate k2"] * x["concentration"] * x["time"] >= 10
    H = depth(x["volume"], L)  # a bigger volume means a deeper ocean, so depth varies with it
    heals = ((L / A) ** 2 / x["horizontal D_h"] <= x["time"]) & ((H / A) ** 2 / x["vertical D_z"] <= x["time"])
    strong = x["chem. factor F"] * delta_at(x["concentration"], x["bias g"], x["rate k2"], x["volume"], x["time"]) >= 2
    win = finishes & heals & strong
    across = (L / A) ** 2 / x["horizontal D_h"] <= x["time"]
    return x, win, {"finishes": finishes.mean(), "heals": heals.mean(), "heals across": across.mean(),
                    "strong enough": strong.mean()}


# README headline numbers (Summary + Part 8 table), checked below against this script's own output.
# (Correction 2026-09-25: were 18% / 18% with one mixing D; see the WORLDS comment.)
README_WIN_PCT = {"Enceladus": 0.43, "Europa": 0.58, "Early Earth ocean": 1.00}
WIN_TOL = 0.03  # percentage points; N=200,000 with a fixed seed makes this tight
README_DZ_SPLIT = {"Enceladus": (0.18, 0.68), "Europa": (0.28, 0.87)}  # Part 8: D_z low half / high half


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
        if name in ("Enceladus", "Europa"):
            # README: "vertical mixing is the deciding input"; horizontal mixing never limits healing.
            assert effects[0][1] == "vertical D_z", (
                f"{name}: most decisive input is {effects[0][1]!r}, not 'vertical D_z' as the README claims")
            lo, hi = README_DZ_SPLIT[name]
            assert abs(effects[0][2] - lo) < WIN_TOL and abs(effects[0][3] - hi) < WIN_TOL, (
                f"{name}: D_z split {effects[0][2]:.0%}/{effects[0][3]:.0%} no longer matches README")
            assert parts["heals across"] == 1.0, f"{name}: horizontal healing now binds; README says it never does"
    # README Part 3 caveat: mixing thresholds for healing within the oldest / youngest plausible age
    for name, lo_across, lo_down, hi_down, ratio in (("Enceladus", 5.6e-4, 1.7e-9, 1.7e-6, 334),
                                                      ("Europa", 2.2e-4, 4.0e-9, 1.6e-7, 1386)):
        V, L, _, _, _, (t0, t1) = WORLDS[name]
        H = depth(V, L)
        need = lambda size, t: (size / A) ** 2 / (t * YR)
        print(f"{name}: heals across if D_h >= {need(L, t0):.1e}; down if D_z >= {need(H, t1):.1e}..{need(H, t0):.1e}"
              f" m^2/s; depth {H / 1e3:.0f} km; vertical threshold {(L / H) ** 2:.0f}x below the old one-D rule")
        for got, want in ((need(L, t0), lo_across), (need(H, t1), lo_down), (need(H, t0), hi_down), ((L / H) ** 2, ratio)):
            assert abs(got / want - 1) < 0.05, f"{name}: threshold {got:.2e} drifted from README's {want:.2e}"
    # sanity: a best-case input set must win, a hopeless one must lose
    assert delta_at(1e-6, 1e-16, 1e-3, 1.3e18, 1e8 * YR) * 0.7 > 2
    assert delta_at(1e-6, 1e-18, 1e-3, 3.0, 1e2 * YR) * 0.7 < 2


if __name__ == "__main__":
    main()
