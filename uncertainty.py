"""Uncertainty sweep: vary every uncertain input at once; which conclusions survive, and which
input decides the answer?

For each world, draw 200,000 random input sets (log-uniform over each plausible range) and
count "physics wins" only if ALL hold:
  finishes      k2*c*tau >= 10                                   (ocean.py time rule)
  one ocean     patchwork heals in time: (L/A)^2/D <= tau        (coarsen3d.py, A = 5.5)
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


WORLDS = {  # volume m^3, size L m, conc range M, mixing range m^2/s, time range yr
    "Enceladus": (2.7e16, np.pi * 232e3, (1e-12, 1e-4), (1e-10, 1e-3), (1e6, 1e9)),
    "Europa": (3.3e18, np.pi * 1460e3, (1e-12, 1e-4), (1e-10, 1e-3), (1e8, 4e9)),
    "Early Earth ocean": (1.332e18, np.pi * 6.371e6, (1e-7, 3e-4), (1e2, 1e4), (1e7, 5e8)),
}


def sweep(V, L, c_rng, d_rng, t_rng):
    x = {"bias g": logu(4e-19, 4e-16), "rate k2": logu(1e-6, 1.0), "concentration": logu(*c_rng),
         "mixing D": logu(*d_rng), "time": logu(*t_rng) * YR, "volume": V * RNG.uniform(0.7, 1.3, N),
         "chem. factor F": RNG.uniform(0.3, 0.7, N)}
    finishes = x["rate k2"] * x["concentration"] * x["time"] >= 10
    heals = (L / A) ** 2 / x["mixing D"] <= x["time"]
    strong = x["chem. factor F"] * delta_at(x["concentration"], x["bias g"], x["rate k2"], x["volume"], x["time"]) >= 2
    win = finishes & heals & strong
    return x, win, {"finishes": finishes.mean(), "heals": heals.mean(), "strong enough": strong.mean()}


# README headline numbers (Summary + Part 8 table), checked below against this script's own output.
README_WIN_PCT = {"Enceladus": 0.18, "Europa": 0.18, "Early Earth ocean": 1.00}
WIN_TOL = 0.03  # percentage points; N=200,000 with a fixed seed makes this tight
README_MIXING_HIGH = {"Enceladus": 0.36, "Europa": 0.37}  # Part 8: mixing D low half 0%, high half 36%/37%


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
            # README: "mixing the deciding input" / "most decisive input: mixing D" for both icy moons.
            assert effects[0][1] == "mixing D", (
                f"{name}: most decisive input is {effects[0][1]!r}, not 'mixing D' as the README claims")
            assert effects[0][2] < 0.01 and abs(effects[0][3] - README_MIXING_HIGH[name]) < WIN_TOL, (
                f"{name}: mixing split {effects[0][2]:.0%}/{effects[0][3]:.0%} no longer matches README")
    # sanity: a best-case input set must win, a hopeless one must lose
    assert delta_at(1e-6, 1e-16, 1e-3, 1.3e18, 1e8 * YR) * 0.7 > 2
    assert delta_at(1e-6, 1e-18, 1e-3, 3.0, 1e2 * YR) * 0.7 < 2


if __name__ == "__main__":
    main()
