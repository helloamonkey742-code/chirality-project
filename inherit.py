"""Can a small pool inherit the ocean's bias?

A pool fills with water carrying a small enantiomeric excess e0 (from the ocean, a
meteorite, ...), then concentrates (evaporation, thermophoresis) until its own
autocatalysis switches on. Model: start AT the transition (lam = 0) with alpha = e0.

Linear theory gives
    Delta = (e0 + g*sqrt(pi/(2*gamma))) / sqrt(eps * sqrt(pi/gamma) / 2)
so with eps = 1/N the pool reliably follows the inherited hand (Delta >= 2) once
    e0 >= e_need = 2 * (pi*k*tau)**0.25 / sqrt(2N)
Then compare e_need with what each source can supply.
"""
import numpy as np
from scipy.stats import norm
from ocean import NA, YR


def simulate(e0, g, eps, gamma, runs=4000, dt=0.01, seed=0):
    rng = np.random.default_rng(seed)
    a = np.full(runs, e0)
    for t in np.arange(0, 1 / gamma, dt):
        a += (gamma * t * a - a**3 + g) * dt + np.sqrt(eps * dt) * rng.standard_normal(runs)
    return (a > 0).mean()


def delta(e0, g, eps, gamma):
    return (e0 + g * np.sqrt(np.pi / (2 * gamma))) / np.sqrt(eps * np.sqrt(np.pi / gamma) / 2)


def e_need(c, V, k2, tau):
    N = c * V * 1000 * NA
    return 2 * (np.pi * k2 * c * tau) ** 0.25 / np.sqrt(2 * N)


POOLS = {  # volume m^3, concentration M, time s  (as in scenarios.py)
    "Warm little pond (1 m × 1 m)": (np.pi, 1e-2, 1e2 * YR),
    "Vent pore (1 mm³)": (1e-9, 1e-1, 1e5 * YR),
    "Lake / lagoon (1 km² × 10 m)": (1e7, 1e-4, 1e4 * YR),
}
SOURCES = {
    "weak force alone (equilibrium, g/2)": 5e-18,
    "meteorite isovaline (Murchison, Glavin & Dworkin 2009)": 0.185,
}


def main():
    eps, gamma = 1e-4, 0.01
    print("check: inherited excess e0 (no bias) -> P(pool follows it)")
    for e0 in [0.0, 0.01, 0.03, 0.06]:
        pred, sim = norm.cdf(delta(e0, 0, eps, gamma)), simulate(e0, 0, eps, gamma)
        print(f"  e0={e0:<6} predicted {pred:.3f}  simulated {sim:.3f}")
        assert abs(pred - sim) < 0.03
    pred, sim = norm.cdf(delta(0.02, 5e-4, eps, gamma)), simulate(0.02, 5e-4, eps, gamma)
    print(f"  e0=0.02 + bias g=5e-4: predicted {pred:.3f}  simulated {sim:.3f}")
    assert abs(pred - sim) < 0.03

    print("\ninherited excess each pool needs to follow reliably (k2 = 1e-3 /M/s):")
    for name, (V, c, tau) in POOLS.items():
        print(f"  {name:30} e0 >= {e_need(c, V, 1e-3, tau):.0e}")
    print("\nwhat sources can supply:")
    for name, e in SOURCES.items():
        print(f"  {name:55} {e:.0e}")


if __name__ == "__main__":
    main()
