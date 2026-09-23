"""Bias vs noise at a chiral symmetry-breaking transition.

Normal form of a Frank/Kondepudi-Nelson autocatalytic system passing slowly
through its bifurcation:

    d(alpha) = (lam(t)*alpha - alpha**3 + g) dt + sqrt(eps) dW,   lam = gamma*t

alpha > 0 means the favoured hand won (L if the parity-violating energy
difference favours L, D if it favours D -- the sign is unknown, see README).

  g      chiral bias      ~ dE_pv / kT            (weak force)
  eps    chiral noise     ~ 1 / N                 (molecules in well-mixed volume)
  gamma  sweep rate       ~ 1 / (k * tau)         (reaction rate x time to cross)

Linearising around lam = 0 gives the selection probability

    P(favoured) = Phi(Delta),  Delta = sqrt(2) * pi**0.25 * g * gamma**-0.25 * eps**-0.5

This script checks that formula against brute-force stochastic runs, then tests
whether temperature noise (random jitter in lam) breaks it.
"""
import numpy as np
from scipy.stats import norm


def delta(g, eps, gamma):
    return np.sqrt(2) * np.pi**0.25 * g * gamma**-0.25 * eps**-0.5


def simulate(g, eps, gamma, runs=4000, dt=0.01, sigma_T=0.0, seed=0):
    """Fraction of runs ending on the favoured hand. sigma_T = temperature noise on lam."""
    rng = np.random.default_rng(seed)
    a = np.zeros(runs)
    for t in np.arange(-1 / gamma, 1 / gamma, dt):
        lam = gamma * t + sigma_T * rng.standard_normal(runs) / np.sqrt(dt)
        a += (lam * a - a**3 + g) * dt + np.sqrt(eps * dt) * rng.standard_normal(runs)
    return (a > 0).mean()


def main():
    eps, gamma = 1e-4, 0.01
    print("check: analytic Phi(Delta) vs simulated P(favoured hand)")
    print(f"{'g':>9} {'Delta':>6} {'analytic':>9} {'simulated':>9}")
    rows = []
    for g in [0, 5e-4, 1e-3, 2e-3, 4e-3]:
        d = delta(g, eps, gamma)
        p_sim = simulate(g, eps, gamma)
        rows.append((norm.cdf(d), p_sim))
        print(f"{g:9.1e} {d:6.2f} {norm.cdf(d):9.3f} {p_sim:9.3f}")
    # 4000 runs -> std error <= 0.008; allow 3 sigma
    assert all(abs(a - s) < 0.03 for a, s in rows), "formula disagrees with simulation"

    print("\ntemperature noise (jitter in lam), g=1e-3:")
    base = simulate(1e-3, eps, gamma)
    for s in [0.0, 0.05, 0.2]:
        print(f"  sigma_T={s:<5} P={simulate(1e-3, eps, gamma, sigma_T=s):.3f}")
    print(f"  (no-temperature-noise baseline {base:.3f})")
    return rows


if __name__ == "__main__":
    main()
