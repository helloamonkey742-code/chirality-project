"""Threat test: does a wobbling environment break the slow-sweep formula?

ocean.py uses tau = the ocean's age as the time to cross the transition, so Delta ~ (k*tau)^(1/4).
That assumes the control parameter lam rises smoothly and only once. A real ocean's temperature,
heat flow and feed wobble (orbital and heat-flux cycles), so lam = gamma*t + w*sin(omega*t).
If the wobble swings lam back below zero, the excess decays and the hand is re-chosen on each
cycle, and only the last crossings count.

Same normal form as sim.py:  d(alpha) = (lam*alpha - alpha^3 + g) dt + sqrt(eps) dW.
Reports Delta_eff = Phi^-1(P favoured) against the no-wobble Delta.
"""
import sys
import numpy as np
from scipy.stats import norm
from scipy.special import i0e
from sim import delta

EPS, GAMMA = 1e-4, 0.01
G = 1.5 / delta(1.0, EPS, GAMMA)  # no-wobble Delta = 1.5 (P = 0.933)
RUNS = 4000


def simulate(w, omega, runs=RUNS, dt=0.01, seed=0):
    rng = np.random.default_rng(seed)
    a = np.zeros(runs)
    t_end = (w + 1) / GAMMA  # lam >= 1 from here on: locked in
    sq = np.sqrt(EPS * dt)
    for t in np.arange(-t_end, t_end, dt):
        lam = GAMMA * t + w * np.sin(omega * t)
        a += (lam * a - a**3 + G) * dt + sq * rng.standard_normal(runs)
    return (a > 0).mean()


def bessel(w, omega):
    """Linear theory, many wobble periods per decision window: Delta * I0(x)/sqrt(I0(2x)), x = w/omega."""
    x = w / omega
    return 1.5 * i0e(x) / np.sqrt(i0e(2 * x))


def main():
    base = simulate(0, 1)
    print(f"no wobble: P = {base:.3f}  (formula {norm.cdf(1.5):.3f})")
    assert abs(base - norm.cdf(1.5)) < 0.03
    print(f"{'omega':>6} {'w':>5} {'w*omega/gamma':>13} {'P':>6} {'Delta_eff':>9} {'ratio':>6} {'bessel':>6}")
    for omega, ws in [(1.0, [1, 3, 10]), (0.2, [0.5, 2, 8]), (0.02, [0.5, 2, 8])]:
        for w in ws:
            p = simulate(w, omega)
            d = norm.ppf(min(max(p, 1e-4), 1 - 1e-4))
            print(f"{omega:6.2f} {w:5.1f} {w * omega / GAMMA:13.1f} {p:6.3f} {d:9.2f} {d / 1.5:6.2f} {bessel(w, omega) / 1.5:6.2f}",
                  flush=True)


if __name__ == "__main__":
    main()
