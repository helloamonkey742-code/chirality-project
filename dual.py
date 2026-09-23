"""The dual explanation: L-amino acids AND D-sugars.

Two coupled chiral systems passing through the transition together:
    alpha: amino acids   (+ = L)
    beta:  sugars        (+ = D)
    d(alpha) = (lam*alpha - alpha**3 + kappa*beta + g_a) dt + sqrt(eps) dW1
    d(beta)  = (lam*beta  - beta**3  + kappa*alpha + g_s) dt + sqrt(eps) dW2
kappa > 0: each hand catalyses the partner hand (L-aa -> D-sugar, as in Pizzarello &
Weber 2004 / Breslow & Cheng 2010; D-RNA -> L-aa, as in Tamura & Schimmel 2004).
g_a = weak-force bias toward L-aa; g_s = bias toward D-sugar (negative if it favours L-sugar).

Prediction (pair mode u = (alpha+beta)/sqrt2 grows first, at lam = -kappa):
  1. matched pair (L-aa with D-sugar, or D-aa with L-sugar) wins once kappa >> sqrt(gamma)
  2. P(L-aa & D-sugar | matched) = Phi(Delta) with bias (g_a + g_s)/sqrt2
     -> the two weak-force biases ADD if they favour L-aa and D-sugar, CANCEL if not.
"""
import numpy as np
from scipy.stats import norm
from sim import delta


def simulate(kappa, g_a, g_s, eps=1e-4, gamma=0.01, runs=4000, dt=0.01, seed=0):
    rng = np.random.default_rng(seed)
    a, b = np.zeros(runs), np.zeros(runs)
    for t in np.arange(-(1 + kappa) / gamma, 1 / gamma, dt):
        lam = gamma * t
        da = (lam * a - a**3 + kappa * b + g_a) * dt
        db = (lam * b - b**3 + kappa * a + g_s) * dt
        a += da + np.sqrt(eps * dt) * rng.standard_normal(runs)
        b += db + np.sqrt(eps * dt) * rng.standard_normal(runs)
    matched = np.sign(a) == np.sign(b)
    natural = (a > 0) & (b > 0)
    return matched.mean(), natural.sum() / max(matched.sum(), 1)


def main():
    eps, gamma = 1e-4, 0.01
    print("1. how often the pair is matched (no bias)")
    for kappa in [0.0, 0.02, 0.05, 0.1, 0.3]:
        m, _ = simulate(kappa, 0.0, 0.0, runs=2000)
        print(f"   kappa={kappa:<5} matched {m:.3f}")
    m0, _ = simulate(0.0, 0.0, 0.0, runs=2000)
    m1, _ = simulate(0.3, 0.0, 0.0, runs=2000)
    assert abs(m0 - 0.5) < 0.05 and m1 > 0.99

    print("\n2. which pair wins (kappa=0.3): biases add or cancel")
    print(f"   {'g_a':>7} {'g_s':>7} {'predicted':>9} {'simulated':>9}")
    for g_a, g_s in [(1e-3, 0.0), (1e-3, 1e-3), (1e-3, -1e-3), (2e-3, -1e-3)]:
        pred = norm.cdf(delta((g_a + g_s) / np.sqrt(2), eps, gamma))
        m, p = simulate(0.3, g_a, g_s)
        print(f"   {g_a:7.0e} {g_s:7.0e} {pred:9.3f} {p:9.3f}   (matched {m:.3f})")
        assert abs(p - pred) < 0.03, "pair-mode prediction fails"


if __name__ == "__main__":
    main()
