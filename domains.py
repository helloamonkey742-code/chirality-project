"""One ocean or a patchwork? Spatial version of sim.py.

A 1D ring of cells, each running the same normal form, coupled by mixing D:

    d(alpha_i) = (lam*alpha_i - alpha_i**3 + g + D*laplacian(alpha)_i) dt + sqrt(eps) dW

Three checks:
  A. domain size at the transition   l = C * sqrt(D) * gamma**-0.25   (Kibble-Zurek-type)
  B. bias-driven wall speed          v = (3/sqrt2) * g * sqrt(D)       (favoured hand invades)
  C. local selection inside a patch  P = Phi(Delta) with N_eff = K * l cells

C decides whether an ocean that breaks into patches still ends up mostly one hand.
"""
import numpy as np
from scipy.stats import norm

PREF = np.sqrt(2) * np.pi**0.25


def lap(a):
    return np.roll(a, 1, -1) + np.roll(a, -1, -1) - 2 * a


def sweep(D, gamma, g=0.0, eps=1e-4, M=2000, runs=16, seed=0):
    """Sweep lam from -1 to +1. Returns final field, shape (runs, M)."""
    rng = np.random.default_rng(seed)
    dt = min(0.05, 0.2 / D)
    a = np.zeros((runs, M))
    for t in np.arange(-1 / gamma, 1 / gamma, dt):
        a += (gamma * t * a - a**3 + g + D * lap(a)) * dt + np.sqrt(eps * dt) * rng.standard_normal(a.shape)
    return a


def domain_length(a):
    walls = (np.sign(a) != np.sign(np.roll(a, 1, -1))).sum()
    return a.size / max(walls, 1)


def wall_speed(g=0.02, D=4.0, M=800, T=400.0):
    dt = min(0.05, 0.2 / D)
    x = np.arange(M)
    a = np.where((x > M // 4) & (x < 3 * M // 4), -1.0, 1.0)   # unfavoured block in the middle
    width0 = (a < 0).sum()
    for _ in range(int(T / dt)):
        a += (a - a**3 + g + D * lap(a)) * dt
    return (width0 - (a < 0).sum()) / 2 / T          # each of the two walls moves inward


def main():
    print("A. domain size at freeze-out, no bias")
    ratios = []
    for gamma in [0.01, 0.0025]:
        for D in [1.0, 4.0, 16.0]:
            l = domain_length(sweep(D, gamma))
            r = l / (np.sqrt(D) * gamma**-0.25)
            ratios.append(r)
            print(f"  gamma={gamma:<7} D={D:<5} domain={l:7.1f} cells   C={r:.2f}")
    C = float(np.mean(ratios))
    print(f"  C = {C:.2f} (spread {np.std(ratios) / C:.0%})")
    assert np.std(ratios) / C < 0.25, "domain size does not follow sqrt(D)*gamma^-1/4"

    print("\nB. wall speed with bias")
    for g, D in [(0.02, 4.0), (0.01, 16.0)]:
        v, v_pred = wall_speed(g, D), 3 / np.sqrt(2) * g * np.sqrt(D)
        print(f"  g={g} D={D}: simulated {v:.4f}  predicted {v_pred:.4f}")
        assert abs(v / v_pred - 1) < 0.1

    print("\nC. local selection inside patches (D=4, gamma=0.01, eps=1e-4)")
    D, gamma, eps = 4.0, 0.01, 1e-4
    l = C * np.sqrt(D) * gamma**-0.25
    Ks = []
    for g in [1e-4, 2e-4, 3e-4]:
        p = (sweep(D, gamma, g=g, eps=eps) > 0).mean()
        # solve Phi(PREF*g*gamma^-1/4*sqrt(K*l/eps)) = p for K
        K = (norm.ppf(p) / (PREF * g * gamma**-0.25)) ** 2 * eps / l
        Ks.append(K)
        print(f"  g={g}: fraction favoured {p:.3f}  -> K={K:.2f}")
    K = float(np.mean(Ks))
    print(f"  K = {K:.2f} (spread {np.std(Ks) / K:.0%})")
    assert np.std(Ks) / K < 0.3, "patch selection not described by N_eff = K*l"
    return C, K


if __name__ == "__main__":
    main()
