"""3D per-patch selection constant K3 (the 3D version of domains.py check C).

Stop the sweep just after the transition (lam = 0.3), before patches coarsen much,
measure patch size l and the fraction on the favoured hand, and solve
    Phi(PREF * g * gamma**-0.25 * sqrt(K3 * l**3 / eps)) = fraction
for K3. Consistency across g values is the check.
"""
import sys

import numpy as np
from scipy.stats import norm
from domains3d import sweep3d, domain_length

PREF = np.sqrt(2) * np.pi**0.25


def main(seed_l=5, seed_g=1):
    D, gamma, eps, lam_end = 1.0, 0.1, 1e-4, 0.3
    l = domain_length(sweep3d(D, gamma, eps=eps, lam_end=lam_end, runs=2, seed=seed_l))
    print(f"patch size at lam={lam_end}: {l:.1f} cells", flush=True)
    Ks = []
    for g in (5e-5, 1e-4, 2e-4):
        p = (sweep3d(D, gamma, g=g, eps=eps, lam_end=lam_end, runs=3, seed=seed_g) > 0).mean()
        K = (norm.ppf(p) / (PREF * g * gamma**-0.25)) ** 2 * eps / l**3
        Ks.append(K)
        print(f"  g={g:.0e}: favoured {p:.3f} -> K3={K:.3f}", flush=True)
    K3 = float(np.mean(Ks))
    print(f"K3 = {K3:.3f} (spread {np.std(Ks) / K3:.0%})")
    assert np.std(Ks) / K3 < 0.3, "3D patch selection not described by N_eff = K3*l^3"
    return l, Ks, K3


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:3]]
    main(*args)
