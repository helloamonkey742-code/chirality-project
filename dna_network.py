"""Asymmetric DNA Sugar - Amino Acid Coupling and Parity Violation (dna_network.py)

Extends dual.py to model asymmetric cross-catalysis between nucleic acids (DNA/RNA)
and peptides, addressing the unsettled PVED sign of sugars (Faglioni et al. 2005/2006):

    d(alpha) = (lam*alpha - alpha**3 + kappa_sa*beta + g_a) dt + sqrt(eps) dW1   [amino acids, + = L]
    d(beta)  = (lam*beta  - beta**3  + kappa_as*alpha + g_s) dt + sqrt(eps) dW2   [sugars/DNA,  + = D]

kappa_sa: RNA/DNA template -> peptide synthesis stereocontrol (Tamura & Schimmel 2004; Ozturk & Sasselov 2025)
kappa_as: Amino acid / peptide -> sugar / nucleoside catalysis (Breslow & Cheng 2010; Hein et al. 2011)

Linearized mode analysis:
  Growing pair eigenvalue: lam_+ = lam + sqrt(kappa_sa * kappa_as) = lam + kappa_bar
  Left eigenvector: w = (sqrt(kappa_as), sqrt(kappa_sa)) / sqrt(kappa_as + kappa_sa)
  Effective projected bias:
      g_eff = (sqrt(kappa_as) * g_a + sqrt(kappa_sa) * g_s) / sqrt(kappa_as + kappa_sa)
  Condition to rescue natural (L-aa, D-sugar) homochirality when g_s < 0 (DNA PVED opposing):
      kappa_as / kappa_sa > ( |g_s| / g_a )^2

Within this model, the case where DNA double-helix PVED opposes the natural form is handled:
as long as amino-acid-to-sugar transmission exceeds the squared bias ratio, L-aa drags DNA sugars into D!
"""
import numpy as np
from scipy.stats import norm
from sim import delta


def simulate_asymmetric(kappa_sa, kappa_as, g_a, g_s, eps=1e-4, gamma=0.01, runs=4000, dt=0.01, seed=0):
    rng = np.random.default_rng(seed)
    kappa_bar = np.sqrt(kappa_sa * kappa_as)
    a, b = np.zeros(runs), np.zeros(runs)
    # Start sweep before the earliest bifurcating mode: lam = -kappa_bar
    t_start = -(1 + kappa_bar) / gamma
    t_end = 1.0 / gamma
    for t in np.arange(t_start, t_end, dt):
        lam = gamma * t
        da = (lam * a - a**3 + kappa_sa * b + g_a) * dt
        db = (lam * b - b**3 + kappa_as * a + g_s) * dt
        a += da + np.sqrt(eps * dt) * rng.standard_normal(runs)
        b += db + np.sqrt(eps * dt) * rng.standard_normal(runs)
    matched = (np.sign(a) == np.sign(b))
    natural = (a > 0) & (b > 0)
    return matched.mean(), natural.sum() / max(matched.sum(), 1)


def theoretical_p(kappa_sa, kappa_as, g_a, g_s, eps=1e-4, gamma=0.01):
    w_a = np.sqrt(kappa_as)
    w_s = np.sqrt(kappa_sa)
    norm_w = np.sqrt(kappa_as + kappa_sa)
    g_eff = (w_a * g_a + w_s * g_s) / norm_w
    d = delta(g_eff, eps, gamma)
    return float(norm.cdf(d))


def main():
    eps, gamma = 1e-4, 0.01
    print("=== Asymmetric DNA-Peptide Coupling Verification ===")
    
    # 1. Verification of symmetric baseline (matches dual.py)
    print("\n1. Symmetric baseline check (kappa_sa = kappa_as = 0.3):")
    p_sim_sym = []
    p_thy_sym = []
    for g_a, g_s in [(1e-3, 0.0), (1e-3, 1e-3), (1e-3, -1e-3), (2e-3, -1e-3)]:
        thy = theoretical_p(0.3, 0.3, g_a, g_s, eps, gamma)
        matched, sim = simulate_asymmetric(0.3, 0.3, g_a, g_s, eps, gamma, runs=3000)
        print(f"   g_a={g_a:7.0e}  g_s={g_s:7.0e}  predicted={thy:.3f}  simulated={sim:.3f}  (matched {matched:.3f})")
        assert abs(sim - thy) < 0.035, f"Symmetric prediction failed: {sim} vs {thy}"
    
    # 2. Asymmetric cross-catalysis: aa -> sugar dominance (kappa_as >> kappa_sa)
    # Testing rescue when sugar PVED opposes natural D-sugar (g_s = -1e-3, g_a = +1e-3)
    print("\n2. Asymmetric coupling: Can amino acid catalysis rescue opposing DNA sugar PVED?")
    print("   Setting g_a = +1e-3 (favours L-aa), g_s = -1e-3 (favours L-sugar, opposing natural D)")
    g_a, g_s = 1e-3, -1e-3
    for kappa_as, kappa_sa in [(0.1, 0.9), (0.3, 0.3), (0.9, 0.1), (1.6, 0.1)]:
        ratio = kappa_as / kappa_sa
        thy = theoretical_p(kappa_sa, kappa_as, g_a, g_s, eps, gamma)
        matched, sim = simulate_asymmetric(kappa_sa, kappa_as, g_a, g_s, eps, gamma, runs=3000)
        print(f"   kappa_as={kappa_as:.2f}  kappa_sa={kappa_sa:.2f} (ratio {ratio:4.1f}) -> predicted={thy:.3f}  simulated={sim:.3f}")
        assert abs(sim - thy) < 0.035, f"Asymmetric prediction failed: {sim} vs {thy}"
    
    # 3. Critical threshold check:
    # When kappa_as / kappa_sa > ( |g_s| / g_a )^2 = 1, P > 0.50 (natural wins).
    # When kappa_as / kappa_sa < 1, P < 0.50 (mirror wins).
    thy_sub = theoretical_p(0.9, 0.1, g_a, g_s, eps, gamma)  # ratio = 0.1 / 0.9 = 0.11 < 1 -> mirror wins
    thy_sup = theoretical_p(0.1, 0.9, g_a, g_s, eps, gamma)  # ratio = 0.9 / 0.1 = 9 > 1 -> natural wins
    assert thy_sub < 0.40, "Expected mirror pair to dominate when template transmission is reversed"
    assert thy_sup > 0.60, "Expected natural pair to dominate when aa catalysis rescues"

    print("\nInternal self-consistency check passed: linear theory agrees with this model's own stochastic runs.")


if __name__ == "__main__":
    main()
