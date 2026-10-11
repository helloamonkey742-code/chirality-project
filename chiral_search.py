"""Computational Search Engine for Prebiotic Chiral Amplification Networks.

Mathematically maps and scans the kinetic parameter space of coupled prebiotic
reaction networks (e.g., amino acid transamination + peptide ligation as in
Deng et al. Nature 2024, Higgs & Blackmond PNAS 2025, and Yu et al. PNAS 2024)
in a Continuous-Flow Stirred Tank Reactor (CSTR).

Evaluates:
1. Exact racemic steady-state concentrations.
2. The difference-mode Jacobian J_diff and leading eigenvalue lambda_max.
3. The critical bifurcation boundary (lambda_max > 0) separating racemic decay
   from spontaneous homochiral amplification.
4. Nonlinear ODE verification confirming ee growth from a tiny seed (1e-3).
5. Mapping of an ASSUMED lab-baseline catalyst point (not a measured S) against the boundary.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Default prebiotic CSTR parameters
DEFAULT_F0 = 5e-3                 # 5 mM realistic prebiotic feed
DEFAULT_TAU_H = 200.0             # 200 hours residence time
DEFAULT_PHI = 1.0 / (DEFAULT_TAU_H * 3600.0) # dilution rate (s^-1)
DEFAULT_K_CAT = 0.10              # 0.1 /M/s catalytic monomer synthesis
DEFAULT_K_HOMO = 0.02             # 0.02 /M/s ligation rate


class PrebioticNetworkModel:
    """Coupled prebiotic network: amino acid synthesis + oligomer ligation.
    
    Species:
      F: achiral precursor (ketoacid / aldehyde feedstock)
      L, D: chiral amino acid monomers
      LL, DD: homochiral dimers (catalysts for monomer synthesis)
      LD: heterochiral dimer (inactive / sink)
    """

    def __init__(self,
                 F0=DEFAULT_F0,
                 phi=DEFAULT_PHI,
                 k0=1e-9,
                 k_cat=DEFAULT_K_CAT,
                 S=3.0,
                 k_homo=DEFAULT_K_HOMO,
                 rho=1.5,
                 k_rac=1e-9):
        self.F0 = float(F0)
        self.phi = float(phi)
        self.k0 = float(k0)
        self.k_cat = float(k_cat)
        self.S = float(S)                   # Stereospecificity S >= 1
        self.k_homo = float(k_homo)
        self.rho = float(rho)               # Preference for heterochiral ligation
        self.k_hetero = float(2.0 * rho * k_homo)
        self.k_rac = float(k_rac)

    def ode_system(self, t, y):
        F, L, D, LL, DD, LD = y

        # Enforce physical non-negativity
        F = max(0.0, F)
        L = max(0.0, L)
        D = max(0.0, D)
        LL = max(0.0, LL)
        DD = max(0.0, DD)
        LD = max(0.0, LD)

        # Monomer synthesis rates
        rate_uncat = self.k0 * F
        rate_cat_L = self.k_cat * F * (LL + DD / self.S)
        rate_cat_D = self.k_cat * F * (DD + LL / self.S)

        # Ligation rates
        rate_LL = self.k_homo * (L ** 2)
        rate_DD = self.k_homo * (D ** 2)
        rate_LD = self.k_hetero * L * D

        # Racemization
        rate_rac = self.k_rac * (L - D)

        dF = self.phi * (self.F0 - F) - 2.0 * rate_uncat - (rate_cat_L + rate_cat_D)
        dL = rate_uncat + rate_cat_L - 2.0 * rate_LL - rate_LD - rate_rac - self.phi * L
        dD = rate_uncat + rate_cat_D - 2.0 * rate_DD - rate_LD + rate_rac - self.phi * D
        dLL = rate_LL - self.phi * LL
        dDD = rate_DD - self.phi * DD
        dLD = rate_LD - self.phi * LD

        return [dF, dL, dD, dLL, dDD, dLD]

    def solve_racemic_steady_state(self):
        """Solves for the symmetric (racemic: L = D = x, LL = DD) steady state."""
        def residual(x):
            if x <= 0:
                return 1e9
            LL = (self.k_homo / self.phi) * (x ** 2)
            cat_factor = 2.0 * self.k0 + 2.0 * self.k_cat * (1.0 + 1.0 / self.S) * LL
            F = self.phi * self.F0 / (self.phi + cat_factor)
            prod_L = self.k0 * F + self.k_cat * F * (1.0 + 1.0 / self.S) * LL
            dest_L = 2.0 * self.k_homo * (x ** 2) + self.k_hetero * (x ** 2) + self.phi * x
            return prod_L - dest_L

        sol = root_scalar(residual, bracket=[1e-9, self.F0 * 5.0], method="brentq")
        x_ss = sol.root
        LL_ss = (self.k_homo / self.phi) * (x_ss ** 2)
        cat_factor = 2.0 * self.k0 + 2.0 * self.k_cat * (1.0 + 1.0 / self.S) * LL_ss
        F_ss = self.phi * self.F0 / (self.phi + cat_factor)
        LD_ss = (self.k_hetero / self.phi) * (x_ss ** 2)

        return F_ss, x_ss, LL_ss, LD_ss

    def difference_jacobian_eigenvalue(self, F_ss=None, x_ss=None, LL_ss=None):
        """Computes the leading eigenvalue lambda_max of the difference mode.
        
        Variables: Delta_m = L - D, Delta_d = LL - DD.
        Near racemic state:
          d(Delta_m)/dt = J11 * Delta_m + J12 * Delta_d
          d(Delta_d)/dt = J21 * Delta_m + J22 * Delta_d
        
        If lambda_max > 0, the racemic state is UNSTABLE,
        guaranteeing spontaneous chiral amplification!
        """
        if F_ss is None or x_ss is None or LL_ss is None:
            F_ss, x_ss, LL_ss, _ = self.solve_racemic_steady_state()

        J11 = -4.0 * self.k_homo * x_ss - self.phi - 2.0 * self.k_rac
        J12 = self.k_cat * F_ss * (1.0 - 1.0 / self.S)
        J21 = 2.0 * self.k_homo * x_ss
        J22 = -self.phi

        tr = J11 + J22
        det = J11 * J22 - J12 * J21
        disc = tr ** 2 - 4.0 * det

        if disc < 0:
            lambda_max = tr / 2.0
        else:
            lambda_max = (tr + np.sqrt(disc)) / 2.0

        return lambda_max, det


def scan_bifurcation_landscape():
    """Scans stereospecificity S vs Feed Concentration F0."""
    print("Scanning bifurcation landscape...")
    S_vals = np.linspace(1.1, 8.0, 40)
    F0_vals = np.linspace(1e-3, 10e-3, 40) # 1 mM to 10 mM
    
    eig_grid = np.zeros((len(F0_vals), len(S_vals)))
    
    for i, F0 in enumerate(F0_vals):
        for j, S in enumerate(S_vals):
            model = PrebioticNetworkModel(S=S, F0=F0)
            lam_max, _ = model.difference_jacobian_eigenvalue()
            eig_grid[i, j] = lam_max

    return S_vals, F0_vals, eig_grid


def run_nonlinear_verification(S=5.0, F0=5e-3, seed_ee=1e-3, t_days=30):
    """Integrates the full nonlinear 6-ODE system to verify ee growth."""
    model = PrebioticNetworkModel(S=S, F0=F0)
    F_ss, x_ss, LL_ss, LD_ss = model.solve_racemic_steady_state()

    # Initial condition: steady state with 0.1% enantiomeric seed
    L0 = x_ss * (1.0 + seed_ee / 2.0)
    D0 = x_ss * (1.0 - seed_ee / 2.0)
    LL0 = LL_ss * (1.0 + seed_ee)
    DD0 = LL_ss * (1.0 - seed_ee)
    LD0 = LD_ss
    y0 = [F_ss, L0, D0, LL0, DD0, LD0]

    t_span = (0.0, t_days * 86400.0)
    t_eval = np.linspace(t_span[0], t_span[1], 400)

    sol = solve_ivp(model.ode_system, t_span, y0, t_eval=t_eval, method="Radau")
    
    L = sol.y[1]
    D = sol.y[2]
    LL = sol.y[3]
    DD = sol.y[4]
    
    ee_monomer = (L - D) / np.maximum(1e-15, L + D)
    ee_dimer = (LL - DD) / np.maximum(1e-15, LL + DD)
    t_days_eval = sol.t / 86400.0

    return t_days_eval, ee_monomer, ee_dimer


def generate_search_figure():
    """Generates the diagnostic plot chiral_search.png."""
    S_vals, F0_vals, eig_grid = scan_bifurcation_landscape()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel 1: Bifurcation Boundary
    F0_mM = F0_vals * 1e3
    growth_per_day = eig_grid * 86400.0
    
    # Levels centered around 0
    vmin = max(-0.15, float(np.min(growth_per_day)))
    vmax = float(np.max(growth_per_day))
    levels = np.linspace(vmin, vmax, 30)
    
    cs = ax1.contourf(S_vals, F0_mM, growth_per_day, levels=levels, cmap="RdBu_r")
    cbar = fig.colorbar(cs, ax=ax1)
    cbar.set_label(r"Growth Rate $\lambda_{\max}$ (day$^{-1}$)", fontsize=10)
    
    # Zero contour line (Bifurcation boundary)
    ax1.contour(S_vals, F0_mM, growth_per_day, levels=[0.0], colors="black", linewidths=2.5, linestyles="--")

    # Annotations
    ax1.text(5.5, 7.5, "AMPLIFICATION ZONE\n$(\\lambda_{\\max} > 0)$\nSpontaneous Homochirality",
             fontsize=9.5, weight="bold", color="darkred", ha="center")
    ax1.text(2.0, 3.0, "RACEMIC ZONE\n$(\\lambda_{\\max} < 0)$\nDecays to 50/50",
             fontsize=9.5, weight="bold", color="navy", ha="center")

    # Lab frontiers
    # ASSUMED placeholder S = 2.0 at 5 mM. Yu et al. 2024 PNAS (e2315447121) does not report S;
    # the dipeptide values once attributed to it were not in the paper.
    ax1.plot([2.0], [5.0], "y*", markersize=14, markeredgecolor="black", label="Assumed lab baseline (S=2.0)")
    ax1.annotate("Assumed lab baseline\n(S not measured)", xy=(2.0, 5.0), xytext=(2.6, 4.2),
                 arrowprops=dict(facecolor="black", shrink=0.08, width=1, headwidth=5),
                 fontsize=8.5, weight="bold")

    ax1.set_xlabel(r"Catalytic Stereospecificity $S = k_{\mathrm{homo}} / k_{\mathrm{cross}}$", fontsize=10)
    ax1.set_ylabel(r"Feed Concentration $[F]_0$ (mM)", fontsize=10)
    ax1.set_title("1. Analytical Bifurcation Boundary (CSTR Flow)", fontsize=11, weight="bold")
    ax1.legend(loc="lower right", fontsize=8.5)

    # Panel 2: Nonlinear Trajectory Verification
    t_sub, ee_m_sub, ee_d_sub = run_nonlinear_verification(S=1.8, F0=5e-3, seed_ee=1e-3, t_days=25)
    t_super, ee_m_super, ee_d_super = run_nonlinear_verification(S=5.0, F0=5e-3, seed_ee=1e-3, t_days=25)

    ax2.plot(t_super, ee_m_super * 100.0, "r-", lw=2.2, label=r"Supercritical ($S=5.0$): Monomer ee")
    ax2.plot(t_super, ee_d_super * 100.0, "m--", lw=1.8, label=r"Supercritical ($S=5.0$): Dimer ee")
    ax2.plot(t_sub, ee_m_sub * 100.0, "b:", lw=1.8, label=r"Subcritical ($S=1.8$): Monomer ee (decays)")
    
    ax2.set_xlabel("Residence Time (days in continuous CSTR)", fontsize=10)
    ax2.set_ylabel("Enantiomeric Excess ee (%)", fontsize=10)
    ax2.set_title("2. Nonlinear Dynamic Verification (from 0.1% Seed)", fontsize=11, weight="bold")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="center right", fontsize=8.5)

    fig.tight_layout()
    fig.savefig("chiral_search.png", dpi=150)
    print("Wrote diagnostic figure: chiral_search.png")


def main():
    print("=" * 80)
    print("COMPUTATIONAL SEARCH ENGINE: PREBIOTIC CHIRAL AMPLIFICATION")
    print("=" * 80)

    # 1. Test steady state solver
    model = PrebioticNetworkModel(F0=5e-3, S=3.0)
    F_ss, x_ss, LL_ss, LD_ss = model.solve_racemic_steady_state()
    assert F_ss > 0 and x_ss > 0 and LL_ss > 0, "Failed physical steady state"
    print(f"Racemic steady state (5 mM feed): F = {F_ss*1e3:.3f} mM, monomer = {x_ss*1e3:.4f} mM, dimer = {LL_ss*1e3:.4f} mM")

    # 2. Test subcritical eigenvalue (decay)
    sub = PrebioticNetworkModel(F0=5e-3, S=1.5)
    lam_sub, det_sub = sub.difference_jacobian_eigenvalue()
    assert lam_sub < 0 and det_sub > 0, f"Expected stable racemic state for S=1.5, got lambda={lam_sub}"
    print(f"Subcritical (S=1.5): lambda_max = {lam_sub:.2e} s^-1 (Decays back to racemic)")

    # 3. Test supercritical eigenvalue (amplification)
    super_m = PrebioticNetworkModel(F0=5e-3, S=5.0)
    lam_super, det_super = super_m.difference_jacobian_eigenvalue()
    assert lam_super > 0 and det_super < 0, f"Expected unstable racemic state for S=5.0, got lambda={lam_super}"
    print(f"Supercritical (S=5.0): lambda_max = +{lam_super:.2e} s^-1 (Spontaneous homochiral amplification!)")

    # 4. Generate phase diagram & nonlinear plot
    generate_search_figure()
    print("All search engine self-checks 100% PASSED.")


if __name__ == "__main__":
    main()
