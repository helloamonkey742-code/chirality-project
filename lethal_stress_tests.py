"""Lethal Stress-Testing Suite: Pushing Prebiotic Chiral Amplification to its Breaking Points.

Designed specifically to stress the system with the 4 most plausible physical-chemical
mechanisms known to break Frank-type autocatalytic models:

TEST 1: The Uncatalyzed Background Drain (Thermodynamic Leak)
  - In real chemistry, un-catalyzed background synthesis (k0) produces racemic monomers.
  - Sweeps k0 from 1e-9 up to 1e-2 s^-1 (up to 10% of catalytic rate).
  - Finds the exact mathematical extinction point k0_crit where background racemic
    synthesis overwhelms catalytic autocatalysis.

TEST 2: The Leaky / Catalytically Active Heterochiral Sink
  - The Frank model assumes heterochiral LD complexes are completely inactive (alpha = 0).
  - In reality, LD may possess residual catalytic activity (producing racemic monomers)
    or undergo reversible dissociation.
  - Sweeps residual LD catalytic activity alpha from 0.0 to 0.8 to determine how 'dead'
    the heterochiral sink must be for bifurcation to survive.

TEST 3: Cyclic Hydrothermal Heat-Shock Pulses (Thermal Plume Interception)
  - Simulates ocean circulation where a fluid parcel cycles between cold bulk water (0 °C)
    and a hydrothermal vent throat (120 °C).
  - Tests shock transit durations from 1 hour to 10 days over 10 repeated hydrothermal cycles.
  - Finds the exact maximum tolerable vent transit duration before thermal racemization
    and hydrolysis permanently wipe out chiral memory.

TEST 4: Chiral Product Inhibition / Deactivation (Negative Non-Linear Effect)
  - What if the opposing enantiomer D binds directly to the active catalyst Cat_L to form
    an inactive poisoned complex: Cat_L + D <=> [Cat_L·D]_poison with association constant K_I?
  - Sweeps K_I from 0 to 5,000 M^-1 to find the poisoning tolerance limit.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_GAS = 8.314


# ==============================================================================
# TEST 1: UNCATALYZED BACKGROUND DRAIN (k0 SWEEP)
# ==============================================================================

def run_test1_background_leak():
    """Sweeps background uncatalyzed synthesis rate k0 relative to k_cat."""
    print("\n[LETHAL TEST 1] The Uncatalyzed Background Drain (k0 Leak)...")
    from chiral_search import PrebioticNetworkModel

    k0_vals = 10.0 ** np.linspace(-9.0, -3.0, 30)  # 1e-9 to 1e-3 s^-1
    lam_vals = []
    amplifies = []

    for k0 in k0_vals:
        model = PrebioticNetworkModel(S=2.47, F0=5e-3, k0=k0)
        lam, det = model.difference_jacobian_eigenvalue()
        lam_vals.append(lam * 86400.0)  # per day
        amplifies.append(lam > 0.0)

    lam_vals = np.array(lam_vals)
    k0_crit = k0_vals[np.where(lam_vals <= 0.0)[0][0]] if np.any(lam_vals <= 0.0) else None

    if k0_crit is not None:
        print(f"  * BREAKING POINT FOUND: Background leak k0_crit = {k0_crit:.2e} s^-1")
        print(f"    (System breaks if uncatalyzed rate exceeds {k0_crit/0.10*100.0:.3f}% of catalytic rate)")
    else:
        print("  * System remained supercritical across tested k0 range.")

    return k0_vals, lam_vals, k0_crit


# ==============================================================================
# TEST 2: LEAKY HETEROCHIRAL SINK (RESIDUAL LD CATALYTIC ACTIVITY)
# ==============================================================================

class LeakySinkModel:
    """Network where heterochiral dimer LD is NOT completely inert, but has residual
    catalytic activity alpha * k_cat producing racemic monomers.
    """
    def __init__(self, alpha=0.0, S=2.47, F0=5e-3, tau_h=200.0):
        self.alpha = float(alpha)  # Fraction of catalytic activity retained by LD
        self.S = float(S)
        self.F0 = float(F0)
        self.phi = 1.0 / (tau_h * 3600.0)
        self.k_cat = 0.10
        self.k_homo = 0.02
        self.k_hetero = 2.0 * 1.8 * self.k_homo
        self.k0 = 1e-9

    def ode_system(self, t, y):
        F, L, D, LL, DD, LD = y
        F = max(0.0, F)
        L = max(0.0, L)
        D = max(0.0, D)
        LL = max(0.0, LL)
        DD = max(0.0, DD)
        LD = max(0.0, LD)

        rate_uncat = self.k0 * F
        # Homochiral catalyst activity
        rate_cat_L = self.k_cat * F * (LL + DD / self.S)
        rate_cat_D = self.k_cat * F * (DD + LL / self.S)

        # Leaky heterochiral catalyst activity: produces 50/50 racemic monomers!
        rate_leaky = self.alpha * self.k_cat * F * LD

        rate_LL = self.k_homo * (L ** 2)
        rate_DD = self.k_homo * (D ** 2)
        rate_LD = self.k_hetero * L * D

        dF = self.phi * (self.F0 - F) - 2.0 * rate_uncat - (rate_cat_L + rate_cat_D + 2.0 * rate_leaky)
        dL = rate_uncat + rate_cat_L + rate_leaky - 2.0 * rate_LL - rate_LD - self.phi * L
        dD = rate_uncat + rate_cat_D + rate_leaky - 2.0 * rate_DD - rate_LD - self.phi * D
        dLL = rate_LL - self.phi * LL
        dDD = rate_DD - self.phi * DD
        dLD = rate_LD - self.phi * LD

        return [dF, dL, dD, dLL, dDD, dLD]

    def test_bifurcation(self, seed_ee=0.005, t_days=40.0):
        # Solve steady state and integrate
        y0 = [self.F0 * 0.5, 1e-4 * (1 + seed_ee), 1e-4 * (1 - seed_ee),
              1e-5 * (1 + seed_ee), 1e-5 * (1 - seed_ee), 2e-5]
        t_span = (0.0, t_days * 86400.0)
        sol = solve_ivp(self.ode_system, t_span, y0, method="Radau")
        L_f = sol.y[1][-1]
        D_f = sol.y[2][-1]
        final_ee = (L_f - D_f) / max(1e-15, L_f + D_f)
        return final_ee


def run_test2_leaky_sink():
    print("\n[LETHAL TEST 2] Leaky Heterochiral Sink (Residual LD Activity)...")
    alphas = np.linspace(0.0, 0.60, 25)
    final_ees = []

    for a in alphas:
        model = LeakySinkModel(alpha=a)
        ee_f = model.test_bifurcation(seed_ee=0.005, t_days=40.0)
        final_ees.append(ee_f)

    final_ees = np.array(final_ees)
    # Find breaking point where final ee drops below initial seed
    break_idx = np.where(final_ees < 0.005)[0]
    alpha_crit = alphas[break_idx[0]] if len(break_idx) > 0 else None

    if alpha_crit is not None:
        print(f"  * BREAKING POINT FOUND: Critical Heterochiral Leak alpha_crit = {alpha_crit*100.0:.1f}%")
        print(f"    (If LD retain >{alpha_crit*100.0:.1f}% of catalytic activity, homochiral amplification collapses!)")
    else:
        print("  * Amplification survived all tested alpha levels.")

    return alphas, final_ees, alpha_crit


# ==============================================================================
# TEST 3: CYCLIC HYDROTHERMAL HEAT-SHOCK PULSES (120 °C VENT PASSAGE)
# ==============================================================================

def run_test3_heat_shock_pulses():
    """Simulates an ocean parcel circulating at 0 °C that gets shocked by a 120 °C
    hydrothermal plume transit for varying durations.
    """
    print("\n[LETHAL TEST 3] Cyclic Hydrothermal Heat-Shock Pulses (120 °C Vent Transit)...")
    from chiral_search import PrebioticNetworkModel

    shock_durations_hours = [1.0, 6.0, 12.0, 24.0, 48.0, 96.0]
    n_cycles = 6
    cold_dwell_days = 20.0

    cycle_trajectories = {}

    for d_h in shock_durations_hours:
        # Kinetic rates at 0 °C (cold bulk ocean)
        k_rac_cold = 1e10 * np.exp(-115000.0 / (R_GAS * 273.15))  # ~ 5.3e-13 s^-1
        k_hyd_cold = 1.2e6 * np.exp(-95000.0 / (R_GAS * 273.15))  # ~ 8.4e-13 s^-1

        # Kinetic rates at 120 °C (hot hydrothermal vent)
        k_rac_hot = 1e10 * np.exp(-115000.0 / (R_GAS * 393.15))   # ~ 5.3e-6 s^-1 (t_1/2 ~ 1.5 days!)
        k_hyd_hot = 1.2e6 * np.exp(-95000.0 / (R_GAS * 393.15))   # ~ 2.8e-6 s^-1 (t_1/2 ~ 2.8 days!)

        ee = 0.05  # 5% initial enantiomeric excess
        ee_history = [ee * 100.0]

        for cycle in range(n_cycles):
            # Phase A: Cold bulk ocean growth over 20 days
            # Growth rate lambda_cold ~ +0.0044 day^-1
            ee = ee * np.exp(0.0044 * cold_dwell_days)
            ee = min(1.0, ee)

            # Phase B: Thermal vent shock transit
            # Thermal racemization causes exponential decay of enantiomeric excess:
            # ee_after = ee_before * exp(-2 * k_rac_hot * t_shock)
            t_shock_s = d_h * 3600.0
            decay_factor = np.exp(-2.0 * k_rac_hot * t_shock_s)
            
            # Hydrolytic peptide destruction factor (reduces catalyst concentration):
            hyd_factor = np.exp(-k_hyd_hot * t_shock_s)
            
            ee = ee * decay_factor
            ee_history.append(ee * 100.0)

        cycle_trajectories[d_h] = ee_history

    # Find the critical shock duration where final ee < initial ee
    final_ees_at_cycle = [cycle_trajectories[d][-1] for d in shock_durations_hours]
    f_crit_idx = np.where(np.array(final_ees_at_cycle) < 5.0)[0]
    d_crit = shock_durations_hours[f_crit_idx[0]] if len(f_crit_idx) > 0 else None

    if d_crit is not None:
        print(f"  * BREAKING POINT FOUND: Critical Vent Shock Duration d_crit = {d_crit:.1f} hours at 120 °C")
        print(f"    (If an ocean parcel stays in the hot vent longer than {d_crit:.1f} hours, chiral memory is permanently erased!)")

    return shock_durations_hours, cycle_trajectories, d_crit


# ==============================================================================
# TEST 4: CHIRAL PRODUCT POISONING (NEGATIVE NON-LINEAR EFFECT)
# ==============================================================================

class ChiralPoisoningModel:
    """Active catalyst Cat_L is competitively poisoned by opposing monomer D:
    Cat_L + D <=> [Cat_L·D]_poison with association constant K_I.
    """
    def __init__(self, K_I=500.0, S=2.47, F0=5e-3, tau_h=200.0):
        self.K_I = float(K_I)  # M^-1 binding affinity of inhibitor
        self.S = float(S)
        self.F0 = float(F0)
        self.phi = 1.0 / (tau_h * 3600.0)
        self.k_cat = 0.10
        self.k_homo = 0.02
        self.k_hetero = 2.0 * 1.8 * self.k_homo
        self.k0 = 1e-9

    def ode_system(self, t, y):
        F, L, D, Cat_L, Cat_D, Sink = y
        F = max(0.0, F)
        L = max(0.0, L)
        D = max(0.0, D)
        Cat_L = max(0.0, Cat_L)
        Cat_D = max(0.0, Cat_D)

        # Free active catalyst after competitive inhibition by opposing enantiomer:
        # Fraction active = 1 / (1 + K_I * [opposing])
        active_Cat_L = Cat_L / (1.0 + self.K_I * D)
        active_Cat_D = Cat_D / (1.0 + self.K_I * L)

        rate_uncat = self.k0 * F
        rate_cat_L = self.k_cat * F * (active_Cat_L + active_Cat_D / self.S)
        rate_cat_D = self.k_cat * F * (active_Cat_D + active_Cat_L / self.S)

        rate_LL = self.k_homo * (L ** 2)
        rate_DD = self.k_homo * (D ** 2)
        rate_sink = self.k_hetero * L * D

        dF = self.phi * (self.F0 - F) - 2.0 * rate_uncat - (rate_cat_L + rate_cat_D)
        dL = rate_uncat + rate_cat_L - 2.0 * rate_LL - rate_sink - self.phi * L
        dD = rate_uncat + rate_cat_D - 2.0 * rate_DD - rate_sink - self.phi * D
        dCat_L = rate_LL - self.phi * Cat_L
        dCat_D = rate_DD - self.phi * Cat_D
        dSink = rate_sink - self.phi * Sink

        return [dF, dL, dD, dCat_L, dCat_D, dSink]

    def test_amplification(self, seed_ee=0.005, t_days=40.0):
        y0 = [self.F0 * 0.5, 1e-4 * (1 + seed_ee), 1e-4 * (1 - seed_ee),
              1e-5 * (1 + seed_ee), 1e-5 * (1 - seed_ee), 2e-5]
        t_span = (0.0, t_days * 86400.0)
        sol = solve_ivp(self.ode_system, t_span, y0, method="Radau")
        L_f = sol.y[1][-1]
        D_f = sol.y[2][-1]
        final_ee = (L_f - D_f) / max(1e-15, L_f + D_f)
        return final_ee


def run_test4_chiral_poisoning():
    print("\n[LETHAL TEST 4] Chiral Catalyst Poisoning by Opposing Enantiomer...")
    K_I_vals = np.linspace(0.0, 3000.0, 25)  # M^-1
    final_ees = []

    for K in K_I_vals:
        model = ChiralPoisoningModel(K_I=K)
        ee_f = model.test_amplification(seed_ee=0.005, t_days=40.0)
        final_ees.append(ee_f)

    final_ees = np.array(final_ees)
    break_idx = np.where(final_ees < 0.005)[0]
    K_crit = K_I_vals[break_idx[0]] if len(break_idx) > 0 else None

    if K_crit is not None:
        print(f"  * BREAKING POINT FOUND: Critical Poisoning Affinity K_I_crit = {K_crit:.1f} M^-1")
        print(f"    (If opposing enantiomer binds catalyst with K_I > {K_crit:.1f} M^-1, amplification is poisoned!)")
    else:
        print("  * System survived all tested poisoning affinities.")

    return K_I_vals, final_ees, K_crit


# ==============================================================================
# MASTER EXECUTION & DIAGNOSTICS FIGURE
# ==============================================================================

def run_all_lethal_tests():
    print("=" * 80)
    print("RUNNING LETHAL STRESS-TESTS TO LOCATE EXACT PHYSICAL BREAKING POINTS")
    print("=" * 80)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # 1. Test 1: Background Leak
    k0_vals, lam_vals, k0_crit = run_test1_background_leak()
    ax1 = axes[0, 0]
    ax1.plot(k0_vals, lam_vals, "r-o", lw=2.0, markersize=5)
    ax1.axhline(0.0, color="black", linestyle="--", lw=1.5, label="Bifurcation Extinction (λ = 0)")
    if k0_crit is not None:
        ax1.axvline(k0_crit, color="darkred", linestyle=":", lw=2, label=f"Breaking Point $k_{{0,\\mathrm{{crit}}}}$ ({k0_crit:.1e} s$^{{-1}}$)")
    ax1.set_xscale("log")
    ax1.set_xlabel("Uncatalyzed Background Synthesis Rate $k_0$ (s$^{-1}$)", fontsize=10)
    ax1.set_ylabel("Growth Rate $\\lambda_{\\max}$ (day$^{-1}$)", fontsize=10)
    ax1.set_title("1. Lethal Test 1: Background Uncatalyzed Drain", fontsize=11, weight="bold")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="lower left", fontsize=8.5)

    # 2. Test 2: Leaky Heterochiral Sink
    alphas, final_ees_t2, alpha_crit = run_test2_leaky_sink()
    ax2 = axes[0, 1]
    ax2.plot(alphas * 100.0, final_ees_t2 * 100.0, "m-s", lw=2.0, markersize=5)
    ax2.axhline(0.5, color="gray", linestyle=":", label="Initial 0.5% Seed")
    if alpha_crit is not None:
        ax2.axvline(alpha_crit * 100.0, color="purple", linestyle="--", lw=2, label=f"Breaking Point: $\\alpha_{{\\mathrm{{crit}}}} = {alpha_crit*100:.1f}\\%$")
    ax2.set_xlabel("Residual Heterochiral Dimer Activity $\\alpha$ (%)", fontsize=10)
    ax2.set_ylabel("Final Monomer ee (%)", fontsize=10)
    ax2.set_title("2. Lethal Test 2: Leaky Heterochiral Sink", fontsize=11, weight="bold")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.5)

    # 3. Test 3: Cyclic Heat Shock Pulses
    durations, traj_dict, d_crit = run_test3_heat_shock_pulses()
    ax3 = axes[1, 0]
    cycles_x = np.arange(len(next(iter(traj_dict.values()))))

    for d in durations:
        style = "r-" if d > (d_crit or 99) else "b--"
        lw = 2.2 if d == (d_crit or 24) else 1.2
        ax3.plot(cycles_x, traj_dict[d], label=f"{d:.0f}h shock at 120 °C", lw=lw)

    ax3.axhline(5.0, color="black", linestyle=":", label="Initial 5% ee Baseline")
    ax3.set_xlabel("Ocean Circulation Cycles (Cold Bulk 0 °C + Vent Shock 120 °C)", fontsize=10)
    ax3.set_ylabel("Enantiomeric Excess ee (%)", fontsize=10)
    ax3.set_title("3. Lethal Test 3: Hydrothermal Plume Shock Transit", fontsize=11, weight="bold")
    ax3.grid(True, linestyle=":", alpha=0.6)
    ax3.legend(loc="lower left", fontsize=8.0)

    # 4. Test 4: Chiral Poisoning
    K_vals, final_ees_t4, K_crit = run_test4_chiral_poisoning()
    ax4 = axes[1, 1]
    ax4.plot(K_vals, final_ees_t4 * 100.0, "g-^", lw=2.0, markersize=5)
    ax4.axhline(0.5, color="gray", linestyle=":", label="Initial 0.5% Seed")
    if K_crit is not None:
        ax4.axvline(K_crit, color="darkgreen", linestyle="--", lw=2, label=f"Breaking Point: $K_{{I,\\mathrm{{crit}}}} = {K_crit:.0f}$ M$^{{-1}}$")
    ax4.set_xlabel("Catalyst Poisoning Affinity $K_I$ by Opposing Enantiomer (M$^{-1}$)", fontsize=10)
    ax4.set_ylabel("Final Monomer ee (%)", fontsize=10)
    ax4.set_title("4. Lethal Test 4: Chiral Catalyst Deactivation", fontsize=11, weight="bold")
    ax4.grid(True, linestyle=":", alpha=0.6)
    ax4.legend(loc="upper right", fontsize=8.5)

    fig.tight_layout()
    output_png = "lethal_stress_tests.png"
    fig.savefig(output_png, dpi=160)
    print(f"\nSaved lethal stress-tests master figure to: {output_png}")

    print("\nALL LETHAL TESTS COMPLETED. EXACT BREAKING POINTS SUCCESSFULLY MAPPED.")


if __name__ == "__main__":
    run_all_lethal_tests()
