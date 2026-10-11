"""Titan Cryogenic Chiral Amplification & Quantum Tunneling Model (94 K) (2026-10-10)

Tests Dr. Chris McKay's (NASA Ames) suggestion:
1. Classical Arrhenius freeze-out at Titan's surface temperature (94 K).
2. Quantum mechanical proton/hydrogen tunneling bypassing classical thermal barriers.
   Refs: H. W. Longo & Remsing 2024 arXiv:2405.19195, "Quantum Tunneling Could Enable
   Proton Transfer Reactions on Titan" (path-integral MD of proton-transfer thermodynamics
   in liquid ethane; it gives NO rate constants and no WKB, so every rate/barrier below is
   ASSUMED, not taken from it). Trixler 2013 Current Organic Chemistry 17:1758-1770,
   doi:10.2174/13852728113179990083 (generic tunnelling review; no rates or T_c values).
3. Temperature scaling of electroweak parity-violating energy difference (PVED): g ~ 1/T.
4. Coherence volume and Frank bifurcation thresholds in Titan's hydrocarbon lakes
   (Kraken Mare: V ~ 1.6e14 m^3, T = 94 K, liquid methane-ethane).

PHYSICAL EQUATIONS:
  - Electroweak bias: g(T) = Delta_E_pv / (k_B * T) = g(273 K) * (273.15 / T)
  - Classical Arrhenius rate: k_arr(T) = A * exp(-E_a / (R * T))
  - Crossover temperature: T_c = (hbar * omega_0) / (2 * pi * k_B) = (hbar * sqrt(V_0 / (m * d^2))) / (2 * pi * k_B)
  - Quantum tunneling transmission (Bell-Eckart / WKB inverted parabola):
      P_tunnel(E=0) = exp(-pi^2 * d / h * sqrt(2 * m * V_0))
  - Effective rate with tunneling: k_eff(T) = k_arr(T) + nu_0 * P_tunnel
  - Coherence parameter: Delta = PREF * g * (k_2 * c * tau)^0.25 * sqrt(N)

PREDICTIONS (written 2026-10-10, before the first run):
P1. Classical freeze-out: At 94 K, for any activation barrier E_a >= 40 kJ/mol,
    classical Arrhenius rates drop below 1e-18 /M/s. Over 1 Gyr, k_arr * c * tau << 1,
    making classical thermal chiral amplification strictly impossible (0% win rate).
P2. Quantum crossover: For hydrogen/proton transfer (m = 1.67e-27 kg, d = 0.6..0.9 A,
    V_0 = 15..35 kJ/mol), the crossover temperature T_c is 120..240 K.
    Since Titan's surface is 94 K (< T_c), all light-atom reactions are in the
    deep quantum tunneling regime, with rates exceeding classical Arrhenius by > 1e15x.
P3. Bifurcation threshold: Tunneling rate constants (1e-10 to 1e-6 /M/s, consistent with
    Longo & Remsing 2024) allow Frank networks in Kraken Mare to cross the bifurcation
    threshold (k_2 * c * tau >= 10) for prebiotic concentrations c >= 3.2 uM over 10^8 yr.
P4. PVED temperature boost: Because g ~ 1/T, g(94 K) = 2.91e-17 (2.91x higher than 273 K),
    reducing the critical coherence volume N_needed by (2.91)^2 = 8.47x.
RESULTS (titan_95k.log, corrected 2026-10-10). Four predictions (P1-P4) were made;
P5 below was added after the run and is not a prediction.
  P1 partly held: at 94 K, Ea = 40 kJ/mol gives 1-Gyr k*c*tau = 0.19 < 1 (no bifurcation), but
     k(94 K, 40 kJ/mol) = 5.9e-14 /M/s, NOT below 1e-18; the 1e-18 bound holds only from
     Ea >= 50 kJ/mol (k = 1.6e-19 /M/s).
  P2 FAILED: predicted T_c = 120-240 K; computed T_c = 78-95 K for d = 0.6-0.9 A. Titan (94 K)
     is near or above T_c (above 3 of 4 cases), not deep in the tunnelling regime. The tunnel/Arrhenius ratio at 94 K
     is 6.2e6x (15 kJ/mol) to 1.2e20x (50 kJ/mol), so the predicted > 1e15x also fails at low barriers.
  P3 held only for an ASSUMED rate: k_tunnel = 1e-8 /M/s is hand-set (Longo & Remsing 2024 gives
     no rate constants). With it, 3.2 uM over 100 Myr gives k*c*tau = 101 >= 10.
  P4 held: g(94 K) = 2.91e-17 (2.91x the 273 K value); molecules needed drop by 8.44x
     (predicted 8.47x, rounding). Kraken Mare volume margin at 1 uM is 78x.
  P5 (post hoc): proton-abstraction tunnelling t_1/2 = 2.5e6 yr (the limiting channel);
     heavy-atom inversion (m = m_C) t_1/2 = 6.5e67 yr. Racemization is slow, not zero.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

# Physical constants
R_GAS = 8.314462618        # J / (mol K)
K_B = 1.380649e-23         # J / K
H_BAR = 1.054571817e-34    # J s
H_PLANCK = 6.62607015e-34  # J s
M_PROTON = 1.6726219e-27   # kg (hydrogen nucleus)
M_CARBON = 1.9926468e-26   # kg
NA = 6.02214076e23         # 1 / mol
YR = 3.15576e7             # s / year
PREF = 1.103               # sqrt(2) * pi^(1/4) from sim.py

# Titan lake parameters (Kraken Mare)
T_TITAN = 94.0             # K
T_EARTH_0C = 273.15        # K
V_KRAKEN = 1.6e14          # m^3 (Lorenz et al. 2014)
AREA_KRAKEN = 5.0e11       # m^2 (500,000 km^2)
DEPTH_KRAKEN = 320.0       # m


def pved_g(T, g_273=1.0e-17):
    """Electroweak bias g = Delta_E_pv / (k_B * T), scaling as 1/T."""
    return g_273 * (T_EARTH_0C / T)


def crossover_temp(V0_kJ, d_angstrom, m=M_PROTON):
    """Crossover temperature T_c below which tunneling dominates Arrhenius kinetics."""
    V0_J = (V0_kJ * 1e3) / NA
    d_m = d_angstrom * 1e-10
    omega_0 = np.sqrt(V0_J / (m * d_m**2))
    return (H_BAR * omega_0) / (2 * np.pi * K_B)


def arrhenius_rate(T, Ea_kJ, A=1.0e9):
    """Classical Arrhenius rate constant (M^-1 s^-1 or s^-1)."""
    return A * np.exp(-(Ea_kJ * 1e3) / (R_GAS * T))


def tunneling_prob(V0_kJ, d_angstrom, m=M_PROTON):
    """WKB/Bell tunneling transmission probability through inverted parabolic barrier."""
    V0_J = (V0_kJ * 1e3) / NA
    d_m = d_angstrom * 1e-10
    exponent = (np.pi**2 * d_m / H_PLANCK) * np.sqrt(2 * m * V0_J)
    return np.exp(-exponent)


def effective_rate(T, Ea_kJ, d_angstrom, m=M_PROTON, A=1.0e9, nu_0=1.0e11):
    """Combined rate: Arrhenius thermal activation + quantum mechanical tunneling."""
    k_arr = arrhenius_rate(T, Ea_kJ, A=A)
    p_tun = tunneling_prob(Ea_kJ, d_angstrom, m=m)
    k_tun = nu_0 * p_tun  # s^-1 (converted to effective bimolecular with encounter prefactor)
    # At low T, k_eff is floored by tunneling
    return k_arr + k_tun


def n_needed_titan(T, k_tau, target_delta=2.0, g_273=1.0e-17):
    """Required coherent molecules for Delta >= target_delta at temperature T."""
    g = pved_g(T, g_273=g_273)
    return (target_delta / (PREF * g * k_tau**0.25))**2


def run_titan_audit():
    print("=" * 80)
    print("TITAN 94 K CRYOGENIC & QUANTUM TUNNELING AUDIT (Dr. Chris McKay Follow-Up)")
    print("=" * 80)

    # 1. PVED Temperature Scaling
    g_titan = pved_g(T_TITAN)
    g_earth = pved_g(T_EARTH_0C)
    g_boost = g_titan / g_earth
    print(f"\n[AUDIT 1] Electroweak PVED Scaling with Temperature:")
    print(f"  * Earth (273.15 K) : g = {g_earth:.3e}")
    print(f"  * Titan (94.0 K)   : g = {g_titan:.3e}")
    print(f"  * Titan PVED boost : {g_boost:.2f}x stronger bias at cryogenic temperatures!")

    # 2. Classical Arrhenius Freeze-Out Check (P1)
    print(f"\n[AUDIT 2] Classical Arrhenius Kinetics at 94 K vs 273 K (P1 check):")
    ea_values = [20.0, 40.0, 50.0, 80.0, 120.0]
    tau_age = 1.0e9 * YR  # 1 Gyr
    c_nom = 1.0e-4        # 100 uM
    for Ea in ea_values:
        k_273 = arrhenius_rate(T_EARTH_0C, Ea)
        k_94 = arrhenius_rate(T_TITAN, Ea)
        kct_94 = k_94 * c_nom * tau_age
        print(f"  * Ea = {Ea:5.1f} kJ/mol | k(273 K) = {k_273:10.2e} /M/s | k(94 K) = {k_94:10.2e} /M/s | 1-Gyr k*c*tau = {kct_94:10.2e}")
    k_40_94 = arrhenius_rate(T_TITAN, 40.0)
    k_50_94 = arrhenius_rate(T_TITAN, 50.0)
    assert k_40_94 * c_nom * tau_age < 1.0, "P1 check failed: at 40 kJ/mol 1-Gyr k*c*tau is below bifurcation threshold"
    assert k_50_94 < 1.0e-18, "P1 check failed: classical rate at 50 kJ/mol should be below 1e-18"
    print("  -> P1 partly held: at Ea = 40 kJ/mol 1-Gyr k*c*tau < 0.2 (frozen), but k(94 K) is 5.9e-14,")
    print("     not < 1e-18 /M/s; that bound holds only from Ea >= 50 kJ/mol.")

    # 3. Quantum Tunneling Crossover Temperature (P2)
    print(f"\n[AUDIT 3] Quantum Tunneling Crossover Temperatures T_c (P2 check):")
    proton_barriers = [
        (15.0, 0.6, "HCN tautomerization / fast proton transfer"),
        (25.0, 0.7, "H-transfer addition (assumed barrier)"),
        (35.0, 0.8, "Prebiotic condensation step"),
        (50.0, 0.9, "Enol-keto isomerization"),
    ]
    tc_list = []
    for V0, d, label in proton_barriers:
        Tc = crossover_temp(V0, d, m=M_PROTON)
        tc_list.append(Tc)
        P_tun = tunneling_prob(V0, d, m=M_PROTON)
        k_arr = arrhenius_rate(T_TITAN, V0)
        k_tun = 1e11 * P_tun
        print(f"  * {label:42}: V0={V0:4.1f} kJ, d={d:3.1f} A -> T_c = {Tc:5.1f} K | Tunnel/Arrhenius ratio = {k_tun / max(k_arr, 1e-99):.1e}")
    p2_held = all(120.0 <= tc <= 240.0 for tc in tc_list)  # the PREDICTED range
    # Regression check on what the model actually computes (recorded 78-95 K), not on P2.
    assert all(75.0 <= tc <= 100.0 for tc in tc_list), "T_c moved from the recorded 78-95 K"
    print(f"  -> P2 {'HELD' if p2_held else 'FAILED'}: predicted T_c = 120-240 K, computed "
          f"T_c = {min(tc_list):.0f}-{max(tc_list):.0f} K. Titan (94 K) is near or above T_c,")
    print("     not deep in the tunnelling regime.")

    # 4. Bifurcation Window in Kraken Mare (P3)
    print(f"\n[AUDIT 4] Tunneling-Driven Bifurcation in Kraken Mare (P3 check):")
    # ASSUMED: k_tunnel is hand-set. Longo & Remsing 2024 gives no rate constants.
    tau_100myr = 1.0e8 * YR
    c_thresh_p3 = 3.2e-6  # 3.2 uM
    k_tunnel_nom = 1.0e-8  # 1e-8 /M/s, ASSUMED
    kct_p3 = k_tunnel_nom * c_thresh_p3 * tau_100myr
    print(f"  * At k_tunnel = {k_tunnel_nom:.1e} /M/s, c = {c_thresh_p3*1e6:.1f} uM, tau = 100 Myr:")
    print(f"    k*c*tau = {kct_p3:.2f} (Threshold to cross bifurcation is k*c*tau >= 10.0)")
    assert kct_p3 >= 10.0, "P3 check failed: bifurcation threshold not crossed"
    print("  -> P3 holds only for the ASSUMED k_tunnel = 1e-8 /M/s (no measured or published rate).")

    # 5. Coherence Volume & Molecule Sizing (P4)
    print(f"\n[AUDIT 5] Coherence Volume & Lake Capacity (P4 check):")
    n_needed_94 = n_needed_titan(T_TITAN, k_tau=10.0, target_delta=2.0)
    n_needed_273 = n_needed_titan(T_EARTH_0C, k_tau=10.0, target_delta=2.0)
    n_ratio = n_needed_273 / n_needed_94
    print(f"  * Molecules needed at 273 K (Earth) : N = {n_needed_273:.3e}")
    print(f"  * Molecules needed at 94 K (Titan)  : N = {n_needed_94:.3e}")
    print(f"  * Coherence reduction factor        : {n_ratio:.2f}x fewer molecules needed at 94 K!")
    assert abs(n_ratio - g_boost**2) < 0.1, "P4 check failed: coherence scaling must match (g_boost)^2"

    vol_needed_1uM = n_needed_94 / (1.0e-6 * 1000 * NA)
    print(f"  * Volume required at 1 uM in Titan methane : {vol_needed_1uM:.2e} m^3 ({vol_needed_1uM/1e9:.3f} km^3)")
    print(f"  * Kraken Mare actual volume                : {V_KRAKEN:.2e} m^3 ({V_KRAKEN/1e9:.0f} km^3)")
    print(f"  * Kraken Mare capacity margin              : {V_KRAKEN / vol_needed_1uM:.0f}x greater than minimum required volume!")
    print("  -> P4 held: the 1/T scaling of g gives the predicted (g_boost)^2 reduction (algebraic).")

    # 6. Racemization Freedom Check (P5)
    print(f"\n[AUDIT 6] Racemization Longevity at 94 K (P5 check):")
    # Alpha-carbon stereocenter inversion involves:
    # 1. Light-atom proton abstraction (if basic solvent present): m = M_PROTON
    # 2. Heavy-atom carbon/substituent planar rehybridization (sp3 -> sp2): m = M_CARBON
    P_rac_H = tunneling_prob(120.0, 1.5, m=M_PROTON)
    P_rac_C = tunneling_prob(120.0, 1.5, m=M_CARBON)
    t_half_H_yr = np.log(2) / (1e11 * P_rac_H * YR) if P_rac_H > 0 else np.inf
    t_half_C_yr = np.log(2) / (1e11 * P_rac_C * YR) if P_rac_C > 0 else np.inf

    print(f"  * Proton abstraction tunneling limit (m=m_p) : P = {P_rac_H:.2e} | t_1/2 = {t_half_H_yr:.2e} yr")
    print(f"  * Heavy stereocenter inversion (m=m_C)       : P = {P_rac_C:.2e} | t_1/2 = {t_half_C_yr:.2e} yr")
    assert t_half_H_yr > 1.0e6, "P5 check failed: even light proton abstraction t_1/2 must exceed 1 Myr at 94 K"
    assert t_half_C_yr > 1.0e30, "P5 check failed: heavy atom inversion must be completely frozen (> 1e30 yr)"
    print(f"  -> P5 (post hoc, not a prediction): the limiting channel is proton abstraction,")
    print(f"     t_1/2 = {t_half_H_yr:.1e} yr; heavy-atom inversion alone is {t_half_C_yr:.1e} yr.")
    print("     Racemization is slow on Myr scales, not frozen out.")

    # Plot comprehensive diagnostic figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11))
    fig.patch.set_facecolor("#ffffff")

    # Panel 1: Arrhenius vs Tunneling Rates across Temperature
    temps = np.linspace(50, 350, 300)
    Ea_ex = 30.0  # kJ/mol
    d_ex = 0.7    # A
    k_arr_curve = [arrhenius_rate(T, Ea_ex) for T in temps]
    P_tun_ex = tunneling_prob(Ea_ex, d_ex)
    k_tun_curve = [1e11 * P_tun_ex] * len(temps)
    k_tot_curve = [arr + 1e11 * P_tun_ex for arr in k_arr_curve]

    ax1.semilogy(temps, k_arr_curve, "r--", linewidth=2.0, label=f"Classical Arrhenius ($E_a={Ea_ex:.0f}$ kJ/mol)")
    ax1.semilogy(temps, k_tun_curve, "b--", linewidth=2.0, label=f"Quantum Tunneling Plateau ($d={d_ex:.1f}$ Å)")
    ax1.semilogy(temps, k_tot_curve, "k-", linewidth=2.5, label="Total Effective Rate $k_{\\mathrm{eff}}$")
    ax1.axvline(T_TITAN, color="darkorange", linestyle=":", linewidth=2.5, label=f"Titan Surface (94 K)")
    ax1.axvline(crossover_temp(Ea_ex, d_ex), color="purple", linestyle="-.", linewidth=2.0, label=f"Crossover $T_c$ ({crossover_temp(Ea_ex, d_ex):.0f} K)")
    ax1.set_xlabel("Temperature (K)", fontsize=11)
    ax1.set_ylabel("Rate Constant ($s^{-1}$)", fontsize=11)
    ax1.set_title("1. Quantum Tunneling Bypass of Arrhenius Freeze-out", fontsize=12, weight="bold")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="lower right", fontsize=9)
    ax1.set_ylim(1e-25, 1e12)

    # Panel 2: PVED Scaling and Coherence Threshold
    T_range = np.linspace(60, 320, 200)
    g_curve = [pved_g(T) for T in T_range]
    n_curve = [n_needed_titan(T, 10.0) for T in T_range]

    ax2_twin = ax2.twinx()
    l1 = ax2.plot(T_range, np.array(g_curve) * 1e17, "g-", linewidth=2.5, label="Electroweak Bias $g(T) \\times 10^{17}$")
    l2 = ax2_twin.semilogy(T_range, n_curve, "m-", linewidth=2.5, label="Molecules Needed $N_{\\mathrm{crit}}$")
    ax2.axvline(T_TITAN, color="darkorange", linestyle=":", linewidth=2.5, label="Titan (94 K)")
    ax2.axvline(T_EARTH_0C, color="blue", linestyle=":", linewidth=2.0, label="Earth (273 K)")
    ax2.set_xlabel("Temperature (K)", fontsize=11)
    ax2.set_ylabel("Electroweak Tilt $g \\times 10^{17}$", color="green", fontsize=11)
    ax2_twin.set_ylabel("Molecules Needed $N_{\\mathrm{crit}}$", color="purple", fontsize=11)
    ax2.set_title("2. Cryogenic PVED Enhancement ($g \\propto 1/T$)", fontsize=12, weight="bold")
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc="upper center", fontsize=9)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Concentration vs Time Bifurcation Phase Diagram at 94 K
    conc_grid = np.logspace(-7, -2, 100)  # 0.1 uM to 10 mM
    time_grid = np.logspace(4, 9, 100)   # 10 kyr to 1 Gyr
    C_mesh, T_mesh = np.meshgrid(conc_grid, time_grid)
    KCT = k_tunnel_nom * C_mesh * (T_mesh * YR)

    cs = ax3.contourf(conc_grid * 1e6, time_grid / 1e6, KCT, levels=[0, 1, 10, 1e2, 1e4, 1e6], cmap="viridis")
    cbar = plt.colorbar(cs, ax=ax3)
    cbar.set_label("Dynamical Progress $k_2 c \\tau$", fontsize=10)
    ax3.contour(conc_grid * 1e6, time_grid / 1e6, KCT, levels=[10.0], colors="red", linewidths=2.5)
    ax3.scatter([10.0], [100.0], color="red", s=100, zorder=5, label="Titan Baseline (10 uM, 100 Myr)")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel("Reactant Concentration ($\\mu$M)", fontsize=11)
    ax3.set_ylabel("Reaction Timescale (Myr)", fontsize=11)
    ax3.set_title("3. Titan 94 K Bifurcation Boundary ($k_2 c \\tau \\geq 10$)", fontsize=12, weight="bold")
    ax3.legend(loc="upper left", fontsize=9)
    ax3.grid(True, alpha=0.3)

    # Panel 4: Comparison of Titan 94 K vs Earth 273 K vs Enceladus 273 K
    environments = ["Earth (273 K)", "Enceladus (273 K)", "Titan Melt (250 K)", "Titan Surface (94 K)"]
    g_vals = [pved_g(273.15)*1e17, pved_g(273.15)*1e17, pved_g(250.0)*1e17, pved_g(94.0)*1e17]
    # ASSUMED illustrative racemization rates (/s); Titan uses this script's own
    # proton-channel tunnelling rate (the limiting channel), not a hand-set 1e-35.
    rac_rates = [1e-12, 1e-12, 1e-13, np.log(2) / (t_half_H_yr * YR)]

    x_pos = np.arange(len(environments))
    width = 0.35

    ax4.bar(x_pos - width/2, g_vals, width, label="Electroweak Bias $g \\times 10^{17}$", color="#2b5c8f")
    ax4_twin = ax4.twinx()
    ax4_twin.bar(x_pos + width/2, [-np.log10(max(r, 1e-35)) for r in rac_rates], width, label="Chiral Stability $-\\log_{10}(k_{\\mathrm{rac}})$", color="#d95f02")

    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(environments, rotation=15, ha="right", fontsize=9)
    ax4.set_ylabel("PVED Tilt ($g \\times 10^{17}$)", color="#2b5c8f", fontsize=11)
    ax4_twin.set_ylabel("Racemization Resistance ($-\\log_{10} k_{\\mathrm{rac}}$)", color="#d95f02", fontsize=11)
    ax4.set_title("4. Environmental Parameter Comparison", fontsize=12, weight="bold")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = "titan_95k_audit.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"\nSaved Titan audit master figure to: {fig_path}")
    print("=" * 80)


if __name__ == "__main__":
    run_titan_audit()
