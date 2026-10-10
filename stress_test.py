"""Stress test: Quantitative modeling of Red Team critical attacks on icy-moon chirality.

Computes:
1. Arrhenius racemization timescales for Alanine, Isovaline, and Aspartic Acid across
   0 °C to 150 °C vs. the layer-sinking timescale H*R / (2*D_z) for Enceladus, Europa,
   and early Earth.
2. Hydrothermal transit destruction: survival of chiral excess during vent circulation.
3. Dimensionless fluid-dynamical ratios between Allen-Cahn curvature wall velocity
   (2*D_z/R) and turbulent eddy shear / chaotic advection.
"""

import numpy as np

# Physical constants
R_GAS = 8.314462618     # J/(mol*K)
YR_TO_SEC = 3.15576e7   # seconds per year

# World geometric parameters (Radius R in m, Depth H in m, Volume V in m^3)
# Sourced from uncertainty.py / ocean.py / Čadek et al. 2016
WORLDS = {
    "Enceladus": {
        "R": 232e3,             # 232 km radius
        "H": 3.99e4,            # ~40 km ocean depth
        "V": 2.7e16,            # 2.7e16 m^3
        "age_min": 1e6,         # 1 Myr (young ocean hypothesis)
        "age_max": 1e9,         # 1 Gyr (ancient ocean hypothesis)
        "u_eddy": 1e-3,         # 1 mm/s typical eddy velocity (Kang et al. 2022)
        "L_eddy": 20e3,         # 20 km Rossby deformation radius / eddy size
        "tau_hydro": 2.5e5 * YR_TO_SEC, # Whole-ocean hydrothermal turnover ~250 kyr (Choblet et al. 2017)
    },
    "Europa": {
        "R": 1460e3,            # 1460 km radius
        "H": 1.23e5,            # ~123 km ocean depth
        "V": 3.3e18,            # 3.3e18 m^3
        "age_min": 1e8,         # 100 Myr
        "age_max": 4e9,         # 4 Gyr
        "u_eddy": 1e-2,         # 1 cm/s typical eddy velocity (Soderlund et al. 2014, 2020)
        "L_eddy": 50e3,         # 50 km eddy scale
        "tau_hydro": 1e7 * YR_TO_SEC, # Hydrothermal circulation ~10 Myr
    },
    "Early Earth": {
        "R": 6.371e6,           # 6371 km radius
        "H": 2.61e3,            # ~2.6 km average depth
        "V": 1.332e18,          # 1.332e18 m^3
        "age_min": 1e7,
        "age_max": 5e8,
        "u_eddy": 5e-2,         # 5 cm/s mesoscale eddy velocity
        "L_eddy": 100e3,        # 100 km deformation radius
        "tau_hydro": 1e6 * YR_TO_SEC, # Whole-ocean vent cycling ~1 Myr
    }
}

# Racemization Arrhenius parameters
# Rate equation for excess decay: de/dt = -2*k_r*e => tau_half = ln(2) / (2*k_r)
# Alanine: Bada 1972, Schroeder & Bada 1976 (Ea = 130 kJ/mol, A = 3.0e13 s^-1)
# Aspartic acid: Bada 1972 (Ea = 115 kJ/mol, A = 2.4e12 s^-1)
# Isovaline: alpha-dialkyl, no alpha-H. Thermal racemization via C-C homolysis
# (Ea ~ 270 kJ/mol, A ~ 1e15 s^-1; Cohen & Chyba 2000, Glavin & Dworkin 2009)
AMINO_ACIDS = {
    "Alanine (alpha-H)": {
        "Ea": 130.0e3,          # J/mol
        "A": 3.0e13,            # s^-1
        "notes": "Representative proteinogenic amino acid; carbanion mechanism"
    },
    "Aspartic acid (fast alpha-H)": {
        "Ea": 115.0e3,          # J/mol
        "A": 2.4e12,            # s^-1
        "notes": "Rapid racemizer due to beta-carboxyl electron withdrawal"
    },
    "Isovaline (alpha-methyl)": {
        "Ea": 270.0e3,          # J/mol (C-C cleavage barrier)
        "A": 1.0e15,            # s^-1
        "notes": "Non-proteinogenic, lacks alpha-H; immune to thermal deprotonation"
    }
}


def k_racemization(Ea, A, T_celsius):
    """Compute first-order enantiomeric interconversion rate constant k_r in s^-1."""
    T_kelvin = T_celsius + 273.15
    return A * np.exp(-Ea / (R_GAS * T_kelvin))


def t_half_racemization(k_r):
    """Half-life of enantiomeric excess in seconds: e(t) = e0 * exp(-2*k_r*t)."""
    return np.log(2.0) / (2.0 * k_r)


def layer_sinking_time(H, R, Dz):
    """Curvature-driven layer-sinking timescale H * R / (2 * Dz) in seconds."""
    return (H * R) / (2.0 * Dz)


def curvature_wall_velocity(Dz, R):
    """Allen-Cahn curvature wall velocity 2 * Dz / R in m/s."""
    return 2.0 * Dz / R


def run_racemization_analysis():
    print("=" * 80)
    print("1. ARRHENIUS RACEMIZATION TIMESCALES ACROSS TEMPERATURE (0 °C - 150 °C)")
    print("=" * 80)
    temps = [0.0, 25.0, 50.0, 75.0, 100.0, 150.0]

    for name, params in AMINO_ACIDS.items():
        print(f"\n--- {name} (Ea = {params['Ea']/1e3:.1f} kJ/mol, A = {params['A']:.1e} s^-1) ---")
        print(f"Notes: {params['notes']}")
        print(f"{'Temp (°C)':>10} {'k_r (s^-1)':>14} {'k_r (yr^-1)':>14} {'Half-life':>18} {'Decay timescale':>18}")
        for t_c in temps:
            kr = k_racemization(params["Ea"], params["A"], t_c)
            t_half = t_half_racemization(kr)
            t_half_yr = t_half / YR_TO_SEC
            tau_decay_yr = (1.0 / (2.0 * kr)) / YR_TO_SEC

            if t_half_yr > 1e12:
                th_str = f"> 1e12 yr"
                td_str = f"> 1e12 yr"
            elif t_half_yr >= 1.0:
                th_str = f"{t_half_yr:10.2e} yr"
                td_str = f"{tau_decay_yr:10.2e} yr"
            elif t_half_yr * 365.25 >= 1.0:
                th_str = f"{t_half_yr * 365.25:10.2f} days"
                td_str = f"{tau_decay_yr * 365.25:10.2f} days"
            elif t_half_yr * 365.25 * 24.0 >= 1.0:
                th_str = f"{t_half_yr * 365.25 * 24.0:10.2f} hours"
                td_str = f"{tau_decay_yr * 365.25 * 24.0:10.2f} hours"
            else:
                th_str = f"{t_half:10.2f} sec"
                td_str = f"{1.0/(2.0*kr):10.2f} sec"

            print(f"{t_c:10.1f} {kr:14.3e} {kr*YR_TO_SEC:14.3e} {th_str:>18} {td_str:>18}")


def run_sinking_vs_racemization_analysis():
    print("\n" + "=" * 80)
    print("2. LAYER-SINKING TIMESCALES (H*R / 2*D_z) VS. RACEMIZATION LIFETIMES")
    print("=" * 80)
    dz_values = [1e-10, 1e-7, 1e-5, 1e-3]

    print(f"\n{'World':>12} {'D_z (m^2/s)':>14} {'v_wall (m/s)':>14} {'v_wall (mm/yr)':>16} {'t_sink (yr)':>16}")
    for world_name, data in WORLDS.items():
        print(f"--- {world_name} (R = {data['R']/1e3:.0f} km, H = {data['H']/1e3:.1f} km) ---")
        for dz in dz_values:
            v_wall = curvature_wall_velocity(dz, data["R"])
            t_sink_sec = layer_sinking_time(data["H"], data["R"], dz)
            t_sink_yr = t_sink_sec / YR_TO_SEC
            v_wall_mm_yr = v_wall * YR_TO_SEC * 1e3
            print(f"{world_name:>12} {dz:14.1e} {v_wall:14.3e} {v_wall_mm_yr:16.3e} {t_sink_yr:16.2e}")

    print("\n--- RATIO: Racemization Half-life / Layer Sinking Time (tau_rac / t_sink) ---")
    print("(Values << 1 mean racemization COMPLETELY ERASES chiral excess before layers can sink!)")
    print(f"{'World':>12} {'D_z (m^2/s)':>12} {'Alanine 0°C':>16} {'Alanine 25°C':>16} {'Aspartic 0°C':>16} {'Aspartic 25°C':>16}")
    for world_name, data in WORLDS.items():
        for dz in [1e-7, 1e-5, 1e-3]:
            t_sink_yr = layer_sinking_time(data["H"], data["R"], dz) / YR_TO_SEC
            
            ala_0 = (t_half_racemization(k_racemization(AMINO_ACIDS["Alanine (alpha-H)"]["Ea"], AMINO_ACIDS["Alanine (alpha-H)"]["A"], 0.0)) / YR_TO_SEC) / t_sink_yr
            ala_25 = (t_half_racemization(k_racemization(AMINO_ACIDS["Alanine (alpha-H)"]["Ea"], AMINO_ACIDS["Alanine (alpha-H)"]["A"], 25.0)) / YR_TO_SEC) / t_sink_yr
            asp_0 = (t_half_racemization(k_racemization(AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["Ea"], AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["A"], 0.0)) / YR_TO_SEC) / t_sink_yr
            asp_25 = (t_half_racemization(k_racemization(AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["Ea"], AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["A"], 25.0)) / YR_TO_SEC) / t_sink_yr

            print(f"{world_name:>12} {dz:12.1e} {ala_0:16.2e} {ala_25:16.2e} {asp_0:16.2e} {asp_25:16.2e}")


def run_hydrothermal_transit_analysis():
    print("\n" + "=" * 80)
    print("3. HYDROTHERMAL VENT PASSAGE ERASURE (100 °C - 300 °C)")
    print("=" * 80)
    print("When ocean water circulates through hydrothermal cores:")
    
    vent_temps = [100.0, 150.0, 200.0, 300.0]
    transit_times_days = [0.1, 1.0, 10.0, 100.0] # 2.4 hours to ~3 months passage through hot zone

    print(f"\nRemaining chiral excess e/e0 after passing through hot vent zone:")
    print(f"{'Vent Temp (°C)':>16} {'Transit Time':>16} {'Alanine e/e0':>18} {'Aspartic e/e0':>18} {'Isovaline e/e0':>18}")
    for temp in vent_temps:
        for t_days in [1.0, 10.0]:
            t_sec = t_days * 86400.0
            
            kr_ala = k_racemization(AMINO_ACIDS["Alanine (alpha-H)"]["Ea"], AMINO_ACIDS["Alanine (alpha-H)"]["A"], temp)
            kr_asp = k_racemization(AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["Ea"], AMINO_ACIDS["Aspartic acid (fast alpha-H)"]["A"], temp)
            kr_iso = k_racemization(AMINO_ACIDS["Isovaline (alpha-methyl)"]["Ea"], AMINO_ACIDS["Isovaline (alpha-methyl)"]["A"], temp)
            
            surv_ala = np.exp(-2.0 * kr_ala * t_sec)
            surv_asp = np.exp(-2.0 * kr_asp * t_sec)
            surv_iso = np.exp(-2.0 * kr_iso * t_sec)
            
            print(f"{temp:16.1f} {f'{t_days:.0f} days':>16} {surv_ala:18.2e} {surv_asp:18.2e} {surv_iso:18.4f}")


def run_fluid_dynamics_analysis():
    print("\n" + "=" * 80)
    print("4. FLUID DYNAMICS: CURVATURE WALL VELOCITY VS. TURBULENT EDDY SHEAR")
    print("=" * 80)
    
    print(f"\n{'World':>12} {'D_z (m^2/s)':>12} {'v_wall (m/s)':>14} {'u_eddy (m/s)':>14} {'Vel Ratio Pi_v':>18} {'Shear Ratio S':>16}")
    
    for world_name, data in WORLDS.items():
        print(f"--- {world_name} (u_eddy = {data['u_eddy']*1e3:.1f} mm/s, L_eddy = {data['L_eddy']/1e3:.0f} km) ---")
        gamma_eddy = data["u_eddy"] / data["L_eddy"] # s^-1
        
        for dz in [1e-7, 1e-5, 1e-3]:
            v_wall = curvature_wall_velocity(dz, data["R"])
            gamma_curv = v_wall / data["H"] # characteristic wall strain rate = 2*Dz / (H*R)
            
            pi_v = v_wall / data["u_eddy"]
            shear_ratio = gamma_eddy / gamma_curv # = u_eddy * H * R / (2 * Dz * L_eddy)
            
            print(f"{world_name:>12} {dz:12.1e} {v_wall:14.3e} {data['u_eddy']:14.3e} {pi_v:18.3e} {shear_ratio:16.2e}")

    print("\n--- DAMKÖHLER NUMBER (Da = k2*c / gamma_eddy) ---")
    print("Compares chemical autocatalytic phase-separation rate to turbulent shear strain rate:")
    print("If Da << 1, turbulent shear tears apart domains faster than autocatalysis can form them!")
    
    concentrations = [1e-9, 1e-6, 1e-3] # 1 nM, 1 uM, 1 mM
    k2_values = [1e-6, 1e-3, 1.0]       # /M/s from project range
    
    for world_name, data in WORLDS.items():
        gamma_eddy = data["u_eddy"] / data["L_eddy"]
        print(f"\n{world_name} (eddy shear rate gamma = {gamma_eddy:.2e} s^-1):")
        print(f"{'Conc c (M)':>12} {'k2 (/M/s)':>12} {'Rate k2*c (s^-1)':>18} {'Da = k2*c / gamma':>20} {'Regime':>24}")
        for c in concentrations:
            for k2 in k2_values:
                rate = k2 * c
                da = rate / gamma_eddy
                if da < 1e-2:
                    regime = "Turbulence shreds"
                elif da < 1.0:
                    regime = "Sub-critical / unformed"
                else:
                    regime = "Autocatalysis survives"
                print(f"{c:12.1e} {k2:12.1e} {rate:18.2e} {da:20.2e} {regime:>24}")


def main():
    run_racemization_analysis()
    run_sinking_vs_racemization_analysis()
    run_hydrothermal_transit_analysis()
    run_fluid_dynamics_analysis()


if __name__ == "__main__":
    main()
