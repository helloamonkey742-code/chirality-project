"""Can a lab see a reaction as slow as nature allows? Sizing the slow experiment.

Seeded assay: start with enantiomeric excess e0; a real amplifier (self-copying AND
mutual suppression) makes |ee| grow as e0*exp(k*t) with k = k2*c, sign following the seed.
Plain self-copying without suppression keeps ee constant, so growth tests R1+R2 together.

Detection: ee measured with noise sigma per sample; call it detected when the rise
exceeds 3*sigma*sqrt(2) (difference of two noisy points). Then
    k2_min = ln(1 + 3*sqrt2*sigma/e0) / (c * T)
Compare with nature's floor from scorecard.py (3e-7 /M/s pond at 10 mM, 3e-9 ocean at 1 uM;
and 3e-6 /M/s for a realistic 1 mM pond over 100 years).

Arm F parameters:
Prebiotic nucleotide precursor (ribo-aminooxazoline RAO) synthesis coupled with chiral
dipeptide catalysis, operated in a fed (flow) regime to prevent thermodynamic relaxation
to racemic equilibrium (R3/R6).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

YR = 3.156e7
SIGMA = 0.002       # ee measurement noise (0.2% absolute) -- baseline assumption
E0 = 0.05           # seed excess, 5%
FLOOR = {"pond floor (10 mM, 100 yr)": 3e-7, "ocean floor (1 uM, 1e8 yr)": 3e-9}
POND_FLOOR_1MM = 10 / (1e-3 * 100 * YR)  # 3.168e-6 /M/s for 1 mM pond over 100 yr

# Arm F specific parameters (prebiotic nucleotide precursor / RAO synthesis + dipeptide catalysis)
ARM_F = {
    "name": "Arm F: Prebiotic RAO synthesis + dipeptide catalysis (fed regime)",
    "c_realistic": 1e-3,            # 1 mM realistic meteoritic / prebiotic concentration
    "c_accelerated": [1e-2, 1e-1],   # 10 mM and 100 mM accelerated tiers
    "e0_seed": 0.20,                # 20% enantiomeric seed, ASSUMED (chosen near Murchison isovaline 18.5%)
    "sigma_per_vial": 0.026,        # ASSUMED per-vial ee noise. 0.026 is the Murchison isovaline
                                    # uncertainty (Glavin & Dworkin 2009), NOT a measured per-vial lab noise.
    "n_vials": 9,                   # 9 independent reactors per condition (Dworkin test 5)
    "flow_rate_uL_h": 5.0,          # 5 uL/h fed flow rate
    "reactor_volume_mL": 2.0,       # 2 mL reactor working volume
    "pH_range": (7.0, 8.5),         # neutral to slightly alkaline pH
    "reactants": {
        "glycolaldehyde": 1e-3,     # 1 mM realistic
        "glyceraldehyde": 1e-3,     # 1 mM realistic
        "2-aminooxazole": 1e-3,     # 1 mM realistic
        "dipeptide_catalyst": 1e-4, # 0.1 mM (10 mol%)
    },
}


def k2_min(c, T, e0=E0, sigma=SIGMA):
    """Smallest detectable rate constant k2 (/M/s) over time T (seconds)."""
    return np.log(1 + 3 * np.sqrt(2) * sigma / e0) / (c * T)


def arm_f_detectable_rate(c, T, e0=ARM_F["e0_seed"], sigma=ARM_F["sigma_per_vial"], n_vials=ARM_F["n_vials"]):
    """Smallest detectable k2 for Arm F accounting for n-vial noise reduction (sigma_eff = sigma / sqrt(n))."""
    sigma_eff = sigma / np.sqrt(n_vials)
    return k2_min(c, T, e0=e0, sigma=sigma_eff)


def arm_f_fed_parameters(flow_rate_uL_h=5.0, reactor_volume_mL=2.0):
    """Hydrodynamic parameters for Arm F's fed (CSTR) regime.
    
    Returns residence time, dilution rate, and washout threshold k2.
    """
    Q_mL_h = flow_rate_uL_h / 1000.0                # mL/h
    tau_res_hours = reactor_volume_mL / Q_mL_h       # hours
    tau_res_days = tau_res_hours / 24.0              # days
    tau_res_s = tau_res_hours * 3600.0               # seconds
    dilution_rate_s = 1.0 / tau_res_s                # 1/s
    return {
        "flow_rate_uL_h": flow_rate_uL_h,
        "reactor_volume_mL": reactor_volume_mL,
        "tau_res_hours": tau_res_hours,
        "tau_res_days": tau_res_days,
        "tau_res_s": tau_res_s,
        "dilution_rate_per_s": dilution_rate_s,
        "dilution_rate_per_day": dilution_rate_s * 86400.0,
    }


def arm_f_asymmetric_rescue_check(g_a, g_s, kappa_as, kappa_sa):
    """Check if asymmetric peptide-to-sugar transmission rescues natural homochirality.
    
    Condition from dna_network.py:
        kappa_as / kappa_sa > (|g_s| / g_a)**2
    Ensures that L-amino acids / dipeptides drive D-sugars / D-nucleotide precursors
    even when sugar PVED opposes the natural form (g_s < 0).
    """
    if kappa_sa <= 0:
        return True
    return (kappa_as / kappa_sa) > (abs(g_s) / g_a)**2


def main():
    # self-check: at k2 = k2_min the seeded ee rises by exactly the detection threshold
    c, T = 0.1, YR
    rise = E0 * (np.exp(k2_min(c, T) * c * T) - 1)
    assert abs(rise - 3 * np.sqrt(2) * SIGMA) < 1e-12

    # Arm F self-checks
    sigma_eff_f = ARM_F["sigma_per_vial"] / np.sqrt(ARM_F["n_vials"])
    k2_f_1yr = arm_f_detectable_rate(ARM_F["c_realistic"], YR)
    rise_f = ARM_F["e0_seed"] * (np.exp(k2_f_1yr * ARM_F["c_realistic"] * YR) - 1)
    assert abs(rise_f - 3 * np.sqrt(2) * sigma_eff_f) < 1e-12

    # Arm F 3-year realistic arm reaches below pond floor
    k2_f_3yr = arm_f_detectable_rate(ARM_F["c_realistic"], 3 * YR)
    assert k2_f_3yr < POND_FLOOR_1MM, f"Arm F 3-year run must beat pond floor: {k2_f_3yr} vs {POND_FLOOR_1MM}"
    assert abs(k2_f_3yr - 1.83e-6) < 0.05e-6

    # Arm F fed hydrodynamic parameters check
    fed = arm_f_fed_parameters(ARM_F["flow_rate_uL_h"], ARM_F["reactor_volume_mL"])
    assert abs(fed["tau_res_hours"] - 400.0) < 1e-6
    assert abs(fed["tau_res_days"] - 16.666666666666668) < 1e-6

    # Arm F asymmetric rescue check (matches dna_network.py)
    assert arm_f_asymmetric_rescue_check(1e-3, -1e-3, 0.9, 0.1) is True
    assert arm_f_asymmetric_rescue_check(1e-3, -1e-3, 0.1, 0.9) is False

    durations = {"1 month": YR / 12, "6 months": YR / 2, "1 year": YR, "3 years": 3 * YR}
    print(f"=== Baseline Detection Limits (seed ee {E0:.0%}, noise {SIGMA:.1%}) ===")
    print(f"{'conc.':>8} " + " ".join(f"{d:>10}" for d in durations))
    for c in (0.01, 0.1, 1.0):
        print(f"{c:>6} M " + " ".join(f"{k2_min(c, T):10.0e}" for T in durations.values()))
    print("nature's floor: " + ", ".join(f"{k} {v:.0e}" for k, v in FLOOR.items()))

    print("\n=== Arm F & Realistic Sizing (seed ee 20%, 9 vials, noise 2.6%/sqrt(9) = 0.87%) ===")
    print(f"{'tier / conc':>22} " + " ".join(f"{d:>10}" for d in durations))
    # 1 mM realistic
    row_1mM = [arm_f_detectable_rate(1e-3, T) for T in durations.values()]
    print(f"{'1 mM (realistic)':>22} " + " ".join(f"{v:10.0e}" for v in row_1mM))
    # 10 mM accelerated
    row_10mM = [arm_f_detectable_rate(1e-2, T) for T in durations.values()]
    print(f"{'10 mM (accelerated)':>22} " + " ".join(f"{v:10.0e}" for v in row_10mM))
    # 0.1 M accelerated
    row_100mM = [arm_f_detectable_rate(1e-1, T) for T in durations.values()]
    print(f"{'0.1 M (accelerated)':>22} " + " ".join(f"{v:10.0e}" for v in row_100mM))
    print(f"nature's floor (1 mM pond, 100 yr): {POND_FLOOR_1MM:.1e} /M/s")

    print("\n=== Arm F Fed Regime Hydrodynamics ===")
    print(f"Flow rate: {fed['flow_rate_uL_h']} uL/h | Reactor volume: {fed['reactor_volume_mL']} mL")
    print(f"Residence time: {fed['tau_res_hours']:.1f} h ({fed['tau_res_days']:.1f} days)")
    print(f"Dilution rate D: {fed['dilution_rate_per_s']:.2e} /s ({fed['dilution_rate_per_day']:.3f} /day)")
    print(f"Washout threshold k2*c > D -> at 1 mM, k2 > {fed['dilution_rate_per_s'] / 1e-3:.2e} /M/s")

    Ts = np.logspace(np.log10(YR / 12), np.log10(5 * YR), 100)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    # Arm F realistic curve (1 mM, 20% seed, 9 vials)
    ax.loglog(Ts / YR, arm_f_detectable_rate(1e-3, Ts), "b-", lw=2, label="1 mM (Arm F realistic, 20% seed, 9 vials)")
    # Accelerated curves
    ax.loglog(Ts / YR, arm_f_detectable_rate(1e-2, Ts), "g--", label="10 mM (accelerated, 20% seed)")
    ax.loglog(Ts / YR, arm_f_detectable_rate(1e-1, Ts), "m--", label="0.1 M (accelerated, 20% seed)")
    # Baseline comparison
    ax.loglog(Ts / YR, k2_min(1.0, Ts), "c:", label="1.0 M (baseline 5% seed, 0.2% noise)")

    # Nature's floors
    ax.axhline(POND_FLOOR_1MM, color="red", ls=":", lw=1.2)
    ax.text(0.09, POND_FLOOR_1MM * 1.2, "pond floor (1 mM, 100 yr: 3e-6)", fontsize=8, color="red")
    ax.axhline(FLOOR["pond floor (10 mM, 100 yr)"], color="0.5", ls="--", lw=0.8)
    ax.text(0.09, FLOOR["pond floor (10 mM, 100 yr)"] * 1.2, "pond floor (10 mM, 100 yr: 3e-7)", fontsize=8, color="0.35")
    ax.axhline(FLOOR["ocean floor (1 uM, 1e8 yr)"], color="0.5", ls="--", lw=0.8)
    ax.text(0.09, FLOOR["ocean floor (1 uM, 1e8 yr)"] * 1.3, "ocean floor (1 uM, 10^8 yr: 3e-9)", fontsize=8, color="0.35")

    ax.set(xlabel="experiment length (years)", ylabel="smallest detectable k2 (/M/s)",
           title="Experimental reach vs nature's floors (including Arm F)")
    ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig("design.png", dpi=150)
    print("wrote design.png")


if __name__ == "__main__":
    main()
