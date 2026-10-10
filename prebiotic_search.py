"""Systematic Computational Search for Prebiotic Chiral Amplification Candidates.

Evaluates two distinct regimes:
1. The Single-Molecule Regime:
   Systematic multi-criterion scoring of known single prebiotic chemical intermediates
   (Ribo-aminooxazoline, Proline, Diketopiperazines, Thiazolidine-4-carboxylic acid,
   Glyceraldehyde) across synthesis plausibility, aqueous stability, and phase/kinetic
   symmetry-breaking mechanisms.

2. The Minimal Prebiotic Network Regime (Combinatorial Peptide Search):
   Physical organic chemistry screening of all 100 dipeptides and 1,000 tripeptides
   derived from the 10 prebiotic amino acids (Gly, Ala, Val, Leu, Ile, Pro, Asp, Glu, Ser, Thr).
   Uses Taft steric parameters (E_s), conformational rigidity, and transition-state steric
   clash models calibrated against Yu et al. 2024 (PNAS) to compute stereoselectivity S,
   cross-referencing against Murchison meteorite prebiotic abundance data to identify
   the Pareto-optimal candidates crossing S_crit >= 2.4.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: SINGLE-MOLECULE PREBIOTIC CANDIDATE AUDIT
# ==============================================================================

SINGLE_MOLECULE_CANDIDATES = [
    {
        "name": "Ribo-aminooxazoline (RAO)",
        "formula": "C4 H6 N2 O2",
        "prebiotic_steps": 1,  # One-pot: glycolaldehyde + glyceraldehyde + cyanamide + phosphate
        "hydrolysis_t12_days_25C": 1e4,  # Highly stable in water up to 100 °C
        "mechanism": "Conglomerate crystallization + Viedma ripening / CISS",
        "phase": "Solid-liquid phase equilibrium (slurry)",
        "amplification_to_100_ee": True,
        "homogeneous_fluid_autocatalytic": False,  # Does not autocatalyze own synthesis in clear solution
        "experimental_ref": "Hein et al. 2011 Nat Chem; Powner et al. 2009 Nature; Ozturk & Sasselov 2023",
        "notes": "Forms chiral conglomerate crystals. Under abrasive grinding or magnetic CISS, deracemizes to 100% ee in water."
    },
    {
        "name": "Cyclo(L-Phe-L-His) [Diketopiperazine / DKP]",
        "formula": "C15 H16 N4 O2",
        "prebiotic_steps": 2,  # Dipeptide condensation -> spontaneous cyclization
        "hydrolysis_t12_days_25C": 1e3,  # Cyclic amide resists hydrolysis
        "mechanism": "Asymmetric organocatalysis (cyanohydrin addition)",
        "phase": "Homogeneous solution",
        "amplification_to_100_ee": False,  # Catalyzes cyanohydrins, but does NOT catalyze own synthesis
        "homogeneous_fluid_autocatalytic": False,
        "experimental_ref": "Inoue et al. 1981; Danger et al. 2012 Chem Eur J",
        "notes": "Excellent stereocontrol for secondary reactions (ee up to 90%), but lacks autocatalytic closure."
    },
    {
        "name": "L-Proline",
        "formula": "C5 H9 N O2",
        "prebiotic_steps": 2,  # Strecker or cyanosulfidic synthesis
        "hydrolysis_t12_days_25C": 1e5,  # Free monomer is stable
        "mechanism": "Organocatalytic aldol / enamine catalysis",
        "phase": "Homogeneous solution",
        "amplification_to_100_ee": False,  # Does not autocatalyze proline synthesis in water
        "homogeneous_fluid_autocatalytic": False,
        "experimental_ref": "List, Lerner, Barbas 2000 JACS; Breslow 2010 PNAS",
        "notes": "Prototypical organocatalyst for sugars, but cannot self-replicate without biological enzyme machinery."
    },
    {
        "name": "Thiazolidine-4-carboxylic acid (TCA)",
        "formula": "C4 H7 N O2 S",
        "prebiotic_steps": 1,  # Spontaneous condensation: Cysteine + Formaldehyde
        "hydrolysis_t12_days_25C": 500,  # Reversible under strong acid/base
        "mechanism": "Spontaneous resolution upon crystallization",
        "phase": "Solid-liquid phase equilibrium",
        "amplification_to_100_ee": False,  # Racemic conglomerate in specific salts only
        "homogeneous_fluid_autocatalytic": False,
        "experimental_ref": "Givens & Schrier 1978; Danger et al. 2013",
        "notes": "Prebiotically abundant from sulfur chemistry, but lacks self-amplifying feedback loop."
    },
    {
        "name": "D-Glyceraldehyde",
        "formula": "C3 H6 O3",
        "prebiotic_steps": 1,  # Formose reaction dimer of glycolaldehyde
        "hydrolysis_t12_days_25C": 30,  # Readily isomerizes and degrades via Lobry de Bruyn-van Ekenstein
        "mechanism": "Formose autocatalysis (weak stereocontrol)",
        "phase": "Homogeneous solution",
        "amplification_to_100_ee": False,  # Racemizes and forms complex mixtures (asphalt)
        "homogeneous_fluid_autocatalytic": True,  # Aldol autocatalysis exists, but unselective
        "experimental_ref": "Breslow 1959; Weber 2005 Orig Life Evol Biosph",
        "notes": "Autocatalytic in formose reaction, but suffers from low selectivity and rapid thermal degradation."
    }
]

# ==============================================================================
# PART 2: MINIMAL PREBIOTIC PEPTIDE NETWORK SEARCH
# ==============================================================================

# 10 Earliest Prebiotic Amino Acids (Murchison meteorite & Miller spark discharge)
# Abundance: relative to Glycine = 1.0 (Parker et al. 2011, Cronin & Chang 1993)
# Taft steric parameter E_s: negative values indicate larger steric bulk near alpha/beta carbon
# Rigid: 1 if conformationally restricted ring/branch, 0 otherwise
AMINO_ACIDS = {
    "Gly": {"name": "Glycine",       "abund": 1.000, "Es": 0.00,  "beta_branch": 0, "ring": 0},
    "Ala": {"name": "Alanine",       "abund": 0.650, "Es": -1.24, "beta_branch": 0, "ring": 0},
    "Val": {"name": "Valine",        "abund": 0.080, "Es": -2.66, "beta_branch": 1, "ring": 0},
    "Leu": {"name": "Leucine",       "abund": 0.040, "Es": -2.17, "beta_branch": 0, "ring": 0},
    "Ile": {"name": "Isoleucine",    "abund": 0.030, "Es": -2.90, "beta_branch": 1, "ring": 0},
    "Pro": {"name": "Proline",       "abund": 0.020, "Es": -2.85, "beta_branch": 0, "ring": 1},
    "Asp": {"name": "Aspartate",     "abund": 0.070, "Es": -1.60, "beta_branch": 0, "ring": 0},
    "Glu": {"name": "Glutamate",     "abund": 0.060, "Es": -1.65, "beta_branch": 0, "ring": 0},
    "Ser": {"name": "Serine",        "abund": 0.015, "Es": -1.50, "beta_branch": 0, "ring": 0},
    "Thr": {"name": "Threonine",     "abund": 0.008, "Es": -2.40, "beta_branch": 1, "ring": 0},
}

# Empirical calibration factors to match experimental dipeptide stereoselectivity:
# Yu et al. 2024 (PNAS):
# Ala-Ala: ee ~ 11% -> S = (1+ee)/(1-ee) = 1.11 / 0.89 = 1.25
# Val-Val: ee ~ 33% -> S = 1.33 / 0.67 = 1.99 ~ 2.0
# Leu-Phe: ee ~ 37% -> S = 1.37 / 0.63 = 2.17 ~ 2.2
# Gly-Gly: S = 1.0 (achiral)
# S_crit for spontaneous bifurcation in CSTR is ~ 2.40 at 5 mM feed, 200 h residence.

R_GAS = 8.314  # J/mol/K
TEMP = 298.15  # 25 °C in Kelvin

ALPHA_ES = 0.42     # steric coefficient (kJ/mol per |Es|)
BETA_BRANCH = 0.55  # extra barrier from beta-carbon branching (kJ/mol)
GAMMA_RING = 0.85   # conformational locking barrier from proline ring (kJ/mol)
POSITION_WEIGHTS_DI = [0.60, 0.40]              # N-terminus vs C-terminus influence
POSITION_WEIGHTS_TRI = [0.45, 0.35, 0.20]       # N-terminus, middle, C-terminus


def compute_dipeptide_properties(aa1, aa2):
    """Computes stereospecificity S and prebiotic abundance for dipeptide aa1-aa2."""
    p1 = AMINO_ACIDS[aa1]
    p2 = AMINO_ACIDS[aa2]

    # Prebiotic joint abundance
    abund = p1["abund"] * p2["abund"]

    # If both are Gly, purely achiral
    if aa1 == "Gly" and aa2 == "Gly":
        return {"seq": "Gly-Gly", "S": 1.00, "ee_pct": 0.0, "abund": abund, "supercritical": False}

    # Transition state free energy differential ddG_ddagger (kJ/mol)
    ddG = (
        POSITION_WEIGHTS_DI[0] * (ALPHA_ES * abs(p1["Es"]) + BETA_BRANCH * p1["beta_branch"] + GAMMA_RING * p1["ring"]) +
        POSITION_WEIGHTS_DI[1] * (ALPHA_ES * abs(p2["Es"]) + BETA_BRANCH * p2["beta_branch"] + GAMMA_RING * p2["ring"])
    )

    # Stereoselectivity S = exp(ddG / RT)
    S = float(np.exp((ddG * 1000.0) / (R_GAS * TEMP)))
    ee_pct = float(100.0 * (S - 1.0) / (S + 1.0))
    supercritical = S >= 2.40

    return {
        "seq": f"{aa1}-{aa2}",
        "S": S,
        "ee_pct": ee_pct,
        "abund": abund,
        "supercritical": supercritical,
        "ddG_kJ": ddG
    }


def compute_tripeptide_properties(aa1, aa2, aa3):
    """Computes stereospecificity S and prebiotic abundance for tripeptide aa1-aa2-aa3."""
    p1 = AMINO_ACIDS[aa1]
    p2 = AMINO_ACIDS[aa2]
    p3 = AMINO_ACIDS[aa3]

    abund = p1["abund"] * p2["abund"] * p3["abund"]

    if aa1 == "Gly" and aa2 == "Gly" and aa3 == "Gly":
        return {"seq": "Gly-Gly-Gly", "S": 1.00, "ee_pct": 0.0, "abund": abund, "supercritical": False}

    # Tripeptides fold into compact secondary turns, enhancing chiral pocket discrimination
    # Synergy bonus if proline provides a rigid turn
    turn_bonus = 0.40 if (p1["ring"] or p2["ring"]) else 0.0

    ddG = (
        POSITION_WEIGHTS_TRI[0] * (ALPHA_ES * abs(p1["Es"]) + BETA_BRANCH * p1["beta_branch"] + GAMMA_RING * p1["ring"]) +
        POSITION_WEIGHTS_TRI[1] * (ALPHA_ES * abs(p2["Es"]) + BETA_BRANCH * p2["beta_branch"] + GAMMA_RING * p2["ring"]) +
        POSITION_WEIGHTS_TRI[2] * (ALPHA_ES * abs(p3["Es"]) + BETA_BRANCH * p3["beta_branch"] + GAMMA_RING * p3["ring"]) +
        turn_bonus
    )

    S = float(np.exp((ddG * 1000.0) / (R_GAS * TEMP)))
    ee_pct = float(100.0 * (S - 1.0) / (S + 1.0))
    supercritical = S >= 2.40

    return {
        "seq": f"{aa1}-{aa2}-{aa3}",
        "S": S,
        "ee_pct": ee_pct,
        "abund": abund,
        "supercritical": supercritical,
        "ddG_kJ": ddG
    }


def run_full_peptide_search():
    """Screens all 100 dipeptides and 1,000 tripeptides."""
    aa_keys = list(AMINO_ACIDS.keys())

    # 1. Screen Dipeptides
    dipeptides = []
    for a1 in aa_keys:
        for a2 in aa_keys:
            dipeptides.append(compute_dipeptide_properties(a1, a2))

    # 2. Screen Tripeptides
    tripeptides = []
    for a1 in aa_keys:
        for a2 in aa_keys:
            for a3 in aa_keys:
                tripeptides.append(compute_tripeptide_properties(a1, a2, a3))

    return dipeptides, tripeptides


def generate_search_visualization(dipeptides, tripeptides):
    """Plots Pareto frontier: Prebiotic Abundance vs Catalytic Stereospecificity S."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # ------------------ PANEL 1: Dipeptide Landscape ------------------
    s_di = [d["S"] for d in dipeptides]
    ab_di = [d["abund"] for d in dipeptides]
    labels_di = [d["seq"] for d in dipeptides]

    ax1.scatter(s_di, ab_di, color="navy", alpha=0.6, edgecolors="none", s=40, label="Dipeptides (100 total)")

    # Threshold line
    ax1.axvline(2.40, color="crimson", linestyle="--", lw=2, label=r"Bifurcation Threshold ($S_{\mathrm{crit}} = 2.40$)")

    # Annotate key benchmarks
    for d in dipeptides:
        if d["seq"] in ["Val-Val", "Leu-Val", "Pro-Val", "Ile-Ile", "Ala-Ala", "Gly-Gly", "Pro-Ile"]:
            offset_y = 1.3 if d["seq"] in ["Ala-Ala", "Gly-Gly"] else 0.7
            ax1.scatter([d["S"]], [d["abund"]], color="red" if d["supercritical"] else "gold",
                        edgecolors="black", s=80, zorder=5)
            ax1.annotate(f"{d['seq']}\n(S={d['S']:.2f})", xy=(d["S"], d["abund"]),
                         xytext=(d["S"] + 0.05, d["abund"] * offset_y),
                         fontsize=8, weight="bold")

    ax1.set_yscale("log")
    ax1.set_xlabel(r"Catalytic Stereospecificity $S = k_{\mathrm{homo}} / k_{\mathrm{cross}}$", fontsize=11)
    ax1.set_ylabel("Relative Prebiotic Abundance (Murchison/Miller proxy)", fontsize=11)
    ax1.set_title("A. Prebiotic Dipeptides: 100 Combinations", fontsize=12, weight="bold")
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(loc="lower left", fontsize=8.5)

    # ------------------ PANEL 2: Tripeptide Landscape ------------------
    s_tri = [t["S"] for t in tripeptides]
    ab_tri = [t["abund"] for t in tripeptides]

    ax2.scatter(s_tri, ab_tri, color="teal", alpha=0.4, edgecolors="none", s=25, label="Tripeptides (1,000 total)")
    ax2.axvline(2.40, color="crimson", linestyle="--", lw=2, label=r"Bifurcation Threshold ($S_{\mathrm{crit}} = 2.40$)")

    # Annotate top candidates crossing the boundary with highest abundance
    super_tri = [t for t in tripeptides if t["supercritical"]]
    super_tri.sort(key=lambda x: x["abund"], reverse=True)

    for t in super_tri[:5]:
        ax2.scatter([t["S"]], [t["abund"]], color="lime", edgecolors="black", s=75, zorder=5)
        ax2.annotate(f"{t['seq']}\n(S={t['S']:.2f})", xy=(t["S"], t["abund"]),
                     xytext=(t["S"] + 0.08, t["abund"] * 0.75),
                     fontsize=8, weight="bold")

    ax2.set_yscale("log")
    ax2.set_xlabel(r"Catalytic Stereospecificity $S = k_{\mathrm{homo}} / k_{\mathrm{cross}}$", fontsize=11)
    ax2.set_ylabel("Relative Prebiotic Abundance (Murchison/Miller proxy)", fontsize=11)
    ax2.set_title("B. Prebiotic Tripeptides: 1,000 Combinations", fontsize=12, weight="bold")
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(loc="lower left", fontsize=8.5)

    fig.tight_layout()
    plot_path = "prebiotic_candidates.png"
    fig.savefig(plot_path, dpi=160)
    print(f"Saved candidate diagnostic figure to: {plot_path}")


def main():
    print("=" * 80)
    print("PREBIOTIC CHIRAL AMPLIFIER SYSTEMATIC SEARCH")
    print("=" * 80)

    # 1. Single-Molecule Candidate Summary
    print("\n--- 1. SINGLE-MOLECULE PREBIOTIC AUDIT ---")
    for cand in SINGLE_MOLECULE_CANDIDATES:
        print(f"\n* Molecule: {cand['name']} [{cand['formula']}]")
        print(f"  - Prebiotic Steps: {cand['prebiotic_steps']}")
        print(f"  - Mechanism: {cand['mechanism']}")
        print(f"  - Phase Regime: {cand['phase']}")
        print(f"  - Reaches 100% ee: {cand['amplification_to_100_ee']}")
        print(f"  - Fluid-phase Autocatalytic: {cand['homogeneous_fluid_autocatalytic']}")
        print(f"  - Key Ref: {cand['experimental_ref']}")
        print(f"  - Notes: {cand['notes']}")

    # 2. Peptide Search
    print("\n--- 2. PEPTIDE COMBINATORIAL SEARCH ---")
    dipeptides, tripeptides = run_full_peptide_search()

    # Sort dipeptides by S
    dipeptides.sort(key=lambda x: x["S"], reverse=True)
    super_di = [d for d in dipeptides if d["supercritical"]]
    print(f"\nTotal Dipeptides Screened: {len(dipeptides)}")
    print(f"Dipeptides Crossing Bifurcation (S >= 2.40): {len(super_di)} / {len(dipeptides)}")
    print("Top 5 Dipeptides by Stereospecificity S:")
    for d in dipeptides[:5]:
        print(f"  {d['seq']}: S = {d['S']:.2f} (ee = {d['ee_pct']:.1f}%), Rel Abundance = {d['abund']:.2e}, Supercritical = {d['supercritical']}")

    # Top dipeptides by Pareto efficiency (Abundance * (S - 1))
    pareto_di = sorted(dipeptides, key=lambda x: x["abund"] * max(0.01, x["S"] - 1.0), reverse=True)
    print("\nTop 5 Prebiotically Viable Dipeptides (Pareto Abundance x Induction):")
    for d in pareto_di[:5]:
        print(f"  {d['seq']}: S = {d['S']:.2f} (ee = {d['ee_pct']:.1f}%), Rel Abundance = {d['abund']:.2e}, Supercritical = {d['supercritical']}")

    # Sort tripeptides by S
    tripeptides.sort(key=lambda x: x["S"], reverse=True)
    super_tri = [t for t in tripeptides if t["supercritical"]]
    print(f"\nTotal Tripeptides Screened: {len(tripeptides)}")
    print(f"Tripeptides Crossing Bifurcation (S >= 2.40): {len(super_tri)} / {len(tripeptides)} ({100.0*len(super_tri)/len(tripeptides):.1f}%)")

    # Top tripeptides by Abundance that ARE supercritical
    super_tri_by_abund = sorted(super_tri, key=lambda x: x["abund"], reverse=True)
    print("\nTop 5 Prebiotically Abundant SUPERCRITICAL Tripeptides (S >= 2.40):")
    for t in super_tri_by_abund[:5]:
        print(f"  {t['seq']}: S = {t['S']:.2f} (ee = {t['ee_pct']:.1f}%), Rel Abundance = {t['abund']:.2e}")

    generate_search_visualization(dipeptides, tripeptides)
    print("\nSearch complete.")


if __name__ == "__main__":
    main()
