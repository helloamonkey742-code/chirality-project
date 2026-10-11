"""Titan at 94 K: replace titan_95k.py's hand-set rate with Longo & Remsing's barriers. (2026-10-11)

titan_95k.py needs k2 >= 1.0e-9 /M/s (3.2 uM, 100 Myr) and hand-sets k2 = 1e-8 /M/s. This script asks what reaction
barrier that bar means at 94 K, and where the one published 94 K calculation sits against it.

Source: Longo & Remsing 2024, arXiv:2405.19195v2, "Quantum Tunneling Could Enable Proton Transfer Reactions on
Titan" (full text read 2026-10-11). Path-integral MD of a model proton-transfer complex (modified Azzouz-Borgis,
donor-acceptor distance fixed at 2.7 A) in liquid ethane at 94 K. Free-energy barrier: classical 18.0 kcal/mol
(~97 kT), quantum 2.8 kcal/mol (~15 kT); by transition-state theory (TST) ln(kq/kc) ~ 82, i.e. ~35 orders of
magnitude. They give barriers and this ratio, not rate constants; they note TST neglects recrossing and that the
model is generic, not a specific Titan molecule.
Viscosity of the saturated liquid at 94.4 K (NIST WebBook, Lemmon/Friend correlations): ethane 1.098 mPa s,
methane 1.73e-4 Pa s at 94.7 K.

Method (standard TST, 1 M standard state): k = (kT/h) exp(-dG/RT). For a bimolecular step dG is the effective
activation free energy including association. Encounter cap (Smoluchowski, equal-size spheres): k_d = 8RT/(3 eta).
A step inside an already-formed complex runs at (kT/h) exp(-dG/RT) per second, so a bimolecular step can be no
faster than k_d.

Not modelled: whether the steps of a real Frank network at 94 K are proton transfers. Autocatalytic bond forming
(C-C, C-N) moves heavy atoms, which tunnel far less; Longo & Remsing's barrier does not apply to them. Solubility
in hydrocarbons is also not addressed (they flag it too). This is a threshold, not a measured Titan rate.

PREDICTION (written 2026-10-11, before the first run):
P1. TST at 94 K reproduces Longo & Remsing: barriers of 15.0 and 96.4 kT (quantum, classical) and an enhancement
    of 10^35.6 (within 1 kT and 0.5 decades of their ~15, ~97 kT and 10^35).
P2. The model's bar k2 >= 1.0e-9 /M/s (100 Myr) equals an effective barrier of <= 9.2 kcal/mol at 94 K; the 4 Gyr
    bar 2.5e-11 /M/s equals <= 9.8 kcal/mol. Both lie between the quantum 2.8 and classical 18.0 kcal/mol, so
    tunnelling, not temperature alone, decides whether a proton-transfer-limited network can run.
P3. A step with the quantum 2.8 kcal/mol barrier runs at ~6e5 /s inside a complex; a bimolecular version is then
    capped by encounter at k_d = 1.9e9 /M/s (ethane) to 1.2e10 /M/s (methane), clearing the 100 Myr bar by >= 18
    decades. titan_95k.py's hand-set 1e-8 is conservative by that measure.
P4. With the classical 18.0 kcal/mol, k2 ~ 3e-30 /M/s: >= 20 decades short of the bar.
P5. titan_95k.py's assumed racemization half-life 2.5e6 yr equals a unimolecular barrier of 11.3 kcal/mol at 94 K.
    The same tunnelling that helps cuts both ways: an alpha-C-H proton abstraction with a 2.8 kcal/mol barrier
    would racemize in ~1 us. So titan_95k.py's numbers assume bond forming <= 9.2 and racemization ~11.3 kcal/mol
    (effective, quantum-corrected): a gap of ~2 kcal/mol between two unmeasured barriers.
"""
import numpy as np

R, KB, H, NA = 8.314462618, 1.380649e-23, 6.62607015e-34, 6.02214076e23
KCAL = 4184.0
YR = 3.15576e7
T = 94.0
DG_Q, DG_C = 2.8, 18.0                  # kcal/mol, Longo & Remsing 2024
ETA = {"ethane": 1.098e-3, "methane": 1.7259e-4}  # Pa s, NIST saturated liquid ~94.5 K
C, TAU = 3.2e-6, {"100 Myr": 1e8 * YR, "4 Gyr": 4e9 * YR}
T_HALF_RAC = 2.5e6 * YR                 # titan_95k.py P5, ASSUMED

nu = KB * T / H                         # TST attempt frequency (/s)
rate = lambda dg: nu * np.exp(-dg * KCAL / (R * T))
barrier = lambda k: R * T * np.log(nu / k) / KCAL  # kcal/mol giving rate k


def main():
    kT = R * T / KCAL
    print(f"T = {T} K, kT = {kT:.4f} kcal/mol, kT/h = {nu:.3e} /s")

    print("\nP1. Longo & Remsing barriers in kT, and the TST enhancement")
    q, c = DG_Q / kT, DG_C / kT
    print(f"  quantum {DG_Q} kcal/mol = {q:.1f} kT; classical {DG_C} kcal/mol = {c:.1f} kT; "
          f"ln(kq/kc) = {c - q:.1f}, 10^{(c - q) / np.log(10):.1f}")
    assert abs(q - 15) < 1 and abs(c - 97) < 1 and abs((c - q) / np.log(10) - 35) < 1

    print("\nP2. The model's bar k2*c*tau >= 10 at c = 3.2 uM, as an effective barrier at 94 K")
    bars = {}
    for name, tau in TAU.items():
        k = 10 / (C * tau)
        bars[name] = barrier(k)
        print(f"  {name:8}: k2 >= {k:.2e} /M/s  <=>  dG <= {bars[name]:.2f} kcal/mol ({bars[name] / kT:.1f} kT)")
    assert DG_Q < min(bars.values()) and max(bars.values()) < DG_C

    print("\nP3. Quantum barrier: rate inside a complex, and the encounter cap in Titan liquids")
    k_in = rate(DG_Q)
    print(f"  in-complex step: {k_in:.2e} /s (half-life {np.log(2) / k_in * 1e6:.2f} us)")
    kd = {s: 8 * R * T / (3 * e) * 1e3 for s, e in ETA.items()}  # m^3/mol/s -> /M/s
    for s, k in kd.items():
        print(f"  encounter cap, {s:7}: k_d = {k:.2e} /M/s, {np.log10(k / (10 / (C * TAU['100 Myr']))):.1f} "
              f"decades above the 100 Myr bar")
    assert all(k > 1e9 for k in kd.values())

    print("\nP4. Classical barrier as a bimolecular TST step")
    kc = rate(DG_C)
    print(f"  k2 = {kc:.1e} /M/s, {np.log10(10 / (C * TAU['100 Myr']) / kc):.1f} decades short of the 100 Myr bar")
    assert kc < 1e-9

    print("\nP5. Racemization: the barrier behind titan_95k.py's assumed half-life")
    rac = barrier(np.log(2) / T_HALF_RAC)
    print(f"  half-life 2.5e6 yr  <=>  dG = {rac:.2f} kcal/mol; a 2.8 kcal/mol abstraction would take "
          f"{np.log(2) / k_in * 1e6:.2f} us")
    print(f"  window: bond forming <= {bars['100 Myr']:.2f}, racemization >= {rac:.2f} kcal/mol "
          f"(width {rac - bars['100 Myr']:.2f} kcal/mol)")
    assert bars["100 Myr"] < rac


if __name__ == "__main__":
    main()
