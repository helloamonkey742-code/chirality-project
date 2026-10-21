"""Arm E at measured rates: replaces the ASSUMED catalyst point of chiral_search.py.

Sources (every number below is CITED with table/figure and units, DERIVED from cited numbers,
or ASSUMED and scanned over a stated range):
  Yu, Daru, Deng & Blackmond 2024, PNAS 121 e2315447121 (open access; main text + SI)  [Yu]
  Deng, Yu & Blackmond 2024, Nature 626 1019 (closed access; PubMed 38418914 abstract only) [Deng]
  Higgs & Blackmond 2025, PNAS 122 e2423683122 (open access; main text + SI)  [Higgs]

Answers, written to arm_e_rates.log:
  (1) which parameters are measured and which remain assumed
  (2) at a 1-5 mM feed, is the fed network past the bifurcation; if not, by what factor must the
      deciding rate rise
  (3) R9: as the flow (feed) rate rises slowly, supercritical pitchfork / slow tip-over / jump
  (4) k2 against EXPERIMENT.md's realistic-arm floor 1.8e-6 /M/s (scorecard R5)

Units: M, s.  Run: uv run ... python arm_e_rates.py
"""
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from chiral_search import PrebioticNetworkModel, DEFAULT_K_CAT

H, mM = 3600.0, 1e-3
FLOOR = 1.8e-6  # /M/s, EXPERIMENT.md "Choice: two tiers", realistic arm (1 mM, 9 vials, 3 years)

# ------------------------------------------------------------------ CITED  [Yu]
# All Yu data: pH 8.5, room temperature, 80% H2O / 20% D2O, peptide catalyst 100 mM.
T1 = {"none": (9.2, 0.2), "Pro-Pro": (31.0, -50.0)}  # Table 1 (5 h): Ala yield %, ee % (+ = L); pyr 200, PM 100 mM
T2 = {"Pro-Pro": (23.0, 15.0), "Pro-Val": (26.0, 11.0), "Pro-Lys": (29.0, 8.4)}  # Table 2: conversion %, ee % L
#   of the remaining Ala (6 h Pro-Pro, 24 h others); rac-Ala 100 mM + PL 100 mM
R0 = {"fwd none": 1.48, "fwd Pro-Pro": 5.29, "rev none": 0.91, "rev Pro-Pro": 4.04}  # SI Table 5 initial
#   rates, mM/h; 100 mM of each reactant (pyr + PM forward, Ala + PL reverse)
S41 = [(2, 13.91, 9.55), (4, 22.06, 14.38), (6, 28.45, 16.21), (8, 32.47, 18.87)]  # SI Table S4-1, L-Pro-Pro
#   reverse: h, conversion %, ee % L
S44 = {50: 16.58, 100: 23.57, 200: 30.50, 400: 29.32}  # SI Table S4-4: 6 h conversion % vs [Pro-Pro] mM
# ------------------------------------------------------------------ CITED  [Deng]
# Abstract: "the reaction favours heterochiral ligation" -> R = k_hetero/k_homo > 1. No rate constant.
# ------------------------------------------------------------------ CITED  [Higgs]
# Destructive mechanism = SI Eq. 5a-c with dimers 4d-g (rates written as k1b, k2b; k_f = K k_b).
# Fig. 3 parameters (C0 = 1, C0 k_homo = 1, R = 4, k1b/k_homo = 10, phi/(C0 k_homo) = 0.1, K = 0.01);
# read by eye from Fig. 3: racemic [L]/C0 = 0.0545 at k2b/k1b = 0, symmetry breaks at k2b/k1b = 6.9.
FIG3 = dict(C0=1.0, kh=1.0, R=4.0, f=0.1, K=0.01, ku=0.0)
FIG3_L0, FIG3_SCRIT = 0.0545, 6.9


# ------------------------------------------------------------------ DERIVED
def kagan_s(c, ee):
    """k_fast/k_slow from conversion c and ee of the remaining substrate (fractions)."""
    return np.log((1 - c) * (1 - ee)) / np.log((1 - c) * (1 + ee))


def minus_background(s_obs, share):
    """Selectivity of the catalysed path alone; share = racemic uncatalysed part of the total rate."""
    return (s_obs / (1 + s_obs) - share / 2) / (1 / (1 + s_obs) - share / 2)


BG = R0["rev none"] / R0["rev Pro-Pro"]                       # uncatalysed share of the reverse rate
S_T2 = minus_background(kagan_s(.23, .15), BG)                 # Table 2, Pro-Pro
S_TS = [minus_background(kagan_s(c / 100, e / 100), BG) for _, c, e in S41]
S_LO, S_MID, S_HI = min(S_TS + [S_T2]), S_T2, max(S_TS + [S_T2])
y, ee = T1["Pro-Pro"]
S_FWD_RAW = (1 + ee / 100) / (1 - ee / 100)                    # k1f/k2f = L/D made by an L,L catalyst
S_FWD = (y * (1 + ee / 100) - T1["none"][0]) / (y * (1 - ee / 100) - T1["none"][0])
KB_SUM = (R0["rev Pro-Pro"] - R0["rev none"]) * mM / H / (0.100 * 0.050)  # k1b + k2b, /M/s ([PL] 100 mM)
KU = R0["rev none"] * mM / H / 0.100                                     # uncatalysed breakdown, 1/s
KF_SUM = (R0["fwd Pro-Pro"] - R0["fwd none"]) * mM / H / (0.100 * 0.100)  # k1f + k2f, /M/s ([PM] 100 mM)
K0F = R0["fwd none"] * mM / H / 0.100                                    # uncatalysed synthesis, 1/s
SAT = (S44[50] / 50) / (S44[100] / 100)          # per-catalyst rate at 50 vs 100 mM (saturation)


def kb(s, F=1.0):
    """(k1b, k2b) at selectivity s, catalysed rate scaled by F."""
    return F * KB_SUM / (1 + s), F * KB_SUM * s / (1 + s)


# ------------------------------------------------------------------ MODEL  [Higgs] Eq. 5 + 4d-g
# plus the measured uncatalysed path ku (racemic, same K). State z = [x, a, A, y, b] with
# x, a = ([L] +- [D])/2 and y, b = ([LL] +- [DD])/2; heterodimers are sinks only.
KEYS = ("C0", "f", "kh", "R", "K", "k1", "k2", "ku")


def rhs(t, z, p, f_of_t=None, g=0.0):
    x, a, A, y, b = z
    C0, f, kh, R, K, k1, k2, ku = (p[k] for k in KEYS)
    f = f_of_t(t) if f_of_t else f
    P, Q, u = (k1 + k2) * y + ku, (k1 - k2) * b, x - K * A
    return [f * (C0 / 2 - x) - (P * u + Q * a) - 2 * kh * (x * x + a * a) - 2 * R * kh * (x * x - a * a),
            f * (g * C0 / 2 - a) - (P * a + Q * u) - 4 * kh * x * a,  # g = feed bias (ramp test only)
            -f * A + 2 * (P * u + Q * a),
            -f * y + kh * (x * x + a * a),
            -f * b + 2 * kh * x * a]


def _parts(x, C0, f, kh, R, K, k1, k2, ku):
    y = kh * x * x / f
    P = (k1 + k2) * y + ku
    A = 2 * P * x / (f + 2 * P * K)
    return y, P, A, x - K * A


def steady(p):
    """Racemic steady state, vectorised bisection (the residual falls monotonically in x)."""
    v = np.broadcast_arrays(*(np.asarray(p[k], float) for k in KEYS))
    C0, f, kh, R = v[:4]
    lo, hi = np.zeros_like(C0), C0 / 2
    for _ in range(80):
        x = (lo + hi) / 2
        _, P, _, u = _parts(x, *v)
        up = f * (C0 / 2 - x) - P * u - 2 * (1 + R) * kh * x * x > 0
        lo, hi = np.where(up, x, lo), np.where(up, hi, x)
    return x, v


def _jac(p):
    """Exact (a, b) block of the Jacobian at the racemic state."""
    x, v = steady(p)
    C0, f, kh, R, K, k1, k2, ku = v
    _, P, _, u = _parts(x, *v)
    return -f - P - 4 * kh * x, (k2 - k1) * u, 2 * kh * x, -f


def lam(p):
    """Leading eigenvalue (1/s); > 0 = symmetry breaks."""
    J11, J12, J21, J22 = _jac(p)
    return (J11 + J22 + np.sqrt(np.maximum((J11 - J22) ** 2 + 4 * J12 * J21, 0))) / 2


def margin(p):
    """J12 J21 / (J11 J22): feedback over damping; > 1 exactly when lam > 0 (J11, J22 < 0)."""
    J11, J12, J21, J22 = _jac(p)
    return J12 * J21 / (J11 * J22)


def racemic_z(p):
    x, v = steady(p)
    y, _, A, _ = _parts(x, *v)
    return np.array([float(x), 0.0, float(A), float(y), 0.0])


def simulate(p, z0, T, f_of_t=None, g=0.0, n=600):
    atol = np.array([1e-13, 1e-40, 1e-13, 1e-13, 1e-40]) * p["C0"]
    s = solve_ivp(rhs, (0, T), z0, args=(p, f_of_t, g), method="LSODA", rtol=1e-9, atol=atol,
                  t_eval=np.linspace(0, T, n))
    return s.t, s.y


# ------------------------------------------------------------------ ASSUMED grid (scanned)
GRID = dict(kh=np.logspace(-9, 0, 19)[:, None, None, None],    # /M/s ligation: [Deng] rate not obtained
            R=np.array([1.0, 2.0, 4.0])[None, :, None, None],  # [Deng] says R > 1; [Higgs] Fig. S1 to 4
            K=np.array([1e-6, 1e-4, 1e-3, 1e-2, 1e-1])[None, None, :, None],  # [PM]/[PL], [Higgs] SI Eq. 3
            # f floor 1e-13 /s is too high for the ocean: its seed optimum needs f*tau >= 5, i.e. f >= 1.6e-15 /s at tau = 1e8 yr
            #   (arm_e_r6.py:58 extends the grid to 1e-19, POST HOC; arm_e_r6.log:156-174)
            f=np.logspace(-13, -3, 41)[None, None, None, :])   # 1/s flow, a design choice
C0S = (1 * mM, 2 * mM, 5 * mM)                                  # prebiotic feed (task)


def best(C0, s, F=1.0, ku=KU):
    """Over the ASSUMED grid: max margin, max lambda, and where the margin peaks."""
    k1, k2 = kb(s, F)
    p = dict(C0=C0, k1=k1, k2=k2, ku=ku, **GRID)
    M, L = margin(p), lam(p)
    i = np.unravel_index(np.argmax(M), M.shape)
    return M[i], L.max(), {k: float(GRID[k].ravel()[j]) for k, j in zip(("kh", "R", "K", "f"), i)}


def f_crit(C0, s, ku=KU, k2_min=0.0):
    """Factor on the catalysed rate (fixed s) that first gives lambda > k2_min C0 somewhere on the grid."""
    ok = lambda F: best(C0, s, F, ku)[1] > k2_min * C0
    if ok(1.0):
        return 1.0
    if not ok(1e8):
        return np.inf
    lo, hi = 0.0, 8.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ok(10 ** mid) else (mid, hi)
    return 10 ** hi


def constructive_lambda_max():
    """chiral_search.py's own (constructive) model at the measured S, scanned; returns max lambda."""
    top, bad = -np.inf, 0
    for S in (S_FWD, S_FWD_RAW):
        for F0 in C0S:
            for phi in np.logspace(-10, -3, 15):
                for kh in np.logspace(-9, 0, 10):
                    for rho in (1.0, 2.0, 4.0):
                        m = PrebioticNetworkModel(F0=F0, phi=phi, k0=K0F / 2, k_cat=KF_SUM * S / (1 + S),
                                                  S=S, k_homo=kh, rho=rho, k_rac=0.0)
                        try:
                            top = max(top, m.difference_jacobian_eigenvalue()[0])
                        except ValueError:
                            bad += 1
    return top, bad


def r9_test(p, f_lo, f_hi, g=1e-9):
    """Static branch, bistability probe, and slow up/down ramp of the flow f across f_c (lower edge)."""
    fc = np.exp(np.mean([np.log(f_lo), np.log(f_hi)]))
    lo, hi = np.log(f_lo), np.log(f_hi)
    for _ in range(60):                                   # lower critical flow
        fc = np.exp((lo + hi) / 2)
        lo, hi = (lo, np.log(fc)) if lam({**p, "f": fc}) > 0 else (np.log(fc), hi)
    fc = np.exp(hi)

    def settle(f, ee0):
        q = {**p, "f": f}
        z = racemic_z(q) * np.array([1, 0, 1, 1, 0])
        z[1], z[4] = ee0 * z[0], ee0 * z[3]
        T = 60 / max(abs(float(lam(q))), 1e-3 * f) + 400 / f
        _, Z = simulate(q, z, T, n=50)
        return Z[1, -1] / Z[0, -1]

    deltas = np.array([0.01, 0.02, 0.04, 0.08])
    branch = np.array([settle(fc * (1 + d), 1e-3) for d in deltas])
    expo = np.polyfit(np.log(deltas), np.log(np.abs(branch)), 1)[0]
    bist = settle(fc / 1.05, 0.9)                         # below f_c, started near homochiral
    f_end = fc * 1.5
    T = 2000 / fc                                         # log-ramp; d ln f/dt = ln(3)/T << fc
    up = lambda t: fc / 2 * 3 ** (t / T)
    t, Z = simulate({**p, "f": fc}, racemic_z({**p, "f": fc / 2}), T, up, g)
    ee_up, f_up = Z[1] / Z[0], up(t)
    ee_static_end = settle(f_end, 1e-3)
    i_end = np.searchsorted(f_up, f_end)
    down = lambda s: fc * 1.5 * 3 ** (-s / T)
    t2, Z2 = simulate({**p, "f": fc}, Z[:, i_end], T, down, g)
    ee_dn, f_dn = Z2[1] / Z2[0], down(t2)
    half_up = f_up[np.argmax(np.abs(ee_up) > 0.5 * abs(ee_static_end))]
    gone_dn = f_dn[np.argmax(np.abs(ee_dn) < 0.01 * abs(ee_static_end))]
    return dict(fc=fc, deltas=deltas, branch=branch, expo=expo, bist=bist, ee_ramp=ee_up[i_end],
                ee_static=ee_static_end, half_up=half_up / fc, gone_dn=gone_dn / fc, g=g)


# ------------------------------------------------------------------ self-check
def self_check():
    rng = np.random.default_rng(1)
    p = dict(C0=5 * mM, f=3e-8, kh=2e-5, R=2.0, K=0.01, k1=2e-5, k2=1.5e-4, ku=KU)
    # (i) the (x, a) form equals Higgs SI Eq. 5a-c / 4d-e written in L, D
    L, D, A, LL, DD = rng.uniform(0.2, 1.5, 5) * mM
    C0, f, kh, R, K, k1, k2, ku = (p[k] for k in KEYS)
    cL, cD = (k1 * LL + k2 * DD + ku) * (L - K * A), (k1 * DD + k2 * LL + ku) * (D - K * A)
    eq5 = [f * (C0 / 2 - L) - cL - 2 * kh * L * L - 2 * R * kh * L * D,
           f * (C0 / 2 - D) - cD - 2 * kh * D * D - 2 * R * kh * L * D,
           -f * A + cL + cD, -f * LL + kh * L * L, -f * DD + kh * D * D]
    mine = rhs(0, [(L + D) / 2, (L - D) / 2, A, (LL + DD) / 2, (LL - DD) / 2], p)
    want = [(eq5[0] + eq5[1]) / 2, (eq5[0] - eq5[1]) / 2, eq5[2], (eq5[3] + eq5[4]) / 2, (eq5[3] - eq5[4]) / 2]
    assert np.allclose(mine, want, rtol=1e-12, atol=0)
    # (ii) racemic state is a fixed point and lam() equals a finite-difference Jacobian eigenvalue
    z = racemic_z(p)
    assert np.allclose(rhs(0, z, p), 0, atol=1e-12 * p["C0"] * p["f"])
    J = np.zeros((2, 2))
    for j, idx in enumerate((1, 4)):
        dz = np.zeros(5); dz[idx] = 1e-6 * z[0]
        J[:, j] = (np.array(rhs(0, z + dz, p)) - rhs(0, z - dz, p))[[1, 4]] / (2 * dz[idx])
    assert np.isclose(max(np.linalg.eigvals(J).real), float(lam(p)), rtol=1e-5)
    # (iii) reproduces Higgs Fig. 3: racemic [L]/C0 at k2b/k1b = 0, and critical k2b/k1b ~ 6.9
    fig = lambda s: {**FIG3, "k1": 10.0, "k2": 10.0 * s}
    assert abs(float(steady(fig(0.0))[0]) - FIG3_L0) < 0.001
    ss = np.linspace(3, 12, 9001)
    s_crit = ss[np.argmax(lam({**FIG3, "k1": 10.0, "k2": 10.0 * ss}) > 0)]
    assert abs(s_crit - FIG3_SCRIT) < 0.2, s_crit
    assert lam(fig(3.0)) < 0                              # Eq. 11 needs k2b > 3 k1b
    # (iv) the sign of lambda predicts growth or decay of a seeded excess in the full nonlinear model
    for s, grows in ((s_crit + 1, True), (s_crit - 1, False)):
        q = fig(s)
        z0 = racemic_z(q); z0[1] = 1e-6 * z0[0]
        _, Z = simulate(q, z0, 400.0, n=3)               # compare t = 200 and 400 (past the transient)
        assert (abs(Z[1, -1]) > abs(Z[1, 1])) == grows
    # (v) measured data: forward and reverse selectivities agree (Higgs SI Eq. 2: k2f/k1f = k2b/k1b)
    assert 0.5 < (1 / S_FWD) / S_MID < 2
    return s_crit


# ------------------------------------------------------------------ report
def main():
    out = open(Path(__file__).with_suffix(".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        out.write(s + "\n")

    s_crit = self_check()
    log("self-check passed: Eq. 5 transcription, Jacobian, Higgs Fig. 3 (critical k2b/k1b "
        f"{s_crit:.2f} here vs {FIG3_SCRIT} read from the figure), nonlinear sign, fwd/rev consistency\n")

    log("(1) PARAMETERS  [D] = derived from cited measurements, [A] = assumed and scanned")
    rows = [
        ("k1b+k2b  L-Pro-Pro, Ala breakdown", f"{KB_SUM:.2e} /M/s", "D", "Yu SI Table 5: (4.04-0.91) mM/h at 100 mM Pro-Pro, 50+50 mM Ala, PL 100 mM"),
        ("s = k2b/k1b  catalysed", f"{S_MID:.1f} (range {S_LO:.1f}-{S_HI:.1f})", "D", f"Yu Table 2 / SI S4-1 via Kagan s; uncatalysed share {BG:.3f} removed (SI Table 5)"),
        ("  raw s before background", ", ".join(f"{kagan_s(c / 100, e / 100):.2f}" for c, e in T2.values()), "D", "Yu Table 2: Pro-Pro, Pro-Val, Pro-Lys"),
        ("ku  uncatalysed Ala breakdown", f"{KU:.2e} /s", "D", "Yu SI Table 5: 0.91 mM/h at 100 mM Ala, PL 100 mM"),
        ("k1f+k2f  L-Pro-Pro, Ala synthesis", f"{KF_SUM:.2e} /M/s", "D", "Yu SI Table 5: (5.29-1.48) mM/h at 100 mM pyr, PM, Pro-Pro"),
        ("S = k1f/k2f  (chiral_search S)", f"{S_FWD:.3f} (raw {S_FWD_RAW:.3f})", "D", "Yu Table 1: L-Pro-Pro gives 50% ee D at 31% yield; 9.2% racemic background"),
        ("  check 1/S vs s (Higgs SI Eq. 2)", f"{1 / S_FWD:.1f} vs {S_MID:.1f}", "D", "same catalyst both ways: consistent"),
        ("k0  uncatalysed synthesis", f"{K0F:.2e} /s", "D", "Yu SI Table 5: 1.48 mM/h at 100 mM pyr, PM"),
        ("rate per catalyst, 50 vs 100 mM", f"{SAT:.2f}x", "D", "Yu SI Table S4-4 (saturation; mM catalyst may be up to ~2x faster per molecule)"),
        ("R = k_hetero/k_homo", "1, 2, 4", "A", "Deng abstract: heterochiral ligation favoured (R > 1); no value obtained"),
        ("k_homo  ligation", "1e-9 - 1 /M/s", "A", "Deng 2024 rates not obtained (closed access)"),
        ("K = [Ala]/[pyr] at equilibrium", "1e-6, 1e-4, 1e-3, 0.01, 0.1", "A", "not measured (Yu); set by [PM]/[PL] (Higgs SI Eq. 3); must be < 1"),
        ("f  flow", "1e-13 - 1e-3 /s", "A", "design choice, optimised"),
        ("C0  racemic Ala feed", "1, 2, 5 mM", "A", "prebiotic feed (task)"),
        ("[PL]", "100 mM, held", "A", "as Yu; Higgs assumes PL, PM in excess; both paths taken first order in PL"),
        ("catalyst = the product's own LL", "behaves as L-Pro-Pro", "A", "Yu tested Pro dipeptides, not Ala-Ala; Higgs SI p.7: one-pot coupling not yet run"),
    ]
    for r in rows:
        log(f"  [{r[2]}] {r[0]:<36} {r[1]:<26} {r[3]}")
    log(f"  Yu's Pro-Pro makes the opposite hand going forward (S < 1): it fits Higgs's DESTRUCTIVE mechanism "
        f"(k2b > k1b, K < 1, racemic Ala fed), not chiral_search's constructive one.")
    log(f"  chiral_search's ASSUMED point (k_cat {DEFAULT_K_CAT} /M/s, S = 2) vs measured "
        f"(k1f {KF_SUM * S_FWD / (1 + S_FWD):.1e} /M/s, S {S_FWD:.2f}).\n")

    log("(2) PAST THE BIFURCATION AT 1-5 mM?")
    top, bad = constructive_lambda_max()
    log(f"  constructive (chiral_search.PrebioticNetworkModel, measured S, k_cat, k0; k_homo, rho, phi scanned): "
        f"max lambda = {top:.2e} /s (= -phi at the slowest flow; {bad} points without a root). Never past: J12 ~ (1 - 1/S) < 0, "
        f"so no rate factor helps; S would have to exceed 1 (the catalyst would have to make its own hand).")
    log("  destructive (Higgs Eq. 5 + measured uncatalysed path), best over the ASSUMED grid. margin =")
    log("  J12 J21/(J11 J22), feedback over damping of the (L-D, LL-DD) mode; past the bifurcation iff > 1:")
    log(f"  {'C0':>5} {'s':>5} {'max margin':>11} {'factor needed':>14} {'margin, no background':>22}   at kh, R, K, f")
    res = {}
    for C0 in C0S:
        for s in (S_LO, S_MID, S_HI):
            M1, _, arg = best(C0, s)
            M0 = best(C0, s, ku=0.0)[0]
            F = f_crit(C0, s)
            res[C0, s] = F
            log(f"  {C0 / mM:4.0f}mM {s:5.1f} {M1:11.4f} {F:14.0f} {M0:22.2f}"
                f"   {arg['kh']:.0e} {arg['R']:.0f} {arg['K']:.0e} {arg['f']:.0e}")
    F5 = res[5 * mM, S_MID]
    bound = KU / (5 * mM / 4 * (kb(S_MID)[1] - 3 * kb(S_MID)[0]))
    assert F5 >= bound   # Eq. 11 bound: y (k2b - 3 k1b) > ku with y <= C0/4
    log(f"  Not past at any feed or assumed setting (best at the most irreversible K and R = 1, the grid edges")
    log(f"  that help). The deciding rate is the catalysed breakdown against the uncatalysed one: at 5 mM, s = {S_MID:.1f},")
    log(f"  the catalysed rate must rise {F5:.0f}x (or the background fall as much); {res[1 * mM, S_MID]:.0f}x at 1 mM.")
    log(f"  Eq. 11 alone gives >= {bound:.0f}x (homodimer <= C0/4 = 1.25 mM vs ku/(k2b-3k1b) = "
        f"{KU / (kb(S_MID)[1] - 3 * kb(S_MID)[0]) / mM:.0f} mM). The rest:")
    q = dict(C0=5 * mM, **{k: v for k, v in best(5 * mM, S_MID, F5)[2].items()}, **dict(zip(("k1", "k2"), kb(S_MID, F5))), ku=KU)
    zq = racemic_z(q)
    log(f"  a faster catalyst also breaks down more of the fed monomer, so at the critical point free Ala is")
    log(f"  {zq[0] / (q['C0'] / 2):.1%} of the feed and the homodimer that carries the excess only {zq[3] / (q['C0'] / 4):.0%} of C0/4. Catalyst")
    log(f"  saturation (SI S4-4) can supply ~{SAT:.1f}-2x. With no uncatalysed path the measured s would suffice (margin > 1).\n")

    log("(3) R9: SLOW RISE OF THE FLOW RATE")
    log("  At measured rates lambda < 0 at every flow on the grid: the racemic state is never lost, so there is")
    log("  nothing to tip over. R9 not met at measured rates.")
    arg = best(5 * mM, S_MID, 2 * F5)[2]
    k1, k2 = kb(S_MID, 2 * F5)
    p = dict(C0=5 * mM, kh=arg["kh"], R=arg["R"], K=arg["K"], k1=k1, k2=k2, ku=KU, f=arg["f"])
    r = r9_test(p, arg["f"] / 1e4, arg["f"])
    log(f"  Hypothetical point (catalysed rate x{2 * F5:.0f}, 5 mM, s {S_MID:.1f}, kh {p['kh']:.0e}, R {p['R']:.0f}, "
        f"K {p['K']:.0e}): lower critical flow f_c = {r['fc']:.2e} /s (residence {1 / r['fc'] / 86400:.0f} d).")
    log(f"  static ee at f = f_c(1+d), d = {', '.join(f'{d:g}' for d in r['deltas'])}: "
        f"{', '.join(f'{e:.3f}' for e in r['branch'])}; ee ~ d^{r['expo']:.2f} (pitchfork: 0.5)")
    log(f"  started at ee 0.9 just below f_c (f_c/1.05): ends at ee {r['bist']:.1e} (no second stable state)")
    log(f"  slow log-ramp f_c/2 -> 1.5 f_c, feed bias g = {r['g']:g} (ASSUMED, sets only the delay): ee reaches half")
    log(f"  its branch value at f = {r['half_up']:.2f} f_c, ends at {r['ee_ramp']:.3f} vs static {r['ee_static']:.3f}; "
        f"ramping back, ee is gone by f = {r['gone_dn']:.2f} f_c")
    sup = 0.4 < r["expo"] < 0.6 and abs(r["bist"]) < 1e-3
    follows = abs(r["ee_ramp"] / r["ee_static"] - 1) < 0.05 and r["gone_dn"] > 0.8
    log(f"  -> {'supercritical pitchfork' if sup else 'NOT supercritical (jump)'}; "
        f"{'the ramp tips over slowly and follows the static branch back and forth (normal form of sim.py / frank.py)' if sup and follows else 'the ramp does not follow the branch'}."
        " R9 would hold IF the rates were high enough.\n")

    log("(4) k2 AND THE 1.8e-6 /M/s FLOOR (R5)")
    L1 = best(1 * mM, S_MID)[1]
    L0 = best(1 * mM, S_MID, ku=0.0)[1]
    ceil = (kb(S_MID)[1] - 3 * kb(S_MID)[0]) / 4
    F_R5 = f_crit(1 * mM, S_MID, k2_min=FLOOR)
    log(f"  k2 = k/c as EXPERIMENT.md (k = lambda, c = C0), 1 mM, s = {S_MID:.1f}, best over the grid:")
    log(f"    measured rates:            {L1 / mM:9.1e} /M/s  (negative: an excess decays; = -f at the slowest flow)")
    log(f"    uncatalysed path removed:  {L0 / mM:9.1e} /M/s  ({FLOOR / (L0 / mM):.0f}x below the floor)")
    log(f"    bound (k2b - 3 k1b)/4:     {ceil:9.1e} /M/s  (needs every monomer in homodimers and no losses;")
    log(f"                                                 range {(kb(S_LO)[1] - 3 * kb(S_LO)[0]) / 4:.1e}-{(kb(S_HI)[1] - 3 * kb(S_HI)[0]) / 4:.1e}; not reachable: the excess lives in")
    log(f"                                                 dimers, which turn over only at the flow rate)")
    log(f"  To reach the 1.8e-6 floor at 1 mM the catalysed rate must rise {F_R5:.0f}x.")
    log(f"  k1b + k2b = {KB_SUM:.1e} /M/s is a conversion rate (EXPERIMENT.md's 6e-6 - 1e-4 estimate), not k2.")
    log(f"  R5 at measured rates: N.\n")

    log("STILL MISSING: ligation k_homo and k_hetero (Deng 2024, closed access); K (or the PM/PL ratio) for")
    log("Ala/pyruvate; the rate at mM catalyst and mM PL (only 50-400 mM measured); whether Ala-Ala (the")
    log("product's own dimer) catalyses like Pro-Pro; and k2 itself, which no paper has measured.")
    out.close()


if __name__ == "__main__":
    main()
