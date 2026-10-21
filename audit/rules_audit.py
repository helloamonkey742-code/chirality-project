"""Independent audit of the six numeric thresholds scorecard.py prints (R5 time rule, R6 sensitivity).

Written blind: ocean.py, inherit.py, scorecard.py were not opened until the comparison stage (see COMPARE below).
Run from the repo root:  /opt/anaconda3/bin/python3 audit/rules_audit.py | tee audit/rules_audit.log

PREDICTIONS FROZEN BEFORE THE MONTE CARLO RAN (Parts A-C of audit/rules_audit_result.txt, first 5493 bytes):
    $ head -c 5493 audit/rules_audit_result.txt | shasum -a 256
    b8c5890fc08ac1616075b02174ccbbcb54590c362241ffc3dad0b9105c807f75  -
    $ date
    Mon Oct 19 23:38:57 PDT 2026      (2026-10-20T06:38:57Z)
This script re-checks that hash, keeps those bytes, and appends everything it computes after them
(Part D onwards), so a rerun replaces only the appended part.

Step 1  blind recompute of R5/R6 from journals/ija_draft.md alone (inputs below carry draft line numbers).
Step 2  Monte Carlo of inheritance: seed e0 at the transition, noise 1/N, slow sweep, find e0 at Delta = 2.
Step 3  COMPARE: input-by-input comparison with the code, filled in only after unblinding.
"""
import hashlib
import math
import sys
import numpy as np
from scipy.stats import norm

RES = "audit/rules_audit_result.txt"
FROZEN_BYTES = 5493
FROZEN_SHA = "b8c5890fc08ac1616075b02174ccbbcb54590c362241ffc3dad0b9105c807f75"
NA, YR = 6.02214076e23, 365.25 * 86400
out = []


def say(s=""):
    print(s)
    out.append(s)


# ---------------------------------------------------------------- Step 1: blind recompute
def k2_min(c, tau_yr, kt=10.0):
    return kt / (c * tau_yr * YR)                       # ija_draft.md:448, k2*c*tau >= 10


def n_mol(c, V, f=1.0):
    return c * V * f * 1000 * NA                        # ija_draft.md:440-442, Eq. 12


def e_need(c, V, kt=10.0, F=1.0, delta=2.0):
    return delta * (math.pi * kt) ** 0.25 / math.sqrt(2 * n_mol(c, V)) / F


SET = {  # name: (c [M], V [m^3], tau [yr], draft receipt)
    "pond": (1e-3, 1.0, 100, "804 (c, tau); 621/623 (1 m^3)"),
    "lake_1e7m3": (1e-3, 1e7, None, "623 (c); 631 (1 km^2 x 10 m); tau unstated (1e4-1e5 yr at 327)"),
    "lake_1km3": (1e-3, 1e9, None, "623"),
    "lake_10km3": (1e-3, 1e10, None, "623"),
    "ocean_Enc": (1e-6, 2.7e16, 1e8, "804 (c, tau); 551 (V)"),
    "ocean_Eur": (1e-6, 3.3e18, 1e8, "804; 551"),
    "ocean_Earth": (1e-6, 1.332e18, 1e8, "804; 551"),
}
PRED_R5 = {"pond": 3.17e-6, "lake_1e4yr": 3.17e-8, "lake_1e5yr": 3.17e-9, "ocean": 3.17e-9}
PRED_R6 = {  # F=1, F=0.7 (frozen block, Part B)
    "pond": (4.31e-12, 6.16e-12), "lake_1e7m3": (1.36e-15, 1.95e-15), "lake_1km3": (1.36e-16, 1.95e-16),
    "lake_10km3": (4.31e-17, 6.16e-17), "ocean_Enc": (8.30e-19, 1.19e-18),
    "ocean_Eur": (7.51e-20, 1.07e-19), "ocean_Earth": (1.18e-19, 1.69e-19)}

R5 = {"pond": k2_min(1e-3, 100), "lake_1e4yr": k2_min(1e-3, 1e4), "lake_1e5yr": k2_min(1e-3, 1e5),
      "ocean": k2_min(1e-6, 1e8)}
R6 = {k: (e_need(c, V), e_need(c, V, F=0.7), e_need(c, V, F=0.3)) for k, (c, V, _, _) in SET.items()}

for k, v in PRED_R5.items():
    assert abs(R5[k] / v - 1) < 0.005, (k, R5[k], v)
for k, (a, b) in PRED_R6.items():
    assert abs(R6[k][0] / a - 1) < 0.005 and abs(R6[k][1] / b - 1) < 0.005, (k, R6[k], a, b)
assert abs(norm.cdf(2) - 0.97725) < 1e-5

say("=== PART D. STEP 1 RECOMPUTE (rules_audit.py; matches the frozen Part B to 0.5%) ===")
for k, v in R5.items():
    say(f"R5 k2_min {k:12s} {v:.3e}")
for k, (a, b, c3) in R6.items():
    say(f"R6 e_need {k:12s} F=1 {a:.3e}  F=0.7 {b:.3e}  F=0.3 {c3:.3e}   N={n_mol(SET[k][0], SET[k][1]):.3e}  [draft {SET[k][3]}]")


# ---------------------------------------------------------------- Step 2: Monte Carlo of inheritance
def mc(N, kt, n=60000, rs=(0.6, 0.8, 1.0, 1.25), seed=0):
    """d(alpha) = (gamma t alpha - alpha^3) dt + sqrt(eps) dW from alpha(0) = r*e_formula; time in 1/k."""
    gam, eps = 1.0 / kt, 1.0 / N
    ef = 2 * (math.pi * kt) ** 0.25 / math.sqrt(2 * N)
    sg = math.sqrt(gam)
    dt = 0.004 / sg
    steps = int(round(5 / sg / dt))
    rng = np.random.default_rng(seed)
    a = np.outer(np.ones(n), np.array(rs) * ef)
    sq, t = math.sqrt(eps * dt), 0.0
    for _ in range(steps):
        a += (gam * t - a * a) * a * dt + (rng.standard_normal(n) * sq)[:, None]  # shared noise across r
        t += dt
    return ef, np.array(rs), (a > 0).mean(0), n


def fit_ratio(rs, P, n):
    """z_i = Phi^-1(P_i) = A*r_i; e0 at Delta = 2 is e_formula * 2/A. Returns (ratio, se, z)."""
    P = np.clip(P, 0.5 / n, 1 - 0.5 / n)
    z = norm.ppf(P)
    var = P * (1 - P) / (n * norm.pdf(z) ** 2)
    w = 1 / var
    A = (w * rs * z).sum() / (w * rs * rs).sum()
    seA = 1 / math.sqrt((w * rs * rs).sum())
    return 2 / A, 2 * seA / A ** 2, z


# self-check of the fit: exact P for a formula that is 10% too optimistic must give ratio 1/1.1
_rs = np.array([0.6, 0.8, 1.0, 1.25])
assert abs(fit_ratio(_rs, norm.cdf(2.2 * _rs), 10**9)[0] - 1 / 1.1) < 1e-6

CFG_LIN = [(1e6, 1e2), (1e8, 1e2), (1e6, 1e4), (1e8, 1e4), (1e10, 1e4), (1e8, 1e6), (1e10, 1e6), (1e6, 10)]
CFG_STRESS = [(1e5, 1e4), (1e4, 1e4)]
say("")
say("=== PART E. STEP 2 MONTE CARLO (n = 60000 paths, seed 0, shared noise over r) ===")
say("N, k*tau, sigma=sqrt(k*tau/N), e_formula, P(r=1), ratio e0*/e_formula +/- se, z at r = 0.6 0.8 1.0 1.25")
rows = {}
for N, kt in CFG_LIN + CFG_STRESS:
    ef, rs, P, n = mc(N, kt)
    ratio, se, z = fit_ratio(rs, P, n)
    rows[(N, kt)] = (ef, P, ratio, se)
    say(f"{N:.0e} {kt:.0e} sigma={math.sqrt(kt / N):.3f} e_f={ef:.3e} P(1)={P[2]:.4f} "
        f"ratio={ratio:.4f}+/-{se:.4f} z=" + " ".join(f"{x:.3f}" for x in z))
    sys.stdout.flush()

lin = [(N, kt) for N, kt in CFG_LIN]
X = np.array([[1, math.log(N), math.log(kt)] for N, kt in lin])
y = np.array([math.log(rows[c][2]) for c in lin])
w = np.array([1 / (rows[c][3] / rows[c][2]) ** 2 for c in lin])
cov = np.linalg.inv(X.T @ (w[:, None] * X))
beta = cov @ X.T @ (w * y)
expN, expK = -0.5 + beta[1], 0.25 + beta[2]
seN, seK = math.sqrt(cov[1, 1]), math.sqrt(cov[2, 2])
chi2 = float((w * (y - X @ beta) ** 2).sum())
say(f"exponent of N : {expN:.4f} +/- {seN:.4f}   (formula -0.5)")
say(f"exponent of k*tau: {expK:.4f} +/- {seK:.4f}   (formula +0.25)")
say(f"fit chi2 = {chi2:.1f} for 5 dof (weighted, MC errors only)")

checks = {
    "C1 P(r=1) = Phi(2) +/-0.004 in all linear configs": all(abs(rows[c][1][2] - norm.cdf(2)) <= 0.004 for c in lin),
    "C2 ratio in 0.95-1.05 in all linear configs": all(0.95 <= rows[c][2] <= 1.05 for c in lin),
    "C3 exponents -0.50+/-0.03 (N), +0.25+/-0.03 (k*tau)": abs(expN + 0.5) <= 0.03 and abs(expK - 0.25) <= 0.03,
    "C4 ratio outside 0.95-1.05 at sigma = 1": not (0.95 <= rows[(1e4, 1e4)][2] <= 1.05),
}
say("")
say("=== PART F. PREDICTIONS C1-C4 ===")
for k, v in checks.items():
    say(f"{'PASS' if v else 'FAIL'}  {k}")
say(f"max |ratio-1| over linear configs = {max(abs(rows[c][2] - 1) for c in lin):.4f}; "
    f"ratio at sigma=0.32: {rows[(1e5, 1e4)][2]:.4f}; at sigma=1: {rows[(1e4, 1e4)][2]:.4f}")

# ---------------------------------------------------------------- Step 3: COMPARE (filled after unblinding)
YR_C, NA_C = 3.156e7, 6.022e23                            # ocean.py:25-26 (draft side: 3.15576e7, 6.02214076e23)
SC = {"pond": (1e-2, 1e2 * YR_C, math.pi),                # scorecard.py:24  (c M, tau s, V m^3)
      "lake": (1e-4, 1e4 * YR_C, 1e7),                    # scorecard.py:25
      "ocean": (1e-6, 1e8 * YR_C, 1.332e18)}              # scorecard.py:26
K2_SC = 1e-3                                              # scorecard.py:33 passes 1e-3 /M/s to e_need
PRINTED_R5 = {"pond": 3e-7, "lake": 3e-7, "ocean": 3e-9}  # scorecard.py output, '/opt/anaconda3/bin/python3 scorecard.py' lines 2-4
PRINTED_R6 = {"pond": 6e-12, "lake": 3e-14, "ocean": 3e-18}
BLIND_REF = {"pond": ("pond", 3.169e-6, "pond"),          # (R6 blind key, R5 blind value, label)
             "lake": ("lake_1e7m3", 3.169e-8, "lake 1 mM, 1e4 yr"),
             "ocean": ("ocean_Earth", 3.169e-9, "ocean")}
WEAK_FORCE = 5e-18                                        # inherit.py:41


def sc_r5(c, tau):                                        # scorecard.py:33
    return 10 / (c * tau)


def sc_r6(c, V, tau):                                     # inherit.py:30-32
    return 2 * (math.pi * K2_SC * c * tau) ** 0.25 / math.sqrt(2 * c * V * 1000 * NA_C)


say("")
say("=== PART G. STEP 3 COMPARE (code opened only after Parts A-F were frozen/run) ===")
say("G1. do the printed numbers reproduce from the code's own inputs? (assert, %.0e as printed)")
for k, (c, tau, V) in SC.items():
    r5, r6 = sc_r5(c, tau), sc_r6(c, V, tau)
    assert float(f"{r5:.0e}") == PRINTED_R5[k] and float(f"{r6:.0e}") == PRINTED_R6[k], (k, r5, r6)
    say(f"  {k:6s} R5 {r5:.4e} -> {r5:.0e}   R6 {r6:.4e} -> {r6:.0e}   k*tau at k2=1e-3: {K2_SC * c * tau:.3e}")

say("G2. R5 input by input, code vs draft (blind)")
for k, (c, tau, V) in SC.items():
    key, ref, lab = BLIND_REF[k]
    cd, td = {"pond": (1e-3, 100), "lake": (1e-3, 1e4), "ocean": (1e-6, 1e8)}[k]
    say(f"  {k:6s} code c={c:.0e} M tau={tau / YR_C:.0e} yr -> {sc_r5(c, tau):.3e} | draft c={cd:.0e} M tau={td:.0e} yr -> "
        f"{ref:.3e} | code/draft = {sc_r5(c, tau) / ref:.3f}")
assert abs(sc_r5(*SC["pond"][:2]) / 3.169e-6 - 0.1) < 0.002 and abs(sc_r5(*SC["ocean"][:2]) / 3.169e-9 - 1) < 0.002

say("G3. R6 code vs blind F=1, split into factors (e_need ~ (k*tau)^0.25 / sqrt(c*V))")
for k, (c, tau, V) in SC.items():
    key = BLIND_REF[k][0]
    cb, Vb, _, _ = SET[key]
    blind = R6[key][0]
    f_kt = (K2_SC * c * tau / 10) ** 0.25
    f_N = math.sqrt(n_mol(cb, Vb) / (c * V * 1000 * NA_C))
    code = sc_r6(c, V, tau)
    assert abs(code / blind / (f_kt * f_N) - 1) < 0.005, (k, code / blind, f_kt * f_N)   # factors reproduce the gap
    say(f"  {k:6s} code {code:.3e} blind(F=1) {blind:.3e} ratio {code / blind:.2f} = k*tau factor {f_kt:.2f} x N factor {f_N:.3f}"
        f"   | code/F: F=0.7 {code / 0.7:.2e}  F=0.3 {code / 0.3:.2e}")
wf = sc_r6(*SC["ocean"][:1], SC["ocean"][2], SC["ocean"][1])
say(f"  ocean e_need {wf:.2e} vs weak-force excess {WEAK_FORCE:.0e} [inherit.py:41]: F=1 ratio {WEAK_FORCE / wf:.2f}, "
    f"F=0.7 {WEAK_FORCE / (wf / 0.7):.2f}, F=0.3 {WEAK_FORCE / (wf / 0.3):.2f}")
say("  pond R6 equals the blind F=0.7 value (6.16e-12) only by coincidence: code 5.77e-12 comes from c*V=0.0314 and k*tau=3.2e4, not from F.")

say("G4. decisions (reasoned from the files, no new numbers):")
say("  scorecard.py:38-40 - R5 letters Y rest on hours-days lab rates vs a floor <= 3e-7; N letters rest on k2 <= 0 [README.md:516-517]. Floor moves x0.1-x10: no R5 letter changes.")
say("  R6 letters rest on 'works from racemic or <= ~1% ee' [scorecard.py:40] vs e_need <= 6e-12: no R6 letter changes.")
say("  aped.py:464 uses floor = 10/tau (no c), so the pond c mismatch does not enter its P10 (pond 3.17e-9 fails, lake and ocean pass, aped.py:74-75).")
say("  draft ija_draft.md:806: lab arm 1.8e-6 'below the pond value' holds for the draft's 3e-6 (1 mM) but not for scorecard's 3e-7 (1e-2 M).")

# ---------------------------------------------------------------- write result file (frozen block kept)
with open(RES, "rb") as f:
    data = f.read()
frozen = data[:FROZEN_BYTES]
assert hashlib.sha256(frozen).hexdigest() == FROZEN_SHA, "frozen predictions were edited"
with open(RES, "wb") as f:
    f.write(frozen + ("\n".join(out) + "\n").encode())
print("\nwrote", RES, "- frozen block sha256 verified")
