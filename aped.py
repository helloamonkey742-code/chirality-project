"""APED: noncatalytic Activation-Polymerization-Epimerization-Depolymerization network
(Plasson, Bersini & Commeyras 2004, PNAS 101:16733, doi:10.1073/pnas.0405293101).
Scorecard row "APED peptides via NCAs (2004-2010)": decides R9 and R5.

MODEL. Closed (residue-conserving), recycled network of chiral amino-acid derivatives, no autocatalyst.
Species L, D, L*, D* (activated, an NCA), LL, LD, DL, DD (dimer-limited; first letter = reacting residue).
15 reactions (each with its L<->D mirror):
  activation X + L <-> L*          deactivation L* <-> L + Y
  polymerization L* + L <-> LL + Y   (rate kP);   L* + D <-> LD + Y   (alpha kP)
  depolymerization LL <-> L + L (kD);             LD <-> L + D   (beta kD)
  epimerization of the first residue DL <-> LL (kE; back gamma kE);  racemization L <-> D (kR)
  pool: ctot = l + d + l* + d* + 2(ll + ld + dl + dd), conserved.
Stereoselectivities alpha (heterochiral/homochiral polymerization), beta (same, depolymerization),
gamma (same, epimerization; also Kf,LD/Kf,LL, the equilibrium ratio). alpha, beta, gamma < 1 favours
homochiral dimers; > 1 favours heterochiral ones (the Danger et al. 2010 sign: DKP formation removes
heterochiral dimers faster and the equilibrium favours them). Back-rates follow microreversibility.

SOURCES. The PNAS paper's full text could not be retrieved (PMC and publisher pages are behind a
reCAPTCHA; not bypassed). The model, rate laws and parameters are taken instead from the authors'
open arXiv versions: Plasson & Bersini 2009, J Phys Chem B 113:3477 (arXiv:0804.4834v3; microreversible
form: Eqs. 30-34 and Table 2, "Analysis of the System", Figs 5-8) and Plasson comment arXiv:0804.3939
(the PNAS model is its irreversible limit). Dimensionless "real system" values alpha = 0.35, beta = 0.2,
gamma = 0.3, eps = eta = 10, a = 1e-8 /s are quoted from Gleiser & Walker (arXiv:0810.5398, Sec. II-III),
not read from Plasson 2004. The 2009 set differs from the PNAS set (kR added, microreversible).
Base parameters (2009): ctot = 2 M, kA = kH = kD = kE = 1 /s, kP = 1 /M/s, kR = 1e-4 /s, alpha = beta =
gamma = 0.1, Kact = 1e-6, Kdim = 1e-3 /M, Kf,x = 1e-9 M, Vy = 1. Reported there: racemic until
x_c = 0.0105 M, a non-racemic branch with ee peaking near x = 0.15 M, racemic again above x = 2.11 M.

DRIVE. x = concentration of the activating fuel X (held constant). Vx = x/Kf,x. The activation step of
Table 2 can be read two ways; both are in the code (p["act"]). Reading A (the first run, registered with
the predictions): forward rate kA [/s] independent of x, reverse kA/(Kact Vx). Reading B (final): forward
kA x with kA = 1 /M/s, reverse kA Kf,x/Kact, the ordinary mass-action law at constant x. B was adopted
after run 1 (see RESULTS: under A Re(lambda) saturates and there is no upper threshold); no numerical
parameter was changed. Both have the same activation equilibrium L*/L = Kact Vx. Thermodynamic
equilibrium (zero drive) is Vx = Vy, i.e. x = Kf,x Vy = 1e-9 M; the affinity is ln(Vx/Vy) = ln(x/1e-9 M).

METHOD. Racemic state: fixed point of the 4-variable symmetric subspace (l = d, l* = d*, ll = dd,
ld = dl), from a monomer-only start, polished by Newton with the pool as a constraint. Difference mode
(D - L, antisymmetric subspace): 4x4 block of the full 8x8 Jacobian, central finite differences
(step 1e-6 max(z_i, 1e-3 max z); rhs is quadratic, so central differences are exact up to rounding) and
an analytic mass-action Jacobian, which must agree. Continuous time, so the racemic
state is stable iff every Re(lambda) < 0 (not |lambda| < 1 as in the discrete-map files); complex
pairs are allowed and tested for.

R9 (scorecard): racemic state stable at low drive, lost as the drive rises. Here: Re(lambda_lead)
< 0 at x -> equilibrium and > 0 above a threshold x_c. Y if both literature parameter sets (homochiral-
favouring Plasson and heterochiral-favouring Danger-like) show a threshold; N if neither; P if only the
Plasson sign does (the sign of alpha, gamma against 1 is the contested input).
R5 (scorecard): k2 >= 10/(c tau). The model is invariant under a common rescaling of time, so only rate
ratios are fixed. sigma = max_x Re(lambda_lead) [1/s] at the base rates; a = apparent activation rate at
the optimum (kA x* under B, kA under A); s* = sigma/a; effective k2 = sigma/ctot, e = k2_eff/kP =
sigma/(kP ctot). Rescaling time so that a = 1e-8 /s (the one real-system rate quoted, a tunable fuel-supply
rate and not a measured one) gives sigma = s* a, compared with the floors sigma >= 10/tau: pond 3.2e-9,
lake 3.2e-11, ocean 3.2e-15 /s (scorecard SETTINGS). Y if all three pass, N if all fail; mixed or
a-dependent (sigma is proportional to a) -> "?" and name the missing number (the measured apparent NCA
activation rate a at each setting's concentration).

PREDICTIONS (dated 2026-10-19, written before any run)
Base set (2009 parameters above, alpha = beta = gamma = 0.1), drive x in M.
(P1) Lower threshold x_c1, where Re(lambda_lead) crosses 0 upward, lies in [0.005, 0.021] M
     (within 2x of the paper's 0.0105 M).
(P2) Upper threshold x_c2, where it crosses back below 0, lies in [1.05, 4.2] M (within 2x of 2.11 M).
(P3) Re(lambda_lead) < 0 for x < x_c1 and x > x_c2, > 0 between; Im(lambda_lead) = 0 at both crossings
     (steady pitchfork) and no complex leading pair for any x in 1e-9..1e4 M.
(P4) Finite-difference and analytic difference-mode Jacobians agree to max|diff|/max|J| < 1e-5.
(P5) Seeded at ee = 1e-6: x = 0.15 M grows to |ee| > 0.5; x = 0.004 M and x = 8 M decay to |ee| < 1e-6.
(P6) No hysteresis (supercritical): integrating down from a homochiral state at x = 0.15 M, ee stays
     > 1e-3 until x is within 10% of x_c1; integrating up, until within 10% of x_c2.
(P7) alpha >= 1 never breaks symmetry: for alpha in {1, 3, 10}, beta in {0.1, 0.2, 1, 10},
     gamma in {0.1, 0.3, 1, 3, 10} (60 sets) Re(lambda_lead) < 0 at every x in 1e-6..1e3 M.
(P8) The Plasson "real system" set (alpha, beta, gamma) = (0.35, 0.2, 0.3) at the same rates is unstable
     in some window of x.
(P9) alpha = 0.1 with beta = gamma = 1 (selective polymerization only) is stable at every x.
(P10) s* = max_x Re(lambda_lead)/kA lies in [0.01, 0.5]; at a = kA = 1e-8 /s the pond floor fails
     (sigma < 3.2e-9 /s) while the lake and ocean floors pass.
(P11) Proposed verdicts: R9 = P (threshold for Plasson's signs, none for Danger-like ones), R5 = ?.

RESULTS (aped.log, 2026-10-19; Ln = line n of aped.log)
MODEL CHANGE AFTER RUN 1. Run 1 used reading A. It crashed in the threshold step (IndexError: no downward
crossing) because Re(lambda_lead) never returns below 0: it saturates at 1.5239e-3 /s for x >= 1e2 M and
x_c1 = 1.1206e-3 M (L5). Scored under A as registered: P1 failed (1.12e-3 is outside [0.005, 0.021]), P2
failed (no upper threshold), P3 failed; P4-P11 were not reached, so 0 of 3 evaluated held. I then adopted
reading B and reran everything; no numerical parameter changed (L2). B is the standard mass-action law,
but it was chosen after seeing that A cannot give an upper threshold, so P1 and P2 holding under B is not
an independent test of the reading. Score below is under B: 8 of 11 held (P10 is a conjunction of two
halves and counts as wrong; its second half held).
(P1) held: x_c1 = 0.010457 M (L36; paper 0.0105).
(P2) held: x_c2 = 2.1133 M (L37; paper 2.11).
(P3) held: Re<0 / >0 / <0 as predicted (L38); the leading eigenvalue is real everywhere on the scan, max
     |Im| = 0 (L39). A complex pair exists but is not leading (x = 0.15 M, L41).
(P4) held: FD vs analytic Jacobian max |dJ|/max|J| = 1.84e-10 (L44).
(P5) held: x = 0.004 decays and x = 8 decays (ee ~ 6e-12), x = 0.15 grows to ee = +0.764; measured early
     rates equal the eigenvalues to 5 digits (L47-L49); pool drift <= 5.2e-14.
(P6) WRONG: the lower threshold shows hysteresis. The down-sweep keeps |ee| > 1e-3 to x = 0.0089349 M =
     0.854 x_c1 (L52), outside the 10% band, with |ee| = 0.0279-0.475 over the last nine points (L55).
     The upper threshold shows none (1.0026 x_c2, L53). The log's "within 10% of thresholds: True" (L54) is
     a one-sided test (x_down <= 1.1 x_c1) that any below-threshold ratio passes; the registered
     two-sided claim failed. Caveat: |lambda| ~ 1e-5 /s near x = 0.01 M (L23) against 5e4 s per step, so
     critical slowing down over ~3 steps could account for part of the persistence; a subcritical branch
     versus sweep lag was not separated (no branch continuation done).
(P7) held: alpha >= 1 sets with an unstable window: 0 of 60 (L83).
(P8) WRONG: (0.35, 0.2, 0.3) is stable at these rates, max Re(lambda) = -1.8971e-4 /s (L84; grid L65).
(P9) held: (0.1, 1, 1) max Re(lambda) = -1.3267e-3 /s, 0 unstable points (L85).
(P10) WRONG as a conjunction. s* = 0.0063166 (a = kA x* = 0.3852 /s; L89) and sigma/kA = 0.0024332 (L89)
     are both below [0.01, 0.5]. Second half held: at a = 1e-8 /s, sigma = 6.32e-11 /s; the pond floor
     3.17e-9 fails, the lake floor 3.17e-11 and ocean floor 3.17e-15 pass (L99-L102).
(P11) held: R9 = P, R5 = ? (L103).
PROPOSED R9 = P. Base set: Re(lambda_lead) = -2.35e-4 /s at equilibrium (L9), above zero between x_c1 =
0.010457 M and x_c2 = 2.1133 M (L36-L37). Plasson signs (alpha, gamma < 1): a window in 5 of 16 sets;
Danger-like signs (alpha, gamma > 1): 0 of 16 (L86); all 60 alpha >= 1 sets are stable (L83); the Plasson
"real-system" set is stable at the 2009 rates (L84). By the registered rule (Y both signs, N neither, P only
the Plasson sign): P. Windows need alpha < 1 and gamma <= 0.3 with beta <= 0.2 (L58-L67).
PROPOSED R5 = ?. sigma is proportional to the unmeasured activation rate a. At a = 1e-8 /s sigma =
6.32e-11 /s: pond fails, lake and ocean pass (L99-L101); a needed 5.0e-7, 5.0e-9, 5.0e-13 /s. Missing
number: the measured apparent NCA activation rate a at each setting's concentration. Second caveat: the
unstable window exists for ctot = 1-1000 M with the other rates fixed and not at ctot = 0.5 (L91-L97);
lower ctot was not scanned, and the scorecard concentrations are 1e-2, 1e-4 and 1e-6 M, so the pass
verdicts use sigma at ctot = 2 M against each setting's own floor.
Runtime: 123.3 s on the first model-B run, 86.0 s on the final rerun (L105, nice -n 10); log otherwise identical.
POST-HOC (2026-10-20, after the run, not registered; changes no verdict). Desfoux, Rossi, Boiteau, Rydzek & Pascal
2025 (Eur. J. Org. Chem. 29(1), doi:10.1002/ejoc.202500786) give one end point, not a rate law: Phe 25 mM +
K3Fe(CN)6 50 mM, pH 6.5, COS for 24 h -> Phe2/Phe3/Phe4 1.4/0.5/0.3%, peptide ureas 47%, peptide hydantoins
16.7% of Phe residues. As first-order a = ln(1/(1-f))/24 h: all products, f = 0.659 -> a = 1.25e-5 /s (all
three pass); free-peptide bonds only, f = 0.0126 -> a = 1.47e-7 /s (pond fails at 0.29x its need, lake and
ocean pass). R5 stays ?: the pond verdict flips on which a counts (ureas and hydantoins are dead ends this
model lacks, and peptide ureas hold some peptide bonds, so the productive a lies between the two); the COS
dose is a lab exposure with no stated concentration; the oxidant is 2 equiv ferricyanide. Settling it needs
a measured a at a stated COS level (Leman et al. 2004, Danger et al. 2006; copies requested 2026-10-20).
POST-HOC (2026-10-20, copies received, not registered; changes no verdict). Same a = ln(1/(1-f))/t. Leman, Orgel
& Ghadiri 2004 (Science 306, 283, doi:10.1126/science.1102722) Table 1 entry 1: Phe 50 mM, COS ~400 mM, borate
pH 9.6, 25 C, 56 h -> dipeptide 6.8%. The % basis is not stated, so f = 0.034-0.068 -> a = 1.7e-7 to 3.5e-7 /s:
pond fails at 0.34-0.70x, lake and ocean pass. Thiocarbamate-start runs (entries 2-4) give 6-7% in 41-65 h, so
cyclization limits the rate and more COS would not lift a. With 100 mM K3Fe(CN)6 (entry 9, 36.3% in 20 min)
a = 1.7e-4 /s, all three pass. Danger, Boiteau, Cottet & Pascal 2006 (JACS 128, 7412, doi:10.1021/ja061339+)
Table 1 entry 3: carbamoyl-Ala 20 mM + Ala 100 mM, 80 C, 53 d -> Ala2 12.7%, Ala3 1.9% of the activator; one NCA
per new bond over 120 mM Ala -> a = 6.1e-9 /s (pond 0.012x, lake 1.2x, ocean pass), but at 80 C, not 25 C, and
the SI time courses were not received. R5 stays ?: every measured uncatalysed a fails the pond and passes the
lake and ocean, the mixed case; only an added oxidant or metal salt passes all three.
CAVEAT ON R5 (added 2026-10-21 from aped_r5_scan2.py; changes no verdict): a measured a is necessary, not sufficient.
sigma = s* a (s* = 0.0063166) is fitted at a/kH = 0.385 (L89) and does not hold at the lab rates for plain COS. In a slice at
x = 1e4 M the difference mode is stable (sigma < 0) below a/kH = 1.6e-3 (kR = 0, kP = 1), 1.0e-2 (kR = 1e-4, kP = 1; this
base set) and 4.1e-4 (kR = 0, kP = 22) (aped_r5_scan2.log:19, 37, 55); the floor a >= 10/(tau s*) is not the whole test.
DO NOT IMPORT: the scan and prints below run at import (~86 s) and there is no __main__ guard. lit_r6.py execs this file
up to its first header print; aped_r5.py, aped_r5_scan.py and aped_r5_scan2.py copy the functions. Keep both working.
"""
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

T0 = time.time()
YR = 3.156e7
SETTINGS = {"Warm little pond": (1e-2, 1e2 * YR),  # c [M], tau [s]; copied from scorecard.py
            "Lake / lagoon": (1e-4, 1e4 * YR),
            "Open ocean": (1e-6, 1e8 * YR)}
N = ["l", "d", "ls", "ds", "ll", "ld", "dl", "dd"]
IX = {n: i for i, n in enumerate(N)}
W = np.array([1, 1, 1, 1, 2, 2, 2, 2.0])  # residues per species
W4 = np.array([2, 2, 4, 4.0])  # same for the symmetric coordinates (l, ls, ll, ld), each doubled
SL, SD = [0, 2, 4, 5], [1, 3, 7, 6]  # L-side species and their D mirrors
BASE = dict(ctot=2.0, kA=1.0, kH=1.0, kP=1.0, kD=1.0, kE=1.0, kR=1e-4, al=0.1, be=0.1, ga=0.1,
            Kact=1e-6, Kdim=1e-3, Kfx=1e-9, Vy=1.0, act="B")
BASEA = dict(BASE, act="A")  # first-run reading, kept only to document why it was dropped


def net(x, p):
    """Mass-action network at fuel x: (stoichiometry S, reactant counts, product counts, kf, kr)."""
    Vx, Vy, Ka, Kd, al, be, ga = x / p["Kfx"], p["Vy"], p["Kact"], p["Kdim"], p["al"], p["be"], p["ga"]
    kA, kH, kP, kD, kE, kR = (p[k] for k in ("kA", "kH", "kP", "kD", "kE", "kR"))
    rx = []
    for h, o in ("ld", "dl"):  # reacting hand h, other hand o; dimer name = reacting residue first
        hh, ho, oh = h + h, h + o, o + h
        fA, rA = (kA, kA / (Ka * Vx)) if p["act"] == "A" else (kA * x, kA * p["Kfx"] / Ka)
        rx += [([h], [h + "s"], fA, rA),  # activation X + h <-> h* ("B": kA [/M/s] x; "A": kA [/s], reverse kA/(Kact Vx))
               ([h + "s"], [h], kH, kH * Ka * Vy),  # deactivation h* <-> h + Y
               ([h + "s", h], [hh], kP, kP * Ka * Vy / Kd),  # polymerization (+Y)
               ([h + "s", o], [ho], al * kP, al * kP * Ka * Vy / (ga * Kd)),
               ([hh], [h, h], kD, kD * Kd),  # depolymerization
               ([ho], [h, o], be * kD, be * ga * kD * Kd),
               ([oh], [hh], kE, ga * kE)]  # epimerization of the first residue
    rx.append((["l"], ["d"], kR, kR))  # racemization
    Rm, Pm = np.zeros((len(rx), 8)), np.zeros((len(rx), 8))
    for j, (a, b, _, _) in enumerate(rx):
        for s in a:
            Rm[j, IX[s]] += 1
        for s in b:
            Pm[j, IX[s]] += 1
    return (Pm - Rm).T, Rm, Pm, np.array([r[2] for r in rx]), np.array([r[3] for r in rx])


def prod(z, m):
    return np.prod(z ** m, axis=1)


def dprod(z, m):  # d(prod z^m)/dz_i, shape (reactions, 8)
    D = np.empty_like(m)
    for i in range(8):
        t = z ** m
        t[:, i] = m[:, i] * z[i] ** np.maximum(m[:, i] - 1, 0)
        D[:, i] = t.prod(1)
    return D


def rhs(z, n):
    S, Rm, Pm, kf, kr = n
    return S @ (kf * prod(z, Rm) - kr * prod(z, Pm))


def jac(z, n):  # analytic
    S, Rm, Pm, kf, kr = n
    return S @ (kf[:, None] * dprod(z, Rm) - kr[:, None] * dprod(z, Pm))


def jac_fd(z, n):  # central differences
    J, zm = np.empty((8, 8)), z.max()
    for i in range(8):
        h = 1e-6 * max(z[i], 1e-3 * zm)
        d = np.zeros(8)
        d[i] = h
        J[:, i] = (rhs(z + d, n) - rhs(z - d, n)) / (2 * h)
    return J


def anti(J):  # difference (D - L) block: perturbation +1 on L-side, -1 on D-side of each pair
    return J[np.ix_(SL, SL)] - J[np.ix_(SL, SD)]


def sym(J):
    return J[np.ix_(SL, SL)] + J[np.ix_(SL, SD)]


def full(y):
    z = np.empty(8)
    z[SL], z[SD] = y, y
    return z


def ee(z, p):
    return (z[0] + z[2] + 2 * z[4] - z[1] - z[3] - 2 * z[7]) / p["ctot"]


def newton(y, n, p):
    """Symmetric fixed point by Newton on the 4 reduced equations, one replaced by the pool constraint."""
    for _ in range(60):
        z = full(y)
        F, J = rhs(z, n)[SL], sym(jac(z, n))
        F[0], J[0] = W4 @ y - p["ctot"], W4
        dy = np.linalg.solve(J, -F)
        t = 1.0
        while np.any(y + t * dy <= 0) and t > 1e-8:
            t /= 2
        y = y + t * dy
        if np.max(np.abs(t * dy) / y) < 1e-10 and t == 1.0:
            return y
    return None


def racemic(x, p, y0=None):
    """Racemic (symmetric-subspace) steady state at fuel x: Newton from y0 if given, else relax the
    4-variable symmetric ODE from a monomer-only pool, then Newton-polish."""
    n = net(x, p)
    if y0 is not None:
        y = newton(np.array(y0, float), n, p)
        if y is not None:
            return y
    y = np.array([p["ctot"] / 2, 1e-30, 1e-30, 1e-30])
    for T in (1e4, 1e7, 1e10):
        sol = solve_ivp(lambda t, v: rhs(full(v), n)[SL], [0, T], y, method="Radau",
                        jac=lambda t, v: sym(jac(full(v), n)), rtol=1e-9, atol=1e-22)
        y = np.maximum(sol.y[:, -1], 1e-30)
        yn = newton(y, n, p)
        if yn is not None:
            return yn
    raise RuntimeError(f"no symmetric fixed point at x={x}")


def lam(x, p, y):
    """Leading (largest real part) eigenvalue of the difference-mode block, analytic Jacobian."""
    ev = np.linalg.eigvals(anti(jac(full(y), net(x, p))))
    return ev[np.argmax(ev.real)]


def fdcheck(x, p, y):
    n, z = net(x, p), full(y)
    Ja, Jf = anti(jac(z, n)), anti(jac_fd(z, n))
    e1, e2 = np.linalg.eigvals(Ja), np.linalg.eigvals(Jf)
    return np.abs(Ja - Jf).max() / np.abs(Ja).max(), abs(e1.real.max() - e2.real.max())


def evolve(z, x, p, T, ts=None):
    n = net(x, p)
    sol = solve_ivp(lambda t, v: rhs(v, n), [0, T], z, method="Radau", jac=lambda t, v: jac(v, n),
                    rtol=1e-12, atol=1e-22, t_eval=ts)
    assert sol.success, sol.message
    return sol


def seeded(x, p, ee0=1e-6):
    """Integrate the full 8-species ODE from the racemic state plus an ee0 seed; return the final state,
    ee(t_end), the analytic lambda and the measured early rate of ee (fit over 1e-5<|ee|<1e-3, or
    t in [0.1 T, 0.6 T] when decaying)."""
    y = racemic(x, p)
    s = lam(x, p, y).real
    z = full(y)
    z[0] += 0.5 * ee0 * p["ctot"]
    z[1] -= 0.5 * ee0 * p["ctot"]
    T = min(60 / s if s > 0 else 12 / -s, 1e9)
    sol = evolve(z, x, p, T, np.linspace(0, T, 2001))
    e = np.array([ee(c, p) for c in sol.y.T])
    w = (np.abs(e) > 1e-5) & (np.abs(e) < 1e-3) if s > 0 else (sol.t > 0.1 * T) & (sol.t < 0.6 * T)
    rate = np.polyfit(sol.t[w], np.log(np.abs(e[w])), 1)[0] if w.sum() > 3 else np.nan
    return sol.y[:, -1], e[-1], s, rate


def sweep(z, xs, p, T=5e4):
    out = []
    for x in xs:
        z = evolve(z, x, p, T).y[:, -1]
        out.append(abs(ee(z, p)))
    return np.array(out)


def scan(p, xs, refine=False):
    """Continuation along xs (integration for the first point). Returns lead eigenvalues (complex)."""
    y, L, ys = None, [], []
    for x in xs:
        y = racemic(x, p, y)
        ys.append(y)
        L.append(lam(x, p, y))
    return np.array(L), ys


def aapp(x, p):  # apparent activation rate [1/s] at fuel x
    return p["kA"] * (1.0 if p["act"] == "A" else x)


def smax(p, xs=np.geomspace(1e-4, 1e4, 65)):
    """Largest Re(lambda) [1/s] over the drive, refined by a bounded search around the best grid point."""
    L, ys = scan(p, xs)
    i = int(np.argmax(L.real))
    if L[i].real <= 0:
        return L[i].real, xs[i]
    f = lambda lx: -lam(10 ** lx, p, racemic(10 ** lx, p, ys[i])).real
    r = minimize_scalar(f, bounds=(np.log10(xs[max(i - 1, 0)]), np.log10(xs[min(i + 1, len(xs) - 1)])),
                        method="bounded", options=dict(xatol=1e-5))
    return -r.fun, 10 ** r.x


def crossing(xs, re, ys, p, i):
    f = lambda x: lam(x, p, racemic(x, p, ys[i - 1])).real
    return brentq(f, xs[i - 1], xs[i], xtol=1e-14, rtol=1e-10)


# ------------------------------------------------------------------------------------------------
print("APED racemic-state stability (aped.py); base = Plasson & Bersini 2009 parameters")
print("base:", {k: v for k, v in BASE.items()})

xs = np.geomspace(1e-9, 1e4, 131)
LA, ysA = scan(BASEA, xs)
upA = [i for i in range(1, 131) if LA[i - 1].real < 0 <= LA[i].real]
xcA = crossing(xs, LA.real, ysA, BASEA, upA[0])
print("\n== 0. first-run model A (activation kA [/s] independent of x; only the reverse rate carries x) ==")
print(f"[A] lower threshold x_c1 = {xcA:.5g} M; downward crossings on 1e-9..1e4 M: "
      f"{sum(LA[i - 1].real > 0 >= LA[i].real for i in range(1, 131))}; "
      f"Re(lambda_lead) at 1e2, 1e3, 1e4 M = {LA[xs >= 1e2][0].real:.6e}, {LA[xs >= 1e3][0].real:.6e}, {LA[-1].real:.6e} /s "
      f"(saturates: x enters only through the reverse activation, so no upper threshold is possible)")

# 1. drive scan of the difference-mode leading eigenvalue (integration from a monomer pool at every x)
ys = [racemic(x, BASE) for x in xs]
Lm = np.array([lam(x, BASE, y) for x, y in zip(xs, ys)])
fdm = [fdcheck(x, BASE, y) for x, y in zip(xs, ys)]
print("\n== 1. leading difference-mode eigenvalue vs drive x (every 5th of 131 points) ==")
print("      x [M]    affinity ln(Vx/Vy)     Re(lambda) [1/s]    Im(lambda)")
for i in range(0, 131, 5):
    print(f"  {xs[i]:9.3e}   {np.log(xs[i] / BASE['Kfx']):8.2f}         {Lm[i].real:+12.5e}     {Lm[i].imag:+.2e}")
up = [i for i in range(1, 131) if Lm[i - 1].real < 0 <= Lm[i].real]
dn = [i for i in range(1, 131) if Lm[i - 1].real > 0 >= Lm[i].real]
xc1 = crossing(xs, Lm.real, ys, BASE, up[0])
xc2 = crossing(xs, Lm.real, ys, BASE, dn[0])
yc1, yc2 = racemic(xc1, BASE, ys[up[0] - 1]), racemic(xc2, BASE, ys[dn[0] - 1])
print(f"[P1] lower threshold x_c1 = {xc1:.5g} M  (paper: 0.0105)   upward crossings on grid: {len(up)}")
print(f"[P2] upper threshold x_c2 = {xc2:.5g} M  (paper: 2.11)     downward crossings on grid: {len(dn)}")
print(f"[P3] sign: Re<0 below x_c1 and above x_c2: "
      f"{bool((Lm.real[xs < xc1] < 0).all() and (Lm.real[xs > xc2] < 0).all())}; "
      f"Re>0 between: {bool((Lm.real[(xs > xc1) & (xs < xc2)] > 0).all())}")
print(f"[P3] Im(lambda_lead) at x_c1 = {lam(xc1, BASE, yc1).imag:.2e}, at x_c2 = {lam(xc2, BASE, yc2).imag:.2e}; "
      f"max |Im| of the leading eigenvalue over the 131-point scan = {np.abs(Lm.imag).max():.2e}")
ev = lambda x, y: np.sort_complex(np.linalg.eigvals(anti(jac(full(y), net(x, BASE)))))
for x in (xc1, 0.15, xc2):
    y = racemic(x, BASE)
    print(f"     all 4 difference-mode eigenvalues at x={x:.4g}: {np.round(ev(x, y), 6)}")
print(f"[P4] finite-difference vs analytic difference-mode Jacobian over the 131 points: max |dJ|/max|J| = "
      f"{max(f[0] for f in fdm):.2e}; max |d Re(lambda_lead)| = {max(f[1] for f in fdm):.2e}")

# 2. seeded integrations of the full 8-species ODE
print("\n== 2. seeded integration, ee(0) = 1e-6 ==")
res = {}
for x in (0.004, 0.15, 8.0):
    zf, ef, s, rate = seeded(x, BASE)
    res[x] = zf
    print(f"[P5] x = {x:g} M: lambda_lead = {s:+.4e}/s, measured early rate of ee = {rate:+.4e}/s, "
          f"final ee = {ef:+.4e}, pool drift = {W @ zf - BASE['ctot']:+.1e}")

# 3. hysteresis: sweep from the homochiral state at x = 0.15 M
print("\n== 3. hysteresis: quasi-static sweeps from the homochiral state at x = 0.15 M (40 steps/decade, 5e4 s each) ==")
dx = 10 ** (1 / 40)
xd = 0.15 / dx ** np.arange(1, 57)
xu = 0.15 * dx ** np.arange(1, 57)
ed, eu = sweep(res[0.15], xd, BASE), sweep(res[0.15], xu, BASE)
x_down, x_up = xd[ed > 1e-3].min(), xu[eu > 1e-3].max()
print(f"     down: |ee| at x = 0.15 {abs(ee(res[0.15], BASE)):.4f}; last x with |ee| > 1e-3: {x_down:.5g} M; "
      f"ratio to x_c1 = {x_down / xc1:.4f}")
print(f"     up:   last x with |ee| > 1e-3: {x_up:.5g} M; ratio to x_c2 = {x_up / xc2:.4f}")
print(f"[P6] hysteresis present (down ratio < 0.95 or up ratio > 1.05): "
      f"{bool(x_down / xc1 < 0.95 or x_up / xc2 > 1.05)}; "
      f"within 10% of thresholds: {bool(x_down <= 1.1 * xc1 and x_up >= xc2 / 1.1)}")
print("     |ee| on the down-sweep near x_c1:", " ".join(f"{a:.3g}" for a in ed[(xd > 0.9 * xc1) & (xd < 1.5 * xc1)]))

# 4. stereoselectivity grid (Plasson signs: alpha, beta, gamma < 1; Danger-like: alpha, gamma > 1)
print("\n== 4. stereoselectivity grid: window of x (1e-6..1e3 M, 8 pts/decade) where Re(lambda_lead) > 0; s = stable everywhere ==")
AL, BE, GA = (0.1, 0.35, 1, 3, 10), (0.1, 0.2, 1, 10), (0.1, 0.3, 1, 3, 10)
xg = np.geomspace(1e-6, 1e3, 73)
grid = {}
for a in AL:
    print(f"  alpha = {a}:   gamma ->" + "".join(f"{g:>16g}" for g in GA))
    for b in BE:
        row = f"    beta = {b:<6g}   "
        for g in GA:
            L, _ = scan(dict(BASE, al=a, be=b, ga=g), xg)
            pos = xg[L.real > 0]
            grid[(a, b, g)] = (L.real.max(), pos)
            row += f"{'s':>16}" if len(pos) == 0 else f"{pos[0]:8.1e}-{pos[-1]:.0e}"[:16].rjust(16)
        print(row)
n7 = sum(len(grid[(a, b, g)][1]) > 0 for a in (1, 3, 10) for b in BE for g in GA)
print(f"[P7] alpha >= 1 sets with an unstable window: {n7} of {3 * len(BE) * len(GA)}")
r8 = grid[(0.35, 0.2, 0.3)]
print(f"[P8] (0.35, 0.2, 0.3): max Re(lambda) = {r8[0]:+.4e}; unstable window "
      f"{'none' if len(r8[1]) == 0 else f'{r8[1][0]:.2e} .. {r8[1][-1]:.2e} M'}")
print(f"[P9] (alpha, beta, gamma) = (0.1, 1, 1): max Re(lambda) = {grid[(0.1, 1, 1)][0]:+.4e}; unstable points "
      f"{len(grid[(0.1, 1, 1)][1])}")
plasson = [k for k in grid if k[0] < 1 and k[2] < 1 and len(grid[k][1])]
danger = [k for k in grid if k[0] > 1 and k[2] > 1 and len(grid[k][1])]
print(f"     Plasson-sign sets (alpha<1, gamma<1) with a threshold: {len(plasson)} of "
      f"{sum(k[0] < 1 and k[2] < 1 for k in grid)}; Danger-like sets (alpha>1, gamma>1) with a threshold: "
      f"{len(danger)} of {sum(k[0] > 1 and k[2] > 1 for k in grid)}")
R9 = "Y" if plasson and danger else "N" if not plasson and not danger else "P"

# 5. R5: growth rate in ocean.py units
print("\n== 5. R5: growth rate of the difference mode ==")
sig, xstar = smax(BASE)
a_star = aapp(xstar, BASE)
s_a, e_eff = sig / a_star, sig / (BASE["kP"] * BASE["ctot"])
print(f"[P10] base set: sigma_max = {sig:.5g} /s at x = {xstar:.4g} M; apparent activation rate there a = kA x = {a_star:.4g} /s; "
      f"s* = sigma/a = {s_a:.5g} (sigma/kH = {sig / BASE['kH']:.5g}); e = k2_eff/kP = sigma/(kP ctot) = {e_eff:.5g}")
print("      sigma_max(ctot) with the other base rates fixed (x in 1e-4..1e4 M):")
for g in (0.5, 1, 2, 5, 20, 100, 1000):
    sg, xm = smax(dict(BASE, ctot=g))
    print(f"        ctot = {g:<6g} sigma_max = {sg:+.4e} /s at x = {xm:.3e} M  (s* = {sg / aapp(xm, BASE):+.4e}, e = {sg / (BASE['kP'] * g):+.4e})")
a_q = 1e-8
passes = []
print(f"      rescale time so that a = {a_q:g} /s (the one quoted real-system rate, a tunable fuel-supply rate; sigma = s* a):")
for name, (c, tau) in SETTINGS.items():
    floor = 10 / tau
    ok = s_a * a_q >= floor
    passes.append(ok)
    print(f"        {name:17s} floor sigma >= 10/tau = {floor:.2e} /s (k2 >= {10 / (c * tau):.2e} /M/s); "
          f"sigma = {s_a * a_q:.2e} /s -> {'pass' if ok else 'fail'}; a needed = {floor / s_a:.2e} /s; "
          f"kP needed = {10 / (e_eff * c * tau):.2e} /M/s")
R5 = "Y" if all(passes) else "N" if not any(passes) else "?"
print(f"[P10] pond floor fails: {not passes[0]}; lake and ocean floors pass: {bool(passes[1] and passes[2])}")
print(f"[P11] proposed R9 = {R9} (Plasson sign: {'threshold' if plasson else 'none'}; Danger-like sign: "
      f"{'threshold' if danger else 'none'}); proposed R5 = {R5}"
      f" (missing number if '?': measured apparent NCA activation rate a at each setting's c)")

# ------------------------------------------------------------------------------------------------
# self-check
x0 = BASE["Kfx"] * BASE["Vy"]  # zero drive: Vx = Vy
y0 = racemic(x0, BASE)
n0, z0 = net(x0, BASE), full(y0)
f0, r0 = n0[3] * prod(z0, n0[1]), n0[4] * prod(z0, n0[2])
assert np.max(np.abs(f0 - r0) / np.maximum(f0, r0)) < 1e-6, "detailed balance fails at zero drive"
l0 = z0[0]
for got, want in ((z0[2] / l0, BASE["Kact"]), (z0[4] / l0 ** 2, BASE["Kdim"]),
                  (z0[5] / l0 ** 2, BASE["ga"] * BASE["Kdim"]), (z0[6] / z0[4], BASE["ga"])):
    assert abs(got / want - 1) < 1e-6, "equilibrium constants not recovered"
E0 = np.linalg.eigvals(anti(jac(z0, n0)))
assert (E0.real < 0).all(), "difference mode not stable at zero drive"
zb = z0.copy()
zb[0] += 0.25 * BASE["ctot"]
zb[1] -= 0.25 * BASE["ctot"]
zrel = evolve(zb, x0, BASE, 40 / -E0.real.max()).y[:, -1]
assert abs(ee(zrel, BASE)) < 1e-6 and np.abs(zrel / z0 - 1).max() < 1e-4, "zero drive does not relax to racemic equilibrium"
rng = np.random.default_rng(1)
S0 = net(0.15, BASE)
assert np.abs(W @ S0[0]).max() == 0, "stoichiometry does not conserve residues"
assert all(abs(W @ rhs(rng.random(8), S0)) < 1e-9 for _ in range(20)), "residue conservation"
assert all(abs(W @ res[x] - BASE["ctot"]) < 1e-9 for x in res), "pool drift in integration"
yy = racemic(0.15, BASE)
assert np.abs(rhs(full(yy), S0)).max() < 1e-10, "racemic state is not a fixed point of the full ODE"
assert max(f[0] for f in fdm) < 1e-5, "finite-difference and analytic Jacobians disagree"
for x in (1e-3, 0.0105, 0.15, 2.11, 50):  # continuation agrees with integration
    yc = racemic(x, BASE, racemic(x * 1.2, BASE))
    assert np.abs(yc / racemic(x, BASE) - 1).max() < 1e-6, "continuation branch differs"
assert xc1 / 3 < 0.0105 < xc1 * 3 and xc2 / 3 < 2.11 < xc2 * 3, "Plasson 2009 thresholds not reproduced within x3"
assert lam(1e-3, BASE, racemic(1e-3, BASE)).real < 0 < lam(0.15, BASE, yy).real, "no instability at 0.15 M"
assert lam(20.0, BASE, racemic(20.0, BASE)).real < 0, "racemic state not regained at high drive"
assert abs(ee(res[0.15], BASE)) > 0.1, "seeded run at 0.15 M did not grow"
print(f"\nruntime {time.time() - T0:.1f} s")
print("\nself-check passed")
