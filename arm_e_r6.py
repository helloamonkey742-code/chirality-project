"""Arm E R6: does the fed transamination + ligation network follow an inherited excess >= e_need? (scorecard.py:60, R6 "?")

R6 rule (scorecard.py:9): "follows an inherited excess >= e_need (inherit.py)"; inherit.py: followed = Delta >= 2.
Scorecard comment (scorecard.py:40): "R6 here = shown to work from racemic or <= ~1% ee".

METHOD (all imported from arm_e_rates.py: Higgs & Blackmond 2025 destructive mechanism, Yu 2024 measured rates, Deng 2024
  ligation rates NOT obtained so kh, R, K, f are scanned on arm_e_rates.GRID)
  settings   scorecard.SETTINGS (c, tau, V); C0 = c; Omega = V*1000*NA; N = c*Omega; e_need = inherit.e_need(c, V, 1e-3, tau).
  odd sector (a, b) = (L-D, LL-DD)/2 of the monomer and homodimer; its 2x2 Jacobian J (= arm_e_rates._jac) and the counting-noise
             diffusion D per Omega are derived from the channels (feed, outflow, L+L, D+D, L+D, L<->A catalysed and not).
  case (a)   measured rates: F = 1 (catalysed breakdown KB_SUM split at s = 6.0, uncatalysed KU). Racemic state stable everywhere
             (arm_e_rates.log:24-36), so the network is a linear pass-through of an inherited excess.
  case (b)   the crossing rate multiplier: F_c(setting) = smallest F with lambda_max > 0 on the grid at C0 = c (bisection as
             arm_e_rates.f_crit); the network is run at F_b = 1.5 F_c (sensitivity 1.1, 3, 10, s = 4 and 9.1, kh up to 1e4).
  feed excess g   feed L:D = (1+g):(1-g). Stable: steady response a = G g x, Delta = G g / s (s = readout noise, Lyapunov equation).
                  Unstable: Delta_b = (F_g (1-E)/lam)/sqrt(l'Dl (1-E^2)/(2 lam Omega)), E = exp(-lam tau), l = left eigenvector.
  seed ee0        pool starts with a = ee0 x and b = ee0 y (best case; monomer-only is a sensitivity), read at tau through
                  exp(J tau) (closed form) with the stationary noise (stable) or the accumulated noise of the unstable mode.
  followed (literal)  Delta >= 2. (strict)  literal AND readout ee at tau >= 1% (and, case b, linear growth from the inherited
                  excess reaches 1% within tau; saturated nonlinear ee >= 1%, checked with arm_e_rates.simulate).
  guard (stable points only): f tau >= 5 and -lambda_slow tau >= 5 (pool at steady state by tau); sensitivity 1, 2, 3, 10.
  subsets    U: whole guarded grid. E: margin >= half the setting's max margin (network as engaged as the grid permits).
  letters    Y: e_min <= e_need in all three settings; P: some; N: none.

PRE-REGISTERED PREDICTIONS (arm_e_r6_predictions.txt, written before the first results run)
  sha256 235e93aa8a0e4e93bd8d68cecd97befa58eae407b22873ddb90d0f616ca3b9ab   Wed Oct 21 00:58:54 PDT 2026 (2026-10-21T07:58:54Z)
  (this script asserts the file still has that hash)
  DEBUG RUNS BEFORE THE FINAL RUN (disclosed; the predictions file was not changed by any of them):
    smoke 1 and 2: crashed in self_check (a shape error in lyap(); a tolerance in the seed check). No results were printed.
    smoke 3: a FULL run (129 s) whose whole output I read. It scored P4-P8 FAILED and P6 failed through the ocean f-grid floor. After
      reading it I made only these changes (none touches a pre-registered definition or a P rule; failed predictions stay failed):
      (1) self-check (iii) integrator LSODA(atol=0) -> DOP853 (an 'lsoda: Illegal input' warning); (2) retention nan-handling;
      (3) added reporting: unstable-only (b) rows, monomer-only-seed letters, strict (b) also requires ee_sat >= 1% at the chosen point,
      5 significant digits in the s table; (4) POST HOC section (5b): f grid extended down to 1e-19 /s, NOT pre-registered, not used to
      score any P. Also: the pre-registered label "EXISTING-RULE letter = feed" is mine; inherit.py:3-5 defines the inherited excess as a
      SEED e0 in the pool, so the seed letters are the closer analogue. Both are reported with both labels.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with matplotlib python arm_e_r6.py > arm_e_r6.log
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arm_e_rates as ar  # noqa: E402
from inherit import e_need  # noqa: E402
from ocean import NA  # noqa: E402
from scorecard import SETTINGS  # noqa: E402

PRED_SHA = "235e93aa8a0e4e93bd8d68cecd97befa58eae407b22873ddb90d0f616ca3b9ab"
G_EXT = dict(ar.GRID, kh=np.logspace(-9, 4, 27)[:, None, None, None])  # sensitivity only, never used for letters
G_FLOW = dict(ar.GRID, f=np.logspace(-19, -3, 65))                       # POST HOC sensitivity: f floor 1e-19
NAMES = list(SETTINGS)


# ------------------------------------------------------------------ odd-sector building blocks
def base(C0, F, s, grid=ar.GRID, ku=ar.KU):
    """Racemic state, Jacobian J, diffusion D (per Omega), eigenvalues l1 >= l2, margin, on the whole grid."""
    k1, k2 = ar.kb(s, F)
    x, v = ar.steady(dict(C0=C0, k1=k1, k2=k2, ku=ku, **grid))
    _, f, kh, Rr, K, k1v, k2v, kuv = v
    y, P, A, u = ar._parts(x, *v)
    J11, J12, J21, J22 = -f - P - 4 * kh * x, (k2v - k1v) * u, 2 * kh * x, -f
    D = dict(aa=f * C0 / 4 + f * x / 2 + P * x / 2 + P * K * A / 2 + 2 * kh * x * x,
             bb=f * y / 2 + kh * x * x / 2, ab=-kh * x * x)
    tr, det = J11 + J22, J11 * J22 - J12 * J21
    l2 = (tr - np.sqrt((J11 - J22) ** 2 + 4 * J12 * J21)) / 2
    return dict(C0=C0, F=F, s=s, k1=k1, k2=k2, f=f, kh=kh, R=Rr, K=K, x=x, y=y, A=A, P=P, J=(J11, J12, J21, J22), D=D,
                m=J12 * J21 / (J11 * J22), l1=det / l2, l2=l2, cat=(k1v + k2v) * y / P)


def lyap(J, D):
    """Sigma/Omega from J S + S J' + D = 0 (2x2, closed 3x3 linear solve, time scaled by |J11| for conditioning)."""
    J11, J12, J21, J22 = J
    z = np.zeros_like(J11)
    sc = -J11
    a, b, c, d = J11 / sc, J12 / sc, J21 / sc, J22 / sc
    M = np.stack([np.stack([2 * a, 2 * b, z], -1), np.stack([c, a + d, b], -1), np.stack([z, 2 * c, 2 * d], -1)], -2)
    rhs = (-np.stack([D["aa"], D["ab"], D["bb"]], -1) / sc[..., None])[..., None]
    S = np.linalg.solve(M, rhs)[..., 0]
    return S[..., 0], S[..., 1], S[..., 2]  # Saa, Sab, Sbb


def feed_stable(b, Om):
    """Steady readout excess per unit feed excess (Ga monomer, Gb dimer) and readout noise sd (sa, sb)."""
    J11, J12, J21, J22 = b["J"]
    h = b["f"] * b["C0"] / 2
    det = J11 * J22 - J12 * J21
    Saa, _, Sbb = lyap(b["J"], b["D"])
    return (-J22 * h / det) / b["x"], (J21 * h / det) / b["y"], np.sqrt(Saa / Om) / b["x"], np.sqrt(Sbb / Om) / b["y"]


def seed_stable(b, tau, mono_only=False):
    """Readout excess at tau per unit seed (a, b) = ee0 (x, y) [or (x, 0)], closed form of exp(J tau) w0."""
    J11, J12, J21, J22 = b["J"]
    x, y, l1, l2 = b["x"], b["y"], b["l1"], b["l2"]
    w = 0.0 if mono_only else y
    e1, e2, dl = np.exp(l1 * tau), np.exp(l2 * tau), l1 - l2
    a = (e1 * ((J11 - l2) * x + J12 * w) - e2 * ((J11 - l1) * x + J12 * w)) / dl
    bb = (e1 * (J21 * x + (J22 - l2) * w) - e2 * (J21 * x + (J22 - l1) * w)) / dl
    return a / x, bb / y


def unstable(b, tau, Om):
    """Unstable mode (lam = l1 > 0): per variant (Delta-only e_min, strict e_min with the 1% growth gate)."""
    J11, J12, J21, J22 = b["J"]
    D, x, y, f, C0, l = b["D"], b["x"], b["y"], b["f"], b["C0"], b["l1"]
    la, lb = J21, l - J11                              # left eigenvector (J21, lam - J11); right is (J12, lam - J11)
    lr = J21 * J12 + lb * lb
    nu2 = la * la * D["aa"] + 2 * la * lb * D["ab"] + lb * lb * D["bb"]
    lt = l * tau
    sd = np.sqrt(nu2 * (-np.expm1(-2 * lt)) / (2 * l * Om))
    gain = np.maximum(J12 / x, lb / y) / lr            # readout ee per unit mode amplitude (better of monomer, dimer)
    out = {}
    for k, c in (("feed", la * (f * C0 / 2) * (-np.expm1(-lt)) / l), ("seed", la * x + lb * y), ("seedm", la * x)):
        d = 2 * sd / c
        gate = np.exp(np.minimum(np.log(0.01) - lt - np.log(c * gain), 700))
        out[k] = (d, np.maximum(d, gate))
    return out


def emins(b, tau, Om, guard=5.0):
    """Smallest followed excess at every grid point, per variant/readout. inf where not applicable."""
    l1, f = b["l1"], b["f"]
    st, un = l1 < 0, l1 > 0
    ok = st & (f * tau >= guard) & (-l1 * tau >= guard)
    out = {}
    with np.errstate(all="ignore"):
        Ga, Gb, sa, sb = feed_stable(b, Om)
        src = {"feed": ((Ga, sa), (Gb, sb)),
               "seed": tuple((q, s_) for q, s_ in zip(seed_stable(b, tau), (sa, sb))),
               "seedm": tuple((q, s_) for q, s_ in zip(seed_stable(b, tau, True), (sa, sb)))}
        for k, pair in src.items():
            for nm, (G, s_) in zip("ab", pair):
                lit = 2 * s_ / np.abs(G)
                out[f"{k}_{nm}"] = np.where(ok, lit, np.inf)
                out[f"{k}S_{nm}"] = np.where(ok, np.maximum(lit, 0.01 / np.abs(G)), np.inf)
        u = unstable(b, tau, Om)
        for k, (d, dS) in u.items():
            out[f"u_{k}"], out[f"u_{k}S"] = np.where(un, d, np.inf), np.where(un, dS, np.inf)
        for k in ("feed", "seed", "seedm"):
            out[f"{k}_lit"] = np.fmin(np.fmin(out[f"{k}_a"], out[f"{k}_b"]), out[f"u_{k}"])
            out[f"{k}_strict"] = np.fmin(np.fmin(out[f"{k}S_a"], out[f"{k}S_b"]), out[f"u_{k}S"])
        out["Ga"], out["Gb"], out["ret_m"] = np.where(st, Ga, np.nan), np.where(st, Gb, np.nan), seed_stable(b, tau, True)[0]
    return {k: np.where(np.isnan(v), np.inf, v) if k not in ("Ga", "Gb", "ret_m") else v for k, v in out.items()}


def pick(arr, mask=None):
    a = np.where(np.isnan(arr), np.inf, arr)
    a = np.where(mask, a, np.inf) if mask is not None else a
    i = np.unravel_index(np.argmin(a), a.shape)
    return float(a[i]), i


def case(name, F, s=6.0, grid=ar.GRID, guard=5.0):
    c, tau, V = SETTINGS[name]
    Om = V * 1000 * NA
    b = base(c, F, s, grid)
    A = emins(b, tau, Om, guard)
    E = b["m"] >= 0.5 * b["m"].max()
    P = {(k, u): pick(A[k], E if u == "E" else None) for k in A if k not in ("Ga", "Gb", "ret_m") for u in ("U", "E")}
    return dict(name=name, c=c, tau=tau, V=V, Om=Om, N=c * Om, en=e_need(c, V, 1e-3, tau), b=b, A=A, P=P, E=E)


def rec(cs, key, sub):
    """JSON record of the optimum: parameters, odd-sector numbers and every e_min component at that point."""
    v, i = cs["P"][(key, sub)]
    b = cs["b"]
    r = {k: float(b[k][i]) for k in ("f", "kh", "R", "K", "x", "y", "A", "P", "m", "l1", "l2", "cat")}
    r.update(C0=b["C0"], F=b["F"], s=b["s"], tau=cs["tau"], Om=cs["Om"], N=cs["N"], e_need=cs["en"], e_min=v, key=key, sub=sub)
    r.update({f"J{n}": float(j[i]) for n, j in zip(("11", "12", "21", "22"), b["J"])})
    r.update({f"D{k}": float(d[i]) for k, d in b["D"].items()})
    r["comp"] = {k: (float(a[i]) if np.isfinite(a[i]) else None) for k, a in cs["A"].items() if k != "ret_m"}
    return r


def fcrit(C0, s, grid=ar.GRID):
    ok = lambda F: base(C0, F, s, grid)["l1"].max() > 0
    if ok(1.0):
        return 1.0
    if not ok(1e8):
        return np.inf
    lo, hi = 0.0, 8.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ok(10 ** mid) else (mid, hi)
    return 10 ** hi


def sat(cs, key, sub):
    """Nonlinear saturated monomer ee from a tiny seed at the optimum (arm_e_rates.simulate)."""
    _, i = cs["P"][(key, sub)]
    b = cs["b"]
    p = dict(C0=b["C0"], f=float(b["f"][i]), kh=float(b["kh"][i]), R=float(b["R"][i]), K=float(b["K"][i]),
             k1=b["k1"], k2=b["k2"], ku=ar.KU)
    z0 = ar.racemic_z(p)
    z0[1] = 1e-8 * z0[0]
    try:
        _, Z = ar.simulate(p, z0, 80 / float(b["l1"][i]), n=200)
        ee = Z[1] / Z[0]
        return float(ee[-1]), bool(abs(ee[-1] - ee[-20]) < 1e-3 * abs(ee[-1]))
    except Exception as e:  # noqa: BLE001
        return float("nan"), str(e)[:40]


# ------------------------------------------------------------------ self-check
def self_check():
    # (i) hand case: no dimers (kh = 0), no background (ku = 0): J = -f I, D_aa = f C0/2, Sigma_aa = C0/4, G = 1, s = 1/sqrt(N)
    grid = dict(kh=np.array([0.0]), R=np.array([1.0]), K=np.array([1e-3]), f=np.array([1e-6]))
    with np.errstate(all="ignore"):  # cat share is 0/0 here (no catalysed breakdown)
        b = base(5e-3, 1.0, 6.0, grid, ku=0.0)
    Om = 1e20
    with np.errstate(all="ignore"):  # dimer readout is 0/0 here (no dimers)
        Ga, _, sa, _ = feed_stable(b, Om)
    assert np.isclose(Ga[0], 1.0, rtol=1e-12) and np.isclose(sa[0], 1 / np.sqrt(Om * 5e-3), rtol=1e-12)
    assert np.isclose(2 * sa[0] / Ga[0], 2 / np.sqrt(Om * 5e-3), rtol=1e-12)
    # (ii) grid objects agree with arm_e_rates: Jacobian, lambda sign, Lyapunov vs scipy, expm vs closed form
    rng = np.random.default_rng(3)
    b = base(5e-3, 1.0, 6.0)
    p = dict(C0=5e-3, k1=ar.kb(6.0)[0], k2=ar.kb(6.0)[1], ku=ar.KU, **ar.GRID)
    assert all(np.allclose(u, w, rtol=1e-12) for u, w in zip(b["J"], ar._jac(p)))
    assert np.allclose(b["m"], ar.margin(p), rtol=1e-12) and np.array_equal(b["l1"] > 0, ar.lam(p) > 0)
    Saa, Sab, Sbb = lyap(b["J"], b["D"])
    for _ in range(5):
        i = tuple(int(rng.integers(0, n)) for n in b["m"].shape)
        J = np.array([[b["J"][0][i], b["J"][1][i]], [b["J"][2][i], b["J"][3][i]]])
        Dm = np.array([[b["D"]["aa"][i], b["D"]["ab"][i]], [b["D"]["ab"][i], b["D"]["bb"][i]]])
        S = solve_continuous_lyapunov(J, -Dm)
        assert np.allclose([Saa[i], Sab[i], Sbb[i]], [S[0, 0], S[0, 1], S[1, 1]], rtol=1e-7)
        tq = 1.5 / abs(float(b["l1"][i]))  # a time at which the slow mode is still alive (e^-1.5)
        w = expm(J * tq) @ np.array([b["x"][i], b["y"][i]])
        q = seed_stable(b, tq)
        assert np.allclose([q[0][i] * b["x"][i], q[1][i] * b["y"][i]], w, rtol=1e-6, atol=1e-12 * (b["x"][i] + b["y"][i])), (i, w)
    # (iii) unstable mode: accumulated variance of l.w equals nu2 (e^{2 lam tau} - 1)/(2 lam), by integrating dSigma/dt
    bu = base(5e-3, 4000.0, 6.0)
    i = np.unravel_index(np.argmax(bu["l1"]), bu["l1"].shape)
    assert bu["l1"][i] > 0
    J = np.array([[bu["J"][0][i], bu["J"][1][i]], [bu["J"][2][i], bu["J"][3][i]]])
    Dm = np.array([[bu["D"]["aa"][i], bu["D"]["ab"][i]], [bu["D"]["ab"][i], bu["D"]["bb"][i]]])
    lam = float(bu["l1"][i])
    ell = np.array([J[1, 0], lam - J[0, 0]])
    assert np.allclose(ell @ J, lam * ell, rtol=1e-9)
    from scipy.integrate import solve_ivp
    T = 2.0 / lam
    sol = solve_ivp(lambda t, S: (J @ S.reshape(2, 2) + S.reshape(2, 2) @ J.T + Dm).ravel(), (0, T), np.zeros(4),
                    rtol=1e-10, atol=1e-40, method="DOP853")
    Sn = sol.y[:, -1].reshape(2, 2)
    assert np.isclose(ell @ Sn @ ell, (ell @ Dm @ ell) * np.expm1(2 * lam * T) / (2 * lam), rtol=1e-4)
    # (iv) fcrit equals arm_e_rates.f_crit at 5 mM
    fc = fcrit(5e-3, 6.0)
    assert np.isclose(fc, ar.f_crit(5e-3, 6.0), rtol=1e-4), (fc, ar.f_crit(5e-3, 6.0))
    return fc


# ------------------------------------------------------------------ report
def letter(vals):
    return "Y" if all(v for v in vals) else ("P" if any(vals) else "N")


def main():
    sha = hashlib.sha256((HERE / "arm_e_r6_predictions.txt").read_bytes()).hexdigest()
    assert sha == PRED_SHA, f"predictions file changed: {sha}"
    t0 = time.time()
    def log(*a):
        print(" ".join(str(x) for x in a), flush=True)

    fc5 = self_check()
    log(f"pre-registration file arm_e_r6_predictions.txt: sha256 {sha}")
    log(f"self-check passed: hand case (kh = 0, ku = 0: G = 1, s = 1/sqrt(N)), J/margin/lambda vs arm_e_rates, Lyapunov vs scipy,"
        f" seed closed form vs expm, unstable-mode variance vs ODE, F_c(5 mM, s = 6) = {fc5:.1f} vs arm_e_rates.f_crit\n")
    xn = lambda v, cs: v * np.sqrt(cs["N"])
    RES, REC, chk = {}, {}, []

    def verdict(tag, ok, text):
        chk.append((tag, ok))
        log(f"  {tag} {'held' if ok else 'FAILED'}: {text}")

    # ---- settings
    log("(0) SETTINGS  C0 = c, Omega = V 1000 NA, N = c Omega, e_need = inherit.e_need(c, V, 1e-3, tau)")
    log(f"    {'setting':17} {'c M':>7} {'tau yr':>8} {'V m3':>9} {'N':>10} {'e_need':>10} {'e_need*sqrtN':>13} {'scorecard':>10}")
    cases_a = {n: case(n, 1.0) for n in NAMES}
    for n, cs in cases_a.items():
        log(f"    {n:17} {cs['c']:7.0e} {cs['tau'] / ar.H / 8766:8.0e} {cs['V']:9.3g} {cs['N']:10.3e} {cs['en']:10.3e}"
            f" {xn(cs['en'], cs):13.2f} {cs['en']:10.0e}")
    verdict("P12", [f"{cs['en']:.0e}" for cs in cases_a.values()] == ["6e-12", "3e-14", "3e-18"],
            "e_need to one digit: " + ", ".join(f"{cs['en']:.0e}" for cs in cases_a.values()))

    # ---- (1) stability and crossing
    log("\n(1) MEASURED RATES (F = 1, s = 6): stability, max margin; CROSSING MULTIPLIER F_c on the grid at C0 = c")
    log(f"    {'setting':17} {'lambda_max /s':>14} {'max margin':>11} {'F_c s=4.0':>11} {'F_c s=6.0':>11} {'F_c s=9.1':>11} {'F_c*C0 mM':>10} {'F_c ext kh':>11}")
    FC = {}
    for n, cs in cases_a.items():
        FC[n] = {s: fcrit(cs["c"], s) for s in (ar.S_LO, 6.0, ar.S_HI)}
        fe = fcrit(cs["c"], 6.0, G_EXT)
        log(f"    {n:17} {cs['b']['l1'].max():14.3e} {cs['b']['m'].max():11.3e} " + " ".join(f"{FC[n][s]:11.4g}" for s in (ar.S_LO, 6.0, ar.S_HI))
            + f" {FC[n][6.0] * cs['c'] * 1e3:10.4g} {fe:11.4g}")
        RES[n] = dict(F_c=FC[n][6.0], F_c_ext=fe)
    ms = [cs["b"]["m"].max() for cs in cases_a.values()]
    verdict("P1", all(cs["b"]["l1"].max() <= 0 for cs in cases_a.values()) and 0.01 <= ms[0] <= 0.08 and 5e-5 <= ms[1] <= 1e-3
            and 1e-8 <= ms[2] <= 1e-5, "lambda_max <= 0 everywhere at F = 1; max margins " + ", ".join(f"{m:.2e}" for m in ms))
    fcs = [FC[n][6.0] for n in NAMES]
    verdict("P2", 150 <= fcs[0] <= 1500 and fcs[1] >= 3e4 and fcs[2] >= 3e6, "F_c(s = 6) = " + ", ".join(f"{v:.4g}" for v in fcs))

    # ---- (2) case (a) feed
    def line(cs, key, sub, extra=""):
        v, i = cs["P"][(key, sub)]
        b = cs["b"]
        return (f"e_min = {v:.3e} (xN {xn(v, cs):9.3g}, /e_need {v / cs['en']:9.3g}) at f = {b['f'][i]:.2e} kh = {b['kh'][i]:.1e} R = {b['R'][i]:.0f}"
                f" K = {b['K'][i]:.0e} cat share {b['cat'][i]:.1e} margin {b['m'][i]:.1e}{extra}")

    log("\n(2) CASE (a) F = 1, FEED excess g (guard 5). e_min = smallest g followed (Delta >= 2); xN = e_min*sqrt(N); letter vs e_need")
    for sub, nm in (("U", "U whole guarded grid"), ("E", "E engaged: margin >= half the max")):
        log(f"  subset {nm}")
        for key, lab in (("feed_a", "monomer readout"), ("feed_b", "dimer readout"), ("feed_lit", "best readout, literal"),
                         ("feed_strict", "best readout, strict (>= 1% out)")):
            ok = [cs["P"][(key, sub)][0] <= cs["en"] for cs in cases_a.values()]
            log(f"    {lab}: letter {letter(ok)} ({''.join('Y' if o else 'N' for o in ok)} for pond, lake, ocean)")
            for n, cs in cases_a.items():
                log(f"      {n:17} " + line(cs, key, sub))
            REC[f"a|{key}|{sub}"] = [rec(cs, key, sub) for cs in cases_a.values()]
    # gains
    log("  gain G = readout excess per feed excess (stable points)")
    gmax = []
    for n, cs in cases_a.items():
        b, A = cs["b"], cs["A"]
        i = np.unravel_index(np.argmax(b["m"]), b["m"].shape)
        gmax.append((np.nanmax(A["Ga"]), A["Ga"][i], np.nanmax(A["Gb"]), A["Gb"][i]))
        log(f"      {n:17} max Ga {gmax[-1][0]:.3f}, Ga at max-margin point {gmax[-1][1]:.4f}; max Gb {gmax[-1][2]:.3f}, Gb at max-margin {gmax[-1][3]:.4f}")
    xs_u = [xn(cs["P"][("feed_a", "U")][0], cs) for cs in cases_a.values()]
    fo = [float(cs["b"]["f"][cs["P"][("feed_a", "U")][1]]) for cs in cases_a.values()]
    verdict("P3", all(0.8 <= v <= 2.2 for v in xs_u) and min(fo) >= 1e-4 and
            letter([cs["P"][("feed_lit", "U")][0] <= cs["en"] for cs in cases_a.values()]) == "Y",
            "U monomer xN = " + ", ".join(f"{v:.3g}" for v in xs_u) + "; f_opt = " + ", ".join(f"{v:.0e}" for v in fo))
    strict_n = letter([cs["P"][("feed_strict", "U")][0] <= cs["en"] for cs in cases_a.values()]) == "N"
    verdict("P4", all(1.0 <= g[0] <= 3.0 and 1.0 <= g[1] <= 1.1 for g in gmax) and strict_n,
            "max Ga " + ", ".join(f"{g[0]:.3f}" for g in gmax) + "; Ga at max margin " + ", ".join(f"{g[1]:.4f}" for g in gmax)
            + f"; strict feed letter {letter([cs['P'][('feed_strict', 'U')][0] <= cs['en'] for cs in cases_a.values()])}")
    xa = [xn(cs["P"][("feed_a", "E")][0], cs) for cs in cases_a.values()]
    xb = [xn(cs["P"][("feed_b", "E")][0], cs) for cs in cases_a.values()]
    le = letter([cs["P"][("feed_lit", "E")][0] <= cs["en"] for cs in cases_a.values()])
    verdict("P5", all(8 <= v <= 400 for v in xa) and all(1 <= v <= 25 for v in xb) and le == "Y",
            "E monomer xN " + ", ".join(f"{v:.3g}" for v in xa) + "; dimer xN " + ", ".join(f"{v:.3g}" for v in xb) + f"; best-readout E letter {le}")

    # ---- (3) case (a) seed + guard sensitivity
    log("\n(3) CASE (a) F = 1, SEED excess ee0 (monomer + dimer both ee0; best readout; literal). Guard = f tau >= G and -lambda_slow tau >= G")
    log(f"    {'guard':>5} " + " ".join(f"{n[:6]+' xN':>11} {'/e_need':>8}" for n in NAMES) + "  letter")
    seedx = {}
    for g in (1, 2, 3, 5, 10):
        cg = {n: case(n, 1.0, guard=float(g)) for n in NAMES}
        vals = [cg[n]["P"][("seed_lit", "U")][0] for n in NAMES]
        seedx[g] = [xn(v, cg[n]) for v, n in zip(vals, NAMES)]
        log(f"    {g:5d} " + " ".join(f"{xn(v, cg[n]):11.4g} {v / cg[n]['en']:8.3g}" for v, n in zip(vals, NAMES))
            + f"  {letter([v <= cg[n]['en'] for v, n in zip(vals, NAMES)])}")
        if g == 5:
            REC["a|seed_lit|U"] = [rec(cg[n], "seed_lit", "U") for n in NAMES]
            REC["a|seed_b|U"] = [rec(cg[n], "seed_b", "U") for n in NAMES]
            REC["a|seed_a|U"] = [rec(cg[n], "seed_a", "U") for n in NAMES]
            REC["a|seed_lit|E"] = [rec(cg[n], "seed_lit", "E") for n in NAMES]
    c5 = cases_a  # F = 1, guard 5
    log("  at guard 5, per readout/subset (inf = no grid point passes the guard)")
    for key, sub in (("seed_a", "U"), ("seed_b", "U"), ("seedm_lit", "U"), ("seed_lit", "E"), ("seedm_lit", "E"), ("seed_strict", "U")):
        log(f"    {key} {sub}: letter " + letter([c5[n]["P"][(key, sub)][0] <= c5[n]["en"] for n in NAMES]))
        for n in NAMES:
            log(f"      {key} {sub} {n:17} " + line(c5[n], key, sub))
    log("  f tau at the best guard-5 seed point (a pool at the guard edge has f tau ~ 5; the f grid floor is 1e-13 /s): " + ", ".join(
        f"{n[:6]} {float(c5[n]['b']['f'][c5[n]['P'][('seed_lit', 'U')][1]]) * c5[n]['tau']:.3g}" for n in NAMES))
    ret = []
    for n in NAMES:
        msk = c5[n]["E"] & np.isfinite(c5[n]["A"]["seed_lit"])
        ret.append(float(np.max(c5[n]["A"]["ret_m"][msk])) if msk.any() else float("nan"))
    log("  monomer-only seed, retained fraction of ee0 at tau, max over engaged guard-passing points (nan = no engaged point passes the guard): "
        + ", ".join(f"{v:.2e}" for v in ret))
    s5 = [xn(c5[n]["P"][("seed_lit", "U")][0], c5[n]) for n in NAMES]
    verdict("P6", all(100 <= v <= 5000 for v in s5) and letter([c5[n]["P"][("seed_lit", "U")][0] <= c5[n]["en"] for n in NAMES]) == "N"
            and all(v < 1e-3 for v in ret) and all(3 <= v <= 40 for v in seedx[1]),  # nan retention => not testable => fails
            "guard-5 xN " + ", ".join(f"{v:.4g}" for v in s5) + "; retention " + ", ".join(f"{v:.1e}" for v in ret)
            + "; guard-1 xN " + ", ".join(f"{v:.3g}" for v in seedx[1]))

    # ---- (4) case (b)
    def bcase(mult, s=6.0, grid=ar.GRID):
        res = {}
        for n in NAMES:
            Fc = fcrit(SETTINGS[n][0], s, grid)
            res[n] = (Fc, case(n, mult * Fc, s, grid) if np.isfinite(Fc) else None)
        return res

    log("\n(4) CASE (b) F = 1.5 F_c (s = 6; F_c = inf => not met). lambda_max, lambda*tau, Delta-only (literal) and strict e_min, saturation")
    bb = bcase(1.5)
    sat_ee = {}
    BOK = {}
    for variant in ("feed", "seed"):
        for kind in ("lit", "strict"):
            key = f"{variant}_{kind}"
            ok = [(bb[n][1] is not None) and bb[n][1]["P"][(key, "U")][0] <= bb[n][1]["en"] for n in NAMES]
            log(f"  {variant} {kind}: letter {letter(ok)} ({''.join('Y' if o else 'N' for o in ok)} for pond, lake, ocean; before the ee_sat >= 1% test)")
            for n in NAMES:
                Fc, cs = bb[n]
                if cs is None:
                    log(f"      {n:17} F_c = inf: no crossing by F = 1e8 on the grid")
                    continue
                v, i = cs["P"][(key, "U")]
                lt = float(cs["b"]["l1"][i] * cs["tau"])
                se = sat(cs, key, "U") if cs["b"]["l1"][i] > 0 and np.isfinite(v) else (float("nan"), "stable point")
                sat_ee[(n, key)] = se
                BOK[(n, key)] = v <= cs["en"] and (kind == "lit" or (np.isfinite(se[0]) and se[0] >= 0.01))
                log(f"      {n:17} F_c {Fc:.4g}, F_b {cs['b']['F']:.4g}; " + line(cs, key, "U",
                    f"; lambda {cs['b']['l1'][i]:.2e} lambda*tau {lt:.3g}; ee_sat {se[0]:.3f} ({se[1]})"))
            REC[f"b|{key}|U"] = [rec(bb[n][1], key, "U") if bb[n][1] is not None else None for n in NAMES]
    log("  UNSTABLE POINTS ONLY (lambda > 0 at the chosen parameters), subset U: e_min xN and /e_need per variant")
    for key in ("u_feed", "u_feedS", "u_seed", "u_seedS", "u_seedm", "u_seedmS"):
        log(f"    {key:9} " + "; ".join(f"{n[:6]} xN {xn(bb[n][1]['P'][(key, 'U')][0], bb[n][1]):.3g} ({bb[n][1]['P'][(key, 'U')][0] / bb[n][1]['en']:.3g}x e_need)"
                                       if bb[n][1] else f"{n[:6]} none" for n in NAMES))
    log("  lambda_max at F_b: " + ", ".join(f"{n} {bb[n][1]['b']['l1'].max():.3e}" if bb[n][1] else f"{n} none" for n in NAMES))
    fin = [n for n in NAMES if bb[n][1] is not None]
    dx = {n: xn(bb[n][1]["P"][("feed_lit", "U")][0], bb[n][1]) for n in fin}
    verdict("P7", bb[NAMES[0]][1] is not None and bb[NAMES[0]][1]["b"]["l1"].max() > 0
            and all(0.5 <= dx[n] <= 100 and bb[n][1]["P"][("feed_lit", "U")][0] <= bb[n][1]["en"] for n in fin),
            "literal feed xN " + ", ".join(f"{n} {dx[n]:.3g}" for n in fin) + (" (F_c = inf: " + ", ".join(n for n in NAMES if n not in fin) + ")" if len(fin) < 3 else ""))
    pond = bb[NAMES[0]][1]
    ok_s = [(bb[n][1] is not None) and bb[n][1]["P"][("feed_strict", "U")][0] <= bb[n][1]["en"] for n in NAMES]
    pi = pond["P"][("feed_strict", "U")][1] if pond else None
    plt_ = float(pond["b"]["l1"][pi] * pond["tau"]) if pond else float("nan")
    verdict("P8", pond is not None and plt_ < 23 and pond["P"][("feed_strict", "U")][0] > 1e3 * pond["en"] and letter(ok_s) in ("N", "P"),
            f"pond strict feed lambda*tau {plt_:.3g}, e_min/e_need {pond['P'][('feed_strict', 'U')][0] / pond['en']:.3g}; strict (b) feed letter {letter(ok_s)}")
    sv = [v[0] for v in sat_ee.values() if np.isfinite(v[0])]
    verdict("P9", len(sv) > 0 and all(0.02 <= abs(v) <= 1 for v in sv), "ee_sat at (b) points: " + ", ".join(f"{v:.3f}" for v in sv))

    # ---- (5) sensitivities
    log("\n(5) SENSITIVITY of (b): multiplier on F_c, s, and kh up to 1e4 (feed variant, subset U). xN = e_min*sqrt(N); strict gate = 1% growth in tau")
    log(f"    {'case':26} {'setting':17} {'F_c':>10} {'lam_max':>10} {'lam*tau':>8} {'lit xN':>9} {'lit/e_need':>10} {'strict xN':>10} {'strict/e_need':>13}")
    for lab, mult, s, grid in (("mult 1.1", 1.1, 6.0, ar.GRID), ("mult 1.5 (base)", 1.5, 6.0, ar.GRID), ("mult 3", 3.0, 6.0, ar.GRID),
                               ("mult 10", 10.0, 6.0, ar.GRID), ("s = 4.0, mult 1.5", 1.5, ar.S_LO, ar.GRID),
                               ("s = 9.1, mult 1.5", 1.5, ar.S_HI, ar.GRID), ("kh to 1e4, mult 1.5", 1.5, 6.0, G_EXT)):
        bs = bcase(mult, s, grid)
        for n in NAMES:
            Fc, cs = bs[n]
            if cs is None:
                log(f"    {lab:26} {n:17} {Fc:10.4g} {'no crossing':>10}")
                continue
            vl, il = cs["P"][("feed_lit", "U")]
            vs, is_ = cs["P"][("feed_strict", "U")]
            log(f"    {lab:26} {n:17} {Fc:10.4g} {cs['b']['l1'].max():10.2e} {cs['b']['l1'][is_] * cs['tau']:8.3g} {xn(vl, cs):9.3g} {vl / cs['en']:10.3g}"
                f" {xn(vs, cs):10.3g} {vs / cs['en']:13.3g}")
    log("\n    (a) sensitivity to s (F = 1, feed, guard 5): xN of U best readout / E best readout / strict U")
    for s in (ar.S_LO, 6.0, ar.S_HI):
        cs_ = {n: case(n, 1.0, s) for n in NAMES}
        log(f"      s = {s:.1f}: " + "; ".join(f"{n} U {xn(cs_[n]['P'][('feed_lit', 'U')][0], cs_[n]):.5g} E {xn(cs_[n]['P'][('feed_lit', 'E')][0], cs_[n]):.5g}"
                                              f" strict {cs_[n]['P'][('feed_strict', 'U')][0]:.5e}" for n in NAMES))

    # ---- (5b) POST HOC: f grid extended down to 1e-19 /s (the ocean seed needs f tau ~ 5, f ~ 1.6e-15, below the 1e-13 floor)
    log("\n(5b) POST HOC (not pre-registered): f grid extended to 1e-19 /s (4 per decade). Only the f floor changes; kh and the rest as the base grid.")
    log(f"    {'quantity':34} " + " ".join(f"{n[:6]:>13}" for n in NAMES))
    for g in (1, 5, 10):
        cg = {n: case(n, 1.0, 6.0, G_FLOW, float(g)) for n in NAMES}
        log(f"    {'(a) seed, best readout, guard ' + str(g) + ' xN':34} " + " ".join(f"{xn(cg[n]['P'][('seed_lit', 'U')][0], cg[n]):13.4g}" for n in NAMES))
        log(f"    {'    ratio to e_need':34} " + " ".join(f"{cg[n]['P'][('seed_lit', 'U')][0] / cg[n]['en']:13.4g}" for n in NAMES))
        if g == 5:
            log(f"    {'    f at the optimum (/s)':34} " + " ".join(f"{float(cg[n]['b']['f'][cg[n]['P'][('seed_lit', 'U')][1]]):13.3e}" for n in NAMES))
            log(f"    {'    letter (literal seed, guard 5)':34} " + letter([cg[n]['P'][('seed_lit', 'U')][0] <= cg[n]['en'] for n in NAMES]))
    bx = bcase(1.5, 6.0, G_FLOW)
    log(f"    {'(b) F_c (s = 6)':34} " + " ".join(f"{bx[n][0]:13.4g}" for n in NAMES))
    for key in ("feed_lit", "feed_strict", "seed_lit", "seed_strict"):
        log(f"    {'(b) ' + key + ' xN':34} " + " ".join(f"{xn(bx[n][1]['P'][(key, 'U')][0], bx[n][1]):13.4g}" if bx[n][1] else f"{'none':>13}" for n in NAMES))
        log(f"    {'    ratio to e_need':34} " + " ".join(f"{bx[n][1]['P'][(key, 'U')][0] / bx[n][1]['en']:13.4g}" if bx[n][1] else f"{'none':>13}" for n in NAMES))

    # ---- letters
    log("\n(6) IMPLIED R6 LETTERS (Y = all three settings met, P = some, N = none), pond / lake / ocean")
    L = lambda key, sub, cd: (letter([cd[n]["P"][(key, sub)][0] <= cd[n]["en"] for n in NAMES]),
                              "".join("Y" if cd[n]["P"][(key, sub)][0] <= cd[n]["en"] else "N" for n in NAMES))
    for lab, key, sub, cd in (("(a) literal  feed  U (EXISTING RULE)", "feed_lit", "U", cases_a), ("(a) literal  feed  E", "feed_lit", "E", cases_a),
                              ("(a) strict   feed  U (scorecard comment)", "feed_strict", "U", cases_a), ("(a) strict   feed  E", "feed_strict", "E", cases_a),
                              ("(a) literal  seed  U (guard 5)", "seed_lit", "U", c5), ("(a) literal  seed  E", "seed_lit", "E", c5),
                              ("(a) strict   seed  U", "seed_strict", "U", c5), ("(a) literal  seed  U, monomer-only seed", "seedm_lit", "U", c5)):
        log(f"    {lab:44} {L(key, sub, cd)[0]}  ({L(key, sub, cd)[1]})")
    bc = {n: bb[n][1] for n in NAMES if bb[n][1] is not None}
    for lab, key in (("(b) literal  feed  U", "feed_lit"), ("(b) strict   feed  U", "feed_strict"), ("(b) literal  seed  U", "seed_lit"),
                     ("(b) strict   seed  U", "seed_strict"), ("(b) literal  seed  U, monomer-only seed", "seedm_lit")):
        ok = [(n in bc) and bc[n]["P"][(key, "U")][0] <= bc[n]["en"] and (not key.endswith("strict") or BOK.get((n, key), False)) for n in NAMES]
        log(f"    {lab:44} {letter(ok)}  ({''.join('Y' if o else 'N' for o in ok)}; F_c = inf counts as not met; strict also needs ee_sat >= 1%)")
    log("  Note: inherit.py's inherited excess e0 is a SEED in the pool at the transition (inherit.py:3-5); the feed variant is the continuous-bias analogue.")

    log("\nPREDICTIONS vs RESULT (P10 in arm_e_r6_analytic.log, P11 in arm_e_r6_langevin.log)")
    for tag, ok in chk:
        log(f"  {tag}: {'held' if ok else 'FAILED'}")
    json.dump(dict(sha=sha, records=REC, F_c={n: {str(s): FC[n][s] for s in FC[n]} for n in NAMES}), open(HERE / "arm_e_r6_results.json", "w"), indent=1, default=str)
    log(f"\nresults saved to arm_e_r6_results.json; runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
