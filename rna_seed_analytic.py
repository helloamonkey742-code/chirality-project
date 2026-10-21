"""Arm G, independent analytic cross-check of rna_copying.py: linearise the per-cycle map about the racemic state.

Model (equations of rna_copying.py's docstring, re-derived here, no code taken from it for the formulas):
    attempt(mr, mw): tot = mr + r mw, p = mr/tot, mu = kap tot, kap = kobs t_cold / 24 mM, N ~ Poisson(mu)
        S = p^n P(N>=n),  Q = sum_{j<=n} p^(j-1) P(N>=j),  E_R = p Q,  E_W = (1-p) Q   (so E_R/E_W = p/(1-p) exactly)
    T_D' = (1-w) T_D (1 + c S_D) + w T_in/2
    m_D' = (1-w) [(1-h) m_D - c (T_D E_R,D + T_L E_W,L)] + w m_in/2        (L: swap D, L)
Racemic fixed point: T(m) = (w T_in/2) / (1 - (1-w)(1 + c S(m))), m from the monomer balance (1-D root, brentq).

Difference mode. Perturb (dT_D, dT_L, dm_D, dm_L) = (a, -a, b, -b). With p0 = 1/(1+r), mu0 = kap m (1+r):
    dp/db = p_b = 2r / (m (1+r)^2),   dmu/db = mu_b = kap (1-r)
    sigma = dS_D/db     = n p0^(n-1) P_n p_b + p0^n pi_n mu_b                    pi_j = e^-mu mu^(j-1)/(j-1)!
    eps   = d(E_R-E_W)_D/db = 2 Q p_b + (2 p0 - 1)(Q_p p_b + Q_mu mu_b),  Q_p = sum (j-1) p0^(j-2) P_j,  Q_mu = sum p0^(j-1) pi_j
    J = [[A, B], [C, Dd]] on (a, b):
        A  = (1-w)(1 + c S) = 1 - w T_in / (2 T)            (< 1 for T_in > 0)
        B  = (1-w) c T sigma
        C  = -(1-w) c (E_R - E_W) = -(1-w) c (2 p0 - 1) Q    ((2p0-1) = (1-r)/(1+r))
        Dd = (1-w) [(1-h) - c T eps]
(a) eigenvalues  lam = [tr +- sqrt(tr^2 - 4 det)]/2, tr = A + Dd, det = A Dd - B C. Excess T_D - T_L decays as lam^t
    (complex pair = damped swing, relaxation -1/ln|lam| cycles).
(b) feed excess e (m_D,in = m_in (1+e)/2): forcing f = (0, w m_in e/2) on (a, b); steady state (I-J)(a, b) = f:
        Delta = (1-A)(1-Dd) - B C = det(I - J);   b = f_b (1-A)/Delta;   a = f_b B/Delta
    product excess ee = (F_D - F_L)/(F_D + F_L) = a/T + sigma b/S  -> gain G = ee/e (linear in e; diverges as Delta -> 0).
    A seeded TEMPLATE excess is not a forcing: it decays at lam, so its final gain is 0 whenever |lam| < 1.
(c) racemic state loses stability when Delta = 0 (lam = 1). For r <= 1: B > 0, C <= 0, Dd < 1, A < 1  =>  Delta > 0, stable
    at every feed (r = 1: C = 0 exactly, lam = A, Dd). r > 1: C > 0, crossing possible only if racemic copying outruns
    washout, c S(m) = w/(1-w) reachable; with S <= p0^n the window is  1 < r < r_max = (c (1-w)/w)^(1/n) - 1.
    Threshold feed m*(r): root of Delta(m_in) = 0; m* -> infinity at both window edges.

Constants (rna_copying.py lines 84-88): k4 3.33e-3 /h, kobs 7.35 /h, 24 mM, T_in 3 uM, t_cold 1 h; n 3, c 0.5, w 0.02.

PRE-REGISTERED PREDICTIONS (formulas only; computed with `python3 rna_seed_analytic.py --predict` BEFORE importing,
running or reading any output of rna_copying.py; the target numbers 0.67 / 0.98 / 10.1 mM were known from the paper text but
no parameter or formula was tuned to them). Saved verbatim in rna_seed_analytic_predictions.txt.
    file sha256 fa679df67e0c2be995e77409fc11c4c4ebfbcd4c8db21aa9c9fc21af7d2359fe, written Mon Oct 19 23:38:59 PDT 2026
    feed 20 mM (n = 3, c = 0.5, w = 0.02):
      r      T uM    m mM    A         B         C          Dd        lam1       lam2      cycles/e-fold
      0.40   287.94  2.700  0.999896  0.005627  -0.204227  0.923279  0.979430  0.943745  48.1
      0.56   260.30  2.821  0.999885  0.005020  -0.145608  0.926637  0.987967  0.938554  82.6
      0.98   205.64  3.245  0.999854  0.003758  -0.006117  0.938220  0.999479  0.938596  1918.2
      1.14   188.28  3.468  0.999841  0.003335  +0.041462  0.943310  1.002190  0.940962  growth (lam > 1)
    (a) leading eigenvalue is real at all four r; |lam| < 1 at r = 0.40, 0.56, 0.98 and > 1 (1.00219) at r = 1.14.
        max |lam| over feeds 0.5-100 mM (400 log points): 0.995479 (r 0.40), 0.996293 (r 0.56), 0.999479 (r 0.98) -> stable at every feed for r <= 0.98.
    (b) linear gain G = ee/e (valid for G e << 1) and exact nonlinear plateau at e = +0.2, feed-filled start:
      r      G@5 mM   G@20 mM   0.2 G@20   nonlinear ee@20 (5000 cycles = steady)   nonlinear ee@5
      0.40    5.928    3.395     0.679     0.675183 (0.675183)              0.854350
      0.56    5.693    5.247     1.049     0.951332 (0.951332)              0.843220
      0.98    5.252  115.054    23.011     0.997446 (0.997446)              0.811560
      1.14    5.135    nan      nan           nan (nan)              0.800823
        r = 1.14 at 20 mM is past threshold (no steady state, gain n/a). Template seeds: final gain 0 for every |lam| < 1.
    (c) r_max = 1.90439; threshold feed m*(1.14) = 10.0878 mM (r = 0.40, 0.56, 0.98: none).
        m*(r) mM: 1.02 14.312, 1.05 11.510, 1.10 10.402, 1.14 10.088, 1.20 9.909, 1.50 10.569, 1.80 13.174, 1.90 15.523.
        r_crit at 10 / 20 / 100 mM feed: 1.1605 / 1.0080 / 1.0007.
    Success criteria fixed in advance: matrix elements and eigenvalues within 1e-6 relative of rna_copying.mode; threshold within 1e-4 mM;
    ee within 1e-9 (same 5000-cycle map) and 1e-6 (steady state); stated 0.67 (r = 0.40 is rounded) and 10.1 mM (0.05 mM) within rounding.
Comparison run (python3 rna_seed_analytic.py > rna_seed_analytic.log) re-evaluates the formulas, asserts they reproduce the
file above (and that its sha256 is the one recorded here), then imports rna_copying and compares.
"""
import cmath
import hashlib
import math
import sys
import numpy as np
from scipy.optimize import brentq, fsolve
from scipy.special import gammainc

K4, KOBS, M_REF, T_IN, T_COLD = 3.33e-3, 7.35, 24.0, 0.003, 1.0
N, C0, W = 3, 0.5, 0.02
H = -math.expm1(-K4 * T_COLD)
KAP = KOBS * T_COLD / M_REF
R_LIST, FEED, FEEDS = (0.40, 0.56, 0.98, 1.14), 20.0, np.geomspace(0.5, 100, 400)
PREREG_SHA256 = "fa679df67e0c2be995e77409fc11c4c4ebfbcd4c8db21aa9c9fc21af7d2359fe"
PRED_FILE = "rna_seed_analytic_predictions.txt"
R_MAX = (C0 * (1 - W) / W) ** (1 / N) - 1


def P(j, mu):
    """P(Poisson(mu) >= j)."""
    return gammainc(j, mu)


def pi(j, mu):
    """Poisson pmf at j-1 = dP(j, mu)/dmu."""
    return math.exp(-mu + (j - 1) * math.log(mu) - math.lgamma(j)) if mu > 0 else 0.0


def sym(r, m):
    """Sums at the racemic point: p0, mu0, S, Q, Q_p, Q_mu and the S derivatives."""
    p, mu = 1 / (1 + r), KAP * (1 + r) * m
    Pj = [P(j, mu) for j in range(1, N + 1)]
    pj = [pi(j, mu) for j in range(1, N + 1)]
    Q = sum(p ** (j - 1) * Pj[j - 1] for j in range(1, N + 1))
    Qp = sum((j - 1) * p ** (j - 2) * Pj[j - 1] for j in range(2, N + 1))
    Qmu = sum(p ** (j - 1) * pj[j - 1] for j in range(1, N + 1))
    return dict(p=p, mu=mu, S=p ** N * Pj[-1], Q=Q, Qp=Qp, Qmu=Qmu,
                Sp=N * p ** (N - 1) * Pj[-1], Smu=p ** N * pj[-1])


def fixed(r, m_in):
    """Racemic fixed point (T, m)."""
    grow = lambda m: (1 - W) * (1 + C0 * sym(r, m)["S"]) - 1
    hi = W * m_in / 2 / (1 - (1 - W) * (1 - H))
    if grow(hi) > 0:
        hi = brentq(grow, 1e-12, hi)
    T = lambda m: W * T_IN / 2 / -grow(m)
    f = lambda m: (1 - W) * ((1 - H) * m - C0 * T(m) * sym(r, m)["Q"]) + W * m_in / 2 - m
    m = brentq(f, 1e-9 * hi, hi * (1 - 1e-12), xtol=1e-14)
    return T(m), m


def lin(r, m_in):
    """A, B, C, Dd and the derivatives sigma, eps at the racemic fixed point."""
    T, m = fixed(r, m_in)
    s = sym(r, m)
    pb, mb = 2 * r / (m * (1 + r) ** 2), KAP * (1 - r)
    sig = s["Sp"] * pb + s["Smu"] * mb
    eps = 2 * s["Q"] * pb + (2 * s["p"] - 1) * (s["Qp"] * pb + s["Qmu"] * mb)
    A = 1 - W * T_IN / (2 * T)
    B = (1 - W) * C0 * T * sig
    C = -(1 - W) * C0 * (2 * s["p"] - 1) * s["Q"]
    D = (1 - W) * ((1 - H) - C0 * T * eps)
    return dict(T=T, m=m, A=A, B=B, C=C, D=D, sig=sig, eps=eps, S=s["S"])


def eig(L):
    tr, det = L["A"] + L["D"], L["A"] * L["D"] - L["B"] * L["C"]
    q = cmath.sqrt(tr * tr - 4 * det)
    return (tr + q) / 2, (tr - q) / 2


def delta(L):
    return (1 - L["A"]) * (1 - L["D"]) - L["B"] * L["C"]


def gain(r, m_in):
    """Linear product-excess gain ee/e for a feed excess e (None when Delta <= 0: no steady state)."""
    L = lin(r, m_in)
    d = delta(L)
    if d <= 0:
        return None
    fb = W * m_in / 2
    a, b = L["B"] * fb / d, fb * (1 - L["A"]) / d
    return a / L["T"] + L["sig"] * b / L["S"]


def mapF(z, m_in, r, e=0.0):
    """Independent re-implementation of the per-cycle map (for the nonlinear steady state and a finite-difference check)."""
    def att(mr, mw):
        tot = mr + r * mw
        p, mu = mr / tot, KAP * tot
        Q = sum(p ** (j - 1) * P(j, mu) for j in range(1, N + 1))
        return p ** N * P(N, mu), p * Q, (1 - p) * Q
    TD, TL, mD, mL = z
    SD, RD, WD = att(mD, mL)
    SL, RL, WL = att(mL, mD)
    z2 = np.array([(1 - W) * TD * (1 + C0 * SD) + W * T_IN / 2, (1 - W) * TL * (1 + C0 * SL) + W * T_IN / 2,
                   (1 - W) * ((1 - H) * mD - C0 * (TD * RD + TL * WL)) + W * m_in * (1 + e) / 2,
                   (1 - W) * ((1 - H) * mL - C0 * (TL * RL + TD * WD)) + W * m_in * (1 - e) / 2])
    return z2, C0 * TD * SD, C0 * TL * SL


def forced(r, m_in, e, cap=400000):
    """Nonlinear forced response from the feed-filled reactor (as in rna_copying's run): iterate the own map until the
    per-cycle relative change is < 1e-15. Returns ee after 5000 cycles, steady-state ee, cycles used."""
    z = np.array([T_IN / 2, T_IN / 2, m_in * (1 + e) / 2, m_in * (1 - e) / 2])
    ee5000 = float("nan")
    for t in range(1, cap + 1):
        z2, FD, FL = mapF(z, m_in, r, e)
        ee = (FD - FL) / (FD + FL)
        if t == 5000:
            ee5000 = ee
        if t > 5000 and np.abs(z2 / z - 1).max() < 1e-15:
            return ee5000, ee, t
        z = z2
    return ee5000, ee, cap


def threshold(r, grid=FEEDS):
    """Feed where Delta crosses zero from + to - (None if no crossing on the grid)."""
    v = [delta(lin(r, m)) for m in grid]
    for i in range(1, len(grid)):
        if v[i - 1] > 0 >= v[i]:
            return brentq(lambda m: delta(lin(r, m)), grid[i - 1], grid[i], xtol=1e-9)
    return None


def r_crit(m_in):
    """Smallest r in (1, r_max) where Delta(m_in) < 0 first appears, by bracketing + brentq (None if none)."""
    rs = np.linspace(1.0005, R_MAX - 1e-3, 400)
    v = [delta(lin(r, m_in)) for r in rs]
    for i in range(1, len(rs)):
        if v[i - 1] > 0 >= v[i]:
            return brentq(lambda r: delta(lin(r, m_in)), rs[i - 1], rs[i], xtol=1e-10)
    return None


def predict():
    out = {}
    for r in R_LIST:
        L = lin(r, FEED)
        l1, l2 = eig(L)
        big = max(abs(l1), abs(l2))
        k = f"r{r:.2f}."
        out |= {k + "T_uM@20": L["T"] * 1e3, k + "m_mM@20": L["m"], k + "A@20": L["A"], k + "B@20": L["B"],
                k + "C@20": L["C"], k + "Dd@20": L["D"], k + "lam1_re@20": l1.real, k + "lam1_im@20": abs(l1.imag),
                k + "lam2_re@20": l2.real, k + "absmax@20": big,
                k + "cycles_per_efold@20": -1 / math.log(big) if big < 1 else float("nan")}
        lams = [eig(lin(r, m)) for m in FEEDS]
        out[k + "absmax_allfeeds"] = max(max(abs(a), abs(b)) for a, b in lams)
        out[k + "maxRe_allfeeds"] = max(max(a.real, b.real) for a, b in lams)
        ms = threshold(r)
        out[k + "thr_feed_mM"] = float("nan") if ms is None else ms
        for m_in in (5.0, FEED):
            g = gain(r, m_in)
            out[k + f"gain_lin@{m_in:g}"] = float("nan") if g is None else g
            out[k + f"ee_lin_e0.2@{m_in:g}"] = float("nan") if g is None else 0.2 * g
            f5, fss, _ = (float("nan"),) * 2 + (0,) if g is None else forced(r, m_in, 0.2)
            out[k + f"ee_nl_5000cyc_e0.2@{m_in:g}"] = f5
            out[k + f"ee_nl_steady_e0.2@{m_in:g}"] = fss
    out["r_max"] = R_MAX
    out["Delta_at_r1.000_20mM"] = delta(lin(1.0, FEED))
    for r in (1.02, 1.05, 1.10, 1.14, 1.20, 1.50, 1.80, 1.90):
        ms = threshold(r, np.geomspace(0.5, 1e5, 600))
        out[f"mstar_r{r:.2f}"] = float("nan") if ms is None else ms
    for m_in in (10.0, 20.0, 100.0):
        rc = r_crit(m_in)
        out[f"r_crit@{m_in:g}"] = float("nan") if rc is None else rc
    return out


def fmt(out):
    return "\n".join(f"{k} = {v:.6g}" for k, v in out.items())


def read_pred():
    txt = open(PRED_FILE).read()
    return {a.strip(): float(b) for a, b in (ln.split(" = ") for ln in txt.splitlines() if ln and ln[0] != "#")}


def close(x, y, rtol=1e-5, atol=1e-9):
    return (math.isnan(x) and math.isnan(y)) or abs(x - y) <= atol + rtol * abs(y)


def sympy_check():
    """Hand-derived sigma, eps vs sympy differentiation of S and E_R-E_W with explicit Poisson tails (n = 3)."""
    import sympy as sp
    mr, mw, r, kap, b, m = sp.symbols("mr mw r kap b m", positive=True)
    tot = mr + r * mw
    p, mu = mr / tot, kap * tot
    tail = lambda j: 1 - sp.exp(-mu) * sum(mu ** k / sp.factorial(k) for k in range(j))
    Q = sum(p ** (j - 1) * tail(j) for j in range(1, N + 1))
    S, Ed, ER, EW = p ** N * tail(N), (2 * p - 1) * Q, p * Q, (1 - p) * Q
    assert sp.simplify(ER / EW - p / (1 - p)) == 0, "E_R/E_W identity"
    sub = {mr: m + b, mw: m - b}
    val = {r: 0.56, kap: KAP, m: 3.7}
    sig = sp.diff(S.subs(sub), b).subs(b, 0).subs(val).evalf(20)
    eps = sp.diff(Ed.subs(sub), b).subs(b, 0).subs(val).evalf(20)
    s = sym(0.56, 3.7)
    pb, mb = 2 * 0.56 / (3.7 * 1.56 ** 2), KAP * 0.44
    sig_h = s["Sp"] * pb + s["Smu"] * mb
    eps_h = 2 * s["Q"] * pb + (2 * s["p"] - 1) * (s["Qp"] * pb + s["Qmu"] * mb)
    assert abs(float(sig) / sig_h - 1) < 1e-10 and abs(float(eps) / eps_h - 1) < 1e-10, "sympy vs hand derivatives"
    return float(sig) / sig_h - 1, float(eps) / eps_h - 1


def fd_jacobian_check(r, m_in):
    """Antisymmetric block of a finite-difference Jacobian of mapF (central, 1e-6 relative) vs analytic A, B, C, Dd."""
    L = lin(r, m_in)
    z = np.array([L["T"], L["T"], L["m"], L["m"]])
    PM = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]]) / math.sqrt(2)
    J = np.empty((4, 4))
    for i in range(4):
        d = np.zeros(4)
        d[i] = 1e-6 * z[i]
        J[:, i] = (mapF(z + d, m_in, r)[0] - mapF(z - d, m_in, r)[0]) / (2 * d[i])
    Jd = PM.T @ J @ PM
    ana = np.array([[L["A"], L["B"]], [L["C"], L["D"]]])
    return np.abs(Jd / ana - 1).max()


def main():
    if "--predict" in sys.argv:
        print(fmt(predict()))
        return
    sha = hashlib.sha256(open(PRED_FILE, "rb").read()).hexdigest()
    print(f"pre-registration file {PRED_FILE}: sha256 {sha}")
    assert sha == PREREG_SHA256, "pre-registration file changed after hashing"
    pred = read_pred()
    live = predict()
    bad = [k for k in pred if not close(live[k], pred[k], rtol=2e-5)]
    assert not bad and set(pred) == set(live), f"formulas no longer reproduce the pre-registered values: {bad}"
    dev = max(abs(live[k] - pred[k]) / abs(pred[k]) for k in pred if pred[k] and not math.isnan(pred[k]))
    print(f"live formulas reproduce all {len(pred)} pre-registered values (max rel deviation {dev:.1e}: the file keeps 6 significant digits)")
    pred = live          # the file is 6-digit text; comparisons below use the same formulas at full precision
    ds, de = sympy_check()
    print(f"sympy check: hand-derived sigma and eps vs sympy differentiation, rel diff {ds:.1e}, {de:.1e}; E_R/E_W = p/(1-p) simplifies to 0")
    for r in (0.40, 1.14):
        assert fd_jacobian_check(r, FEED) < 1e-6, "analytic block vs finite-difference of own map"
    print("analytic A, B, C, Dd match a finite-difference Jacobian of the independent map (r = 0.40, 1.14, 20 mM) to < 1e-6 relative")

    print("\n(a)(b) PREDICTIONS (feed 20 mM; ee gain = product excess per unit feed excess e, linear; nl = exact forced steady state)")
    hdr = f"{'r':>5} {'T uM':>7} {'m mM':>7} {'A':>8} {'B':>8} {'C':>8} {'Dd':>8} {'lam1':>17} {'|lam2|':>7} {'cyc/e-fold':>10} " \
          f"{'G@20':>7} {'ee_lin':>7} {'ee_nl':>7} {'G@5':>7}"
    print(hdr)
    for r in R_LIST:
        k = f"r{r:.2f}."
        l1 = complex(pred[k + "lam1_re@20"], pred[k + "lam1_im@20"])
        print(f"{r:5.2f} {pred[k + 'T_uM@20']:7.3f} {pred[k + 'm_mM@20']:7.3f} {pred[k + 'A@20']:8.5f} {pred[k + 'B@20']:8.5f} "
              f"{pred[k + 'C@20']:8.5f} {pred[k + 'Dd@20']:8.5f} {l1.real:8.5f}{l1.imag:+8.5f}j {abs(pred[k + 'lam2_re@20']):7.4f} "
              f"{pred[k + 'cycles_per_efold@20']:10.2f} {pred[k + 'gain_lin@20']:7.3f} {pred[k + 'ee_lin_e0.2@20']:7.4f} "
              f"{pred[k + 'ee_nl_steady_e0.2@20']:7.4f} {pred[k + 'gain_lin@5']:7.3f}")

    # structural proof for r <= 1 on a grid
    for r in (0.02, 0.1, 0.4, 0.56, 0.98, 1.0):
        for m in FEEDS[::20]:
            L = lin(r, m)
            assert L["A"] < 1 and L["D"] < 1 and L["B"] > 0 and L["C"] <= 0 and delta(L) > 0, (r, m, L)
    assert abs(lin(1.0, FEED)["C"]) == 0 or abs(lin(1.0, FEED)["C"]) < 1e-15
    print("\n(c) r <= 1 (r = 0.02-1.0, 20 feeds each, 0.5-100 mM): A < 1, Dd < 1, B > 0, C <= 0, Delta > 0 everywhere "
          "(lam = 1 impossible); r = 1: C = 0 exactly, lam = (A, Dd)")
    print(f"    window edges: r = 1 (C = 0) and r_max = (c(1-w)/w)^(1/n) - 1 = {R_MAX:.4f}")
    print("    threshold feed m*(r) (Delta = 0), feed grid 0.5-1e5 mM:")
    for r in (1.02, 1.05, 1.10, 1.14, 1.20, 1.50, 1.80, 1.90):
        print(f"      r = {r:.2f}: m* = {pred[f'mstar_r{r:.2f}']:.4g} mM")
    rs = np.linspace(1.01, 1.89, 89)
    ms = [threshold(r, np.geomspace(0.5, 1e5, 600)) for r in rs]
    ms = [(m, r) for m, r in zip(ms, rs) if m is not None]
    print(f"    minimum of m*(r) on r = 1.01-1.89: {min(ms)[0]:.3f} mM at r = {min(ms)[1]:.2f} (m* -> infinity at both window edges)")
    for m_in in (10.0, 20.0, 100.0):
        print(f"    r_crit at feed {m_in:g} mM: {pred[f'r_crit@{m_in:g}']:.4f}")

    # ---- comparison with rna_copying.py (imported only now) ----
    import rna_copying as rc
    consts = dict(K4=K4, KOBS=KOBS, M_REF=M_REF, T_IN=T_IN, T_COLD=T_COLD)
    assert all(getattr(rc, a) == v for a, v in consts.items()), "constants differ from rna_copying.py"
    assert (rc.LO["n"], rc.LO["c"], rc.LO["w"]) == (N, C0, W) and abs(rc.H - H) < 1e-15 and abs(rc.A / rc.M_REF - KAP) < 1e-15
    print("\nCOMPARISON with rna_copying.py (imported). tolerances (pre-registered): matrix/eig 1e-6 rel, "
          "threshold 1e-4 mM abs, ee 1e-9 abs (5000 cycles), 1e-6 abs (steady)")
    print(f"{'r':>5} {'quantity':>22} {'analytic':>12} {'rna_copying':>12} {'diff':>10} {'tol':>8}  ok")
    fails, rows = [], []

    def cmp(r, name, a, b, tol, rel=True):
        d = abs(a - b) / abs(b) if rel and b else abs(a - b)
        ok = d <= tol
        rows.append(ok)
        if not ok:
            fails.append((r, name, a, b, d))
        print(f"{r:5.2f} {name:>22} {a:12.7f} {b:12.7f} {d:10.2e} {tol:8.0e}  {'yes' if ok else 'NO'}")

    for r in R_LIST:
        k = dict(rc.LO, r=r)
        lam, Jd = rc.mode(rc.fixed_point(FEED, **k), FEED, k)
        for name, v, key in (("A", Jd[0, 0], "A@20"), ("B", Jd[0, 1], "B@20"), ("C", Jd[1, 0], "C@20"), ("Dd", Jd[1, 1], "Dd@20")):
            cmp(r, name, pred[f"r{r:.2f}.{key}"], v, 1e-6)
        cmp(r, "|lam| max", pred[f"r{r:.2f}.absmax@20"], abs(lam).max(), 1e-6)
        ms_rc, mods, dd = rc.threshold(k)
        a_ms = pred[f"r{r:.2f}.thr_feed_mM"]
        if r < 1:
            ok = ms_rc is None and math.isnan(a_ms)
            rows.append(ok)
            fails += [] if ok else [(r, "stable all feeds", a_ms, ms_rc, 0)]
            print(f"{r:5.2f} {'stable at all feeds':>22} {'yes' if math.isnan(a_ms) else a_ms:>12} {'yes' if ms_rc is None else ms_rc:>12} "
                  f"{'':>10} {'':>8}  {'yes' if ok else 'NO'}")
            grid = np.geomspace(0.5, 100, 40)
            mine = max(max(abs(x) for x in eig(lin(r, m))) for m in grid)
            cmp(r, "max|lam| 40-pt grid", mine, mods, 1e-6)
            cmp(r, "max|lam| (mine, 400pt)", pred[f"r{r:.2f}.absmax_allfeeds"], mods, 2e-3)
        else:
            cmp(r, "threshold feed mM", a_ms, ms_rc, 1e-4, rel=False)
            print(f"      stated in the paper / rna_copying.py output: 10.1 mM -> analytic {a_ms:.3f}, |analytic - 10.1| = {abs(a_ms - 10.1):.3f}")
            assert abs(a_ms - 10.1) < 0.05, "threshold differs from the stated 10.1 mM beyond rounding"
        if r < 1:
            z0 = np.array([rc.T_IN / 2, rc.T_IN / 2, FEED * 1.2 / 2, FEED * 0.8 / 2])
            ee_run = rc.run(z0, FEED, k, 5000, e=0.2)[1][-1]
            cmp(r, "ee @5000 cycles", pred[f"r{r:.2f}.ee_nl_5000cyc_e0.2@20"], ee_run, 1e-9, rel=False)
            ee_long = rc.run(z0, FEED, k, 400000, e=0.2)[1][-1]
            cmp(r, "ee steady (400k cyc)", pred[f"r{r:.2f}.ee_nl_steady_e0.2@20"], ee_long, 1e-6, rel=False)
            print(f"      linear prediction {pred[f'r{r:.2f}.ee_lin_e0.2@20']:.4f}; run {ee_run:.4f} at 5000 cycles, {ee_long:.4f} steady "
                  f"(linear valid only if gain*0.2 << 1)")
    # stated: 0.67 at r = 0.40, +20% feed
    ee40 = pred["r0.40.ee_nl_5000cyc_e0.2@20"]
    print(f"\nstated plateau 0.67 at r = 0.40, +20% feed: analytic nonlinear {ee40:.4f}, linear {pred['r0.40.ee_lin_e0.2@20']:.4f}")
    lo, hi = (forced(r, FEED, 0.2)[0] for r in (0.395, 0.405))
    r67 = brentq(lambda r: forced(r, FEED, 0.2)[0] - 0.67, 0.38, 0.41, xtol=1e-6)
    print(f"    sensitivity: r = 0.395 -> {lo:.4f}, r = 0.405 -> {hi:.4f}; ee = 0.670 at r = {r67:.4f}, i.e. the stated 0.67 is "
          f"consistent with r = 0.40 as rounded (|ee(0.40) - 0.67| = {abs(ee40 - 0.67):.4f})")
    assert lo < 0.67 < hi, "stated 0.67 outside the plateau range for r = 0.40 +- 0.005"
    # window edge and onset
    k = dict(rc.HI)
    cmp(1.2, "threshold feed mM (HI)", threshold(1.2), rc.threshold(k)[0], 1e-4, rel=False)
    cmp(1.2, "m* docstring 9.9", threshold(1.2), 9.9, 6e-3)
    k = dict(rc.LO, r=1.2)
    mon = brentq(lambda m: C0 * sym(1.2, m)["S"] - W / (1 - W), 1e-9, 1e6)
    cmp(1.2, "T_in->0 onset mM", 2 * mon * (1 - (1 - W) * (1 - H)) / W, rc.onset(k), 1e-6)
    cmp(1.2, "onset docstring 8.3", 2 * mon * (1 - (1 - W) * (1 - H)) / W, 8.3, 6e-3)
    cmp(R_MAX, "r_max vs rna_copying", R_MAX, (0.5 * 0.98 / 0.02) ** (1 / 3) - 1, 1e-12)
    assert rc.onset(dict(rc.LO, r=R_MAX + 0.01)) is None and rc.onset(dict(rc.LO, r=R_MAX - 0.05)) is not None

    print(f"\nSUMMARY: {sum(rows)} of {len(rows)} comparisons within tolerance; mismatches: {fails if fails else 'none'}")
    assert not fails, "mismatch with rna_copying.py"
    print("self-check passed")


if __name__ == "__main__":
    main()
