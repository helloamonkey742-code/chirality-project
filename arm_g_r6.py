"""Arm G, R6 for 'Nonenzymatic RNA copying, 2-AI' (scorecard.py:61): feed seed (a), template seed (b), racemic-pool poisoning (c).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with matplotlib python arm_g_r6.py > arm_g_r6.log
(matplotlib is only needed because scorecard.py imports inherit.py -> ocean.py -> matplotlib.)

Scorecard R6 = "follows an inherited excess >= e_need" (scorecard.py line 9; e_need from inherit.py), here "shown to work from racemic or <= ~1% ee".

(a) FEED SEED. Letters are NOT re-derived: they are read from rna_seed_scan.log (Delta = ee_final/s >= 2, stationary linear-noise s,
    Omega = N_setting/m_in per mM; ALL LETTERS line printed with its line number). This script re-derives the gains/|lam| it needs and
    adds a supplement (block S): feed = the setting's own concentration, with T_in fixed or scaled with the feed.
(b) TEMPLATE SEED. Racemic feed, excess ee0 only in T_D vs T_L at the racemic fixed point (the primer-template duplex is the model's unit;
    there is no free primer pool). Linear response of the difference mode: [x, y](t) = Jd^t [x0, 0], x0 = sqrt(2) T* ee0 (eigen-expansion,
    checked against matrix_power and the nonlinear map). rho(t) = ee_T(t)/ee0. LOST/KEPT/GROWN at tau; e_min^T = 2 s/|g_T(tau)|;
    letters weak (Delta >= 2, as (a)) and strong (weak AND not LOST). Cycle length is unmeasured: 1 h and 24 h bracket n_cyc = tau/L.
(c) POISONING. Per-attempt S (racemic pool vs pure pool at the same total monomer, vs same right-hand level), onset feeds, steady
    reactor (feed excess e = 0 vs e = 1), yield Y and template amplification Q, ee-vs-yield curve, survival of a feed excess.

PRE-REGISTERED PREDICTIONS (arm_g_r6_predictions.txt, saved before the first results run; informed by rna_seed_scan.log, arm_g_rates.log
and closed-form hand estimates, not blind to them; blind to this script's output)
  sha256 79bea991b1ec5e91a3ab008b4fee5c85bece239e2e153392edf458286d8c64ce   2026-10-21 00:35:05 PDT (2026-10-21T07:35:05Z)
  (this script asserts the file still has that hash)
Smoke runs before the registered run: none were made of this script's quantities. The only command run earlier was `python scorecard.py`
(import error: no matplotlib). Anything added after the first results run is labelled 'NOT pre-registered' in the output.
"""
import hashlib
import time
import warnings

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

import rna_copying as rc
from ocean import NA, YR
from rna_seed_scan import linear, readout, template_seed
from scorecard import SETTINGS, spec

PRED_SHA = "79bea991b1ec5e91a3ab008b4fee5c85bece239e2e153392edf458286d8c64ce"
RS, FEEDS = (0.40, 0.56, 0.98), (5.0, 20.0, 200.0)
LEN_H = {"1 h": 1.0, "24 h": 24.0}                       # cycle length brackets (unmeasured)
HPY = YR / 3600.0                                        # hours per year as ocean.YR
LOG = "rna_seed_scan.log"
SQ2 = np.sqrt(2.0)


# ---------------------------------------------------------------- (b) linear template-seed response
def slopes(z, m_in, k):
    """Readout slope wrt the difference coordinates [x, y] (same finite differences as rna_seed_scan.linear)."""
    out = np.empty(2)
    for i in range(2):
        d = 1e-6 * z[2 * i]
        out[i] = (readout(z + d * rc.PM[:, i], m_in, k) - readout(z - d * rc.PM[:, i], m_in, k)) / (2 * d)
    return out


def seed_block(m_in, k, noise=True):
    """Eigen-expansion of the template-seed response. rho(t) = x(t)/x0 = sum a_rho vals^t; g_T(t) = readout(t)/ee0 = sum a_g vals^t."""
    z = rc.fixed_point(m_in, **k)
    _, Jd = rc.mode(z, m_in, k)
    vals, V = np.linalg.eig(Jd)
    Vi = np.linalg.inv(V)
    a_rho = V[0, :] * Vi[:, 0]
    a_g = SQ2 * z[0] * (slopes(z, m_in, k) @ V) * Vi[:, 0]
    assert abs(a_rho.sum() - 1) < 1e-9, "rho(0) must be 1"
    su2 = linear(m_in, k)[5] if noise else np.nan
    return dict(z=z, Jd=Jd, vals=vals, a_rho=a_rho, a_g=a_g, su2=su2, lam=abs(vals).max(), cplx=vals.dtype.kind == "c")


def exact(a, vals, t):
    """Re sum a_k vals_k^t (complex logs; fine for huge t, underflows to 0)."""
    with np.errstate(all="ignore"):
        return float((a * np.exp(t * np.log(vals.astype(complex)))).sum().real)


def logenv(a, vals, t):
    """log10 of sum |a_k| |vals_k|^t: an upper bound on |sum a_k vals_k^t| (so 'envelope < 0.5' proves LOST)."""
    with np.errstate(all="ignore"):
        l = np.log10(np.abs(a)) + t * np.log10(np.abs(vals))
    m = l.max()
    return m + np.log10(np.sum(10 ** (l - m)))


def rho_max(Jd, steps=100000):
    """max_{t >= 1} x(t)/x0 for [x, y](0) = [1, 0] (transient growth check); stops once the state has decayed below 1e-30."""
    v, m = np.array([1.0, 0.0]), -np.inf
    for _ in range(steps):
        v = Jd @ v
        m = max(m, v[0])
        if np.abs(v).max() < 1e-30:
            break
    return m


def classify(b, ncyc, gmax):
    if b["lam"] >= 1 or gmax > 1 + 1e-9:
        return "GROWN"
    if logenv(b["a_rho"], b["vals"], ncyc) < np.log10(0.5):
        return "LOST"
    return "KEPT" if exact(b["a_rho"], b["vals"], ncyc) >= 0.5 else "LOST"


def letter(flags):
    return "Y" if all(flags) else "P" if any(flags) else "N"


def be_cycles(b):
    """cycles for the envelope of rho to fall to 0.5: leading term only."""
    i = np.argmax(abs(b["vals"]))
    amp = abs(b["a_rho"][i]) * (2 if b["cplx"] else 1)
    return np.log(amp / 0.5) / -np.log(abs(b["vals"][i]))


# ---------------------------------------------------------------- (c) pure vs racemic
def s_pure(m_total, n=3):
    """P(full length) for a pure pool: p = 1, mu = A m_total/M_REF."""
    return float(poisson.sf(n - 1, rc.A * m_total / rc.M_REF))


def onset_pure(k):
    """Feed at which one-hand copying just keeps up with washout: c S(m) = w/(1-w), m = w m_in/(1-(1-w)(1-h))."""
    m = brentq(lambda m: k["c"] * s_pure(m, k["n"]) - k["w"] / (1 - k["w"]), 1e-9, 1e6)
    return m * (1 - (1 - k["w"]) * (1 - rc.H)) / k["w"]


def reactor(m_in, k, e, block=20000, maxcyc=400000):
    """Steady state of rna_copying.step with feed excess e by direct iteration from the feed-filled reactor."""
    z = np.array([rc.T_IN / 2, rc.T_IN / 2, m_in * (1 + e) / 2, m_in * (1 - e) / 2])
    cyc = 0
    while True:
        for _ in range(block):
            z = rc.step(z, m_in, e=e, **k)[0]
        cyc += block
        z2, XD, XL = rc.step(z, m_in, e=e, **k)
        conv = bool(np.all(np.abs(z2 - z) <= 1e-9 * np.maximum(np.abs(z), 1e-6)))
        if conv or cyc >= maxcyc:
            return dict(z=z2, XD=XD, XL=XL, conv=conv, cyc=cyc, Y=XD + XL, ee=(XD - XL) / (XD + XL) if XD + XL > 0 else 0.0)


def set_tin(tin):
    old, rc.T_IN = rc.T_IN, tin
    return old


def fmt_e(x):
    return "none(>1)" if x > 0 else f"1e{x:+.1f}"


def main():
    t0 = time.time()
    assert hashlib.sha256(open("arm_g_r6_predictions.txt", "rb").read()).hexdigest() == PRED_SHA, "predictions file changed"
    names = list(SETTINGS)
    need = {n: e for n, _, e in spec()}
    N = {n: c * V * 1000 * NA for n, (c, tau, V) in SETTINGS.items()}
    tau_h = {n: SETTINGS[n][1] / 3600.0 for n in names}
    print("arm_g_r6: R6 for nonenzymatic RNA copying (rna_copying.py), feed seed vs template seed, racemic-pool poisoning")
    print("settings: " + "; ".join(f"{n}: c = {SETTINGS[n][0]:g} M, tau = {SETTINGS[n][1] / YR:.0e} yr = {tau_h[n]:.2e} h, N = {N[n]:.2e}, e_need = {need[n]:.1e}" for n in names))
    print(f"model: rna_copying with n = 3, c = 0.5, w = 0.02, T_in = {rc.T_IN * 1e3:g} uM, kobs = {rc.KOBS}/h, h = {rc.H:.4e}; cycle length L = {', '.join(LEN_H)} (unmeasured)")

    # ---- self-checks (hand-computed) ----
    mu = rc.A * 15.0 / rc.M_REF                                           # racemic, m = 10 mM each, r = 0.5: tot = 15 mM
    S_hand = (2 / 3) ** 3 * (1 - np.exp(-mu) * (1 + mu + mu ** 2 / 2))     # P(N>=3) = 1 - e^-mu (1 + mu + mu^2/2)
    S_code = rc.attempt(10.0, 10.0, 3, 0.5)[0]
    assert abs(S_code - S_hand) < 1e-12 and abs(S_hand - 0.24791) < 2e-5, (S_code, S_hand)
    assert abs(s_pure(10.0) - poisson.sf(2, rc.A * 10 / rc.M_REF)) < 1e-15
    k0 = dict(rc.LO, r=0.40)
    b0 = seed_block(20.0, k0)
    tr, det = np.trace(b0["Jd"]), np.linalg.det(b0["Jd"])
    lam_q = max(abs((tr + np.sqrt(complex(tr * tr - 4 * det))) / 2), abs((tr - np.sqrt(complex(tr * tr - 4 * det))) / 2))
    assert abs(lam_q - b0["lam"]) < 1e-12 and abs(b0["lam"] - 0.97943) < 5e-5, "lam does not match rna_seed_scan.log Table 1"
    for t in (1, 7, 100):
        assert abs(exact(b0["a_rho"], b0["vals"], t) - np.linalg.matrix_power(b0["Jd"], t)[0, 0]) < 1e-9, "eigen-expansion vs matrix_power"
    z0 = np.array([rc.T_IN / 2, rc.T_IN / 2, 20.0 * 1.2 / 2, 20.0 * 0.8 / 2])
    plateau = rc.run(z0, 20.0, dict(rc.LO, r=0.396), 5000, e=0.2)[1][-1]
    assert abs(plateau - 0.668) < 1e-3, plateau
    lin_err = []                                                           # linear template response vs nonlinear map, ee0 = 1e-4, t = 100
    for r in RS:
        k = dict(rc.LO, r=r)
        b = seed_block(20.0, k, noise=False)
        z = b["z"] * [1 + 1e-4, 1 - 1e-4, 1, 1]
        for _ in range(100):
            z = rc.step(z, 20.0, **k)[0]
        nl = (z[0] - z[1]) / (z[0] + z[1]) / 1e-4
        lin_err.append(nl / exact(b["a_rho"], b["vals"], 100) - 1)
        assert abs(lin_err[-1]) < 1e-3, ("linear vs nonlinear template response", r, nl, exact(b["a_rho"], b["vals"], 100))
    print(f"\nself-checks passed: S(10,10; r = .5) = {S_code:.6f} (hand {S_hand:.6f}); lam r = .4 by quadratic formula = {lam_q:.6f} (log 0.97943); "
          f"eigen-expansion = matrix_power at t = 1, 7, 100; plateau r = .396, +20% feed, 5000 cycles = {plateau:.4f} (arm_g_rates.log:40: 0.668); "
          f"linear vs nonlinear template response at t = 100, ee0 = 1e-4: " + ", ".join(f"{x:+.1e}" for x in lin_err))

    # ---- (a) provenance ----
    print("\n=== (a) FEED SEED: letters from rna_seed_scan.log (not re-derived here) ===")
    lines = open(LOG).read().splitlines()
    t2 = False
    for i, ln in enumerate(lines, 1):
        t2 = True if ln.startswith("TABLE 2") else False if ln.startswith("TABLE 3") else t2
        if ln.startswith("ALL LETTERS") or "letters at 20 mM" in ln or (t2 and ln[:5].strip() in ("0.40", "0.56", "0.98")):
            print(f"  {LOG}:{i}: {ln.strip()}")
    print("  definition used there: followed = Delta = ee_final/s >= 2; e_min = 2 s/G; Y if e_min <= e_need in all three settings.")
    print("  soundness: a signal-to-noise test of a bounded response (gain G = 2.4-115, racemic state stable) to a continuously supplied excess,")
    print("  against an e_need that comes from a bifurcation ramp (inherit.py). It says the readout tracks the feed excess above counting noise.")
    print("  It does not say the system selects one hand; the product excess stays below 1 for small seeds (see (c) ee-vs-yield).")

    # ---- (b) template seed ----
    print("\n=== (b) TEMPLATE SEED (racemic feed; excess ee0 in T_D vs T_L only) ===")
    blk, gmax = {}, {}
    for r in RS:
        for f in FEEDS:
            k = dict(rc.LO, r=r)
            blk[r, f] = seed_block(f, k)
            gmax[r, f] = rho_max(blk[r, f]["Jd"])
    print("difference-mode block at the racemic state: |lam| (leading), complex pair?, a_rho (weight of the leading mode in rho), max rho(t >= 1), T* (uM)")
    print(f"{'r':>5} {'feed':>5} {'|lam|':>9} {'cplx':>5} {'|a_rho|':>8} {'max rho(t>=1)':>14} {'T* uM':>8} {'log10 s_pond':>13}")
    for (r, f), b in blk.items():
        i = np.argmax(abs(b["vals"]))
        s = np.sqrt(b["su2"] / (N[names[0]] / f))
        print(f"{r:5.2f} {f:5g} {b['lam']:9.5f} {str(b['cplx']):>5} {abs(b['a_rho'][i]):8.4f} {gmax[r, f]:14.6f} {b['z'][0] * 1e3:8.2f} {np.log10(s):13.2f}")
    print("\nOUTCOME at tau (LOST: rho(tau) < 0.5; KEPT >= 0.5; GROWN: rho > 1 at some t) and e_min^T = 2 s/|g_T(tau)|; 20 mM feed")
    print(f"{'r':>5} {'L':>5} " + "".join(f"| {n.split()[-1]:>6} log10 ncyc  log10 rho_env  outcome  log10 e_min^T  log10 e_need " for n in names) + "| weak strong")
    cell, res20 = {}, {}
    for r in RS:
        b = blk[r, 20.0]
        for Ln, Lh in LEN_H.items():
            row, wk, st = "", [], []
            for n in names:
                nc = tau_h[n] / Lh
                out = classify(b, nc, gmax[r, 20.0])
                s = np.sqrt(b["su2"] / (N[n] / 20.0))
                lge = np.log10(2 * s) - logenv(b["a_g"], b["vals"], nc)
                lgr = logenv(b["a_rho"], b["vals"], nc)
                ok = lge <= np.log10(need[n])
                wk.append(ok)
                st.append(ok and out != "LOST")
                cell[r, 20.0, Ln, n] = (out, lge, lgr)
                row += f"| {'':>6} {np.log10(nc):9.2f}  {lgr:13.3g}  {out:>7}  {fmt_e(lge):>13}  {np.log10(need[n]):12.1f} "
            res20[r, Ln] = (letter(wk), letter(st))
            print(f"{r:5.2f} {Ln:>5} {row}|  {res20[r, Ln][0]:>3}   {res20[r, Ln][1]:>3}")
    # all 54 cells
    tally = {"LOST": 0, "KEPT": 0, "GROWN": 0}
    for r in RS:
        for f in FEEDS:
            for Ln, Lh in LEN_H.items():
                for n in names:
                    tally[classify(blk[r, f], tau_h[n] / Lh, gmax[r, f])] += 1
    print(f"\nall 54 cells (3 r x 3 feeds x 3 settings x 2 cycle lengths): {tally}")
    fin = [(r, Ln, n) for r in RS for Ln in LEN_H for n in names if cell[r, 20.0, Ln, n][1] <= 0]
    print(f"20 mM cells with e_min^T < 1: {fin if fin else 'none'}")
    for (r, Ln, n) in fin:
        print(f"    r = {r}, L = {Ln}, {n}: e_min^T = 10^{cell[r, 20.0, Ln, n][1]:.2f} = {10 ** cell[r, 20.0, Ln, n][1]:.1e}; e_need = {need[n]:.1e}")

    print("\nBREAK-EVEN cycle length for a template seed to be KEPT (rho_env = 0.5) over tau, if the per-cycle parameters (w = 0.02 ...) were unchanged")
    print("  (a scale, not a model result: w, h would change with the cycle length)")
    print(f"{'r':>5} {'feed':>5} {'cycles to 0.5':>14} " + " ".join(f"{n.split()[-1]:>14}" for n in names))
    be = {}
    for (r, f), b in blk.items():
        nb = be_cycles(b)
        txt = []
        for n in names:
            Lbe = tau_h[n] / nb                                            # hours
            be[r, f, n] = Lbe
            txt.append(f"{Lbe / 24 / 365.25:11.3g} yr" if Lbe > 24 * 365.25 * 0.5 else f"{Lbe / 24:11.3g} d ")
        print(f"{r:5.2f} {f:5g} {nb:14.1f} " + " ".join(f"{t:>14}" for t in txt))

    print("\nNONLINEAR CHECK (20 mM, 20000 direct iterations of rna_copying.step; |ee_T| = max over the last 200 cycles)")
    nl = {}
    print(f"{'r':>5} " + " ".join(f"ee0 = {e0:<5}" for e0 in (0.1, 0.5, 0.99)) + "   |ee_T|/ee0")
    for r in RS:
        k = dict(rc.LO, r=r)
        nl[r] = [template_seed(20.0, k, e0) for e0 in (0.1, 0.5, 0.99)]
        print(f"{r:5.2f} " + " ".join(f"{v:11.2e}" for v in nl[r]) + "   " + " ".join(f"{v / e0:9.1e}" for v, e0 in zip(nl[r], (0.1, 0.5, 0.99))))

    print("\nPOSITIVE CONTROL (can the classifier say GROWN?): r = 1.141 (arm_g_rates.log:44: Fig. 9 reading A, lam 1.00221, m* 10.08 mM), 20 mM")
    kc = dict(rc.LO, r=1.141)
    bc = seed_block(20.0, kc, noise=False)
    zc = bc["z"] * [1 + 1e-3, 1 - 1e-3, 1, 1]
    grown_at, ee_t = None, 0.0
    for i in range(200000):
        zc = rc.step(zc, 20.0, **kc)[0]
        ee_t = abs(zc[0] - zc[1]) / (zc[0] + zc[1])
        if ee_t > 0.5:
            grown_at = i + 1
            break
    lab7 = classify(bc, 1e5, rho_max(bc["Jd"], 100000))
    print(f"  leading |lam| = {bc['lam']:.5f}; seed ee0 = 1e-3 reaches |ee_T| > 0.5 after {grown_at} cycles (|ee_T| = {ee_t:.3f}); classifier: {lab7}")

    # ---- (c) poisoning ----
    print("\n=== (c) RACEMIC-POOL POISONING INSIDE THE MODEL ===")
    print("C1. per-attempt full-length probability S: racemic pool (m each hand) vs pure pool with the same TOTAL monomer (2m) vs same right-hand level (m)")
    print("    'fresh' = m = m_in/2 at 20 mM feed; 'reactor' = m = steady racemic monomer level at 20 mM feed. r = 0: wrong hand is an inert spectator.")
    print(f"{'r':>5} {'pool':>8} {'m (mM)':>7} {'mu rac':>7} {'p':>6} {'S rac':>8} {'S pure(2m)':>11} {'rac/pure':>9} {'S(m,0)':>8} {'rac/(m,0)':>10}")
    c1 = {}
    for r in (0.0, 0.06, 0.40, 0.56, 0.98):
        k = dict(rc.LO, r=r)
        for pool, m in (("fresh", 10.0), ("reactor", rc.fixed_point(20.0, **k)[2])):
            S = rc.attempt(m, m, 3, r)[0]
            Sp, S1 = s_pure(2 * m), rc.attempt(m, 0.0, 3, r)[0]
            c1[r, pool] = (m, S, Sp, S1)
            print(f"{r:5.2f} {pool:>8} {m:7.3f} {rc.A * m * (1 + r) / rc.M_REF:7.3f} {1 / (1 + r):6.3f} {S:8.4f} {Sp:11.4f} {S / Sp:9.3f} {S1:8.4f} {S / S1:10.3f}")
    print("    in the reactor the racemic S is pinned by washout (c S = w/(1-w) = 0.0408 wherever T* >> T_in); poisoning shows up as the monomer level")
    print("    and feed needed to reach it, not as a lower S.")

    print("\nC2. onset feed (copying just keeps up with washout, c S = w/(1-w)): racemic pool vs pure pool, same total feed")
    print(f"{'r':>5} {'racemic mM':>11} {'pure mM':>9} {'ratio':>6}")
    po = onset_pure(rc.LO)
    assert abs(po - rc.onset(dict(rc.LO, r=0.0)) / 2) < 1e-6 * po, "pure onset formula vs rna_copying.onset(r = 0)/2"
    on = {}
    for r in (0.06,) + RS:
        on[r] = rc.onset(dict(rc.LO, r=r))
        print(f"{r:5.2f} {on[r]:11.3f} {po:9.3f} {on[r] / po:6.2f}")
    print(f"    pure onset = rna_copying.onset(r = 0)/2 = {po:.3f} mM (asserted equal to the one-hand formula)")

    print("\nC3. steady reactor, same total feed: racemic (e = 0) vs pure (e = 1). Q = T_D*/(T_in/2); Y = full-length copies per cycle, both hands (uM/cycle)")
    print(f"{'r':>5} {'feed':>5} | {'T* rac uM':>9} {'Q rac':>7} {'Y rac':>9} | {'T_D* pure':>9} {'Q pure':>7} {'Y pure':>9} {'max|m_L|':>9} | {'Y rac/Y pure':>12} conv")
    c3 = {}
    for r in RS:
        k = dict(rc.LO, r=r)
        for f in (5.0, 20.0):
            a, p = reactor(f, k, 0.0), reactor(f, k, 1.0)
            c3[r, f] = (a, p)
            print(f"{r:5.2f} {f:5g} | {a['z'][0] * 1e3:9.2f} {a['z'][0] / (rc.T_IN / 2):7.2f} {a['Y'] * 1e3:9.3f} | {p['z'][0] * 1e3:9.2f} "
                  f"{p['z'][0] / (rc.T_IN / 2):7.2f} {p['Y'] * 1e3:9.3f} {abs(p['z'][3]):9.1e} | {a['Y'] / p['Y']:12.3f} {a['conv'] and p['conv']}")

    print("\nC4. ee-vs-yield at 20 mM: feed excess e from racemic to pure; product excess ee, yield relative to the pure pool, gain = ee/e")
    print(f"{'r':>5} {'e':>6} {'ee':>9} {'gain':>8} {'Y/Y(e=1)':>9} {'T_D* uM':>9} {'T_L* uM':>9} conv")
    c4 = {}
    for r in RS:
        k = dict(rc.LO, r=r)
        for e in (0.0, 0.01, 0.1, 0.2, 0.5, 0.9, 1.0):
            x = c4[r, e] = reactor(20.0, k, e) if e not in (0.0, 1.0) else c3[r, 20.0][0 if e == 0 else 1]
            print(f"{r:5.2f} {e:6.2f} {x['ee']:9.4f} {(x['ee'] / e if e else float('nan')):8.3f} {x['Y'] / c3[r, 20.0][1]['Y']:9.3f} "
                  f"{x['z'][0] * 1e3:9.2f} {x['z'][1] * 1e3:9.2f} {x['conv']}")
    g = {r: c4[r, 0.01]['ee'] / 0.01 for r in RS}
    print("    gain at e = 0.01: " + ", ".join(f"r = {r}: {g[r]:.3g}" for r in RS) + "  (rna_seed_scan.log Table 1: 3.4, 5.25, 79)")

    print("\nC5. poison-off control: gain (linear, from rna_seed_scan.linear) with the wrong hand inert (r = 0) vs measured r, 20 mM")
    for r in (0.0, 0.06) + RS:
        zz, lam, G = linear(20.0, dict(rc.LO, r=r))[:3]
        print(f"    r = {r:4.2f}: gain {G:7.3f}, |lam| {lam:.5f}, T* {zz[0] * 1e3:7.2f} uM")

    # ---- block S: setting-consistent supplement ----
    print("\n=== S. SETTING-CONSISTENT SUPPLEMENT (feed = the setting's own concentration; Omega = N_setting/m_in; NOT part of the earlier definition) ===")
    print("    rna_seed_scan.py used the setting's N with a 5-200 mM feed. Here m_in = pond 10 mM, lake 0.1 mM, ocean 1e-3 mM.")
    print("    T_in (assumed template feed, 3 uM at the experiment's 24 mM) either fixed, or scaled T_in = 3 uM x m_in/20 mM.")
    print(f"{'r':>5} {'T_in':>7} " + "".join(f"| {n.split()[-1]:>6} m_in(mM)  Q-1    gain  e_min/e_need " for n in names) + "| letter")
    sup = {}
    for conv in ("fixed", "scaled"):
        for r in RS:
            row, rat, extra = "", [], []
            for n in names:
                m_in = SETTINGS[n][0] * 1000
                old = set_tin(3e-3 if conv == "fixed" else 3e-3 * m_in / 20.0)
                try:
                    z, lam, G, _, _, su2 = linear(m_in, dict(rc.LO, r=r))
                    Q = z[0] / (rc.T_IN / 2)
                finally:
                    rc.T_IN = old
                s = np.sqrt(su2 / (N[n] / m_in))
                ratio = 2 * s / G / need[n]
                sup[conv, r, n] = dict(Q=Q, G=G, ratio=ratio, lam=lam)
                rat.append(ratio <= 1)
                row += f"| {'':>6} {m_in:9.3g} {Q - 1:7.1e} {G:7.3f} {ratio:12.3g} "
            sup[conv, r] = letter(rat)
            print(f"{r:5.2f} {conv:>7} {row}|   {sup[conv, r]}")

    # ---- P13 sensitivity (template seed) ----
    print("\n=== SENSITIVITY of the template-seed outcome at (pond, 24 h) = 36,525 cycles, 20 mM, one parameter at a time ===")
    print(f"{'r':>5} {'change':>12} {'|lam|':>9} {'cplx':>5} {'outcome':>8} {'log10 rho_env':>14} {'T* uM':>8}")
    sens_bad, sens_skip = [], []
    nc_short = tau_h[names[0]] / 24.0
    for r in RS:
        for lab, ch in (("base", {}), ("w=0.005", dict(w=0.005)), ("w=0.05", dict(w=0.05)), ("c=0.1", dict(c=0.1)), ("c=1.0", dict(c=1.0)),
                        ("T_in=0.3uM", dict(tin=3e-4)), ("T_in=30uM", dict(tin=3e-2))):
            tin = ch.pop("tin", None)
            k = dict(rc.LO, r=r, **ch)
            old = set_tin(tin) if tin else rc.T_IN
            try:
                b = seed_block(20.0, k, noise=False)
                gm = rho_max(b["Jd"])
            except Exception:
                sens_skip.append((r, lab))
                print(f"{r:5.2f} {lab:>12}   no fixed point / no racemic state")
                continue
            finally:
                rc.T_IN = old
            if b["lam"] >= 1:
                sens_skip.append((r, lab))
                print(f"{r:5.2f} {lab:>12} {b['lam']:9.5f} {str(b['cplx']):>5}  racemic state unstable (not scored)")
                continue
            out = classify(b, nc_short, gm)
            if out != "LOST":
                sens_bad.append((r, lab, out))
            print(f"{r:5.2f} {lab:>12} {b['lam']:9.5f} {str(b['cplx']):>5} {out:>8} {logenv(b['a_rho'], b['vals'], nc_short):14.3g} {b['z'][0] * 1e3:8.2f}")

    # ---- predictions vs result ----
    print("\n=== PREDICTIONS vs RESULT (arm_g_r6_predictions.txt) ===")
    ck = lambda ok: "held" if ok else "FAILED"
    P = {}
    ref = {(0.40, 5.0): 0.99059, (0.40, 20.0): 0.97943, (0.40, 200.0): 0.97845, (0.56, 5.0): 0.99003, (0.56, 20.0): 0.98797,
           (0.56, 200.0): 0.98551, (0.98, 5.0): 0.98860, (0.98, 20.0): 0.99948, (0.98, 200.0): 0.99941}
    P[1] = all(abs(blk[key].get("lam") - v) < 5e-5 for key, v in ref.items())
    print(f"  P1 |lam| matches rna_seed_scan.log Table 1 to 5e-5 in 9 cells: {ck(P[1])} (max |diff| {max(abs(blk[key]['lam'] - v) for key, v in ref.items()):.1e})")
    P[2] = all(v <= 1 + 1e-9 for v in gmax.values())
    print(f"  P2 no transient growth of a template seed (max rho(t >= 1) <= 1): {ck(P[2])} (max {max(gmax.values()):.6f})")
    P[3] = tally == {"LOST": 54, "KEPT": 0, "GROWN": 0}
    print(f"  P3 54/54 LOST: {ck(P[3])} ({tally})")
    rng = {0.40: (2.2 * 365.25, 3.9 * 365.25), 0.56: (1.3 * 365.25, 2.3 * 365.25), 0.98: (20.0, 35.0)}
    p4 = {r: be[r, 20.0, names[0]] / 24 for r in RS}
    P[4] = all(rng[r][0] <= p4[r] <= rng[r][1] for r in RS)
    print(f"  P4 break-even cycle length, pond, 20 mM: " + "; ".join(f"r = {r}: {p4[r]:.1f} d = {p4[r] / 365.25:.2f} yr (pred {rng[r][0]:.0f}-{rng[r][1]:.0f} d)" for r in RS) + f": {ck(P[4])}")
    one = [(r, Ln, n) for r in RS for Ln in LEN_H for n in names if cell[r, 20.0, Ln, n][1] <= 0]
    p5a = one == [(0.98, "24 h", names[0])] and 1e-5 <= 10 ** cell[0.98, 20.0, "24 h", names[0]][1] <= 1
    p5b = all(cell[r, 20.0, Ln, n][1] > np.log10(need[n]) for r in RS for Ln in LEN_H for n in names) and all(res20[r, Ln][0] == "N" for r in RS for Ln in LEN_H)
    p5c = all(res20[r, Ln][1] == "N" for r in RS for Ln in LEN_H)
    P[5] = p5a and p5b and p5c
    print(f"  P5 a (only pond/24 h/r = 0.98 has e_min^T < 1): {ck(p5a)}; b (weak letter N in 6/6, e_min^T > e_need everywhere): {ck(p5b)}; c (strong letter N in 6/6): {ck(p5c)}")
    P[6] = all(v / e0 < 1e-3 for r in RS for v, e0 in zip(nl[r], (0.1, 0.5, 0.99)))
    print(f"  P6 nonlinear decay |ee_T(20000)|/ee0 < 1e-3 in 9/9: {ck(P[6])}")
    P[7] = bc["lam"] > 1 and grown_at is not None and lab7 == "GROWN"
    print(f"  P7 positive control r = 1.141: lam > 1, seed reaches 0.5, classifier GROWN: {ck(P[7])}")
    r8 = {0.40: (5.5, 7.5), 0.56: (5.8, 7.9), 0.98: (6.5, 9.0)}
    ratios = [on[r] / po for r in RS]
    P[8] = (all(r8[r][0] <= on[r] <= r8[r][1] for r in RS) and 2.5 <= po <= 3.3 and all(1.8 <= x <= 3.2 for x in ratios)
            and ratios[0] < ratios[1] < ratios[2])
    print(f"  P8 onset feeds (mM): racemic " + ", ".join(f"{on[r]:.2f}" for r in RS) + f" (pred 6.3, 6.6, 7.6); pure {po:.3f} (pred 2.86); ratios " + ", ".join(f"{x:.2f}" for x in ratios) + f": {ck(P[8])}")
    qr = {r: c3[r, 5.0][0]["z"][0] / (rc.T_IN / 2) for r in RS}
    qp = {r: c3[r, 5.0][1]["z"][0] / (rc.T_IN / 2) for r in RS}
    P[9] = all(q < 3 for q in qr.values()) and all(30 <= q <= 250 for q in qp.values()) and max(qp.values()) / min(qp.values()) - 1 < 0.01
    print(f"  P9 5 mM: racemic Q " + ", ".join(f"{qr[r]:.2f}" for r in RS) + " (< 3); pure Q " + ", ".join(f"{qp[r]:.1f}" for r in RS) + f" (30-250, r-independent to 1%): {ck(P[9])}")
    yr = {r: c3[r, 20.0][0]["Y"] / c3[r, 20.0][1]["Y"] for r in RS}
    P[10] = all(0.3 <= v <= 0.9 for v in yr.values()) and yr[0.40] > yr[0.56] > yr[0.98]
    print(f"  P10 Y_rac/Y_pure at 20 mM: " + ", ".join(f"{yr[r]:.3f}" for r in RS) + f" (pred 0.61, 0.55, 0.44; range 0.3-0.9, decreasing): {ck(P[10])}")
    lg = {0.40: 3.4, 0.56: 5.25, 0.98: 79.0}
    gains = {(r, f): (blk[r, f] and linear(f, dict(rc.LO, r=r))[2]) for r in RS for f in FEEDS}
    P[11] = (all(abs(g[r] / lg[r] - 1) < 0.05 for r in RS) and all(v >= 2 for v in gains.values())
             and all(abs(c4[r, 0.01]["Y"] / c4[r, 0.0]["Y"] - 1) < 0.02 for r in RS))
    print(f"  P11 gain at e = 0.01: " + ", ".join(f"{g[r]:.3g}" for r in RS) + " (log 3.4, 5.25, 79, +-5%); linear gain >= 2 in 9 cells (min "
          f"{min(gains.values()):.3f}); Y(0.01)/Y(0) - 1 = " + ", ".join(f"{c4[r, 0.01]['Y'] / c4[r, 0.0]['Y'] - 1:+.1e}" for r in RS) + f": {ck(P[11])}")
    a_ = all(sup[c_, r, n]["Q"] - 1 < 1e-3 and 2.9 <= sup[c_, r, n]["G"] <= 3.1 for c_ in ("fixed", "scaled") for r in RS for n in names[1:]) \
        and all(sup[c_, r, names[0]]["Q"] > 3 for c_ in ("fixed", "scaled") for r in RS)
    b_ = all(0.05 <= sup["fixed", r, names[1]]["ratio"] <= 0.5 and 0.001 <= sup["fixed", r, names[2]]["ratio"] <= 0.02 for r in RS) and all(sup["fixed", r] == "Y" for r in RS)
    c_ = all(1 <= sup["scaled", r, names[1]]["ratio"] <= 5 and 0.3 <= sup["scaled", r, names[2]]["ratio"] <= 1.0 for r in RS) and all(sup["scaled", r] != "Y" for r in RS)
    P[12] = a_ and b_ and c_
    print(f"  P12 a (no copying at lake/ocean, gain ~ 3, pond Q > 3): {ck(a_)}; b (T_in fixed: ratios in range, Y): {ck(b_)}; c (T_in scaled: lake ratio 1-5, ocean 0.3-1, not Y): {ck(c_)}")
    P[13] = not sens_bad
    print(f"  P13 template seed LOST at (pond, 24 h) in every scored sensitivity cell: {ck(P[13])} (not LOST: {sens_bad if sens_bad else 'none'}; not scored: {sens_skip if sens_skip else 'none'})")
    # verdicts
    V1 = all(x >= 1.8 for x in ratios)
    V2 = all(c3[r, 20.0][0]["z"][0] / (rc.T_IN / 2) >= 50 for r in RS)
    V3 = all(v >= 2 for v in gains.values())
    print("\nVERDICT (c): V1 racemic onset >= 1.8 x pure onset: " + str(V1) + "; V2 racemic Q >= 50 at 20 mM: " + str(V2) + " (Q = "
          + ", ".join(f"{c3[r, 20.0][0]['z'][0] / (rc.T_IN / 2):.0f}" for r in RS) + "); V3 feed-excess gain >= 2 in all 9 cells: " + str(V3))
    print("LETTERS (b) at 20 mM, weak / strong: " + "; ".join(f"r = {r}, L = {Ln}: {res20[r, Ln][0]}/{res20[r, Ln][1]}" for r in RS for Ln in LEN_H))
    print("LETTERS (a) supplement (block S): " + "; ".join(f"{c_} T_in, r = {r}: {sup[c_, r]}" for c_ in ("fixed", "scaled") for r in RS))
    print(f"\nself-check passed; runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    warnings.simplefilter("default")
    main()
