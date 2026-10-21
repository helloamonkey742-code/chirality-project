"""Arm G (rna_copying.py): does the fed RNA-copying reactor follow an inherited excess as small as each setting's e_need?

Scorecard R6 = "follows an inherited excess >= e_need" (scorecard.py line 9). inherit.py: "follows" = Delta >= 2, noise eps = 1/N.

METHOD
  seed ee0   enantiomeric excess of the monomer FEED, m_in (1 +/- ee0)/2: the inherited excess (rna_copying.run(e = ee0)).
  ee_final   steady-state product excess (D-L)/(D+L) of full-length copies, the rna_copying.run readout. gain = ee_final/ee0.
             ee0 >= 1e-5: nonlinear fixed point of rna_copying.step (continuation from the linear solution, attractor
             checked). ee0 < 1e-5: exact linear response rho (I - J)^-1 b, because 1 + 1e-18 == 1 in float64. The two agree
             at 1e-5 (asserted) and the gain is constant below it.
  noise      rna_copying.py HAS a noise term (ramp(): chemical-Langevin counting noise on the two template concentrations,
             variance (X + w (T + X) + w T_in/2)/Omega per cycle; monomers are mM pools). Used unchanged. inherit.py's
             convention eps = 1/N is the same thing with Omega = molecules per concentration unit, so
             Omega = N_setting / m_in per mM, N_setting = c V 1000 N_A (scorecard.SETTINGS, as e_need). Stationary readout
             variance s^2 from the discrete Lyapunov equation of the linearised map (stable r < 1: no bifurcation to ramp
             through, so inherit's ramp factor (pi k c tau)^(1/4) does not enter). Delta = ee_final / s, P(win) = Phi(Delta)
             for a snapshot of the readout; followed = Delta >= 2. The linear-noise s is checked against a Langevin
             ensemble (ramp()'s noise) at a small Omega; s scales exactly as Omega^-1/2.
  letter     Y: e_min <= e_need in all three settings; P: in some; N: in none (scorecard.py line 9).
  also       seeding the TEMPLATES (rna_copying.seed_run semantics, no feed excess) for contrast.

PRE-REGISTERED PREDICTIONS (rna_seed_scan_predictions.txt, saved before the first run)
  sha256 14c5f4ea46961666d2b488f4dd90d634121afc82527fd027a9ab82e02fabb801   2026-10-19 23:37:34 PDT (2026-10-20T06:37:34Z)
  (this script asserts the file still has that hash)
  NOT pre-registered, added after the first run: TABLE 4 (sensitivity to the assumed self-copying efficiency c = 0.1-1, which
  is a scanned model parameter, not a measured one, and to adding monomer counting noise), and the per-feed split of P1.
  The base-case tables (1-3), predictions and letters above are unchanged by it.
------------------------------------------------------------------------------------------------------------------------
PREDICTIONS for rna_seed_scan.py (written before the first run; no gain, noise or e_min number from this script had been computed)

Definitions (what is predicted)
  seed ee0   = enantiomeric excess of the monomer FEED, m_in (1 +/- ee0)/2 (the "inherited excess"; rna_copying.run(e=ee0)).
  ee_final   = steady-state product excess (D-L)/(D+L) of full-length copies (rna_copying.run readout), deterministic.
  gain       = ee_final / ee0.
  noise      = chemical-Langevin counting noise on the two template concentrations (the variance used in rna_copying.ramp),
               system size Omega = N_setting / m_in molecules per mM, N_setting = c V 1000 N_A from scorecard.SETTINGS
               (pond 1.9e25, lake 6.0e29, ocean 8.0e38); linear-noise stationary variance s^2 of the readout.
  followed   = Delta = ee_final / s >= 2.   e_min = smallest seed followed.
  cells      = r in {0.40, 0.56, 0.98} x feed in {5, 20, 200} mM; LO's n = 3, c = 0.5, w = 0.02.

P1  The racemic state is stable (leading difference-mode |lambda| < 1) in all 9 cells; lambda grows with r.
P2  Gain is constant (linear regime) for ee0 <= 1e-4 to 3 digits; at ee0 = 0.1 it is below the linear gain at r = 0.98.
P3  Linear gain at 20 mM: r = 0.40: 3.5 (range 3-5); r = 0.56: 5.5 (4-9); r = 0.98: 150 (50-600).
    At 5 and 200 mM the gain is within a factor 5 of the 20 mM value for each r.
P4  e_min at 20 mM, pond: r = 0.40: 1e-10 (3e-11 - 3e-10); r = 0.56: 7e-11 (2e-11 - 2e-10); r = 0.98: 1.5e-11 (3e-12 - 1e-10).
    Lake = pond / 178, ocean = pond / 6.5e6 (square roots of the N ratios; exact by construction).
P5  Against scorecard e_need (pond 6e-12, lake 3e-14, ocean 3e-18):
    r = 0.40 and 0.56 at 20 mM: e_min > e_need in all three settings (misses by 5x to 50x).
    r = 0.98 at 20 mM: ocean meets e_need (e_min ~ 2e-18), pond and lake miss (50% confidence on the ocean).
    Implied R6 letter (Y = all three met, P = some, N = none): N at r = 0.40 and 0.56, P at r = 0.98.
P6  Seeding the TEMPLATES instead (rna_copying.seed_run semantics, no feed excess): the excess decays (final |ee_T| < ee0, gain ~ 0) for every r.
P7  Self-check: run(e = 0.2, 5000 cycles, 20 mM) at r = 0.396 reproduces arm_g_rates.log:40 (0.668); at r = 0.40 it is 0.67 +/- 0.01.
------------------------------------------------------------------------------------------------------------------------
"""
import hashlib
import time
import warnings

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import fsolve
from scipy.stats import norm

import rna_copying as rc
from ocean import NA
from scorecard import SETTINGS, spec

PRED_SHA = "14c5f4ea46961666d2b488f4dd90d634121afc82527fd027a9ab82e02fabb801"
RS, FEEDS = (0.40, 0.56, 0.98), (5.0, 20.0, 200.0)
SEEDS = [10.0 ** -i for i in range(1, 19)]
NL_MIN = 1e-5                       # nonlinear solve for ee0 >= this, linear response below


def readout(z, m_in, k):
    _, FD, FL = rc.step(z, m_in, **k)
    return (FD - FL) / (FD + FL)


def jac(z, m_in, k, e):
    J = np.empty((4, 4))
    for i in range(4):
        d = np.zeros(4)
        d[i] = 1e-6 * z[i]
        J[:, i] = (rc.step(z + d, m_in, e=e, **k)[0] - rc.step(z - d, m_in, e=e, **k)[0]) / (2 * d[i])
    return J


def linear(m_in, k, mono=False):
    """Racemic state, leading difference eigenvalue, linear gain, y* per unit feed excess, readout slope, noise s^2 * Omega.
    mono=True adds monomer counting noise (sensitivity only): at steady state in + out events = w m_in per species, independent of T."""
    z = rc.fixed_point(m_in, **k)
    ev, Jd = rc.mode(z, m_in, k)
    h = 1e-3                                                     # e enters step() linearly, so this is exact
    b = rc.PM.T @ (rc.step(z, m_in, e=h, **k)[0] - rc.step(z, m_in, e=-h, **k)[0]) / (2 * h)
    rho = np.empty(2)
    for i in range(2):
        d = 1e-6 * z[2 * i]
        rho[i] = (readout(z + d * rc.PM[:, i], m_in, k) - readout(z - d * rc.PM[:, i], m_in, k)) / (2 * d)
    ystar = np.linalg.solve(np.eye(2) - Jd, b)
    S = rc.attempt(z[2], z[3], k["n"], k["r"])[0]
    X, w = k["c"] * z[0] * S, k["w"]
    v = X + w * (z[0] + X) + w * rc.T_IN / 2                     # per-cycle variance x Omega of each template count
    Sig = solve_discrete_lyapunov(Jd, np.diag([v, w * m_in if mono else 0.0])).real    # noise only on (T_D - T_L)/sqrt 2
    return z, abs(ev).max(), rho @ ystar, ystar, rho, rho @ Sig @ rho


def steady(zg, m_in, k, e):
    f = lambda z: rc.step(z, m_in, e=e, **k)[0] - z
    with warnings.catch_warnings():                  # 'xtol too small' = already at machine precision; residual asserted below
        warnings.simplefilter("ignore", RuntimeWarning)
        z = fsolve(f, zg, xtol=1e-14)
    assert np.abs(f(z) / z).max() < 1e-9 and (z > 0).all(), "steady state not found"
    assert abs(np.linalg.eigvals(jac(z, m_in, k, e))).max() < 1, "steady state not an attractor"
    return z


def langevin(z, m_in, k, e0, omega, runs, burn, seed=0):
    """Ensemble of ramp()'s noise at fixed feed, started at the deterministic steady state z; readout after burn cycles."""
    rng = np.random.default_rng(seed)
    z, w = z[:, None] * np.ones(runs), k["w"]
    for _ in range(burn):
        z2, XD, XL = rc.step(z, m_in, e=e0, **k)
        for j, X in ((0, XD), (1, XL)):
            var = (X + w * (z[j] + X) + w * rc.T_IN / 2) / omega
            z2[j] = np.maximum(z2[j] + np.sqrt(var) * rng.standard_normal(runs), 0)
        z = z2
    return readout(z, m_in, k)


def template_seed(m_in, k, ee0, cycles=20000):
    """rna_copying.seed_run semantics on the 4-variable map: excess ee0 in T_D vs T_L, no feed excess."""
    z = rc.fixed_point(m_in, **k) * [1 + ee0, 1 - ee0, 1, 1]
    tail = []
    for i in range(cycles):
        z = rc.step(z, m_in, **k)[0]
        if i >= cycles - 200:
            tail.append(abs(z[0] - z[1]) / (z[0] + z[1]))
    return max(tail)


def main():
    t0 = time.time()
    assert hashlib.sha256(open("rna_seed_scan_predictions.txt", "rb").read()).hexdigest() == PRED_SHA, "predictions file changed"
    names = list(SETTINGS)
    need = {n: e for n, _, e in spec()}
    N = {n: c * V * 1000 * NA for n, (c, tau, V) in SETTINGS.items()}
    print("Arm G seed scan: feed excess ee0 -> product excess; does it follow ee0 >= e_need under counting noise at N?")
    print("settings: " + "; ".join(f"{n}: N = {N[n]:.2e}, e_need = {need[n]:.0e}" for n in names))
    print(f"model: rna_copying LO with r = {RS}, n = 3, c = 0.5, w = 0.02, feed {FEEDS} mM (20 mM = main runs, 5 mM = below m*, 200 mM = activator scan)")

    # ---- self-check 1: the existing plateau (arm_g_rates.log:40, +20% feed, 5000 cycles, feed-filled reactor) ----
    plateau = {}
    for r in (0.396, 0.40):
        k = dict(rc.LO, r=r)
        z0 = np.array([rc.T_IN / 2, rc.T_IN / 2, 20.0 * 1.2 / 2, 20.0 * 0.8 / 2])
        plateau[r] = rc.run(z0, 20.0, k, 5000, e=0.2)[1][-1]
    print(f"\nself-check: +20% feed, 20 mM, 5000 cycles: r = 0.396 -> {plateau[0.396]:.4f} (arm_g_rates.log:40: 0.668), "
          f"r = 0.40 -> {plateau[0.40]:.4f} (0.67)")
    assert abs(plateau[0.396] - 0.668) < 1e-3 and abs(plateau[0.40] - 0.67) < 0.01, "plateau not reproduced"

    # ---- scan ----
    cells, gain = {}, {}
    for r in RS:
        for m_in in FEEDS:
            k = dict(rc.LO, r=r)
            z, lam, G, ystar, rho, su2 = linear(m_in, k)
            assert lam < 1, f"racemic state unstable at r={r}, feed {m_in}"
            g, zg = {}, z + rc.PM @ (ystar * NL_MIN)
            for e0 in sorted(s for s in SEEDS if s >= NL_MIN):          # continuation upward from the linear regime
                zg = steady(zg, m_in, k, e0)
                g[e0] = readout(zg, m_in, k) / e0
            assert abs(g[NL_MIN] / G - 1) < 1e-3, "nonlinear and linear gain disagree at 1e-5"
            gain[r, m_in] = {e0: g.get(e0, G) for e0 in SEEDS}
            cells[r, m_in] = dict(lam=lam, G=G, su2=su2, T=z[0], k=k, z=z, ystar=ystar)
    # nonlinear solve vs direct iteration of the map (independent of fsolve)
    for r, m_in, e0, cyc in ((0.40, 20.0, 0.1, 3000), (0.98, 20.0, 0.1, 120000)):
        k = dict(rc.LO, r=r)
        zi, ee = rc.run(cells[r, m_in]["z"], m_in, k, cyc, e=e0)
        print(f"check: direct iteration r = {r}, {m_in:g} mM, ee0 = {e0}, {cyc} cycles: gain {ee[-1] / e0:.6f} vs fixed point {gain[r, m_in][e0]:.6f}")
        assert abs(ee[-1] / e0 / gain[r, m_in][e0] - 1) < 1e-5, "direct iteration disagrees with fixed point"

    print("\nTABLE 1. deterministic gain ee_final/ee0 (ee0 >= 1e-5 nonlinear fixed point; below: exact linear response)")
    print(f"{'r':>5} {'feed':>5} {'|lam|':>8} {'T* uM':>7} " + " ".join(f"{e0:>8.0e}" for e0 in SEEDS))
    for (r, m_in), c in cells.items():
        print(f"{r:5.2f} {m_in:5g} {c['lam']:8.5f} {c['T'] * 1e3:7.2f} " + " ".join(f"{gain[r, m_in][e0]:8.3g}" for e0 in SEEDS))
    print("    gain is the product excess per unit inherited excess; it is bounded (racemic stable, r < 1), no amplification to 1")

    print("\nTABLE 2. noise s of the readout (stationary, linear-noise), Delta = ee_final/s, smallest seed followed (Delta >= 2)")
    print(f"    e_min = 2 s/G exactly (linear regime); grid = smallest of 1e-1 ... 1e-18 with Delta >= 2; 'x need' = e_min/e_need")
    hdr = "".join(f"{n.split()[-1]:>10} e_min {'x need':>8} {'grid':>6} " for n in names)
    print(f"{'r':>5} {'feed':>5} {'G':>8}  {hdr}  letter")
    print("    (grid ends at 1e-18: '<=1e-18' means e_min is below the grid)")
    letters = {}
    for (r, m_in), c in cells.items():
        out, ok = [], []
        for n in names:
            s = np.sqrt(c["su2"] / (N[n] / m_in))
            c[n] = dict(s=s, emin=2 * s / c["G"], dneed=c["G"] * need[n] / s)
            grid = [e0 for e0 in SEEDS if gain[r, m_in][e0] * e0 / s >= 2]
            c[n]["grid"] = min(grid) if grid else None
            ok.append(c[n]["emin"] <= need[n])
            out.append(f"{'':>10} {c[n]['emin']:7.1e} {c[n]['emin'] / need[n]:8.2g} {'-' if not grid else ('<=' if min(grid) == SEEDS[-1] else '') + f'{min(grid):.0e}':>6} ")
        letters[r, m_in] = "Y" if all(ok) else "P" if any(ok) else "N"
        print(f"{r:5.2f} {m_in:5g} {c['G']:8.3g}  {''.join(out)}  {letters[r, m_in]}")

    print("\nTABLE 3. seeded hand at e0 = e_need: Delta and P(win) = Phi(Delta) (a coin flip is 0.5, followed reliably is >= 0.977)")
    print(f"{'r':>5} {'feed':>5}  " + "  ".join(f"{n.split()[-1]:>8} Delta  P(win)" for n in names))
    for (r, m_in), c in cells.items():
        print(f"{r:5.2f} {m_in:5g}  " + "  ".join(f"{'':>8} {c[n]['dneed']:5.2f}  {norm.cdf(c[n]['dneed']):6.3f}" for n in names))

    # ---- sensitivity: the assumed self-copying efficiency c, and monomer counting noise (added after the first run) ----
    print("\nTABLE 4. sensitivity of e_min/e_need (pond | ocean; lake equals pond) to the assumed c (docstring: scanned 0.1-1, not measured)")
    print("    and to adding monomer counting noise (independent of template noise; variance w m_in/Omega per species). T* = racemic template level.")
    print(f"{'r':>5} {'feed':>5} {'c':>5} {'T* uM':>8} {'G':>7} {'|lam|':>8}  {'template-only: pond':>20} {'ocean':>8} {'L':>2}  {'+monomer noise: pond':>21} {'ocean':>8} {'L':>2}")
    letter = lambda rat: "Y" if all(x <= 1 for x in rat) else "P" if any(x <= 1 for x in rat) else "N"
    sens, bad = {}, []
    for r in RS:
        for m_in in FEEDS:
            for c in (0.1, 0.2, 0.5, 1.0):
                k, line = dict(rc.LO, r=r, c=c), f"{r:5.2f} {m_in:5g} {c:5.2f} "
                try:
                    res = [linear(m_in, k, mono) for mono in (False, True)]
                except ValueError:
                    bad.append((r, m_in, c))
                    continue
                if res[0][1] >= 1:
                    bad.append((r, m_in, c))
                    continue
                z, lam, G = res[0][:3]
                rat = [[2 * np.sqrt(x[5] / (N[n] / m_in)) / x[2] / need[n] for n in names] for x in res]
                sens[r, m_in, c] = (letter(rat[0]), letter(rat[1]), rat)
                print(line + f"{z[0] * 1e3:8.2f} {G:7.3g} {lam:8.5f}  {rat[0][0]:20.3g} {rat[0][2]:8.2g} {letter(rat[0]):>2}  {rat[1][0]:21.3g} {rat[1][2]:8.2g} {letter(rat[1]):>2}")
    print(f"    cells with no stable racemic state or no fixed point (skipped): {bad if bad else 'none'}")
    for c in (0.1, 0.2, 0.5, 1.0):
        row = [v for (r, f, cc), v in sens.items() if cc == c]
        print(f"    c = {c}: letters over r x feed, template-only {''.join(v[0] for v in row)}, +monomer noise {''.join(v[1] for v in row)}")

    # ---- validation of the linear-noise s with a Langevin ensemble (ramp()'s noise) at a small Omega ----
    print("\nVALIDATION. Langevin ensemble vs linear-noise theory, 20 mM, e0 = 1e-4, Omega chosen so that Delta_LNA = 1")
    print(f"{'r':>5} {'Omega/mM':>10} {'burn':>6} {'mean sim':>10} {'mean th':>10} {'sd sim':>10} {'sd th':>10} {'P sim':>6} {'P th':>6}")
    for r in RS:
        c, e0, runs = cells[r, 20.0], 1e-4, 3000
        omega = c["su2"] / (c["G"] * e0) ** 2
        burn = int(4 / (1 - c["lam"]))
        zs = steady(c["z"] + rc.PM @ (c["ystar"] * e0), 20.0, c["k"], e0)
        x = langevin(zs, 20.0, c["k"], e0, omega, runs, burn)
        sd_th = np.sqrt(c["su2"] / omega)
        print(f"{r:5.2f} {omega:10.3g} {burn:6d} {x.mean():10.3e} {c['G'] * e0:10.3e} {x.std():10.3e} {sd_th:10.3e} "
              f"{(x > 0).mean():6.3f} {norm.cdf(1.0):6.3f}")
        assert abs(x.mean() - c["G"] * e0) < 4 * sd_th / np.sqrt(runs) and abs(x.std() / sd_th - 1) < 0.08
        assert abs((x > 0).mean() - norm.cdf(1.0)) < 0.03, "linear-noise Delta does not match the Langevin ensemble"

    # ---- templates seeded instead of the feed (seed_run semantics) ----
    print("\nTEMPLATE SEED (seed_run semantics: excess ee0 in T_D vs T_L, no feed excess), 20 mM, |ee_T| over last 200 of 20000 cycles")
    ts = {}
    for r in RS:
        ts[r] = [template_seed(20.0, dict(rc.LO, r=r), e0) for e0 in (1e-1, 1e-3)]
        print(f"  r = {r:4.2f}: ee0 = 0.1 -> {ts[r][0]:.2e}, ee0 = 1e-3 -> {ts[r][1]:.2e}   (< ee0: decays)")
        assert ts[r][0] < 1e-1 and ts[r][1] < 1e-3, "template seed should decay for r < 1"

    # ---- predictions vs result ----
    print("\nPREDICTIONS vs RESULT")
    G20 = {r: cells[r, 20.0]["G"] for r in RS}
    chk = lambda ok: "held" if ok else "FAILED"
    mono_up = {f: cells[0.40, f]["lam"] < cells[0.56, f]["lam"] < cells[0.98, f]["lam"] for f in FEEDS}
    print(f"  P1 all 9 cells stable (|lam| < 1, asserted): held; |lam| rises with r: {chk(all(mono_up.values()))} "
          f"(per feed {', '.join(f'{f:g} mM: ' + ('yes' if mono_up[f] else 'NO') for f in FEEDS)}; at 20 mM: "
          + ", ".join(f"{cells[r, 20.0]['lam']:.5f}" for r in RS) + ")")
    lin_ok = all(abs(gain[r, f][e0] / cells[r, f]["G"] - 1) < 1e-3 for r in RS for f in FEEDS for e0 in SEEDS if e0 <= 1e-4)
    print(f"  P2 gain constant to 3 digits for ee0 <= 1e-4: {chk(lin_ok)}; at ee0 = 0.1, r = 0.98, 20 mM: gain {gain[0.98, 20.0][0.1]:.3g} vs linear {G20[0.98]:.3g}: {chk(gain[0.98, 20.0][0.1] < G20[0.98])}")
    rng_ = {0.40: (3, 5), 0.56: (4, 9), 0.98: (50, 600)}
    print("  P3 linear gain at 20 mM: " + "; ".join(f"r = {r}: {G20[r]:.3g} (pred {rng_[r][0]}-{rng_[r][1]}: {chk(rng_[r][0] <= G20[r] <= rng_[r][1])})" for r in RS))
    f5 = all(0.2 <= cells[r, f]["G"] / G20[r] <= 5 for r in RS for f in (5.0, 200.0))
    print(f"     5 and 200 mM within a factor 5 of 20 mM: {chk(f5)} (ratios " + ", ".join(f"{cells[r, f]['G'] / G20[r]:.2g}" for r in RS for f in (5.0, 200.0)) + ")")
    p4 = {0.40: (3e-11, 3e-10), 0.56: (2e-11, 2e-10), 0.98: (3e-12, 1e-10)}
    pond = names[0]
    print("  P4 e_min, pond, 20 mM: " + "; ".join(f"r = {r}: {cells[r, 20.0][pond]['emin']:.1e} (pred {p4[r][0]:.0e}-{p4[r][1]:.0e}: {chk(p4[r][0] <= cells[r, 20.0][pond]['emin'] <= p4[r][1])})" for r in RS))
    miss = all(cells[r, 20.0][n]["emin"] > need[n] for r in (0.40, 0.56) for n in names)
    oc = cells[0.98, 20.0][names[2]]["emin"] <= need[names[2]]
    print(f"  P5 r = 0.40, 0.56 miss e_need in all three at 20 mM: {chk(miss)}; r = 0.98 ocean meets it: {chk(oc)}; "
          f"r = 0.98 pond and lake miss: {chk(all(cells[0.98, 20.0][n]['emin'] > need[n] for n in names[:2]))}")
    print("     letters at 20 mM: " + ", ".join(f"r = {r}: {letters[r, 20.0]}" for r in RS) + "  (pred N, N, P)")
    print(f"  P6 template seed decays for every r: {chk(all(ts[r][0] < 1e-1 and ts[r][1] < 1e-3 for r in RS))}")
    print(f"  P7 plateau reproduced: {chk(True)} (self-check asserted above)")
    print("\nALL LETTERS at the base c = 0.5, template-only noise (r/feed): " + ", ".join(f"{r}/{f:g}: {letters[r, f]}" for r in RS for f in FEEDS))
    print(f"\nself-check passed; runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
