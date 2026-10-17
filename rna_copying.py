"""Arm G: nonenzymatic template-directed RNA copying (2-aminoimidazole-activated
nucleotides) in a fed, well-mixed reactor with temperature cycling.

Per-cycle map (one cycle = cold copying phase, then a hot phase that separates strands):
    each primer-template duplex attempts a copy with probability c   (self-copying
        efficiency: NOT shown experimentally, so scanned)
    an attempt adds monomers one at a time: N ~ Poisson(mu) addition events, mu set by the
        measured extension rate; each event is right-handed with probability
        p = m_R / (m_R + r m_W); a wrong-handed monomer at the 3' end ends the strand
        (Joyce 1984; Essex et al. 2026); n right additions in a row = a full-length copy
    full-length copies become templates after the hot phase; partial/terminated ones do not
    T_D' = (1-w) T_D (1 + c S_D) + w T_in/2
    m_D' = (1-w) [(1-h) m_D - c (T_D E_R,D + T_L E_W,L)] + w m_in (1+e)/2      (L: swap D, L)
    S = p^n P(N>=n);  E_R, E_W = right / wrong monomers used per attempt.

Exact identity: E_R/E_W = p/(1-p) for any n, any truncation (each event before the stop is
drawn the same way). At 50/50 monomer this is 1/r. So a template uses more of its OWN hand's
monomer than the other's iff r < 1.

Difference mode (x = T_D - T_L, y = m_D - m_L) at the racemic fixed point:
    J_diff = [[A, B], [C, Dd]],  B > 0 (own-monomer surplus helps),
    C = -(1-w) c (E_R - E_W)  ->  C < 0 for r < 1, C > 0 for r > 1
    loses stability (eigenvalue +1) when (1-A)(1-Dd) = B C.
So r > 1 (cross-competition beats self-competition, Lotka-Volterra) is NECESSARY.
For r < 1 each hand eats its own monomer: coexistence, the racemic state is stable at every feed.

Inputs.
MEASURED (Walton & Szostak 2017, Biochemistry 56:5739, doi:10.1021/acs.biochem.7b00792):
    k4   = 3.33e-3 /h  activated monomer (Cp*) hydrolysis, Table 1 (upper limit, t1/2 > 8.7 d)
    kobs = 5.0-9.7 /h  primer extension with 24 mM Cp*, Fig. 8a (midpoint 7.35 used)
    T_in = 3 uM        template concentration used there (as a FEED level it is assumed)
QUALITATIVE ONLY (Essex, Jia, Zhang & Szostak 2026, ChemRxiv, doi:10.26434/chemrxiv.15005970/v1):
    L/L bridged dinucleotides extend D primer-templates "much more slowly" than D/D, and an
    L nucleotide at a D 3' end is a chain terminator -> r < 1. No number is given, so r is scanned.
ASSUMED (scanned or stated):
    r      wrong/right addition rate ratio     0.1 (measured direction) and 1.2; scan 0.1-2
    c      self-copying efficiency per cycle  0.5; scan 0.1-1
    n      nucleotides to full length          3; scan 3-6
    w      exchange (wash) fraction per cycle  0.02; 0.005 and 0.05 also run
    t_cold 1 h copying phase per cycle (h = 1 - exp(-k4 t), mu = kobs t m_tot / 24 mM)
    extension first order in monomer (real chemistry goes through the bridged dinucleotide)
    noise  template counting noise only (monomers are mM pools)

RESULT (rna_copying.log):
 (a) +20% D and +20% L give exactly equal and opposite product excess. It rises over cycles to a
     bounded plateau (r = 0.1: 0.35 from a 0.20 feed excess, after a damped swing), and is amplified
     to ~1 only past m* with r > 1.
 (b) The racemic state is lost only if r > 1 AND c p^n P(N>=n) > w/(1-w) (racemic copying outpaces
     washout), so 1 < r < (c(1-w)/w)^(1/n) - 1 is necessary. n = 3, c = 0.5, w = 0.02, r = 1.2:
     m* = 9.9 mM feed (T_in -> 0 formula 8.3 mM). r <= 1: stable at every feed and c scanned.
 (c) Inside that window: supercritical pitchfork in feed, no hysteresis, slow tip-over; the
     normal form holds (slow sweeps within 0.02) with F_RNA = 0.58 vs frank.py's 0.70.
     For r < 1, the direction the chemistry points, there is nothing to cross: R9 not met.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import gammainc
from scipy.stats import norm
from sim import delta
import frank

K4, KOBS, M_REF, T_IN = 3.33e-3, 7.35, 24.0, 0.003      # /h, /h, mM, mM (measured, see above)
T_COLD = 1.0                                            # h, ASSUMED
H, A = 1 - np.exp(-K4 * T_COLD), KOBS * T_COLD
LO = dict(n=3, r=0.1, c=0.5, w=0.02)                    # r < 1: direction the chemistry points
HI = dict(n=3, r=1.2, c=0.5, w=0.02)                    # r > 1: hypothetical
PM = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]]) / np.sqrt(2)   # difference subspace


def attempt(mr, mw, n, r):
    """Full-length probability S and right/wrong monomers used, per copy attempt."""
    mr, mw = np.maximum(mr, 1e-12), np.maximum(mw, 1e-12)
    tot = mr + r * mw
    p, mu = mr / tot, A * tot / M_REF
    er = ew = 0.0
    pj = 1.0
    for j in range(1, n + 1):
        t = pj * gammainc(j, mu)              # event j happens and events 1..j-1 were right
        er, ew, pj = er + t * p, ew + t * (1 - p), pj * p
    return pj * gammainc(n, mu), er, ew


def step(z, m_in, n, r, c, w, e=0.0, g=0.0):
    """One cycle. g = parity-violating bias on every D addition (+g) and L addition (-g)."""
    TD, TL, mD, mL = z
    SD, RD, WD = attempt((1 + g) * mD, (1 - g) * mL, n, r)
    SL, RL, WL = attempt((1 - g) * mL, (1 + g) * mD, n, r)
    z2 = np.array([(1 - w) * TD * (1 + c * SD) + w * T_IN / 2,
                   (1 - w) * TL * (1 + c * SL) + w * T_IN / 2,
                   (1 - w) * ((1 - H) * mD - c * (TD * RD + TL * WL)) + w * m_in * (1 + e) / 2,
                   (1 - w) * ((1 - H) * mL - c * (TL * RL + TD * WD)) + w * m_in * (1 - e) / 2])
    return z2, c * TD * SD, c * TL * SL


def fixed_point(m_in, n, r, c, w):
    """Racemic fixed point [T, T, m, m]."""
    S = lambda m: attempt(m, m, n, r)[0]
    T = lambda m: w * T_IN / 2 / (1 - (1 - w) * (1 + c * S(m)))
    hi = w * m_in / 2 / (1 - (1 - w) * (1 - H))           # no-copying monomer level
    grow = lambda m: (1 - w) * (1 + c * S(m)) - 1
    if grow(hi) > 0:                                      # copying outruns washout: T unbounded above here
        hi = brentq(grow, 1e-12, hi)
    f = lambda m: (1 - w) * ((1 - H) * m - c * T(m) * sum(attempt(m, m, n, r)[1:])) + w * m_in / 2 - m
    m = brentq(f, 1e-9 * hi, hi * (1 - 1e-12), xtol=1e-14)
    return np.array([T(m), T(m), m, m])


def mode(z, m_in, k):
    """Difference-mode block J_diff = [[A, B], [C, Dd]] at a racemic state, and its eigenvalues."""
    J = np.empty((4, 4))
    for i in range(4):
        d = np.zeros(4)
        d[i] = 1e-6 * z[i]
        J[:, i] = (step(z + d, m_in, **k)[0] - step(z - d, m_in, **k)[0]) / (2 * d[i])
    Jd = PM.T @ J @ PM
    return np.linalg.eigvals(Jd), Jd


def onset(k):
    """T_in -> 0 boundary: racemic copying just keeps up with washout, c S(m) = w/(1-w)."""
    n, r, c, w = k["n"], k["r"], k["c"], k["w"]
    if c * (1 + r) ** -n <= w / (1 - w):
        return None
    m = brentq(lambda m: c * attempt(m, m, n, r)[0] - w / (1 - w), 1e-9, 1e6)
    return 2 * m * (1 - (1 - w) * (1 - H)) / w


def threshold(k, grid=np.geomspace(0.5, 100, 40)):
    """Computed boundary: feed m_in (mM) where the leading difference eigenvalue crosses +1.
    Also returns max |eigenvalue| and min Dd over the scan (Dd > 0: no per-cycle monomer overshoot)."""
    blocks = [mode(fixed_point(m, **k), m, k) for m in grid]
    lead = lambda m: mode(fixed_point(m, **k), m, k)[0].real.max() - 1
    vals = [b[0].real.max() - 1 for b in blocks]
    mods, dd = max(abs(b[0]).max() for b in blocks), min(b[1][1, 1] for b in blocks)
    for i in range(1, len(grid)):
        if vals[i - 1] < 0 < vals[i]:
            return brentq(lead, grid[i - 1], grid[i], xtol=1e-6), mods, dd
    return None, mods, dd


def run(z, m_in, k, cycles, e=0.0):
    ee = []
    for _ in range(cycles):
        z, FD, FL = step(z, m_in, e=e, **k)
        ee.append((FD - FL) / (FD + FL))
    return z, np.array(ee)


def ramp(g, omega, tau, k, mstar, runs=2000, seed=0):
    """Chemical-Langevin ramp of m_in from 0.8 to 1.6 m*. Returns P(D wins) and normal-form prediction."""
    rng = np.random.default_rng(seed)
    w = k["w"]
    z = fixed_point(0.8 * mstar, **k)[:, None] * np.ones(runs)
    feed = mstar * (0.8 + 0.8 * np.arange(tau) / tau)
    Tbar, mbar = np.empty(tau), np.empty(tau)
    for i in range(tau):
        Tbar[i], mbar[i] = np.median(z[:2].mean(0)), z[2:].mean()
        z2, XD, XL = step(z, feed[i], g=g, **k)
        for j, X in ((0, XD), (1, XL)):                 # copies + washout + feed, counted per template
            var = (X + w * (z[j] + X) + w * T_IN / 2) / omega
            z2[j] = np.maximum(z2[j] + np.sqrt(var) * rng.standard_normal(runs), 0)
        z = z2
    pick = np.arange(0, tau, max(tau // 400, 1))
    lam = np.array([mode(np.array([Tbar[i], Tbar[i], mbar[i], mbar[i]]), feed[i], k)[0].real.max() for i in pick])
    i = pick[np.argmax(lam >= 1)]
    s = max(tau // 50, 2)                                 # slope over +-2% of the ramp, as frank.py
    sel = (pick >= i - s) & (pick <= i + s)
    gamma = np.polyfit(pick[sel], lam[sel], 1)[0]
    zc = np.array([Tbar[i], Tbar[i], mbar[i], mbar[i]])
    lam, Jd = mode(zc, feed[i], k)
    vals, V = np.linalg.eig(Jd)
    lu, U = np.linalg.eig(Jd.T)
    v, u = PM @ V[:, np.argmax(vals.real)].real, PM @ U[:, np.argmax(lu.real)].real
    u = u / (u @ v)                                       # left eigenvector, u.v = 1
    b = (step(zc, feed[i], g=g, **k)[0] - step(zc, feed[i], g=-g, **k)[0]) / 2
    _, X, _ = step(zc, feed[i], **k)
    var = (X + w * (zc[0] + X) + w * T_IN / 2) / omega
    d_real = delta(u @ b, (u[0] ** 2 + u[1] ** 2) * var, gamma)
    S = [attempt((1 + s_) * zc[2], (1 - s_) * zc[3], k["n"], k["r"])[0] for s_ in (g, -g)]
    g_strand, kc = (S[0] - S[1]) / (S[0] + S[1]), (1 - w) * k["c"] * (S[0] + S[1]) / 2
    d_simple = delta(g_strand * kc, kc / (2 * zc[0] * omega), gamma)
    ratio, amp = (d_real / d_simple, g_strand / g) if g else (np.nan, np.nan)
    return (z[0] > z[1]).mean(), norm.cdf(d_real), ratio, amp, zc[0]


def main():
    print("Arm G: nonenzymatic RNA copying (2-AI monomers), fed reactor, temperature cycling")
    print(f"measured: k4 = {K4} /h, kobs = {KOBS} /h (midpoint of 5.0-9.7) at 24 mM, T_in = {T_IN * 1e3:.0f} uM (Walton & Szostak 2017)")
    print(f"-> h = {H:.4f} per cycle, mu = {A:.2f} addition events per 1 h cold phase at 24 mM")
    print("assumed: r, c, n, w, 1 h cold phase, first-order extension, template-only noise (see docstring)")
    print(f"LO = {LO}   HI = {HI}")
    p, n, r = 0.7, 4, 0.37
    _, er, ew = attempt(p, (1 - p) / r, n, r)            # p = m_R/(m_R + r m_W) = 0.7
    print(f"\nidentity E_R/E_W = p/(1-p): {er / ew:.12f} vs {p / (1 - p):.12f} (truncated Poisson, n = {n})")

    print("\n(a) seed-following. Readout: excess of full-length D over L product per cycle, (D-L)/(D+L)")
    print("    reactor starts filled with feed; e = monomer feed excess")
    cyc = [1, 10, 30, 100, 300, 1000, 5000]
    print(f"{'case':>22} {'feed':>7} " + " ".join(f"{c:>7}" for c in cyc))
    sym = True
    for name, k, m_in in (("r=0.1, 20 mM", LO, 20.0), ("r=1.2, 5 mM (< m*)", HI, 5.0), ("r=1.2, 20 mM (> m*)", HI, 20.0)):
        out = {}
        for e, lab in ((0.2, "+20% D"), (-0.2, "+20% L"), (0.0, "0%")):
            z0 = np.array([T_IN / 2, T_IN / 2, m_in * (1 + e) / 2, m_in * (1 - e) / 2])
            out[e] = run(z0, m_in, k, cyc[-1], e=e)[1]
            print(f"{name:>22} {lab:>7} " + " ".join(f"{out[e][c - 1]:7.3f}" for c in cyc))
        sym &= np.allclose(out[0.2], -out[-0.2], atol=1e-12) and np.all(out[0.0] == 0)
    print(f"    +20% D and +20% L give equal and opposite excess at every cycle: {sym}")
    print("    r = 0.1: excess rises (with a damped cycle-to-cycle swing) to a bounded plateau set by the feed")
    print("    r = 1.2 above m*: the seed is amplified to near-homochiral; below m*: bounded, like r < 1")
    print("    0% stays exactly 0 here (no noise); with noise, (c) shows who wins")

    print("\n(b) racemic stability. J_diff = [[A, B], [C, Dd]] at the racemic fixed point, m_in = 20 mM")
    print("    lost when (1-A)(1-Dd) = B*C;  C = -(1-w) c (E_R - E_W), E_R/E_W = 1/r  ->  needs r > 1")
    print("    T_in -> 0 boundary (onset): c * S(m*) = w/(1-w), S = p^n P(N>=n), p = 1/(1+r),")
    print("        mu = kobs t_cold m (1+r) / 24 mM, m = w m_in / 2 / (1 - (1-w)(1-h))")
    print("    necessary window (S <= p^n):  1 < r < (c (1-w)/w)^(1/n) - 1")
    print(f"{'n':>3} {'r':>5} {'c':>5} {'w':>6} {'A':>8} {'B':>8} {'C':>8} {'Dd':>7} {'r_max':>6} "
          f"{'onset':>7} {'m* (mM)':>8} {'max|lam|':>8}")
    rows = [dict(LO, r=r) for r in (0.1, 0.5, 0.9, 1.0, 1.1, 1.2, 1.5, 1.8, 2.0)]
    rows += [dict(HI, c=c) for c in (0.1, 0.2, 0.3, 1.0)] + [dict(HI, c=1.0, n=n) for n in (4, 5, 6)]
    rows += [dict(HI, w=w) for w in (0.005, 0.05)] + [dict(LO, c=c) for c in (0.1, 1.0)]
    res, dd_min = {}, np.inf
    for k in rows:
        _, Jd = mode(fixed_point(20.0, **k), 20.0, k)
        ms, mods, dd = threshold(k)
        dd_min = min(dd_min, dd)
        res[tuple(k.values())] = ms
        rmax = (k["c"] * (1 - k["w"]) / k["w"]) ** (1 / k["n"]) - 1
        on = onset(k)
        print(f"{k['n']:3d} {k['r']:5.2f} {k['c']:5.2f} {k['w']:6.3f} {Jd[0, 0]:8.5f} {Jd[0, 1]:8.5f} {Jd[1, 0]:8.4f} "
              f"{Jd[1, 1]:7.4f} {rmax:6.2f} {'-' if on is None else f'{on:7.2f}':>7} "
              f"{'none' if ms is None else f'{ms:8.2f}':>8} {mods:8.4f}")
    print(f"    feed scanned 0.5-100 mM; min Dd over all scans = {dd_min:.3f} (> 0: per-cycle map valid)")
    print("    r <= 1: no crossing at any feed or c. r > 1: crossing only inside the window, just above")
    print("    the onset (offset by the fed templates T_in, the k0 of frank.py).")
    lo_stable = all(v is None for key, v in res.items() if key[1] <= 1)
    mstar = res[tuple(HI.values())]

    print(f"\n(c) R9 and the normal form. r = 1.2 (HI): m* = {mstar:.3f} mM")
    print(f"{'m_in':>7} " + " ".join(f"{m:>15}" for m in ("lam r=0.1", "lam r=1.2")))
    for m in (2, 5, 8, 9.5, 10.5, 12, 20, 50, 100):
        lams = [mode(fixed_point(m, **k), m, k)[0] for k in (LO, HI)]
        print(f"{m:7.1f} " + " ".join(f"{l[np.argmax(l.real)].real:9.5f}{'(c)' if l.imag.any() else '   ':>6}" for l in lams))
    print("    (c) = complex pair: damped cycle-to-cycle D/L swing, not symmetry breaking")
    print("    steady state just above m* (deterministic, tiny seed): supercritical if ee ~ sqrt(d)")
    for d in (0.005, 0.01, 0.02, 0.05):
        m = mstar * (1 + d)
        ee = run(fixed_point(m, **HI) * np.array([1.001, 0.999, 1, 1]), m, HI, 200000)[1][-1]
        print(f"      d = (m_in - m*)/m* = {d:5.3f}: ee = {ee:.3f}, ee/sqrt(d) = {ee / np.sqrt(d):.2f}")
    z = run(fixed_point(2 * mstar, **HI) * np.array([1.5, 0.5, 1, 1]), 2 * mstar, HI, 60000)[0]
    ee_back = run(z, 0.97 * mstar, HI, 200000)[1][-1]
    print(f"    homochiral state at 2 m*, feed dropped to 0.97 m*: ee -> {ee_back:.1e} (no hysteresis)")
    print("    growth of the difference mode above m*: lam - 1 <= 0.006 per cycle -> tips over in ~10^2-10^3 cycles")

    print("\n    chemical-Langevin ramp of m_in from 0.8 m* to 1.6 m* over tau cycles, bias g on every addition")
    print("    prediction: P = Phi(Delta), Delta from sim.delta with g_a = u.b, eps_a = u Sigma u, gamma = d(lam)/dt")
    print(f"{'g':>7} {'Omega':>7} {'tau':>6}  {'predicted':>9} {'simulated':>9} {'F_RNA':>6} {'g_strand/g':>10}")
    diffs, F = [], []
    for g, om, tau, runs in ((0.0, 1e5, 20000, 2000), (1e-3, 1e5, 2000, 2000), (1e-3, 1e5, 20000, 2000),
                             (5e-4, 4e5, 20000, 2000), (2e-3, 1e5, 20000, 2000), (1e-3, 1e5, 80000, 1000)):
        P, pred, f, amp, T = ramp(g, om, tau, HI, mstar, runs=runs)
        if tau >= 20000:                                  # slow sweeps only, as in frank.py
            diffs.append(abs(P - pred))
        if g:
            F.append(f)
        print(f"{g:7.4f} {om:7.0e} {tau:6d}  {pred:9.3f} {P:9.3f} {f:6.3f} {amp:10.2f}   (T at crossing {T * 1e3:.1f} uM)")
    x = np.sqrt(frank.K0 / frank.KI)
    Ff = np.sqrt(x / (frank.K0 + 2 * x))
    print(f"    slow sweeps: max |simulated - predicted| = {max(diffs):.3f}; fast sweep selects MORE (as in frank.py)")
    print(f"    F_RNA = Delta_real / Delta_simple = {np.mean(F):.2f}  vs frank.py F = {Ff:.2f}")
    print(f"      (Delta_simple: bias = per-strand copy bias, noise 1/N strands, same sweep, as frank.py/ocean.py)")
    print(f"      -> required concentrations rise by F^(-4/3) = {np.mean(F) ** (-4 / 3):.2f}x (frank.py: {Ff ** (-4 / 3):.2f}x)")
    print("      a per-nucleotide bias g becomes a per-strand bias ~2n(1-p) g = 3.1 g (n = 3, p = 1/2.2)")

    print("\nANSWERS")
    print(" (a) yes: +20% D and +20% L give equal and opposite product excess that rises with cycle number")
    print("     (r = 0.1 after a damped swing) to a bounded plateau, larger as m* is approached;")
    print("     it goes to ~1 only past m* with r > 1.")
    print(" (b) the racemic state is lost only if r > 1 (each hand uses MORE of the other hand's monomer,")
    print("     i.e. wrong-hand addition faster than right-hand) AND racemic copying outpaces washout:")
    print("     c p^n P(N>=n) > w/(1-w). Window 1 < r < (c(1-w)/w)^(1/n) - 1; m* from the table.")
    print("     Measured chemistry points to r < 1 (Essex 2026: L/L extends D primers much more slowly) ->")
    print("     racemic stable at every feed and c scanned: seed-following only.")
    print(" (c) R9: met only for r > 1 inside the window (supercritical pitchfork, slow tip-over, the")
    print(f"     normal form holds, F_RNA = {np.mean(F):.2f} vs 0.70). For r < 1: no instability to cross -> R9 not met.")
    assert (abs(er / ew - p / (1 - p)) < 1e-12 and sym and lo_stable and mstar is not None
            and max(diffs) < 0.03), "self-check failed"
    print("\nself-check passed")


if __name__ == "__main__":
    main()
