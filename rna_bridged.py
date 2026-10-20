"""Arm G with bridged-dinucleotide copying: can heterochiral bridges make 50/50 unstable at r < 1?

rna_copying.py extends primers first order in free monomer, so every addition draws from one pool and
E_R/E_W = p/(1-p): for r < 1 each template uses more of its own hand, C < 0 and 50/50 is stable at every feed.
Real 2-AI copying goes through imidazolium-bridged dinucleotides N*N (Walton & Szostak 2017): only the bridge
binds the template; the primer adds the +1 nucleotide and the +2 leaves as activated monomer N*. A heterochiral
bridge ties one D and one L together, so a D template adding through L(+1)/D(+2) eats an L and hands back a D.

State z = [T_D, T_L, m_D, m_L, B_DD, B_LL, B_DL, B_LD, a]. B_XY = bridge with hand X at +1 and Y at +2 as seen on a
D template; D*L and L*D are different molecules (heteronucleotide A*C, the case of Essex Fig. 4); the mirror swaps
B_DL and B_LD. a = free 2-AI. Per cycle (one cold phase, t = 1 h, as rna_copying.py), one ODE, fixed-step RK4:
  chemistry   2 N* -> N*N + AI  k1 [N*]^2, split by hand: k1 fo m_D^2 (DD), k1 fh m_D m_L (DL and LD each), k1 fo m_L^2,
                                fo = 2 hom/(het+hom), fh = 2 het/(het+hom)  (the weights of arm_g_rates.r_bridged)
              N*N + AI -> 2 N*  k2      N*N -> Np + N* (either end, 1/2 each)  k3      N* -> Np + AI  k4
  copying     concurrent: c T_i primers of hand-i templates, q_k = fraction still extending with k right (k < n),
              rate per primer by species s = (+1, +2) hand relative to its template, competitive Michaelis-Menten,
                v_s = g kmax_s [B_s]/Km_s / (1 + sum_j [B_j]/Km_j + (m_D + m_L)/Ki)
              q_k gains right events (s = own +1) from q_(k-1), loses all events; S_i = primers reaching n right.
              Each event removes bridge s for good and returns its +2 to the free monomer pool.
  T' = (1-w) T (1 + c S) + w T_in/2;  pools' = (1-w) pools(t) + w feed.
With pools held fixed this is rna_copying.attempt exactly (S = p^n P(N >= n), N ~ Poisson(sum v t)).
bridged=False (k1 = 0, extension first order in free monomer at rna_copying's rate, frozen at start-of-cycle pools as
there) reproduces rna_copying.step and its D-L block (self-check: to < 1e-3, leading |lam| to < 1e-5; the residue is
hydrolysis running alongside copying here instead of before it, O(k4 t x copying use / pool)).

D-L block: the mirror splits J; P projects on x = T_D - T_L, y = m_D - m_L, b1 = B_DD - B_LL, b2 = B_DL - B_LD (each
/ sqrt 2). At the racemic fixed point the script reports the 4x4 block's leading eigenvalue (> 1 in modulus: 50/50
lost; negative real = per-cycle overshoot) and two cross terms:
  C_m   = dy'/dx           free monomer, the C of rna_copying's J_diff = [[A, B], [C, Dd]]
  C_tot = d(y + 2 b1)'/dx  chiral inventory, free + bridged (D*L and L*D hold one of each)
An event spends only its +1 nucleotide and moves its +2 between pools, so while copying uses a small part of the
pools within the hour C_tot = -(1-w) c (E_R - E_W) as in rna_copying.py, negative iff r_fp < 1,
r_fp = (v_LD + v_LL)/(v_DD + v_DL) at the fixed point. When copying depletes bridges within the hour (g >> 1),
E_R, E_W of both hands respond to x as well, and r_fp drifts from r because depletion is uneven across species.

Inputs.
MEASURED, Walton & Szostak 2017, Biochemistry 56:5739, doi:10.1021/acs.biochem.7b00792 (C*C, Tris pH 8.3, 90 mM Mg):
    k1 = 4.49e-3 /h/mM  Table 1, Fig. 1 (2 Cp* -> Cp*pC + AI, d[Cp*pC]/dt = k1 [Cp*]^2; plot slope 2 k1)
    k2 = 0.238 /h/mM    Table 1, Fig. 2 (Cp*pC + AI -> 2 Cp*)
    k3 = 0.171 /h       Table 1, Fig. 3 (Cp*pC -> Cp + Cp*)
    k4 = 3.33e-3 /h     Table 1, Fig. 4 (Cp* -> Cp + AI, upper limit) = rna_copying.K4
    Ki = 24.7 mM        Fig. 5 (Cp inhibition of Cp*pC extension; Cp* taken as equal, as there)
MEASURED, Essex, Jia, Zhang & Szostak 2026, ChemRxiv, doi:10.26434/chemrxiv.15005970/v1, imported from arm_g_rates.py:
    AC    kmax, Km of the four A*C hand species (Fig. 4)  -> v_s; at 50/50 r0 = ag.r_bridged = 0.396
    FORM  hetero/homo bridge formation from 1:1 *C (Fig. S3) -> fh, fo (ratio only; absolute 2.4e-3 is ~k1/2)
    HYD   A*C hydrolysis D/L = D/D (Fig. S6) -> one k3 for every hand pair (ratio only; absolute 0.088 /h vs 0.171)
    T_in, t_cold, n = 3, w = 0.02 and rc.attempt's statistics as rna_copying.py.
ASSUMED (scanned or stated):
    r      0.06-1.2 scanned by scaling kmax of the +1-wrong species (L/D, L/L on a D template) by r/r0
    c      0.1-1 scanned; feed 0.5-100 mM scanned (16 points)
    g      copying speed-up, all kmax x g: 1 = measured; 10 and 30 scanned to reach copying faster than washout
    k1-k3 of C*C used for every bridge and hand pair (k2 hand-independent); Walton and Essex buffers differ
    no free 2-AI in the feed (it builds up to ~ free monomer); inactive Np not tracked (inhibits like N*, ignored)
    C*C on a GG template (one hetero molecule binding both ways) doubles the hetero weight: r = 0.76, in the scan
"""
import numpy as np
from scipy.optimize import brentq, fsolve

import arm_g_rates as ag
import rna_copying as rc

K1, K2, K3, KI = 4.49e-3, 0.238, 0.171, 24.7          # /h/mM, /h/mM, /h, mM  (Walton & Szostak 2017)
HET, HOM = ag.FORM["hetero"][0], ag.FORM["homo"][0]
FH, FO = 2 * HET / (HET + HOM), 2 * HOM / (HET + HOM)
SP = ("DD", "DL", "LD", "LL")                         # (+1, +2) hand as seen by a D template
KMAX, KMS = (np.array([ag.AC[s][i] for s in SP])[:, None] for i in (0, 2))
R0 = ag.r_bridged(ag.eff(ag.AC), HET, HOM)
R_CC = ag.r_bridged(ag.eff(ag.AC), 2 * HET, HOM)      # homonucleotide bridge: hetero molecule binds both ways
ON = ([2, 4, 5, 3], [3, 5, 4, 2])                     # pool index (z[2:]) of species SP as seen by a D / L template
NSUB = 40                                             # RK4 steps per 1 h cold phase, x g // 10 (stiff when fast)
P = np.kron(np.eye(5, 4), [[1], [-1]])[:9] / 2 ** .5  # (9, 4): D-L subspace [x, y, b1, b2] of z
SYM = [0, 0, 1, 1, 2, 2, 3, 3, 4]                     # racemic z from u = [T, m, B_homo, B_het, a]


def chem(y, k1):
    """d/dt of pools y = [m_D, m_L, B_DD, B_LL, B_DL, B_LD, a] from bridge chemistry alone."""
    mD, mL, dd, ll, dl, ld, a = y
    fdd, fll, fx, rev = k1 * FO * mD * mD, k1 * FO * mL * mL, k1 * FH * mD * mL, K2 * a
    back = (rev + K3 / 2) * (dl + ld)                 # each hetero bridge returns one of each (rev) or either (hyd)
    return np.array([-rc.K4 * mD - 2 * fdd - 2 * fx + (2 * rev + K3) * dd + back,
                     -rc.K4 * mL - 2 * fll - 2 * fx + (2 * rev + K3) * ll + back,
                     fdd - (K3 + rev) * dd, fll - (K3 + rev) * ll, fx - (K3 + rev) * dl, fx - (K3 + rev) * ld,
                     rc.K4 * (mD + mL) + fdd + fll + 2 * fx - rev * (dd + ll + dl + ld)])


def rates(y, i, r, g=1.0, bridged=True):
    """Extension rate (/h) of one active hand-i primer by species SP; first order in monomer if not bridged."""
    if bridged:
        B = y[ON[i]] / KMS
        return g * KMAX * np.array([1, 1, r / R0, r / R0])[:, None] * B / (1 + B.sum(0) + (y[0] + y[1]) / KI)
    v = np.zeros_like(y[:4])
    v[0], v[2] = y[i], r * y[1 - i]                   # right (own monomer), wrong (other monomer)
    return rc.A / rc.T_COLD / rc.M_REF * v


def rhs(Y, cT, n, r, g, bridged, y0):
    """Pools, then per hand: q_0..q_(n-1) and S. Not bridged: rates frozen at start pools y0 (rna_copying)."""
    y, out = Y[:7], []
    d = chem(y, K1 if bridged else 0.0)
    for i in (0, 1):
        q = Y[7 + i * (n + 1): 7 + i * (n + 1) + n]
        v = rates(y if bridged else y0, i, r, g, bridged)
        f = cT[i] * q.sum(0) * v                      # bridges (monomers) used per h by species
        if bridged:
            d[ON[i]] -= f
            d[i] += f[0] + f[2]                       # +2 of own hand (D/D, L/D on a D template)
            d[1 - i] += f[1] + f[3]
        else:
            d[i] -= f[0]
            d[1 - i] -= f[2]
        R = v[0] + v[1]
        dq = -v.sum(0) * q
        dq[1:] += R * q[:-1]
        out += [dq, R * q[-1:]]
    return np.vstack([d] + out)


def step(z, m_in, n, r, c, w, e=0.0, bridged=True, g=1.0):
    """One cycle; z is (9,) or (9, K) for K states at once. Returns z', full-length D and L copies made."""
    z = np.asarray(z, float)
    flat, z = z.ndim == 1, z.reshape(9, -1)
    q0 = np.zeros((n + 1, z.shape[1]))
    q0[0] = 1
    ns = NSUB * max(1, int(g) // 10)
    Y, cT, h = np.vstack([z[2:], q0, q0]), c * z[:2], rc.T_COLD / ns
    for _ in range(ns):
        k1 = rhs(Y, cT, n, r, g, bridged, z[2:])
        k2 = rhs(Y + h / 2 * k1, cT, n, r, g, bridged, z[2:])
        k3 = rhs(Y + h / 2 * k2, cT, n, r, g, bridged, z[2:])
        Y = Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + rhs(Y + h * k3, cT, n, r, g, bridged, z[2:]))
    S = Y[[7 + n, 8 + 2 * n]]
    feed = np.zeros_like(Y[:7])
    feed[0], feed[1] = m_in * (1 + e) / 2, m_in * (1 - e) / 2   # no free 2-AI fed
    z2 = np.vstack([(1 - w) * z[:2] * (1 + c * S) + w * rc.T_IN / 2, (1 - w) * Y[:7] + w * feed])
    F = c * z[:2] * S
    return (z2[:, 0], *F[:, 0]) if flat else (z2, *F)


def mode(z, m_in, k):
    """Eigenvalues of the D-L block Jd = P^T J P (central differences, all 18 states in one call), and Jd."""
    d = np.where(z > 0, 1e-6 * z, 1e-9)
    Z = step(np.hstack([z[:, None] + np.diag(d), z[:, None] - np.diag(d)]), m_in, **k)[0]
    Jd = P.T @ ((Z[:, :9] - Z[:, 9:]) / (2 * d)) @ P
    return np.linalg.eigvals(Jd), Jd


def fixed(m_in, k, u0):
    """Racemic fixed point from guess u0 = [T, m, B_homo, B_het, a], Newton in logs."""
    f = lambda x: step(np.exp(x)[SYM], m_in, **k)[0][[0, 2, 4, 6, 8]] / np.exp(x) - 1
    x = fsolve(f, np.log(u0), xtol=1e-13)
    if abs(f(x)).max() > 1e-9:                         # poor guess: relax by iterating the map, then polish
        z = np.exp(x)[SYM]
        for _ in range(3000):
            z = step(z, m_in, **k)[0]
        x = fsolve(f, np.log(z[[0, 2, 4, 6, 8]]), xtol=1e-13)
    assert abs(f(x)).max() < 1e-9, f"no fixed point at {m_in} mM, {k}"
    return np.exp(x)[SYM]


def scan_row(k, grid):
    """Fixed point and D-L block along the feed grid, continued up from the no-copying pools."""
    m = k["w"] * grid[0] / 2 / (k["w"] + rc.K4)
    a = 2 * m * np.sqrt(K3 * K1 / K2 / k["w"])
    b = K1 * m * m / (K3 + K2 * a + k["w"])
    out, u = [], np.array([rc.T_IN / 2, m, FO * b, FH * b, a])
    for m in grid:
        z = fixed(m, k, u)
        u = z[[0, 2, 4, 6, 8]]
        out.append((z, *mode(z, m, k)))
    return out


def lead(ev):
    return ev[np.argmax(abs(ev))]


def crossing(k, grid, rows):
    """First feed where the leading |eigenvalue| crosses 1 (bisection between grid points)."""
    gap = [abs(ev).max() - 1 for _, ev, _ in rows]
    for i in range(1, len(grid)):
        if gap[i - 1] < 0 < gap[i]:
            u = rows[i][0][[0, 2, 4, 6, 8]]
            return brentq(lambda m: abs(mode(fixed(m, k, u), m, k)[0]).max() - 1, grid[i - 1], grid[i], xtol=1e-3)
    return None


def r_fp(z, r, g):
    v = rates(z[2:, None], 0, r, g)[:, 0]
    return (v[2] + v[3]) / (v[0] + v[1])


def run(z, m_in, k, cycles, e):
    """Feed excess e (array: several runs at once); product excess (D-L)/(D+L) each cycle."""
    ee = []
    for _ in range(cycles):
        z, FD, FL = step(z, m_in, e=e, **k)
        ee.append((FD - FL) / (FD + FL))
    return z, np.array(ee)


def main():
    print("Arm G, bridged-dinucleotide copying (2-AI), fed reactor, 1 h cold phase per cycle, n = 3, w = 0.02")
    print(f"measured: k1 {K1} /h/mM, k2 {K2} /h/mM, k3 {K3} /h, k4 {rc.K4} /h, Ki {KI} mM (Walton & Szostak 2017)")
    print(f"          A*C kmax, Km by hand (Essex Fig. 4); hetero/homo formation {HET:.1e}/{HOM:.1e} -> fh {FH:.3f}, fo {FO:.3f}")
    print(f"          r0 = ag.r_bridged = {R0:.3f} (heteronucleotide bridges); C*C homonucleotide bridge: r = {R_CC:.3f}")
    print("assumed:  r scanned by kmax(+1 wrong) x r/r0; c scanned; copying speed-up g (1 = measured); no 2-AI fed")

    # self-check 1: bridging off = rna_copying.step / mode
    k4 = dict(n=3, r=0.5, c=0.5, w=0.02)
    zt = np.array([3e-4, 2e-4, 2.0, 3.0, 0, 0, 0, 0, 1.0])
    d_step = abs(step(zt, 20.0, e=0.1, bridged=False, **k4)[0][:4] / rc.step(zt[:4], 20.0, e=0.1, **k4)[0] - 1).max()
    z4 = rc.fixed_point(20.0, **k4)
    z9 = np.r_[z4, 0, 0, 0, 0, 1.0]
    d_fix = abs(step(z9, 20.0, bridged=False, **k4)[0][:4] / z4 - 1).max()
    ev9, Jd9 = mode(z9, 20.0, dict(k4, bridged=False))
    ev4, Jd4 = rc.mode(z4, 20.0, k4)
    d_J, d_lam = abs(Jd9[:2, :2] - Jd4).max() / abs(Jd4).max(), abs(abs(ev9).max() - abs(ev4).max())
    print(f"\nbridging off vs rna_copying.py (r = 0.5, c = 0.5, 20 mM): step rel. diff {d_step:.1e}, fixed point {d_fix:.1e},"
          f"\n    [[A, B], [C, Dd]] rel. diff {d_J:.1e}, leading |lam| {abs(ev9).max():.7f} vs {abs(ev4).max():.7f}")

    grid = np.geomspace(0.5, 100, 16)
    rs = sorted([0.06, 0.1, 0.2, R0, 0.6, R_CC, 0.9, 1.0, 1.1, 1.2])
    print("\nD-L block (4x4) at the racemic fixed point, feed 0.5-100 mM (16 points)")
    print("  lam = leading eigenvalue at the feed where |lam| is largest; C_m, C_tot = cross terms (docstring) over the feed")
    print("  range; r_fp at 100 mM; m* = feed where |lam| crosses 1, bridged / rna_copying.threshold (first order, g = 1)")
    print(f"{'g':>4} {'r':>6} {'c':>4} {'lam':>9} {'C_m min':>9} {'C_m max':>9} {'C_tot min':>9} {'C_tot max':>9} "
          f"{'r_fp':>6} {'m* bridged':>10} {'m* first':>9}")
    res, keep = {}, {}
    for g in (1, 10, 30):
        for r in rs:
            for c in (0.1, 0.3, 0.5, 1.0):
                k = dict(n=3, r=r, c=c, w=0.02, g=g)
                rows = scan_row(k, grid)
                cm = [Jd[1, 0] for _, _, Jd in rows]
                ct = [Jd[1, 0] + 2 * Jd[2, 0] for _, _, Jd in rows]
                lam = lead(max((ev for _, ev, _ in rows), key=lambda ev: abs(ev).max())).real
                ms = crossing(k, grid, rows)
                m1 = rc.threshold(dict(n=3, r=r, c=c, w=0.02))[0] if g == 1 else None
                rf = [r_fp(z, r, g) for z, _, _ in rows]
                res[g, r, c] = (lam, min(cm), max(cm), min(ct), max(ct), ms, min(rf), max(rf))
                keep[g, r, c] = rows
                print(f"{g:4d} {r:6.3f} {c:4.1f} {lam:9.5f} {min(cm):9.2e} {max(cm):9.2e} {min(ct):9.2e} {max(ct):9.2e} "
                      f"{rf[-1]:6.3f} {'none' if ms is None else f'{ms:.2f}':>10} "
                      f"{'' if g > 1 else 'none' if m1 is None else f'{m1:.2f}':>9}", flush=True)

    print("\nracemic operating point, r0, c = 0.5 (mM per hand or species): g, feed, free m, B_homo, B_het, 2-AI, T (uM),"
          " r_fp, lam")
    for g in (1, 30):
        for m, (z, ev, _) in list(zip(grid, keep[g, R0, 0.5]))[::3]:
            print(f"  {g:3d} {m:6.2f} {z[2]:8.3f} {z[4]:8.4f} {z[6]:8.4f} {z[8]:8.3f} {z[0] * 1e3:9.2f} {r_fp(z, R0, g):6.3f}"
                  f" {lead(ev).real:9.5f}")

    print("\n+20% seed (feed excess +0.2 / -0.2 / 0, EXPERIMENT.md arm G), 20 mM feed, c = 0.5, reactor filled with feed")
    cyc = [1, 10, 30, 100, 300, 1000, 5000]
    print(f"{'case':>26} {'feed':>7} " + " ".join(f"{x:>7}" for x in cyc))
    e, sym, plat = np.array([0.2, -0.2, 0.0]), True, {}
    for g, r in ((1, R0), (1, 1.2), (30, R0), (30, 1.2)):
        z0 = np.zeros((9, 3))
        z0[:2], z0[2], z0[3] = rc.T_IN / 2, 10 * (1 + e), 10 * (1 - e)
        ee = run(z0, 20.0, dict(n=3, c=0.5, w=0.02, r=r, g=g), cyc[-1], e)[1]
        for j, lab in enumerate(("+20% D", "+20% L", "0%")):
            print(f"{f'bridged g = {g}, r = {r:.3f}':>26} {lab:>7} " + " ".join(f"{ee[x - 1, j]:7.3f}" for x in cyc))
        if g == 1:
            ee1 = rc.run(np.array([rc.T_IN / 2] * 2 + [12.0, 8.0]), 20.0, dict(n=3, r=r, c=0.5, w=0.02), cyc[-1], e=0.2)[1]
            print(f"{f'first order r = {r:.3f}':>26} {'+20% D':>7} " + " ".join(f"{ee1[x - 1]:7.3f}" for x in cyc))
        sym &= np.allclose(ee[:, 0], -ee[:, 1], atol=1e-9) and np.all(abs(ee[:, 2]) < 1e-12)
        plat[g, r] = ee[-1, 0]
    print(f"    +20% D and +20% L equal and opposite at every cycle: {sym}")

    # self-check 2: RK4 converged at the stiffest point (fast copying, depleted bridges)
    k = dict(n=3, r=R0, c=1.0, w=0.02, g=30)
    z = keep[30, R0, 1.0][-1][0]
    global NSUB
    a = step(z, 100.0, **k)[0]
    NSUB *= 4
    d_rk = abs(step(z, 100.0, **k)[0] / a - 1).max()
    NSUB //= 4
    print(f"RK4 {3 * NSUB} vs {12 * NSUB} steps per cold phase at g = 30, r0, c = 1, 100 mM: rel. diff {d_rk:.1e}")

    lost = {key: v[5] is not None or abs(v[0]) > 1 for key, v in res.items()}
    print("\nANSWER")
    for g in (1, 10, 30):
        sub = {key: v for key, v in res.items() if key[0] == g}
        low = [key for key in sub if key[1] < 1 and lost[key]]
        print(f"  g = {g:3d}: C_m > 0 at every r < 1: {all(v[1] > 0 for key, v in sub.items() if key[1] < 1)};"
              f"  C_tot < 0 at every r < 1: {all(v[4] < 0 for key, v in sub.items() if key[1] < 1)};"
              f"  r_fp {min(v[6] for v in sub.values()):.3f}-{max(v[7] for v in sub.values()):.3f}")
        print(f"           50/50 lost at r < 1: " + (", ".join(f"(r {r:.2f}, c {c:.1f}, m* {sub[g, r, c][5]}, r_fp "
              f"{sub[g, r, c][7]:.2f})" for _, r, c in low) or "nowhere") +
              f";  at r >= 1: {sum(lost[key] for key in sub if key[1] >= 1)} of {sum(key[1] >= 1 for key in sub)} (r, c)")
    print(f"  measured rates (g = 1, r0 = {R0:.2f} and {R_CC:.2f}): max |lam| = "
          f"{max(abs(res[key][0]) for key in res if key[0] == 1 and key[1] in (R0, R_CC)):.5f}")
    print("  +20% seed plateau (5000 cycles): " + ", ".join(f"g {g}, r {r:.2f}: {v:.3f}" for (g, r), v in plat.items()))

    assert d_step < 1e-3 and d_fix < 1e-3 and d_J < 1e-3 and d_lam < 1e-5 and sym and d_rk < 1e-5, "self-check failed"
    print("\nself-check passed")


if __name__ == "__main__":
    main()
