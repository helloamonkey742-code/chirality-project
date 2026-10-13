"""front2d: does gradient_ramp.py's 1D front result survive in 2D and 3D, and do real oceans cross the transition as a
front? (H21, H22 of journals/ija_draft.md; 2026-10-12.)

PART A - the front in 2D and 3D (H21). Same local rule as gradient_ramp.py / domains.py:
    d(alpha) = (lam alpha - alpha^3 + g + D lap(alpha)) dt + sqrt(eps / h^dim) dW,  lam(x, t) = clip(gamma (t - x / v), -1, 1)
on an nx x W (2D) or nx x W x W (3D) grid: no-flux in x (x = 0 is the wall the transition reaches first), periodic
across (a planar front of width W). gamma = 0.01, eps = 1e-4 per cell, D = 4 (cheaper than gradient_ramp's D = 16):
xi = sqrt(D) gamma^-1/4 = 6.32 cells, V0 = sqrt(D) gamma^1/4 = 0.632 cells per time, u = v / V0 (front advance per
selection window gamma^-1/2, in coherence lengths). Airy layer l_A = 2.2393 (D v / gamma)^(1/3) = 2.24 u^(1/3) xi.
L = 8 xi = 51 cells unless stated. For each width the bias g is set so Delta_total = 2 at L = 8 xi (whole volume).
Outcome: majority hand at the end (score 1 / 0.5 / 0) and whether the whole grid is one hand.
Compared: Phi(Delta_total) (whole volume) and Phi(Delta_A) (volume l_A x W^(dim-1)). l_eff per unit width is read
off P: l_eff = (Phi^-1(P) / Delta(g, eps = 1))^2 / W^(dim-1), the same inversion as gradient_ramp.l_eff.
Why the 2D answer is not obvious: along the front the transverse modes of the Airy layer are almost degenerate
(gap D (2 pi / W)^2), so the layer is itself a 1D transition along y and can break into patches of size ~xi. If
patches decide independently, the majority still follows the summed bias (kz_patches.py: c = 0.98 in 2D), so the
expectation is N = l_A x W, with the one-hand fraction falling as W grows.

PART B - per-world drive uniformity (H22). Physical assumptions, each labelled:
 F1 (drive): lam crosses 0 because the reacting concentration builds up (k = k2 c rises at the sweep rate 1/tau). It
    is supplied locally (vents) and spread by horizontal mixing D_h. The spread in crossing time across the ocean is
    the mixing time t_mix = (L/A)^2 / D_h, A = 5.5 (sweep2's healing-across time), in both regimes: if t_mix > tau the
    drive is a front that travels at the mixing speed; if t_mix < tau the residual heterogeneity is
    d(lam) ~ t_mix / tau, i.e. again a crossing spread of ~t_mix. So the front speed is v = L / t_mix = A^2 D_h / L.
 F2 (hand mixing): the hand field mixes with the same D_h, so in physical units xi = sqrt(D_h/k) (k tau)^(1/4),
    V0 = sqrt(D_h k) (k tau)^(-1/4), and u = v / V0 = A^2 xi / L.
 F3 (front layer): a front spans the depth of the deciding layer and a length L across the ocean (one vent region).
    Its deciding volume is l_A L H f_old, l_A = 2.24 (D_h v tau / k)^(1/3), so f_new = f_old min(1, pi l_A / (4 L))
    (ocean area 4 L^2 / pi). Equivalent: f_new / f_old = 0.0582 u^(4/3).
 F4 (rule): u < u* -> N from f_new; u >= u* -> whole-volume N kept (as multiwall). u* = 3 (gradient_ramp's
    (1, 3]; Part A tests it in 2D). Sensitivities: u* = 1; v x A and v / A (the order-one choices in F1).
 F5: arrival, healing, wall and recursion conditions are unchanged (the front-set hand invades behind the drive front,
    which crosses the ocean in t_mix, the same as sweep2's "across" time).
 F6: a globally set drive (cooling through the shell) is uniform: front fraction 0, shares = multiwall main.
Implementation: multiwall.run4's main loop is repeated with the same draws (sweep2.chains, default_rng([seed, i])),
and ocean.delta_at is wrapped in the sweep2 / sweep3 / multiwall namespaces so every N gets the per-draw factor
f_new / f_old. Nothing is edited. With the front rule off (factor 1) the shares must equal multiwall main.

PREDICTIONS (written 2026-10-12, before any result-producing run):
 P1. Uniform limit (2 x 2 cells, v = infinity): P = Phi(Delta_total) within 3 SE (asserted in check()).
 P2. 1D limit (W = 1 cell) at this D: for u <= 1, l_eff / l_A lies in [0.67, 1.5] (gradient_ramp's P5 at D = 4).
 P3. 2D slow front (u <= 0.3, W = 2 xi and 8 xi): l_eff per unit width / l_A lies in [0.5, 1.5]. The front layer
     acts as one N across the whole width: N = l_A x W.
 P4. 2D L-independence (W = 2 xi, u = 0.3): P at L = 4, 8, 16 xi agrees within 2.5 combined SE, and at 16 xi lies
     below Phi(Delta_total) by more than 5 SE.
 P5. 2D crossover: u*_2D lies in (1, 3], as in 1D. At W = 2 xi and 8 xi, l_eff / l_A <= 1.5 at u = 1 and > 1.5 at
     u = 3 (L = 8 xi); and at W = 2 xi, u = 3, P(L = 16 xi) - P(L = 4 xi) > 2 combined SE (whole-volume growth
     returns above u*).
 P6. 2D fast limit (u = infinity, W = 2 xi and 8 xi): l_eff / L in [0.64, 1.3] (P = Phi(c Delta_total), c in
     [0.8, 1.14], as kz_patches' 2D c = 0.98).
 P7. Patches along the front: at W = 8 xi and u <= 0.3 the one-hand fraction is lower than at W = 1 cell by >= 0.1.
 P8. Convergence (W = 2 xi, L = 8 xi, u = 0.3): dt/2 and dx/2 each change P by less than 2 combined SE.
 P9. 3D (W = 2 xi x 2 xi): l_eff per unit cross-section / l_A in [0.5, 1.5] at u = 0.3 and 1, > 1.5 at u = 3;
     l_eff / L in [0.64, 1.3] at u = infinity.
 P10. Regression: front rule off reproduces multiwall.run4(seed, full=False) main shares exactly (equal floats) for
     all 15 worlds and seed 0, and multiwall.log Section 4 printed to 0.1.
 P11. Front fraction (u < 3, all draws, F1-F2) is <= 30 % in every world and <= 15 % for Early Earth; it orders
     Enceladus < Europa < Ganymede (smaller oceans mix across faster).
 P12. Among draws that reach one hand (no lid), the front fraction is >= the all-draw fraction in every world.
 P13. Front rule on (u* = 3): every world's P_pved (no lid and final) falls by 0-5 points; P_none is unchanged
     exactly (the rule only moves draws from P_pved to P_chance); Early Earth stays >= 99 %; the final ordering
     Ganymede > Europa > Enceladus is unchanged. The headlines do not change.
 P14. Sensitivity: v / A (5.5x slower front) lowers P_pved by <= 10 points in every world.
 P15. Seed 1: front-rule shares within 0.5 points of seed 0 in every world.

Usage: python front2d.py [W1] [W2] [W8a] [W8b] [L] [conv] [3Da] [3Db] [worlds]   (check() always runs last)
"""
import sys
import time

import numpy as np
from scipy.stats import norm

from gradient_ramp import EPS, GAMMA, KAPPA_A, lap_n
from kz_patches import pse
from sim import delta

D, DT = 4.0, 0.05
XI = np.sqrt(D) * GAMMA ** -0.25             # 6.32 cells
V0 = np.sqrt(D) * GAMMA ** 0.25              # 0.632 cells per unit time
L8 = round(8 * XI)                           # 51 cells
US = (0.1, 0.3, 1.0, 3.0, 10.0, np.inf)


def lap(a):
    out = lap_n(a)                                           # no-flux along x (axis 1)
    for ax in range(2, a.ndim):                              # periodic across the front
        out += np.roll(a, 1, ax) + np.roll(a, -1, ax) - 2 * a
    return out


def front(nx, w, dim, u, runs, seed, g, dt=DT, h=1.0):
    """nx x w^(dim-1) cells of size h. Returns (score, one hand) per run."""
    Dl = D / h ** 2
    assert Dl * dt <= 1 / (2 * dim) + 1e-12, "Euler stability D dt/dx^2 <= 1/(2 dim)"
    rng = np.random.default_rng(seed)
    tc = ((np.arange(nx) + 0.5) * h / (u * V0)).reshape((nx,) + (1,) * (dim - 1))
    t0, t1 = -1 / GAMMA, tc.max() + 1 / GAMMA
    a = np.zeros((runs, nx) + (w,) * (dim - 1), np.float32)
    noise, dtf, Dl, gf = np.float32(np.sqrt(EPS / h ** dim * dt)), np.float32(dt), np.float32(Dl), np.float32(g)
    for i in range(int(np.ceil((t1 - t0) / dt))):
        lam = np.clip(GAMMA * (t0 + i * dt - tc), -1, 1).astype(np.float32)
        a += (lam * a - a ** 3 + gf + Dl * lap(a)) * dtf + noise * rng.standard_normal(a.shape, dtype=np.float32)
    s = np.sign(a).reshape(runs, -1).mean(1)
    return (s > 0) + 0.5 * (s == 0), np.abs(s) == 1


def g_for(w, dim):
    return 2.0 / delta(1.0, EPS / (L8 * w ** (dim - 1)), GAMMA)      # Delta_total = 2 at L = 8 xi


def point(Lx, w, dim, u, runs, seed, dt=DT, h=1.0):
    t0 = time.time()
    L, g, area = round(Lx * XI), g_for(w, dim), w ** (dim - 1)
    sc, one = front(int(round(L / h)), int(round(w / h)), dim, u, runs, seed, g, dt, h)
    p, se = pse(sc)
    lA = min(L, KAPPA_A * (D * u * V0 / GAMMA) ** (1 / 3))
    k = delta(g, EPS, GAMMA)
    z = norm.ppf(p) if 0.5 < p < 1 else np.nan
    le, lse = (z / k) ** 2 / area, 2 * z / k ** 2 / norm.pdf(z) * se / area
    P = lambda vol: norm.cdf(delta(g, EPS / vol, GAMMA))
    print(f"  {dim}D L={Lx:<3g}xi ({L:3d}) W={w:3d} u={u:<4g} n={runs} dt={dt:g} h={h:g}  P={p:.3f}+-{se:.3f} | "
          f"Phi(D_tot)={P(L * area):.3f} Phi(D_A)={P(lA * area):.3f} | l_eff={le:6.1f}+-{lse:5.1f} l_A={lA:5.1f} "
          f"ratio={le / lA:5.2f}+-{lse / lA:4.2f} l_eff/L={le / L:4.2f} | 1hand {one.mean():.3f}  [{time.time() - t0:.0f}s]",
          flush=True)
    return p, se


def check():
    """Uniform fast-sweep limit: 2 x 2 cells, v = inf, must give sim.py's Phi(Delta)."""
    runs, g = 4000, 1.0 / delta(1.0, EPS / 4, GAMMA)
    sc, _ = front(2, 2, 2, np.inf, runs, seed=7, g=g)
    p, se = pse(sc)
    want = norm.cdf(delta(g, EPS / 4, GAMMA))
    print(f"\ncheck: uniform limit (2x2 cells, v=inf) P={p:.4f}+-{se:.4f}, sim.py Phi(Delta)={want:.4f}, "
          f"diff {(p - want) / se:+.2f} SE")
    assert abs(p - want) < 3 * se, "uniform limit does not reproduce sim.py Phi(Delta)"


# ---------------------------------------------------------------- Part B: per-world front rule (wrapper, no edits)
def worlds():
    import multiwall as M
    import ocean
    import sweep2 as S
    import sweep3 as S3
    import uncertainty as U

    PHI = [1.0]
    wrapped = lambda c, g, k2, V, tau, f=1.0: ocean.delta_at(c, g, k2, V, tau, f * PHI[0])
    for mod in (S, S3, M):
        mod.delta_at = wrapped

    def factor(x, L, ustar, vmul):
        c, k2, Dh, t = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time"))
        k = k2 * c
        v = vmul * U.A ** 2 * Dh / L                                     # F1
        u = v / (np.sqrt(Dh * k) * (k * t) ** -0.25)                     # F2
        lA = KAPPA_A * (Dh * v * t / k) ** (1 / 3)                       # F3
        return np.where(u < ustar, np.minimum(1.0, np.pi * lA / (4 * L)), 1.0), u

    def shares(seed, W, ustar=3.0, vmul=1.0, on=True):
        out = {}
        for i, (name, (x, _, a, salt)) in enumerate(W.items()):
            rng = np.random.default_rng([seed, i])                       # multiwall.run4's extra draws, same order
            ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N),
                  "rho": 10 ** rng.uniform(0, 2, U.N)}
            none = np.zeros(U.N, bool)
            fac, u = factor(x, a[1], ustar, vmul)
            PHI[0] = fac if on else 1.0
            res = [M.multi(x, a[1], a[4], ex, l)[:2] for l in (none, salt)]
            PHI[0] = 1.0
            out[name] = {"main": [S3.split(*r) for r in res], "u": u, "reached": res[0][1], "front": u < ustar}
        return out

    t0 = time.time()
    print(f"\n=== Part B: per-world front rule (multiwall Section 4 with N -> l_A layer where u < u*), N = {U.N} ===")
    _, W0 = S.chains(0)
    off = shares(0, W0, on=False)
    ref = M.run4(0, full=False)
    exact = all(off[n]["main"] == ref[n]["main"] for n in ref)
    log = {}
    for line in open("multiwall.log"):
        if "||" in line and line[:18].strip() in ref:
            log[line[:18].strip()] = line.split("||")[1].replace("|", "/").split("/")
    printed = all([f"{v:.1f}" for s in off[n]["main"] for v in s] == [t.strip() for t in log[n]] for n in ref)
    print(f"regression, front rule off vs multiwall.run4(0) main: {'EXACT' if exact else 'FAILED'}; "
          f"vs multiwall.log Section 4 (printed 0.1, 15 worlds): {'EXACT' if printed and len(log) == 15 else 'FAILED'}")
    assert exact and printed and len(log) == 15, "front rule off does not reproduce multiwall main"

    main = shares(0, W0)
    var = {"u*=1": shares(0, W0, ustar=1.0), "v x A": shares(0, W0, vmul=U.A), "v / A": shares(0, W0, vmul=1 / U.A)}
    for n in main:
        assert all(main[n]["main"][l][2] == off[n]["main"][l][2] for l in (0, 1)), f"{n}: P_none moved"
    order = sorted(main, key=lambda n: -off[n]["main"][1][0])
    print("\nDrive uniformity (F1-F2, u = A^2 xi / L): % of draws with u < 3 (all | reaching one hand, no lid) | "
          "u quartiles (all draws) | % with u < 1 | t_mix > tau %")
    for n in order:
        x, a, q = W0[n][0], W0[n][2], main[n]
        tmix = (a[1] / U.A) ** 2 / x["horizontal D_h"] > x["time"]
        r = q["reached"]
        print(f"{n:18} {100 * q['front'].mean():5.1f} | {100 * q['front'][r].mean():5.1f} | "
              + " ".join(f"{v:8.3g}" for v in np.percentile(q["u"], (25, 50, 75)))
              + f" | {100 * (q['u'] < 1).mean():5.1f} | {100 * tmix.mean():5.1f}")
    print("\nShares, %: P_pved / P_chance / P_none, no lid | final")
    print(f"{'world':18} {'front rule off (= multiwall main)':>38} || {'front rule on, u* = 3':>36}")
    for n in order:
        print(f"{n:18} " + " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in off[n]["main"]) + "  ||  "
              + " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in main[n]["main"]))
    print("\nP_pved change vs multiwall main, points, no lid / final: main (u*=3) | " + " | ".join(var))
    for n in order:
        d = lambda r: " / ".join(f"{r[n]['main'][l][0] - off[n]['main'][l][0]:+5.2f}" for l in (0, 1))
        print(f"{n:18} {d(main)} | " + " | ".join(d(var[k]) for k in var))
    _, W1 = S.chains(1)
    s1 = shares(1, W1)
    dmax = max(abs(s1[n]["main"][l][j] - main[n]["main"][l][j]) for n in main for l in (0, 1) for j in range(3))
    print(f"\nseed 1, front rule on: max |seed 1 - seed 0| over worlds, splits, lid = {dmax:.2f} points  "
          f"[{time.time() - t0:.0f}s]")


def main():
    secs = sys.argv[1:] or ["W1", "W2", "W8a", "W8b", "L", "conv", "3Da", "3Db", "worlds"]
    t0 = time.time()
    print(f"gamma={GAMMA} eps={EPS}/cell D={D} dt={DT}: xi={XI:.2f} cells, V0={V0:.3f} cells/time, L=8xi={L8} cells; "
          f"sections {secs}", flush=True)
    w2, w8 = round(2 * XI), round(8 * XI)
    if "W1" in secs:
        print("\n1D limit (W = 1 cell)")
        for j, u in enumerate(US):
            point(8, 1, 2, u, 400, 100 + j)
    if "W2" in secs:
        print(f"\n2D, W = 2 xi = {w2} cells")
        for j, u in enumerate(US):
            point(8, w2, 2, u, 400, 200 + j)
    if "W8a" in secs:
        print(f"\n2D, W = 8 xi = {w8} cells (a)")
        point(8, w8, 2, 0.1, 400, 300)
    if "W8b" in secs:
        print(f"\n2D, W = 8 xi = {w8} cells (b)")
        for j, u in enumerate(US[1:]):
            point(8, w8, 2, u, 400, 301 + j)
    if "L" in secs:
        print(f"\n2D size scan, W = 2 xi, g fixed (Delta_total = 1.41 / 2.83 at 4 / 16 xi)")
        for j, (Lx, u) in enumerate(((4, 0.3), (16, 0.3), (4, 3.0), (16, 3.0))):
            point(Lx, w2, 2, u, 400, 400 + j)
    if "conv" in secs:
        print("\nconvergence at W = 2 xi, L = 8 xi, u = 0.3 (compare with the W2 row u = 0.3)")
        point(8, w2, 2, 0.3, 400, 500, dt=DT / 2)
        point(8, w2, 2, 0.3, 400, 501, dt=DT / 4, h=0.5)
    if "3Da" in secs:
        print(f"\n3D, W = 2 xi x 2 xi = {w2} x {w2} cells (a)")
        for j, u in enumerate((0.3, np.inf)):
            point(8, w2, 3, u, 250, 600 + j, dt=0.04)
    if "3Db" in secs:
        print(f"\n3D, W = 2 xi x 2 xi = {w2} x {w2} cells (b)")
        for j, u in enumerate((1.0, 3.0)):
            point(8, w2, 3, u, 250, 610 + j, dt=0.04)
    if "worlds" in secs:
        worlds()
    check()
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
