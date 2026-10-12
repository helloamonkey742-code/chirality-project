"""Objection 2 (non-uniform sweep / first crossing): if the transition (lam > 0) is crossed first in a small region
(near a vent, under a thinning shell) and that region then invades the rest, is the effective N the first-crossing
volume rather than the whole ocean?

1D line of L cells, no-flux ends (x = 0 is a shore/ice wall where the transition arrives first), same local rule as
domains.py with a transition FRONT:
    d(alpha_i) = (lam_i*alpha_i - alpha_i**3 + g + D*lap(alpha)_i) dt + sqrt(eps) dW_i,
    lam(x, t) = clip(gamma * (t - x / v_front), -1, 1)
Each cell sees the same sweep as sim.py (-1 -> +1 at rate gamma); only the start time varies along x.
gamma = 0.01, eps = 1e-4 per cell, D = 16: coherence length xi = sqrt(D) gamma^-1/4 = 12.6 cells, front-speed scale
V0 = sqrt(D gamma^1/2) = 1.26 cells per time. u = v_front / V0 = 0.03 .. 30 and infinity (uniform sweep).
L = 4, 16, 64 xi. One bias g for every run, chosen so Delta_total = 4 at L = 16 xi (2 at 4 xi, 8 at 64 xi).
Outcome: majority hand at the end; also whether the line is a single hand, and whether the final majority equals
the hand at x = 0 when lam(0) first reached 0.4 ("inheritance").

Predictions compared:  (i)  Phi(Delta_total), N = all L cells;
 (ii)  Phi(Delta_first), N = cells crossed within one selection window, l_first = min(L, v_front * gamma^-1/2);
 (iii) Phi(Delta_Airy), N = the first unstable mode. Theory (written before the runs): near the wall the linear
 operator D d2/dx2 + lam(x, t) with lam' = -gamma/v has a Neumann Airy ground state psi = Ai(x/w + a1'),
 a1' = -1.0188, of width w = (D v / gamma)^(1/3). Projecting the noise and the bias on psi gives an effective length
 l_A = (int psi)^2 / int psi^2 = 2.239 w. The ground state is separated from the next mode by 2.23 D/w^2 in lam,
 i.e. 2.23 u^(-2/3) selection windows, so it alone sets the sign when u < ~1. In general the effective length
 l_eff = eps * (Phi^-1(P) / Delta(g, eps=1))^2 is read off the simulated P.
 Invasion: the selected hand's edge lags the transition front by d = v^3 / (4 D gamma) (Fisher speed 2 sqrt(D lam)
 = v). Water in the lag amplifies its own noise by exp(u^4/32) before it is invaded; independent selection starts
 when that reaches 1 / sigma_0, sigma_0 ~ sqrt(eps gamma^-1/2 / xi) ~ 0.009, so u* = (32 ln(1/sigma_0))^(1/4) ~ 3.5.

PREDICTIONS (written 2026-10-12, before any result-producing run):
 P1. Uniform limit (v = infinity, L = 3 cells): P = Phi(Delta_total) within 3 SE (asserted in check()).
 P2. Fast front (u = 30; L / v <= 0.53 gamma^-1/2 at L = 16 xi): P equals the v = infinity value within 2 combined SE
     at L = 4 and 16 xi.
 P3. Slow front (u <= 0.3): the final majority equals the hand at x = 0 in >= 95% of runs and >= 90% of runs end
     as a single hand, at every L.
 P4. Slow front (u <= 0.3): P does not depend on L (4, 16, 64 xi agree within 2.5 combined SE) and lies below
     Phi(Delta_total) at L = 16 xi by more than 5 SE. The whole-volume rule fails here.
 P5. Slow front (u <= 0.3): l_eff / l_A lies in [0.67, 1.5]; the task's first-crossing length v gamma^-1/2
     UNDERestimates l_eff by more than 2x at u <= 0.1.
 P6. Crossover: at L = 64 xi the single-hand fraction falls through 50% at u* in [1, 10] (point estimate 3.5).
     Below u*, the effective N is a front-set length that does not grow with L; above it, the 1D patch result of
     kz_patches.py applies with the whole L.
 P7. Convergence: halving dt and halving dx at L = 16 xi, u = 1 change P by less than 2 combined SE.

Usage: python gradient_ramp.py [L4] [L16] [L64a] [L64b] [conv]   (default: all; check always runs last)
"""
import sys
import time

import numpy as np
from scipy.stats import norm

from kz_patches import pse
from sim import delta

GAMMA, EPS, D, DT = 0.01, 1e-4, 16.0, 0.02
XI = np.sqrt(D) * GAMMA ** -0.25             # 12.65 cells
V0 = np.sqrt(D) * GAMMA ** 0.25              # 1.265 cells per unit time
KAPPA_A = 2.2393                             # (int Ai)^2 / int Ai^2 on [a1', inf), a1' = first zero of Ai'
G = 4.0 / delta(1.0, EPS / round(16 * XI), GAMMA)      # Delta_total = 4 at L = 16 xi
US = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, np.inf)


def lap_n(a):
    p = np.concatenate([a[:, :1], a, a[:, -1:]], axis=1)    # no-flux ends
    return p[:, 2:] + p[:, :-2] - 2 * a


def ramp(n, u, runs, seed, g=G, dt=DT, h=1.0):
    """n cells of size h. Returns (score, single hand, inherited) per run."""
    Dl = D / h ** 2
    assert Dl * dt <= 0.5 + 1e-12, "Euler stability D dt/dx^2 <= 1/2"
    rng = np.random.default_rng(seed)
    x = (np.arange(n) + 0.5) * h
    tc = x / (u * V0)                                         # u = inf -> 0
    t0, t1 = -1 / GAMMA, tc.max() + 1 / GAMMA
    a = np.zeros((runs, n), np.float32)
    noise, dtf, Dl, gf = np.float32(np.sqrt(EPS / h * dt)), np.float32(dt), np.float32(Dl), np.float32(g)
    early = None
    for i in range(int(np.ceil((t1 - t0) / dt))):
        lam = np.clip(GAMMA * (t0 + i * dt - tc), -1, 1).astype(np.float32)
        a += (lam * a - a ** 3 + gf + Dl * lap_n(a)) * dtf + noise * rng.standard_normal(a.shape, dtype=np.float32)
        if early is None and lam[0] >= 0.4:
            early = np.sign(a[:, 0])
    s = np.sign(a).mean(1)
    return (s > 0) + 0.5 * (s == 0), np.abs(s) == 1, np.sign(s) == early


def lengths(L, u):
    v = u * V0
    return min(L, v * GAMMA ** -0.5), min(L, KAPPA_A * (D * v / GAMMA) ** (1 / 3))


def l_eff(p, se):
    k = delta(G, EPS, GAMMA)                                  # Delta per sqrt(cell)
    if not 0.5 < p < 1:
        return np.nan, np.nan
    z = norm.ppf(p)
    return (z / k) ** 2, 2 * z / k ** 2 / norm.pdf(z) * se


def point(Lx, u, runs, seed, **kw):
    t0 = time.time()
    h = kw.get("h", 1.0)
    L = round(Lx * XI)
    sc, one, inh = ramp(int(round(L / h)), u, runs, seed, **kw)
    p, se = pse(sc)
    lf, la = lengths(L, u)
    le, lse = l_eff(p, se)
    P = lambda ell: norm.cdf(delta(G, EPS / ell, GAMMA))
    print(f"  L={Lx:<3g}xi ({L:4d}) u={u:<5g} n={runs}  P={p:.3f}+-{se:.3f} | Phi(D_tot)={P(L):.3f} "
          f"Phi(D_first)={P(lf):.3f} Phi(D_Airy)={P(la):.3f} | l_eff={le:6.1f}+-{lse:5.1f} l_first={lf:6.1f} "
          f"l_A={la:5.1f} | 1hand {one.mean():.3f} inherit {inh.mean():.3f}  [{time.time() - t0:.0f}s]", flush=True)
    return p, se, one.mean(), inh.mean()


def sweep_L(Lx, us, runs, seed0):
    print(f"\nL = {Lx} xi = {round(Lx * XI)} cells, Delta_total = {delta(G, EPS / round(Lx * XI), GAMMA):.2f}", flush=True)
    return {u: point(Lx, u, runs, seed0 + j) for j, u in enumerate(us)}


def check():
    """Uniform fast-sweep limit: 3 cells, v = inf, must give sim.py's Phi(Delta)."""
    runs = 4000
    g = 1.0 / delta(1.0, EPS / 3, GAMMA)
    sc, _, _ = ramp(3, np.inf, runs, seed=7, g=g)
    p, se = pse(sc)
    want = norm.cdf(delta(g, EPS / 3, GAMMA))
    print(f"\ncheck: uniform limit (3 cells, v=inf) P={p:.4f}+-{se:.4f}, sim.py Phi(Delta)={want:.4f}, "
          f"diff {(p - want) / se:+.2f} SE")
    assert abs(p - want) < 3 * se, "uniform limit does not reproduce sim.py Phi(Delta)"


def main():
    secs = sys.argv[1:] or ["L4", "L16", "L64a", "L64b", "conv"]
    t0 = time.time()
    print(f"gamma={GAMMA} eps={EPS}/cell D={D} dt={DT}: xi={XI:.2f} cells, V0={V0:.3f} cells/time, g={G:.3e}")
    if "L4" in secs:
        sweep_L(4, US, 600, 400)
    if "L16" in secs:
        sweep_L(16, US, 600, 500)
    if "L64a" in secs:
        sweep_L(64, (0.3,), 400, 600)
    if "L64b" in secs:
        sweep_L(64, (1.0, 3.0, 10.0, np.inf), 400, 610)
    if "conv" in secs:
        print("\nconvergence at L=16xi, u=1 (compare with the L16 row u=1)")
        point(16, 1.0, 600, 700, dt=DT / 2)
        point(16, 1.0, 600, 701, dt=DT / 4, h=0.5)
    check()
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
