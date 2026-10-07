"""Polar gap: does a hole in the stratified lid (polar convection, as in Earth's deep-water formation) free a pinned wall?

two_zone.py (1D, radius only) pinned the wall under a stratified layer that covers the whole sphere. Real lids can have
holes: Earth's thermocline outcrops at the poles; Zeng & Jansen 2024 find polar FREEZING on Enceladus (brine rejection
-> convection) unless the ice shell is tidally heated, while Ames et al. 2025 find polar MELTING (stratified) at steady state.
Same rule as two_zone.py on an axisymmetric hemisphere (radius r, colatitude th; both poles have the gap, no flux at
the equator by symmetry):
    da/dt = a - a^3 + h + (1/r^2) d_r(r^2 D_r d_r a) + (1/(r^2 sin th)) d_th(sin th D_th d_th a)
Stratified lid: top 0.4 of the ocean, at colatitude th > th_gap. It cuts VERTICAL mixing only (D_r = D_hi/ratio);
sideways mixing stays D_hi (stratification damps vertical motion, not horizontal).

PREDICTIONS (written before any run, 2026-10-06):
 (P1) No lid: + takes the whole ocean. Lid with no gap: wall pinned (two_zone.py).
 (P2) With a gap, + sinks through the gap to the floor, leaving a deep + column at the pole and a - pocket under the lid.
      The pocket's side wall (height h ~ 0.65, tension ~ sqrt(D_hi)) can move toward the equator only if that shrinks the
      total wall energy: the lid wall above the pocket (tension ~ sqrt(D_hi/ratio)) loses area 2 pi R^2 sin th dth while
      the side wall gains 2 pi R h cos th dth. The pocket shrinks only if
          tan(th_gap) > (h / R) * sqrt(D_side / D_lid)                       ("critical gap" th_c)
      R = 4, ratio 10: th_c ~ 27 deg -> gap 15 deg stuck, 35 deg + wins.
      R = 4, ratio 100: th_c ~ 58 deg -> gap 45 deg stuck, 70 deg + wins.
      (Ceiling: energy argument for a static wall, ignores the junction shape; expect th_c within ~10 deg.)
 (P3) Bias: a hand-favouring bias h adds a pressure 2h that can push the wall through the lid. Needed: 2h ~ (jump in
      tension)/(wall width) -> h_c ~ 0.05-0.3 in model units (two_zone.py geometry, ratio 10, lid on top).
      The real parity-violating bias is g ~ 1e-19..1e-16 (h ~ g in these units), so it cannot unpin anything.

RESULTS (polar_gap.log, 2026-10-06): P1, P2, P3 held.
 Controls: no lid takes over at t=105 (dt/2: 102); lid without gap pinned 0.09 above the lid's lower edge (two_zone: 0.08).
 Critical gap: ratio 10 between 25 and 30 deg (predicted 27); ratio 100 between 55 and 70 deg (predicted 58).
 Below it the + water in the gap does not even reach the floor (narrow gap = costly side walls). Above it + wins everywhere.
 Bias to push through a full lid: h_c = 0.200-0.206 (predicted 0.05-0.3). Real g <= 4e-16: ~15 orders too small.
 Ceiling: "STUCK" = + share moved < 2e-4 in 500 time units; a very slow creep would be misread as stuck.

What P2 means for real oceans (real units, sigma ~ sqrt(D) across each wall): tan(th_c) = (h/R) sqrt(D_h,deep / D_z,lid).
"""
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix, identity
from scipy.sparse.linalg import splu

import two_zone

R, H, D_HI, LAYER, DT = 4.0, 1.0, 1e-2, 0.4, 0.1
NR, NTH = 50, 314                                   # dr = 0.02, R*dth = 0.02 (wall width sqrt(2 D_HI) = 0.14)


def run(ratio, gap_deg, lid=True, dt=DT, t_max=4000.0, W=0.0, d_hi=D_HI, grid=1.0, top=0.0):
    """gap_deg: polar gap half-width in degrees of colatitude (0 = no gap). W: overturning current (currents.py),
    psi = -W R^2 sin(pi (r - r_in)/H) sin^2(th) cos(th): sinks at the pole at speed ~2W, rises at low latitudes.
    d_hi, grid: thinner walls (realistic.py) with grid refined by `grid`; top: fast mixed layer above the lid. Returns (outcome, t, plus share, pole depth, side angle)."""
    r_in = R - H
    nr, nth = int(round(NR * grid)), int(round(NTH * grid))
    slow_t = D_HI / d_hi                                          # stuck test window scales with the diffusion time
    dr, dth = H / nr, (np.pi / 2) / nth
    rc = r_in + (np.arange(nr) + 0.5) * dr
    tc = (np.arange(nth) + 0.5) * dth
    rf, tf = r_in + np.arange(1, nr) * dr, np.arange(1, nth) * dth
    w = 0.02
    s = lambda x, wid: 0.5 * (1 + np.tanh(x / wid))
    gap = np.radians(gap_deg)

    def slow(r, th):                                              # 1 inside the lid, 0 outside
        if not lid:
            return np.zeros(np.broadcast(r, th).shape)
        on = s(th - gap, w / R) if gap > 0 else np.ones_like(th)
        return s(r - (R - LAYER), w) * (s((R - top) - r, w) if top else 1) * on

    lnDr = np.log(d_hi) - slow(rf[:, None], tc[None, :]) * np.log(ratio)  # radial faces (NR-1, NTH)
    Dr = np.exp(lnDr)
    I, J = np.meshgrid(np.arange(nr), np.arange(nth), indexing="ij")
    k = lambda i, j: i * nth + j
    rows, cols, vals = [], [], []

    def couple(ka, kb, ca, cb):                                   # flux between cells a and b, coefficients seen from a, b
        rows.extend([ka, ka, kb, kb]); cols.extend([kb, ka, ka, kb]); vals.extend([ca, -ca, cb, -cb])

    i, j = I[:-1], J[:-1]                                         # radial neighbours
    couple(k(i, j).ravel(), k(i + 1, j).ravel(),
           (Dr * rf[:, None] ** 2 / (rc[:-1, None] ** 2 * dr ** 2)).ravel(),
           (Dr * rf[:, None] ** 2 / (rc[1:, None] ** 2 * dr ** 2)).ravel())
    i, j = I[:, :-1], J[:, :-1]                                   # colatitude neighbours, D_th = d_hi everywhere
    couple(k(i, j).ravel(), k(i, j + 1).ravel(),
           (d_hi * np.sin(tf)[None, :] / (rc[:, None] ** 2 * np.sin(tc[:-1])[None, :] * dth ** 2)).ravel(),
           (d_hi * np.sin(tf)[None, :] / (rc[:, None] ** 2 * np.sin(tc[1:])[None, :] * dth ** 2)).ravel())
    if W:                                                         # centred advection, face fluxes from psi at corners
        psi = lambda r, th: -W * R ** 2 * np.sin(np.pi * (r - r_in) / H) * np.sin(th) ** 2 * np.cos(th)
        V = rc[:, None] ** 2 * np.sin(tc)[None, :] * dr * dth
        tb = np.arange(nth + 1) * dth
        Fr = psi(rf[:, None], tb[None, 1:]) - psi(rf[:, None], tb[None, :-1])          # outward, (NR-1, NTH)
        rb = r_in + np.arange(nr + 1) * dr
        Ft = -(psi(rb[1:, None], tf[None, :]) - psi(rb[:-1, None], tf[None, :]))       # equatorward, (NR, NTH-1)
        for (ka, kb, F, Va, Vb) in ((k(I[:-1], J[:-1]), k(I[1:], J[1:]), Fr, V[:-1], V[1:]),
                                    (k(I[:, :-1], J[:, :-1]), k(I[:, 1:], J[:, 1:]), Ft, V[:, :-1], V[:, 1:])):
            ka, kb, F, Va, Vb = ka.ravel(), kb.ravel(), F.ravel(), Va.ravel(), Vb.ravel()
            rows.extend([ka, ka, kb, kb]); cols.extend([ka, kb, ka, kb])
            vals.extend([-F / (2 * Va), -F / (2 * Va), F / (2 * Vb), F / (2 * Vb)])
    n = nr * nth
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsc()
    lu = splu((identity(n, format="csc") - dt * A).tocsc())

    r0 = R - 0.25
    Dc = np.exp(np.log(d_hi) - slow(rc[:, None], tc[None, :]) * np.log(ratio))
    a = np.tanh((rc[:, None] - r0) / np.sqrt(2 * Dc)).ravel()
    vol = (rc[:, None] ** 2 * np.sin(tc)[None, :]).ravel()
    i_mid = np.searchsorted(rc, r_in + 0.3)                       # side wall measured 0.3 above the floor
    t, hist = 0.0, []
    while t < t_max:
        a = lu.solve(a + dt * (a - a ** 3))
        t += dt
        if round(t / dt) % 50 == 0:
            frac = vol[a > 0].sum() / vol.sum()
            hist.append((t, frac))
            if a.min() > 0:
                return "TAKEOVER", t, 1.0, H, 90.0
            if a.max() < 0:
                return "DEEP HAND WON", t, 0.0, 0.0, 0.0
            old = [f for tt, f in hist if tt <= t - 500 * slow_t]
            if t > 1500 * slow_t and old and abs(frac - old[-1]) < 2e-4:
                break
    g = a.reshape(nr, nth)
    pole = R - rc[np.flatnonzero(g[:, 0] > 0)[0]] if (g[:, 0] > 0).any() else 0.0
    side = np.degrees(tc[np.flatnonzero(g[i_mid] > 0)[-1]]) if (g[i_mid] > 0).any() else 0.0
    return ("STUCK" if t < t_max else "unfinished"), t, hist[-1][1], pole, side


def show(label, ratio, gap, **kw):
    t0 = time.time()
    out, t, frac, pole, side = run(ratio, gap, **kw)
    print(f"  {label:30s} ratio {ratio:<4g} gap {gap:3g} deg -> {out:10s} t={t:6.0f}  + share {frac:.3f}  "
          f"+ reaches depth {pole:.2f} at the pole, side wall at {side:4.1f} deg  [{time.time() - t0:.0f}s]", flush=True)
    return out


def main():
    t0 = time.time()
    sec = sys.argv[1:] or ["1", "2", "3"]
    if "1" in sec:
        print(f"(1) controls: R={R}, depth {H}, lid {LAYER} thick, D_hi={D_HI}, grid {NR}x{NTH}, dt={DT}")
        assert show("no lid", 1, 0, lid=False) == "TAKEOVER"
        assert show("no lid, dt/2", 1, 0, lid=False, dt=DT / 2) == "TAKEOVER"
        assert show("lid, no gap", 10, 0) == "STUCK"
    if "2" in sec:
        print("\n(2) polar gap; predicted critical gap 27 deg (ratio 10), 58 deg (ratio 100)")
        for ratio, gap, want in ((10, 15, "STUCK"), (10, 20, "STUCK"), (10, 25, "STUCK"), (10, 30, "TAKEOVER"),
                                 (10, 35, "TAKEOVER"), (100, 45, "STUCK"), (100, 55, "STUCK"), (100, 70, "TAKEOVER")):
            assert show("lid with polar gap", ratio, gap, t_max=16000) == want, f"ratio {ratio} gap {gap}: outcome changed"
    if "3" in sec:
        print("\n(3) bias needed to push the wall through a full lid (two_zone.py 1D, ratio 10, lid on top)")
        lo, hi = 0.0, 0.4
        for _ in range(6):
            mid = (lo + hi) / 2
            done = two_zone.run(10, "top", h=mid)[3]
            print(f"  h = {mid:.4f}: {'pushed through' if done else 'pinned'}", flush=True)
            lo, hi = (lo, mid) if done else (mid, hi)
        print(f"  critical bias h_c between {lo:.4f} and {hi:.4f}; real parity bias g ~ 1e-19..1e-16")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
