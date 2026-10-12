"""Objection 1 (Kibble-Zurek / patch sum): does the WHOLE-volume Delta still predict the hand when the ocean
first breaks into many independent patches?

Same field model as domains.py / coarsen3d.py (Eq. 14 of journals/ija_draft.md), periodic grid, cell size 1:
    d(alpha_i) = (lam*alpha_i - alpha_i**3 + g + D*lap(alpha)_i) dt + sqrt(eps) dW_i,   lam = gamma*t, -1 -> +1
eps is the noise PER CELL, so the uniform mode (the spatial mean) feels eps_total = eps / (number of cells) and the
whole-system prediction is  P = Phi(Delta_total),  Delta_total = sim.delta(g, eps / cells, gamma).
gamma = 0.01, eps = 1e-4 per cell throughout. Coherence length xi = sqrt(D) * gamma**-1/4 (domains.py check A:
the 1D domain length at freeze-out is C_KZ = 12.81 xi, run_all.log). The "final hand" is the majority hand (sign of
the mean of sign(alpha)); in 1D walls barely move, so the majority IS what perfect coarsening would deliver.

(A1) 1D, D = 4 (xi = 6.3 cells), L = 0.5 .. 512 xi (M = L / (12.8 xi) up to 40 patches), Delta_total = 0.5, 1,
     1.5, 2.5 (g set from L), 600 runs per point. Fit P_sim = Phi(c * Delta_total) by maximum likelihood per L.
(A2) 2D, D = 1 (xi = 3.2 cells), n = 8 .. 64 cells per side, Delta_total = 1, 2, 400 runs per point. Majority
     recorded at lam = 0.4 (just after selection), at the end of the sweep, and after holding lam = 1 for t = 500
     (curvature-driven coarsening, as coarsen3d.py).
(B)  2D coarsening from a random patchwork (coarsen3d.patchwork, 2D), g = 0, no noise, n = 64, patch filter width
     w = 2 and 4. sigma_m = std of the volume fraction of a g=0 patchwork (exact, from the lattice correlation);
     M_eff = 1 / (4 sigma_m^2), so m in units of 1/sqrt(M_eff) is 2x. Start with imbalance m = x * sigma_m exactly
     (threshold at the quantile), x = 0 .. 3, 300 runs; report P(final majority = favoured) at t = 150, 300, 600.

Theory written down before the runs. The patch-sum argument (M equal, independent, fully +-1 patches; coarsening
keeps the majority) gives c = sqrt(2/pi) = 0.798. But the field at freeze-out is a smooth Gaussian field, not a set of
boxes: the mean of sign(field) has mean sqrt(2/pi)*mu/sigma and variance (1/V) int (2/pi) arcsin(rho), while the
field mean has variance (1/V) int rho. Hence c = [int rho / int arcsin(rho)]^(1/2), independent of the patch size.
For Gaussian rho: c = 0.914 (1D), 0.958 (2D), 0.979 (3D); boxes give arcsin(1) = pi/2 and recover 0.798.

PREDICTIONS (written 2026-10-12, before any result-producing run):
 P1. Uniform limit (1D, L = 3 cells ~ 0.5 xi, all non-uniform modes damped): P = Phi(Delta_total) within 3 SE
     (asserted in check()).
 P2. 1D: the sqrt(N) scaling survives patching. For L >= 128 xi (M >= 10), P_sim = Phi(c * Delta_total) with one c:
     c(128 xi) and c(512 xi) agree within 2 combined SE.
 P3. 1D plateau: c = 0.91 +- 0.05 (Gaussian-field sign-majority theory), and the patch-sum value 0.80 is excluded by
     more than 3 SE at L = 512 xi.
 P4. 1D: c(L) falls from ~1 to the plateau; c >= 0.97 for L <= 2 xi.
 P5. 1D: the majority at lam = 0.4 equals the majority at the end of the sweep in >= 99% of runs (walls frozen).
 P6. 2D, majority at lam = 0.4: c_early = 0.96 +- 0.05 at n = 64.
 P7. 2D, after coarsening (t = 500 at lam = 1): c_final < c_early, and c_final lies in [0.6, 0.9] at n = 64.
 P8. (B) P(final favoured | x) is monotone in x, equals 0.5 at x = 0 within 2 SE, and is >= 0.95 at x = 3. A probit
     fit P = Phi(kappa * x) gives kappa in [0.7, 3]. The w = 2 (M_eff ~ 117) and w = 4 (M_eff ~ 29) curves agree
     within 2 combined SE at every x (they collapse in units of 1/sqrt(M)).
 P9. Consistency of A2 and B: c_final = c_early * kappa / sqrt(1 + kappa^2) within +- 0.1 (coarsening acts on the
     post-selection imbalance as B's majority curve).
 P10. Convergence: halving dt (1D at L = 128 xi; 2D at n = 32) and halving dx (1D at L = 128 xi; 2D at n = 16)
     each change P by less than 2 combined SE.

Usage: python kz_patches.py [1d] [2d] [imb2] [imb4] [check]   (default: all; check always runs last)

RESULTS (kz_patches.log, 2026-10-12):
 Four sections run in parallel, each exit 0; check() passed in every section (P = 0.8478 +- 0.0057 vs 0.8413,
 +1.13 SE). Wall times 1407 s (1d), 1245 s (2d), 964 s (imb2), 961 s (imb4); the 1d and 2d sections ran over the
 20-minute budget under 9-way parallel load.
 P1 HELD. Uniform limit +1.13 SE.
 P2 HELD. c(128 xi) = 0.876 +- 0.029, c(512 xi) = 0.882 +- 0.029 (end of sweep); diff 0.15 combined SE.
 P3 FAILED (second clause). Plateau c = 0.879 +- 0.021 (128 and 512 xi pooled; theory 0.914, -1.7 SE) lies in
    0.91 +- 0.05, but at 512 xi alone 0.80 is excluded by only 2.90 SE, short of the 3 SE required. (Pooling the
    two plateau points gives 3.95 SE; that pooling was not pre-registered.)
 P4 FAILED. c(0.5 xi) = 0.986 +- 0.033 but c(2 xi) = 0.936 +- 0.031 < 0.97. c(L) is not monotone either:
    c(8 xi) = 1.015 +- 0.034, c(32 xi) = 0.941 +- 0.031.
 P5 HELD. early == end majority in 1.0000, 1.0000, 0.9988, 0.9921, 0.9900, 0.9908 of runs (0.5 .. 512 xi);
    the 128 xi value sits exactly on the 0.99 threshold.
 P6 FAILED. c_early(n = 64) = 1.031 +- 0.053, outside 0.96 +- 0.05; it is consistent with 0.958 at 1.4 SE.
 P7 FAILED. c_final = 0.981 +- 0.050 < c_early = 1.031, but c_final is not in [0.6, 0.9]. Coarsening during the
    hold changed nothing: P_end = 0.843 / 0.973 vs P_final = 0.844 / 0.973 (Delta = 1 / 2), from the same runs,
    while the single-hand fraction rose 0.15 -> 0.81.
 P8 FAILED (collapse clause). Monotone in x, P(x = 0) = 0.483 +- 0.029 (w = 2) and 0.533 +- 0.029 (w = 4), P(3) = 1.000,
    kappa = 1.45 +- 0.07 (w = 2) and 2.12 +- 0.12 (w = 4), both in [0.7, 3]. But the curves do not collapse in units of
    1/sqrt(M_eff): at x = 0.5, P = 0.790 +- 0.024 vs 0.915 +- 0.016 (4.3 combined SE). P is the same at t = 150, 300 and
    600, so the whole loss happens in the first relaxation of the thresholded patchwork, not in late coarsening.
 P9 AMBIGUOUS. kappa (w) was not specified. c_early * kappa/sqrt(1 + kappa^2) = 0.849 (w = 2, off by 0.13: fails)
    or 0.933 (w = 4, off by 0.05: holds) against c_final = 0.981. The mechanism behind P9 (coarsening reduces c) is
    contradicted by P7: c_end = 0.979 and c_final = 0.981.
 P10 HELD. 1D dt/2 -1.03 SE, 1D dx/2 +0.23 SE, 2D dt/2 -1.70 SE, 2D dx/2 +1.29 SE (all < 2).
 Net: the whole-volume Delta survives patching with a factor c on Delta. In 1D, c = 0.88 +- 0.02, flat from 10 to 40
 patches. In 2D, c = 0.98 +- 0.05 after coarsening. The box patch-sum value 0.80 is a lower bound, not the answer.
"""
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import norm

from domains import lap as lap1
from sim import delta

GAMMA, EPS = 0.01, 1e-4
C_KZ = 12.81                                   # domains.py check A, run_all.log
C_PATCH, C_GAUSS = np.sqrt(2 / np.pi), {1: 0.914, 2: 0.958}


def xi(D):
    return np.sqrt(D) * GAMMA ** -0.25


def g_for(Delta, volume):
    """Bias giving Delta_total for a system of `volume` base cells (sim.delta is linear in g)."""
    return Delta / delta(1.0, EPS / volume, GAMMA)


def lap2(a):
    return sum(np.roll(a, s, ax) for ax in (-2, -1) for s in (1, -1)) - 4 * a


def majority(a):
    s = np.sign(a).reshape(len(a), -1).mean(1)
    return (s > 0) + 0.5 * (s == 0), np.abs(s) == 1          # score (ties count 1/2), single hand


def sweep(n, dim, D, g, dt, seed, h=1.0, hold=0.0):
    """Periodic grid, n cells per side of size h (base units). lam = gamma*t from -1 to +1, then held at 1 for `hold`.
    g: per-run bias, shape (runs,). Noise per cell EPS / h**dim, so the noise per base volume is EPS for any h.
    Returns {'early': lam=0.4, 'end': lam=1, 'final': after hold} -> (score, single)."""
    Dl = D / h ** 2
    assert Dl * dt <= 1 / (2 * dim) + 1e-12, "Euler stability D dt/dx^2 <= 1/(2 dim)"
    rng = np.random.default_rng(seed)
    runs = len(g)
    a = np.zeros((runs,) + (n,) * dim, np.float32)
    gb = np.asarray(g, np.float32).reshape((runs,) + (1,) * dim)
    lap = lap1 if dim == 1 else lap2
    noise, dtf, Dl = np.float32(np.sqrt(EPS / h ** dim * dt)), np.float32(dt), np.float32(Dl)
    nsweep, nhold = int(round(2 / GAMMA / dt)), int(round(hold / dt))
    out = {}
    for i in range(nsweep + nhold):
        lam = np.float32(min(1.0, -1.0 + GAMMA * dt * i))
        a += (lam * a - a ** 3 + gb + Dl * lap(a)) * dtf + noise * rng.standard_normal(a.shape, dtype=np.float32)
        if "early" not in out and lam >= 0.4:
            out["early"] = majority(a)
        if i == nsweep - 1:
            out["end"] = majority(a)
    out["final"] = majority(a)
    return out


def fit_probit(x, k, n):
    """MLE of c in P = Phi(c x) from k successes of n; SE from the curvature of the log-likelihood."""
    x, k, n = map(np.asarray, (x, k, n))
    nll = lambda c: -np.sum(k * norm.logcdf(c * x) + (n - k) * norm.logsf(c * x))
    c = minimize_scalar(nll, bounds=(0.02, 5), method="bounded").x
    e = 1e-3
    return c, 1 / np.sqrt(max((nll(c + e) - 2 * nll(c) + nll(c - e)) / e ** 2, 1e-12))


def pse(score):
    p = float(np.mean(score))
    return p, np.sqrt(max(p * (1 - p), 0.25 / len(score)) / len(score))   # floor keeps SE > 0 at p = 0 or 1


def batch(n, dim, D, deltas, runs, dt, seed, h=1.0, hold=0.0):
    """One sweep with all Deltas stacked along the run axis. Returns {key: [(Delta, score, single), ...]}."""
    vol = (n * h) ** dim
    g = np.repeat([g_for(d, vol) for d in deltas], runs)
    out = sweep(n, dim, D, g, dt, seed, h=h, hold=hold)
    return {k: [(d, s[j * runs:(j + 1) * runs], o[j * runs:(j + 1) * runs]) for j, d in enumerate(deltas)]
            for k, (s, o) in out.items()}


def report(label, rows, keys):
    """Print one table block; return {key: (c, se)} fitted over the Deltas."""
    fits = {}
    for key in keys:
        ds = [d for d, _, _ in rows[key]]
        ks = [s.sum() for _, s, _ in rows[key]]
        ns = [len(s) for _, s, _ in rows[key]]
        c, se = fit_probit(ds, ks, ns)
        fits[key] = (c, se)
        cells = "  ".join(f"D={d:<4g}P={pse(s)[0]:.3f}+-{pse(s)[1]:.3f}(1hand {o.mean():.2f})" for d, s, o in rows[key])
        print(f"  {label:22s} {key:5s} n={ns[0]}/pt  {cells}  ->  c = {c:.3f} +- {se:.3f}", flush=True)
    return fits


def theory_line(deltas):
    return "  predictions          " + "  ".join(
        f"D={d:<4g}Phi(D)={norm.cdf(d):.3f} Phi(.914D)={norm.cdf(C_GAUSS[1] * d):.3f} Phi(.798D)={norm.cdf(C_PATCH * d):.3f}"
        for d in deltas)


def sec_1d():
    D, runs, dt, deltas = 4.0, 600, 0.05, (0.5, 1.0, 1.5, 2.5)
    print(f"(A1) 1D periodic, D={D}, xi={xi(D):.2f} cells, KZ domain {C_KZ * xi(D):.0f} cells, gamma={GAMMA}, "
          f"eps={EPS}/cell, dt={dt}", flush=True)
    print(theory_line(deltas))
    fits = {}
    for j, Lx in enumerate((0.5, 2, 8, 32, 128, 512)):
        t0 = time.time()
        n = int(round(Lx * xi(D)))
        rows = batch(n, 1, D, deltas, runs, dt, seed=100 + j)
        agree = np.mean(np.concatenate([e[1] == f[1] for e, f in zip(rows["early"], rows["end"])]))
        fits[Lx] = report(f"L={Lx:g}xi ({n} cells, M={Lx / C_KZ:.2f})", rows, ("early", "end"))
        print(f"      early==end majority in {agree:.4f} of runs  [{time.time() - t0:.0f}s]", flush=True)
        if Lx == 128:
            base = rows["end"][1]
    print("  c(L) at end of sweep: " + "  ".join(f"{Lx:g}xi:{c:.3f}+-{se:.3f}" for Lx, ((c, se)) in
                                                 ((k, v["end"]) for k, v in fits.items())))
    # P10 convergence at L = 128 xi, Delta = 1
    n = int(round(128 * xi(D)))
    p0, s0 = pse(base[1])
    for label, kw in (("dt/2", dict(dt=dt / 2)), ("dx/2 (dt/4)", dict(dt=dt / 4, h=0.5))):
        h = kw.get("h", 1.0)
        r = batch(int(round(n / h)), 1, D, (1.0,), runs, kw["dt"], seed=190, h=h)["end"][0][1]
        p, s = pse(r)
        print(f"  convergence L=128xi Delta=1 {label:12s} P={p:.3f}+-{s:.3f} vs base {p0:.3f}+-{s0:.3f}: "
              f"diff {(p - p0) / np.hypot(s, s0):+.2f} combined SE", flush=True)
    return fits


def sec_2d():
    D, runs, dt, hold, deltas = 1.0, 400, 0.1, 500.0, (1.0, 2.0)
    print(f"\n(A2) 2D periodic, D={D}, xi={xi(D):.2f} cells, dt={dt}, hold lam=1 for t={hold:g}", flush=True)
    print("  predictions          " + "  ".join(f"D={d:<4g}Phi(D)={norm.cdf(d):.3f} Phi(.958D)={norm.cdf(.958 * d):.3f}"
                                              for d in deltas))
    base = {}
    for j, n in enumerate((8, 16, 32, 64)):
        t0 = time.time()
        rows = batch(n, 2, D, deltas, runs, dt, seed=200 + j, hold=hold)
        report(f"n={n} (L={n / xi(D):.1f}xi)", rows, ("early", "end", "final"))
        print(f"      [{time.time() - t0:.0f}s]", flush=True)
        base[n] = rows["final"][0][1]
    for label, n, kw in (("dt/2", 32, dict(dt=dt / 2)), ("dx/2 (dt/4)", 16, dict(dt=dt / 4, h=0.5))):
        h = kw.get("h", 1.0)
        r = batch(int(round(n / h)), 2, D, (1.0,), runs, kw["dt"], seed=290 + n, h=h, hold=hold)["final"][0][1]
        p, s = pse(r)
        p0, s0 = pse(base[n])
        print(f"  convergence n={n} Delta=1 final {label:12s} P={p:.3f}+-{s:.3f} vs base {p0:.3f}+-{s0:.3f}: "
              f"diff {(p - p0) / np.hypot(s, s0):+.2f} combined SE", flush=True)


def patch_filter(n, w):
    k = np.fft.fftfreq(n)
    return np.exp(-2 * (np.pi * w) ** 2 * (k[:, None] ** 2 + k[None, :] ** 2))   # coarsen3d.patchwork, in 2D


def sigma_m(n, w):
    """Exact std of the +hand volume fraction of sign(z), z the filtered noise: Var = (1/n^2) sum_r arcsin(rho_r)/(2 pi)."""
    rho = np.fft.ifft2(patch_filter(n, w) ** 2).real
    rho /= rho[0, 0]
    return float(np.sqrt(np.arcsin(np.clip(rho, -1, 1)).sum() / (2 * np.pi) / n ** 2))


def sec_imb(w, n=64, runs=300, dt=0.1, xs=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0), times=(150, 300, 600)):
    sm = sigma_m(n, w)
    rng = np.random.default_rng(300 + w)
    H = patch_filter(n, w)
    zchk = np.fft.ifft2(np.fft.fft2(rng.standard_normal((400, n, n))) * H).real
    print(f"\n(B) 2D coarsening from a patchwork, n={n}, w={w}: sigma_m={sm:.4f} (sampled {(zchk > 0).mean((1, 2)).std():.4f}),"
          f" M_eff={1 / (4 * sm * sm):.0f}; D=1, dt={dt}, g=0, no noise", flush=True)
    res = {}
    for x in xs:
        t0 = time.time()
        z = np.fft.ifft2(np.fft.fft2(rng.standard_normal((runs, n, n))) * H).real
        f = 0.5 + x * sm
        a = np.where(z > np.quantile(z, 1 - f, axis=(1, 2), keepdims=True), 1.0, -1.0).astype(np.float32)
        m0 = (a > 0).mean((1, 2)) - 0.5
        t, rec = 0.0, []
        for T in times:
            while t < T - 1e-9:
                a += (a - a ** 3 + lap2(a)) * np.float32(dt)
                t += dt
            rec.append(majority(a))
        res[x] = rec
        print(f"  x={x:<4g} m={m0.mean():.4f} (m*sqrt(M)={x / 2:.2f}): " + "  ".join(
            f"t={T}: P={pse(s)[0]:.3f}+-{pse(s)[1]:.3f} 1hand {o.mean():.2f}" for T, (s, o) in zip(times, rec))
              + f"  [{time.time() - t0:.0f}s]", flush=True)
    for i, T in enumerate(times):
        xx = [x for x in xs if x > 0]
        c, se = fit_probit(xx, [res[x][i][0].sum() for x in xx], [runs] * len(xx))
        print(f"  probit fit t={T}: kappa = {c:.3f} +- {se:.3f}  (P = Phi(kappa x), x>0 points)")
    return res


def check():
    """Uniform fast-sweep limit: 3 cells, D=4 (non-uniform modes at -3D, never unstable) must give sim.py's Phi(Delta)."""
    runs, Delta = 4000, 1.0
    s = batch(3, 1, 4.0, (Delta,), runs, 0.05, seed=7)["end"][0][1]
    p, se = pse(s)
    want = norm.cdf(delta(g_for(Delta, 3), EPS / 3, GAMMA))
    print(f"\ncheck: uniform limit (3 cells, D=4) P={p:.4f}+-{se:.4f}, sim.py Phi(Delta)={want:.4f}, "
          f"diff {(p - want) / se:+.2f} SE")
    assert abs(p - want) < 3 * se, "uniform limit does not reproduce sim.py Phi(Delta)"


def main():
    secs = sys.argv[1:] or ["1d", "2d", "imb2", "imb4"]
    t0 = time.time()
    if "1d" in secs:
        sec_1d()
    if "2d" in secs:
        sec_2d()
    for w in (2, 4):
        if f"imb{w}" in secs:
            sec_imb(w)
    check()
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
