"""h13_sphere: does a shallow ocean on a SPHERE coarsen sideways to one hand in the time the sweep assumes? (H13)

Model: the normal form da/dt = a - a^3 + D lap(a) (coarsen3d.py) on the surface of a sphere of radius R (a shallow
ocean mixed over its depth is a 2D surface), with a flat periodic square as the control (the aniso3d.py slab,
whose two flat walls persisted in 4 of 8 runs). Solver: spherical harmonics on a Gauss-Legendre grid (flat square:
FFT), SBDF2 time stepping (diffusion implicit, a - a^3 explicit), dt = 0.25.
Wall length P = (3 sqrt(2D)/4) int |grad a|^2 dA (tanh wall); domain length l = 2 Area / P.
One hand = minority area fraction < 1%.  Sweep's term: tau_sw = (pi R / 5.5)^2 / D (uncertainty.py, "across").
Seeds: run k of a batch uses seed0 + k.  Start: smoothed noise (width 4) cut at area fraction f0 (+ is "favoured").

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python h13_sphere.py check
     ... python h13_sphere.py random kind=sphere size=16 f0=0.5 n=150 seed=0 procs=4 > h13_sphere_<name>.log
     ... python h13_sphere.py special size=16 > h13_sphere_special.log

Pre-registration: h13_sphere_predictions.txt, written Wed Oct 21 00:37:01 PDT 2026, before any results run,
  sha256 b316fc86801e185cedeecc3d3f6689cd4f5eaaf8fbba78fbed5ccdc5667373ae (asserted below by need_pred()).
Smoke runs before that (disclosed; none ran a random patchwork to the end):
  - cap-law convergence tests (offset circle, R = 16 and 32; dt 0.5/0.25/0.125; lmax 48/64/80), timing at R = 16, 32;
  - two short random patchworks (seed 999, f0 = 0.5, R = 16 and 32, stopped at t = 100): l(100) = 63.53 and 62.30.
"""
import hashlib
import sys
import time
import warnings
from functools import lru_cache
from multiprocessing import Pool

import numpy as np
from scipy.special import sph_harm_y

A_SW = 5.5  # uncertainty.py: A
DT = 0.25
PRED = "h13_sphere_predictions.txt"
PRED_SHA, PRED_DATE = "b316fc86801e185cedeecc3d3f6689cd4f5eaaf8fbba78fbed5ccdc5667373ae", "Wed Oct 21 00:37:01 PDT 2026"


class Sphere:
    def __init__(self, R, lmax=None):
        self.R, self.M = R, int(lmax or 4 * R) + 1
        self.nlat = self.M + self.M % 2  # even, so nlon = 2 nlat is divisible by 4 (exact quarter-turn symmetry)
        self.nlon = 2 * self.nlat
        x, w = np.polynomial.legendre.leggauss(self.nlat)
        self.theta, self.phi = np.arccos(x), 2 * np.pi * np.arange(self.nlon) / self.nlon
        P = np.zeros((self.M, self.M, self.nlat))  # P[m, l, j] = Y_l^m(theta_j, 0), zero for l < m
        for m in range(self.M):
            P[m, m:] = sph_harm_y(np.arange(m, self.M)[:, None], m, self.theta[None, :], 0.0).real
        self.Pw, self.PT = P * w, np.ascontiguousarray(P.transpose(0, 2, 1))
        ell, m = np.arange(self.M), np.arange(self.M)
        self.valid = (ell[None, :] >= m[:, None])[..., None]
        self.ll1 = (ell * (ell + 1.0))[None, :, None] * np.ones((self.M, 1, 1))
        self.k2 = self.ll1 / R**2  # minus the Laplacian's eigenvalue
        self.wm = np.where(m == 0, 1.0, 2.0)[:, None, None]  # m and -m
        self.W = np.outer(w, np.ones(self.nlon)) * R**2 * 2 * np.pi / self.nlon  # area of each grid cell
        self.area = 4 * np.pi * R**2

    def fwd(self, a):  # grid -> coefficients (m, l, re/im)
        X = np.fft.rfft(a, axis=1)[:, : self.M] * (2 * np.pi / self.nlon)
        return self.Pw @ np.stack([X.real, X.imag], -1).transpose(1, 0, 2)

    def inv(self, c):
        F = self.PT @ c
        X = np.zeros((self.nlat, self.nlon // 2 + 1), complex)
        X[:, : self.M] = (F[..., 0] + 1j * F[..., 1]).T * self.nlon
        return np.fft.irfft(X, self.nlon, axis=1)

    def energy(self, c):  # int |grad a|^2 dA = sum l(l+1) |c_lm|^2  (R cancels)
        return (self.wm * self.ll1 * c**2).sum()

    def noise(self, rng):  # isotropic random field, in coefficients
        c = rng.standard_normal((self.M, self.M, 2)) * self.valid
        c[0, :, 1] = 0
        return c


class Torus:
    def __init__(self, L, n=None):
        self.L, self.n, self.area = L, n or 2 * int(round(L / 1.6)), L * L  # dx ~ 0.8
        k = 2 * np.pi * np.fft.fftfreq(self.n, L / self.n)
        kr = 2 * np.pi * np.fft.rfftfreq(self.n, L / self.n)
        self.k2 = k[:, None] ** 2 + kr[None, :] ** 2
        j = np.arange(self.n // 2 + 1)
        self.wk = np.where((j == 0) | (j == self.n // 2), 1.0, 2.0)[None, :]
        self.W = np.full((self.n, self.n), L * L / self.n**2)

    def fwd(self, a):
        return np.fft.rfft2(a)

    def inv(self, c):
        return np.fft.irfft2(c, s=(self.n, self.n))

    def energy(self, c):
        return self.L**2 / self.n**4 * (self.wk * self.k2 * np.abs(c) ** 2).sum()

    def noise(self, rng):
        return self.fwd(rng.standard_normal((self.n, self.n)))


@lru_cache(None)
def grid(kind, size):
    return Sphere(size) if kind == "sphere" else Torus(size)


def tau(kind, size, D=1.0):
    return ((np.pi * size if kind == "sphere" else size) / A_SW) ** 2 / D


def patchwork(g, f0, rng, width=4.0):
    f = g.inv(g.noise(rng) * np.exp(-g.k2 * width**2 / 2))
    o = np.argsort(f, axis=None)[::-1]
    cw = np.cumsum(g.W.ravel()[o]) / g.area
    return np.where(f > f.ravel()[o][np.searchsorted(cw, f0)], 1.0, -1.0)


def evolve(g, a0, D=1.0, dt=DT, tmax=1e9, rec=2.0, stop=0.01, probe=None):
    """SBDF2 steps of da/dt = a - a^3 + D lap(a). Returns t, domain length l, minority fraction, T1, winner."""
    lam, pw = D * g.k2, 0.75 * np.sqrt(2 * D)
    c = g.fwd(a0)
    a = g.inv(c)
    Nm = g.fwd(a - a**3)
    cm, c = c, (c + dt * Nm) / (1 + dt * lam)  # one Euler step to start
    out = dict(t=[], ell=[], mino=[], probe=[], T1=np.nan, winner=0)
    nrec, step = max(1, round(rec / dt)), 1
    while step * dt <= tmax:
        a = g.inv(c)
        if step % nrec == 0:
            fp = (g.W * (a > 0)).sum() / g.area
            out["t"].append(step * dt), out["mino"].append(min(fp, 1 - fp))
            out["ell"].append(2 * g.area / (pw * g.energy(c)))
            if probe:
                out["probe"].append(probe(c))
            out["winner"] = 1 if fp > 0.5 else -1
            if out["mino"][-1] < stop:
                out["T1"] = step * dt
                break
        N = g.fwd(a - a**3)
        cm, c, Nm = c, (4 * c - cm + 2 * dt * (2 * N - Nm)) / (3 + 2 * dt * lam), N
        step += 1
    return {k: np.array(v) if isinstance(v, list) else v for k, v in out.items()}


def cap_ic(g, theta0, D=1.0):  # + inside a cap of colatitude theta0, tanh wall
    return np.tanh(g.R * (theta0 - g.theta[:, None] + 0 * g.phi[None, :]) / np.sqrt(2 * D))


def wavy_ic(g, eps, m, D=1.0):  # wall at colatitude pi/2 + eps cos(m phi)
    th = g.theta[:, None] + 0 * g.phi[None, :]
    return np.tanh(g.R * (np.pi / 2 + eps * np.cos(m * g.phi[None, :]) - th) / np.sqrt(2 * D))


def cap_T1(R, u0, D=1.0):  # sharp-interface cap law: cos(theta) = u0 exp(D t / R^2); minority 1% <=> cos(theta) = 0.98
    return R**2 / D * np.log(0.98 / u0)


def check():
    g = Sphere(2.0, 24)
    G = 2 * np.pi * (g.Pw @ g.PT)  # harmonics orthonormal on the grid
    assert abs(G - np.eye(g.M)[None] * g.valid).max() < 1e-10
    rng = np.random.default_rng(1)
    c = g.noise(rng)
    assert abs(g.fwd(g.inv(c)) - c).max() < 1e-12  # round trip
    x, y, z = [f(g.theta)[:, None] * h(g.phi)[None, :] for f, h in
               ((np.sin, np.cos), (np.sin, np.sin), (np.cos, lambda p: np.ones_like(p)))]
    s2 = x**2 + y**2
    for l, h, gr2 in ((5, ((x + 1j * y) ** 5).real, 25 * s2**4), (3, x * y * z, y**2 * z**2 + x**2 * z**2 + x**2 * y**2)):
        E = g.energy(g.fwd(h))
        E_ref = (g.W / g.R**2 * (gr2 - l * l * h**2)).sum()  # |grad_S h|^2 = |grad h|^2 - (l h)^2 on the unit sphere
        rq = E / (g.W * h**2).sum() * g.R**2
        print(f"l={l}: energy {E:.8f} vs polynomial calculus {E_ref:.8f}; Rayleigh quotient * R^2 = {rq:.8f} (l(l+1) = {l * (l + 1)})")
        assert abs(E / E_ref - 1) < 1e-9 and abs(rq / (l * (l + 1)) - 1) < 1e-9
    t = Torus(20.0)
    h = np.cos(2 * np.pi * 3 * np.arange(t.n) / t.n)[:, None] * np.ones((1, t.n))
    rq = t.energy(t.fwd(h)) / (t.W * h**2).sum()
    print(f"torus: Rayleigh quotient {rq:.8f} vs (2 pi 3 / L)^2 = {(2 * np.pi * 3 / 20) ** 2:.8f}")
    assert abs(rq / (2 * np.pi * 3 / 20) ** 2 - 1) < 1e-9
    g = grid("sphere", 16)
    th0 = 1.2
    r = evolve(g, cap_ic(g, th0), tmax=600)
    print(f"cap law R=16 theta0={th0}: T1 {r['T1']:.1f} vs sharp-interface {cap_T1(16, np.cos(th0)):.1f}")
    assert abs(r["T1"] / cap_T1(16, np.cos(th0)) - 1) < 0.05


def run_one(job):
    kind, size, f0, seed, tmax, stop = job
    g = grid(kind, size)
    r = evolve(g, patchwork(g, f0, np.random.default_rng(seed)), tmax=tmax, stop=stop)
    return seed, r


def need_pred():
    h = hashlib.sha256(open(PRED, "rb").read()).hexdigest()
    assert PRED_SHA and h == PRED_SHA, f"{PRED} does not match the hash recorded in this script: {h}"


def random_job(kind="sphere", size="16", f0="0.5", n="100", seed="0", procs="4", tmaxtau="30", tmax=None, stop="0.01",
               fitlo="30"):
    need_pred()
    size, f0, n, seed, procs, stop = int(size), float(f0), int(n), int(seed), int(procs), float(stop)
    ta = tau(kind, size)
    tm = float(tmax) if tmax else float(tmaxtau) * ta
    Lc = np.pi * size if kind == "sphere" else size
    print(f"# {kind} size={size} f0={f0} n={n} seed0={seed} tmax={tm:.0f} (= {tm / ta:.1f} tau_sw) tau_sw={ta:.2f} R^2/D={size**2:.0f}")
    t0, res = time.time(), []
    with Pool(procs) as p:
        for s, r in p.imap_unordered(run_one, [(kind, size, f0, seed + k, tm, stop) for k in range(n)]):
            res.append((s, r))
            print(f"run seed={s} T1={r['T1']:.1f} T1/tau_sw={r['T1'] / ta:.3f} winner={r['winner']:+d} final_minority={r['mino'][-1]:.4f}", flush=True)
    res.sort(key=lambda x: x[0])
    T = np.array([r["T1"] for _, r in res])
    Tn = np.where(np.isnan(T), np.inf, T) / ta
    q = lambda p: np.percentile(Tn, p, method="higher")
    print(f"SUMMARY {kind} size={size} f0={f0}: n={n} reached one hand {np.isfinite(Tn).sum()}/{n} by {tm / ta:.0f} tau_sw; "
          f"T1/tau_sw median {q(50):.2f} p90 {q(90):.2f} p99 {q(99):.2f} max {Tn.max():.2f} mean(finished) {Tn[np.isfinite(Tn)].mean():.2f}")
    if kind == "sphere":
        W = np.array([r["winner"] for _, r in res])
        print(f"SUMMARY favoured (+) hand wins {np.mean(W > 0):.3f} ({(W > 0).sum()}/{n})")
        ok = np.sort(Tn[np.isfinite(Tn)]) * ta / size**2  # in units of R^2/D
        i = np.arange(len(ok))
        sel = ok >= np.median(ok)
        if sel.sum() > 5:
            s = np.polyfit(ok[sel], np.log(1 - (i[sel] + 0.5) / n), 1)[0]
            print(f"SUMMARY tail: survival e-folds in {-1 / s:.3f} R^2/D (n_tail={sel.sum()}); R^2/D = {size**2 / ta:.3f} tau_sw")
    L = max(len(r["t"]) for _, r in res)
    ell = np.full((len(res), L), np.nan)
    for k, (_, r) in enumerate(res):
        ell[k, : len(r["ell"])] = r["ell"]
    t = np.arange(1, L + 1) * 2.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gm = np.exp(np.nanmean(np.log(ell), axis=0))
        alive = np.isfinite(ell).mean(0) > 0.5
    hi = np.argmax(~alive | (gm > 0.3 * Lc)) if (~alive | (gm > 0.3 * Lc)).any() else L - 1
    sel = (t >= float(fitlo)) & (np.arange(L) <= hi)
    for tt in (10, 20, 40, 80, 160, 320):
        k = tt // 2 - 1
        if k < L and alive[k]:
            print(f"  l({tt}) = {gm[k]:.2f}  l/sqrt(Dt) = {gm[k] / np.sqrt(tt):.3f}")
    if sel.sum() > 5:
        e, lnA = np.polyfit(np.log(t[sel]), np.log(gm[sel]), 1)
        A2 = np.mean(gm[sel] / np.sqrt(t[sel]))
        print(f"SUMMARY exponent: l ~ t^{e:.3f} over t in [{t[sel][0]:.0f}, {t[sel][-1]:.0f}] (l from {gm[sel][0]:.1f} to {gm[sel][-1]:.1f}, "
              f"limit 0.3*{Lc:.0f} = {0.3 * Lc:.1f}); A = l/sqrt(Dt) at fixed 0.5: {A2:.3f}")
    print(f"# wall time {time.time() - t0:.0f} s")


def special_job(size="16"):
    need_pred()
    R = int(size)
    g = grid("sphere", R)
    ta, R2 = tau("sphere", R), R**2
    print(f"# special R={R}: tau_sw={ta:.2f}, R^2/D={R2}")
    for u0 in (0.3, 0.1, 0.01, 0.001):  # offset circle: cos(theta0) = u0
        r = evolve(g, cap_ic(g, np.arccos(u0)), tmax=60 * R2)
        print(f"offset circle u0={u0}: T1 {r['T1']:.1f} = {r['T1'] / ta:.2f} tau_sw; sharp-interface cap law {cap_T1(R, u0):.1f}; ratio {r['T1'] / cap_T1(R, u0):.3f}")
    for eps in (0.05, 0.3):  # symmetric wavy circle (mean(a) = 0 by symmetry), m = 2
        tm = 30 * ta
        r = evolve(g, wavy_ic(g, eps, 2), tmax=tm, stop=0.01, probe=lambda c: np.sqrt((c[2] ** 2).sum()))
        t, p = r["t"], r["probe"]
        sel = (t >= 0.3 * R2) & (t <= 1.5 * R2)
        s = np.polyfit(t[sel], np.log(p[sel]), 1)[0]
        k = np.searchsorted(t, 1.5 * R2)
        print(f"wavy m=2 eps={eps}: m=2 amplitude decay rate {-s * R2:.3f} D/R^2 (linear theory 3); amplitude(1.5 R^2/D)/amplitude(0.3 R^2/D) = "
              f"{p[k] / p[np.searchsorted(t, 0.3 * R2)]:.4f}; minority at {t[-1]:.0f} ({t[-1] / ta:.1f} tau_sw) = {r['mino'][-1]:.4f}; T1 = {r['T1']}")


if __name__ == "__main__":
    job = sys.argv[1]
    kw = dict(a.split("=") for a in sys.argv[2:])
    {"check": check, "random": random_job, "special": special_job}[job](**kw)
    if job == "check":
        print("self-check passed")
