"""h14b_advect: how much slower is wall removal under eddies much larger than the wall width? (open item 14 / H14,
journals/ija_draft.md:250-254; E86 at :187; Section 2.3 "Eddy diffusion" at :479; 2026-10-21.)

Question. advect2d.py measured wall removal at ONE eddy size above the wall, l/w = 3.6 (advect2d.log:49-51: of the 29% of
runs still walled at freeze-out, 14 / 25 / 50 / 80% were single-domain 10 / 20 / 40 / 80 time units later; Eq. 17 with
D_eff = 8 gives 16.9). The sibling h14_eddy.py (h14_eddy.log:60-77 for l/w, :81-175 for the break-even factors) asks which
uniform slowdown f of the sideways coarsening time would lose 5% of each world's {l > w} draws, and finds that l/w reaches
4-9 (median) and 35-178 (90th percentile) at U = 0.02 m/s. This script measures the slowdown at l/w = 10 and 30 at the
same D_eff = 8, plus the exponent p in slowdown ~ (l/w)^p, and tests the sibling's stress powers p = 1.1 and p = 2.

Model: advect2d.py unchanged in physics (see its docstring): 2D Allen-Cahn field, d(alpha) = (lam alpha - alpha^3 + g
+ D_mol lap(alpha)) dt + noise, lam = clip(-1 + gamma t, ., 1), gamma = 0.01, eps = 1e-4 per unit area, D_mol = 1, k = 1,
kicked alternating random-phase sine shear (cubic Lagrange kick), D_eddy = A^2 / (8 Tk) = 7, D_eff = 8, w = sqrt(8) = 2.83,
kappa = A / l = 1.469 as in advect2d's B arm. Eddy size l = L / (2 pi n) (advect2d's convention; ambiguity up to pi against
the sibling's Rhines-type l). Same-geometry scaling: Tk = kappa^2 l^2 / 56 rounded to a multiple of 0.2, A = sqrt(56 Tk).
Slowdown = KM median (S50) or 80th percentile (S80) of the time from freeze-out to single domain among runs that are
walled at freeze-out, divided by t17(L) = (L / 5.5)^2 / D_eff (Eq. 17 at ija_draft.md:493, A_coarsen = 5.5).
Cases A-K, thresholds P1-P12 and the decision rule are in h14b_advect_predictions.txt (asserted below).

PREDICTIONS: h14b_advect_predictions.txt (v1, 2026-10-21 01:07 PDT), sha256
d077b49b367a6a7160696d767a2a19dc97497974dc3e61611264144d83aaa869 (asserted below).

Changes made after the hash and before any results run (cost estimate; no outcome seen; thresholds untouched):
 (1) Kernel. A fused numba step replaces the numpy expression (identical Euler update; noise-off max |diff| vs the numpy
     path is asserted in self_check). Cases A and B draw numpy Gaussian noise exactly as advect2d does. All other cases draw
     in-kernel xorshift32 noise (uniform, unit variance, per-row state) because numpy noise is 60-80% of the step cost; a
     case B2 (= B with the xorshift noise, 192 runs) is added as the noise-distribution control (not in the predictions).
 (2) Run counts: case I has 18 runs (six jobs of 3) instead of 16; F, I, G are split into jobs of equal size and pooled.
 (3) Runs that are single-domain are dropped from the batch after the check at which they become single (their wall
     length is 0 from then on). This changes the random stream after freeze-out relative to advect2d, so even case A is a
     new realisation after the sweep, not a replay.
 (4) Pre-registered break-even f (h14_eddy.log) is the first ladder f (FL, h14_eddy.py:53) with >= 5% lost, so the largest
     ladder f known to be safe is the entry before it; the verdict table gives the pre-registered rule (S < f_be) and the
     stricter rule (S < previous ladder f).
 (5) Censored quantiles (not spelled out in the predictions): an undefined (not reached within the hold) S50/S80 makes
     P2, P3, P4, P5, P11 UNDECIDED unless the censoring floor already breaks the range; in P10 an undefined S80 bound is
     "unbounded" (never COVERED) and an undefined lower bound of S50 at I is replaced by the censoring floor.
Disclosure of what was seen before the hash: advect2d.py/.log, the h14_eddy files, coarsen_slow.py, and a scratch kernel
timing/equivalence pilot (no model statistic).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python h14b_advect.py \
       > h14b_advect.log      (4 worker processes; self-check first; "report" re-reads h14b_advect_raw.pkl; "smoke" = self-check)
"""
import hashlib
import os
import pickle
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numba as nb
import numpy as np
from scipy.ndimage import label

import kz_patches as K

HERE = os.path.dirname(os.path.abspath(__file__))
PRED_SHA = "d077b49b367a6a7160696d767a2a19dc97497974dc3e61611264144d83aaa869"
assert hashlib.sha256(open(os.path.join(HERE, "h14b_advect_predictions.txt"), "rb").read()).hexdigest() == PRED_SHA

GAMMA, EPS = K.GAMMA, K.EPS
D_MOL, D_EFF, DELTAS = 1.0, 8.0, (0.5, 1.0, 2.0)
W_WALL = float(np.sqrt(D_EFF))                       # sqrt(D_eff / k), k = 1
L_B = 64.0 / (2 * np.pi)                             # advect2d B arm eddy size
KAPPA2 = 8 * (D_EFF - D_MOL) * 4.0 / L_B ** 2        # kappa^2 of the B arm (Tk = 4): A^2 / l^2
AMIN, T17_A = 8.0, 5.5
PKL = os.path.join(HERE, "h14b_advect_raw.pkl")
GUARD_S = 40 * 60                                    # pre-registered wall-clock guard
NWORK = 4


def t17(Lbox):
    return (Lbox / T17_A) ** 2 / D_EFF


def make_flow(n, Lbox):
    l = Lbox / (2 * np.pi * n)
    Tk = round(KAPPA2 * l ** 2 / (8 * (D_EFF - D_MOL)) / 0.2) * 0.2
    return dict(n=n, Tk=float(Tk), A=float(np.sqrt(8 * (D_EFF - D_MOL) * Tk)), l=float(l))


# ----------------------------------------------------------------------------------------------- kernels
@nb.njit
def rd_g(a, b, nz, lam, g, Dl, dt, nzamp):
    """One explicit Euler step of reaction + diffusion + noise for a batch; noise array nz is supplied (Gaussian)."""
    R, N, _ = a.shape
    for r in range(R):
        gr = g[r]
        for i in range(N):
            ip = i + 1 if i + 1 < N else 0
            im = i - 1 if i > 0 else N - 1
            for j in range(N):
                jp = j + 1 if j + 1 < N else 0
                jm = j - 1 if j > 0 else N - 1
                x = a[r, i, j]
                lap = a[r, ip, j] + a[r, im, j] + a[r, i, jp] + a[r, i, jm] - 4 * x
                b[r, i, j] = x + (lam * x - x * x * x + gr + Dl * lap) * dt + nzamp * nz[r, i, j]


M32 = np.uint64(0xFFFFFFFF)


@nb.njit
def xs_next(s):
    s ^= (s << np.uint64(13)) & M32
    s ^= s >> np.uint64(17)
    s ^= (s << np.uint64(5)) & M32
    return s


@nb.njit
def rd_x(a, b, st, lam, g, Dl, dt, nzamp):
    """Same step; the noise is uniform, mean 0, variance 1 (xorshift32, one state per row, advanced along the row)."""
    R, N, _ = a.shape
    for r in range(R):
        gr = g[r]
        for i in range(N):
            ip = i + 1 if i + 1 < N else 0
            im = i - 1 if i > 0 else N - 1
            s = st[r, i]
            for j in range(N):
                jp = j + 1 if j + 1 < N else 0
                jm = j - 1 if j > 0 else N - 1
                s = xs_next(s)
                u = (np.float32(s >> np.uint64(8)) * np.float32(1.0 / 16777216.0) - np.float32(0.5)) * np.float32(3.4641016)
                x = a[r, i, j]
                lap = a[r, ip, j] + a[r, im, j] + a[r, i, jp] + a[r, i, jm] - 4 * x
                b[r, i, j] = x + (lam * x - x * x * x + gr + Dl * lap) * dt + nzamp * u
            st[r, i] = s


@nb.njit
def xs_fill(st, n):
    out = np.empty((n, st.shape[1]), np.float32)
    for k in range(n):
        for i in range(st.shape[1]):
            s = xs_next(st[0, i])
            st[0, i] = s
            out[k, i] = (np.float32(s >> np.uint64(8)) * np.float32(1.0 / 16777216.0) - np.float32(0.5)) * np.float32(3.4641016)
    return out


# ----------------------------------------------------------------------------------------------- advect2d pieces
def kick(a, s, h, axis):
    """advect2d.kick: shear displacement, cubic Lagrange, circulant (conserves the mean)."""
    N = a.shape[-1]
    b = a if axis == 2 else a.swapaxes(1, 2)
    sc = s / h
    m = np.floor(sc)
    u = (1 - (sc - m)).astype(np.float32)[..., None]
    base = np.arange(N)[None, None, :] - m.astype(np.intp)[:, :, None]
    ws = ((-2, -u * (u - 1) * (u - 2) / 6), (-1, (u + 1) * (u - 1) * (u - 2) / 2),
          (0, -(u + 1) * u * (u - 2) / 2), (1, (u + 1) * u * (u - 1) / 6))
    out = sum(w * np.take_along_axis(b, (base + off) % N, axis=2) for off, w in ws)
    return np.ascontiguousarray(out if axis == 2 else out.swapaxes(1, 2))


def shear(rng, flow, runs, N, h, axis, Lbox):
    phi = rng.uniform(0, 2 * np.pi, runs)[:, None]
    coord = (np.arange(N) + 0.5) * h
    return flow["A"] * np.sin(2 * np.pi * flow["n"] * coord[None, :] / Lbox + phi), axis


def coarse(a, f):
    R, N = a.shape[0], a.shape[1]
    return a if f == 1 else a.reshape(R, N // f, f, N // f, f).mean((2, 4))


def n_domains(s, amin=AMIN):
    """advect2d.n_domains: periodic 4-connected domains of both signs with area >= amin cells."""
    tot = 0
    for val in (True, False):
        lab, nl = label(s == val)
        if nl == 0:
            continue
        par = list(range(nl + 1))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for p, q in list(zip(lab[0], lab[-1])) + list(zip(lab[:, 0], lab[:, -1])):
            if p and q:
                par[find(p)] = find(q)
        area = np.bincount(lab.ravel(), minlength=nl + 1)
        roots = {}
        for i in range(1, nl + 1):
            r = find(i)
            roots[r] = roots.get(r, 0) + area[i]
        tot += sum(v >= amin for v in roots.values())
    return tot


def wall_len(s, Lbox):
    return ((s != np.roll(s, 1, 1)).sum((1, 2)) + (s != np.roll(s, 1, 2)).sum((1, 2))) / Lbox


# ----------------------------------------------------------------------------------------------- one case (one job)
def run_case(sp, deadline):
    """Sweep lam -1 -> 1 then hold. Runs that become single-domain are dropped from the batch. Returns per-run results."""
    t_start = time.time()
    Lbox, h, dt, D = sp["L"], sp["h"], sp["dt"], sp["D"]
    gamma = sp.get("gamma", GAMMA)
    N = int(round(Lbox / h))
    assert abs(N * h - Lbox) < 1e-9 and D * dt / h ** 2 <= 0.25 + 1e-9, "grid / Euler stability"
    rng = np.random.default_rng(sp["seed"])
    S = sp["S"]
    runs = S * len(DELTAS)
    flow = sp["flow"]
    g = np.repeat([K.g_for(d, Lbox ** 2) for d in DELTAS], S).astype(np.float32)
    a = np.zeros((runs, N, N), np.float32)
    b = np.empty_like(a)
    nz = np.empty_like(a) if sp["noise"] == "g" else None
    st = rng.integers(1, 2 ** 32 - 1, (runs, N), dtype=np.uint64) if sp["noise"] == "x" else None
    nzamp, dtf, Dl = np.float32(np.sqrt(EPS / h ** 2 * dt)), np.float32(dt), np.float32(D / h ** 2)
    nsweep = int(round(2 / gamma / dt))
    ke = int(round(flow["Tk"] / dt)) if flow else 0
    if flow:
        assert abs(ke * dt - flow["Tk"]) < 1e-6
    f = max(1, int(round(1 / h)))
    cstep = int(round(sp["cad"] / dt))
    nh = int(sp["hold"] / (cstep * dt) + 1e-9)               # number of post-freeze-out cadence checks
    idx = np.arange(runs)
    t_single = np.full(runs, np.inf)
    tobs = np.zeros(runs)
    tgrid, wall_hist = [0.0], []
    res = dict(spec=sp, runs=runs, deltas=np.repeat(DELTAS, S), cut=False)
    kn, j_done, done = 0, 0, False
    for i in range(nsweep + nh * cstep):
        if flow and i % ke == 0:
            sh, ax = shear(rng, flow, len(idx), N, h, 2 if kn % 2 == 0 else 1, Lbox)
            a = kick(a, sh, h, ax)
            b = np.empty_like(a) if b.shape != a.shape else b
            kn += 1
        lam = np.float32(min(1.0, -1.0 + gamma * dt * i))
        if sp["noise"] == "g":
            rng.standard_normal(out=nz, dtype=np.float32)
            rd_g(a, b, nz, lam, g, Dl, dtf, nzamp)
        else:
            rd_x(a, b, st, lam, g, Dl, dtf, nzamp)
        a, b = b, a
        if i < nsweep - 1 or (i - (nsweep - 1)) % cstep:
            continue
        j = (i - (nsweep - 1)) // cstep
        c = coarse(a, f) > 0
        one = c.all((1, 2)) | (~c).all((1, 2))
        wl = np.zeros(runs)
        wl[idx] = wall_len(c, Lbox)
        if j == 0:
            res["wall0"] = wl.copy()
            res["ndom0"] = np.array([n_domains(x) if sp.get("ndom", True) else 0 for x in c])
            res["walled_at_t0"] = np.zeros(runs, bool)
            res["walled_at_t0"][idx[~one]] = True
        else:
            tgrid.append(j * cstep * dt)
            wall_hist.append(wl)
        t_single[idx[one]] = j * cstep * dt
        tobs[idx[one]] = j * cstep * dt
        tobs[idx[~one]] = j * cstep * dt
        keep = ~one
        if not keep.all():
            a, g, idx = np.ascontiguousarray(a[keep]), g[keep], idx[keep]
            b = np.empty_like(a)
            if nz is not None:
                nz = np.empty_like(a)
            if st is not None:
                st = np.ascontiguousarray(st[keep])
        if len(idx) == 0:
            break
        if time.time() > deadline:
            res["cut"] = True
            break
    W = np.zeros((len(tgrid), runs))
    W[0] = res["wall0"]
    for jj, wl in enumerate(wall_hist):
        W[jj + 1] = wl
    res.update(t_single=t_single, tobs=tobs, tgrid=np.array(tgrid), W=W, active_end=len(idx), elapsed=time.time() - t_start)
    assert np.isfinite(a).all(), "NaN/inf in field"
    return res


def job(args):
    name, sp, deadline = args
    return name, run_case(sp, deadline)


# ----------------------------------------------------------------------------------------------- case matrix
def cases():
    """name -> spec. Seeds: A uses advect2d's B-arm seed 702."""
    def sp(L, n, h, dt, S, hold_t17, noise, seed, D=D_MOL, hold=None, flowless=False):
        fl = None if flowless else make_flow(n, L)
        return dict(L=float(L), n=n, h=h, dt=dt, S=S, D=D, flow=fl, noise=noise, seed=seed,
                    hold=hold if hold is not None else float(np.ceil(hold_t17 * t17(L) / 2) * 2),
                    cad=max(2.0, round(t17(L) / 60 / 0.2) * 0.2))
    c = {}
    c["A"] = sp(64, 1, 0.5, 0.05, 64, 0, "g", 702, hold=160.0)
    c["B"] = sp(64, 1, 1.0, 0.2, 64, 12, "g", 1702)
    c["B2"] = sp(64, 1, 1.0, 0.2, 64, 12, "x", 2702)
    c["C"] = sp(128, 2, 1.0, 0.2, 32, 12, "x", 1703)
    c["D"] = sp(256, 4, 1.0, 0.2, 8, 8, "x", 1704)
    c["E"] = sp(384, 6, 1.0, 0.2, 3, 6, "x", 1705)
    for k in range(3):
        c[f"F#{k}"] = sp(178, 1, 0.5, 0.05, 4, 8, "x", 1710 + k)
    c["G"] = sp(178, 1, 1.0, 0.2, 32, 8, "x", 1720)
    c["H"] = sp(356, 2, 1.0, 0.2, 4, 6, "x", 1721)
    for k in range(6):
        c[f"I#{k}"] = sp(533, 1, 1.0, 0.2, 1, 8, "x", 1730 + k)
    c["J"] = sp(64, 1, 1.0, 0.2, 32, 16, "x", 1740, flowless=True, D=D_MOL)
    c["K"] = sp(178, 1, 1.0, 0.2, 8, 16, "x", 1741, flowless=True, D=D_MOL)
    for k, v in c.items():
        v["name"] = k
    return c


ORDER = ["I#0", "I#1", "I#2", "I#3", "I#4", "I#5", "A", "F#0", "F#1", "F#2", "E", "H", "D", "G", "K", "C", "B", "B2", "J"]


# ----------------------------------------------------------------------------------------------- statistics
def km_quantiles(T, ev, qs):
    """Kaplan-Meier times at which a fraction q has cleared (first t with S(t) <= 1 - q); inf if never reached."""
    if len(T) == 0:
        return [np.inf] * len(qs)
    o = np.lexsort((~ev, T))                    # events before censorings at equal times
    T, ev = T[o], ev[o]
    n = len(T)
    S, out, k = 1.0, [np.inf] * len(qs), 0
    for i in range(n):
        if ev[i]:
            S *= 1 - 1.0 / (n - i)
            for m, q in enumerate(qs):
                if out[m] == np.inf and S <= 1 - q + 1e-12:
                    out[m] = T[i]
    return out


def half_time(tgrid, W, idx):
    mw = W[:, idx].mean(1)
    ref = mw[0]
    if ref <= 0:
        return 0.0
    below = np.nonzero(mw <= 0.5 * ref)[0]
    if len(below) == 0:
        return np.inf
    j = below[0]
    if j == 0:
        return 0.0
    f = (mw[j - 1] - 0.5 * ref) / (mw[j - 1] - mw[j])
    return tgrid[j - 1] + f * (tgrid[j] - tgrid[j - 1])


def pool(parts):
    """Merge jobs of one case (same grids); truncate histories to the shortest."""
    nt = min(len(p["tgrid"]) for p in parts)
    r = dict(spec=parts[0]["spec"], cut=any(p["cut"] for p in parts), elapsed=sum(p["elapsed"] for p in parts))
    for k in ("t_single", "tobs", "wall0", "ndom0", "walled_at_t0", "deltas"):
        r[k] = np.concatenate([p[k] for p in parts])
    r["tgrid"] = parts[0]["tgrid"][:nt]
    r["W"] = np.concatenate([p["W"][:nt] for p in parts], axis=1)
    r["runs"] = len(r["t_single"])
    return r


class Stat:
    """Point estimates and bootstrap (resampling runs) for one pooled case; times in units of t17(L)."""

    def __init__(self, c, B=2000, seed=0):
        self.c, self.t17 = c, t17(c["spec"]["L"])
        walled = c["walled_at_t0"]
        self.R, self.nw = c["runs"], int(walled.sum())
        self.pt = self._est(np.arange(self.R))
        rng = np.random.default_rng(seed)
        self.bs = np.array([self._est(rng.integers(0, self.R, self.R)) for _ in range(B)])

    def _est(self, idx):
        c = self.c
        w = c["walled_at_t0"][idx]
        T = c["tobs"][idx][w]
        ev = np.isfinite(c["t_single"][idx][w])
        q50, q80 = km_quantiles(T, ev, (0.5, 0.8)) if w.sum() >= 2 else (np.inf, np.inf)
        return np.array([q50 / self.t17, q80 / self.t17, half_time(c["tgrid"], c["W"], idx) / self.t17, w.mean(),
                         c["wall0"][idx].mean(), (np.isfinite(c["t_single"][idx]) & (c["t_single"][idx] <= 80)).mean()])

    def cens_floor(self):
        """Longest observation time among walled runs / t17: any undefined quantile exceeds this."""
        w = self.c["walled_at_t0"]
        return (self.c["tobs"][w].max() / self.t17) if w.any() else np.nan

    def get(self, k):
        """2.5, 16, 84, 97.5 percentiles of the bootstrap replicates (nearest rank; inf kept as inf)."""
        v = np.sort(self.bs[:, k])
        return np.array([v[min(len(v) - 1, int(np.ceil(q / 100 * (len(v) - 1))))] for q in (2.5, 16, 84, 97.5)])


def fmt(x, nd=2):
    return "undefined" if not np.isfinite(x) else f"{x:.{nd}f}"


def line(st, k, name):
    pt, pc = st.pt[k], st.get(k)
    s = f"{fmt(pt)}" if np.isfinite(pt) else f"undefined (censored; > {st.cens_floor():.2f})"
    return f"{name} {s}  [68% {fmt(pc[1])}-{fmt(pc[2])}; 95% {fmt(pc[0])}-{fmt(pc[3])}]"


def fit_p(lw, vals):
    x, y = np.log(lw), np.log(vals)
    if not np.isfinite(y).all():
        return np.nan
    xm = x - x.mean()
    return float((xm * (y - y.mean())).sum() / (xm ** 2).sum())


# ----------------------------------------------------------------------------------------------- sibling log parsing
def parse_sibling():
    """l/w median and p90 (section 4) and break-even f per world, with line numbers in h14_eddy.log."""
    rows, cur = {}, None
    p4 = re.compile(r"^(\S[\w ]*?)\s+U 0\.02 n\s+\d+ >3\.6:\s+[\d.]+%\s+med\s+([\d.]+)\s+p90\s+([\d.]+) \| U 1 n\s+\d+ >3\.6:\s+[\d.]+%\s+med\s+([\d.]+)\s+p90\s+([\d.]+)")
    p5 = re.compile(r"^(\w[\w ]*?): non-monotone draws")
    pu = re.compile(r"^\s+l at (U 0\.02|U 1|lit hi)\s+subset")
    pb = re.compile(r"break-even f \(5% of subset lost\): main no lid (\S+), main final ([^;]+);")
    lab = None
    for ln, txt in enumerate(open(os.path.join(HERE, "h14_eddy.log")).read().splitlines(), 1):
        m = p4.match(txt)
        if m:
            rows.setdefault(m.group(1), {})["lw"] = dict(zip(("U0.02", "U1"), ((float(m.group(2)), float(m.group(3))), (float(m.group(4)), float(m.group(5))))))
            rows[m.group(1)]["lw_line"] = ln
        m = p5.match(txt)
        if m:
            cur = m.group(1)
            continue
        m = pu.match(txt)
        if m:
            lab = m.group(1).replace(" ", "")
            continue
        m = pb.search(txt)
        if m and cur and lab:
            f = [None if v.strip() == "None" else float(v) for v in (m.group(1), m.group(2))]
            rows.setdefault(cur, {}).setdefault("f5", {})[lab] = (f, ln)
    return rows


FL = (1.0, 2.0, 3.0, 4.0, 6.0, 10.0, 14.0, 20.0, 30.0, 50.0, 100.0, 200.0, 500.0, 1e3, 1e4, 1e6, 1e10)


def prev_ladder(f):
    return FL[FL.index(f) - 1]


# ----------------------------------------------------------------------------------------------- self-check
def self_check():
    t0 = time.time()
    rng = np.random.default_rng(5)
    # (a) D_eddy = A^2 / (8 Tk) = D_eff - D_mol for every flow used
    for sp in cases().values():
        if sp["flow"]:
            fl = sp["flow"]
            assert abs(fl["A"] ** 2 / (8 * fl["Tk"]) - (D_EFF - D_MOL)) < 1e-9, fl
    assert abs(make_flow(1, 64)["Tk"] - 4.0) < 1e-9 and abs(make_flow(1, 64)["A"] - np.sqrt(224)) < 1e-9
    # (b) noise-off numba step == advect2d's numpy expression (with a kick)
    N, R, steps = 96, 3, 60
    a0 = (rng.standard_normal((R, N, N)) * 0.3).astype(np.float32)
    fl = dict(n=1, Tk=1.0, A=10.0)
    g = np.array([0.0, 1e-3, -1e-3], np.float32)
    h, dt, D = 0.5, 0.05, 1.0
    f32 = np.float32
    ref, c1, c2 = a0.copy(), a0.copy(), a0.copy()
    b1, b2 = np.empty_like(a0), np.empty_like(a0)
    nz0 = np.zeros_like(a0)
    st = rng.integers(1, 2 ** 32 - 1, (R, N), dtype=np.uint64)
    r2 = np.random.default_rng(9)
    for i in range(steps):
        if i % 20 == 0:
            sh, ax = shear(r2, fl, R, N, h, 2 if (i // 20) % 2 == 0 else 1, N * h)
            ref, c1, c2 = kick(ref, sh, h, ax), kick(c1, sh, h, ax), kick(c2, sh, h, ax)
        lam = f32(min(1.0, -1.0 + 0.01 * dt * i * 40))
        ref = ref + (lam * ref - ref ** 3 + g.reshape(R, 1, 1) + f32(D / h ** 2) * K.lap2(ref)) * f32(dt)
        rd_g(c1, b1, nz0, lam, g, f32(D / h ** 2), f32(dt), f32(0.0))
        c1, b1 = b1, c1
        rd_x(c2, b2, st, lam, g, f32(D / h ** 2), f32(dt), f32(0.0))
        c2, b2 = b2, c2
    d1, d2 = float(np.abs(ref - c1).max()), float(np.abs(ref - c2).max())
    assert d1 < 1e-4 and d2 < 1e-4, (d1, d2)
    # (c) xorshift noise: mean 0, variance 1, no lag-1 correlation
    o = xs_fill(rng.integers(1, 2 ** 32 - 1, (1, 4096), dtype=np.uint64), 400)
    assert abs(o.mean()) < 0.01 and abs(o.var() - 1) < 0.02, (o.mean(), o.var())
    assert abs((o[1:] * o[:-1]).mean()) < 0.01 and abs((o[:, 1:] * o[:, :-1]).mean()) < 0.01
    # (d) Kaplan-Meier and half-time on known inputs
    assert km_quantiles(np.array([1.0, 2.0, 3.0]), np.array([True] * 3), (0.5, 0.8)) == [2.0, 3.0]
    assert km_quantiles(np.array([1.0, 2.0, 3.0]), np.array([True, False, False]), (0.5,)) == [np.inf]
    tg, W = np.array([0.0, 1.0, 2.0]), np.array([[4.0, 4.0], [2.0, 2.0], [0.0, 0.0]])
    assert abs(half_time(tg, W, np.arange(2)) - 1.0) < 1e-12
    # (e) a tiny run: shapes, drop logic, history consistency (both noise kernels)
    for nzk in ("g", "x"):
        sp = dict(L=20.0, n=1, h=1.0, dt=0.2, S=2, D=1.0, flow=dict(n=1, Tk=1.0, A=3.0, l=3.0), noise=nzk, seed=3, hold=10.0,
                  cad=2.0, gamma=0.5, name="smoke")
        r = run_case(sp, time.time() + 600)
        assert r["W"].shape == (len(r["tgrid"]), 6) and np.isfinite(r["W"]).all()
        assert (np.isfinite(r["t_single"]) | (r["tobs"] > 0) | (r["tobs"] == 0)).all()
        assert (r["W"][0] == r["wall0"]).all() and r["walled_at_t0"].sum() == (r["t_single"] != 0).sum()
        cleared = np.isfinite(r["t_single"])
        assert (r["W"][-1][cleared & (r["t_single"] < r["tgrid"][-1])] == 0).all()
    # (f) sibling log parses
    sib = parse_sibling()
    assert len(sib) == 15 and all("lw" in v and "f5" in v for v in sib.values()), sorted(sib)
    print(f"self-check passed ({time.time() - t0:.0f}s): D_eddy = 7 for all flows; numba vs numpy noise-off max |diff| "
          f"{d1:.1e} (Gaussian kernel) / {d2:.1e} (xorshift kernel); xorshift mean {o.mean():+.4f} var {o.var():.4f}; "
          f"KM/half-time known cases; tiny runs consistent; 15 worlds parsed from h14_eddy.log", flush=True)


# ----------------------------------------------------------------------------------------------- report
def report(res, elapsed):
    parts = {}
    for k, v in res.items():
        parts.setdefault(k.split("#")[0], []).append(v)
    P = {k: pool(v) for k, v in parts.items()}
    st = {k: Stat(P[k], seed=sum(map(ord, k))) for k in P}
    lw = {k: (P[k]["spec"]["L"] / (2 * np.pi * P[k]["spec"]["n"])) / W_WALL for k in P if P[k]["spec"]["flow"]}
    print("\nCASES  (grid cells N, box L, eddy count n, l/w, Tk, A, t17(L), runs, walled at freeze-out, hold [t.u.], cadence, elapsed [s])")
    for k in sorted(P):
        sp = P[k]["spec"]
        fl = sp["flow"]
        print(f"  {k:3s} L={sp['L']:g} N={int(round(sp['L'] / sp['h']))} h={sp['h']:g} dt={sp['dt']:g} noise={sp['noise']} "
              + (f"n={sp['n']} l/w={lw[k]:.2f} Tk={fl['Tk']:g} A={fl['A']:.2f} " if fl else "no flow, D=1 ")
              + f"t17={t17(sp['L']):.1f} runs={P[k]['runs']} walled={st[k].nw} ({st[k].pt[3]:.3f}) hold={sp['hold']:g} cad={sp['cad']:g} "
              f"cens_at_end={int(np.isinf(P[k]['t_single'][P[k]['walled_at_t0']]).sum())} cut={P[k]['cut']} cpu={P[k]['elapsed']:.0f}")
    print("\nSLOWDOWN  (time since freeze-out / t17(L); Kaplan-Meier over runs walled at freeze-out; bootstrap over runs, 2000x)")
    for k in sorted(P):
        s = st[k]
        print(f"  {k:3s} {line(s, 0, 'S50')} | {line(s, 1, 'S80')} | {line(s, 2, 'S_half')}")
        print(f"      walled fraction {s.pt[3]:.3f}  wall0 {s.pt[4]:.3f}  single by t_after=80: {s.pt[5]:.3f}  mean domains at t0 {P[k]['ndom0'].mean():.2f}")

    def ratio_ci(a, b, k=0):
        va, vb = st[a].bs[:, k], st[b].bs[:, k]
        ok = np.isfinite(va) & np.isfinite(vb)
        r = va[ok] / vb[ok]
        return (st[a].pt[k] / st[b].pt[k], np.percentile(r, [2.5, 97.5]) if len(r) > 100 else (np.nan, np.nan))

    # exponent
    ser = [k for k in ("B", "G", "I") if k in st]
    lws = np.array([lw[k] for k in ser])
    pfit = {}
    for kk, nm in ((0, "S50"), (1, "S80"), (2, "S_half")):
        vals = np.array([st[k].pt[kk] for k in ser])
        bsv = np.stack([st[k].bs[:, kk] for k in ser], 1)
        pb = np.array([fit_p(lws, v) for v in bsv])
        fin = np.isfinite(pb)
        pfit[nm] = (fit_p(lws, vals), np.percentile(pb[fin], [2.5, 16, 84, 97.5]) if fin.sum() > 100 else None, 1 - fin.mean())
    print("\nEXPONENT  p in slowdown ~ (l/w)^p over cases " + ", ".join(f"{k} (l/w {lw[k]:.1f})" for k in ser))
    for nm, (p, ci, un) in pfit.items():
        print(f"  {nm:6s} p = {fmt(p)}  " + (f"[68% {ci[1]:.2f}..{ci[2]:.2f}; 95% {ci[0]:.2f}..{ci[3]:.2f}]" if ci is not None else "(bootstrap mostly undefined)")
              + f"  share of bootstrap replicates undefined {un:.3f}")
    for p in (1.1, 2.0):
        print(f"  sibling stress p = {p}: slowdown at 3.6 / 10 / 30 would be {3.6 ** p:.1f} / {10 ** p:.1f} / {30 ** p:.0f} "
              f"(p = 1.1 is calibrated so 3.6^1.08 = 4)")
    # n dependence
    print("\nN-DEPENDENCE at l/w = 3.6 (same flow, bigger box): S50 / S80 / S_half by n; and at l/w = 10")
    for ks in (("B", "C", "D", "E"), ("G", "H")):
        for k in ks:
            if k in st:
                print(f"  {k} n={P[k]['spec']['n']} L={P[k]['spec']['L']:g} runs={P[k]['runs']} walled={st[k].nw}  {line(st[k], 0, 'S50')} | {line(st[k], 1, 'S80')} | {line(st[k], 2, 'S_half')}")
    if "D" in st and "B" in st:
        r, ci = ratio_ci("D", "B")
        print(f"  S50(D)/S50(B) = {fmt(r)}  [95% {fmt(ci[0])}-{fmt(ci[1])}]")

    # predictions
    print("\nPREDICTIONS (thresholds from h14b_advect_predictions.txt)")
    out = []

    def verdict(n, ok, txt):
        ok = ok if isinstance(ok, str) else bool(ok)
        out.append((n, ok))
        print(f"P{n} {ok if isinstance(ok, str) else ('HELD' if ok else 'FAILED')}  {txt}")
    a = st["A"]
    a_ok = [0.20 <= a.pt[3] <= 0.38, 0.90 <= a.pt[5] <= 0.98, 0.6 <= a.pt[4] <= 1.3, 1.6 <= a.pt[0] <= 3.6, 3.0 <= a.pt[1] <= 7.5]
    verdict(1, all(a_ok), f"A: walled {a.pt[3]:.3f} (log 0.292), single by 80 {a.pt[5]:.3f} (log 0.943), wall0 {a.pt[4]:.2f} (log 0.92), "
                          f"S50 {fmt(a.pt[0])} (derived 2.4), S80 {fmt(a.pt[1])} (derived 4.7); per-clause {a_ok}")
    r2, ci2 = ratio_ci("B", "A")
    d_w = abs(st["B"].pt[3] - st["A"].pt[3])
    verdict(2, ((0.67 <= r2 <= 1.5) and d_w <= 0.12) if np.isfinite(r2) else "UNDECIDED", f"S50(B)/S50(A) = {fmt(r2)} [95% {fmt(ci2[0])}-{fmt(ci2[1])}], walled B {st['B'].pt[3]:.3f} vs A {a.pt[3]:.3f} (diff {d_w:.3f})")
    if "B2" in st:
        r2b, ci2b = ratio_ci("B2", "B")
        print(f"   (control, not pre-registered) S50(B2)/S50(B) = {fmt(r2b)} [95% {fmt(ci2b[0])}-{fmt(ci2b[1])}], walled B2 {st['B2'].pt[3]:.3f} vs B {st['B'].pt[3]:.3f}")
    r3, ci3 = ratio_ci("G", "F")
    verdict(3, (0.67 <= r3 <= 1.5) if np.isfinite(r3) else "UNDECIDED", f"S50(G)/S50(F) = {fmt(r3)} [95% {fmt(ci3[0])}-{fmt(ci3[1])}]; S80(G)/S80(F) = {fmt(ratio_ci('G', 'F', 1)[0])}")
    f_ = st["F"]

    def inr(s_, k, lo, hi):
        """range check that treats a censored (undefined) quantile as 'above the censoring floor'."""
        v = s_.pt[k]
        if np.isfinite(v):
            return lo <= v <= hi
        return False if s_.cens_floor() >= hi else "UNDECIDED"
    r4 = [inr(f_, 0, 1.8, 6.0), inr(f_, 1, 3, 12)]
    verdict(4, False if False in r4 else ("UNDECIDED" if "UNDECIDED" in r4 else True), f"F (l/w {lw['F']:.1f}): S50 {fmt(f_.pt[0])} (guess 3.1), S80 {fmt(f_.pt[1])} (guess 5); G: S50 {fmt(st['G'].pt[0])}, S80 {fmt(st['G'].pt[1])}")
    i_ = st["I"]
    if np.isfinite(i_.pt[0]):
        s80_ok = (np.isfinite(i_.pt[1]) and 2.5 <= i_.pt[1] <= 16) or (not np.isfinite(i_.pt[1]) and i_.cens_floor() >= 2.5)
        verdict(5, (1.5 <= i_.pt[0] <= 8) and s80_ok, f"I (l/w {lw['I']:.1f}): S50 {fmt(i_.pt[0])} (guess 3.8), S80 {fmt(i_.pt[1]) if np.isfinite(i_.pt[1]) else 'censored, > ' + fmt(i_.cens_floor())} (guess 6)")
    else:
        verdict(5, False if i_.cens_floor() >= 8 else "UNDECIDED", f"I: S50 undefined (censored; > {i_.cens_floor():.2f} t17); pre-registered range [1.5, 8] cannot be tested above {i_.cens_floor():.1f}")
    p50, ci50, un50 = pfit["S50"]
    if np.isfinite(p50) and ci50 is not None:
        verdict(6, (-0.1 <= p50 <= 0.5) and ci50[3] < 1.1, f"p(S50) = {p50:.2f}, 97.5% bound {ci50[3]:.2f}; p = 1.1 {'excluded' if ci50[3] < 1.1 else 'NOT excluded'}; p(S80) = {fmt(pfit['S80'][0])}")
    else:
        verdict(6, "UNDECIDED", f"p undefined: some S50 is censored; lower bound on I gives p >= {fit_p(lws, np.array([st[k].pt[0] if np.isfinite(st[k].pt[0]) else st[k].cens_floor() for k in ser])):.2f}")
    verdict(7, st["G"].pt[3] >= 0.5 and st["F"].pt[3] >= 0.5 and st["I"].pt[3] >= 0.9,
            f"walled fractions G {st['G'].pt[3]:.2f}, F {st['F'].pt[3]:.2f}, I {st['I'].pt[3]:.2f}")
    if "D" in st:
        rr = st["D"].pt[0] / st["B"].pt[0]
        verdict(8, (rr <= 0.75) if np.isfinite(rr) else "UNDECIDED", f"S50(D, n=4)/S50(B, n=1) = {fmt(rr)}; D walled {st['D'].nw}/{st['D'].R}")
    mol = []
    for kf, km_ in (("B", "J"), ("G", "K")):
        sm = st[km_].pt[0]
        mol.append((kf, km_, sm, st[km_].cens_floor()))
    ok9 = []
    for kf, km_, sm, fl_ in mol:
        lo = sm if np.isfinite(sm) else fl_
        ok9.append(st[kf].pt[0] <= 0.7 * lo)
    s_j = st["J"].pt[0]
    j_ok = (7 <= s_j <= 20) if np.isfinite(s_j) else True
    verdict(9, j_ok and all(ok9), "S50_mol J (L=64) " + (fmt(s_j) if np.isfinite(s_j) else f"censored > {st['J'].cens_floor():.1f}") +
            ", K (L=178) " + (fmt(st['K'].pt[0]) if np.isfinite(st['K'].pt[0]) else f"censored > {st['K'].cens_floor():.1f}") +
            f"; S50 flow/mol bounds ok {ok9}")
    ok11 = []
    for k in ("B", "G", "I"):
        s50, sh = st[k].pt[0], st[k].pt[2]
        ok11.append(bool(np.isfinite(s50) and np.isfinite(sh) and 0.5 <= sh / s50 <= 2))
    verdict(11, all(ok11) if all(np.isfinite(st[k].pt[0]) and np.isfinite(st[k].pt[2]) for k in ("B", "G", "I")) else "UNDECIDED", "S_half/S50: " + " ".join(f"{k} {fmt(st[k].pt[2] / st[k].pt[0])}" for k in ("B", "G", "I")))
    cut = any(v["cut"] for v in res.values())
    verdict(12, (elapsed < GUARD_S) and not cut, f"total {elapsed:.0f} s, {NWORK} workers, cut by guard: {cut}")

    # worlds
    sib = parse_sibling()
    hi_cases = [k for k in ("B", "G", "I") if k in st]
    s80_hi = max(st[k].get(1)[3] for k in hi_cases)
    s50_lo_I = st["I"].get(0)[0]
    if not np.isfinite(s50_lo_I):
        s50_lo_I = st["I"].cens_floor()          # even the lower bound is censored: S50 exceeds the longest observation
    p_hi = max(0.0, pfit["S80"][1][3]) if (pfit["S80"][1] is not None and pfit["S80"][2] <= 0.025) else np.inf
    p_lo = max(0.0, pfit["S50"][1][0]) if pfit["S50"][1] is not None else 0.0
    print(f"\nWORLD VERDICTS (pre-registered rule P10). bound_hi(lw) = {fmt(s80_hi)} x max(1, lw/30)^{fmt(p_hi)}; "
          f"bound_lo(lw) = {fmt(s50_lo_I)} x max(1, lw/30)^{fmt(p_lo)}  (undefined = unbounded: quantile not reached within the hold)")
    print("  f_be = smaller non-None of the main no-lid / main final break-even f (h14_eddy.log line cited); strict = previous ladder f")
    ncov = {"U0.02": 0, "U1": 0}
    for wname in sib:
        s_ = sib[wname]
        for U in ("U0.02", "U1"):
            fs, ln = s_["f5"][U]
            fe = min(f for f in fs if f is not None)
            med, p90 = s_["lw"][U]
            bh = s80_hi * max(1.0, p90 / 30) ** p_hi if np.isfinite(s80_hi) else np.inf
            bl = s50_lo_I * max(1.0, med / 30) ** p_lo
            v = "COVERED" if bh < fe else ("NOT COVERED" if bl >= fe else "UNDECIDED")
            strict = "COVERED" if bh < prev_ladder(fe) else ("NOT COVERED" if bl >= fe else "UNDECIDED")
            ncov[U] += v == "COVERED"
            print(f"  {wname:18s} {U:6s} lw med {med:5.1f} p90 {p90:6.1f} (h14_eddy.log:{s_['lw_line']}) f_be {fe:g} (:{ln}) "
                  f"bound_hi {fmt(bh, 1)} bound_lo {fmt(bl, 1)} -> {v}; strict (f<{prev_ladder(fe):g}) -> {strict}")
    ok10 = ncov["U0.02"] == 15 and ncov["U1"] >= 13
    verdict(10, ok10, f"COVERED at U 0.02: {ncov['U0.02']}/15; at U 1: {ncov['U1']}/15 (needs 15 and >= 13)")
    held = [n for n, ok in out if ok is True]
    print("\nSUMMARY  held: " + ", ".join(f"P{n}" for n in sorted(held)) + " | failed: " + ", ".join(f"P{n}" for n, ok in sorted(out) if ok is False)
          + " | undecided: " + ", ".join(f"P{n}" for n, ok in sorted(out) if isinstance(ok, str)))
    print("DECISION (pre-registered): " + ("close item 14 as measured to l/w = 30" if (ok10 and any(n == 6 and ok is True for n, ok in out) and all(ok is True for n, ok in out if n in (2, 3)))
                                           else "narrow / keep open: see which of P2, P3, P6, P10 failed or are undecided"))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    t0 = time.time()
    if mode == "report":
        res = pickle.load(open(PKL, "rb"))
        report(res["res"], res["elapsed"])
        return
    self_check()
    if mode == "smoke":
        return
    cs = cases()
    print(f"\nL-scalings: kappa^2 = {KAPPA2:.3f}; w = {W_WALL:.3f}; t17: " + ", ".join(f"L={L}: {t17(L):.1f}" for L in (64, 128, 178, 256, 356, 384, 533)), flush=True)
    for k in ORDER:
        if k in cs:
            s = cs[k]
            print(f"  job {k}: L={s['L']:g} N={int(round(s['L'] / s['h']))} h={s['h']:g} dt={s['dt']:g} runs={s['S'] * 3} noise={s['noise']} hold={s['hold']:g} "
                  + (f"n={s['n']} l/w={s['flow']['l'] / W_WALL:.2f} Tk={s['flow']['Tk']:g} A={s['flow']['A']:.2f}" if s["flow"] else "no flow"), flush=True)
    deadline = t0 + GUARD_S
    res = {}
    with ProcessPoolExecutor(NWORK) as ex:
        futs = [ex.submit(job, (k, cs[k], deadline)) for k in ORDER]
        for fu in as_completed(futs):
            k, r = fu.result()
            res[k] = r
            print(f"  job done: {k} [{r['elapsed']:.0f}s] walled {int(r['walled_at_t0'].sum())}/{r['runs']} cleared {int(np.isfinite(r['t_single'][r['walled_at_t0']]).sum())} "
                  f"cut={r['cut']} (elapsed {time.time() - t0:.0f}s)", flush=True)
    elapsed = time.time() - t0
    pickle.dump(dict(res=res, elapsed=elapsed), open(PKL, "wb"))
    report(res, elapsed)


if __name__ == "__main__":
    main()
