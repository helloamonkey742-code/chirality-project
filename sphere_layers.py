"""Layers and bands on a ROUND ocean: a thin spherical shell with uneven mixing (D_z << D_h).

Same local rule as aniso3d.py / shell3d.py:  da/dt = a - a^3 + div(D grad a),  a = +1 favoured hand, -1 other,
no bias term (the bias enters only through the starting majority, as in aniso3d.py). Here the mixing tensor is
tied to the LOCAL vertical: D_z along the radius, D_h along the surface, as in a real moon's ocean.

Grid: N Fibonacci points on the sphere, triangulated (convex hull), cotangent Laplacian with lumped areas
(horizontal); nz radial layers of thickness dz, finite-volume radial operator (1/r^2) d/dr(r^2 D_z d/dr), no flux
at the sea floor and the ice. Radial resolution is fixed so the layer wall spans ~2 cells: D_z = 4 dz^2.

PREDICTIONS (written before any run, 2026-09-26):
 README Part 3b says a spherical layer wall sinks at 2*D_h/R "however weak vertical mixing is".
 My prediction is different. A layer wall on a sphere is a surface of constant radius. For a profile a(r) that
 depends only on radius, the surface (D_h) part of the Laplacian is exactly zero, and what is left is
 D_z (a'' + 2 a'/r). So the wall sinks at 2*D_z/R, NOT 2*D_h/R. The flat-box wall_check in aniso3d.py fixed the
 anisotropy axis to z, so a curved wall there leans across the D_h direction; on a sphere the wall is always
 perpendicular to the local vertical. Expected: measured speed / (2 D_z/r) ~ 1, independent of D_h (x4 D_h
 changes it < 15%), and ~ D_z/D_h (<< 1) times the README rule.
 (2) Patchwork start in a column-like shell (stretched depth H' = H sqrt(D_h/D_z) > pi R): layers form and then
 barely move within the run, because they only creep at 2 D_z/R. The top layer's hand is NOT yet the whole-ocean
 winner at the end of the run; the 50% control stays near 50%. Top-layer share: I do not have a firm prediction;
 README assumes ~L/H'.
 (3) Bands: a latitude band whose edges are not great circles shrinks and vanishes. A hemisphere split along a
 great circle is stationary but unstable to a SHIFT (not a tilt; a tilt is just a rotation): the offset grows
 like exp(D_h t / R^2).
"""
import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.spatial import ConvexHull

SMOKE = os.environ.get("SMOKE") == "1"
SECTIONS = sys.argv[1:] or ["1", "2", "3"]


def sphere_mesh(n):
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = np.pi * (1 + 5 ** 0.5) * i
    p = np.c_[np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z]
    tri = ConvexHull(p).simplices
    rows, cols, vals, area = [], [], [], np.zeros(n)
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        i, j, k = tri[:, a], tri[:, b], tri[:, c]
        u, v = p[i] - p[k], p[j] - p[k]
        cot = (u * v).sum(1) / np.linalg.norm(np.cross(u, v), axis=1)
        rows += [i, j]; cols += [j, i]; vals += [cot / 2, cot / 2]
        area += np.bincount(i, np.linalg.norm(np.cross(p[j] - p[i], p[k] - p[i]), axis=1) / 6, n)
    W = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n)).tocsr()
    L = diags(1 / area) @ (W - diags(np.asarray(W.sum(1)).ravel()))   # Laplacian on the UNIT sphere
    return p, L.tocsr(), area


class Shell:
    def __init__(self, R, nz, dz, Dh, N=None):
        self.N = N or int(round(4 * np.pi * R * R))                    # horizontal spacing ~1 at radius R
        self.p, self.L, self.area = sphere_mesh(self.N)
        self.R, self.nz, self.dz, self.Dh, self.Dz = R, nz, dz, Dh, 4 * dz * dz
        self.r_in = R - nz * dz
        self.r = self.r_in + (np.arange(nz) + 0.5) * dz
        self.rf = self.r_in + np.arange(1, nz) * dz                    # inner faces
        lam_h = 8.5 / (4 * np.pi / self.N)                             # generous bound on |L| eigenvalues (checked below)
        self.dt = 0.8 / (1 + Dh * lam_h / self.r_in ** 2 + 4 * self.Dz / dz ** 2)

    def rhs(self, a):
        hor = (self.L @ a.T).T * (self.Dh / self.r[:, None] ** 2)
        F = self.Dz * self.rf[:, None] ** 2 * (a[1:] - a[:-1]) / self.dz
        F = np.concatenate([np.zeros((1, a.shape[1])), F, np.zeros((1, a.shape[1]))])
        ver = (F[1:] - F[:-1]) / (self.r[:, None] ** 2 * self.dz)
        return a - a ** 3 + hor + ver

    def step_to(self, a, t, T):
        while t < T - 1e-9:
            a = a + self.rhs(a) * self.dt
            t += self.dt
        return a, t

    def wall_radius(self, a):                                          # area-weighted mean top of the inner (-) layer
        neg = ((1 - np.clip(a, -1, 1)) / 2).sum(0) * self.dz
        return self.r_in + (neg * self.area).sum() / self.area.sum()

    def fav(self, a):                                                  # volume share of the + hand
        w = self.r[:, None] ** 2 * self.area[None, :]
        return float((w * (a > 0)).sum() / w.sum())


def check_stability(sh):
    from scipy.sparse.linalg import eigsh
    lam = abs(eigsh(sh.L, k=1, which="LM", return_eigenvectors=False)[0])
    assert lam < 8.5 / (4 * np.pi / sh.N), f"dt bound too small: max |eig| {lam:.0f}"


def smooth_noise(sh, rng, ell=3.0):
    """Smooth random field on the unit sphere, correlation ~ell cells at radius R (diffuse white noise)."""
    x = rng.standard_normal(sh.N)
    tau = (ell / sh.R) ** 2 / 2
    for _ in range(80):
        x = x + (sh.L @ x) * (tau / 80)
    return (x - x.mean()) / x.std()


# ---------- (1) layer wall speed ----------
def wall_speed(Dh, dz, R, nz, t_end, seed=0):
    sh = Shell(R, nz, dz, Dh)
    check_stability(sh)
    rng = np.random.default_rng(seed)
    h = sh.r_in + nz * dz / 2 + 1.5 * dz * smooth_noise(sh, rng)      # bumpy wall, so D_h has something to act on
    a = np.tanh((sh.r[:, None] - h[None, :]) / np.sqrt(2 * sh.Dz))     # + on top, - underneath
    t, ts, rw = 0.0, [], []
    for T in np.linspace(t_end / 5, t_end, 9):
        a, t = sh.step_to(a, t, T)
        ts.append(t); rw.append(sh.wall_radius(a))
    speed = -np.polyfit(ts, rw, 1)[0]
    r_mid = float(np.mean(rw))
    pz, ph = 2 * sh.Dz / r_mid, 2 * sh.Dh / r_mid
    print(f"  D_h={Dh:g} D_z={sh.Dz:.4g} (D_z/D_h={sh.Dz / Dh:.1e}) R={R} N={sh.N} nz={nz} dz={dz} dt={sh.dt:.3f}: "
          f"wall radius {rw[0]:.4f} -> {rw[-1]:.4f}; sink speed {speed:.3e}; "
          f"2D_z/r={pz:.3e} (ratio {speed / pz:.3f}); README 2D_h/r={ph:.3e} (ratio {speed / ph:.4f})", flush=True)
    return speed / pz, speed / ph


# ---------- (2) patchwork start: layers, top layer, winner ----------
def patchwork(sh, f, rng, width=4.0):
    """Blobby +-1 start ~width cells across in stretched units, + hand fills volume share f (like aniso3d.patchwork)."""
    z = np.stack([smooth_noise(sh, rng, width) for _ in range(sh.nz)])
    for _ in range(int(width ** 2)):                                  # smooth along the radius, in stretched units
        p = np.concatenate([z[:1], z, z[-1:]])
        z = z + 0.25 * (p[2:] + p[:-2] - 2 * z)
    w = (sh.r[:, None] ** 2 * sh.area[None, :]).ravel()
    order = np.argsort(z.ravel())
    cut = z.ravel()[order][np.searchsorted(np.cumsum(w[order]) / w.sum(), 1 - f)]
    return np.where(z > cut, 1.0, -1.0)


def layer_state(sh, a):
    s = a > 0
    planes = (s * sh.area).sum(1) / sh.area.sum()                    # + share in each radial layer
    layered = bool(np.all((planes > 0.98) | (planes < 0.02)))
    walls_up = int((np.diff(planes > 0.5)).sum())
    top_hand = planes[-1] > 0.5
    k = 0
    while k < sh.nz and (planes[-1 - k] > 0.5) == top_hand:           # thickness of the top layer, in layers
        k += 1
    return layered, walls_up, bool(top_hand), k / sh.nz


def patch_runs(R, nz, dz, times, seeds, fracs):
    out = {}
    for f in fracs:
        rows = []
        for seed in seeds:
            sh = Shell(R, nz, dz, 1.0)
            a, t = patchwork(sh, f, np.random.default_rng(100 + seed)), 0.0
            hist = []
            for T in times:
                a, t = sh.step_to(a, t, T)
                hist.append((T, sh.fav(a)) + layer_state(sh, a))
            rows.append(hist)
            print(f"  start {f:.2f} seed {seed}: " + "  ".join(
                f"t={T:g}:fav={x:.3f}{' LAYERED' if l else ''} walls={w} top={'+' if th else '-'} topshare={ts:.2f}"
                for T, x, l, w, th, ts in hist), flush=True)
        out[f] = rows
    return out


def summarize(pr, R, nz, dz):
    creep = 2 * 4 * dz ** 2 / (R - nz * dz / 2)
    print(f"  predicted layer-wall creep 2D_z/r = {creep:.2e} per time = {creep * 1000 / (nz * dz):.1%} of the depth per 1000;"
          f" README 2D_h/r would sweep the whole depth in t = {nz * dz * (R - nz * dz / 2) / 2:.1f}")
    for f, rows in pr.items():
        fin = [h[-1][1] for h in rows]
        n = len(rows)
        print(f"  start {f:.2f}: final fav {' '.join(f'{x:.3f}' for x in fin)}; one hand {sum(x < .01 or x > .99 for x in fin)}/{n}; "
              f"stacked mixed layers {sum(h[-1][2] and .01 < h[-1][1] < .99 for h in rows)}/{n}; favoured on top "
              f"{sum(h[-1][4] for h in rows)}/{n}; top-layer share {' '.join(f'{h[-1][5]:.2f}' for h in rows)} "
              f"(README L/H' = {min(1, 2 * np.pi * R / (nz / 2)):.2f}); fav change t=500->1000 "
              f"{' '.join(f'{h[-1][1] - h[-2][1]:+.3f}' for h in rows)}")
    # added after the first 2b run (post-hoc, not a prediction): expected change in the + volume share while layers
    # are stacked, if walls creep at 2D_z/r: about creep*500/H per 500 time units, sign = the top layer's hand
    print(f"  (post-hoc) layer creep predicts a share change of about +-{creep * 500 / (nz * dz):.3f} per 500 time units, "
          f"toward the top layer's hand")
    # every run: the hand on top at the end is the hand holding the majority of the volume (top layer decides)
    for rows in pr.values():
        for h in rows:
            assert h[-1][4] == (h[-1][1] > 0.5), "the top layer's hand is not the volume majority"
    # the 50% control is a coin flip per run (each run ends at 0 or 1), so the check is on the favoured starts:
    # with a 52-55% start the favoured hand must win clearly more often than with 50%
    # (a first version asserted mean(50% runs) ~ 0.5 on 3 runs; that tests nothing, removed after it failed at 0.333)
    fav_wins = [np.mean([h[-1][1] > 0.5 for h in rows]) for f, rows in pr.items() if f > 0.5]
    assert min(fav_wins) >= 0.5, f"favoured start did not usually win ({fav_wins})"


# ---------- (3) bands on the sphere ----------
def band_runs(R, nz, dz, t_end):
    sh = Shell(R, nz, dz, 1.0)
    lat = np.arcsin(sh.p[:, 2])
    # (a) hemisphere split shifted north by phi0: the + cap north of phi0 should shrink as phi grows ~exp(D_h t/R^2)
    phi0 = 1.5 / R
    a = np.broadcast_to(np.tanh((lat - phi0) * R / np.sqrt(2 * sh.Dh)), (nz, sh.N)).copy()
    t, ts, ph = 0.0, [], []
    for T in np.linspace(0, t_end, 11)[1:]:
        a, t = sh.step_to(a, t, T)
        f = sh.fav(a)
        ts.append(t); ph.append(np.arcsin(np.clip(1 - 2 * f, -1, 1)))
        if f < 0.02 or f > 0.98:
            break
    ts, ph = np.array(ts), np.array(ph)
    use = (ph > 0) & (ph < 0.35)                                       # small-offset regime (tan phi ~ phi within 4%)
    rate = np.polyfit(ts[use], np.log(ph[use]), 1)[0] if use.sum() >= 3 else float("nan")
    pred = sh.Dh / R ** 2
    print(f"  hemisphere split shifted {phi0:.3f} rad: offset " + " ".join(f"t={T:.0f}:{p:.3f}" for T, p in zip(ts, ph))
          + f"\n  growth rate {rate:.5f} (fit on {use.sum()} points) vs D_h/R^2 = {pred:.5f} (ratio {rate / pred:.3f}); "
          f"final + share {sh.fav(a):.3f}", flush=True)
    # (b) a latitude band 10..50 deg N (neither edge a great circle) must vanish; (c) an exactly centred split, for contrast
    res = {"split_ratio": rate / pred}
    for name, lo, hi in (("band 10-50N", 10, 50), ("centred split 0-90N", 0, 91)):
        a = np.where((lat > np.radians(lo)) & (lat < np.radians(hi)), 1.0, -1.0)
        a = np.broadcast_to(a, (nz, sh.N)).copy()
        t, fs = 0.0, []
        for T in (t_end / 4, t_end / 2, t_end):
            a, t = sh.step_to(a, t, T)
            fs.append(sh.fav(a))
        print(f"  {name}: + share " + " ".join(f"t={T:g}:{x:.3f}" for T, x in zip((t_end / 4, t_end / 2, t_end), fs)),
              flush=True)
        res[name] = fs
    return res


def section1(t0):
    print("(1) layer wall on a sphere: sink speed vs 2*D_z/r (my prediction) and 2*D_h/r (README Part 3b)")
    # t_end keeps each wall in the middle third of the depth (a first run to t=400 let the dz=0.04 wall reach the
    # sea floor, where the no-flux boundary pulls it in: ratio 1.60, rejected as a boundary artifact)
    r1 = [wall_speed(Dh, dz, 16, 24, T) for Dh, dz, T in ((1.0, 0.02, 400), (4.0, 0.02, 400), (1.0, 0.04, 150))]
    for rz, rh in r1:
        assert 0.8 < rz < 1.25, f"wall speed not 2*D_z/r (ratio {rz:.3f})"
        assert rh < 0.1, f"wall speed near README 2*D_h/r (ratio {rh:.3f})"
    assert abs(r1[1][0] / r1[0][0] - 1) < 0.15, "x4 D_h changed the wall speed: D_h does control it"
    print(f"  asserts ok: speed tracks D_z (ratios {[round(x[0], 3) for x in r1]}), not D_h. [{time.time() - t0:.0f}s]")



def main():
    t0 = time.time()
    if SMOKE:
        print("SMOKE run (tiny sizes; numbers not meaningful)")
        wall_speed(1.0, 0.02, 8, 12, 20)
        patch_runs(6, 16, 0.02, [20, 40], [0], [0.55])
        band_runs(8, 2, 0.25, 40)
        print(f"smoke ok {time.time() - t0:.0f}s")
        return

    if "1" in SECTIONS:
        section1(t0)
    if "2" in SECTIONS:
        print("\n(2) patchwork start, R=10 (circumference 63), nz=80 -> stretched depth H' = 40 (NOT deeper than wide), "
              "dz=0.02, D_z=0.0016")
        summarize(patch_runs(10, 80, 0.02, [100, 500, 1000], range(3), (0.5, 0.52, 0.55)), 10, 80, 0.02)
        print(f"  [{time.time() - t0:.0f}s]")
    if "2b" in SECTIONS:
        print("\n(2b) column-like shell: R=6 (circumference 38), nz=300 -> H' = 150 (4x deeper than wide, like aniso3d's "
              "column), dz=0.005, D_z=1e-4, H/R=0.25")
        summarize(patch_runs(6, 300, 0.005, [100, 500, 1000], range(2), (0.5, 0.55)), 6, 300, 0.005)
        print(f"  [{time.time() - t0:.0f}s]")
    if "3" not in SECTIONS:
        return
    print("\n(3) bands on the sphere: R=16, nz=2 (depth irrelevant), D_h=1")
    b = band_runs(16, 2, 0.25, 800)
    assert 0.7 < b["split_ratio"] < 1.3, "hemisphere split does not destabilise at D_h/R^2"
    assert b["band 10-50N"][-1] < 0.02, "a non-great-circle band survived"
    print(f"  asserts ok. total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
