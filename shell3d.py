"""Coarsening in a thin ocean shell instead of a cube.

Real icy-moon oceans are ~10-100x wider than deep. Model: 256 x 256 x 8 slab, periodic sideways,
closed top and bottom (no flux). Once patches are wider than the ocean is deep, the ocean acts 2D.
Checks: (1) patches still grow like sqrt(t); (2) the majority hand still takes over.
"""
import numpy as np

NX, NZ = 256, 8          # growth test; the majority test uses NX = 512 so many patches fit


def lap(a):
    out = sum(np.roll(a, s, ax) for ax in (1, 2) for s in (1, -1)) - 4 * a
    p = np.concatenate([a[..., :1], a, a[..., -1:]], axis=3)      # no-flux top/bottom
    return out + p[..., 2:] + p[..., :-2] - 2 * a


def sideways_length(a):
    s = np.sign(a)
    walls = sum((s != np.roll(s, 1, ax)).sum() for ax in (1, 2))
    return 2 * a.size / max(walls, 1)


def patchwork(f, rng, width=4, nx=NX):
    """Saturated +-1 patches ~width cells across; favoured hand fills fraction f."""
    kx, kz = np.fft.fftfreq(nx), np.fft.fftfreq(NZ)
    kk = kx[:, None, None]**2 + kx[None, :, None]**2 + kz[None, None, :]**2
    z = np.fft.ifftn(np.fft.fftn(rng.standard_normal((nx, nx, NZ))) * np.exp(-2 * (np.pi * width)**2 * kk)).real
    return np.where(z > np.quantile(z, 1 - f), 1.0, -1.0).astype(np.float32)[None]


def run(f, times, D=1.0, seed=1, nx=NX):
    rng = np.random.default_rng(seed)
    a = patchwork(f, rng, nx=nx)
    dt, t, out = 0.05, 0.0, []
    for T in times:
        while t < T:
            a = a + (a - a**3 + D * lap(a)) * dt
            t += dt
        out.append((T, sideways_length(a), float((a > 0).mean())))
    return out


def main():
    rows = run(0.5, [25, 50, 100, 200, 400])
    ts, ls = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
    late = ls > NZ                                                    # patches wider than the ocean is deep
    slope = np.polyfit(np.log(ts[late]), np.log(ls[late]), 1)[0]
    print("f=0.5: " + "  ".join(f"t={T}:l={l:.0f}" for T, l, _ in rows) + f"   exponent (l > depth) {slope:.2f}")
    assert abs(slope - 0.5) < 0.12, "thin-shell growth is not sqrt(t)"
    print("\nmajority test: 512 x 512 x 8, 4 seeds each, favoured fraction at t=1 -> t=200 (mean +- sd)")
    for f in (0.5, 0.52, 0.55):
        start, end = [], []
        for seed in range(4):
            rows = run(f, [1, 200], seed=seed, nx=512)
            start.append(rows[0][2]); end.append(rows[1][2])
        print(f"  f={f}: {np.mean(start):.3f} -> {np.mean(end):.3f} +- {np.std(end):.3f}   "
              f"(change {np.mean(end) - np.mean(start):+.3f})", flush=True)
        if f == 0.5:
            assert abs(np.mean(end) - 0.5) < 0.05, "control drifted: finite-size noise too large"

if __name__ == "__main__":
    main()
