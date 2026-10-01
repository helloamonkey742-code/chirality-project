"""Grid-resolution check for shell3d.py (asked for by A. Brandenburg, 2026-09-30).

Numerical scheme used by every spatial script (domains3d, coarsen3d, shell3d, aniso3d):
  space: second-order centred finite differences (7-point Laplacian), error O(dx^2)
  time:  explicit forward Euler, error O(dt); stable for D*dt/dx^2 <= 1/6 in 3D
A domain wall is only ~sqrt(2D) = 1.4 cells wide at dx = 1, D = 1, so walls are barely resolved
and 8 cells of depth is coarse. Test: SAME physical ocean (256 x 256 x 8 units, D = 1), grid spacing
dx = 1 (256x256x8, as in shell3d.py) vs dx = 1/2 (512x512x16), dt scaled with dx^2 so D*dt/dx^2 is
fixed; plus dx = 1 with the small dt, which isolates the time-step error from the grid error.
If the patch sizes and majority takeover agree across grids, the shell3d.py results are converged.
"""
import sys
import numpy as np

LX, LZ = 256, 8          # physical box, units where D = 1 and the reaction rate = 1


def lap(a, dx):
    out = sum(np.roll(a, s, ax) for ax in (0, 1) for s in (1, -1)) - 4 * a
    p = np.concatenate([a[..., :1], a, a[..., -1:]], axis=2)      # no-flux top/bottom
    return (out + p[..., 2:] + p[..., :-2] - 2 * a) / np.float32(dx * dx)


def patchwork(f, nx, nz, dx, rng, width=4.0):
    """Saturated +-1 patches ~width physical units across; favoured hand fills fraction f."""
    kx, kz = np.fft.fftfreq(nx, dx), np.fft.fftfreq(nz, dx)
    kk = kx[:, None, None]**2 + kx[None, :, None]**2 + kz[None, None, :]**2
    z = np.fft.ifftn(np.fft.fftn(rng.standard_normal((nx, nx, nz))) * np.exp(-2 * (np.pi * width)**2 * kk)).real
    return np.where(z > np.quantile(z, 1 - f), 1.0, -1.0).astype(np.float32)


def sideways_length(a, dx):
    s = a > 0
    walls = sum((s != np.roll(s, 1, ax)).sum() for ax in (0, 1))
    return 2 * a.size / max(walls, 1) * dx                        # physical units


def run(f, times, dx, dt, seed, lx=LX):
    nx, nz = round(lx / dx), round(LZ / dx)
    a = patchwork(f, nx, nz, dx, np.random.default_rng(seed))
    t, out, dtf = 0.0, [], np.float32(dt)
    for T in times:
        while t < T - 1e-9:
            a = a + (a - a**3 + lap(a, dx)) * dtf
            t += dt
        out.append((T, sideways_length(a, dx), float((a > 0).mean())))
    return out


def main(which):
    grids = {"dx=1 dt=0.05 (shell3d.py)": (1.0, 0.05), "dx=1 dt=0.0125": (1.0, 0.0125),
             "dx=1/2 dt=0.0125 (512x512x16)": (0.5, 0.0125)}
    if which == "growth":
        times = [25, 50, 100, 200, 400]
        res = {}
        for name, (dx, dt) in grids.items():
            rows = run(0.5, times, dx, dt, seed=1)
            ls = np.array([r[1] for r in rows])
            late = ls > LZ
            slope = np.polyfit(np.log(np.array(times)[late]), np.log(ls[late]), 1)[0]
            res[name] = (ls, slope)
            print(f"{name:32} " + " ".join(f"t={T}:l={l:.0f}" for T, l, _ in rows) + f"  exponent {slope:.2f}", flush=True)
        base, fine = res["dx=1 dt=0.05 (shell3d.py)"][0], res["dx=1/2 dt=0.0125 (512x512x16)"][0]
        print(f"patch size, fine / coarse grid: " + " ".join(f"{r:.2f}" for r in fine / base))
        assert np.all(np.abs(fine / base - 1) < 0.1), "patch sizes not converged in grid spacing"
        assert abs(res["dx=1 dt=0.05 (shell3d.py)"][1] - res["dx=1/2 dt=0.0125 (512x512x16)"][1]) < 0.05
    else:  # majority takeover, physical box 128 x 128 x 8 so the fine grid (256x256x16) stays affordable
        for name, (dx, dt) in grids.items():
            label = name.replace("512x512x16", "256x256x16")    # this mode's box is half as wide
            for f in (0.5, 0.55):
                ends = [run(f, [1, 200], dx, dt, seed=s, lx=128) for s in range(2)]
                print(f"{label:32} f={f}: " + "  ".join(f"{r[0][2]:.3f}->{r[1][2]:.3f}" for r in ends), flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "growth")
