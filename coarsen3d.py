"""3D patches keep growing after the choice (curvature-driven coarsening).

After the transition (lam = 1, no bias) run the 3D field and check:
  D3. patch size grows as  l = A * sqrt(D * t)
  E3. the majority hand takes over: start with fraction f of the favoured hand in
      random small patches; the minority disappears (non-conserved dynamics).
If both hold, a patchwork ocean heals within t ~ (L/A)**2 / D, and the winner is the
ocean-wide majority -- whose excess is decided by ALL the ocean's molecules together.
"""
import numpy as np
from domains3d import domain_length


def step(a, D, dt):
    lap = sum(np.roll(a, s, ax) for ax in (1, 2, 3) for s in (1, -1)) - 6 * a
    return a + (a - a**3 + np.float32(D) * lap) * np.float32(dt)


def patchwork(f, n, runs, width, rng):
    """Saturated +-1 patches of size ~width; the favoured hand fills fraction f."""
    k = np.fft.fftfreq(n)
    kk = k[:, None, None]**2 + k[None, :, None]**2 + k[None, None, :]**2
    z = np.fft.ifftn(np.fft.fftn(rng.standard_normal((runs, n, n, n)), axes=(1, 2, 3))
                     * np.exp(-2 * (np.pi * width)**2 * kk), axes=(1, 2, 3)).real
    cut = np.quantile(z, 1 - f, axis=(1, 2, 3), keepdims=True)
    return np.where(z > cut, 1.0, -1.0).astype(np.float32)


def coarsen(D, f, times, n=128, runs=1, seed=0, width=0):
    rng = np.random.default_rng(seed)
    if width:
        a = patchwork(f, n, runs, width, rng)
    else:
        a = np.where(rng.random((runs, n, n, n)) < f, 1.0, -1.0).astype(np.float32)
    dt, t, out = min(0.05, 0.1 / D), 0.0, []
    for T in times:
        while t < T:
            a = step(a, D, dt); t += dt
        out.append((T, domain_length(a), float((a > 0).mean())))
    return out


def main():
    print("D3. patch growth, f = 0.5", flush=True)
    for D in (1.0, 2.25):
        rows = coarsen(D, 0.5, [10, 20, 40, 80])
        ts, ls = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
        slope = np.polyfit(np.log(ts), np.log(ls), 1)[0]
        A = np.mean(ls / np.sqrt(D * ts))
        print(f"  D={D}: " + "  ".join(f"t={T:g}:l={l:.1f}" for T, l, _ in rows) + f"  exponent {slope:.2f}  A={A:.2f}", flush=True)
        assert abs(slope - 0.5) < 0.1, "not sqrt(t) growth"

    print("\nE3. pre-formed patches (~4 cells): does the majority take over? (n=128)", flush=True)
    for f in (0.5, 0.52, 0.55):
        rows = coarsen(1.0, f, [1, 50, 200, 400], n=128, runs=1, seed=1, width=4)
        print("  f={}: ".format(f) + "  ".join(f"t={T}: favoured {p:.3f}" for T, _, p in rows), flush=True)
        if f > 0.5:
            assert rows[-1][2] > rows[0][2] + 0.02, "majority did not grow"


if __name__ == "__main__":
    main()
