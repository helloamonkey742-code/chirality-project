"""3D version of the transition sweep in domains.py (periodic 3D grid).

Result: the 1D patch-size law fails in 3D because curved patch walls keep moving after
the choice (curvature-driven coarsening) -- tested directly in coarsen3d.py.
"""
import numpy as np


def sweep3d(D, gamma, g=0.0, eps=1e-4, n=128, runs=2, seed=0, lam_end=1.0):
    rng = np.random.default_rng(seed)
    dt = np.float32(min(0.05, 0.1 / D))          # 3D stability: D*dt < 1/6
    a = np.zeros((runs, n, n, n), np.float32)
    noise = np.float32(np.sqrt(eps * dt))
    for t in np.arange(-1 / gamma, lam_end / gamma, dt):
        lap = sum(np.roll(a, s, ax) for ax in (1, 2, 3) for s in (1, -1)) - 6 * a
        a += (np.float32(gamma * t) * a - a**3 + np.float32(g) + np.float32(D) * lap) * dt
        a += noise * rng.standard_normal(a.shape, dtype=np.float32)
    return a


def domain_length(a):
    s = np.sign(a)
    walls = sum((s != np.roll(s, 1, ax)).sum() for ax in (1, 2, 3))
    return 3 * a.size / max(walls, 1)


def law_spread(ratios):
    """Coefficient of variation of the l/(sqrt(D)*gamma^-1/4) ratios (domains.py's C)."""
    ratios = np.asarray(ratios, dtype=float)
    return np.std(ratios) / np.mean(ratios)


def main():
    """Patch size at the end of the sweep. In 3D this does NOT follow the 1D law
    (spread ~31%): patches keep growing after the choice. See coarsen3d.py."""
    ratios = []
    for D, gamma in [(1.0, 0.1), (1.0, 0.025), (2.25, 0.1), (2.25, 0.025)]:
        l = domain_length(sweep3d(D, gamma))
        r = l / (np.sqrt(D) * gamma**-0.25)
        ratios.append(r)
        print(f"D={D:<5} gamma={gamma:<6} domain={l:6.1f} cells  l/(sqrt(D)*gamma^-1/4)={r:.2f}", flush=True)
    spread = law_spread(ratios)
    print(f"spread {spread:.0%}")
    # domains.py asserts spread < 0.25 for the law to "hold" in 1D; here we assert
    # the opposite -- if this ever holds in 3D, the "fails in 3D" claim is wrong.
    assert spread > 0.25, "1D patch law unexpectedly held in 3D (spread <= 25%)"


if __name__ == "__main__":
    main()
