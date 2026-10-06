"""Simulation check for earth.py: does the bigger hemisphere really win, or does the field pin the wall?

Earth's field gives the two hemispheres opposite local biases +h / -h. Each side picks its own hand,
leaving one wall at the bias boundary. Curvature pushes the wall into the smaller side (majority wins,
Part 3), but the bias on the smaller side pushes back. Axisymmetric reaction-diffusion on a sphere
(time in 1/k, length in R, eps = sqrt(D/k)/R = wall width / radius):

    da/dt = eps^2 (1/sin th) d/dth (sin th da/dth) + a - a^3 + h(th),   h = +b (th < th_b), -b (th > th_b)

The + side covers a share (1+A)/2 of the sphere (cos th_b = -A). Analytic guess: a bias b moves a flat
wall at (3/sqrt2)*eps*b, curvature at eps^2*cot(th_b) = eps^2*A/sqrt(1-A^2), so the wall is pinned when
    b > b_crit = (sqrt2/3) * eps * A / sqrt(1-A^2)
Measured: final + share for b = (0.25..4) x b_crit.
"""
import numpy as np


def b_crit(eps, A):
    return np.sqrt(2) / 3 * eps * A / np.sqrt(1 - A**2)


def run(eps, A, bs, t_max):
    n = int(np.ceil(np.pi / (eps / 6)))
    dth = np.pi / n
    th = (np.arange(n) + 0.5) * dth
    face = np.sin(np.arange(n + 1) * dth)            # sin at faces; 0 at both poles -> no flux
    coef = eps**2 / (np.sin(th) * dth**2)
    th_b = np.arccos(-A)
    h = np.where(th < th_b, 1.0, -1.0)[None, :] * np.asarray(bs)[:, None]
    a = np.sign(h)                                   # each side has picked its own hand
    area = np.sin(th) * dth / 2                      # share of the sphere per cell (sums to 1)
    dt = min(0.05, 0.2 * dth**2 / eps**2)
    for _ in range(int(t_max / dt)):
        flux = np.zeros((len(bs), n + 1))
        flux[:, 1:-1] = face[1:-1] * np.diff(a, axis=1)
        a = a + dt * (coef * np.diff(flux, axis=1) + a - a**3 + h)
    return ((a > 0) * area).sum(axis=1)               # final + share


def main():
    mult = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 8.0, 16.0]
    print(f"{'eps':>5} {'A':>4} {'start':>6}  final + share at b/b_crit = " + " ".join(f"{m:>5}" for m in mult))
    fits = []
    for eps, A in [(0.2, 0.3), (0.1, 0.3), (0.2, 0.1), (0.1, 0.1), (0.05, 0.3), (0.05, 0.1)]:
        t_max = 20 / (eps**2 * A)
        share = run(eps, A, [m * b_crit(eps, A) for m in mult], t_max)
        start = (1 + A) / 2
        print(f"{eps:5.2f} {A:4.1f} {start:6.2f}  " + " " * 30 + " ".join(f"{s:5.2f}" for s in share), flush=True)
        won = [m for m, s in zip(mult, share) if s > 0.99]
        pinned = [m for m, s in zip(mult, share) if s < 0.99]
        fits.append((max(won) if won else None, min(pinned) if pinned else None))
    print("threshold bracket (b/b_crit): bigger side wins up to / pinned from:", fits)


if __name__ == "__main__":
    main()
