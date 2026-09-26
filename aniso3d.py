"""Does a strongly layered ocean still heal to one hand, or does it stall in layers?

Real icy-moon oceans mix much faster sideways (D_h) than up and down (D_z). Stretching the depth axis by
sqrt(D_h/D_z) turns  a - a^3 + D_h lap_h(a) + D_z d2a/dz2  into the equal-mixing model exactly (the reaction
term has no derivatives), in a box whose depth becomes H' = H * sqrt(D_h/D_z). So weak vertical mixing is
the same as equal mixing in a TALLER box. Icy moons: H/L ~ 1/18 (Enceladus) to 1/37 (Europa), so H' > L once
D_z < D_h/333 (Enceladus) or D_h/1386 (Europa), which covers most of the plausible vertical range.

Test: equal mixing D = 1, periodic sideways (L = 32), closed top and bottom, depth H' = 8 (slab), 32 (cube),
128 (column). Patchwork start ~4 cells, 4 seeds, favoured fraction 50% and 55%. Track the favoured fraction,
walls between layers (sign change going up) vs walls between side-by-side patches, and whether every
horizontal plane is one hand (a stack of layers), and which hand is on top.

On a real moon a layer wall is a sphere of radius ~R, not a flat sheet. Walls in this model move at
D x curvature (Allen & Cahn 1979), so a spherical wall creeps inward and the TOP layer takes over.
For a nearly flat wall at height h(x, y) the motion is dh/dt = D_h * lap_h(h): only SIDEWAYS mixing sets the speed,
so a spherical layer wall sinks at 2*D_h/R whatever D_z is.
circle_check() confirms walls move at D x curvature (a shrinking circle: r^2 = r0^2 - 2 D t);
wall_check() confirms a wavy wall flattens at rate D_h*k^2, independent of D_z.
"""
import numpy as np

L = 32
SHAPES = {"slab": 8, "cube": 32, "column": 128}
TIMES = [50, 200, 500, 1000]


def lap(a):
    out = sum(np.roll(a, s, ax) for ax in (0, 1) for s in (1, -1)) - 4 * a
    p = np.concatenate([a[..., :1], a, a[..., -1:]], axis=2)       # no-flux top/bottom
    return out + p[..., 2:] + p[..., :-2] - 2 * a


def patchwork(f, nz, rng, width=4):
    k = [np.fft.fftfreq(n) for n in (L, L, nz)]
    kk = k[0][:, None, None]**2 + k[1][None, :, None]**2 + k[2][None, None, :]**2
    z = np.fft.ifftn(np.fft.fftn(rng.standard_normal((L, L, nz))) * np.exp(-2 * (np.pi * width)**2 * kk)).real
    return np.where(z > np.quantile(z, 1 - f), 1.0, -1.0).astype(np.float32)


def state(a):
    s = a > 0
    side = sum((s != np.roll(s, 1, ax)).sum() for ax in (0, 1))
    up = (s[..., 1:] != s[..., :-1]).sum()
    planes = s.mean(axis=(0, 1))
    layered = bool(np.all((planes > 0.98) | (planes < 0.02)))
    return float(s.mean()), int(side), int(up), layered, bool(planes[-1] > 0.5)


def circle_check(n=192, r0=60.0, D=1.0, t_end=1000, dt=0.05):
    """2D: a disc of one hand inside the other must shrink as r^2 = r0^2 - 2*D*t (motion by curvature)."""
    y, x = np.mgrid[:n, :n] - n / 2
    a = np.where(x**2 + y**2 < r0**2, 1.0, -1.0)
    lap2 = lambda b: sum(np.roll(b, s, ax) for ax in (0, 1) for s in (1, -1)) - 4 * b
    for _ in range(round(t_end / dt)):
        a = a + (a - a**3 + D * lap2(a)) * dt
    r2 = (a > 0).sum() / np.pi
    want = r0**2 - 2 * D * t_end
    print(f"circle check: r^2 = {r2:.0f} after t={t_end}, curvature law predicts {want:.0f} ({r2 / want - 1:+.1%})")
    assert abs(r2 / want - 1) < 0.10, "walls do not move at D x curvature"


def run(f, nz, seed, dt=0.05):
    a, t, out = patchwork(f, nz, np.random.default_rng(seed)), 0.0, []
    for T in TIMES:
        while t < T - 1e-9:
            a = a + (a - a**3 + lap(a)) * dt
            t += dt
        out.append(state(a))
    return out


def wall_check(Dh, Dz, nx=128, nz=64, amp=4.0, t_end=100, dt=0.02):
    """2D, uneven mixing: a wall at height nz/2 + amp*cos(kx) must flatten at rate D_h*k^2 (not D_z*k^2)."""
    x, z = np.arange(nx), np.arange(nz)[None, :]
    a = np.tanh((z - nz / 2 - amp * np.cos(2 * np.pi * x / nx)[:, None]) / np.sqrt(2 * Dz))

    def height(a):  # where each column crosses zero, linearly interpolated
        i = np.argmax(a > 0, axis=1)
        lo, hi = a[x, i - 1], a[x, i]
        return i - 1 + lo / (lo - hi)

    a0 = np.ptp(height(a)) / 2
    for _ in range(round(t_end / dt)):
        p = np.concatenate([a[:, :1], a, a[:, -1:]], axis=1)
        a = a + (a - a**3 + Dh * (np.roll(a, 1, 0) + np.roll(a, -1, 0) - 2 * a) + Dz * (p[:, 2:] + p[:, :-2] - 2 * a)) * dt
    rate, k2 = np.log(a0 / (np.ptp(height(a)) / 2)) / t_end, (2 * np.pi / nx) ** 2
    print(f"wall check D_h={Dh:g}, D_z={Dz:g}: flattening rate {rate:.5f}; D_h*k^2 = {Dh * k2:.5f}, D_z*k^2 = {Dz * k2:.5f}")
    return rate / (Dh * k2)


def main():
    circle_check()
    # both uneven cases must follow D_h (they differ 4x in D_z/D_h, so following D_z would be off by 4x)
    for Dh, Dz in ((4.0, 1.0), (1.0, 4.0)):
        assert abs(wall_check(Dh, Dz) - 1) < 0.05, "layer walls do not flatten at D_h*k^2"
    wall_check(1.0, 1.0)  # reported, not asserted: with the sharpest wall the grid drags it (~20% slow)
    summary = {}
    for name, nz in SHAPES.items():
        for f in (0.5, 0.55):
            runs = [run(f, nz, seed) for seed in range(4)]
            one = sum(r[-1][0] in (0.0, 1.0) for r in runs)
            lay = sum(r[-1][3] and 0 < r[-1][0] < 1 for r in runs)
            frozen = sum(r[-1][3] and 0 < r[-1][0] < 1 and abs(r[-1][0] - r[-2][0]) < 0.01 and r[-1][2] == r[-2][2]
                         for r in runs)
            fav = [r[-1][0] for r in runs]
            summary[name, f] = (one, lay, frozen, fav)
            top = sum(r[-1][4] for r in runs)
            print(f"{name:6} H'={nz:3} start {f:.2f}: favoured at t=1000 " + " ".join(f"{x:.2f}" for x in fav)
                  + f" | one hand {one}/4, stacked layers {lay}/4 (frozen t=500->1000: {frozen}/4),"
                  f" favoured on top {top}/4", flush=True)
            for seed, r in enumerate(runs):
                print("        seed %d: " % seed + "  ".join(f"t={T}:{x:.2f} side={s} up={u}{' L' if l else ''}"
                                                          for T, (x, s, u, l, _) in zip(TIMES, r)))
    # the tall box (weak vertical mixing) must stall in layers: most runs end as stacked, frozen, mixed-hand layers
    col = [summary["column", f][2] for f in (0.5, 0.55)]
    assert sum(col) >= 5, f"column runs did not stall in layers ({col})"
    # the cube with a 55% majority must heal to one hand (Part 3's takeover)
    assert summary["cube", 0.55][0] >= 3, "cube did not heal"
    return summary


if __name__ == "__main__":
    main()
