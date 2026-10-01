"""Pictures of coarsening (asked for by A. Brandenburg): what the patchwork looks like as it heals or layers.

Row 1: 3D cube (96^3, equal mixing D = 1), 55% favoured start, horizontal slice at t = 1, 20, 100, 400.
Row 2: tall "column" box (32 x 32 x 128, = weak vertical mixing after stretching, aniso3d.py), vertical
       slice at the same times: the patchwork freezes into stacked layers of opposite hands.
Right: patch size l(t) in the cube vs the sqrt(t) law (fit t = 20-100, before patches span the box).
The column is also run on to t = 5000 to check the layers stay frozen (long-run check).
Same model and scheme as coarsen3d.py / aniso3d.py (explicit Euler, 7-point Laplacian).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from coarsen3d import step, patchwork, coarsen
from domains3d import domain_length
import aniso3d

TIMES = [1, 20, 100, 400]


def cube_frames(n=96, f=0.55, seed=1):
    a = patchwork(f, n, 1, 4, np.random.default_rng(seed))
    t, frames, ls = 0.0, [], []
    for T in TIMES:
        while t < T - 1e-9:
            a = step(a, 1.0, 0.05); t += 0.05
        frames.append(a[0, :, :, 0].copy()); ls.append((T, domain_length(a), float((a > 0).mean())))
    return frames, ls


def column_frames(f=0.55, seed=1, t_long=5000):
    """seed 1 is one of aniso3d.py's frozen-layer runs; also run it on to t_long (12x longer) to check it stays frozen."""
    a = aniso3d.patchwork(f, 128, np.random.default_rng(seed))
    t, frames = 0.0, []
    for T in TIMES + [t_long]:
        while t < T - 1e-9:
            a = a + (a - a**3 + aniso3d.lap(a)) * 0.05; t += 0.05
        frames.append(a[:, 0, :].T.copy())
    return frames, [float((f > 0).mean()) for f in frames]


def main():
    cube, ls = cube_frames()
    col, col_fav = column_frames()
    fig = plt.figure(figsize=(12, 6.2))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 1.4])
    for i, T in enumerate(TIMES):
        ax = fig.add_subplot(gs[0, i]); ax.imshow(cube[i], cmap="coolwarm_r", vmin=-1, vmax=1)
        ax.set_title(f"cube, t = {T}\nfavoured {ls[i][2]:.0%}", fontsize=9); ax.axis("off")
        j = i if i < 3 else 4                           # last panel: the t = 5000 long run
        ax = fig.add_subplot(gs[1, i]); ax.imshow(col[j], cmap="coolwarm_r", vmin=-1, vmax=1, origin="lower", aspect=0.25)
        ax.set_title(f"layered, side view\nt = {(TIMES + [5000])[j]}, favoured {col_fav[j]:.0%}", fontsize=9); ax.axis("off")
    ax = fig.add_subplot(gs[:, 4])
    grow = coarsen(1.0, 0.5, [10, 20, 40, 80], n=128)       # coarsen3d.py's D3 growth run (random start, f = 0.5)
    ts, l = np.array([r[0] for r in grow]), np.array([r[1] for r in grow])
    ax.loglog(ts, l, "o-", label="measured, 128³ cube")
    ax.loglog(ts, 5.5 * np.sqrt(ts), "k--", lw=0.8, label="l = 5.5·√(D t)")
    ax.set_xticks([10, 20, 40, 80], ["10", "20", "40", "80"]); ax.minorticks_off()
    ax.set(xlabel="time (reaction units)", ylabel="patch size (cells)", title="patch growth (D = 1)")
    ax.legend(fontsize=8)
    fig.suptitle("Blue = favoured hand, red = other. Top: equal mixing heals to the majority. "
                 "Bottom: weak vertical mixing freezes into layers (unchanged to t = 5000).", fontsize=10)
    fig.tight_layout()
    fig.savefig("coarsening.png", dpi=150)
    slope = np.polyfit(np.log(ts), np.log(l), 1)[0]
    print("cube: " + "  ".join(f"t={T}: l={x:.1f} favoured {p:.3f}" for T, x, p in ls) + f"  exponent {slope:.2f}")
    print("column (side slice): favoured " + "  ".join(f"t={T}: {p:.3f}" for T, p in zip(TIMES + [5000], col_fav))
          + "; wrote coarsening.png")
    assert abs(col_fav[-1] - col_fav[-2]) < 0.01, "layers did not stay frozen over a 12x longer run"
    assert ls[-1][2] > ls[0][2] + 0.1, "majority did not grow in the cube"
    assert 0 < col_fav[-1] < 1, "column healed instead of layering"
    assert 0.35 < slope < 0.6


if __name__ == "__main__":
    main()
