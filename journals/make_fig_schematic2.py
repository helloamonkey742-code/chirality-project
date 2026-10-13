"""Figure 1, redrawn (H26): fig_schematic2.png. Copy of make_fig_schematic.py (fig_schematic.png is untouched) with
panel (c) over the sourced amino-acid g range 4.4e-18 .. 2.64e-16 (sweep2.log g_min / g_max; Table 1) instead of
4e-19 .. 4e-16, asserted against sweep2.log together with the N needed at both ends (k tau = 10), and 300 dpi.

(a) Tilted pitchfork: fixed points of lam*a - a^3 + g = 0 (Eq. 1 without noise), g exaggerated so it is visible.
(b) Ocean column cartoon: ice shell, stratified lid, convecting interior, D_h / D_z, a wall between hands. Not to scale.
(c) Molecules needed, N = (2 / (sqrt2 pi^1/4 g (k tau)^1/4))^2 (review.py: n_needed), over the sourced range of g,
    against N = c V 1000 N_A for a 1 m^3 pond at 1 M and the two moon oceans at 1 uM (volumes from Table 1).
No simulated or measured data are plotted. Run: /opt/anaconda3/bin/python journals/make_fig_schematic2.py
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from review import n_needed          # noqa: E402  (the scale rule used in the paper)
from ocean import NA                 # noqa: E402

# self-checks against the numbers quoted in ija_draft.md Sec. 3.1
assert f"{n_needed(1e-17, 10):.1e}" == "3.6e+33"
assert f"{1.0 * 1.0 * 1000 * NA:.0e}" == "6e+26"          # 1 m^3 pond at 1 M
G_RANGE = (4.4e-18, 2.64e-16)                               # sweep2.log g_min, g_max
with open(os.path.join(os.path.dirname(HERE), "sweep2.log")) as fh:
    LOG = fh.read()
for key, g0 in zip(("g_min", "g_max"), G_RANGE):
    m = re.search(rf"^{key}\s+(\S+)\s+(\S+)", LOG, re.M)
    assert float(m.group(1)) == g0 and f"{n_needed(g0, 10):.2e}" == m.group(2), (key, m.groups(), n_needed(g0, 10))

plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "font.family": "DejaVu Sans"})
fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1, 1.1, 1]})

# ---------- (a) tilted pitchfork ----------
G = 0.05                                     # exaggerated bias, for visibility only (real g ~ 1e-17)
a = np.linspace(-1.35, 1.35, 4000)   # even count: skips a = 0
lam = a**2 - G / a                           # every fixed point of lam*a - a^3 + g = 0, solved for lam
stable = lam - 3 * a**2 < 0
for mask, ls in ((stable, "-"), (~stable, "--")):
    y = np.where(mask & (lam > -1) & (lam < 1.6), a, np.nan)
    axa.plot(np.where(np.isnan(y), np.nan, lam), y, ls, color="k", lw=1.1)
lam0 = np.linspace(0, 1.6, 200)              # untilted pitchfork (g = 0), for comparison
for s in (1, -1):
    axa.plot(lam0, s * np.sqrt(lam0), color="0.6", lw=0.7, ls=":")
axa.plot([-1, 0], [0, 0], color="0.6", lw=0.7, ls=":")
axa.axvspan(-0.25, 0.25, color="0.88", lw=0)
axa.text(0, -1.27, "outcome set\nhere", ha="center", va="bottom", fontsize=6.5)
axa.text(-0.95, 0.85, "favoured hand\n(tilted by g)", ha="left", fontsize=6.5)
axa.text(0.95, -0.6, "other hand", ha="left", fontsize=6.5)
axa.annotate("", xy=(-0.8, G / 0.8 + 0.02), xytext=(-0.8, 0.0),
             arrowprops=dict(arrowstyle="->", lw=0.7))
axa.text(-0.95, -0.28, r"offset $\approx g/|\lambda|$", fontsize=6.5)
axa.set_xlim(-1, 1.6)
axa.set_ylim(-1.3, 1.3)
axa.set_xlabel(r"control parameter $\lambda = \gamma t$")
axa.set_ylabel(r"excess $\alpha = (L-D)/(L+D)$")
axa.set_title(r"(a) fixed points, $g$ exaggerated", fontsize=8, loc="left")

# ---------- (b) ocean column (cartoon, not to scale) ----------
ax = axb
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")
ax.add_patch(Rectangle((0, 8.4), 10, 1.4, facecolor="0.92", edgecolor="k", hatch="///", lw=0.6))
ax.text(5, 9.1, "ice shell", ha="center", va="center", fontsize=7, backgroundcolor="0.92")
ax.add_patch(Rectangle((0, 7.0), 10, 1.4, facecolor="#cfe0ee", edgecolor="none"))
ax.text(9.8, 7.7, r"stratified lid ($S<S^*$): low $D_z$", ha="right", va="center", fontsize=6.3)
ax.add_patch(Rectangle((0, 0.6), 10, 6.4, facecolor="#eef4f9", edgecolor="none"))
ax.add_patch(Rectangle((0, 0), 10, 0.6, facecolor="0.55", edgecolor="none"))
ax.text(5, 0.3, "sea floor", ha="center", va="center", fontsize=6.3, color="w")
ax.add_patch(Rectangle((0, 0.6), 4.6, 6.4, facecolor="#f4c7a8", alpha=0.55, edgecolor="none"))
ax.add_patch(Rectangle((4.6, 0.6), 5.4, 6.4, facecolor="#a8c8f4", alpha=0.55, edgecolor="none"))
ax.plot([4.6, 4.6], [0.6, 7.0], color="k", lw=1.4)
ax.text(2.3, 5.9, r"$\alpha \approx +1$" + "\n(L-dominated)", ha="center", fontsize=6.5)
ax.text(7.3, 5.9, r"$\alpha \approx -1$" + "\n(D-dominated)", ha="center", fontsize=6.5)
ax.annotate("wall, width" + "\n" + r"$\sim\sqrt{D/k}$", xy=(4.6, 4.6), xytext=(5.0, 4.2), fontsize=6.3,
            arrowprops=dict(arrowstyle="->", lw=0.6))
for cx in (1.6, 7.6):                        # convection cells
    ax.add_patch(FancyArrowPatch((cx - 0.8, 3.0), (cx + 0.8, 3.0), connectionstyle="arc3,rad=-1.0",
                                 arrowstyle="->", mutation_scale=7, lw=0.7, color="0.3"))
    ax.add_patch(FancyArrowPatch((cx + 0.8, 3.0), (cx - 0.8, 3.0), connectionstyle="arc3,rad=-1.0",
                                 arrowstyle="->", mutation_scale=7, lw=0.7, color="0.3"))
ax.text(2.3, 1.0, "convecting\ninterior", ha="center", fontsize=6.3)
ax.annotate("", xy=(3.8, 4.6), xytext=(0.8, 4.6), arrowprops=dict(arrowstyle="<->", lw=0.8))
ax.text(2.3, 4.75, r"$D_h$", ha="center", fontsize=7.5)
ax.annotate("", xy=(9.4, 5.4), xytext=(9.4, 3.9), arrowprops=dict(arrowstyle="<->", lw=0.8))
ax.text(9.25, 4.65, r"$D_z$", ha="right", va="center", fontsize=7.5)
ax.annotate("", xy=(-0.25, 7.0), xytext=(-0.25, 0.6), arrowprops=dict(arrowstyle="<->", lw=0.6),
            annotation_clip=False)
ax.text(-0.45, 3.8, r"$H$", ha="right", va="center", fontsize=7.5)
ax.set_title("(b) ocean column (not to scale)", fontsize=8, loc="left")

# ---------- (c) molecules needed vs molecules available ----------
g = np.logspace(np.log10(G_RANGE[0]), np.log10(G_RANGE[1]), 200)   # sourced range of g (sweep2.log)
axc.loglog(g, n_needed(g, 10), color="k", lw=1.2, label=r"needed, $k\tau=10$")
axc.loglog(g, n_needed(g, 1e5), color="k", lw=0.8, ls="--", label=r"needed, $k\tau=10^5$")
axc.plot([1e-17], [n_needed(1e-17, 10)], "o", color="k", ms=3)
axc.annotate(r"$3.6\times10^{33}$", xy=(1e-17, n_needed(1e-17, 10)), xytext=(1.5e-17, 2e35), fontsize=6.5,
             arrowprops=dict(arrowstyle="->", lw=0.5))
avail = [("Europa ocean, 1 µM", 3.3e18 * 1e-6, "#3a6ea5"),
         ("Enceladus ocean, 1 µM", 2.7e16 * 1e-6, "#6a9bd1"),
         ("1 m³ pond, 1 M", 1.0 * 1.0, "#b5651d")]
for name, cv, col in avail:                  # N = c V 1000 N_A, fully mixed (f = 1)
    n = cv * 1000 * NA
    axc.axhline(n, color=col, lw=0.9)
    axc.text(G_RANGE[0] * 1.15, n * 1.8, name, fontsize=6.3, color=col)
axc.set_xlim(*G_RANGE)
axc.set_ylim(1e26, 1e41)
axc.set_xlabel(r"bias $g = \Delta E_{pv}/kT$")
axc.set_ylabel("molecules N")
axc.legend(loc="lower left", bbox_to_anchor=(0, 0.1), fontsize=6.3, frameon=False)
axc.set_title(r"(c) scale rule, $\Delta = 2$", fontsize=8, loc="left")

fig.tight_layout(w_pad=1.0)
out = os.path.join(HERE, "fig_schematic2.png")
fig.savefig(out, dpi=300)
print("wrote", out)
