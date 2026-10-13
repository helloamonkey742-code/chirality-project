"""Figure 3, redrawn (H26): threshold2.png. Same curves as ocean.py's threshold.png (ocean.c_size / c_time,
k2 = 1e-3 /M/s), with the grey band moved to the sourced amino-acid g range 4.4e-18 .. 2.64e-16 (0.10-6.2 aeV;
sweep2.py, sweep2.log g_min / g_max lines, journals/new_analyses_sweep2.md) and no anthropomorphic title.
Self-check: the band edges match sweep2.log, and each plotted curve matches ocean.py at the band edges.
Run: /opt/anaconda3/bin/python make_fig_threshold2.py   (does not touch threshold.png)
"""
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ocean import BODIES, c_size, c_time, delta_at, DELTA_TARGET

G_BAND = (4.4e-18, 2.64e-16)        # copied from sweep2.log; asserted below
K2 = 1e-3

log = open("sweep2.log").read()
g_log = tuple(float(re.search(rf"^{k}\s+([0-9.e+-]+)", log, re.M).group(1)) for k in ("g_min", "g_max"))
assert g_log == G_BAND, (g_log, G_BAND)

gs = np.logspace(-19, -15, 400)
curves = {name: np.maximum(c_size(gs, K2, V, tau), c_time(K2, tau)) for name, (V, tau) in BODIES.items()}
for name, (V, tau) in BODIES.items():   # plotted line = max(size, time) limit; size limit gives Delta = 2
    for g in G_BAND:
        c = max(c_size(g, K2, V, tau), c_time(K2, tau))
        assert np.isclose(np.interp(np.log10(g), np.log10(gs), np.log10(curves[name])), np.log10(c), atol=1e-3)
        assert delta_at(c, g, K2, V, tau) >= DELTA_TARGET * (1 - 1e-9)

fig, ax = plt.subplots(figsize=(7, 4.5))
for name, y in curves.items():
    ax.loglog(gs, y, label=name)
ax.axvspan(*G_BAND, color="0.9", zorder=0, label="published amino-acid PVED range (0.10–6.2 aeV)")
for c_ref, lab in [(1e-9, "1 nM"), (1e-6, "1 µM")]:
    ax.axhline(c_ref, color="0.5", lw=0.8, ls="--")
    ax.text(1.2e-19, c_ref * 1.4, lab, color="0.4", fontsize=8)
ax.set_xlabel("parity-violating bias  g = ΔE/kT")
ax.set_ylabel("minimum concentration (M)")
ax.set_title("Minimum concentration for PVED selection  (k2 = 1e-3 /M/s)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("threshold2.png", dpi=300)
print("self-check passed: band", G_BAND, "= sweep2.log; wrote threshold2.png (300 dpi)")
