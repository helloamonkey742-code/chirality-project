"""RNAAS figure: weak-force win rate per ocean world, values copied from ocean_worlds.log."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (world, no-lid %, final %, ocean status) from ocean_worlds.log
ROWS = [
    ("Umbriel", 38.7, 38.7, "possible"), ("Ariel", 37.0, 37.0, "possible"),
    ("Oberon", 36.5, 36.5, "possible"), ("Titania", 36.4, 36.4, "possible"),
    ("Pluto", 30.2, 30.2, "disputed"), ("Callisto", 28.5, 28.5, "likely"),
    ("Ganymede", 28.3, 28.3, "confirmed"), ("Triton", 27.7, 27.7, "likely"),
    ("Titan", 24.3, 24.3, "disputed"), ("Europa", 28.6, 21.7, "confirmed"),
    ("Dione", 33.6, 21.0, "possible"), ("Miranda", 26.9, 8.9, "past only"),
    ("Mimas", 14.7, 4.7, "likely"), ("Enceladus", 21.8, 0.0, "confirmed"),
]
COLOR = {"confirmed": "#1f4e79", "likely": "#4a86b8", "disputed": "#9aa5b1",
         "possible": "#c9d3dd", "past only": "#e3c9a8"}

fig, ax = plt.subplots(figsize=(6.5, 4.6))
y = range(len(ROWS))[::-1]
for yi, (w, nolid, final, st) in zip(y, ROWS):
    ax.barh(yi, nolid, color="none", edgecolor="#888", hatch="///", lw=0.6)
    ax.barh(yi, final, color=COLOR[st], edgecolor="#333", lw=0.6)
    ax.text(final + 0.6, yi, f"{final:.0f}%", va="center", fontsize=8,
            bbox=dict(fc="white", ec="none", pad=0.5))
ax.set_yticks(list(y), [r[0] for r in ROWS])
ax.set_xlabel("Draws in which the weak-force-favoured hand wins (%)")
ax.set_xlim(0, 45)
handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="#333") for c in COLOR.values()]
handles.append(plt.Rectangle((0, 0), 1, 1, fc="none", ec="#888", hatch="///"))
ax.legend(handles, [f"ocean {k}" for k in COLOR] + ["if no stratified lid"],
          fontsize=7, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig("rnaas_fig.png", dpi=300)
