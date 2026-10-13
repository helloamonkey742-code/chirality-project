"""Figure 6, redrawn (H26): shares2.png. Success share per Table 3 world from multiwall.log Section 4
("P_pved %, no lid / final" table): each bar runs from the cascade x10 column (low end) to the multiwall column
(main), with the lid test (solid) and, where a lid can form, without it (hatched). Replaces the base-model
rnaas_fig.png (make_rnaas_fig.py, ocean_worlds.log), which is left untouched.
Self-check: every number below is re-read from multiwall.log and must match to 0.1 point.
Run: /opt/anaconda3/bin/python make_fig_shares2.py
"""
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (world, status, no lid: x10 low, main; final: x10 low, main), copied from multiwall.log Section 4; Table 3 order
ROWS = [
    ("Ganymede", "confirmed", 33.5, 44.1, 33.5, 44.1),
    ("Europa", "confirmed", 29.0, 39.6, 22.0, 30.1),
    ("Enceladus", "confirmed", 19.9, 28.3, 0.0, 0.0),
    ("Callisto", "likely", 33.0, 43.5, 33.0, 43.5),
    ("Triton", "likely", 26.5, 37.0, 26.5, 37.0),
    ("Mimas", "likely", 11.6, 17.3, 3.7, 5.6),
    ("Titan", "disputed", 25.5, 35.8, 25.5, 35.8),
    ("Pluto", "disputed", 30.6, 41.1, 30.6, 41.1),
    ("Umbriel", "possible", 43.6, 53.7, 43.6, 53.7),
    ("Ariel", "possible", 39.7, 49.6, 39.7, 49.6),
    ("Oberon", "possible", 40.8, 51.4, 40.8, 51.4),
    ("Titania", "possible", 40.9, 51.3, 40.9, 51.3),
    ("Dione", "possible", 32.2, 42.0, 20.1, 26.2),
    ("Miranda", "past only", 21.3, 27.1, 7.0, 9.0),
]
COLOR = {"confirmed": "#1f4e79", "likely": "#4a86b8", "disputed": "#9aa5b1",
         "possible": "#c9d3dd", "past only": "#e3c9a8"}


def check():
    """Parse 'world  sweep3 / sweep3 | bound | multiwall | cascade x10 | sum rule' lines of Section 4."""
    sec = open("multiwall.log").read().split("=== Section 4")[1]
    block = sec.split("P_pved %, no lid / final")[1].split("\n\n")[0]
    num = r"\s*([0-9.]+)\s*/\s*([0-9.]+)\s*"
    pat = re.compile(rf"^(\S[\w ]*?)\s{{2,}}{num}\|{num}\|{num}\|{num}\|{num}$", re.M)
    got = {m.group(1): [float(x) for x in m.groups()[1:]] for m in pat.finditer(block)}
    assert len(got) == 15, sorted(got)
    for w, _, nl_lo, nl_hi, f_lo, f_hi in ROWS:
        v = got[w]          # [s3 nl, s3 f, bound nl, bound f, mw nl, mw f, x10 nl, x10 f, sum nl, sum f]
        assert (nl_lo, nl_hi, f_lo, f_hi) == (v[6], v[4], v[7], v[5]), (w, v)
    return got


def r(x):
    return int(x + 0.5)      # round half up, as Table 3 (26.5 -> 27), not Python's half-even


check()
fig, ax = plt.subplots(figsize=(6.5, 5.4))
y = list(range(len(ROWS)))[::-1]
for yi, (w, st, nl_lo, nl_hi, f_lo, f_hi) in zip(y, ROWS):
    lid = (nl_lo, nl_hi) != (f_lo, f_hi)
    if lid:     # no-lid range, hatched, just above the solid final range
        ax.barh(yi + 0.2, nl_hi - nl_lo, left=nl_lo, height=0.32, color="none", edgecolor="#777", hatch="///", lw=0.6)
    h, off = (0.32, -0.2) if lid else (0.6, 0.0)
    ax.barh(yi + off, max(f_hi - f_lo, 0.4), left=f_lo, height=h, color=COLOR[st], edgecolor="#333", lw=0.6)
    lab = f"{r(f_lo)}–{r(f_hi)}%" if f_hi > 0 else "0%"
    ax.text(max(f_hi, nl_hi if lid else 0) + 0.8, yi, lab + (f" ({r(nl_lo)}–{r(nl_hi)}%)" if lid else ""),
            va="center", fontsize=7)
ax.set_yticks(y, [r[0] for r in ROWS])
ax.set_xlabel("Draws in which the PVED-favoured hand is selected (%)")
ax.set_xlim(0, 70)
handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="#333") for c in COLOR.values()]
handles.append(plt.Rectangle((0, 0), 1, 1, fc="none", ec="#777", hatch="///"))
ax.legend(handles, [f"ocean {k}" for k in COLOR] + ["no stratified lid (in brackets)"],
          fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
ax.text(0.99, 1.01, "bar: wall times ×10 (left end) to as computed (right end)", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=7, color="0.3")
fig.tight_layout()
fig.savefig("shares2.png", dpi=300)
print(f"self-check passed: {len(ROWS)} worlds match multiwall.log Section 4; wrote shares2.png (300 dpi)")
