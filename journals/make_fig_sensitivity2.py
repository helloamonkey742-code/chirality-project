"""Figure 5 (rerun): tornado plot of each input's effect on the PVED share, current pipeline.

Built only from ../rerun_base.log section 1 (no reruns), main cascade time. Same style as make_fig_sensitivity.py:
each input's draws are split at the median; the bar runs from the share in the low half to the share in the high
half; the dashed line is the share over all draws. Self-checks: every parsed line re-formats to the log line exactly,
the stated effect is high minus low, and the regression line of section 0 gives the same all-draw share.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "..", "rerun_base.log")
OUT = os.path.join(HERE, "fig_sensitivity2.png")

HEAD = re.compile(r"^(.+) \((main|cascade x10|x10)\): PVED-set hand in (\d+\.\d)% of draws, no lid$")
ROW = re.compile(r"^  (.{15}) low half +(\d+\.\d)%  high half +(\d+\.\d)%   \(effect +([+-]\d+\.\d)\)$")
NAMES = ["Enceladus", "Europa", "Early Earth ocean"]


def parse(path):
    lines = open(path).read().splitlines()
    worlds, cur = {}, None
    for raw in lines:
        m = HEAD.match(raw)
        if m:
            cur = (m.group(1), m.group(2))
            worlds[cur] = {"base": float(m.group(3)), "rows": []}
            continue
        m = ROW.match(raw)
        if m and cur:
            name, lo, hi, eff = m.group(1).strip(), float(m.group(2)), float(m.group(3)), float(m.group(4))
            assert raw == f"  {name:15} low half {lo:5.1f}%  high half {hi:5.1f}%   (effect {eff:+5.1f})", raw
            assert abs((hi - lo) - eff) <= 0.11, raw
            worlds[cur]["rows"].append((name, lo, hi))
        elif cur and not raw.startswith("  "):
            cur = None
    # regression lines of section 0: "<world> <main no lid> / <final> (...multiwall.log...)"
    for n in NAMES:
        reg = [l for l in lines if l.startswith(f"{n:18} ") and "multiwall.log:" in l]
        assert len(reg) == 1, n
        assert float(reg[0][19:].split("/")[0]) == worlds[n, "main"]["base"], (n, reg[0])
    return worlds


def main():
    w = parse(LOG)
    assert all((n, "main") in w and len(w[n, "main"]["rows"]) == 11 for n in NAMES), list(w)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=False)
    for ax, n in zip(axes, NAMES):
        d = w[n, "main"]
        rows = sorted(d["rows"], key=lambda r: abs(r[2] - r[1]))  # largest effect on top
        for i, (name, lo, hi) in enumerate(rows):
            left, right = min(lo, hi), max(lo, hi)
            ax.barh(i, right - left, left=left, color="#3b6ea5" if hi >= lo else "#b5533c", height=0.6)
            ax.text(right + 1.5, i, "%.1f→%.1f" % (lo, hi), va="center", fontsize=7)
        ax.axvline(d["base"], color="k", lw=0.8, ls="--")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows], fontsize=8)
        ax.set_xlim(0, 115)
        ax.set_xlabel("PVED share (%)", fontsize=8)
        ax.set_title("%s (all inputs: %.1f%%)" % (n.replace(" ocean", ""), d["base"]), fontsize=9)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(OUT, dpi=300)
    print("wrote", OUT, "; self-checks passed for", sum(len(v["rows"]) for v in w.values()), "log lines")


if __name__ == "__main__":
    main()
