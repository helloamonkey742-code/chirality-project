"""Figure 5: tornado plot of each input's effect on the success share.

Built only from ../uncertainty.log (no reruns). For each world, each input's
range is split in half; the bar runs from the success share in the low half
to the share in the high half, and the vertical line is the full-range share.
Self-check: every parsed line is re-formatted and must equal the log line.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "..", "uncertainty.log")
OUT = os.path.join(HERE, "fig_sensitivity.png")

HEAD = re.compile(r"^(.+): physics wins in (\d+)% of plausible input space$")
ROW = re.compile(r"^  (.{15}) low half +(\d+)%  high half +(\d+)%   \(effect ([+-]\d+)%\)$")


def parse(path):
    worlds, cur = {}, None
    for raw in open(path).read().splitlines():
        m = HEAD.match(raw)
        if m:
            cur = m.group(1)
            worlds[cur] = {"base": int(m.group(2)), "rows": []}
            continue
        m = ROW.match(raw)
        if m and cur:
            name, lo, hi, eff = m.group(1).strip(), int(m.group(2)), int(m.group(3)), m.group(4)
            # self-check 1: re-formatted line equals the log line exactly
            assert raw == "  %-15s low half %4d%%  high half %4d%%   (effect %s%%)" % (name, lo, hi, eff), raw
            # self-check 2: the stated effect is high minus low (±1 for rounding, "-0" allowed)
            assert abs((hi - lo) - int(eff)) <= 1, raw
            worlds[cur]["rows"].append((name, lo, hi))
    return worlds


def main():
    w = parse(LOG)
    names = ["Enceladus", "Europa", "Early Earth ocean"]
    assert list(w) == names, list(w)
    # self-check 3: values quoted in the paper (Sections 3.3, 2.6)
    assert [w[n]["base"] for n in names] == [22, 29, 75]
    dz = {n: dict((r[0], r[1:]) for r in w[n]["rows"])["vertical D_z"] for n in names}
    assert dz == {"Enceladus": (0, 43), "Europa": (0, 57), "Early Earth ocean": (52, 97)}, dz
    assert all(len(w[n]["rows"]) == 8 for n in names)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=False)
    for ax, n in zip(axes, names):
        rows = sorted(w[n]["rows"], key=lambda r: abs(r[2] - r[1]))  # largest effect on top
        for i, (name, lo, hi) in enumerate(rows):
            left, right = min(lo, hi), max(lo, hi)
            ax.barh(i, right - left, left=left, color="#3b6ea5" if hi >= lo else "#b5533c", height=0.6)
            ax.text(right + 1.5, i, "%d→%d" % (lo, hi), va="center", fontsize=7)
        ax.axvline(w[n]["base"], color="k", lw=0.8, ls="--")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows], fontsize=8)
        ax.set_xlim(0, 115)
        ax.set_xlabel("success share (%)", fontsize=8)
        ax.set_title("%s (all inputs: %d%%)" % (n.replace(" ocean", ""), w[n]["base"]), fontsize=9)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(OUT, dpi=300)
    print("wrote", OUT, "; self-checks passed for", sum(len(w[n]["rows"]) for n in names), "log lines")


if __name__ == "__main__":
    main()
