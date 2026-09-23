"""Where was the hand decided? Open ocean vs lake vs pond vs vent pore.

Same selection formula as ocean.py, applied to each origin-of-life setting (small
settings are treated as well mixed). If life began in one isolated setting, that
setting's P(favoured) is the chance the weak force, not luck, picked the hand.

Concentrations of the reacting molecules are rough ranges (see README); only the
vent-pore concentration factor (Baaske et al. 2007) and pond size (Pearce et al. 2017)
come directly from papers. Everything else is a labelled assumption.
"""
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ocean import NA, YR, delta_at

# name: (volume m^3, (c_low, c_high) M, time available s)
SETTINGS = {
    "Open ocean": (1.332e18, (1e-7, 1e-4), 1e8 * YR),
    "Lake / lagoon (1 km² × 10 m)": (1e7, (1e-5, 1e-2), 1e4 * YR),
    "Warm little pond (1 m × 1 m)": (np.pi, (1e-6, 1e-2), 1e2 * YR),
    "Vent pore (1 mm³)": (1e-9, (1e-4, 1.0), 1e5 * YR),   # ocean c x 1e3..1e8, capped at 1 M
}


def p_favoured(c, g, k2, V, tau):
    if k2 * c * tau < 10:
        return np.nan                                    # chemistry never finishes
    return norm.cdf(delta_at(c, g, k2, V, tau))


def main():
    assert p_favoured(1e-6, 1e-17, 1e-3, *SETTINGS["Open ocean"][::2]) > 0.999
    assert p_favoured(1e-2, 1e-17, 1e-3, *SETTINGS["Warm little pond (1 m × 1 m)"][::2]) < 0.51

    print("P(weak force picks the hand), g = 1e-17 (1e-16 in brackets); '-' = never finishes\n")
    print(f"{'setting':30} {'conc.':>7}  {'k2=1e-6':>15} {'k2=1e-3':>15} {'k2=1':>15}")
    best = {}
    for name, (V, (lo, hi), tau) in SETTINGS.items():
        for c in (lo, hi):
            cells = []
            for k2 in (1e-6, 1e-3, 1.0):
                p, p16 = p_favoured(c, 1e-17, k2, V, tau), p_favoured(c, 1e-16, k2, V, tau)
                cells.append("-" if np.isnan(p) else f"{p:.3f} ({p16:.3f})")
                if not np.isnan(p16):
                    best[name] = max(best.get(name, 0.5), p16)
            print(f"{name:30} {c:7.0e}  " + " ".join(f"{x:>15}" for x in cells))

    fig, ax = plt.subplots(figsize=(7, 3.2))
    names = list(best)
    ax.barh(names, [best[n] - 0.5 for n in names], left=0.5, color=["#2a9d8f" if best[n] > 0.977 else "#aaa" for n in names])
    for i, n in enumerate(names):
        ax.text(max(best[n], 0.5) + 0.005 if best[n] < 0.9 else 0.99, i, f"{best[n]:.3f}" + (" (coin flip)" if best[n] < 0.55 else ""),
                va="center", ha="left" if best[n] < 0.9 else "right", fontsize=8, color="0.2" if best[n] < 0.9 else "w")
    ax.axvline(0.977, color="0.3", ls="--", lw=0.8)
    ax.text(0.975, -0.6, "physics wins (97.7%)", ha="right", fontsize=8, color="0.3")
    ax.set_xlim(0.5, 1.0)
    ax.set_xlabel("best-case chance the weak force picks the hand (g = 1e-16)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig("scenarios.png", dpi=150)
    print("\nwrote scenarios.png")


if __name__ == "__main__":
    main()
