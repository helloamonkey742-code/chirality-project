"""Can a lab see a reaction as slow as nature allows? Sizing the slow experiment.

Seeded assay: start with enantiomeric excess e0; a real amplifier (self-copying AND
mutual suppression) makes |ee| grow as e0*exp(k*t) with k = k2*c, sign following the seed.
Plain self-copying without suppression keeps ee constant, so growth tests R1+R2 together.

Detection: ee measured with noise sigma per sample; call it detected when the rise
exceeds 3*sigma*sqrt(2) (difference of two noisy points). Then
    k2_min = ln(1 + 3*sqrt2*sigma/e0) / (c * T)
Compare with nature's floor from scorecard.py (3e-7 /M/s pond, 3e-9 ocean).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

YR = 3.156e7
SIGMA = 0.002       # ee measurement noise (0.2% absolute) -- assumption, check with the lab's method
E0 = 0.05           # seed excess, 5%
FLOOR = {"pond floor": 3e-7, "ocean floor": 3e-9}


def k2_min(c, T, e0=E0, sigma=SIGMA):
    return np.log(1 + 3 * np.sqrt(2) * sigma / e0) / (c * T)


def main():
    # self-check: at k2 = k2_min the seeded ee rises by exactly the detection threshold
    c, T = 0.1, YR
    rise = E0 * (np.exp(k2_min(c, T) * c * T) - 1)
    assert abs(rise - 3 * np.sqrt(2) * SIGMA) < 1e-12

    durations = {"1 month": YR / 12, "6 months": YR / 2, "1 year": YR, "3 years": 3 * YR}
    print(f"smallest detectable k2 (/M/s), seed ee {E0:.0%}, measurement noise {SIGMA:.1%}")
    print(f"{'conc.':>8} " + " ".join(f"{d:>10}" for d in durations))
    for c in (0.01, 0.1, 1.0):
        print(f"{c:>6} M " + " ".join(f"{k2_min(c, T):10.0e}" for T in durations.values()))
    print("nature's floor: " + ", ".join(f"{k} {v:.0e}" for k, v in FLOOR.items()))

    Ts = np.logspace(np.log10(YR / 12), np.log10(5 * YR), 100)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for c in (0.01, 0.1, 1.0):
        ax.loglog(Ts / YR, k2_min(c, Ts), label=f"{c:g} M")
    for lab, v in FLOOR.items():
        ax.axhline(v, color="0.5", ls="--", lw=0.8)
        ax.text(0.09, v * 1.3, lab, fontsize=8, color="0.35")
    ax.set(xlabel="experiment length (years)", ylabel="smallest detectable k2 (/M/s)",
           title="What a slow experiment can see (5% seed, 0.2% ee noise)")
    ax.legend(title="concentration", fontsize=8)
    fig.tight_layout()
    fig.savefig("design.png", dpi=150)
    print("wrote design.png")


if __name__ == "__main__":
    main()
