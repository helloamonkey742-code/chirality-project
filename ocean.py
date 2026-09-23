"""Could the weak force pick life's hand in an icy-moon ocean?

Plugs real ocean sizes into the selection formula verified in sim.py:

    Delta = sqrt(2) * pi**0.25 * g * (k*tau)**0.25 * N**0.5,   P(favoured) = Phi(Delta)

with
    g   = dE_pv / kT                      parity-violating bias (unknown sign, swept)
    k   = k2 * c                          autocatalysis rate; dilute -> slow
    tau = time to cross the transition    capped by the ocean's age
    N   = c * V * f * 1000 * N_A          molecules in the well-mixed volume

Two things must both hold:
  size:  Delta >= 2   (favoured hand wins >= 97.7% of the time)
  time:  k*tau >= 10  (the chemistry actually finishes; formula needs a slow sweep)

Delta grows like c**0.75, so each constraint gives a minimum concentration.
Order-of-magnitude only: O(1) prefactors in g, eps and gamma are dropped.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NA = 6.022e23
YR = 3.156e7
PREF = np.sqrt(2) * np.pi**0.25
DELTA_TARGET = 2.0

# volume m^3, time available s -- sources in README
BODIES = {
    "Small lake (1 km² × 10 m)": (1e7, 1e4 * YR),
    "Enceladus ocean": (2.7e16, 1e8 * YR),   # derived from Cadek 2016 geometry; age debated 1e6-1e9 yr
    "Earth ocean": (1.332e18, 1e8 * YR),
    "Europa ocean": (3.3e18, 1e9 * YR),      # ~2.5x Earth; long-lived
}


def c_size(g, k2, V, tau, f=1.0):
    """Min concentration (M) for Delta >= target."""
    return (DELTA_TARGET / (PREF * g * (k2 * tau) ** 0.25 * (V * f * 1000 * NA) ** 0.5)) ** (4 / 3)


def c_time(k2, tau):
    """Min concentration (M) for the reaction to finish: k2*c*tau >= 10."""
    return 10 / (k2 * tau)


def delta_at(c, g, k2, V, tau, f=1.0):
    return PREF * g * (k2 * c * tau) ** 0.25 * (c * V * f * 1000 * NA) ** 0.5


def main():
    # self-check: c_size round-trips to Delta = target
    V, tau = BODIES["Enceladus ocean"]
    c = c_size(1e-17, 1e-3, V, tau)
    assert abs(delta_at(c, 1e-17, 1e-3, V, tau) - DELTA_TARGET) < 1e-9

    print("Minimum concentration (M) for the weak force to pick the hand")
    print("= max(size limit, time limit). Mixing only a fraction f of the ocean")
    print("acts like a weaker bias g*sqrt(f): 1% mixed == 10x smaller g.\n")
    for k2 in [1e-6, 1e-3, 1.0]:
        print(f"autocatalysis rate k2 = {k2:g} /M/s")
        print(f"  {'body':28} {'g=1e-18':>9} {'g=1e-17':>9} {'g=1e-16':>9} {'time limit':>11}")
        for name, (V, tau) in BODIES.items():
            cs = [c_size(g, k2, V, tau) for g in (1e-18, 1e-17, 1e-16)]
            ct = c_time(k2, tau)
            print(f"  {name:28} " + " ".join(f"{max(x, ct):9.1e}" for x in cs) + f" {ct:11.1e}")
        print()

    # figure: required concentration vs bias, k2 = 1e-3
    k2 = 1e-3
    gs = np.logspace(-19, -15, 100)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name, (V, tau) in BODIES.items():
        ax.loglog(gs, np.maximum(c_size(gs, k2, V, tau), c_time(k2, tau)), label=name)
    ax.axvspan(4e-19, 4e-16, color="0.9", zorder=0, label="modern PVED range (fJ–pJ/mol)")
    for c_ref, lab in [(1e-9, "1 nM"), (1e-6, "1 µM")]:
        ax.axhline(c_ref, color="0.5", lw=0.8, ls="--")
        ax.text(1.2e-19, c_ref * 1.4, lab, color="0.4", fontsize=8)
    ax.set_xlabel("parity-violating bias  g = ΔE/kT")
    ax.set_ylabel("minimum concentration (M)")
    ax.set_title("Where physics, not chance, picks life's hand  (k2 = 1e-3 /M/s)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("threshold.png", dpi=150)
    print("wrote threshold.png")


if __name__ == "__main__":
    main()
