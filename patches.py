"""One ocean or a patchwork -- applied to Enceladus and Europa.

Uses the laws checked in domains.py (constants C, K measured there):
  patch size at the transition   l = C * sqrt(D/k) * (k*tau)**0.25
  molecules deciding one patch   N = K * c * (patch volume) * 1000 * N_A
  wall speed (favoured invades)  v = (3/sqrt2) * g * sqrt(D*k)
with k = k2*c. If l >= ocean size the whole ocean decides together (ocean.py case).

Closed-form result: a "one ocean" window (coherent AND finishes in time) exists only if
    D >= sqrt(10) * L**2 / (C**2 * tau)
"""
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ocean import NA, YR, PREF

C, K = 12.8, 0.61          # from domains.py (1D; 3D prefactors may differ by O(1))

# ocean volume m^3, horizontal size L m (half circumference), thickness H m, age s
MOONS = {
    "Enceladus": (2.7e16, np.pi * 232e3, 40e3, 1e8 * YR),
    "Europa": (3.3e18, np.pi * 1460e3, 100e3, 1e9 * YR),
}
D_MOLECULAR, D_EDDY = 1e-9, 5e-5   # m^2/s; eddy value used in Zeng & Jansen 2021 Enceladus simulations
D_RANGE = (1e-10, 1e-3)            # their plausible range: very poorly constrained


def favoured_fraction(c, D, g, k2, V, L, H, tau):
    """Expected fraction of the ocean on the favoured hand; NaN if chemistry can't finish."""
    k = k2 * c
    l = C * np.sqrt(D / k) * (k * tau) ** 0.25
    vol = np.where(l >= L, V, np.where(l >= H, l**2 * H, l**3))
    n = np.where(l >= L, 1.0, K) * c * vol * 1000 * NA
    frac = norm.cdf(PREF * g * (k * tau) ** 0.25 * np.sqrt(n))
    return np.where(k * tau >= 10, frac, np.nan), l


def d_min(L, tau):
    return np.sqrt(10) * L**2 / (C**2 * tau)


def main():
    g, k2 = 1e-17, 1e-3
    for name, (V, L, H, tau) in MOONS.items():
        print(f"{name}: one-ocean window needs mixing D >= {d_min(L, tau):.1e} m^2/s "
              f"(molecular {D_MOLECULAR:.0e}, modelled eddy {D_EDDY:.0e})")
        for D in (D_MOLECULAR, D_EDDY):
            for c in (1e-12, 1e-9, 1e-6):
                f, l = favoured_fraction(c, D, g, k2, V, L, H, tau)
                wall_T = L / (3 / np.sqrt(2) * g * np.sqrt(D * k2 * c)) / YR
                print(f"   D={D:.0e} c={c:.0e} M: patch {float(l):9.2e} m  favoured {float(f):.3f}"
                      f"  wall takeover {wall_T:.0e} yr")
    # sanity: very strong mixing + enough molecules -> whole ocean, physics wins
    V, L, H, tau = MOONS["Europa"]
    assert favoured_fraction(1e-9, 1.0, g, k2, V, L, H, tau)[0] > 0.97

    # regime map for Enceladus
    V, L, H, tau = MOONS["Enceladus"]
    Ds, cs = np.logspace(-10, 0, 200), np.logspace(-15, -3, 200)
    DD, CC = np.meshgrid(Ds, cs)
    F, Lp = favoured_fraction(CC, DD, g, k2, V, L, H, tau)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    m = ax.pcolormesh(DD, CC, F, cmap="viridis", vmin=0.5, vmax=1.0, shading="auto")
    fig.colorbar(m, label="fraction of ocean on the favoured hand")
    ax.contour(DD, CC, Lp >= L, levels=[0.5], colors="w", linewidths=1.2)
    ax.contour(DD, CC, np.nan_to_num(F), levels=[0.977], colors="orange", linewidths=1.2)
    for d, lab in [(D_MOLECULAR, "molecular diffusion"), (D_EDDY, "modelled mixing (Zeng & Jansen)")]:
        ax.axvline(d, color="0.8", ls="--", lw=0.8)
        ax.text(d * 1.3, 2e-15, lab, color="0.45", fontsize=8, rotation=90, va="bottom")
    ax.text(4e-9, 2e-13, "too dilute: never finishes", color="0.3", fontsize=8)
    ax.set(xscale="log", yscale="log", xlabel="mixing D (m²/s)", ylabel="concentration (M)",
           title="Enceladus: one ocean or a patchwork?  (g=1e-17, k2=1e-3 /M/s)")
    ax.text(4e-9, 2e-15, "right of white line: whole ocean acts as one\ninside orange line: favoured hand ≥ 97.7%",
            fontsize=8, color="0.2")
    fig.tight_layout()
    fig.savefig("patches.png", dpi=150)
    print("wrote patches.png")


if __name__ == "__main__":
    main()
