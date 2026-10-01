"""Sketch for the Methods section: the landscape the chiral excess alpha moves in, U = -lam*a^2/2 + a^4/4 - g*a.
Before the transition (lam < 0) there is one valley at 50/50; after it (lam > 0) two valleys, one per hand.
The weak-force bias g tilts the landscape very slightly; noise (~1/sqrt(N)) kicks the ball around.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a = np.linspace(-1.6, 1.6, 400)
fig, axs = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
for ax, lam, title in zip(axs, (-1, 0, 1), ("before: one valley (50/50)", "at the transition: flat", "after: two valleys, one per hand")):
    ax.plot(a, -lam * a**2 / 2 + a**4 / 4 - 0.08 * a, "k")
    ax.set(title=title, xlabel="excess α  (−1 all D … +1 all L)", ylim=(-0.5, 1.2), yticks=[])
    ax.annotate("", xy=(0.25, -0.05 if lam >= 0 else 0.25), xytext=(-0.25, -0.05 if lam >= 0 else 0.25),
                arrowprops=dict(arrowstyle="->", color="tab:blue")) if lam == 0 else None
axs[0].set_ylabel("landscape U(α)")
axs[1].text(0, 0.35, "bias g tilts it\n(exaggerated here);\nnoise ~ 1/√N shakes it", ha="center", fontsize=8, color="tab:blue")
fig.tight_layout()
fig.savefig("sketch.png", dpi=150)
print("wrote sketch.png")
