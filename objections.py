"""Two published objections to the weak force deciding, tested in the model (2026-10-09)

1. Lente (PCCP 9:6134, 2007; review: Symmetry 2:767, 2010). In a closed pot where self-copying starts with the first
   chiral molecules, those few molecules decide, so the weak force cannot pick the hand even at ocean scale: copying
   of second order would need more than 1e50 molecules (a whole Europa ocean at 1 mol/L holds ~2e45). This project's
   chemistry is a different class (frank.py): open (fed, with outflow), the 50/50 mixture stable until a slowly rising
   drive tips it, so the whole population chooses at once and the bias adds up as sqrt(N). Lente's review notes this
   open class was computed only for 100-1000 molecules. Test: rerun frank.py's simulation from 1e4 to 1e20 molecules,
   shrinking the bias g so that g*sqrt(Omega) stays 0.3. If the sqrt(N) gain holds at large N, every size gives
   frank.py's formula value. (Past ~1e20, double precision starts to lose a bias this small; 1e20 is still 13
   decades short of an ocean, so this tests the trend, not the ocean itself.)
2. Gleiser & Walker (Int. J. Astrobiol. 11:287, 2012). Impacts or eruptions can wipe the chosen hand out and make it
   be chosen again ("punctuated chirality"), so local events, not the weak force, set the final hand. A wipe-out
   re-rolls a die that the weak force loads the same way every time, but the re-roll may be fast, while uncertainty.py
   lets the drive rise over the world's whole age (slow-passage factor (k*tau)^(1/4)).
   Test: the last choice made (a) as fast as the time rule allows, k*tau = 10, or (b) suddenly, with no slow-passage
   gain at all: each molecule tilted by the weak force and the majority wins, Delta = g*sqrt(N)/2 (Lente's equilibrium
   count). Every draw and every other rule is uncertainty.py's (mixing and healing still get the world's whole age;
   the chemistry factor F is kept, which is conservative for (b)); its headline numbers are reproduced first.

PREDICTION (written 2026-10-09, before the first run):
P1. Every size from 1e4 to 1e20 gives a favoured-hand share within 0.03 of the formula (0.737), with no trend:
    largest minus smallest share < 0.05.
P2. Control: no bias at 1e20 gives 0.50 +/- 0.03 (rounding adds no hand of its own).
P3. Fast last choice (k*tau = 10): Enceladus 13-20% (now 22%), Europa 24-28.5% (now 29%), Earth 73.5-75%.
P4. Sudden last choice (Delta = g*sqrt(N)/2): Enceladus 6-15%, Europa 17-26%, Earth 69-75%.
P5. Earth > Europa > Enceladus in both rows, and neither takes a moon to 0: the ocean's size keeps some advantage
    even with no slow-passage gain.

RESULTS (objections.log, 2026-10-09): all five held.
 P1 held: 0.733, 0.733, 0.762, 0.742, 0.753 from 1e4 to 1e20 molecules (formula 0.733-0.737); largest minus smallest
    0.029. The three largest sizes sit 0.01-0.03 above the formula (1-3 standard errors), never below.
 P2 held: no bias at 1e20 gives 0.493.
 P3 held: fast last choice: Enceladus 17.3% (-4.5), Europa 26.0% (-2.7), Earth 74.7% (0.0).
 P4 held: sudden last choice: Enceladus 11.2% (-10.6), Europa 20.7% (-7.9), Earth 73.2% (-1.5).
 P5 held: Earth > Europa > Enceladus in both rows; no moon reaches 0.
 So (1) in the open class the favoured-hand share depends only on g*sqrt(N) across 16 decades of size: the bias keeps
 adding up as sqrt(N). Lente's limit is for the closed class and does not carry over; which class early chemistry
 belonged to is an assumption, not something this tests. (2) If every wipe-out ended in a sudden re-choice, the moons
 would fall from 22% / 29% to 11% / 21% and Earth from 75% to 73%: wipe-outs cost the moons a third to a half of
 their chance but do not remove it, as long as each ocean still mixes as one. Among the published wins the slow rise
 multiplies Delta by a median x15 (Enceladus), x18 (Europa), x36 (Earth).
"""
import numpy as np

import frank
import uncertainty as U
from ocean import PREF

SIZES = (1e4, 1e8, 1e12, 1e16, 1e20)
orig = U.delta_at


def sizes():
    """frank.run at each size with g*sqrt(Omega) = 0.3 and tau = 800 (frank.log's middle row is the 1e4 case)."""
    rows = [(om, 0.3 / np.sqrt(om), *frank.run(0.3 / np.sqrt(om), om, 800, seed=i)) for i, om in enumerate(SIZES)]
    return [(om, g, frank.predicted(g, om, cross), p) for om, g, p, cross in rows], frank.run(0.0, 1e20, 800, seed=9)[0]


def sweeps(prefactor=None):
    """uncertainty.py's three sweeps; prefactor replaces sqrt2*pi^(1/4)*(k2*c*tau)^(1/4) in Delta (None: as published)."""
    def d(c, g, k2, V, tau, f=1.0):
        return orig(c, g, k2, V, tau, f) * prefactor / (PREF * (k2 * c * tau) ** 0.25)
    U.delta_at = orig if prefactor is None else d
    try:
        U.RNG = np.random.default_rng(0)
        return {w: U.sweep(*a)[:2] for w, a in U.WORLDS.items()}
    finally:
        U.delta_at = orig


def main():
    rows, p0 = sizes()
    print("1. Lente: does the open class keep its sqrt(N) gain at large N? (frank.py, g*sqrt(Omega) = 0.3, tau = 800)")
    print(f"{'molecules':>10} {'bias g':>8} {'formula':>8} {'simulated':>9}")
    for om, g, pred, p in rows:
        print(f"{om:10.0e} {g:8.1e} {pred:8.3f} {p:9.3f}")
    print(f"{1e20:10.0e} {0.0:8.1e} {0.5:8.3f} {p0:9.3f}   (control: no bias)")
    assert abs(rows[0][3] - 0.733) < 6e-4, "no longer reproduces frank.log's 0.733"

    res = {"as published (uncertainty.py)": sweeps(), "fast last choice (k*tau = 10)": sweeps(PREF * 10 ** 0.25),
           "sudden last choice (g*sqrt(N)/2)": sweeps(0.5)}
    base, fast, sudden = res.values()
    for w in U.WORLDS:
        assert abs(base[w][1].mean() - U.README_WIN_PCT[w]) < U.WIN_TOL, (w, base[w][1].mean())
        # every win needs k*tau >= 10, so a faster choice can only shrink Delta: each win set holds the next
        assert not (fast[w][1] & ~base[w][1]).any() and not (sudden[w][1] & ~fast[w][1]).any(), w
    print(f"\n2. Gleiser & Walker: the last choice made fast\n{'share of draws where physics wins':34}"
          + "".join(f"{w:>20}" for w in U.WORLDS))
    for lab, r in res.items():
        print(f"{lab:34}" + "".join(f"{r[w][1].mean():11.1%} ({(r[w][1].mean() - base[w][1].mean()) * 100:+5.1f})"
                                    for w in U.WORLDS))
    gain = {w: np.median(((x["rate k2"] * x["concentration"] * x["time"])[win] / 10) ** 0.25) for w, (x, win) in base.items()}
    print("median slow-passage gain among published wins, (k*tau/10)^(1/4): "
          + ", ".join(f"{w} x{v:.0f}" for w, v in gain.items()))


if __name__ == "__main__":
    main()
