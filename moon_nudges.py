"""Local chiral pushes: how strong can they be before they, not the weak force, pick an ocean's hand? (2026-10-09)

The weak force pushes every molecule the same way. Most other chiral pushes have hands set locally, at random:
magnetite magnetized by the first molecules made on it ("avalanche magnetization", Ozturk et al. 2023, Nat. Commun.
14:6351), chiral crystal faces, a local field. If a share f of the deciding molecules form on M such sources of equal
output, with excess ee and random hands, the ocean (one well-mixed population, as in earth.py) is pushed by
    g_loc = f * ee * nu,   nu = (sum of M random signs) / M,   ee 0.1..0.6 as in earth.py
Sources of unequal output act like M_eff = (sum w)^2 / sum w^2 equal ones. Three cases, each added to the weak force
(which favours L here; its sign is unsettled, README Part 5):
  one source    M = 1: one source sets the push for the whole ocean, its hand at random (nu = +-1). The worst case.
  many sources  M_eff >> 1: nu is Gaussian with spread 1/sqrt(M), so the break-even is f*(M) = f*_G * sqrt(M) exactly,
                where f*_G is the value at spread 1. Only f*_G is simulated; the sqrt(M) columns are arithmetic.
  field         earth.py's planet-wide field acting on magnetite: nu = A * R * s (hemisphere imbalance A 0.01..0.5,
                reversal averaging R with P_rev 1e5..1e8 yr, polarity s = +-1; all as in earth.py). Of these moons
                only Ganymede has a steady field of its own (Kivelson et al. 1996, Nature 384:537), and its ocean may
                sit on high-pressure ice instead of rock (Vance et al. 2014), in which case its f is 0. The row is
                shown for every world for comparison.
Break-even f*: the share at which the weak force's lead (L wins minus D wins) falls to half its value without the
push. Every other draw is uncertainty.py's: its three worlds swept in its order from seed 0, then Ganymede
(jupiter_moons.py's 100 km ocean). Wall pinning (earth.py) is left out: in the field case it needs f above an
analytic floor of 3e-8 (Ganymede) to 7e-7 (Earth), checked below to sit far above every break-even; walls between
local sources are not modelled. This is a property of the model's push arithmetic, not a measurement of any mineral.

PREDICTION (written 2026-10-09, before the first run):
P1. With no push, every world gives its published share within 0.03 (Enceladus 22%, Europa 29%, Early Earth 75%,
    Ganymede 28.5%) and no D wins.
P2. One source: f* is 1e-16..1e-15 for Enceladus, 5e-17..5e-16 for Europa and Ganymede (within x2 of each other) and
    2e-17..8e-17 for Early Earth, whose draws are nearly all strong enough, so there f* ~ median g / median ee
    = 1.3e-17 / 0.35 = 3.6e-17. Order: Enceladus >= Europa >= Early Earth.
P3. Many sources: f*_G is 1.2-2x the one-source value in every world (median |nu| 0.67 instead of 1). For all three
    moons, f*(1e3) and f*(1e20) both come out below the rough figures given in chat on 2026-10-09 (1e-13 and 2e-5),
    which were scaled from Earth's field limit, a push already shrunk by hemisphere imbalance.
P4. Field: f* is 10-80x the one-source value in every world (A and R shrink the push); Ganymede 5e-16..4e-14; Early
    Earth 2e-16..1.5e-15 by the lead-halving rule, and within x1.3 of earth_breakeven.py's 1.9e-15 by its L < 50% rule.
"""
import numpy as np

import uncertainty as U
from breakeven import cross
from earth import logu
from jupiter_moons import args
from ocean import YR

WORLDS = {**U.WORLDS, "Ganymede": args(2634e3, 100e3)}  # Ganymede last: the three published sweeps keep their draws
GRID = 10 ** np.arange(-19, -9.99, 0.25)
KINDS = ("one source", "many sources", "field")
MS = (1e3, 1e6, 1e12, 1e20)
orig = U.delta_at


def shares(kind, f):
    """{world: [L-win share, D-win share]} with the weak force (toward L) plus a push of share f of the given kind."""
    def biased(c, g, k2, V, tau, frac=1.0):
        n = len(tau)
        rng = np.random.default_rng(1)  # identical pushes on every call: both hands, both layer sizes, every world
        ee = rng.uniform(0.1, 0.6, n)
        if kind == "one source":
            nu = rng.choice([-1, 1], n)
        elif kind == "many sources":
            nu = rng.standard_normal(n)  # spread 1; M sources: f*(M) = f*_G * sqrt(M)
        else:  # earth.py's field: A * R * s
            R = np.minimum(1, 0.53 * np.sqrt(logu(rng, 1e5, 1e8, n) * YR / np.sqrt(tau / (k2 * c))))
            nu = logu(rng, 0.01, 0.5, n) * R * rng.choice([-1, 1], n)
        return orig(c, 1.0, k2, V, tau, frac) * (g + f * ee * nu) * hand

    out = {w: [] for w in WORLDS}
    U.delta_at = biased
    try:
        for hand in (1, -1):  # count decisive L, then decisive D
            U.RNG = np.random.default_rng(0)
            for w, a in WORLDS.items():
                out[w].append(U.sweep(*a)[1].mean())
    finally:
        U.delta_at = orig
    return out


def fmt(c):
    return "   never " if c is None else f"{c[0]:8.1e}" + ("*" if c[2] > 1 else " ")


def main():
    base = shares("one source", 0.0)
    for w, (L, D) in base.items():
        assert abs(L - U.README_WIN_PCT.get(w, 0.285)) < U.WIN_TOL and D == 0, (w, L, D)  # 0.285: jupiter_moons.log
    print("weak force alone (favours L), L wins: " + ", ".join(f"{w} {L:.1%}" for w, (L, _) in base.items()))
    res = {(k, w): [] for k in KINDS for w in WORLDS}
    for k in KINDS:
        for f in GRID:
            for w, LD in shares(k, f).items():
                res[k, w].append(LD)
    res = {kw: np.array(v) for kw, v in res.items()}
    lead = {kw: (v[:, 0] - v[:, 1]) / base[kw[1]][0] for kw, v in res.items()}

    print("\nweak force's lead left (L wins minus D wins, as a share of its value with no push), whole decades of f")
    print(f"{'f':>7} | " + " | ".join(f"{k:^31}" for k in KINDS))
    print(f"{'':>7} | " + " | ".join(" ".join(f"{w.split()[0][:7]:>7}" for w in WORLDS) for _ in KINDS))
    for i in range(0, len(GRID), 4):
        print(f"{GRID[i]:7.0e} | " + " | ".join(" ".join(f"{lead[k, w][i]:7.2f}" for w in WORLDS) for k in KINDS))

    star = {kw: cross(GRID, v) for kw, v in lead.items()}
    print("\nbreak-even share f* (the lead halves; * = curve crosses 0.5 more than once). Many sources: f*_G * sqrt(M)")
    print(f"{'':18} {'one src':>9}{'f*_G':>10}" + "".join(f"{'M=' + format(M, '.0e'):>9}" for M in MS)
          + f"{'field':>10}   f*_G/one  field/one")
    for w in WORLDS:
        one, g, fd = (star[k, w] for k in KINDS)
        print(f"{w:18} {fmt(one)} {fmt(g)}" + "".join(f"{g[0] * np.sqrt(M):9.1e}" for M in MS)
              + f" {fmt(fd)}   {g[0] / one[0]:8.2f} {fd[0] / one[0]:10.1f}")

    e = cross(GRID, res["field", "Early Earth ocean"][:, 0])
    print(f"\nEarly Earth, field: L below 50% (earth_breakeven.py's rule) for f above {e[0]:.2e} "
          "(earth_breakeven.py: 1.9e-15)")
    assert abs(np.log10(e[0] / 1.9e-15)) < np.log10(2), "field case no longer matches earth_breakeven.py"
    print("pinning floor, field case (earth.py's rule at the thinnest wall: lowest D_h, k2*c = c_max):")
    for w, a in WORLDS.items():
        eps = np.sqrt(a[3][0] / a[2][1]) / (a[1] / np.pi)
        floor = 0.47 * eps * (0.01 + 3 * eps) / 0.6
        print(f"  {w:18} {floor:.1e}  ({floor / star['field', w][0]:.0e} x the field break-even)")
        assert floor > 100 * star["field", w][0], w


if __name__ == "__main__":
    main()
