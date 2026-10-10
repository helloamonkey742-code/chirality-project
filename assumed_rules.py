"""Two assumed rules and a bigger weak force: do the ocean results lean on them? (2026-10-09)

uncertainty.py's moon results rest on two rules nobody has measured, and on a bias range from the literature:
(1) Top-layer share. When an ocean freezes into stacked layers, only the top layer decides, holding a share
    top = L/H' of the molecules (ASSUMED from the 3-5 layers seen in aniso3d.py). sphere_layers.py measured 0.52
    and 0.07 where the rule gives 0.25, i.e. x2.1 and x0.28. Run: x0.1, x0.28, x2.1, and the whole ocean deciding.
(2) The bar. Physics "wins" when Delta >= 2 (one hand chosen with ~98% odds). Run: Delta >= 1 and Delta >= 3.
(3) Bias size. g = 4e-19..4e-16 kT (about 0.01-10 aeV) sits at or below the bottom of the "aeV to feV" range Quack
    2022 (Chem. Sci., PMC9491092) gives for molecules of light atoms; molecules with S, Cl or Br reach ~200 aeV.
    Run: g x10, x100, and unlimited (the bias never limits). Delta is linear in g, so bias xk is the same as bar 2/k.
Every other draw and rule, and the seed, are uncertainty.py's; its headline numbers are reproduced exactly first.

PREDICTION (written 2026-10-09, before the first run; changes in points from 22/29/75):
P1. No variant moves any world by more than 5 points; Earth > Europa > Enceladus in every row.
P2. Top share x0.1: Enceladus -1 to -4, Europa -0.5 to -2.5. x0.28: Enceladus -0.5 to -2.5, Europa 0 to -1.5.
    x2.1: Enceladus +0.3 to +2, Europa 0 to +1. Whole ocean: Enceladus +1 to +3.5, Europa 0 to +1.5.
P3. Bar 1: Enceladus +0.5 to +2.5, Europa 0 to +1.5. Bar 3: Enceladus -0.3 to -1.5, Europa 0 to -1.
P4. Earth stays within 75 +/- 1 in every row.
(P1-P4 cover the six rule rows, top share and bar; P5, also written before the first run, covers the bias rows.)
P5. Bias x10 / x100 / unlimited: Enceladus 24.5-27.5 / 26-29.5 / 27-31%; Europa 29.5-31.5 / 30-32 / 30-32.5%;
    Earth within 75 +/- 1 and Earth > Europa > Enceladus in all three. Even unlimited, the moons stay near 30%,
    because the reaction must still finish and the layers must still merge in time.

RESULTS (assumed_rules.log, 2026-10-09): all five held.
 P1 held: the six rule rows move a moon by at most 2.2 points (Enceladus, share x0.1); the order never changes.
 P2 held: x0.1 -2.2 / -1.4; x0.28 -1.2 / -0.7; x2.1 +0.6 / +0.3; whole ocean +2.1 / +0.9 (Enceladus / Europa).
 P3 held: bar 1 +1.3 / +0.8; bar 3 -0.8 / -0.5.
 P4 held: Earth is 74.7% in every row, the bias rows included: the bias size never limits Earth.
 P5 held: bias x10 25.7% / 30.7%, x100 27.8% / 31.2%, unlimited 28.6% / 31.2% (Enceladus / Europa).
 So the moon results do not lean on the two unmeasured rules (at most ~2 points). A 10-100x larger weak force, as
 molecules with heavier atoms would give, adds 4-6 points to Enceladus and 2-3 to Europa; beyond that, the reaction
 finishing and the layers merging in time cap the moons near 29% and 31%.
"""
import numpy as np

import uncertainty as U

orig = U.delta_at
ROWS = [  # label, top-layer share multiplier m, Delta multiplier k (= bar 2/k)
    ("as published (uncertainty.py)", 1.0, 1.0),
    ("top-layer share x0.1", 0.1, 1.0),
    ("top-layer share x0.28 (sphere_layers run)", 0.28, 1.0),
    ("top-layer share x2.1 (sphere_layers run)", 2.1, 1.0),
    ("whole ocean decides", np.inf, 1.0),
    ("bar Delta >= 1", 1.0, 2.0),
    ("bar Delta >= 3", 1.0, 2 / 3),
    ("weak-force bias x10", 1.0, 10.0),
    ("weak-force bias x100", 1.0, 100.0),
    ("bias unlimited (never limits)", 1.0, np.inf),
]


def run(m=1.0, k=1.0):
    """uncertainty.py's three sweeps with the top-layer share times m (capped at 1) and Delta times k."""
    def d(c, g, k2, V, tau, f=1.0):
        if not np.isscalar(f):  # only the layered case passes an array: the deciding top-layer share
            f = np.minimum(1.0, m * f)
        return orig(c, g, k2, V, tau, f) * k
    U.delta_at = d
    try:
        U.RNG = np.random.default_rng(0)
        return {w: U.sweep(*a)[1].mean() for w, a in U.WORLDS.items()}
    finally:
        U.delta_at = orig


def main():
    U.RNG = np.random.default_rng(0)
    plain = {w: U.sweep(*a)[1].mean() for w, a in U.WORLDS.items()}
    res = {lab: run(m, k) for lab, m, k in ROWS}
    base = res[ROWS[0][0]]
    assert base == plain, "the wrapper changed uncertainty.py's result"
    for w in U.WORLDS:
        assert abs(base[w] - U.README_WIN_PCT[w]) < U.WIN_TOL, (w, base[w])
        for rows, i in (([r for r in ROWS if r[2] == 1.0], 1), ([r for r in ROWS if r[1] == 1.0], 2)):
            s = [res[r[0]][w] for r in sorted(rows, key=lambda r: r[i])]  # same draws: the win set can only grow
            assert np.all(np.diff(s) >= 0), (w, s)
    print(f"{'share of draws where physics wins':42}" + "".join(f"{w:>20}" for w in U.WORLDS))
    for lab, _, _ in ROWS:
        print(f"{lab:42}" + "".join(f"{res[lab][w]:11.1%} ({(res[lab][w] - base[w]) * 100:+5.1f})" for w in U.WORLDS))


if __name__ == "__main__":
    main()
