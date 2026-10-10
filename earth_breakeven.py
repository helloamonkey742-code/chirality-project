"""Earth's magnetic field: how small must its share be for the weak force's hand to survive? (2026-10-09)

earth.py adds Earth's one unique chiral nudge, magnetite magnetized by the field (Ozturk 2023), to the weak force and
gets a coin flip (L 35%, D 35%, two hands pinned in 7%), because f_mag, the share of the deciding molecules made on
magnetized magnetite, was ASSUMED to be 1e-8..1e-2, and any share in that range swamps the weak force. Nobody has
measured f_mag. This fixes it at each value from 1e-20 to 1e-2 (every other draw exactly as in earth.py, weak force
favouring L) and finds where L falls below 50%: how small the field's share would have to be for the weak force to
keep deciding Earth's hand. Repeated with the reversal period P_rev fixed at 1e4..1e8 yr (the Archean rate is
unknown; 1e4 lies outside earth.py's range). This is a property of earth.py's model, not a measurement of real
magnetite. earth.py and earth.log are unchanged.

PREDICTION (written 2026-10-09, before the first run):
P1. With P_rev drawn as in earth.py, L >= 50% only for f_mag below a value between 3e-16 and 1e-14 (best guess 2e-15),
    6-8 decades below earth.py's lowest assumed f_mag (1e-8).
P2. At f_mag = 1e-14, L is 38-46%. For every f_mag >= 1e-13, L - D < 2 points.
P3. The pinned share rises monotonically with f_mag, reaching 30-50% at f_mag = 1e-2, where L and D are each 15-30%.
P4. Fixing P_rev anywhere in 1e5-1e8 yr changes the break-even f_mag by less than x2.5; at 1e4 yr it is x1.5-x6 above the 1e8 yr value.

RESULTS (earth_breakeven.log, 2026-10-09): all four held.
 P1 held: L >= 50% only for f_mag below 1.9e-15 (P_rev drawn), 6.7 decades below earth.py's lowest assumed 1e-8.
 P2 held: at 1e-14, L 42.3%; from 1e-13 up, L - D is at most 0.6 points.
 P3 held: pinned share 0 up to 3.2e-6, then 0.3% at 1e-5 rising to 39.4% at 1e-2, where L 22.4% and D 22.3%.
 P4 held: P_rev fixed at 1e5..1e8 yr gives 2.7e-15..1.7e-15 (within x1.6 of each other, x1.4 of the drawn value);
    1e4 yr gives 5.8e-15, x3.4 above the 1e8 value.
 So earth.py's coin flip does not hinge on its assumed 1e-8 lower bound: from about 1e-13 up the field already
 splits Earth evenly (L 37.3%, D 37.3%; the other 25.4% are draws where one ocean never forms). For the weak force
 to keep deciding Earth in most of input space, fewer than about 2 in 10^15 of the deciding molecules could have
 formed on magnetized magnetite. (With f_mag = 1e-30 Earth is 74.6% L, not uncertainty.py's 74.7%: earth.py draws
 Earth first from seed 0, uncertainty.py draws it third.)
"""
import numpy as np

import uncertainty as u
import earth as E
from breakeven import cross
from ocean import YR

EARTH = u.WORLDS["Early Earth ocean"]
GRID = 10 ** np.arange(-20, -1.99, 0.5)
P_FIXED = (1e4, 1e5, 1e6, 1e7, 1e8)


def run(f_mag=None, p_rev=None):
    """earth.run(1, 1, 0) (weak force toward L, plus the field), with f_mag and/or P_rev overwritten after drawing so
    every other draw stays identical. Returns the L-win share, D-win share and pinned share."""
    pinned = np.zeros(u.N, bool)

    def biased(c, g, k2, V, tau, f=1.0):
        n = len(tau)
        rng = np.random.default_rng(1)  # earth.py's draw order: P_rev, f_mag, ee, A, s
        P, fm = E.logu(rng, 1e5, 1e8, n), E.logu(rng, 1e-8, 1e-2, n)
        if p_rev is not None:
            P[:] = p_rev
        if f_mag is not None:
            fm[:] = f_mag
        ee, A, s = rng.uniform(0.1, 0.6, n), E.logu(rng, 0.01, 0.5, n), rng.choice([-1, 1], n)
        R = np.minimum(1, 0.53 * np.sqrt(P * YR / np.sqrt(tau / (k2 * c))))
        eps = np.sqrt(E.DH / (k2 * c)) / (EARTH[1] / np.pi)
        pinned[:] = (eps < 1) & (fm * ee > 0.47 * eps * (A + 3 * eps))
        return E.orig(c, 1.0, k2, V, tau, f) * np.where(pinned, 0.0, g + fm * ee * A * R * s) * hand

    out = []
    u.delta_at = biased
    try:
        for hand in (1, -1):  # count decisive L, then decisive D
            u.RNG = np.random.default_rng(0)
            out.append(u.sweep(*EARTH)[1].mean())
    finally:
        u.delta_at = E.orig
    return out[0], out[1], pinned.mean()


def main():
    L0, D0, p0 = run()
    assert [L0, D0] == E.run(1, 1, 0) and p0 == E.PINNED.mean(), "this copy no longer reproduces earth.py"
    L, D, p = run(f_mag=1e-30)
    assert abs(L - u.README_WIN_PCT["Early Earth ocean"]) < u.WIN_TOL and D == 0 and p == 0, (L, D, p)
    print(f"earth.py as published (f_mag drawn 1e-8..1e-2): L {L0:.1%}  D {D0:.1%}  pinned {p0:.1%}")
    print(f"field share ~0 (f_mag = 1e-30):                L {L:.1%}  D {D:.1%}  pinned {p:.1%}\n")
    print(f"{'f_mag':>8} {'L wins':>7} {'D wins':>7} {'no hand':>8} {'pinned':>7}  | L wins with P_rev fixed at (yr):"
          + "".join(f"{P:7.0e}" for P in P_FIXED))
    rows = []
    for fm in GRID:
        L, D, p = run(f_mag=fm)
        if fm < 7e-7:  # analytic floor: pinning needs f_mag*0.6 > 0.47*eps_min*(0.01 + 3*eps_min), f_mag > 7.3e-7
            assert p == 0, f"pinned at f_mag {fm:.1e}"
        rows.append([L, D, p] + [run(f_mag=fm, p_rev=P)[0] for P in P_FIXED])
        r = rows[-1]
        print(f"{fm:8.1e} {r[0]:7.1%} {r[1]:7.1%} {1 - r[0] - r[1]:8.1%} {r[2]:7.1%}  |"
              + " " * 33 + "".join(f"{v:7.1%}" for v in r[3:]))
    rows = np.array(rows)
    print("\nL >= 50% only for f_mag below (interpolated in log f_mag):")
    for j, lab in [(0, "P_rev drawn 1e5..1e8 yr")] + [(3 + i, f"P_rev fixed at {P:.0e} yr") for i, P in enumerate(P_FIXED)]:
        c = cross(GRID, rows[:, j])
        print(f"  {lab:24} " + ("never crosses" if c is None else f"{c[0]:.2e}" + (f"  (crosses {c[2]}x)" if c[2] > 1 else "")))


if __name__ == "__main__":
    main()
