"""How bad is moon_nudges.py's result? Real counts of random-handed patches, and the one push with a shared hand. (2026-10-10)

moon_nudges.py found that a local chiral push whose hand is set at random outvotes the weak force unless fewer than
~5 in 10^17 of the deciding molecules form on it (one source), the limit growing as sqrt(M) with M independent
sources. This script asks how many independent sources a real seafloor has, and how strong the one push that shares
a hand over a whole hemisphere (a magnetic field acting on magnetite) can be on each moon. Pure arithmetic on
moon_nudges.log; no new simulation.

Part A, random hands. If the hand changes every d (the size over which one hand holds: a crystal, a magnetite grain or
film), a rock-water contact area S holds M = S/d^2 independent patches. Unequal output (lognormal, spread sigma) acts
like M * exp(-sigma^2) equal ones. moon_nudges.py's f*(M) = f*_G * sqrt(M) then gives the share of the deciding
molecules that may form on such surfaces before the weak force's lead halves:
    f_allowed = f*_G * sqrt(S) * exp(-sigma^2 / 2) / d
and, turned round, the largest patch size that keeps a share f safe: d_c = f*_G * sqrt(S) * exp(-sigma^2 / 2) / f.
S = (seafloor area) * (depth h of rock that water reaches) * 6(1 - phi)/a (grains of size a, porosity phi). Ranges
(ASSUMED, chosen wide): h 1 m (a sediment skin) to 10 km (fractured rock), a 1 mm to 10 um, phi 0.3 to 0.2.
Enceladus instead uses its whole porous core (radius 190 km, porosity 0.2-0.3; Choblet et al. 2017, Nat. Astron.
1:841). Ganymede's rock is at r ~ 1800 km and its ocean may not touch it at all (Vance et al. 2014).

Part B, a shared hand. A uniform or dipole field magnetizes magnetite the same way over a hemisphere. earth.py assumes
the grains are fully aligned (a = 1). In weaker fields alignment falls: for single-domain grains a = tanh(x),
x = Ms V B / kT (Ms = 4.8e5 A/m, T = 280 K). Steady fields at the rock: Enceladus, Saturn's dipole at 3.95 Saturn radii
(2.1e-5 T / 3.95^3); Europa, Jupiter's at 9.38 Jupiter radii (4.2e-4 T / 9.38^3); Ganymede, its own (719 nT at the
surface equator, Kivelson et al. 1996, Nature 384:537), x (2634/1800)^3 at its rock, x2 at the poles; Earth 25-65 uT.
Break-even for the moons: worst case f*_one / (a * A * R) with A = 0.5 (earth.py's largest hemisphere imbalance) and
R = 1 (Jupiter's and Saturn's fields are not known to reverse); and moon_nudges.py's field value / a (earth.py's A, R
spread). Grain volumes: 30 nm and 100 nm cubes (equant single domain) and 1e-20 m^3 (elongated single domain, the
largest moment a single-domain grain can carry). Compared with Early Earth's field break-even, 6.0e-16 (same rule).

Not modelled: neighbouring patches sharing a hand (that raises d to the hand's correlation length, which is the real
unknown); a magnetite grain with both its top and bottom faces exposed presenting both hands (Ozturk et al. used a
film on a substrate); low temperature (Titan) raising g. This is the model's arithmetic with assumed geology, not a
measurement of any seafloor.

PREDICTION (written 2026-10-10, before the first run):
P1. Read back from moon_nudges.log, the one-source break-evens are 7.3e-17, 5.2e-17, 3.9e-17, 4.7e-17 and f*_G
    1.3e-16, 9.4e-17, 7.4e-17, 8.7e-17 (Enceladus, Europa, Early Earth, Ganymede).
P2. With all deciding molecules on such surfaces (f = 1) and equal output, the critical patch size d_c lies between
    0.01 um (smallest S) and 100 um (largest S) in every world; Enceladus, with its porous core, has the narrowest
    range, 1-20 um. sigma = 2 shrinks d_c by e^2 = 7.4x.
P3. At grain-sized patches (d = 100 um, equal output) the allowed share runs from below 1e-3 (smallest S) to above 0.1
    (largest S) on every world: the model alone cannot say whether random-handed surfaces are a problem; the answer
    depends on how much rock the water touches and how far one hand extends.
P4. Field, worst case (A = 0.5, R = 1, 100 nm grains): Enceladus 2e-15..5e-15 and Europa 1e-15..3e-15, i.e. 2-8x
    Early Earth's 6.0e-16; Ganymede 1.5e-16..4e-16, below Earth's. With earth.py's A, R spread, Enceladus and Europa
    allow 1e-14..5e-14. Elongated grains (1e-20 m^3) bring every moon to within x2 of Earth or below it.

RESULTS (nudge_realism.log, 2026-10-10): P1 and P4 held; P2 held but for Early Earth's top end; P3 failed in part.
 P1 held: the log reads back exactly (asserted).
 P2 held but for one end: with f = 1 and equal output d_c is Enceladus 1.4-15 um, Europa 0.031-33 um, Ganymede
    0.036-38 um, Early Earth 0.11-120 um (above the predicted 100 um). sigma = 2 shrinks each by 7.4x. At f = 0.01,
    d_c is 100x larger: 3 um to 12 mm.
 P3 failed in part: at d = 100 um the smallest-S share is 3e-4 (Europa), 4e-4 (Ganymede), 1e-3 (Early Earth) and
    1e-2 (Enceladus, whose porous core gives it far more surface), so Enceladus and Earth are not below 1e-3; the
    largest-S share is 0.2-1 everywhere. The point held: the model does not decide it; geology does.
 P4 held: worst case (A = 0.5, R = 1, 100 nm grains) Enceladus 3.5e-15 (5.8x Earth's 6.0e-16), Europa 1.6e-15
    (2.7x), Ganymede 1.9e-16..3.4e-16 (0.3-0.6x); with earth.py's A, R spread 2.8e-14 and 2.1e-14. Elongated grains:
    3.7e-16, 1.9e-16, 9.4e-17, all at or below Earth.
 So random-handed surfaces are not fatal: they average out if one hand holds over less than about 1-40 um (all
 chemistry on surfaces) or 0.1-4 mm (1% on surfaces) on the moons; what is unknown is how far one hand extends. The
 push with a shared hand is the real limit: below ~1e-15 of the deciding molecules on field-aligned magnetite on
 Enceladus and Europa (a few times looser than Earth, because their parent planets' fields are ~100x weaker),
 tighter than Earth on Ganymede, and no better than Earth if the magnetite is elongated single-domain grains.
"""
import re

import numpy as np

WORLDS = ("Enceladus", "Europa", "Early Earth ocean", "Ganymede")
SEAFLOOR_R = {"Europa": 1440e3, "Early Earth ocean": 6371e3, "Ganymede": 1800e3}  # radius of the rock surface (m)
D = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)  # patch sizes (m)
KB, MS, T = 1.380649e-23, 4.8e5, 280.0
VOLS = (27e-24, 1e-21, 1e-20)  # 30 nm cube, 100 nm cube, elongated single domain (m^3)
FIELDS = {  # steady field at the rock (T): (low, high)
    "Enceladus": (2.1e-5 / 3.95**3,) * 2,
    "Europa": (4.2e-4 / 9.38**3,) * 2,
    "Ganymede": (719e-9 * (2634 / 1800) ** 3, 2 * 719e-9 * (2634 / 1800) ** 3),
    "Early Earth ocean": (25e-6, 65e-6),
}
EARTH_FIELD = 6.0e-16


def read_log(path="moon_nudges.log"):
    """{world: (one source, f*_G, field)} from moon_nudges.log's break-even table."""
    out = {}
    for line in open(path):
        for w in WORLDS:
            if line.startswith(w) and re.match(r"\s+\d", line[len(w):]):
                v = [float(x) for x in line[len(w):].split()]
                out[w] = (v[0], v[1], v[6])
    return out


def area(w):
    """(low, high) rock-water contact area S in m^2."""
    sv = lambda phi, a: 6 * (1 - phi) / a  # grain surface per m^3 of rock
    if w == "Enceladus":
        V = 4 / 3 * np.pi * 190e3**3
        return V * sv(0.3, 1e-3), V * sv(0.2, 1e-5)
    A = 4 * np.pi * SEAFLOOR_R[w] ** 2
    return A * 1 * sv(0.3, 1e-3), A * 1e4 * sv(0.2, 1e-5)


def main():
    B = read_log()
    assert [B[w][0] for w in WORLDS] == [7.3e-17, 5.2e-17, 3.9e-17, 4.7e-17], B
    assert [B[w][1] for w in WORLDS] == [1.3e-16, 9.4e-17, 7.4e-17, 8.7e-17], B

    print("PART A: random hands. S = rock-water contact area; d_c = largest single-hand patch that keeps share f safe")
    print(f"{'':18} {'S low':>9} {'S high':>9} | d_c, f=1 (um)  sigma=0     sigma=2 | d_c, f=0.01 (mm), sigma=0")
    for w in WORLDS:
        fg, S = B[w][1], area(w)
        dc = lambda s, f, sig=0: fg * np.sqrt(s) * np.exp(-sig**2 / 2) / f
        assert abs(fg * np.sqrt(S[0]) / dc(S[0], 1) - 1) < 1e-12  # f_allowed(d_c) = 1
        print(f"{w:18} {S[0]:9.1e} {S[1]:9.1e} | {dc(S[0], 1)*1e6:9.2g} .. {dc(S[1], 1)*1e6:<9.2g}"
              f"{dc(S[0], 1, 2)*1e6:7.2g} .. {dc(S[1], 1, 2)*1e6:<7.2g} | {dc(S[0], .01)*1e3:8.2g} .. {dc(S[1], .01)*1e3:.2g}")
    print("\nallowed share f (capped at 1) by patch size d, equal output; low S / high S")
    print(f"{'':18} " + " ".join(f"{f'd={d*1e6:g} um':>17}" for d in D))
    for w in WORLDS:
        fg, S = B[w][1], area(w)
        cell = lambda d: "/".join(f"{min(1, fg * np.sqrt(s) / d):.0e}" for s in S)
        print(f"{w:18} " + " ".join(f"{cell(d):>17}" for d in D))

    print("\nPART B: a shared hand (field on magnetite). a = tanh(Ms V B / kT); earth.py assumes a = 1")
    print(f"{'':18} {'B (nT)':>13} | a: 30 nm   100 nm   elong | worst f* (A=.5, R=1): 100 nm   elong "
          f"| x Earth 6.0e-16 | f*/a, A,R spread: 100 nm")
    for w in ("Enceladus", "Europa", "Ganymede"):
        one, _, field = B[w]
        a = [[np.tanh(MS * v * b / (KB * T)) for b in FIELDS[w]] for v in VOLS]  # [vol][low B, high B]
        assert abs(a[0][0] / (MS * VOLS[0] * FIELDS[w][0] / (KB * T)) - 1) < 0.01  # weak field: a ~ x
        worst = [[one / (x * 0.5) for x in row] for row in a]
        rng = lambda p: f"{min(p):.1e}" if np.isclose(min(p), max(p), atol=0) else f"{min(p):.1e}..{max(p):.1e}"
        print(f"{w:18} {FIELDS[w][0]*1e9:5.0f}..{FIELDS[w][1]*1e9:<5.0f} | "
              + " ".join(f"{rng(r):>7}" if np.isclose(min(r), max(r), atol=0) else f"{min(r):.2f}-{max(r):.2f}" for r in a)
              + f" | {rng(worst[1]):>18} {rng(worst[2]):>17} | {rng([x / EARTH_FIELD for x in worst[1]]):>15}"
              + f" | {rng([field / x for x in a[1]])}")
    aE = [np.tanh(MS * VOLS[1] * b / (KB * T)) for b in FIELDS["Early Earth ocean"]]
    print(f"{'Early Earth ocean':18} {25000:5.0f}..{65000:<5.0f} | 100 nm a = {min(aE):.3f}..{max(aE):.3f} "
          f"(earth.py's a = 1 holds for 100 nm grains); field break-even {B['Early Earth ocean'][2]:.1e}")


if __name__ == "__main__":
    main()
