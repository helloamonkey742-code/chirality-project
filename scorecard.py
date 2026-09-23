"""Spec sheet for the missing self-amplifying reaction, and how known candidates score.

Spec (from the models in this folder):
  R1 self-amplification            each hand makes more of itself        (sim.py: growth term)
  R2 hands suppress each other     needed to tip to one side (Frank 1953) (sim.py: saturation term)
  R3 driven / open                 energy or feed keeps it from racemizing
  R4 water, prebiotic inputs       works from simple early-Earth molecules
  R5 fast enough                   k2 >= 10/(c*tau)                       (ocean.py time rule)
  R6 sensitive enough              follows an inherited excess >= e_need  (inherit.py)
  R7 biological product            amino acid, sugar or nucleotide (or passes its hand to one)
"""
import numpy as np
from inherit import e_need
from ocean import YR

SETTINGS = {  # concentration M, time s, volume m^3
    "Warm little pond": (1e-2, 1e2 * YR, np.pi),
    "Lake / lagoon": (1e-4, 1e4 * YR, 1e7),
    "Open ocean": (1e-6, 1e8 * YR, 1.332e18),
}


def spec():
    rows = []
    for name, (c, tau, V) in SETTINGS.items():
        rows.append((name, 10 / (c * tau), e_need(c, V, 1e-3, tau)))
    return rows


# Filled from the literature (see README Part 6); Y = meets, P = partly, N = fails, ? = unknown.
# R5: every candidate runs in hours-days in the lab, far faster than the 3e-7 /M/s floor -> Y,
# except where no rate was found (?). R6 here = shown to work from racemic or <= ~1% ee.
def _m(s):
    return dict(zip(["R1", "R2", "R3", "R4", "R5", "R6", "R7"], s.split()))


CANDIDATES = {
    "Soai reaction (1995; 2003)":            _m("Y P P N Y Y N"),  # from 5e-7 ee; organozinc, toluene
    "Viedma ripening, NaClO3 (2005)":        _m("Y P Y P Y Y N"),  # racemic -> 100%, grinding, inorganic
    "Attrition deracemization, aa deriv (2008)": _m("Y Y Y P Y Y P"),  # needs derivative + racemizer
    "Ghadiri peptide replicator (2001)":     _m("Y Y P N Y Y Y"),  # designed 32-mer peptides
    "Joyce RNA cross-inhibition (1984)":     _m("N Y P Y ? N Y"),  # inhibition only
    "Blackmond eutectic (2006)":             _m("N N N P Y P Y"),  # amplifies existing ee; CHCl3
    "Magnetite crystallization (2023)":      _m("N N P Y Y Y Y"),  # racemic -> ~60% ee, then crystals
}


def main():
    print("Spec: minimum rate and sensitivity per setting")
    for name, k2_min, e in spec():
        print(f"  {name:18} k2 >= {k2_min:.0e} /M/s    follows inherited excess >= {e:.0e}")
    if CANDIDATES:
        keys = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]
        print(f"\n{'candidate':42} " + " ".join(keys) + "  score")
        for name, marks in CANDIDATES.items():
            score = sum({"Y": 1, "P": 0.5}.get(marks[k], 0) for k in keys)
            print(f"{name:42} " + "  ".join(marks[k] for k in keys) + f"   {score:.1f}/7")


if __name__ == "__main__":
    main()
