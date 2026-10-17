"""Spec sheet for the missing self-amplifying reaction, and how known candidates score.

Spec (from the models in this folder):
  R1 self-amplification            each hand makes more of itself        (sim.py: growth term)
  R2 hands suppress each other     needed to tip to one side (Frank 1953) (sim.py: saturation term)
  R3 driven / open                 energy or feed keeps it from racemizing
  R4 water, prebiotic inputs       works from simple early-Earth molecules
  R5 fast enough                   k2 >= 10/(c*tau)                       (ocean.py time rule)
  R6 sensitive enough              follows an inherited excess >= e_need  (inherit.py)
  R7 biological product            amino acid, sugar or nucleotide (or passes its hand to one)
  R8 many molecules decide         homogeneous solution, molecule-sized voting units; no crystal,
                                   surface or aggregate carries the decision (only ~1-300 units
                                   decide in crystallizations: crystal_check.py, blind_units.py)
  R9 tips over slowly              racemic state stable at low drive, lost as the drive rises
                                   slowly (the sweep of Eq. 11; frank.py)
R8 is a gate: a bias of g ~ 1e-17 cannot steer a decision made by a few units, so candidates
are ranked by R8 first (Y, then P, then N), then by score.
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
# except where no rate was found (?), and N where a model at measured rates gives an excess that
# decays (k2 = k/c <= 0; arm_e_rates.py, rna_copying.py). R6 here = shown to work from racemic or <= ~1% ee.
# R9 is N where the racemic state is unstable from the start as run (batch, grinding), and N
# where there is no amplification (no R1), or a model at measured rates finds the racemic state
# stable at every feed, and so no instability to cross.
KEYS = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9"]


def _m(s):
    return dict(zip(KEYS, s.split()))


CANDIDATES = {
    "Soai reaction (1995; 2003)":            _m("Y P P N Y Y N P N"),  # from 5e-7 ee; organozinc, toluene; <= 3.6e9 decide (soai_check.py S4)
    "Viedma ripening, NaClO3 (2005)":        _m("Y P Y P Y Y N N N"),  # racemic -> 100%, grinding, inorganic
    "Attrition deracemization, aa deriv (2008)": _m("Y Y Y P Y Y P N N"),  # needs derivative + racemizer
    "Ghadiri peptide replicator (2001)":     _m("Y Y P N Y Y Y Y ?"),  # designed 32-mer peptides
    "Joyce RNA cross-inhibition (1984)":     _m("N Y P Y ? N Y Y N"),  # inhibition only
    "Blackmond eutectic (2006)":             _m("N N N P Y P Y N N"),  # amplifies existing ee; CHCl3
    "Magnetite crystallization (2023)":      _m("N N P Y Y Y Y N N"),  # racemic -> ~60% ee, then crystals
    # Added 2026-10-17 (abstracts checked in OpenAlex/PubMed):
    "Transamination + ligation network, fed (2025)": _m("P P Y Y N ? Y Y N"),  # Higgs & Blackmond 2025 model; Yu 2024, Deng 2024 steps; arm_e_rates.py: at Yu's rates the excess decays at 1-5 mM, catalysed rate 735-3658x short
    "Nonenzymatic RNA copying, 2-AI (2026)":  _m("P Y P Y N ? Y Y N"),  # Essex 2026 preprint: L-nt terminates D-RNA copying; no sustained replication yet; rna_copying.py: racemic stable for r < 1 (seed-following only)
    "APED peptides via NCAs (2004-2010)":     _m("P P P P ? P Y Y ?"),  # Plasson 2004 model; Danger 2010: DKP removes heterochiral faster, equilibrium favours it; Hitz & Luisi 2004: from 20% ee
}


def main():
    print("Spec: minimum rate and sensitivity per setting")
    for name, k2_min, e in spec():
        print(f"  {name:18} k2 >= {k2_min:.0e} /M/s    follows inherited excess >= {e:.0e}")
    gate = {"Y": 0, "P": 1, "N": 2, "?": 1}
    score = {n: sum({"Y": 1, "P": 0.5}.get(m[k], 0) for k in KEYS) for n, m in CANDIDATES.items()}
    ranked = sorted(CANDIDATES, key=lambda n: (gate[CANDIDATES[n]["R8"]], -score[n]))
    print(f"\n{'candidate (R8 gate first, then score)':47} " + " ".join(KEYS) + "  score")
    for name in ranked:
        print(f"{name:47} " + "  ".join(CANDIDATES[name][k] for k in KEYS) + f"   {score[name]:.1f}/9")
    return ranked

if __name__ == "__main__":
    ranked = main()
    # every crystal or surface route (R8 = N) ranks below every solution route
    assert all(CANDIDATES[n]["R8"] == "N" for n in ranked[ranked.index("Attrition deracemization, aa deriv (2008)"):])
    assert ranked[0] == "Ghadiri peptide replicator (2001)"
