"""Spec sheet for the missing self-amplifying reaction, and how known candidates score.

Spec (from the models in this folder):
  R1 self-amplification            each hand makes more of itself        (sim.py: growth term)
  R2 hands suppress each other     needed to tip to one side (Frank 1953) (sim.py: saturation term)
  R3 driven / open                 energy or feed keeps it from racemizing
  R4 water, prebiotic inputs       works from simple early-Earth molecules
  R5 fast enough                   k2 >= 10/(c*tau)                       (ocean.py time rule)
  R6 sensitive enough              follows an excess seeded in the pool down to e_need (inherit.py)
                                   (Delta >= 2; passing a feed excess through does not count)
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
# decays (k2 = k/c <= 0; arm_e_rates.py, rna_copying.py). R6: an excess seeded in the pool at
# the transition (alpha = e0, inherit.py) is followed down to e_need (Delta >= 2). Following a feed excess does not
# count: a stable racemic state passes it through, and Delta >= 2 then holds at a catalysed share of 1.7e-17
# (arm_e_r6.log:22-24). Literature rows were lettered under the earlier reading ("shown to work from racemic or
# <= ~1% ee"); re-checked under this rule (lit_r6.py, journals/new_analyses_lit_r6.md), all eight are N: a hand
# decided among N_d units cannot follow a seed below ~2/sqrt(N_d), so following e_need needs N_d >= 1.2e23 (pond),
# and every literature route is decided by crystals, a surface or <= 1e23 molecules, or its seed decays or does not
# set the hand.
# R9 is N where the racemic state is unstable from the start as run (batch, grinding), and N
# where there is no amplification (no R1), or a model at measured rates finds the racemic state
# stable at every feed, and so no instability to cross.
KEYS = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9"]


def _m(s):
    return dict(zip(KEYS, s.split()))


CANDIDATES = {
    "Soai reaction (1995; 2003)":            _m("Y P P N Y N N P N"),  # from 5e-7 ee; organozinc, toluene; <= 3.6e9 decide (soai_check.py S4); R6 = N (lit_r6): far below 1.2e23
    "Viedma ripening, NaClO3 (2005)":        _m("Y P Y P Y N N N N"),  # racemic -> 100%, grinding, inorganic; R6 = N (lit_r6): crystals decide
    "Attrition deracemization, aa deriv (2008)": _m("Y Y Y P Y N P N N"),  # needs derivative + racemizer; R6 = N (lit_r6): crystals decide (Noorduin 2008 JACS seeds 2-3% ee; every racemic start went to R, put down to a trace chiral impurity, and 0.1 mol% phenylglycine sets the hand: Noorduin thesis 2010 ch. 2)
    "Ghadiri peptide replicator (2001)":     _m("Y Y P N Y N Y Y N"),  # designed 32-mer peptides; R6 = N (lit_r6): no run follows a seed or amplifies ee; R9 = N inferred, not modelled: Saghatelian 2001 own-hand-only catalysis, no D/L cross-inhibition (no growth order given), so racemic stays stable (as arm G, r < 1)
    "Joyce RNA cross-inhibition (1984)":     _m("N Y P Y Y N Y Y N"),  # inhibition only; R6 = N (lit_r6): no self-copying, nothing to seed; R5 joyce_r5.py: Table 1, 0 C, 7 d, k2(D) >= 1.7e-5, k2(DL) >= 1.7e-6 /M/s (lower bounds)
    "Blackmond eutectic (2006)":             _m("N N N P Y N Y N N"),  # amplifies existing ee; proline in DMSO or water (Klussmann 2006 SI); R6 = N (lit_r6): crystals carry the excess, overall ee only partitioned (Klussmann 2006 SI)
    "Magnetite crystallization (2023)":      _m("N N P Y Y N Y N N"),  # racemic -> up to ~60% ee (best run; repeats ~35%), then crystals; R6 = N (lit_r6): the field-set magnetization sets the hand
    # Added 2026-10-17 (abstracts checked in OpenAlex/PubMed):
    "Transamination + ligation network, fed (2025)": _m("P P Y Y N N Y Y N"),  # Higgs & Blackmond 2025 model; Yu 2024, Deng 2024 steps; arm_e_rates.py: at Yu's rates the excess decays at 1-5 mM, catalysed rate 735-3658x short; R6 = N (arm_e_r6.py): a seed in the pool is followed only down to 44 x e_need (pond, lake; arm_e_r6.log:67) and a feed excess is passed through (monomer gain <= 2.45), so 1% output needs 0.2-0.4% (arm_e_r6.log:34-36, 55-57); its feed Delta >= 2 "Y" is a counting limit, met at a catalysed share of 1.7e-17 (arm_e_r6.log:22-24); at 1.5 x the crossing multiplier F_c (369, 3.9e4, 3.0e7 x measured) a feed excess is followed and a seed is not (1.25-1.85 x e_need; arm_e_r6.log:99-114)
    "Nonenzymatic RNA copying, 2-AI (2026)":  _m("P Y P Y N N Y Y N"),  # Essex 2026 preprint: L-nt terminates D-RNA copying; no sustained replication yet; rna_copying.py: racemic stable for r < 1 (seed-following only); Essex rates give r = 0.40 (arm_g_rates.py); gel fit r ~0.56 (fig_s8_r_fit.py), shared activator does not tip it at p <= 1; R6 = N: a template seed is lost in 54 of 54 cells (arm_g_r6.py; arm_g_r6.log:46; |lambda| 0.979-0.9995 per cycle at 20 mM, rna_seed_scan.log:105); a 1% feed excess is passed through, to 3.4-79% at 20 mM (arm_g_r6.log:129) against a gain of 3.0 with no copying (arm_g_r6.log:142-147), which R6 does not count; rna_seed_scan.py's all-Y (feed excess, Delta >= 2; rna_seed_scan.log:115) is a counting limit (the same test passes arm E at a catalysed share of 1.7e-17, arm_e_r6.log:22-24; registered letters were N, N, P, rna_seed_scan.log:111)
    "APED peptides via NCAs (2004-2010)":     _m("P P P P ? N Y Y P"),  # APED mechanism of Plasson 2004; R6 = N (lit_r6.log:4-8): racemic stable at every setting's ctot, so a seed decays; Danger 2010: DKP removes heterochiral faster, equilibrium favours it; Hitz & Luisi 2004: from 20% ee; aped.py (Plasson & Bersini 2009 reversible model with racemization, not PNAS 2004): R9 = P, racemic lost only with Plasson's signs (0/16 Danger-like sets) and ctot >= ~1 M; PNAS 2004 model (audit/aped_audit.py) has no racemization and no lower drive threshold, lost above ctot 66.9 mM at its Fig. 5 rates; R5 ? (aped_r5.py: NCA activation rate only bracketed, plain-COS end points 1.7e-7 to 1.9e-5 /s against a pond need of 5.0e-7; lake and ocean pass at the lab COS dose; in the aped.py model the difference mode is also stable below a/kH ~ 2e-3 to 1e-2, so sigma = s* a is not valid there, and the letters hold only under aped.py's time rescaling, aped_r5.log:38)
}


def main():
    print("Spec: minimum rate and sensitivity per setting")
    for name, k2_min, e in spec():
        print(f"  {name:18} k2 >= {k2_min:.0e} /M/s    follows a pool seed down to {e:.0e}")
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
    # Ghadiri and Joyce tie for first at 5.5 (lit_r6.log:33); the order between them is dict order only
    assert set(ranked[:2]) == {"Ghadiri peptide replicator (2001)", "Joyce RNA cross-inhibition (1984)"}
