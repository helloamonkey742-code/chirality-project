# The slow experiment: searching for self-amplifying handedness at nature's speed

## Aim
Test whether simple early-Earth chemistry in water can **amplify a small excess of one hand**
when given months to years instead of days. The spec sheet (README Part 6) shows that
speed is never the problem for nature, yet lab searches usually stop after days.

## Hypothesis (stated before any data)
At least one prebiotic, water-based peptide-forming system shows **seed-following growth of
enantiomeric excess** (ee) with an effective rate k2 ≥ 10⁻⁷ /M/s, which is slow enough to be missed
by short experiments but fast enough to matter in a pond (floor: 3×10⁻⁷ /M/s).

## Why a *seeded* design
A reaction that only copies itself keeps ee constant, because both hands grow equally. ee *grows*
only if the hands also **suppress each other**. So watching ee grow from a small seed tests
self-copying and mutual suppression together (spec R1 + R2), which is exactly the missing
combination. Seeding is also far faster than waiting for symmetry to break from pure 50/50.
From 50/50, noise in a 50 µL, 1 M vial would need about 18 e-foldings of growth to become visible, versus
under 1 from a 5% seed.

## How long and how concentrated (`design.py`)
Smallest detectable rate, with a 5% seed and 0.2% ee measurement noise (assumed; confirm
with the lab's method):

| concentration | 1 month | 6 months | 1 year | 3 years |
|---|---|---|---|---|
| 0.01 M | 6×10⁻⁶ | 1×10⁻⁶ | 5×10⁻⁷ | 2×10⁻⁷ |
| **0.1 M** | 6×10⁻⁷ | **1×10⁻⁷** ✓ pond | 5×10⁻⁸ | 2×10⁻⁸ |
| 1 M | 6×10⁻⁸ | 1×10⁻⁸ | 5×10⁻⁹ | **2×10⁻⁹** ✓ ocean |

![design](design.png)

**Choice:** run at **0.1 M for 12 months**. That covers the pond floor with margin by month 2.
A 1 M sub-arm, run for 1–3 years where solubility allows, reaches the ocean floor.

## Arms (chemistries to test)
Chosen from the scorecard for prebiotic inputs (R4) and biological products (R7):

| arm | system | why | source |
|---|---|---|---|
| A | amino acids + carbonyl sulfide (volcanic gas) → peptides, in water | prebiotic peptide formation; peptides can prefer same-hand partners | Leman, Orgel & Ghadiri 2004 |
| B | amino acids + hydroxy acids, **wet-dry cycling** (e.g. daily) | pond-like driving; makes polypeptides | Forsythe et al. 2015 |
| C | prebiotic cysteine peptides that **catalyse peptide joining** in neutral water | the closest known thing to self-copying in prebiotic water | Foden et al. 2020 |
| D | RNA precursor (RAO) on **magnetite**, then crystallization | non-autocatalytic route that already breaks symmetry | Ozturk et al. 2023 |
| + | **Positive control:** Viedma grinding of a known conglomerate | proves the pipeline detects real amplification | Viedma 2005 |

Arm C is the priority: it's the only one with a catalytic loop in prebiotic water.

## Conditions per arm
For each arm, three seeds × 3 replicates:
- **+5% L seed**
- **+5% D seed** (the mirror control, the most important control)
- **0% (racemic)**

Plus three negative controls (3 replicates each): no activator / no cycling; sterile-filtered;
and the seed amino acid alone (to measure plain racemization).
That's 9 + 9 = 18 vials per arm, plus the positive control. Temperature 25 °C, with an optional
40 °C set to speed things up. Racemization is then negligible over a year for most amino
acids, and the "seed alone" vials measure it directly.

## Sampling
Days 0, 3, 7, 14, 30, 60, 90, 120, 180, 270, 365 (roughly evenly spaced on a log scale).
Measure ee of the monomers (and of short peptides where possible) by **chiral GC-MS or HPLC**
on derivatized samples. Keep analysts **blind** to vial labels.

## Decision rule (fixed in advance)
Call an arm an **amplifier** only if **all four** hold:
1. |ee| rises by > 0.85% absolute over its start in both seeded sets (3σ√2 at σ = 0.2%).
2. The **L-seeded and D-seeded vials change by equal and opposite amounts** (within noise).
3. The racemic vials stay at 0 within noise.
4. Sterile and no-activator controls show no change.

Then fit `ee(t) = e0·exp(k·t)` to get k, and k2 = k / c. Compare k2 with the floors
(3×10⁻⁷ pond, 3×10⁻⁹ ocean).

## Contamination: the main way this goes wrong
Biology is left-handed, so any microbe, skin trace or reagent impurity adds **L**.
That would look like "amplification" in the L-seeded vials. Safeguards:
- **The D-seed mirror arm:** real amplification is symmetric; contamination pushes both sets toward L.
- Sterile filtration, sealed vials, and a vial with no amino acid (it should show zero amino acids).
- If only the L arm "amplifies", report it as contamination, not a result.

## What each outcome means
| result | meaning for the project |
|---|---|
| An arm amplifies at k2 ≥ pond floor | First prebiotic, water-based amplifier; pond seeding by meteorites becomes a complete story |
| An arm amplifies but slower than the floor | Real, but too slow to matter in a pond; ocean-scale only |
| Nothing amplifies in 1 year at 0.1 M | Rules out these chemistries down to k2 ≈ 5×10⁻⁸. That's still a publishable upper limit, and it tightens the spec |
| Only the positive control amplifies | Pipeline works; the gap in the scorecard is real |

## Limits and practical notes
- This needs a chemistry lab with chiral GC-MS or HPLC. It's a proposal for a mentor or a summer lab, not a home experiment.
  Carbonyl sulfide is a toxic gas and needs fume-hood handling (arm A).
- The 0.2% noise figure is an assumption. If the lab's method is noisier, rates scale up:
  the smallest detectable k2 grows roughly in proportion to the noise.
- A null result limits only these specific chemistries, not all possible ones.
- The rate model is the single-exponential early phase from `sim.py`. Late saturation isn't used for the fit.

## References
- Leman, Orgel & Ghadiri 2004, *Science*, doi:10.1126/science.1102722
- Forsythe et al. 2015, *Angew. Chem.*, doi:10.1002/anie.201503792
- Foden et al. 2020, *Science*, doi:10.1126/science.abd5680
- Ozturk et al. 2023, *Sci. Adv.*, doi:10.1126/sciadv.adg8274
- Viedma 2005, *PRL* 94:065504, doi:10.1103/PhysRevLett.94.065504
- Frank 1953, *Biochim. Biophys. Acta* 11:459
