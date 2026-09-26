# Could the weak force pick life's hand in an icy-moon ocean?

**DRAFT — not submitted; for student + mentor review; every number traces to README.md or EXPERIMENT.md**

## Abstract

Life uses left-handed amino acids and right-handed sugars, and no one knows why. One hypothesis (Kondepudi & Nelson 1985) is that the weak nuclear force gives mirror-image ("chiral") molecules a tiny energy difference — the parity-violating energy difference, or PVED — and that a large, slow, self-amplifying chemical system could turn that bias into a full one-handed outcome. This project builds and checks the math for when that bias beats random chance, using stochastic simulations in 1D and 3D and a real reaction scheme, then applies the results to Enceladus, Europa, and early Earth's oceans. Headline results: size is not the limit for icy-moon oceans (Part 1); in 3D a patchwork of hands heals and the majority wins; icy-moon oceans, which mix fast sideways but slowly vertically, first freeze into stacked layers of opposite hands, but on a round moon the curved layer boundaries sink and the top layer takes over within ~10⁴–10⁵ years, its hand set by its own majority (Parts 2–3b); small ponds can inherit a decided hand from a larger body of water or from meteorites, even though the weak force alone cannot supply enough bias directly (Part 4b); amino acids and sugars are forced to pair up by mutual catalysis, but the sugar bias's sign is unsettled and could cancel the amino-acid bias (Part 5); no known prebiotic reaction has all the needed properties (Part 6), so a 12-month, mirror-controlled wet-lab experiment has been designed but not yet run (EXPERIMENT.md). A full uncertainty sweep found physics wins in 100% of plausible cases for early Earth's open ocean, 55% for Enceladus and 78% for Europa, with concentration the most decisive input, but only 3% and 7% if the layers never merged (Part 8; corrected 2026-09-25 from 18% for both, then from 43%/58% before layering was modelled; see §3). Every computational step is complete and self-checked; the wet-lab step is the only piece not yet done.

## 1 Introduction

Almost all life uses left-handed ("L") amino acids in proteins and right-handed ("D") sugars in DNA/RNA, even though both mirror-image forms are chemically identical alone. One idea, from Kondepudi & Nelson (1985), is that the weak force makes left- and right-handed molecules very slightly different in energy (PVED). That bias is astronomically small, but a large enough, slow enough, self-amplifying chemical system (each hand helps make more of itself while suppressing the other) could in principle let it beat chance.

No one had applied this size-and-time test to icy-moon oceans before this project; a scan of all 435 papers citing Kondepudi & Nelson 1985 (the paper that proposed the idea) or Brandenburg & Multamäki 2004 (closest prior spatial work, no weak force) found none that do (README, Novelty check). The question: is an icy-moon or early-Earth ocean large and mixed enough for the weak-force bias to reliably beat chance, or is the outcome decided by luck? This is a student computational project. Every number below is either produced and self-checked by a script in the project folder, or a cited literature number (flagged where used), with the source README "Part" noted in brackets.

## 2 Methods

- **`sim.py`**: simulates a chiral system crossing its symmetry-breaking point, derives `P = Φ(Δ)`, `Δ = √2·π^¼·g·(kτ)^¼·√N`, checked against 4,000 brute-force runs/point to ±0.005 (Part 1).
- **`ocean.py`**: applies the formula to real ocean volumes/ages; needs a size condition (Δ ≥ 2) and a time condition (reaction finishes) (Part 1).
- **`domains.py`, `patches.py`**: 1D spatial mixing between patches; derive/check three laws (patch size, wall speed, per-patch outcome) and apply to the moons (Part 2).
- **`domains3d.py`, `coarsen3d.py`**: repeat on a 128³ 3D grid; the 1D patch-size law fails, so the real 3D growth/takeover behavior is measured instead (Part 3).
- **`k3.py`**: measures the per-patch 3D selection constant directly, rerun across extra seeds for spread (Part 3).
- **`scenarios.py`**: applies the formula to four origin-of-life settings (Part 4).
- **`inherit.py`**: pool inheriting a small chiral excess; checks a threshold formula against 4,000 runs/case × 5 cases (Part 4b).
- **`dual.py`**: coupled amino-acid/sugar systems helping each other's hand, checked against 4,000 runs/case (Part 5).
- **`scorecard.py`**: scores real candidate prebiotic reactions against seven derived requirements (Part 6).
- **`frank.py`**: simulates a real Frank-type reaction network in a flow reactor and measures a correction factor F vs. the simplified formula (Part 7). **`frank2.py`** repeats this with a second reaction scheme as a robustness check (Part 8).
- **`uncertainty.py`**: 200,000 random input draws/world (bias, rate, concentration, mixing, time, volume, F) (Part 8).
- **`shell3d.py`**: 3D coarsening in a thin, wide slab (512×512×8), 4 seeds (Part 8).

Every script prints/asserts its own checks. Full pipeline: `./run_all.sh` (fast, ~3 min) or `./run_all.sh --full` (adds slow 3D/real-chemistry runs).

## 3 Results

### 3.1 Is an icy-moon ocean big enough? (Part 1)

At a representative bias g = 1e-17 and rate k2 = 1e-3 /M/s, the minimum reactant concentration needed for physics to win is ~2 mM for a small lake, but only ~50 pM for Enceladus, ~4 pM for Earth, and ~1 pM for Europa. Size is not the limit for icy-moon oceans across the modern PVED range. Time becomes the limit when chemistry is slow: at k2 = 1e-6 /M/s the time condition sets the answer for every ocean (~3 nM for Enceladus); temperature jitter barely matters, shifting the win probability by <0.01. If only a fraction f of the ocean is well mixed, the effective bias scales as g·√f — 1% mixing on Enceladus behaves like a bias ten times smaller.

### 3.2 One ocean, or a patchwork? (Part 2)

Adding 1D spatial mixing, three laws were derived and checked to a few percent. A single ocean-wide outcome needs D ≥ 3×10⁻⁶ m²/s (Enceladus) or 1×10⁻⁵ m²/s (Europa). Molecular diffusion (~10⁻⁹ m²/s) is far too slow for either, but the mixing rate modeled by Zeng & Jansen (2021), 5×10⁻⁵ m²/s, clears both. With only diffusion, oceans fragment into km-scale patches and the favoured hand gets just 50–58% of Enceladus (up to 79% of Europa at 1 µM) — a mostly chance-decided, mixed-hand ocean; with modelled mixing, physics wins ocean-wide down to ~1 nM (Enceladus) or ~1 pM (Europa).

### 3.3 The 3D correction — patches heal themselves (Part 3)

In full 3D, the 1D patch-size law fails (31% spread). Curved 3D patch walls straighten themselves out, so small patches shrink even with no bias: size grows as `l = A·√(D·t)`, exponents 0.45/0.49 (theory 0.5), A = 5.48/5.60 at two mixing strengths (a 2026-09-23 correction fixed a rounded "5.5 at both"; conclusion unchanged). Starting from a 55% majority, it grows to 66% (t=50), 81% (t=200), 93% (t=400) — standard phase-ordering physics. Healing within the moon's age needs D ≥ 6×10⁻⁶ m²/s (Enceladus) or 2×10⁻⁵ m²/s (Europa); diffusion alone would take 6×10¹¹–2×10¹³ years, modelled mixing heals in ~1×10⁷ (Enceladus) or ~4×10⁸ years (Europa). *Correction (2026-09-25):* the thresholds above apply one mixing value over the ocean's width, and the 10⁻¹⁰–10⁻³ m²/s range (and the 5×10⁻⁵ "modelled" value) are Zeng & Jansen's **vertical** diffusivity (their Sec. II.2). Treating the directions separately: sideways healing needs only ≥ ~6×10⁻⁴ (Enceladus) or ~2×10⁻⁴ m²/s (Europa), about 180–450× below the ~0.1 m²/s sideways eddy mixing estimated by Zhang, Kang & Marshall (2024) (50–140× below 0.03, the low end of their runs), so it never limits. Healing over the depth (≈40 km Enceladus, ≈120 km Europa) needs vertical mixing ≥ ~2×10⁻⁹–2×10⁻⁶ m²/s (Enceladus) or ~4×10⁻⁹–2×10⁻⁷ m²/s (Europa), which still sits inside the vertical range: *Update (same day, `aniso3d.py`):* that top-to-bottom rule is the wrong picture for most of the range. Weak vertical mixing is exactly equivalent to equal mixing in a taller ocean (stretch depth by √(D_h/D_z)); the stretched ocean is deeper than wide for 92–96% of the plausible range. In such tall boxes the patchwork froze into 3–5 stacked layers in 7/8 runs (unchanged from t = 500 to 1000), while a cube healed to one hand in 4/4 runs from a 55% start. On a round moon, layer boundaries are spheres: a wavy boundary was measured to flatten at D_h·k² regardless of D_z (0.00967 vs 0.00964 and 0.00239 vs 0.00241), so a spherical boundary sinks at 2D_h/R and the top layer takes over within ≤1.5×10⁴ yr (Enceladus) or ≤2.8×10⁵ yr (Europa). The top layer decides with its own molecules (a share L/H' of the ocean, ASSUMED); the favoured hand ended on top in 3/4 runs from 55% vs 2/4 from 50% (MEASURED, small sample). A direct per-patch check (`k3.py`) found an effective decision volume of 60.2 cells, consistent to 1% across three bias strengths (a correction fixed a rounded "≈57 cells"). Three more seeds gave K3 = 3.698, 3.481, 3.652 (original 3.642), a 2.3% spread across all four (MEASURED).

For early Earth (1.33×10¹⁸ m³, mean depth 3.7 km), a one-ocean outcome needs D ≥ 5×10⁻⁴–2×10⁻² m²/s — cleared by horizontal eddies (~10³ m²/s, Abernathey & Marshall 2013) and even by vertical mixing alone (1.1×10⁻⁵ m²/s, Ledwell et al. 1993: 94% win at 1 nM, 100% at 1 µM). Earth's own left-handedness can't confirm any of this, though — all life shares one ancestor, so any mechanism (chance included) ends with one hand; a second, independent ocean (Enceladus, Europa) is the real test.

### 3.4 Where was the hand decided? (Part 4)

Open ocean gives physics ~100% odds; lake/lagoon 51–100%; warm little pond 50–52% (coin flip); vent pore 50% (coin flip). The tension: settings chemists favour for concentrating molecules (ponds, vent pores) are exactly where physics loses, since they hold too few molecules total. The concentration inputs are the weakest part of this: Stribling & Miller (1987)'s ~3×10⁻⁴ M early-ocean amino-acid estimate is the only sourced number, with no lower estimate found in the literature — treat it as contested, not consensus.

### 3.5 Can small pools inherit a hand? (Part 4b)

A pool that fills with water already carrying a small excess can reliably follow it above a threshold, checked against 4,000 runs/case across 5 cases. Required inherited excess: 3×10⁻¹⁴ (lake/lagoon), 6×10⁻¹² (pond), 1×10⁻⁶ (vent pore). The weak force alone supplies only ~5×10⁻¹⁸ at equilibrium — far too small — but an ocean that already amplified the bias can supply up to ~1, and meteorites (Murchison L-isovaline, Glavin & Dworkin 2009, 18.5% excess) supply enough by 5–13 orders of magnitude. This makes meteorite seeding the simplest explanation, though isovaline is not itself a protein amino acid and its excess would need catalytic pass-through (shown for isovaline → D-sugar, Pizzarello & Weber 2004).

### 3.6 Why L-amino acids pair with D-sugars (Part 5)

Two coupled chiral systems, each helping the partner's hand, were checked against 4,000 runs/case. Real chemistry links them both ways: L-amino acids catalyse D-sugar excess (Pizzarello & Weber 2004; Breslow & Cheng 2010); D-RNA prefers L-amino acids ~4× (Tamura & Schimmel 2004). Without coupling the hands match only 50% of the time; with κ ≈ √γ, 99.9%. The winning pair depends on the *sum* of the two biases — predicted/simulated probabilities of L-amino-acid + D-sugar ranged from 0.500 (biases cancel) to 0.800 (biases add) across four tested combinations. Coupling means physics only makes one combined decision — but the sugar PVED sign is unsettled, and if it actually favours L-sugars (opposite natural D), the biases would partly cancel, pushing back toward a coin flip. A 2005 whole-DNA-helix calculation (Faglioni et al.) concluded the weak force does *not* favour nature's helices, per its abstract. "Not the natural helix" is weaker than a strong push toward the mirror helix; the full paper was inaccessible (paywalled), so this is checked against the abstract only.

### 3.7 The missing reaction (Part 6)

Seven requirements for a self-amplifying prebiotic chiral reaction were derived (self-copying; mutual suppression; driven away from 50/50; works in water from simple molecules; fast enough; works from a tiny excess; biologically relevant product). Scoring seven real candidates, the best (Noorduin 2008 attrition deracemization) scores 6.0/7, but **no known system meets self-copying + suppression + prebiotic plausibility + relevant product together.** Speed is never the limiting factor — natural chemistry has 10⁵–10⁸ years, so even extremely slow reactions qualify. This gap motivates the wet-lab experiment in Section 5.

### 3.8 Does real chemistry follow the formula? (Part 7)

A full Frank-type network, simulated with molecular counting noise from every reaction step, follows the simplified formula within ~1 percentage point in the slow (natural) limit (2.5 points near saturation); a fast sweep selects *more* strongly than predicted, so the formula is conservative. The real network's selection strength is lower by a conversion factor F ≈ 0.7, because every reaction adding/removing a chiral molecule adds noise, not just self-copying ones — so every concentration in Parts 1–4 should rise by about 1.6×. No conclusion changes, since margins there were 10×–10⁶×, except the Europa mixing margin (not concentration-based).

### 3.9 Robustness checks (Part 8)

Varying every uncertain input at once (200,000 draws/world), physics wins in **100%** of plausible cases for early Earth's open ocean, **55%** for Enceladus and **78%** for Europa on a round moon where layers merge (3% and 7% if they never merged). Concentration is the most decisive input (21%/57% win rate in the low half of its range vs. 89%/99% in the high half); bias g (+18 to +23 points), rate k2 (+14) and vertical mixing (+10 to +14) come next; faster sideways mixing lowers the rate slightly (−6 to −8, a thinner deciding layer); time 4–8 points; F and volume ≤2. *(Corrections 2026-09-25: first reported as 18% for both moons, from a run that used the vertical-mixing range for sideways healing; then 43%/58% before layering was modelled.)* A second reaction scheme with back-reactions matches the formula within 1.7 points and gives F = 0.71 (vs. 0.70). A thin, wide "shell" shape closer to a real ocean's proportions still shows the majority taking over, just more slowly (growth exponent 0.41 vs. 0.45–0.49 in a cube).

## 4 Limitations

- All results are order-of-magnitude estimates; the F ≈ 0.7 real-chemistry conversion (Part 7) is checked for only two reaction schemes; other schemes could give a smaller F.
- Part 2's spatial laws are 1D; Part 3's 3D correction uses finite 128³ boxes, and the 55%→93% late-time takeover is shown only qualitatively (patches reach the box edge).
- Real ocean turbulence isn't simple diffusion; icy-moon vertical mixing (Zeng & Jansen 2021) is a model input, uncertain across 10⁻¹⁰–10⁻³ m²/s — the single biggest reason the icy-moon answer isn't settled (Parts 2–3). Sideways mixing (~0.1 m²/s, Zhang et al. 2024) never limits healing. The icy-moon win rates rest on layers merging because the moon is round, which is checked in 2D pieces, not in a full 3D spherical-shell simulation; the deciding top layer's thickness is an estimate.
- τ (how slowly conditions change) is set to the ocean's age, the most favourable case possible.
- No known prebiotic reaction performs the needed self-amplifying chirality; the Soai reaction does it in the lab but isn't prebiotic; rate k2 is a swept assumption, not measured (Part 6).
- PVED has never been measured directly; its sign for amino acids in water is unsettled (the "conformation problem", 2009), and the sugar PVED sign is also unsettled, with some evidence pointing the "wrong" way, which would partly cancel the amino-acid bias (Part 5).
- The Enceladus ocean volume (2.7×10¹⁶ m³) is derived indirectly from ice-shell/core geometry (Čadek et al. 2016), self-checked by `ocean.py` against that paper's implied range (2.45–2.93×10¹⁶ m³); the ocean's age is separately debated (1 Myr–1 Gyr).
- Physics wins for Enceladus/Europa in 55%/78% of the plausible input range (Part 8; first reported as 18%, then 43%/58%, corrected 2026-09-25) — likely-ish, not settled, and it collapses to 3%/7% if layers never merge; concentration is the key unknown.
- Meteorite seeding (Part 4b) predicts the same left-handed outcome throughout the Solar System as the weak-force hypothesis — finding left-handed life on an icy moon can't by itself distinguish the two; only life from another star system could.

## 5 Planned experiment

Because no known prebiotic reaction satisfies the Part 6 requirements, EXPERIMENT.md pre-registers (written before any data collected) a 12-month wet-lab search for one. Hypothesis, stated in advance: at least one prebiotic, water-based peptide-forming system will show seed-following growth of enantiomeric excess (ee) with an effective rate k2 ≥ 10⁻⁷ /M/s — slow enough to be missed by typical short experiments, fast enough to matter for a pond (floor rate 3×10⁻⁷ /M/s).

The design is seeded, not starting from 50/50, because a purely self-copying system keeps ee constant — ee only grows if the hands also suppress each other, the exact combination being tested. Based on `design.py`'s detectability calculations, the chosen conditions are 0.1 M for 12 months (covers the pond-floor threshold with margin by month 2); a parallel 1 M arm for 1–3 years, where solubility allows, would reach the ocean-floor threshold.

Four candidate chemistries are tested, chosen from the Part 6 scorecard for prebiotic plausibility and a biologically relevant product: (A) amino acids + volcanic carbonyl sulfide forming peptides in water (Leman, Orgel & Ghadiri 2004); (B) amino acids + hydroxy acids under wet-dry cycling (Forsythe et al. 2015); (C) prebiotic cysteine peptides catalysing their own joining in neutral water (Foden et al. 2020) — the priority arm, the only one with a genuine catalytic loop in prebiotic water; and (D) an RNA precursor on magnetite followed by crystallization, non-autocatalytic (Ozturk et al. 2023). A positive control (Viedma 2005 grinding) checks the pipeline can detect real amplification at all.

Each arm gets triplicate +5% L, +5% D (mirror control — most important), and 0% (racemic) seedings, plus three negative controls (no activator/no cycling; sterile-filtered; seed amino acid alone, to measure plain racemization), sampled at eleven timepoints over a year by chiral GC-MS or HPLC with blinded analysts. A decision rule fixed in advance calls an arm an "amplifier" only if: ee rises >0.85% absolute in both seeded sets; L- and D-seeded vials change by equal and opposite amounts; racemic vials stay at zero; and sterile/no-activator controls show no change. The main risk flagged is contamination (biology and lab workers are left-handed); the D-seed mirror arm is the main safeguard, since contamination pushes both seeded sets toward L while real amplification pushes them apart symmetrically.

**As of this draft, this experiment has not been run.** It requires a chemistry lab with chiral GC-MS or HPLC access — a proposal for a mentor or summer lab, not a home experiment — and carbonyl sulfide (arm A) needs fume-hood handling.

## 6 Reproducibility

All computational results were produced by the Python scripts in the project folder, each printing/asserting its own self-checks (numpy, scipy, matplotlib required):

```
./run_all.sh          # fast subset, ~3 minutes
./run_all.sh --full   # adds the slow 3D and real-chemistry runs
```

Individual scripts (sim.py, ocean.py, domains.py, patches.py, domains3d.py, coarsen3d.py, scenarios.py, dual.py, inherit.py, scorecard.py, design.py, k3.py, frank.py, frank2.py, shell3d.py, uncertainty.py) can also be run one at a time; approximate runtimes are listed in README.md's "Run" section.

## References

(Only sources cited in README.md/EXPERIMENT.md, in the form given there; fields README omitted are marked "[details in README]".)

- Kondepudi & Nelson 1985, *Nature* 314:438, doi:10.1038/314438a0
- Quack 2002, "How important is parity violation for molecular and biomolecular chirality?" [details in README]
- "Electroweak parity-violating energy shifts of amino acids: the conformation problem", 2009 [details in README]
- "Molecular homochirality and the PVED: a critique with new proposals", 2007 [details in README]
- "The limited roles of autocatalysis and enantiomeric cross-inhibition in achieving homochirality in dilute systems", 2019, doi:10.1007/s11084-019-09579-4
- "Chiral selectivity vs. noise in spontaneous mirror symmetry breaking", 2023, doi:10.1039/d3cp03311b
- Čadek et al. 2016, *Geophysical Research Letters* 43, doi:10.1002/2016GL068634
- Postberg et al. 2018, *Nature*, doi:10.1038/s41586-018-0246-4
- Anderson et al. 1998, *Science* 281:2019
- Brandenburg & Multamäki 2004, *Int. J. Astrobiology*, doi:10.1017/s1473550404001983
- Sandars 2005, *Int. J. Astrobiol.*, doi:10.1017/s1473550405002338
- Gleiser & Walker 2008, "Punctuated chirality", *OLEB*, doi:10.1007/s11084-008-9147-0
- Jafarpour, Biancalani & Goldenfeld 2017, *PRE* 95:032407, doi:10.1103/PhysRevE.95.032407
- Zeng & Jansen 2021, *PSJ* 2:151, doi:10.3847/PSJ/ac1114
- Zhang, Kang & Marshall 2024, *Science Advances* 10, doi:10.1126/sciadv.adn6857
- Allen & Cahn 1979, *Acta Metallurgica* 27:1085, doi:10.1016/0001-6160(79)90196-2 (boundaries move at D × curvature)
- Charette & Smith 2010, *Oceanography* 23(2), doi:10.5670/oceanog.2010.09
- Ledwell, Watson & Law 1993, *Nature* 364:701, doi:10.1038/364701a0
- Abernathey & Marshall 2013, *JGR Oceans* 118:901
- Bray 1994, "Theory of phase-ordering kinetics", *Adv. Phys.* 43:357
- Stribling & Miller 1987, *OLEB* 17:261, doi:10.1007/BF02386466
- Pearce et al. 2017, *PNAS* 114:11327, doi:10.1073/pnas.1710339114
- Baaske et al. 2007, *PNAS* 104:9346, doi:10.1073/pnas.0609592104
- Toner & Catling 2020, *PNAS* 117:883, doi:10.1073/pnas.1916109117
- Pizzarello & Weber 2004, *Science* 303:1151, doi:10.1126/science.1093057
- Breslow & Cheng 2010, *PNAS* 107:5723, doi:10.1073/pnas.1001639107
- Hein, Tse & Blackmond 2011, *Nat. Chem.* 3:704, doi:10.1038/nchem.1108
- Tamura & Schimmel 2004, *Science* 305:1253, doi:10.1126/science.1099141
- Ozturk et al. 2023, *Sci. Adv.* 9:eadg8274, doi:10.1126/sciadv.adg8274
- Faglioni, D'Agostino, Cadioli & Lazzeretti 2005, "Parity violation energy of biomolecules II: DNA", *Chem. Phys. Lett.* 407:522, doi:10.1016/j.cplett.2005.04.009
- Glavin & Dworkin 2009, *PNAS* 106:5487, doi:10.1073/pnas.0811618106
- Frank 1953, *Biochim. Biophys. Acta* 11:459
- Soai et al. 1995, *Nature* 378:767, doi:10.1038/378767a0; Sato et al. 2003, *Angew. Chem.* 42:315, doi:10.1002/anie.200390105
- Viedma 2005, *PRL* 94:065504, doi:10.1103/PhysRevLett.94.065504
- Noorduin et al. 2008, *Angew. Chem.* 47:6445, doi:10.1002/anie.200801846
- Saghatelian et al. 2001, *Nature* 409:797, doi:10.1038/35057238
- Joyce et al. 1984, *Nature* 310:602, doi:10.1038/310602a0
- Klussmann et al. 2006, *Nature* 441:621, doi:10.1038/nature04780
- Leman, Orgel & Ghadiri 2004, *Science*, doi:10.1126/science.1102722
- Forsythe et al. 2015, *Angew. Chem.*, doi:10.1002/anie.201503792
- Foden et al. 2020, *Science*, doi:10.1126/science.abd5680
