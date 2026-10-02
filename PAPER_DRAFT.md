# Could the weak force pick life's hand in an icy-moon ocean?

**Saptarshi Ghosh**  
Cupertino High School, Cupertino, California, USA

*Preprint, September 2026. Code and data: https://github.com/helloamonkey742-code/chirality-project*

## Abstract

Life uses left-handed amino acids and right-handed sugars, and no one knows why. One old idea (Kondepudi & Nelson 1985) is that the weak nuclear force makes one hand very slightly lower in energy (the parity-violating energy difference, PVED), and that a large, slow, self-amplifying chemical system could turn that tiny bias into a one-handed outcome. We ask when the bias beats random chance, using stochastic simulations checked against a formula, and apply the answer to the subsurface oceans of Enceladus and Europa and to early Earth's ocean.

A rule of thumb comes out (§3.0): about 10³³–10³⁴ molecules, a few billion moles, must take part in the decision *together*. That means kilometre-scale volumes of well-mixed water at micromolar to millimolar concentrations. Ponds and single rock pores are many orders of magnitude too small. Icy-moon oceans are big enough, so the answer turns on how well they mix. Their salinity sets how strongly they are stratified, and a stratified ocean mixes slowly up and down. The patchwork of hands then freezes into stacked layers that merge only as fast as vertical mixing allows (§3.3; pictures in Fig. 2). Varying every uncertain input at once, physics wins in 22% (Enceladus), 29% (Europa) and 75% (early Earth) of plausible cases (§3.9).

The racemic amino acids in samples returned from asteroid Bennu, plausibly a fragment of a wet parent body, are a strong test (§3.10). Every amplifier in our assumed rate range would also have run on Bennu's parent body. If an amplified excess lasted forever, Bennu would rule out nearly all of them and cut physics' odds to 2–15%. But amplification is a non-equilibrium process. In a closed rock whose heat runs out, the same reaction network (run with its making steps reversible) rises to a large excess and then falls back to 50/50. In an ocean that vents keep driving, the excess holds. With that rule, a racemic Bennu fits the model, and the odds stay at 22% / 29% / 75% (15–23% / 21–30% / 55–75% across the cases tested). The idea makes a testable prediction: amino acids from a driven ocean, such as Enceladus' plume, should carry a lasting excess, while amplified excesses of molecules that can racemize should fade in extinct closed parent bodies. Meteorites fit this: their surviving excesses are in amino acids that cannot racemize or that locked into crystals (§3.10). Bennu's racemic isovaline, which cannot racemize, is not explained by the rule and remains open. Connected compartments, such as the pores of a vent mound, help only in proportion to the volume they connect within the decision time (§3.11). Long runs show no late drift or flipping at realistic molecule numbers (§3.12). No measured prebiotic reaction yet amplifies handedness this way, though a recently proposed network built from two measured reactions could (Higgs & Blackmond 2025), so we pre-register a wet-lab search sized for meteoritic concentrations (≤1 mM) and published measurement precision (§5). All code is public and self-checking (§6). The stakes are simple: if the weak force chose life's hand, it did so in an ocean, not a pond, and an excess that lasts in Enceladus' plume would be the first sign that a law of physics, not luck, set the handedness of life.

## 1 Introduction

Almost all life uses left-handed ("L") amino acids in proteins and right-handed ("D") sugars in DNA/RNA, even though both mirror-image forms are chemically identical alone. One idea, from Kondepudi & Nelson (1985), is that the weak force makes left- and right-handed molecules very slightly different in energy (PVED). That bias is astronomically small, but a large enough, slow enough, self-amplifying chemical system (each hand helps make more of itself while suppressing the other) could in principle let it beat chance.

No one had applied this size-and-time test to icy-moon oceans before this project; a scan of all 435 papers citing Kondepudi & Nelson 1985 (the paper that proposed the idea) or Brandenburg & Multamäki 2004 (closest prior spatial work, no weak force) found none that do (README, Novelty check). The closest work combines the PVED with metal-catalysed self-amplifying synthesis of sugars and amino acids in an early-Earth ocean (Cowan & Furnstahl 2022; Cowan 2023), or tests PVED selection in noisy open systems (Hochberg et al. 2022) or under a different weak-force bias, supernova neutrinos (Jannat, Shakeri & Shahbazi 2026, preprint). None of these asks how many molecules must decide together, or applies the question to Enceladus or Europa. The question: is an icy-moon or early-Earth ocean large and mixed enough for the weak-force bias to reliably beat chance, or is the outcome decided by luck? This is a student computational project. Every number below is either produced and self-checked by a script in the project folder, or a cited literature number (flagged where used), with the source README "Part" noted in brackets.

## 2 Methods

### 2.1 The core model (`sim.py`)

Let α = (L − D)/(L + D) be the excess of one hand. Near the point where a self-amplifying reaction "switches on", any such reaction reduces to the same simple equation (the normal form; Fig. 1):

$$d\alpha = \left(\lambda(t)\,\alpha - \alpha^3 + g\right)dt + \sqrt{\varepsilon}\,dW, \qquad \lambda = \gamma t .$$

- λ is how far past the switch-on point the system is. Before it (λ < 0) the only stable state is 50/50; after it (λ > 0) there are two stable states, all-L and all-D.
- g = ΔE_pv / kT is the weak-force bias (ΔE_pv is the PVED per molecule; its sign for real amino acids is unsettled).
- ε ≈ 1/N is the noise from counting individual molecules, N being the number of molecules that react together.
- γ = 1/(kτ) is how quickly conditions sweep through the switch-on point, where k = k₂c is the self-copying rate (rate constant k₂ times concentration c) and τ is the time the sweep takes. Time is measured in units of 1/k.

Linearising around λ = 0 gives the chance that the weak-force-favoured hand wins:

$$P = \Phi(\Delta), \qquad \Delta = \sqrt{2}\,\pi^{1/4}\, g\, \gamma^{-1/4}\, \varepsilon^{-1/2} = \sqrt{2}\,\pi^{1/4}\, g\,(k\tau)^{1/4}\sqrt{N},$$

where Φ is the normal cumulative distribution. Δ = 2 means a 97.7% win. The formula is checked against 4,000 brute-force runs per point (agreement ±0.005, Part 1).

![Fig. 1. The landscape the excess α moves in. The weak force tilts it very slightly; molecular noise shakes it. Which valley α falls into is decided in a short window around the switch-on point.](sketch.png)

### 2.2 Oceans (`ocean.py`, `scenarios.py`, `uncertainty.py`)

For a body of water, N = c · V · f · 1000 · N_A, where V is the volume (m³), f the fraction that is well mixed on the decision time, and N_A Avogadro's number. Physics "wins" only if two conditions both hold: the size condition Δ ≥ 2, and the time condition k₂cτ ≥ 10 (the reaction actually finishes). `uncertainty.py` draws 200,000 input sets per world, log-uniform over each plausible range: bias g, rate k₂, concentration c, horizontal and vertical mixing, time, volume, and the real-chemistry factor F (§2.4).

### 2.3 Space and mixing (`domains.py`, `domains3d.py`, `coarsen3d.py`, `k3.py`, `shell3d.py`, `aniso3d.py`, `sphere_layers.py`)

With mixing, α becomes a field α(x, t):

$$\partial_t \alpha = \lambda\alpha - \alpha^3 + g + D_h\nabla_h^2\alpha + D_z\,\partial_z^2\alpha + \text{noise}.$$

D_h and D_z are the horizontal and vertical eddy diffusivities. Rescaling depth by √(D_h/D_z) turns this exactly into the equal-mixing equation in a box that is taller by that factor. So weak vertical mixing acts like a deep, narrow ocean. After the decision, walls between regions of opposite hand move at D × curvature (Allen & Cahn 1979). Closed patches shrink, so patch size grows as l = A√(Dt), and the majority hand takes over (Fig. 2). Flat walls between horizontal layers do not move.

**Numerical scheme.** All spatial runs use explicit Euler (Euler–Maruyama when noise is included), which is first order in time, together with second-order centred differences (the 7-point Laplacian) in space. Walls are periodic sideways; top and bottom are closed (no flux). Units are 1/k for time and √(D/k) for length, with grid spacing 1 and time step 0.05. That satisfies the 3D stability limit D·dt/dx² ≤ 1/6. The resolution and time-step checks are in §3.12.

### 2.4 Real chemistry (`frank.py`, `frank2.py`)

The normal form is checked against a full Frank-type network in a flow reactor. Here A is the achiral feedstock (for example an amino-acid precursor), and L and D are the two hands of the product:

| reaction | rate | role |
|---|---|---|
| feed → A | f(a₀ − a), a₀ ramped slowly upward | supply; the ramp sweeps the system through switch-on |
| A → L, A → D | k₀a each | slow uncatalysed production (unbiased) |
| A + L → 2L | k(1 + g)aL | self-copying; the weak force enters as ±g |
| A + D → 2D | k(1 − g)aD | |
| L + D → waste | k_i LD | mutual destruction ("cross-inhibition") |
| L, D → out | k_d L, k_d D | outflow |

The two ingredients that make it amplify are self-copying and mutual destruction. Every reaction is simulated with its own counting noise (the chemical Langevin equation). At the switch-on point the network maps onto the normal form: g_eff = 2gxka, ε = (2k₀a + 2x(ka + k_d))/Ω, and γ = k·da/dt, where x is the concentration of each hand and Ω the reactor size. The extra noise from non-copying reactions lowers Δ by a factor F ≈ 0.7 compared with §2.2. `frank2.py` adds reverse reactions as a second scheme.

### 2.5 Other scripts

- `inherit.py`: a pool that inherits a small excess (Part 4b).
- `dual.py`: coupled amino-acid and sugar systems (Part 5).
- `scorecard.py`: real reactions scored against seven requirements (Part 6).
- `design.py`: experiment sizing (§5).
- `review.py`: the checks added after expert feedback (§3.0, §3.10–3.12, §5).
- `coarsen_fig.py`, `sketch.py`: the figures.
- `shell_res.py`: the grid-resolution test.
- `quench.py`: the excess frozen in when a closed parent body loses its water (§3.10), and whether an inherited excess survives a closed network (§3.5).

Every script asserts its own checks.

## 3 Results

### 3.0 A common-sense rule: how many molecules must decide together? (`review.py` R1)

Setting Δ = 2 in the formula gives the number of molecules that must react as one well-mixed pool: $N \ge \left(2 / (\sqrt{2}\,\pi^{1/4} g\,(k\tau)^{1/4})\right)^2$. For g = 10⁻¹⁷ (the middle of the modern PVED range) and a reaction that only just finishes (kτ = 10), N ≥ 3.6×10³³ molecules, or 6×10⁹ mol. At 1 µM that is a well-mixed cube of water 18 km on a side; at 1 mM, 1.8 km; at 0.1 M, 0.4 km. A 1 m³ pond at 1 M holds only 6×10²⁶ molecules, about 7 orders of magnitude too few.

The PVED is a fixed energy, so the bias g = ΔE/kT is larger in the cold. The needed N scales as T²: 3.6×10³³ at 0 °C, 4.3×10³³ at 27 °C, 6.7×10³³ at 100 °C. Cold helps only by a factor of about 2. Time helps only weakly, as $(k\tau)^{1/4}$: a sweep 10,000× slower lowers N by 100×. The practical message is that concentration alone cannot rescue a small setting. Only ocean-scale water that is mixed within the decision time can pool enough molecules.

### 3.1 Is an icy-moon ocean big enough? (Part 1)

At a representative bias g = 1e-17 and rate k2 = 1e-3 /M/s, the minimum reactant concentration needed for physics to win is ~2 mM for a small lake, but only ~50 pM for Enceladus, ~4 pM for Earth, and ~1 pM for Europa. Size is not the limit for icy-moon oceans across the modern PVED range. Time becomes the limit when chemistry is slow: at k2 = 1e-6 /M/s the time condition sets the answer for every ocean (~3 nM for Enceladus); temperature jitter barely matters, shifting the win probability by <0.01. If only a fraction f of the ocean is well mixed, the effective bias scales as g·√f — 1% mixing on Enceladus behaves like a bias ten times smaller.

### 3.2 One ocean, or a patchwork? (Part 2)

Adding 1D spatial mixing, three laws were derived and checked to a few percent. A single ocean-wide outcome needs D ≥ 3×10⁻⁶ m²/s (Enceladus) or 1×10⁻⁵ m²/s (Europa). Molecular diffusion (~10⁻⁹ m²/s) is far too slow for either, but the mixing rate modeled by Zeng & Jansen (2021), 5×10⁻⁵ m²/s, clears both. With only diffusion, oceans fragment into km-scale patches and the favoured hand gets just 50–58% of Enceladus (up to 79% of Europa at 1 µM) — a mostly chance-decided, mixed-hand ocean; with modelled mixing, physics wins ocean-wide down to ~1 nM (Enceladus) or ~1 pM (Europa).

### 3.3 The 3D correction — patches heal themselves (Part 3)

In full 3D, the 1D patch-size law fails (31% spread). Curved 3D patch walls straighten themselves out, so small patches shrink even with no bias: size grows as `l = A·√(D·t)`, exponents 0.45/0.49 (theory 0.5), A = 5.48/5.60 at two mixing strengths (a 2026-09-23 correction fixed a rounded "5.5 at both"; conclusion unchanged). Starting from a 55% majority, it grows to 66% (t=50), 81% (t=200), 93% (t=400) — standard phase-ordering physics. Healing within the moon's age needs D ≥ 6×10⁻⁶ m²/s (Enceladus) or 2×10⁻⁵ m²/s (Europa); diffusion alone would take 6×10¹¹–2×10¹³ years, modelled mixing heals in ~1×10⁷ (Enceladus) or ~4×10⁸ years (Europa). *Correction (2026-09-25):* the thresholds above apply one mixing value over the ocean's width, and the 10⁻¹⁰–10⁻³ m²/s range (and the 5×10⁻⁵ "modelled" value) are Zeng & Jansen's **vertical** diffusivity (their Sec. II.2). Treating the directions separately: sideways healing needs only ≥ ~6×10⁻⁴ (Enceladus) or ~2×10⁻⁴ m²/s (Europa), about 180–450× below the ~0.1 m²/s sideways eddy mixing estimated by Zhang, Kang & Marshall (2024) (50–140× below 0.03, the low end of their runs), so it never limits. Healing over the depth (≈40 km Enceladus, ≈120 km Europa) needs vertical mixing ≥ ~2×10⁻⁹–2×10⁻⁶ m²/s (Enceladus) or ~4×10⁻⁹–2×10⁻⁷ m²/s (Europa), which still sits inside the vertical range: *Update (same day, `aniso3d.py`):* that top-to-bottom rule is the wrong picture for most of the range. Weak vertical mixing is exactly equivalent to equal mixing in a taller ocean (stretch depth by √(D_h/D_z)); the stretched ocean is deeper than wide for 92–96% of the plausible range. In such tall boxes the patchwork froze into 3–5 stacked layers in 7/8 runs (unchanged from t = 500 to 1000), while a cube healed to one hand in 4/4 runs from a 55% start. On a round moon, layer boundaries are spheres: a wavy boundary was measured to flatten at D_h·k² regardless of D_z (0.00967 vs 0.00964 and 0.00239 vs 0.00241), which holds in a flat box. *Correction (2026-09-26):* on a sphere "down" follows the radius, so a constant-depth boundary feels only vertical mixing and sinks at 2·D_z/R, not 2D_h/R as first written (which gave takeover within ≤1.5×10⁴ yr Enceladus, ≤2.8×10⁵ yr Europa). `sphere_layers.py` measured sink speeds of 0.991/0.991/0.984 × 2·D_z/r with 0.0% change when D_h ×4; takeover now takes up to 1.5×10¹² yr (Enceladus) or 2.8×10¹³ yr (Europa) at D_z = 10⁻¹⁰ m²/s, and the merge-in-time condition holds in 33%/34% of draws (MEASURED). The hemisphere split is still unstable (rate 1.017× predicted; an off-equator band vanished), because of a shift off the great circle (a tilt is just a rotation, so cannot matter). The top layer decides with its own molecules (a share L/H' = 0.25 of the ocean, ASSUMED; a column-like sphere run measured 0.52 and 0.07 in 2 runs, neither supporting nor ruling it out); the favoured hand ended on top in 3/4 runs from 55% vs 2/4 from 50% (MEASURED, small sample). A direct per-patch check (`k3.py`) found an effective decision volume of 60.2 cells, consistent to 1% across three bias strengths (a correction fixed a rounded "≈57 cells"). Three more seeds gave K3 = 3.698, 3.481, 3.652 (original 3.642), a 2.3% spread across all four (MEASURED). Fig. 2 shows both behaviours: a cube healing to the majority, and a tall (weakly mixed) box freezing into layers that do not change from t = 400 to t = 5000.

![Fig. 2. Coarsening. Top: slices of a 96³ cube starting from a 55% majority; patches round off and shrink, and the majority reaches 99% by t = 400. Bottom: side view of a tall box (equivalent to weak vertical mixing); the patchwork freezes into flat layers of opposite hands, unchanged to t = 5000. Right: patch size grows as l = 5.5√(Dt) (exponent 0.45). `coarsen_fig.py`.](coarsening.png)

For early Earth (1.33×10¹⁸ m³, mean depth 3.7 km), a one-ocean outcome needs D ≥ 5×10⁻⁴–2×10⁻² m²/s — cleared by horizontal eddies (~10³ m²/s, Abernathey & Marshall 2013) and even by vertical mixing alone (1.1×10⁻⁵ m²/s, Ledwell et al. 1993: 94% win at 1 nM, 100% at 1 µM). Earth's own left-handedness can't confirm any of this, though — all life shares one ancestor, so any mechanism (chance included) ends with one hand; a second, independent ocean (Enceladus, Europa) is the real test.

### 3.4 Where was the hand decided? (Part 4)

Open ocean gives physics ~100% odds; lake/lagoon 51–100%; warm little pond 50–52% (coin flip); vent pore 50% (coin flip). The tension: settings chemists favour for concentrating molecules (ponds, vent pores) are exactly where physics loses, since they hold too few molecules total. The concentration inputs are the weakest part of this: Stribling & Miller (1987)'s ~3×10⁻⁴ M early-ocean amino-acid estimate is the only sourced number, with no lower estimate found in the literature — treat it as contested, not consensus.

### 3.5 Can small pools inherit a hand? (Part 4b)

A pool that fills with water already carrying a small excess can reliably follow it above a threshold, checked against 4,000 runs/case across 5 cases. Required inherited excess: 3×10⁻¹⁴ (lake/lagoon), 6×10⁻¹² (pond), 1×10⁻⁶ (vent pore). The weak force alone supplies only ~5×10⁻¹⁸ at equilibrium — far too small — but an ocean that already amplified the bias can supply up to ~1, and meteorites (Murchison L-isovaline, Glavin & Dworkin 2009, 18.5% excess) supply enough by 5–13 orders of magnitude. This made meteorite seeding look like the simplest explanation, though isovaline is not itself a protein amino acid and its excess would need catalytic pass-through (shown for isovaline → D-sugar, Pizzarello & Weber 2004). *Update:* samples returned from Bennu contain chiral non-protein amino acids that are racemic or nearly so (Glavin et al. 2025), and Ryugu's are also racemic (Furusho et al. 2024). The Bennu authors conclude that delivered molecules may not explain life's handedness. So meteoritic excesses are not universal; they vary with parent-body history (Elsila et al. 2016; Glavin et al. 2020a). Not every meteoritic excess was made on the parent body, either: Cooper & Rios (2016) report D excesses in sugar acids that may have formed very early and survived, so some excesses are *inherited*, not amplified in place. The model adds one constraint here (MODEL INFERENCE, `quench.py` Q4): a closed reversible network erases an inherited 20% excess just as it erases a 0.1% one, only 1.2–1.6× later. An inherited excess therefore survives only where no reversible amplifying network acts on it, or where the water left before the network relaxed. One caveat is that the excess a pond needs (6×10⁻¹²) is some 10 orders of magnitude below what any lab can measure (about ±1%). "Racemic within error" therefore neither supplies nor rules out a pond-sized seed.

### 3.6 Why L-amino acids pair with D-sugars (Part 5)

Two coupled chiral systems, each helping the partner's hand, were checked against 4,000 runs/case. Real chemistry links them both ways: L-amino acids catalyse D-sugar excess (Pizzarello & Weber 2004; Breslow & Cheng 2010); D-RNA prefers L-amino acids ~4× in one aminoacylation system (Tamura & Schimmel 2004), and in loop-closing ligation of aminoacyl-RNA, D-RNA favours L-amino acids while mirror-image L-RNA favours D (Kim et al. 2025). The link is not universal: self-aminoacylating D-ribozymes can favour either hand (Kenchel et al. 2024). Without coupling the hands match only 50% of the time; with κ ≈ √γ, 99.9%. The winning pair depends on the *sum* of the two biases — predicted/simulated probabilities of L-amino-acid + D-sugar ranged from 0.500 (biases cancel) to 0.800 (biases add) across four tested combinations. Coupling means physics only makes one combined decision — but the sugar PVED sign is unsettled, and if it actually favours L-sugars (opposite natural D), the biases would partly cancel, pushing back toward a coin flip. A 2005 whole-DNA-helix calculation (Faglioni et al.) concluded the weak force does *not* favour nature's helices, per its abstract. "Not the natural helix" is weaker than a strong push toward the mirror helix; the full paper was inaccessible (paywalled), so this is checked against the abstract only.

### 3.7 The missing reaction (Part 6)

Seven requirements for a self-amplifying prebiotic chiral reaction were derived (self-copying; mutual suppression; driven away from 50/50; works in water from simple molecules; fast enough; works from a tiny excess; biologically relevant product). Scoring seven real candidates, the best (Noorduin 2008 attrition deracemization) scores 6.0/7, but **no known system meets self-copying + suppression + prebiotic plausibility + relevant product together.** Speed is never the limiting factor — natural chemistry has 10⁵–10⁸ years, so even extremely slow reactions qualify. This gap motivates the wet-lab experiment in Section 5.

### 3.8 Does real chemistry follow the formula? (Part 7)

The network (§2.4, table) is the classic Frank (1953) scheme: an achiral feedstock A turns into L or D, each hand copies itself from A, and the two hands destroy each other on contact. Self-copying makes whichever hand leads grow faster, and mutual destruction removes the minority. Together they turn a small lead into a one-handed outcome. The network, simulated with molecular counting noise from every reaction step, follows the simplified formula within ~1 percentage point in the slow (natural) limit (2.5 points near saturation); a fast sweep selects *more* strongly than predicted, so the formula is conservative. The real network's selection strength is lower by a conversion factor F ≈ 0.7, because every reaction adding/removing a chiral molecule adds noise, not just self-copying ones — so every concentration in Parts 1–4 should rise by about 1.6×. No conclusion changes, since margins there were 10×–10⁶×, except the Europa mixing margin (not concentration-based).

### 3.9 Robustness checks (Part 8)

Varying every uncertain input at once (200,000 draws/world), physics wins in **75%** of plausible cases for early Earth's open ocean, **22%** for Enceladus and **29%** for Europa on a round moon (3% and 7% if layers never merged); the merge-in-time condition holds in 33%/34% of draws. Vertical mixing is the most decisive input (0%/0%/52% win rate in the low half of its range vs. 43%/57%/97% in the high half, Enceladus/Europa/Earth); concentration (+21/+10) and time (+15/+11; +32 for Earth) come next, then bias g (+8/+4) and rate k2 (+5); faster sideways mixing lowers the rate slightly (−2 to −9); F and volume ≤2 for the moons. *(Corrections 2026-09-25: first reported as 18% for both moons, from a run that used the vertical-mixing range for sideways healing; then 43%/58% before layering was modelled. Correction 2026-09-26: then 100%/55%/78% with concentration deciding (21%/57% vs. 89%/99%), bias g +18 to +23, k2 +14, vertical mixing +10 to +14, from a run that let layer boundaries sink at 2·D_h/R instead of 2·D_z/R.)* A second reaction scheme with back-reactions matches the formula within 1.7 points and gives F = 0.71 (vs. 0.70). A thin, wide "shell" shape closer to a real ocean's proportions still shows the majority taking over, just more slowly (growth exponent 0.41 vs. 0.45–0.49 in a cube).

### 3.10 Reconciling the model with racemic Bennu and Ryugu (`review.py` R2)

Bennu is plausibly a fragment of a wet parent body. Its chiral non-protein amino acids, including isovaline, which cannot racemize after it forms, are racemic or nearly so (Glavin et al. 2025). This is a direct test of the model, because the model assumes that some amplifying chemistry exists. The inputs below are ASSUMED ranges, chosen to bracket published values:

- one chiral amino acid (D + L) at 10–300 nmol/g of rock. The range spans Bennu's total amino acid content of 70 nmol/g (Glavin et al. 2025) and the 244 nmol/g of isovaline alone in the racemic CR chondrite EET 92042 (Glavin & Dworkin 2009); CR chondrites hold 17–3,300 nmol/g of amino acids in total (Aponte et al. 2020);
- a water-to-rock mass ratio of 0.1–1, giving 10⁻⁵ to 3×10⁻³ M in the fluid. Bennu's evaporites were modelled at 0.5 (0.5–1 allowed; McCoy et al. 2025), and the values cited for CI chondrites (0.5–1) and CM1 chondrites (0.2–0.7) fall inside this range (Zega et al. 2025);
- aqueous alteration lasting 1–10 Myr. Estimates of how long liquid water lasted on the CM parent body run as short as 10²–10⁴ yr (Glavin & Dworkin 2009); `quench.py` Q5 tests 10²–10⁷ yr.

Then c·τ on Bennu's parent body was 3×10⁸ to 9×10¹¹ M·s. Any amplifier with k₂ ≥ 10⁻¹¹ to 3×10⁻⁸ /M/s would have finished there. With only pore-water diffusion, the finished patches would be 1–3 km across, so each gram-sized sample would sit inside one patch and show a large excess of one hand or the other. The racemic samples say no such amplifier finished there.

Suppose the same chemistry, with the same rate constant, operated on Bennu's parent body and in the icy-moon oceans:

- **Rates in the original range (10⁻⁶ to 1 /M/s):** every draw would have finished on Bennu. Taken at face value, Bennu rules out the whole range.
- **Rates extended down to 10⁻¹² /M/s:** 23% of draws keep Bennu racemic. Among those, physics wins in only **2% (Enceladus), 4% (Europa) and 15% (early Earth)** of cases. An amplifier slow enough to stay idle on Bennu is usually too slow to finish in an icy moon as well. Only 1–5% of icy-moon draws have a larger c·τ than Bennu's maximum.

Taken at face value this is the strongest constraint in the paper. But it assumes that an amplified excess, once made, lasts. That is only true while something keeps driving the system.

**Closed rock vs open ocean (`review.py` R6).** Amplification is a non-equilibrium process. At equilibrium, both hands are equal (apart from the 10⁻¹⁷ PVED), so a closed system must end up racemic. This is known: in a closed, reversible Frank network the excess can rise for a while, a "chiral excursion", but the final state is racemic (Blanco, Stich & Hochberg 2011; Stich et al. 2016), and Higgs & Blackmond (2025) find the same for a prebiotic network (§4). What is new here is the comparison of Bennu's parent body with driven icy-moon oceans. To reproduce the result with this paper's own network, we took the `frank2.py` network, made the steps that make each hand (plain and self-copying, from feedstock A) reversible with rates that obey detailed balance, and ran it closed (no feed, no outflow). Detailed balance requires both routes from A to L to share one equilibrium constant (K = 2/s, where s is the factor that slows the reverse steps); otherwise the closed network could cycle A → L → A forever with no energy source. Each D step has the same rates as its L step, so equilibrium favours neither hand. The step in which L and D destroy each other stays one-way, turning them into a waste product W. Starting from a 0.1% seed:

- with each reverse rate constant half its forward one (s = 1, so K = 2), the excess passes a tenth of its peak after about 8 reaction times (1/k₂c), peaks at 70%, and falls back below a tenth of its peak after about 7×10³;
- with reverse steps 100× slower (s = 0.01, so K = 200), the excess passes a tenth of its peak after about 0.6 reaction times, peaks at 100%, and falls back below a tenth of its peak after about 2×10⁶;
- the same network kept open (fed and flushed) holds 99.8% for as long as it is driven.

Bennu's parent body was closed, and its liquid water lasted only while short-lived radioactive heat did. Its excess therefore stays at zero if the excess never rose (k₂cτ below the rise time) or if it rose and then faded (k₂cτ above the fade time). The second branch removes most of the conflict:

| k₂ range | excess lasts (reaction times) | P(win given a racemic Bennu): Enceladus / Europa / Earth |
|---|---|---|
| 10⁻⁶ – 1 | 8 to 7×10³ | 22% / 29% / 75% (unchanged) |
| 10⁻⁶ – 1 | 0.6 to 2×10⁶ | 23% / 30% / 75% |
| 10⁻¹² – 1 | 8 to 7×10³ | 15% / 21% / 57% |
| 10⁻¹² – 1 | 0.6 to 2×10⁶ | 17% / 22% / 55% |

This resolution has two assumptions. First, an icy-moon ocean must be continuously driven; Enceladus' measured plume H₂ points to ongoing water–rock reactions (Waite et al. 2017), and the plume also carries HCN, an amino-acid feedstock, alongside chemical disequilibrium that could power a reaction network (Peter et al. 2024), and phosphate from the ocean (Postberg et al. 2023). H₂ is a reactive fuel, so seeing it today means it is still being made; whether that supply has been steady over long times is not known (Glein & Zolotov 2020). The evidence leans toward a long-lived drive. Tidal heating of a porous core can keep hydrothermal activity going for tens of Myr to Gyr (Choblet et al. 2017). Serpentinization alone can supply the observed H₂ for at most about 500 Myr (Daval et al. 2022), but radiolysis can run for Gyr, and tidal grinding of rock could add bursts lasting up to Myr (Öberg, Magnabosco & Tosca 2026, preprint). Even 500 Myr is far longer than the decision times in this model. Second, the chemical drive on Bennu's parent body must have stopped while liquid water remained. Bennu's parent body held a brine that later evaporated or froze, and closed basins on Earth are the best analogue for it, since it would have received little or no fluid input while that happened (McCoy et al. 2025). So "closed" here cannot mean "had no water". It means the water stopped being fed. Bennu's minerals record alteration at about 25 °C by a fluid that probably evolved from neutral to alkaline. Whether that alteration was a closed system is debated: studies of CI chondrites argue for one, while carbonate veins in Bennu's boulders and fluid-mobile elements in its samples suggest an open one (Zega et al. 2025). The distinction that matters is how long the drive lasts compared with the water. On the parent body, both the heat (from short-lived ²⁶Al) and the fresh rock for water–rock reactions ran out within about 10 Myr, and Bennu's material is almost completely altered, which suggests the rock reactions ran to completion before the water was gone. Enceladus is still tidally heated and still producing H₂ today. If instead the drive and the water ended together, any excess would have frozen in, and for those cases the conflict returns. Bennu's veins and evaporite minerals (McCoy et al. 2025) show that its fluid moved and dried out. A ~2 mm Bennu grain preserves a sulfur-rich solvent front that stopped partway through the rock when the fluid was lost (Connolly et al. 2025, Figs 3–4), and without a solvent, reactions would have been slow and local (ASSUMED). That points to the drive and the water ending together, with any excess frozen in. The next paragraph tests that case directly.

**Freezing before equilibrium (`quench.py`).** Instead of the step rule above (racemic if k₂cτ falls before the rise time or after the fade time), read the excess off the full closed-network curve at the moment the water left, k₂cτ. "Racemic" here means |ee| < 5.2%, twice the 2.6% standard error of the Murchison isovaline value (20 measurements; Glavin & Dworkin 2009). Bennu-like draws (inputs as above) freeze as follows (MEASURED, `quench.log`):

| reverse steps | k₂ range | frozen racemic | 5–20% (the meteoritic isovaline band) | > 20% |
|---|---|---|---|---|
| half of forward (s = 1) | 10⁻⁶ – 1 | 97% | 1.2% | 1.9% |
| half of forward (s = 1) | 10⁻¹² – 1 | 73% | 9.0% | 18% |
| 1/200 of forward (s = 0.01) | 10⁻⁶ – 1 | 64% | 1.2% | 35% |
| 1/200 of forward (s = 0.01) | 10⁻¹² – 1 | 43% | 4.9% | 52% |

Three things follow. (1) *Given* a racemic Bennu, physics' odds match the step rule within 2 points in every case (22%/29%/75% and 23%/30%/75% for the original k₂ range, 16%/21%/57% and 17%/23%/57% for the extended range), so the table above stands. (2) But a racemic Bennu is the likely outcome mainly when reverse steps are fast (73–97% of draws). With slow reverse steps, 35–52% of Bennu-like draws freeze in an excess above 20%, so a racemic Bennu is about 1.5–1.7× less likely (43–64%). Measuring the reverse rates therefore decides how much Bennu counts against the model. (3) The closed network rarely freezes into a modest 5–20% excess (1.2–9% of draws): only 0.6–1.1 of the ~9 decades of k₂cτ on the simulated curve land there. The 5–20% range is used only as a reference scale taken from meteoritic isovaline (L-isovaline 0–20% across 17 meteorites; Dworkin et al. 2024). Isovaline cannot racemize and is not a reactant in this network, so this is not a test against the isovaline data; it only says that modest frozen excesses of a network product would need fine-tuned timing. Caveats: the seed is fixed at 0.1% and the network is the `frank2.py` scheme. The waste W is not counted as amino acid, and in a closed run it takes almost everything: while the excess is 5.2% or more, at most 0.3% (fast reverse steps) or 3.4% (slow) of the chiral material is still free, whereas 45 of Bennu's 70 nmol/g is free (Glavin et al. 2025). Once the excess has faded, less than 0.001% is still free. Bennu's free amino acids therefore cannot be what is left of this network; they would have to have formed outside the amplifier, which makes them a less direct test of it. If W were a bound form that released its L and D on hydrolysis, the hydrolysed total would never show more than 0.2% excess. `quench.py` Q5 tests how much the table rests on its other choices (MEASURED, `quench.log`). For the extended k₂ range the racemic share moves by about 2 points at most when the time unit changes (reading c as the whole 200-unit feedstock budget, or as the free pool at peak excess) or when Bennu's own concentration (7×10⁻⁵–7×10⁻⁴ M) and water lasting 10²–10⁷ yr replace the assumed inputs. For the original range it is not robust: the same changes move it anywhere from 26% to 100%. Letting the waste slowly release L and D (W → L + D at 10⁻⁶ per reaction time) only raises the racemic share (to 52–100%), and judging "racemic" by twice Bennu's own isovaline error (12.4%) raises it by less than 1 point for the original k₂ range and by 3–6 points for the extended range. The main open question in this section is now the reverse rates, not the timing. The resolution also makes a prediction that can be checked: samples from a driven ocean, such as Enceladus' plume, should carry a lasting excess, while amplified excesses of molecules that can racemize should fade in extinct, closed parent bodies. A racemic Enceladus would count against the idea. Without an amplifier, models expect the Enceladus ocean to be racemic, because racemization outpaces production (Steel, Dávila & McKay 2017) and even a biological excess would be erased within 10²–10⁴ years of ocean transport (Higgins et al. 2026, preprint). So any lasting excess in the plume would be a clean sign that something amplifies it. Capillary-electrophoresis instruments being built for ocean worlds can separate enantiomers in salty samples (Mora et al. 2022), so this is a practical measurement. Europa Clipper's dust analyser SUDA, a time-of-flight mass spectrometer that will look for organic molecules (Kempf et al. 2025), cannot make it, because mirror-image molecules have the same mass; the measurement needs a lander or returned samples with a chiral separation step. Radiation does not rule it out: amino acids survive radiolysis anywhere on Enceladus' surface, while on Europa samples would need to come from about 20 cm down (Pavlov et al. 2024).

**What the meteorites say about this rule.** Thermodynamics already requires it: a lasting excess needs a constant supply of free energy, so only a system held away from equilibrium, by a flow of matter or of energy, can keep one (Blackmond 2004; Plasson et al. 2007; Blanco et al. 2013; Ribó et al. 2017). A system closed to matter can keep one too if heat keeps driving it, for example through temperature cycling or through compartments at different temperatures linked by flow, as in hydrothermal vents (Ribó et al. 2013; Laurent et al. 2024). In this paper, "closed" means driven by neither. The meteorites fit if the rule is read as applying to excesses that stay dissolved and can racemize. In Murchison and Murray, the α-H amino acids, which racemize by losing their α-hydrogen, show no excess, and the excesses survive only in α-methyl amino acids such as isovaline, which cannot racemize that way (Pizzarello & Cronin 2000; Bada 1972; Cohen & Chyba 2000). Radiation from radioactive decay inside the parent body could still have racemized about 5% of the isovaline, which would only lower an excess (Glavin & Dworkin 2009). The large L-aspartic and L-glutamic excesses in Tagish Lake (43–59%) are in amino acids that crystallize as separate L and D crystals, while alanine, which cannot, is racemic. Glavin et al. (2012) explain them as amplification during alteration with the excess locked into solids, which takes it out of the reversible solution chemistry the rule describes. The pattern also matches a crystallization model: the amino acids with an excess in Tagish Lake, and the ones that are racemic, are the ones Blackmond's solubility-based model predicts (Klussmann et al. 2006 [confirm this is the intended model]). Only aspartic acid has good isotope data there. Tagish Lake is therefore better read as crystallization than as a driven network. Two points remain open. First, isovaline cannot racemize, so the fade rule does not explain Bennu's racemic isovaline. It says instead that Bennu's parent body either received no starting bias in isovaline or never amplified one, which is a statement about its history and not about this model. The Bennu team itself finds the racemic isovaline unexpected and offers no mechanism (Glavin et al. 2025; confirmed by Baczynski et al. 2026). Others take the racemic amino acids in pristine asteroid samples as support for a terrestrial origin of life's handedness (Ozturk & Sasselov 2025). Second, the Murchison isovaline excess grows with alteration (Glavin & Dworkin 2009), and the least-altered CM chondrite yet studied has racemic isovaline (Glavin et al. 2020b), although across 42 chondrite samples that link is weak (Aponte et al. 2025), so some amplification did happen in a closed body. The rule only says such an excess cannot last in solution, and solid or non-racemizable products escape that. If the rule is still set aside, the Bennu conflict returns and physics' odds fall back to 2%/4%/15%. Both results are reported.

If the closed-system argument fails, the other escape routes all say that the chemistry differed:

- the amplifier needed something Bennu's fluid lacked (longer peptides, mineral surfaces, warmer or longer-lived water);
- Bennu's water was mostly pore water in mud, not a free ocean;
- the amino-acid concentration in Bennu's fluid is not the concentration of the amplifier's reactants.

Each escape route is testable, but none is shown here. Temperature differences between the two settings are not modelled.

### 3.11 Compartments: can a network of connected pools beat one pond? (`review.py` R3)

Several authors suggest that life began in connected compartments, such as the pores of a hydrothermal mound (Martin & Russell 2003; Milner-White & Russell 2005; Russell 2006), acting together as one "reactor". We simulated 16 compartments, each with 1/16 of the molecules, exchanging material at rate q (normal form, g = 10⁻³, γ = 0.01):

| exchange rate q × decision time | P(majority favoured) | runs ending all one hand |
|---|---|---|
| 0 (isolated) | 0.675 | 0% |
| 0.1 | 0.686 | 0% |
| 1 | 0.724 | 53% |
| 10 | 0.714 | 100% |

A single pool holding all the molecules gives P = 0.724, and one isolated compartment gives 0.559. So compartments that exchange faster than the decision time (√(τ/k)) act as one pool of their combined volume, and afterwards they heal to one hand. Isolated compartments vote by majority but stay a mix of hands forever.

In a real vent mound (k₂ = 10⁻³ /M/s, c = 1 mM, 10⁵ yr), the decision takes about 56 years. Pores 1 mm apart exchange by diffusion 2×10⁶ times faster than that, so they are fully connected. Diffusion reaches only about 1.3 m in 56 years, however, so a diffusion-connected mound behaves like a ~2 m³ pond: P = 0.500, a coin flip. If fluid flow connects a 100 m mound (10⁶ m³), P = 0.73; a flow-connected cubic kilometre gives P = 1.000. Compartments do beat ponds, but only by the volume that fluid flow links together within the decision time. That makes a flowing vent system a serious candidate, consistent with §3.0.

### 3.12 Long runs and numerical resolution (`review.py` R4, `coarsen_fig.py`, `shell_res.py`)

*Does the decided hand drift later?* We continued the sweep with the reaction held at saturation for 100× longer than the decision itself. The share of runs on the favoured hand was 0.730 at the end of the sweep, 0.730 after 10³ more time units, and 0.730 after 10⁴. Once decided, a pool flips only by a rare noise-driven jump. We measured the flip rate at high noise and found it matches Kramers' rate (√2/2π)·e^(−1/2ε) within 5%: 1.8×10⁻², 7.7×10⁻³, 3.4×10⁻³ and 1.4×10⁻³ per unit time at ε = 0.2, 0.15, 0.12 and 0.1. With ε = 1/N the rate falls as e^(−N/2). A 1 µm³ droplet at 1 µM (about 600 molecules) can flip; a pond or ocean never does. The frozen layers in Fig. 2 show no change between t = 400 and t = 5000 (a 12× longer run).

*Grid resolution* (suggested by A. Brandenburg). The thin-shell runs (§3.9) use 8 grid points in depth, and a wall is only about 1.4 cells wide. We reran the same physical box (256 × 256 × 8 units) with grid spacing 1 (256 × 256 × 8 points) and 1/2 (512 × 512 × 16 points), with the time step scaled by dx². A run with grid spacing 1 and a 4× smaller time step gave identical patch sizes (30, 37, 47, 68, 91 at t = 25–400), so the time step is converged. On the fine grid, patch sizes were 31, 39, 48, 64 and 93, within 5% of the coarse grid at every time, and the growth exponent was 0.39 against 0.41. The coarse grid therefore resolves the coarsening, and the 8-point depth does not change the results; this is expected because the scheme is second order in space. We also ran a smaller box (128 × 128 × 8 units) with 2 seeds per case to check majority takeover. The coarse and fine grids agree in direction (a 55% majority grows on both), but runs from exact 50/50 vary too much between seeds to say more (`shell_res.py`).

## 4 Limitations

- All results are order-of-magnitude estimates; the F ≈ 0.7 real-chemistry conversion (Part 7) is checked for only two reaction schemes; other schemes could give a smaller F.
- The second scheme (`frank2.py`) gives its plain and self-copying making steps different equilibrium constants (0.005 and 2), which violates detailed balance. In its open, fed reactor this acts as a hidden extra drive (a futile A → L → A cycle, part of the wasted turnover that scheme measures). The closed-system test (§3.10) therefore uses one constant for both routes. Part 7 was not rerun with consistent constants.
- Part 2's spatial laws are 1D; Part 3's 3D correction uses finite 128³ boxes, and the 55%→93% late-time takeover is shown only qualitatively (patches reach the box edge).
- Real ocean turbulence isn't simple diffusion; icy-moon vertical mixing (Zeng & Jansen 2021) is a model input, uncertain across 10⁻¹⁰–10⁻³ m²/s — the single biggest reason the icy-moon answer isn't settled (Parts 2–3). Sideways mixing (~0.1 m²/s, Zhang et al. 2024) never limits healing. Ocean circulation models agree that salinity sets how stratified these oceans are (Kang et al. 2022a; Zeng & Jansen 2024), that the energy for vertical mixing (tides, libration) is very uncertain (Jansen et al. 2023), and that in a stratified Enceladus ocean material takes at least hundreds of years to travel from the floor to the ice (Kang et al. 2022b; Ames et al. 2025). An unstratified ocean heated from below is the fast-mixing end (Bire et al. 2022), and even with convection driven from below, material takes tens of years or more to travel from the seafloor to the ice (Zhang et al. 2025). None of these gives a vertical diffusivity directly, so the full range is kept. One preprint bounds the vertical diffusivity from the shape of Enceladus' ice shell: below 10⁻³ m²/s, the top of the range used here, unless the ocean's salinity is around 10 psu (Kang & Zhang 2026, preprint). On a round moon layers merge only as fast as vertical mixing allows (checked on a spherical shell, `sphere_layers.py`), so vertical mixing is again the key unknown; the deciding top layer's thickness is an estimate.
- τ (how slowly conditions change) is set to the ocean's age, the most favourable case possible.
- No measured prebiotic reaction performs the needed self-amplifying chirality; the Soai reaction does it in the lab but isn't prebiotic. The closest prebiotic candidates are peptide-ligation reactions that amplify an excess (Deng, Yu & Blackmond 2024) and a model network combining them with amino acid synthesis that breaks symmetry on its own (Higgs & Blackmond 2025). That model also finds the excess lasts only in an open, fed system and returns to racemic in a closed one, the same result as `review.py` R6. Its rates are given as ratios, so k2 is still a swept assumption, not measured (Part 6). The nearest lab-studied step is a peptide-catalysed transamination that makes enantioenriched alanine (Yu et al. 2024). Its conversions in the reverse direction (their Table 2) imply rough second-order rates of 6×10⁻⁶–1×10⁻⁴ /M/s depending on the catalyst (our estimate; corrected 2026-10-01 from "roughly 10⁻⁴ /M/s", which came from that reverse reaction). This is a conversion rate, not the amplification rate k2, but it is inside the swept range 10⁻⁶–1 /M/s. Lab tests of the Soai reaction need biases far larger than the PVED to set the outcome (Hawbaker & Blackmond 2019), as expected for flask-sized volumes (§3.0). Amino acids also break down in hot water: alanine lasts about 80 years at 100 °C but billions of years at 0 °C (our half-lives from the Arrhenius fits of Truong et al. 2019, extrapolated below their 386–613 K data). The cold bulk ocean is therefore not limiting, but near vents a slow amplifier must work on molecules that are continuously made. Expected Enceladus amino-acid concentrations, ≥10 nM without life and up to ~90 µM with it (Steel, Dávila & McKay 2017), fall inside the model's range.
- PVED has never been measured directly, and whether it or chance chose life's hand "remains completely open" (Quack, Seyfang & Wichmann 2022; experiments are still being built: Darquié et al. 2010; Cournol et al. 2019; review: Quack & Wichmann 2026), and calculations show it depends strongly on phase and geometry (Berger & Quack 2000; Bast et al. 2011); its sign for amino acids in water is unsettled (the "conformation problem", MacDermott et al. 2009a; that work also finds the L-zwitterion of alanine PVED-stabilized in its lowest-energy conformation), though the L-forms of the Murchison α-methyl amino acids are calculated to be PVED-stabilized (MacDermott et al. 2009b), and the sugar PVED sign is also unsettled, with some evidence pointing the "wrong" way, which would partly cancel the amino-acid bias (Part 5): for natural sugars the favoured hand depends on ring pucker, D for the DNA pucker but L for the RNA pucker (Tranter et al. 1992), DNA double helices are not favoured (Faglioni et al. 2005), and meteoritic sugars are near-racemic (Leyva et al. 2026).
- **Molecules flipping hand.** A chiral molecule can in principle tunnel to its mirror form. For amino acids the barrier makes this far slower than the age of the universe, and the PVED is larger than the tunnelling splitting, so each molecule stays in one hand (Quack 2002). The flipping that does happen is chemical: racemization by loss and return of the α-hydrogen. It is random and goes both ways at the same rate apart from the PVED tilt, so it adds noise rather than a bias, and it is the reverse-reaction route behind the closed-system fade (§3.10). A simulation that adds random flipping at rate r to the model (`review.py` R7) confirms this: the chance the weak force wins only drifts toward 50% (72% → 67% as r goes from 0 to 0.45 of the amplifier rate), and if flipping outpaces amplification no hand is chosen at all.
- **A possible second weak-force bias (beta decay).** Beta decay is a weak-force process that emits electrons spinning mostly one way. Vester, Ulbricht & Krauch (1959) proposed that these destroy one hand faster than the other. In the lab, low-energy spin-polarized electrons break up one enantiomer of bromocamphor more often than the other, by a few parts in 10⁴ (Dreiling & Gay 2014). An icy moon's rocky core contains radioactive ⁴⁰K, U and Th, so its water receives such electrons. The model does not include this. It would add to g, with the same sign in every ocean, but its sign for amino acids and its size in an ocean (dose, shielding by water, fraction of molecules hit) are unknown, so it is listed as a possible extra bias and not used. The experiments are mixed: spin-polarized electrons break bonds faster in one enantiomer on a magnetic surface (Rosenberg et al. 2008), but tests with positrons and beta sources found little or no effect (Gidley et al. 1982; Bonner 1991). A simulation (`review.py` R8) shows what it would take: with an asymmetry of 3 × 10⁻⁴, the beta effect matches the PVED once about one molecule in 10¹³ is destroyed by polarized electrons per reaction time. For an ocean at the edge (PVED alone gives 84%), that raises the odds to 98% if beta pushes the same way, or cuts them to 50% if it pushes the other way.
- The Enceladus ocean volume (2.7×10¹⁶ m³) is derived indirectly from ice-shell/core geometry (Čadek et al. 2016), self-checked by `ocean.py` against that paper's implied range (2.45–2.93×10¹⁶ m³); the ocean's age is separately debated (1 Myr–1 Gyr).
- Physics wins for Enceladus/Europa in 22%/29% of the plausible input range, early Earth 75% (Part 8) — not settled, and it falls to 3%/7% if layers never merge; vertical mixing is the key unknown. *History: first reported as 18%, then 43%/58% (corrected 2026-09-25), then 55%/78% (and 100% for Earth) with concentration as the key unknown; corrected 2026-09-26 because layer boundaries on a sphere sink at 2·D_z/R, not 2·D_h/R (`sphere_layers.py`).*
- **Bennu (§3.10).** If an amplified excess lasted forever, Bennu's racemic samples would cut physics' win rate to 2%/4%/15%. The closed-rock result removes most of that conflict, but it depends on Bennu's parent body behaving as a closed system and on its reverse reaction rates, which are unmeasured. If the water left before the network relaxed (`quench.py`), physics' odds given a racemic Bennu do not drop, but a racemic Bennu is about 1.5–1.7× less likely with slow reverse steps than with fast ones, and for the original k₂ range how likely it is depends on how model time maps onto k₂cτ (`quench.py` Q5). So the reverse rates decide how much Bennu counts against the model. Meteorite excesses survive only in non-racemizable or crystallized amino acids, which fits the rule, but it does not cover Bennu's isovaline, which cannot racemize (§3.10). The closed network also leaves almost no free amino acids once its excess has faded, unlike Bennu (45 of 70 nmol/g free), so Bennu's measured amino acids would have to come from outside the amplifier. Both the 22%/29%/75% and the 2%/4%/15% results should be treated as possible. All Bennu inputs (concentration, duration, same k₂ in both settings) are assumed.
- Concentrations: in Glavin & Dworkin (2009), the meteorites with an isovaline excess hold little isovaline (Murchison about 20 nmol/g, Orgueil about 0.7 nmol/g), from about 5× to over 300× less than the racemic CR chondrites QUE 99177 and EET 92042 (95 and 244 nmol/g). That is ≤ ~0.2 mM in parent-body fluid at a water-to-rock ratio of 0.1–1. CR chondrites hold 17–3,300 nmol/g of amino acids in total (Aponte et al. 2020). The pond (10⁻² M) and vent-pore (up to 1 M) values in Part 4 are upper bounds. Lowering them only strengthens the conclusion that small settings are coin flips.
- The mixing model uses simple eddy diffusion. Salinity-driven stratification enters only through the assumed vertical-diffusivity range (10⁻¹⁰ to 10⁻³ m²/s); no ocean circulation is simulated.
- Meteorite seeding (Part 4b) predicts the same left-handed outcome throughout the Solar System as the weak-force hypothesis — finding left-handed life on an icy moon can't by itself distinguish the two; only life from another star system could.

## 5 Planned experiment

Because no known prebiotic reaction satisfies the Part 6 requirements, EXPERIMENT.md pre-registers a wet-lab search for one; the pre-registration was written before any data were collected. Hypothesis: at least one prebiotic, water-based peptide-forming system will show seed-following growth of enantiomeric excess (ee).

**Revised after expert feedback.** The first design used 0.1 M and assumed 0.2% ee precision. Both were unrealistic:

- Meteoritic isovaline reaches about 20 nmol/g where it carries an excess (Murchison) and 244 nmol/g where it is racemic (EET 92042; Glavin & Dworkin 2009): about 0.2–2 mM in parent-body fluid at a water-to-rock ratio of 0.1, so 1 mM is realistic.
- Published ee uncertainties are ±0.01–1.5% by GC-MS and 1.2–7.2% by LC-MS. Glavin & Dworkin (2009) give 2.6% as the standard error of 20 measurements of Murchison isovaline (their racemic standard read −2.3 ± 1.3%, n = 14). The per-vial errors used below (0.2%, 1%, 2.6%) are assumptions within the published range; a single measurement scatters more than a 20-measurement standard error.

The smallest detectable rate after one year (`review.py` R5) is k₂,min = ln(1 + 3√2 σ/e₀)/(cT), where σ is the ee noise per sample, e₀ the seed excess, c the concentration and T the run length:

| concentration | seed | σ = 0.2% | σ = 1% | σ = 2.6% | 2.6%, 9 vials |
|---|---|---|---|---|---|
| 1 mM | 5% | 5×10⁻⁶ | 2×10⁻⁵ | 4×10⁻⁵ | 2×10⁻⁵ |
| 1 mM | 20% | 1×10⁻⁶ | 6×10⁻⁶ | 1×10⁻⁵ | 5×10⁻⁶ |
| 0.1 M | 5% | 5×10⁻⁸ | 2×10⁻⁷ | 4×10⁻⁷ | 2×10⁻⁷ |

At realistic concentrations the pond floor is 3×10⁻⁶ /M/s (1 mM, 100 yr), and the ocean floor is 3×10⁻⁹ (1 µM, 10⁸ yr). No lab run can reach the ocean floor at natural concentration. The revised design therefore has two tiers:

1. **Realistic arm:** 1 mM, a 20% seed (similar to Murchison isovaline), 3 years, and 9 independent vials per condition at 2.6% per vial. This reaches 1.8×10⁻⁶ /M/s, below the pond floor. Repeat injections of one vial are not replicates (Dworkin et al. 2024, test 5); with only 3 vials the same arm reaches 2.9×10⁻⁶, just 8% below the floor, so the plan uses 9 vials. It asks directly whether amplification happens at meteoritic concentrations.
2. **Accelerated arm:** 0.01 M and 0.1 M, stated plainly as *not* natural conditions. It exists to find any amplifier at all and to measure how the ee growth rate scales with concentration. That scaling (the rate law) is what the model needs in order to extrapolate to 1 µM.

Five candidate chemistries are tested, chosen from the Part 6 scorecard:

- (A) amino acids + carbonyl sulfide (Leman, Orgel & Ghadiri 2004);
- (B) amino acids + hydroxy acids under wet-dry cycling (Forsythe et al. 2015);
- (C) cysteine peptides catalysing their own joining (Foden et al. 2020);
- (D) an RNA precursor on magnetite (Ozturk et al. 2023);
- (E) peptide-catalysed transamination coupled to peptide ligation, run fed rather than closed (Yu et al. 2024; Deng, Yu & Blackmond 2024; Higgs & Blackmond 2025).

Arms E and C are the priority. E is a prebiotic network with a published model prediction of symmetry breaking (Higgs & Blackmond 2025), and a closed vial is predicted to return to 50/50 (`review.py` R6), so it runs in a slowly fed vessel. C is the only other catalytic loop in prebiotic water.

A Viedma (2005) grinding positive control checks that the pipeline can detect real amplification.

Each arm gets +L, +D (mirror control) and racemic seedings in 9 independent vials each (3 in the accelerated arm), plus these negative controls:

- no activator;
- sterile-filtered;
- seed alone, to measure racemization.

Analysts are blinded. An arm counts as an amplifier only if:

- both seeded sets rise by more than the detection threshold;
- the L- and D-seeded sets change by equal and opposite amounts;
- the racemic vials stay at zero;
- the controls stay flat.

The D-seed mirror arm is the main safeguard against contamination by biological (L) amino acids. The realistic arm would be proposed to local groups with chiral GC-MS (NASA Ames: G. Cooper, G. Chaban; San José State: A. Rios) rather than as funded outside work.

**As of this draft, this experiment has not been run.**

## 6 Reproducibility

All computational results were produced by the Python scripts in the project folder, each printing and asserting its own self-checks (numpy, scipy and matplotlib are required). The code is public at https://github.com/helloamonkey742-code/chirality-project (MIT licence) so that every number can be reproduced:

```
./run_all.sh          # fast subset, ~3 minutes
./run_all.sh --full   # adds the slow 3D and real-chemistry runs
```

Individual scripts (sim.py, ocean.py, domains.py, patches.py, domains3d.py, coarsen3d.py, scenarios.py, dual.py, inherit.py, scorecard.py, design.py, k3.py, frank.py, frank2.py, shell3d.py, aniso3d.py, sphere_layers.py, uncertainty.py, review.py, coarsen_fig.py, shell_res.py, quench.py) can also be run one at a time; approximate runtimes are listed in README.md's "Run" section.

## 7 Conclusions

- **A clear threshold.** For the weak force to beat chance, about 10³³–10³⁴ molecules must decide together. Ponds and rock pores fall short by many orders of magnitude; icy-moon oceans and early Earth's ocean clear it. If the weak force chose life's hand, it chose in an ocean.
- **Real odds, one main unknown.** Physics wins in 75% of plausible early-Earth cases and 22–29% for Enceladus and Europa. Most of the remaining uncertainty is a single number, how fast these oceans mix vertically, which ocean models are already working to pin down (§4).
- **The meteorite puzzle has an answer.** An amplified excess fades in a closed rock and holds in a driven ocean (§3.10). That fits racemic Bennu and the pattern of which meteorite amino acids keep their excess, with one open point (Bennu's isovaline).
- **It can be tested, soon and cheaply.** A bench-top search for the missing amplifier is sized for realistic concentrations and one year of runs (§5), and instruments built for ocean worlds can already measure handedness in salty samples (Mora et al. 2022). A lasting L-excess in Enceladus plume amino acids would support the idea; a racemic plume would count against it.
- **Why it matters.** A yes would tie the handedness of life to a fundamental force, and would predict that life anywhere in the Universe uses the same hand as ours.

## References

- Kondepudi & Nelson 1985, *Nature* 314:438, doi:10.1038/314438a0
- Quack 2002, "How important is parity violation for molecular and biomolecular chirality?", *Angew. Chem. Int. Ed.* 41:4618, doi:10.1002/anie.200290005
- MacDermott, Fu, Hyde, Nakatsuka et al. 2009, "Electroweak parity-violating energy shifts of amino acids: the conformation problem", *OLEB* 39:407, doi:10.1007/s11084-009-9161-x (2009a)
- MacDermott et al. 2009, "Parity-violating energy shifts of Murchison L-amino acids are consistent with an electroweak origin of meteorite L-enantiomeric excesses", *OLEB* 39:459, doi:10.1007/s11084-009-9162-9 (2009b)
- Chandrasekhar 2008, "Molecular homochirality and the parity-violating energy difference. A critique with new proposals", *Chirality* 20:84, doi:10.1002/chir.20502
- Brandenburg 2019, "The limited roles of autocatalysis and enantiomeric cross-inhibition in achieving homochirality in dilute systems", *OLEB* 49:49, doi:10.1007/s11084-019-09579-4
- Hochberg, Buhse, Micheau & Ribó 2023, "Chiral selectivity vs. noise in spontaneous mirror symmetry breaking", *PCCP* 25:31583, doi:10.1039/d3cp03311b
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
- Glavin & Dworkin 2009, "Enrichment of the amino acid L-isovaline by aqueous alteration on CI and CM meteorite parent bodies", *PNAS* 106:5487–5492, doi:10.1073/pnas.0811618106
- Connolly et al. 2025, "An overview of the petrography and petrology of particles from aggregate sample from asteroid Bennu", *Meteorit. Planet. Sci.* 60:979–996, doi:10.1111/maps.14335
- Dworkin, Elsila, Glavin, Aponte, McLain, Simkus, Graham & Parker 2024, "Verification of chiral asymmetry in meteoritic organics", 55th Lunar and Planetary Science Conference, Abstract #2645
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
- Glavin et al. 2025, "Abundant ammonia and nitrogen-rich soluble organic matter in samples from asteroid (101955) Bennu", *Nature Astronomy* 9:199–210, doi:10.1038/s41550-024-02472-9
- Furusho et al. 2024, "Enantioselective three-dimensional HPLC determination of amino acids in the Hayabusa2 returned samples from the asteroid Ryugu", *J. Chromatogr. Open*, doi:10.1016/j.jcoa.2024.100134
- Glavin et al. 2020a, "The search for chiral asymmetry as a potential biosignature in our Solar System", *Chem. Rev.* 120:4660, doi:10.1021/acs.chemrev.9b00474
- Elsila et al. 2016, "Meteoritic amino acids: diversity in compositions reflects parent body histories", *ACS Cent. Sci.* 2:370, doi:10.1021/acscentsci.6b00074
- Glavin et al. 2012, "Unusual nonterrestrial L-proteinogenic amino acid excesses in the Tagish Lake meteorite", *Meteorit. Planet. Sci.* 47:1347, doi:10.1111/j.1945-5100.2012.01400.x
- Martin, W. & Russell, M. J. 2003, "On the origins of cells: a hypothesis for the evolutionary transitions from abiotic geochemistry to chemoautotrophic prokaryotes, and from prokaryotes to nucleated cells", *Phil. Trans. R. Soc. Lond. B* 358(1429), 59–85, doi:10.1098/rstb.2002.1183 (the hydrothermal-mound compartment hypothesis; metadata checked on Crossref 2026-10-01)
- Russell, M. J. 2006, "First life", *American Scientist* 94, 32, doi:10.1511/2006.57.32 (popular review of the hydrothermal-mound compartment idea; earlier drafts cited it as "Russell, *Scientific American*"; metadata checked on Crossref 2026-10-01)
- Milner-White, E. J. & Russell, M. J. 2005, "Sites for phosphates and iron-sulfur thiolates in the first membranes: 3 to 6 residue anion-binding motifs (nests)", *Orig. Life Evol. Biosph.* 35(1), 19–27, doi:10.1007/s11084-005-4582-7 (suggested by A. Brandenburg; metadata checked on Crossref 2026-10-01; its title shows it is about peptide anion-binding motifs in the first membranes, so Martin & Russell 2003 is cited alongside it for the compartment idea; abstract not retrieved)
- Kramers 1940, *Physica* 7:284 (escape rate over a barrier)
- Waite et al. 2017, "Cassini finds molecular hydrogen in the Enceladus plume: evidence for hydrothermal processes", *Science* 356:155, doi:10.1126/science.aai8703
- Glein & Zolotov 2020, "Hydrogen, hydrocarbons, and habitability across the Solar System", *Elements* 16:47–52, doi:10.2138/gselements.16.1.47
- Mora, Kok, Noell & Willis 2022, "Detection of biosignatures by capillary electrophoresis mass spectrometry in the presence of salts relevant to ocean worlds missions", *Astrobiology* 22:914, doi:10.1089/ast.2021.0091
- Blackmond 2004, "Asymmetric autocatalysis and its implications for the origin of homochirality", *PNAS* 101:5732, doi:10.1073/pnas.0308363101 (closed systems and equilibrium)
- Vester, Ulbricht & Krauch 1959, "Optische Aktivität und die Paritätsverletzung im β-Zerfall", *Naturwissenschaften* 46:68, doi:10.1007/BF00599091
- Dreiling & Gay 2014, "Chirally sensitive electron-induced molecular breakup and the Vester-Ulbricht hypothesis", *Phys. Rev. Lett.* 113:118103, doi:10.1103/PhysRevLett.113.118103
- Higgs & Blackmond 2025, "Autocatalytic symmetry breaking and chiral amplification in a feedback network combining amino acid synthesis and ligation", *PNAS* 122:e2423683122, doi:10.1073/pnas.2423683122
- Deng, Yu & Blackmond 2024, "Symmetry breaking and chiral amplification in prebiotic ligation reactions", *Nature* 626:1019, doi:10.1038/s41586-024-07059-y
- Choblet et al. 2017, "Powering prolonged hydrothermal activity inside Enceladus", *Nature Astronomy* 1:841, doi:10.1038/s41550-017-0289-8
- Daval, Choblet, Sotin & Guyot 2022, "Theoretical considerations on the characteristic timescales of hydrogen generation by serpentinization reactions on Enceladus", *JGR Planets* 127, doi:10.1029/2021JE006995
- Öberg, Magnabosco & Tosca 2026, "Tidal rock grinding as a source of H₂ on Enceladus", arXiv:2606.16860 (preprint)
- Plasson, Kondepudi, Bersini, Commeyras & Asakura 2007, "Emergence of homochirality in far-from-equilibrium systems: mechanisms and role in prebiotic chemistry", *Chirality* 19:589, doi:10.1002/chir.20440
- Pizzarello & Cronin 2000, "Non-racemic amino acids in the Murray and Murchison meteorites", *Geochim. Cosmochim. Acta* 64:329, doi:10.1016/S0016-7037(99)00280-X
- Peter et al. 2024, "Detection of HCN and diverse redox chemistry in the plume of Enceladus", *Nature Astronomy* 8:164–173, doi:10.1038/s41550-023-02160-0
- Postberg et al. 2023, "Detection of phosphates originating from Enceladus's ocean", *Nature* 618:489–493, doi:10.1038/s41586-023-05987-9
- Blanco, Ribó, Crusats, El-Hachemi, Moyano & Hochberg 2013, "Mirror symmetry breaking with limited enantioselective autocatalysis and temperature gradients: a stability survey", *Phys. Chem. Chem. Phys.* 15:1546–1556, doi:10.1039/c2cp43488a
- Ribó, Hochberg, Crusats, El-Hachemi & Moyano 2017, "Spontaneous mirror symmetry breaking and origin of biological homochirality", *J. R. Soc. Interface* 14:20170699, doi:10.1098/rsif.2017.0699
- Hawbaker & Blackmond 2019, "Energy threshold for chiral symmetry breaking in molecular self-replication", *Nature Chemistry* 11:957–962, doi:10.1038/s41557-019-0321-y
- Bada 1972, "Kinetics of racemization of amino acids as a function of pH", *J. Am. Chem. Soc.* 94:1371–1373, doi:10.1021/ja00759a064
- Cohen & Chyba 2000, "Racemization of meteoritic amino acids", *Icarus* 145:272–281, doi:10.1006/icar.1999.6328
- Glavin et al. 2020b, "Abundant extraterrestrial amino acids in the primitive CM carbonaceous chondrite Asuka 12236", *Meteoritics & Planetary Science* 55:1979–2006, doi:10.1111/maps.13560
- Kang, Mittal, Bire, Campin & Marshall 2022a, "How does salinity shape ocean circulation and ice geometry on Enceladus and other icy satellites?", *Science Advances* 8:eabm4665, doi:10.1126/sciadv.abm4665
- Kang, Marshall, Mittal & Bire 2022b, "Ocean dynamics and tracer transport over the south pole geysers of Enceladus", *MNRAS* 517:3485–3494, doi:10.1093/mnras/stac2882
- Zeng & Jansen 2024, "The effect of salinity on ocean circulation and ice-ocean interaction on Enceladus", *Planetary Science Journal* 5:13, doi:10.3847/psj/ad0cba
- Jansen, Kang, Kite & Zhang 2023, "Energetic constraints on ocean circulations of icy ocean worlds", *Planetary Science Journal* 4:117, doi:10.3847/psj/acda95
- Ames, Ferreira, Czaja & Masters 2025, "Ocean stratification impedes particulate transport to the plumes of Enceladus", *Communications Earth & Environment* 6:63, doi:10.1038/s43247-025-02036-3
- Bire, Kang, Ramadhan, Campin & Marshall 2022, "Exploring ocean circulation on icy moons heated from below", *JGR Planets* 127:e2021JE007025, doi:10.1029/2021je007025
- Truong, Monroe, Glein, Anbar & Lunine 2019, "Decomposition of amino acids in water with application to in-situ measurements of Enceladus, Europa and other hydrothermally active icy ocean worlds", *Icarus* 329:140–147, doi:10.1016/j.icarus.2019.04.009
- Darquié et al. 2010, "Progress toward the first observation of parity violation in chiral molecules by high-resolution laser spectroscopy", *Chirality* 22:870–884, doi:10.1002/chir.20911
- Cournol et al. 2019, "A new experiment to test parity symmetry in cold chiral molecules using vibrational spectroscopy", *Quantum Electronics* 49:288–292, doi:10.1070/QEL16880
- Quack & Wichmann 2026, "Molecular chirality: from structure to the quantum dynamics of tunnelling, parity violation, a molecular quantum switch and the possible astrophysical detection of homochirality as a signature of extraterrestrial life", *CHIMIA* 80:389–400, doi:10.2533/chimia.2026.389
- Berger & Quack 2000, "Electroweak quantum chemistry of alanine: parity violation in gas and condensed phases", *ChemPhysChem* 1:57–60, doi:10.1002/1439-7641(20000804)1:1<57::AID-CPHC57>3.0.CO;2-J
- Bast et al. 2011, "Analysis of parity violation in chiral molecules", *Phys. Chem. Chem. Phys.* 13:864–876, doi:10.1039/C0CP01483D
- Rosenberg, Abu Haija & Ryan 2008, "Chiral-selective chemistry induced by spin-polarized secondary electrons from a magnetic substrate", *Phys. Rev. Lett.* 101:178301, doi:10.1103/PhysRevLett.101.178301
- Gidley, Rich, Van House & Zitzewitz 1982, "β decay and the origins of biological chirality: experimental results", *Nature* 297:639–643, doi:10.1038/297639a0
- Bonner 1991, "The origin and amplification of biomolecular chirality", *Orig. Life Evol. Biosph.* 21:59, doi:10.1007/BF01809580
- Cooper & Rios 2016, "Enantiomer excesses of rare and common sugar derivatives in carbonaceous meteorites", *PNAS* 113(24):E3322–E3331, doi:10.1073/pnas.1603030113
- Cowan & Furnstahl 2022, "Origin of chirality in the molecules of life", *ACS Earth Space Chem.* 6:2575–2581, doi:10.1021/acsearthspacechem.2c00032
- Cowan 2023, "Influence of the weak nuclear force on metal-promoted autocatalytic Strecker synthesis of amino acids", *Life* 14:66, doi:10.3390/life14010066
- Hochberg, Buhse, Micheau et al. 2022, "Resilience of parity-violation-induced chiral selectivity to nonequilibrium temperature fluctuations in open systems", *Phys. Rev. Research* 4:033183, doi:10.1103/PhysRevResearch.4.033183
- Jannat, Shakeri & Shahbazi 2026, "Supernova neutrinos and the origin of biomolecular homochirality", arXiv:2607.12813 (preprint)
- Quack, Seyfang & Wichmann 2022, "Perspectives on parity violation in chiral molecules: theory, spectroscopic experiment and biomolecular homochirality", *Chem. Sci.* 13:10598–10643, doi:10.1039/d2sc01323a
- Tranter, MacDermott, Overill & Speers 1992, "Computational studies of the electroweak origin of biomolecular handedness in natural sugars", *Proc. R. Soc. Lond. A* 436:603–615, doi:10.1098/rspa.1992.0037
- Leyva et al. 2026, "Abiotic sugar enantiomers in the CI carbonaceous chondrite Orgueil", *Nature Communications* 17, doi:10.1038/s41467-026-68709-5
- Aponte, McLain, Saeedi et al. 2025, "Challenges and opportunities in using amino acids to decode carbonaceous chondrite and asteroid parent body processes", *Astrobiology* 25:437–449, doi:10.1089/ast.2025.0017
- Baczynski et al. 2026, "Multiple formation pathways for amino acids in the early Solar System based on carbon and nitrogen isotopes in asteroid Bennu samples", *PNAS* 123:e2517723123, doi:10.1073/pnas.2517723123
- McCoy et al. 2025, "An evaporite sequence from ancient brine recorded in Bennu samples", *Nature* 637:1072–1077, doi:10.1038/s41586-024-08495-6
- Steel, Dávila & McKay 2017, "Abiotic and biotic formation of amino acids in the Enceladus ocean", *Astrobiology* 17:862–875, doi:10.1089/ast.2017.1673
- Higgins, Chen et al. 2026, "A framework for evaluating biosignature potential against the abiotic baseline on ocean worlds", arXiv:2605.15337 (preprint)
- Yu, Darù, Deng et al. 2024, "Prebiotic access to enantioenriched amino acids via peptide-mediated transamination reactions", *PNAS* 121:e2315447121, doi:10.1073/pnas.2315447121
- Aponte, Elsila, Hein et al. 2020, "Analysis of amino acids, hydroxy acids, and amines in CR chondrites", *Meteorit. Planet. Sci.* 55:2422–2439, doi:10.1111/maps.13586
- Zega et al. 2025, "Mineralogical evidence for hydrothermal alteration of Bennu samples", *Nature Geoscience* 18:832–839, doi:10.1038/s41561-025-01741-0
- Blanco, Stich & Hochberg 2011, "Temporary mirror symmetry breaking and chiral excursions in open and closed systems", *Chem. Phys. Lett.* 505:140–147, doi:10.1016/j.cplett.2011.02.032
- Stich, Ribó, Blackmond & Hochberg 2016, "Necessary conditions for the emergence of homochirality via autocatalytic self-replication", *J. Chem. Phys.* 145:074111, doi:10.1063/1.4961021
- Ribó, Crusats, El-Hachemi, Moyano, Blanco & Hochberg 2013, "Spontaneous mirror symmetry breaking in the limited enantioselective autocatalysis model: abyssal hydrothermal vents as scenario for the emergence of chirality in prebiotic chemistry", *Astrobiology* 13:132–142, doi:10.1089/ast.2012.0904
- Laurent, Göppel, Lacoste & Gerland 2024, "Emergence of homochirality via template-directed ligation in an RNA reactor", *PRX Life* 2:013015, doi:10.1103/PRXLife.2.013015
- Ozturk & Sasselov 2025, "Life's homochirality: across a prebiotic network", *PNAS* 122:e2505126122, doi:10.1073/pnas.2505126122
- Kempf et al. 2025, "SUDA: a SUrface Dust Analyser for compositional mapping of the Galilean moon Europa", *Space Sci. Rev.* 221:10, doi:10.1007/s11214-025-01134-0
- Pavlov, McLain, Glavin, Elsila, Dworkin, House & Zhang 2024, "Radiolytic effects on biological and abiotic amino acids in shallow subsurface ices on Europa and Enceladus", *Astrobiology* 24:698–709, doi:10.1089/ast.2023.0120
- Kim, Todisco, Radakovic & Szostak 2025, "Stereoselectivity of aminoacyl-RNA loop-closing ligation", *J. Am. Chem. Soc.* 147:19539–19546, doi:10.1021/jacs.4c16905
- Zhang, Bire, Wang, Nath, Ramadhan, Kang & Marshall 2025, "Long transit time from the seafloor to the ice shell on Enceladus", *MNRAS* 541:859–871, doi:10.1093/mnras/staf1008
- Kang & Zhang 2026, "Subsurface ocean salinity and dissipation rate inferred from Enceladus ice shell morphology", arXiv:2603.22602 (preprint)
- Kenchel et al. 2024, "Prebiotic chiral transfer from self-aminoacylating ribozymes may favor either handedness", *Nature Communications* 15:7980, doi:10.1038/s41467-024-52362-x
