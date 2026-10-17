# The volume experiment: does the Soai reaction count every molecule?

## Aim
Test whether the seed needed to steer the Soai reaction grows with reaction volume, as the counting
rule predicts for one well-mixed pool. The rule (Eqs. 10–11 of `journals/ija_draft.md`) says a seed of
m net molecules among N_d deciding molecules wins with probability

  P = Φ(m / √N_d),  Δ = m / √N_d.

Eq. 13 takes every molecule as an independent vote, so its N is a minimum (Section 2.2.4). Section 3.7
ends by noting that no one has yet tested whether a solution reaction reaches that best case. This experiment
tests it.

**m\*, not m50.** In this model P = 0.5 only at m = 0, so the steering threshold is defined as the seed
that gives Δ = 1 (P = 0.84): m\* = √N_d. The brief called this "m50".

## Hypotheses (stated before any data)
At a fixed concentration, every volume V has the same chemistry. Only the number of molecules that can
vote together changes.

| hypothesis | deciding pool | prediction | slope β of ln m\* on ln V | m\*(10 mL) / m\*(0.1 mL) |
|---|---|---|---|---|
| **H1** one well-mixed pool | N_d = c\*·V | m\* ∝ √V | **0.5** | **10** |
| **H0** local or early events | N_d set by a few early units, independent of V | m\* flat | **0** | **1** |

Two intermediate outcomes are also named in advance. A pool set at the walls (N_d ∝ area ∝ V^⅔) gives
β = ⅓. A pool limited by mixing gives 0 < β < 0.5 that changes with stirring. Arms S and M test these.

**Anchor.** Singleton & Vo (2003) found that about 60,000 seed molecules controlled the hand.
`soai_check.py` S4 gives N_d ≤ 3.6×10⁹ at Δ = 1, so m\* ≈ 6×10⁴ at their scale.
- ASSUMPTION: their scale is 1 mL. Their flask volume is not in the open text.
- Stage 1 measures m\* at 1 mL before any ladder is fixed, so the anchor only sets the starting bracket.
- At the anchor, c\* ≈ 6×10⁻¹² M: one molecule in about 10¹⁰ of the flask's aldehyde votes. A pool this small is plausibly
  the product present when the hand is set, not the whole flask.

**Absorbed constants.** The rule's prefactor a (P = Φ(a·ε·√N)) is taken as 1. Any a ≠ 1 that does not
depend on V is absorbed into m\* and leaves β unchanged. So does any loss in seed delivery shared by all
volumes, and any voting unit larger than one molecule (see Confounds).

## Why a seeded titration, and why volume
Unseeded runs only show 50/50 (`soai_check.py` S1: 46 of 91 S), which every hypothesis predicts. A seed
ladder measures P(seed hand wins) against m, and the probit fit returns m\* at each volume. Volume is the
one variable that changes N_d under H1 and leaves it alone under H0. A 100× span moves m\* by 10× under
H1, more than Stage 1's calibration error (1.91×, below).

Prior hint: in the uncatalysed aldol (Mauksch et al. 2010), RMS ee fell with volume with a log–log slope of −0.33,
against −0.5 if each molecule votes (`blind_solution.py` A2, failed). That study was unseeded and small, so it
does not decide the question.

## Chemistry and conditions
- **Reaction:** 2-(tert-butylethynyl)pyrimidine-5-carbaldehyde (C₁₁H₁₂N₂O, MW 188) + iPr₂Zn → the
  5-pyrimidyl alkanol (C₁₄H₂₀N₂O, MW 232; `soai_check.py` 1c).
- **Planned concentrations:** aldehyde 0.05 M, iPr₂Zn 0.10 M (2 equiv), dry toluene, 0 °C.
  - The lab may substitute its own validated Soai conditions before Stage 0.
  - After Stage 0 starts, conditions are frozen.
- **Volumes:** 0.1, 1 and 10 mL, a 100× span.
  - Silanized borosilicate vials of one shape and one glass lot.
  - Equal fill fraction, magnetic stirring, one stir-bar shape scaled to the vessel.
- **Addition:** iPr₂Zn is added by syringe pump, subsurface, over the same fixed time at every volume.
  An internal thermocouple logs temperature, and the excursion must stay ≤ 1 °C.
- **Seed:** enantiopure (≥ 99% ee) (S)- or (R)-alkanol.
  - Made by gravimetric serial dilution in dry toluene, in pre-conditioned silanized vials.
  - Delivered as one fixed aliquot (5 µL) at every volume, 1 min before iPr₂Zn addition begins.
  - At the anchor, 6×10⁴ molecules in 5 µL is about 2×10⁻¹⁴ M.
  - Each day, one dilution series feeds every volume and level. Shared levels (4.2×10⁴ at 1 and 10 mL;
    1.68×10⁵ at 0.1 and 1 mL) come from the same vial.
- **Quench and work-up:**
  - The product zinc alkoxide precipitates selectively at higher conversion, which changes the solution ee
    (Buono, Iwamura & Blackmond 2004).
  - A precipitate is also a nucleation-limited local event, the kind H0 describes. Every run is therefore
    quenched at a fixed time t_q, before turbidity.
  - Use the whole vessel (the lab's usual Soai quench and extraction), so any solid and the solution are
    scored together.
- **Readout:**
  - The hand is the sign of the alkanol's ee by chiral HPLC, using a method with baseline resolution of the
    enantiomers, run on racemic and enantiopure standards each day.
  - |ee| < 2% is recorded as **undecided** and is excluded from the fit. Every undecided run is reported.
  - Conversion at t_q is recorded for every run.
- **Reagent use:** about 7 g of aldehyde and about 75 mmol of iPr₂Zn over all 209 runs. The 63 runs at
  10 mL use about 94 mg of aldehyde each.

## Stages and arms
Ladder values (molecules) come from `soai_volume_power.py`. After Stage 1b, each ladder entry is multiplied
by m̂_ref / 60,000, where m̂_ref is the m\* that Stage 1b measures at 1 mL. The power does not depend on that
scale.

| stage / arm | what | runs |
|---|---|---|
| **Stage 0** setup | unseeded runs at each volume: temperature log, turbidity onset t_turb (laser and photodiode, or video), mixing time (dye tracer in an identical vessel, at the Arm V and Arm M stir rates). Set **t_q** = 0.5 × the earliest t_turb at any volume, one value for all volumes | 3 per volume (9) |
| **Stage 1a** bracket at 1 mL | seeds 6×10², 6×10³, 6×10⁴, 6×10⁵, 6×10⁶ molecules, 2 S + 2 R each; a rough probit fit gives m̂_1a | 20 |
| **Stage 1b** calibration at 1 mL | the 1 mL ladder × m̂_1a / 6×10⁴, 6 runs per level (3 S + 3 R); a probit fit gives m̂_ref, to within a factor of 1.91 (95%) | 24 |
| **Arm V** main test | 3 volumes × 4 levels × 6 runs (3 S + 3 R): 0.1 mL 1.33, 3.09, 7.21, 16.8 ×10⁴; 1 mL 4.20, 6.67, 10.6, 16.8 ×10⁴; 10 mL 4.20, 9.79, 22.8, 53.1 ×10⁴ (× m̂_ref / 6×10⁴) | 72 |
| **Arm 0** chance zone | unseeded, the same t_q | 8 per volume (24) |
| **Arm P** positive control | 6×10⁷ × m̂_ref / 6×10⁴ molecules (≥ 300 m\* under either hypothesis), 2 S + 2 R per volume | 12 |
| **Arm M** mixing | the 10 mL ladder at a stir rate giving ≥ 3× the Stage 0 mixing time; m\* compared with Arm V at 10 mL | 24 |
| **Arm S** surface | the 1 mL ladder with silanized glass beads that raise the wetted area 20× at the same 1 mL of liquid; compared with Arm V at 1 mL | 24 |
| | **total** | **209** |

**Order and blinding.**
- The Arm V, 0, P, M and S runs are put in one randomized run sheet, in day blocks. Each block mixes volumes,
  seed hands and levels. The sheet is committed to git before the first Arm V run.
- A third person codes the vials. The HPLC analyst does not know a vial's volume, seed hand or level.
- The key is merged only after every Arm V run is scored.
- Stage 0 and Stage 1 data are not pooled into the test. There is no interim look at the slope.

## Analysis and decision rule (fixed in advance)
**Model.** One probit fit over all 72 Arm V runs:

  P(product S) = Φ(σ·m·exp(−θ₀ − β·x_j) + b_j),  x_j = ln(V_j / 1 mL).

- σ = +1 for an S seed and −1 for an R seed.
- m\*_j = exp(θ₀ + β·x_j).
- b_j is a stray bias at each volume, in Δ units, positive toward S.

With mirror seeds, b_j is orthogonal to (θ₀, β) at b = 0. The power calculation leaves b_j out; the
margin at n = 6 (0.93 against 0.80) covers the three extra parameters.

**Tests.**
- Likelihood-ratio tests of β = 0 and of β = 0.5, each at α = 0.05 (χ²₁ > 3.84), profiling over θ₀ and b_j.
- β̂ is reported with its 95% profile-likelihood interval in every case.
- The script's simulated test size is 0.059, slightly above 0.05. Report that too.

**Power** (`soai_volume_power.log`, 4000 simulations per row):

| runs per level | runs per volume | P(reject flat \| √V true) | P(reject √V \| flat true) | P(reject flat \| flat true) |
|---|---|---|---|---|
| 2 | 8 | 0.514 | 0.527 | 0.067 |
| 4 | 16 | 0.799 | 0.814 | 0.061 |
| **6** | **24** | **0.932** | **0.937** | 0.059 |

- If the true m\* at 1 mL is 0.5× or 2× the calibrated value, the power is 0.916/0.906 or 0.927/0.940.
- Stage 1b's 1.91× error lies inside that range.
- A side-arm ladder of 24 runs detects a change in m\* of 3.7× at 1 mL and 5.6× at 10 mL (80% power, α = 0.05).

**Decision rule.**
1. Reject flat, keep √V → **one well-mixed pool** (H1). The solution reaction reaches the best case of
   Eq. 13 at these scales.
2. Reject √V, keep flat → **local or early events decide** (H0).
3. Reject both → neither model. Read β̂:
   - 0 < β̂ < 0.5 together with Arm M positive → a pool limited by mixing.
   - β̂ near ⅓ together with Arm S positive → a pool set at the walls.
   - β̂ > 0.5 → a bias that grows with V, or seed loss to the walls in the small vessels. Check Arm 0, b_j
     and Arm S.
4. Reject neither → **inconclusive** at this n. Report the interval. No runs are added afterward.

**Validity checks.** Arm V is interpreted only if all of these hold:
- Arm 0 gives 7–17 of 24 for one hand.
- The fitted b_j are jointly consistent with 0 (χ²₃ at α = 0.05).
- Arm P gives 12 of 12 following the seed.
- Mean conversion at t_q agrees within 10% across the three volumes.

## Stop rules
- **Arm 0 outside 7–17 of 24, or a significant b_j:** stop. Purify the reagents, change the lot, and restart
  Arm V. The pooled test catches a stray Δ = 1 bias with probability 0.926.
- **Turbidity before t_q** in more than 20% of runs at any volume, or **undecided** runs (|ee| < 2%) above 10%
  at any volume: stop and reset t_q (earlier for turbidity, later for undecided). If no t_q satisfies both,
  the conditions fail: change the solvent or temperature and restart from Stage 0. A stopped Arm V is reported,
  not pooled.
- **Stage 1a bracket fails** (every level follows the seed, or none beats chance): extend the bracket 10× in
  the needed direction, at most twice. Then stop and report m\* as outside 6×10⁰–6×10⁸.
- **Arm P fails** (any run against the seed): stop, because seed delivery is broken.
- **Mixing time > 10% of t_q** at the Arm V stir rate at any volume: fix the stirring before Stage 1.
- **Kinetics mismatch** (conversion at t_q differs by more than 10% across volumes, or the temperature
  excursion exceeds 1 °C): stop and fix heat removal or addition before continuing.
- **Safety:** iPr₂Zn is pyrophoric. Handle it under argon (glovebox or Schlenk line), and quench the
  10 mL runs (1 mmol iPr₂Zn) slowly.

## Confounds, each with a control
| confound | how it would fake or hide a slope | control |
|---|---|---|
| **Chiral impurities** (Singleton & Vo 2002: 48 unseeded runs, all optically active, traced to optically active impurities) | at a fixed concentration, the impurity count grows ∝ V. Its Δ grows as √V under H1 and as V under H0, so the bias b_j depends on V | one purified lot of aldehyde, toluene and iPr₂Zn for every run; acid-washed, silanized, oven-dried glass; mirror seeds, which cancel b_j in m\* to first order; b_j fitted at each volume; Arm 0 (pooled 50/50 test, power 0.926 against a Δ = 1 bias) |
| **Oligomers as voting units** (Athavale, Simon, Houk & Denmark 2020: oligomeric zinc alkoxide aggregates) | if the unit is a k-mer, m\* = √(k·N_d). This changes the absolute m\*, not β, as long as k is the same at every volume | fixed concentration and temperature at all volumes, so speciation is the same; the seed is added as the alkanol and converted in situ the same way everywhere; the kinetics-match stop rule. A concentration series would measure k, but it is not part of this test |
| **Mixing-driven stochasticity** (Valkai & Horváth 2018: in the arsenous acid–periodate "crazy clock" reaction, not a Soai system, random reaction times came from imperfect mixing and dead volumes) | if the 10 mL vessel does not mix within the time the hand is set, smaller local pools decide there, so β falls below 0.5 | mixing time measured in Stage 0 and required to be ≤ 10% of t_q; subsurface syringe-pump addition with matched timing; Arm M (≥ 3× slower mixing at 10 mL detects a change ≥ 5.6×) |
| **Surface-to-volume ratio** | S/V falls 4.6× from 0.1 to 10 mL. Seed lost to the walls raises the apparent m\* in small vessels, biasing β down. A pool set at the walls gives β = ⅓ | one glass lot, silanized; equal fill fraction and shape; pre-conditioned dilution vials; Arm S (20× the wetted area at 1 mL; a wall pool predicts √20 = 4.5×, and the arm detects ≥ 3.7×) |
| **Heat removal** | larger vessels lose heat more slowly, so the exotherm warms them more and changes the kinetics | internal thermocouple, excursion ≤ 1 °C, the same addition time at every volume, the kinetics-match stop rule |
| **Seed-delivery loss** | a loss factor κ shared by all volumes cancels in β. A loss that depends on concentration does not cancel, because the ladders use different dilutions | one dilution series per day for all volumes; the levels shared between volumes come from the same vial; Arm P checks that delivery works at all |

## What each outcome means
| result | meaning for the paper (`journals/ija_draft.md`) |
|---|---|
| H1: β̂ ≈ 0.5, flat rejected | first evidence that a solution reaction pools as Eq. 13 assumes, at least up to 10 mL. Section 3.7's "still untested" becomes a measured case, and the scale rule's N stops being only a minimum for this system |
| H0: β̂ ≈ 0, √V rejected | the hand is set by a few early units even in solution, as in the crystallizations of Section 3.7. ε ≫ 1/N, so Eq. 13's N is a strict minimum and the true N is larger |
| 0 < β̂ < 0.5 with Arm M or Arm S positive | pooling is partial and set by mixing or by the walls. That scaling replaces V in Eq. 13 for a real ocean |
| inconclusive | report the interval on β; that interval sets the n a second lab would need |
| validity check fails | no result on β; report which check failed (impurity, delivery, kinetics) |

## If the Soai reaction fails
The Soai reaction fails here if no workable t_q exists, or if Stage 1 cannot find m\*. The same design
(seed ladders over a 100× volume span) then moves to one of these runners-up:
- **Allylboration of 1,2,3-triazolic aldehydes with triallylborane** (Mikhailov, Gurskii, Kurbangalieva &
  Gridnev 2024; the brief called this "Gridnev 2024").
  - The paper reports spontaneous chirality generation and computes three amplification routes.
  - Caution: its seeded runs with the (+)-R alcohol still gave both hands (7 S : 3 R and 6 S : 4 R). The seed may
    not steer well enough for a titration. Run a Stage 1a bracket first.
- **The catalyst-free aldol** (Valero, Ribó & Moyano 2014).
  - ee appears early in about half the runs, but few runs are still non-racemic at 24 h.
  - The ee is transient, so t_q would sit at the early maximum. Read the paper's time course before planning.

## Limits and practical notes
- This needs an air-free synthesis lab with chiral HPLC. It is a proposal for a mentor or collaborator,
  not something to run at home.
- The 1 mL anchor is an assumption. Stage 1 replaces it with a measurement, and the robustness runs cover a
  residual error of 2×.
- Three volumes give one slope. A fourth volume would test curvature but is not in the power budget.
- Seed amounts of 10⁴–10⁶ molecules (20–900 zmol) cannot be checked by any assay. The design relies on
  shared dilution series, overlapping levels and Arm P, not on measuring what was delivered.
- A result holds for this aldehyde, these conditions and 0.1–10 mL. It says nothing directly about water, or
  about a mixed ocean.
- The analysis is pre-registered in `soai_volume_power.py`. Its self-check asserts:
  - the anchor;
  - that the ladders cover both predictions;
  - that the fit recovers β = 0 and β = 0.5;
  - power ≥ 0.80 both ways;
  - a test size of 0.02–0.09;
  - Stage 1b's calibration within 2×;
  - Arm 0's power ≥ 0.9.

## References
- Singleton & Vo 2002, *J. Am. Chem. Soc.* 124:10010, doi:10.1021/ja027129o
- Singleton & Vo 2003, *Org. Lett.* 5:4337, doi:10.1021/ol035605p
- Buono, Iwamura & Blackmond 2004, *Angew. Chem. Int. Ed.* 43:2099, doi:10.1002/anie.200353086
- Athavale, Simon, Houk & Denmark 2020, "Demystifying the asymmetry-amplifying, autocatalytic behaviour of the Soai reaction through structural, mechanistic and computational studies", *Nat. Chem.* 12:412, doi:10.1038/s41557-020-0421-8
- Valkai & Horváth 2018, "Imperfect mixing as a dominant factor leading to stochastic behavior: a new system exhibiting crazy clock behavior", *Phys. Chem. Chem. Phys.* 20:14145, doi:10.1039/c8cp01156g
- Mauksch, Wei, Freund, Zamfir & Tsogoeva 2010, *Orig. Life Evol. Biosph.* 40:79, doi:10.1007/s11084-009-9177-2
- Mikhailov, Gurskii, Kurbangalieva & Gridnev 2024, "Exploring border conditions for spontaneous emergence of chirality in allylboration of 1,2,3-triazolic aldehydes", *Int. J. Mol. Sci.* 25:11273, doi:10.3390/ijms252011273
- Valero, Ribó & Moyano 2014, *Chem. Eur. J.* 20:17395, doi:10.1002/chem.201404497
