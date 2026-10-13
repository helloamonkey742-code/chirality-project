# Vertical mixing from the tidal energy budget (tidal_mix.py, 2026-10-12)

Script: `tidal_mix.py` (predictions committed before the run, d6854b1). Log: `tidal_mix.log` (N = 200000, 10 s wall).
Idea: Prof. Francis Nimmo (UCSC), after Munk & Wunsch (1998). No existing script, log or md was edited.
`journals/ija_draft.md` was only read. One crash-only smoke run (TIDAL_N=2000, output to /private/tmp, not kept) was
made after the predictions commit. Its numbers are not used, and no code changed between it and the full run.

## The method in plain words

A stratified ocean stores potential energy, because the light water sits on top of the heavy water. Mixing it means
lifting heavy water up, and that work has to be paid for by some energy source. In an ice-covered ocean the main
mechanical source is tidal power. So the time to mix the ocean is roughly the stored energy divided by the tidal
power that goes into mixing.

- For a linear density profile, the stored energy per unit area is PE = ρN²H³/12. Here N² = Δb/H is the
  stratification and Δb is the buoyancy difference from top to bottom. The derivation is in the script docstring.
- Only a fraction Γ of the tidal power P goes into mixing: Γ = 0.2, range 0.1–0.3 (Osborn 1980). The rest becomes
  heat.
- τ_mix = PE·A/(ΓP), where A is the ocean area. The diffusivity is K_z = Γε/N² = ΓP/(ρAΔb), and it does not depend on
  depth.
- Cross-check: K_z·τ_mix = H²/12 exactly, and it holds to 1.000000 for every world.
- The band runs from K_low (Γ = 0.1 with the lowest published P) to K_high (Γ = 0.3 with the highest), with
  K_mid = 0.2·P_mid.
- Two layers are estimated separately.
  - **Interior.** Δb = g·2×10⁻⁶. This is the Europa density contrast of Ashkenazy & Tziperman (2021), applied to every
    world and flagged.
  - **Lid.** This is the stable layer under the ice of the worlds that have an S* lid. Δb = g·|α|·ΔT = g·4×10⁻⁵·2 K
    (Zeng & Jansen 2021). The lid thickness follows their scaling: 30 m at κ = 10⁻⁷, rising with K, adjusted for the
    world's heat flux and capped at the ocean depth.
- The band was then compared with the D_z prior: the effective prior is log-uniform from 10⁻⁹ to 10⁻³ m² s⁻¹, and the
  raw one from 10⁻¹⁰.
- The success shares were rerun with D_z restricted to each world's interior band. The model is `multiwall.multi`,
  imported unedited, with seeds 0–3 and seedfix's world-0 stream.

## Inputs (every DOI checked on Crossref 2026-10-12, titles matched)

| Source | DOI | Where | Used for |
|---|---|---|---|
| Chen, Nimmo & Glatzmaier 2014, Icarus 229 | 10.1016/j.icarus.2013.10.024 | Table 3 (p. 21); Table 2 (p. 15); Sec. 5.1 (p. 20) | Ocean eccentricity + obliquity dissipation (free-surface ocean, no shell); radiogenic and solid-tide heat for the total; g; ×4 obliquity for a decoupled shell |
| Hay & Matsuyama 2019, Icarus 319 | 10.1016/j.icarus.2018.09.019 | Table 3 (accepted MS p. 20) | Ocean dissipation with an ice shell (P_mid where it exists); Triton thick/thin shell 5.71×10¹⁰ / 2.89×10¹¹ W |
| Jansen et al. 2023, PSJ 4:117 | 10.3847/psj/acda95 | Sec. 4; Sec. 2 | Europa 10⁹–10¹¹ W; Enceladus 10–10⁴ W; Earth 2 TW and ε 1.5×10⁻⁹; Europa floor Q 0.02–0.1 W m⁻² |
| Rovira-Navarro et al. 2023, PSJ (via J23) | 10.3847/psj/acae9a | J23 Sec. 4 | Enceladus "most likely" 10–10⁶ W |
| Zeng & Jansen 2021, PSJ 2:151 | 10.3847/PSJ/ac1114 | Sec. 2.1, 2.2 | Lid α, ΔT, thickness scaling, Q = 0.03 W m⁻²; Enceladus E 10–10⁸ W with libration as the upper bound; κ 3×10⁻¹⁰–3×10⁻³ |
| Choblet et al. 2017, Nat. Astron. (via ZJ21) | 10.1038/s41550-017-0289-8 | ZJ21 Sec. 2.1 | Enceladus ~20 GW behind Q = 0.03 |
| Wilson & Kerswell 2018, EPSL (via ZJ21) | 10.1016/j.epsl.2018.08.012 | ZJ21 Sec. 2.2 | Libration O(0.1 GW), second-hand |
| Ashkenazy & Tziperman 2021, Nat. Commun. 12:6376 | 10.1038/s41467-021-26710-0 | rotation-scaling paragraph | Δρ/ρ = 2×10⁻⁶, ρ = 10³, c_p = 4×10³, q = 0.05 W m⁻² |
| Osborn 1980, JPO 10 | 10.1175/1520-0485(1980)010<0083:EOTLRO>2.0.CO;2 | title only | Γ ≈ 0.2 (as used by ZJ21 and J23) |
| Munk & Wunsch 1998, DSR I 45 | 10.1016/S0967-0637(98)00070-3 | title only | the energy-budget idea |
| Wunsch & Ferrari 2004, ARFM 36 (via J23) | 10.1146/annurev.fluid.36.050802.122121 | J23 Sec. 4 | Earth ε |

These were read but gave no usable number: Matsuyama et al. 2018 (10.1016/j.icarus.2018.04.013, values only in
figures), Beuthe 2016 (10.1016/j.icarus.2016.08.009) and Rovira-Navarro et al. 2019 (10.1016/j.icarus.2018.11.010).

**Rule for P.**
- Low is the smallest published value or range end.
- Mid is HM19 (with a shell) where it exists, otherwise C14 ocean eccentricity + obliquity. Triton uses the geometric
  mean of HM19's two shell cases.
- High is the largest published value, including C14's eccentricity + 4 × obliquity.
- Libration is counted only for Enceladus, because only ZJ21 gives a value.
- Total heat (C14 radiogenic + solid tides + ocean) is printed as a ceiling and not used in the band.

**Ocean power P (W) and the band.** K_z is in m² s⁻¹. "% in band" is the share of the effective prior that falls
inside the band, and "% >10⁻³" is the share of the band that lies above the prior.

| World | P low | P mid | P high | K_low | K_mid | K_high | τ_mix mid (yr) | τ_mix range (yr) | % in band | % >10⁻³ | Flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Enceladus | 1.0e1 | 1.3e1 | 1.0e8 | 6.7e-9 | 1.7e-8 | 0.20 | 2.4e8 | 21 – 6.3e8 | 86.2 | 30.8 | high end is libration |
| Europa | 1.0e9 | 6.2e9 | 1.0e11 | 1.4e-3 | 1.8e-2 | 0.43 | 2.3e3 | 94 – 2.8e4 | 0 | 100 | |
| Ganymede | 1.3e9 | 2.8e9 | 5.2e9 | 5.2e-4 | 2.2e-3 | 6.3e-3 | 1.2e4 | 4.2e3 – 5.1e4 | 4.7 | 73.9 | |
| Titan | 1.1e10 | 1.2e10 | 4.4e10 | 4.9e-3 | 1.1e-2 | 5.9e-2 | 1.0e4 | 1.8e3 – 2.2e4 | 0 | 100 | |
| Callisto | 3.8e9 | 3.8e9 | 1.5e10 | 2.1e-3 | 4.2e-3 | 2.5e-2 | 6.3e3 | 1.1e3 – 1.3e4 | 0 | 100 | single source |
| Triton | 5.7e10 | 1.3e11 | 2.9e11 | 0.16 | 0.72 | 2.4 | 83 | 25 – 373 | 0 | 100 | HM19 thin shell above C14 total |
| Titania | 1.3e5 | 1.3e5 | 5.2e5 | 2.2e-6 | 4.4e-6 | 2.6e-5 | 5.4e5 | 9.1e4 – 1.1e6 | 18.0 | 0 | single source |
| Oberon | 3.3e6 | 3.3e6 | 1.3e7 | 6.5e-5 | 1.3e-4 | 7.8e-4 | 1.8e4 | 3.1e3 – 3.7e4 | 18.0 | 0 | single source |
| Umbriel | 2.8e3 | 2.8e3 | 8.2e3 | 1.4e-7 | 2.8e-7 | 1.2e-6 | 3.8e6 | 8.5e5 – 7.6e6 | 15.8 | 0 | single source |
| Ariel | 1.4e3 | 1.4e3 | 1.6e3 | 6.2e-8 | 1.2e-7 | 2.1e-7 | 1.9e7 | 1.2e7 – 3.8e7 | 8.7 | 0 | single source |
| Dione | 6.7e3 | 1.2e4 | 1.3e4 | 3.7e-7 | 1.3e-6 | 2.2e-6 | 8.6e6 | 5.1e6 – 3.1e7 | 12.9 | 0 | |
| Miranda | 2.6e5 | 2.6e5 | 1.0e6 | 2.4e-4 | 4.7e-4 | 2.8e-3 | 5.6e4 | 9.4e3 – 1.1e5 | 10.5 | 41.8 | single source |
| Mimas | 2.6e7 | 2.6e7 | 3.4e7 | 4.2e-2 | 8.3e-2 | 0.16 | 64 | 33 – 128 | 0 | 100 | C14 bracketed value; no shell |
| Pluto | — | — | — | — | — | — | — | — | — | — | not sourced: full prior kept |
| Early Earth | — | — | — | — | — | — | — | — | — | — | Earth calibration only: own prior kept |

**Percentile of K_low / K_mid / K_high in the effective prior.**
- Enceladus: 14 / 21 / 100.
- Titania: 56 / 61 / 74.
- Oberon: 80 / 85 / 98.
- Umbriel: 36 / 41 / 52.
- Ariel: 30 / 35 / 39.
- Dione: 43 / 52 / 56.
- Miranda: 90 / 95 / 100.
- Ganymede: 95 / 100 / 100.
- Europa, Titan, Callisto, Triton and Mimas: 100 throughout.

**Timescales.** The decision time 1/k_rac is 3.47 Myr.
- τ_mix at mid is shorter than the decision time for Europa, Ganymede, Titan, Callisto, Triton, Mimas, Miranda, Oberon
  and Titania. The ratios run from 2×10⁻⁵ (Mimas) to 0.16 (Titania).
- It is longer for Enceladus (70×), Ariel (5.5×), Dione (2.5×) and Umbriel (1.1×).
- Compared with the youngest model ocean age, every sourced world mixes in under a fifth of it, except Enceladus
  (242× longer).

**The lid (worlds with an S* lid).** All of P is assumed to work on the lid, as in the ZJ21 convention.

| World | Q (W m⁻²) | K_low | K_mid | K_high | H_lid lo / mid / hi (m) | τ lo / mid / hi (yr) | Q source |
|---|---|---|---|---|---|---|---|
| Enceladus | 0.030 | 1.7e-10 | 4.4e-10 | 5.0e-3 | 30 / 30 / 4.0e4 | 1.4e4 / 5.5e3 / 835 | ZJ21 Sec. 2.1 |
| Europa | 0.050 | 3.6e-5 | 4.4e-4 | 1.1e-2 | 6.4e3 / 8.0e4 / 1.2e5 | 3.1e3 / 3.8e4 / 3.8e3 | AT21 q |
| Dione | 0.003 | 9.2e-9 | 3.2e-8 | 5.5e-8 | 328 | 3.1e4 / 8.8e3 / 5.2e3 | C14 total / A (derived) |
| Miranda | 0.002 | 5.9e-6 | 1.2e-5 | 7.1e-5 | 3.2e4 / 6.5e4 / 1.0e5 | 4.7e5 / 9.4e5 / 3.7e5 | C14 total / A (derived) |
| Mimas | 1.541 | 1.0e-3 | 2.1e-3 | 4.0e-3 | 6.1e3 / 1.2e4 / 2.4e4 | 94 / 188 / 362 | C14 total / A (derived) |

Check against ZJ21: the Enceladus lid at Γ = 0.2 over E = 10–10⁸ W gives 3.4×10⁻¹⁰ to 3.4×10⁻³ m² s⁻¹. ZJ21 give
3×10⁻¹⁰ to 3×10⁻³, so the match is within a factor of 1.12 at both ends.

## Shares (P_pved / P_chance / P_none %, seed 0, no lid | final)

The full-prior run reproduced multiwall main, with seedfix's Enceladus row and seedfix's seed 1–3 lines, **exactly**
for all 15 worlds. Across seeds 1–3 the band P_pved stays within 0.37 points of seed 0 for every world (Umbriel is
the largest).

| World | Full prior, no lid | Full prior, final | Interior band, no lid | Interior band, final | Δ P_pved no lid / final |
|---|---|---|---|---|---|
| Early Earth ocean | 99.9 / 0.0 / 0.1 | 99.9 / 0.0 / 0.1 | (full prior kept) | (full prior kept) | 0 / 0 |
| Umbriel | 53.7 / 14.9 / 31.4 | same | 63.7 / 23.1 / 13.2 | same | +10.1 / +10.1 |
| Oberon | 51.4 / 11.7 / 36.9 | same | 75.5 / 10.5 / 14.0 | same | +24.1 / +24.1 |
| Titania | 51.3 / 11.6 / 37.1 | same | 73.1 / 16.2 / 10.7 | same | +21.8 / +21.8 |
| Ariel | 49.6 / 13.5 / 36.9 | same | 44.3 / 17.8 / 37.9 | same | −5.3 / −5.3 |
| Ganymede | 44.1 / 3.4 / 52.5 | same | 80.8 / 2.6 / 16.6 | same | +36.6 / +36.6 |
| Callisto | 43.5 / 3.8 / 52.7 | same | 79.8 / 2.1 / 18.1 | same | +36.3 / +36.3 |
| Dione | 42.0 / 10.5 / 47.5 | 26.2 / 6.6 / 67.2 | 51.6 / 17.1 / 31.3 | 32.2 / 10.7 / 57.1 | +9.6 / +6.0 |
| Pluto | 41.1 / 6.1 / 52.8 | same | (full prior kept) | (full prior kept) | 0 / 0 |
| Europa | 39.6 / 4.9 / 55.5 | 30.1 / 3.7 / 66.2 | 79.5 / 3.2 / 17.3 | 60.4 / 2.4 / 37.1 | +39.9 / +30.4 |
| Triton | 37.0 / 4.7 / 58.3 | same | 78.9 / 2.2 / 18.9 | same | +41.9 / +41.9 |
| Titan | 35.8 / 2.5 / 61.7 | same | 80.8 / 1.7 / 17.5 | same | +45.0 / +45.0 |
| Miranda | 27.1 / 8.4 / 64.5 | 9.0 / 2.7 / 88.3 | 68.8 / 16.2 / 15.0 | 22.7 / 5.3 / 72.0 | +41.7 / +13.7 |
| Enceladus | 26.4 / 9.8 / 63.9 | 0.0 / 0.0 / 100.0 | 41.3 / 12.2 / 46.5 | 0.0 / 0.0 / 100.0 | +15.0 / 0.0 |
| Mimas | 17.3 / 6.1 / 76.6 | 5.6 / 2.0 / 92.5 | 63.4 / 13.7 / 22.8 | 20.4 / 4.5 / 75.1 | +46.1 / +14.8 |

What the restriction does:
- It raises P_pved for 12 of the 13 sourced worlds. Most worlds reach 75–81% without a lid.
- The gain comes almost entirely from P_none falling. Draws with D_z too small for the walls to merge in time are
  removed.
- Ariel is the exception (−5.3). Its band, 6×10⁻⁸ to 2×10⁻⁷, sits in the part of the prior where walls are slow.
- With a lid, Enceladus stays at 0. The band moves interior mixing only, and the lid rule in `multiwall` is
  unchanged.

## Predictions: 10 held, 1 failed (full list in the tidal_mix.py docstring)

**P1–P9 and P11 held.** In order, these were:
- an exact regression;
- the predicted sign of K_mid against 10⁻³;
- Europa K_low of 1.4×10⁻³;
- ZJ21's lid range reproduced;
- the lid τ overlapping ZJ21's;
- τ against the decision time;
- τ against the ocean age;
- under 20% of the prior inside the band (86% for Enceladus);
- higher shares for the six high-K worlds;
- seed stability within 0.5 points.

**P10 failed.** I predicted lower shares for Ariel, Umbriel and Dione. Ariel fell (−5.3), but Umbriel (+10.1) and
Dione (+9.6) rose. Their bands, at about 10⁻⁷ to 10⁻⁶, already clear the low part of the prior where walls cannot
merge in time. Only Ariel's band does not.

## Nimmo's caveat

Earth's deep stratification is kept up by a continuous supply of cold, dense water sinking at the poles, and Nimmo
noted that this is unlikely in the outer solar system. Its effect is different for the two layers.
- **Interior.** Nothing obvious renews the stratification. The real Δb is then likely smaller than the Europa value
  used here, which would make K larger and τ shorter for the same power.
- **Not a lower bound, though.** The estimate assumes the tidal power is spent evenly through the water. If it is
  spent in boundary layers at the ice or the seafloor (J23 Sec. 4), interior K is smaller. So the interior K is an
  order-of-magnitude estimate, not a bound.
- **Lid.** The stable layer under the ice of fresh, cold oceans (α < 0) is made by the heat flux itself. It needs no
  outside water source, so the caveat does not weaken it.
- **Effect on the model.** The caveat therefore supports keeping the lid and the interior separate, which is what the
  model does.

## Flags

- **Pluto: not sourced.** Pluto is in neither C14 nor HM19, so it keeps the full prior.
- **Early Earth: not sourced.** J23 give Earth's ε but no N was read, so Early Earth keeps its own prior.
- **Single source.** Callisto, Titania, Oberon, Umbriel, Ariel and Miranda have only C14 (low = mid), so their bands
  are narrow because of the sourcing, not because the physics is known.
- **C14 has no ice shell.** It models a free-surface ocean. Mimas's obliquity value is bracketed in the paper as
  outside the model's validity; it is used here and flagged.
- **Triton.** HM19's thin-shell value, 2.9×10¹¹ W, is above C14's total heat (1.5×10¹¹ W).
- **Δρ/ρ.** The value 2×10⁻⁶ comes from one Europa paper and is applied to every world.
- **Derived lid heat flux.** For Dione, Miranda and Mimas, Q is derived as C14 total / area, not published.
- **Extrapolation above 10⁻³.** Band ends above 10⁻³ are not clipped, so for six worlds and part of Enceladus and
  Miranda the shares run the model past the D_z range it was built for.
- **Γ.** The Osborn (1980) and Munk & Wunsch (1998) full texts were not read. Γ ≈ 0.2 is taken as used by ZJ21 and
  J23.

## Reply to Prof. Nimmo (draft)

Thank you for the energy-budget suggestion. I divided the stored potential energy by the published ocean tidal power
(Chen et al. 2014; Hay & Matsuyama 2019, with a mixing efficiency of 0.2). This puts the interior mixing time at
tens to tens of thousands of years for Europa, Ganymede, Titan, Callisto and Triton, and at
millions of years or more for Enceladus, Ariel, Umbriel and Dione. For most of the larger moons the implied
diffusivity is above the top of the range I had been sampling (10⁻³ m² s⁻¹). Restricting the sampling to each moon's
energy-budget range raises the chance that the weak-force bias sets the hand from about 35–50% to 75–80% for most of
them, but it leaves Enceladus at zero once its cold, fresh lid is included. Your point about cold-water production
seems to me to separate the two layers: the interior has nothing to keep it stratified, while the lid under the ice
is maintained by the heat flux itself. Would you consider that split reasonable?

## Proposed draft sentences (NOT applied to journals/ija_draft.md)

1. Section 2.6, after the D_z prior: "An energy budget gives a physical check on this range. Dividing the potential
   energy of the stratification, ρN²H³/12, by the ocean tidal power (Chen et al., 2014; Hay and Matsuyama, 2019),
   with a mixing efficiency of 0.2 (Osborn, 1980; Munk and Wunsch, 1998), gives interior diffusivities above
   10⁻³ m² s⁻¹ for Europa, Ganymede, Titan, Callisto, Triton and Mimas, and 10⁻⁸ to 5×10⁻⁴ m² s⁻¹ (mid estimates) for Enceladus and
   the smaller moons (`tidal_mix.py`)."
2. Section 3, sensitivity: "With D_z restricted to each moon's energy-budget range, the success share without a lid
   rises for 12 of 13 sourced worlds, to 75–81% for most, and falls for Ariel (49.6% to 44.3%). Enceladus with a lid
   stays at 0%."
3. Section 4, limitations: "The energy-budget range assumes tidal power is dissipated throughout the ocean and that the
   interior density contrast is Europa's (Ashkenazy and Tziperman, 2021). Earth's deep stratification is renewed by
   cold-water production, which is unlikely in icy-moon oceans (F. Nimmo, pers. comm.). This would raise interior
   mixing, while dissipation concentrated in boundary layers would lower it."
4. Section 4: "For Pluto no ocean tidal power has been published, and its full D_z range is kept."
