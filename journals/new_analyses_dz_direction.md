# Does vertical mixing rise with depth? Literature constraints for H25 (2026-10-12)

Question (H25, `ija_draft.md` line 201; assumption 7, line 535): every reversal result samples D_z as a smooth
exponential that rises with depth, bottom/top ratio ρ log-uniform on 1–100 (`sweep3.py:157`,
`rho = 10 ** rng.uniform(0, 2)`). Is there an observational or modelling constraint on the direction for Europa,
Enceladus, Ganymede, Titan, with Earth for comparison?

Short answer: rising with depth is supported for Earth's ocean below the surface mixed layer (by measurement).
For the four moons, no observation exists. The published models give three outcomes:

- a uniformly convective interior (ρ ≈ 1);
- a stratified layer under the ice over a convective bottom (a step rise, which is the lid case);
- dissipation in boundary layers at the ice, the seafloor, or both (either direction).

The one-sided 1–100 prior is not supported for the moons.

## What ρ ≤ 1 does in the code (checked, not rerun)

`sweep3.py:outcome3` pins a wall only if `lr >= 4*H/r_w` with `lr = ln(rho)`. For ρ ≤ 1, lr ≤ 0, so no wall is ever
pinned. `sweep3.py:check` asserts `w_f == w_e & ~pinned`, so an unpinned draw gets exactly the **sweep2 step e**
outcome (all additions a–e on, D_z uniform). That is not the step-0 "uniform-D_z base model" the draft names in
assumption 7.

| World (P_pved %) | step 0 base (`sweep2.log`) | step e = ρ ≤ 1 (`sweep2.log`) | current ρ = 1–100 (`multiwall.log:541–556`; low = cascade ×10 or sweep3 main, high = multiwall) |
|---|---|---|---|
| Early Earth | 74.7 | 74.7 | 67.0–99.9 |
| Ganymede | 28.3 | 29.9 | 33.5–44.1 |
| Titan | 24.3 | 25.5 | 25.5–35.8 |
| Europa, no lid / lid | 28.6 / 21.7 | 30.4 / 23.1 | 29.0–39.6 / 22.0–30.1 |
| Enceladus, no lid / lid | 21.8 / 0.0 | 24.5 / 0.0 | 19.9–28.3 / 0.0 |

So "if D_z is strongest near the top, the base model applies" is approximately right, but the exact fallback is
step e. The step-0 figures understate it by up to 2.7 points (Enceladus, no lid).

## Method

- Sources were found through Crossref bibliographic search. Each DOI was checked with
  `curl -s https://api.crossref.org/works/<doi>`. The title matched for every DOI cited here (JSON in
  `/private/tmp/claude-501/dzdl/`).
- Abstracts and full texts came from publisher pages only: nature.com, journals.ametsoc.org, iopscience,
  the Annual Reviews CC-BY PDF, and the Crossref abstract field. Full texts are in `/private/tmp/claude-501/dzdl/ft/`
  and were read with `python3 -I`.
- These were blocked and not bypassed: Wiley/AGU full text (403), science.org (403), and NASA ADS abstract pages
  (405). Elsevier/Icarus papers have no Crossref abstract.
- Each claim is marked:
  - **[FT]**: read in the full text;
  - **[AB]**: abstract only;
  - **[2nd]**: reported in another paper's full text, not checked at source;
  - **[T]**: Crossref title only. Never used as evidence.

## Table of constraints

ρ > 1 means D_z rises with depth. "Uniform" means a convective interior with ρ ≈ 1 outside thin boundary layers.

| World | Constraint on D_z with depth | Mechanism | Source and DOI | Verified how | Strength of evidence |
|---|---|---|---|---|---|
| Earth | **Up** (below mixed layer). Mean K: 0.3×10⁻⁴ m² s⁻¹ from mixed layer to 1000 m; 4.3×10⁻⁴ below 1000 m. Ratio ≈ 14. | Breaking internal waves from tides over rough topography and from winds. Dissipation is elevated near rough bottoms. | Waterhouse et al. 2014, 10.1175/JPO-D-13-0104.1 | Crossref OK; [FT] (Table 2 and text) | Strong. Global compilation of more than 5200 microstructure profiles. |
| Earth | **Up over rough topography; ≈ uniform over smooth plains.** Weak (~0.1×10⁻⁴) at all depths over smooth abyssal plains. Above 5×10⁻⁴ in the bottom 150 m over the Mid-Atlantic Ridge. | Topographic internal-wave breaking | Polzin et al. 1997, 10.1126/science.276.5309.93 | Crossref OK; [AB] | Strong (direct microstructure), regional |
| Earth | **Up, ratio ~100 near the bottom.** 2–4 cm² s⁻¹ at 500 m above abyssal hills; ~10 cm² s⁻¹ nearer the bottom; ~0.1 cm² s⁻¹ typical of the interior. | Tidal flow over rough bathymetry | Ledwell et al. 2000, 10.1038/35003164 | Crossref OK; [AB] (nature.com page) | Strong (tracer release plus microstructure), one basin |
| Earth | **Up** in tropical and subtropical waters | Finescale shear and strain inference | Kunze et al. 2006, 10.1175/JPO2926.1 | Crossref OK; [AB] | Moderate (indirect method) |
| Earth | **Up toward the bottom.** Mixing is "more vigorous close to the ocean bottom boundary than above". The model profile decreases toward the surface. | Bottom-intensified mixing over sloping topography | Ferrari et al. 2016, 10.1175/JPO-D-15-0244.1 | Crossref OK; [AB] and [FT] (figure captions) | Strong as synthesis; model setup, not new data |
| Earth | **Down at the very top.** A distinct, strongly mixed surface boundary layer sits above the interior. This makes the full-depth profile non-monotonic (high, then low in the pycnocline, then rising to the bottom). | Wind and buoyancy-driven surface boundary layer | Large et al. 1994, 10.1029/94RG01872; pycnocline value from Ledwell et al. 1993 (already cited) | Crossref OK; [AB] | Strong qualitatively. Magnitude not taken from the abstract. |
| Enceladus, Europa | **Uniform.** Bottom heating gives "an unstratified, turbulent, geostrophically controlled ocean" with Taylor columns from the seafloor to the ice. | Rotating convection | Bire et al. 2022, 10.1029/2021JE007025 | Crossref OK; [AB] (Crossref) | Moderate (models only) |
| Enceladus | **Uniform at high salinity; step up at low salinity.** High S: "virtually unstratified, with convection throughout". Below ~20 g/kg: a stratified upper layer (diffusion only) over a convective bottom layer. Turbulent κ_z estimated at 3×10⁻¹⁰–3×10⁻³. | Sign of thermal expansion under ice sets whether convection reaches the ice | Zeng and Jansen 2021 (already cited); Zeng and Jansen 2024 (already cited) | Crossref OK; [FT] | Moderate (models), well matched to the draft's lid test |
| Enceladus | **Step up at low salinity; can stratify throughout.** An inversely stratified layer under the ice over a convective, nearly homogeneous interior. At κ_z = 10⁻³ and S ≤ 10 g/kg the whole ocean stratifies. κ_z "may exhibit spatial and temporal variations … which may be concentrated within boundary layers at the ocean top and bottom". κ_z is held uniform. | Inverse stratification; tidal and libration dissipation in boundary layers | Ames et al. 2025 (already cited), citing Simmons et al. 2004, 10.1016/S1463-5003(03)00011-8 | Crossref OK (both); Ames [FT], Simmons [T] | Moderate. Explicitly says the profile is unknown. |
| Enceladus | **Non-uniform, with fresh water on top.** An overturning circulation with "a shallow freshwater lens in the polar regions" and "significant vertical and horizontal stratification". | Freezing and melting at the ice interface | Lobo et al. 2021, 10.1038/s41561-021-00706-3 | Crossref OK; [AB] (nature.com) | Moderate (one idealized model) |
| Enceladus, Europa | **Non-monotonic.** Double-diffusive layers inhibit transport from the seafloor to the ice. | Double-diffusive staircases | Wong et al. 2022, 10.1029/2022JE007316 | Crossref OK; [AB] | Moderate (models). Soderlund et al. 2024 [FT]: layers may be transient; vertical transport "would be reduced while such layers are present". |
| Enceladus, Europa | **Top, bottom or both; direction unknown.** Bottom-heated flows are weak and rotation-dominated, so "dissipation is likely to be limited mostly to turbulent boundary layers near the seafloor and the ice–ocean interface". It is "not clear" whether libration gives modest dissipation in the ice–ocean boundary layer or strong mixing through the whole depth. | Energy budget: tides, libration, convection | Jansen et al. 2023, 10.3847/PSJ/acda95 | Crossref OK; [FT] | Strong as a statement of ignorance. This is the most direct answer to H25. |
| Enceladus | **Bulk turbulence possible.** Libration-driven elliptical instability persists in shell geometry and is extrapolated to subsurface oceans. | Libration | Lemasquerier et al. 2017, 10.1002/2017JE005340 | Crossref OK; [AB] | Weak for direction (lab and DNS; no depth profile) |
| Enceladus | **No profile given.** Libration dominates linear dissipation (~0.001 GW); inertial-wave resonances can raise it. | Tides and libration, inertial waves | Rekier et al. 2019, 10.1029/2019JE005988 | Crossref OK; [AB] | Weak for direction |
| Europa | **Uniform convective; weak salinity stratification.** The bottom is only 0.01 °C warmer than the top ("well mixed"). A thin, cold, fresh, stable surface layer exists at some latitudes. Constant κ_v = 10⁻⁴ is used because KPP does not apply. | Convection driven by heat, ice melting and salinity | Ashkenazy and Tziperman 2021, 10.1038/s41467-021-26710-0 | Crossref OK; [FT] | Moderate (one GCM) |
| Europa | **Convection throughout**, three-dimensional and more vigorous at low latitudes | Rotating thermal convection | Soderlund et al. 2014, 10.1038/ngeo2021 | Crossref OK; [AB] (nature.com) | Moderate (models) |
| Europa | **Layering possible** | Double-diffusion | Vance and Brown 2005, 10.1016/j.icarus.2005.06.005 | Crossref OK; [T] | Not used |
| All four moons | **Convection throughout; regime differs.** Enceladus and possibly Titan are in a strongly rotating regime. Europa is moderate. Ganymede spans several regimes. | Rotating convection | Soderlund 2019, 10.1029/2018GL081880 | Crossref OK; [AB] | Moderate (models). No D_z profile. |
| All four moons | **Partially stratified possible; a top "stratosphere".** Ocean models usually assume "vigorous convection that precludes the formation of a density gradient", but oceans "could be partially stratified". Melosh et al. (2004) [2nd] proposed a buoyant "stratosphere" under the ice, which turbulence may erode. Libration flows "might be located in the ocean boundary layers … or in the ocean bulk". Basal heat-flux anomalies, including from high-pressure ice layers (Ganymede, Titan), reach the top efficiently. | Review of all of the above | Soderlund et al. 2024, 10.1146/annurev-marine-040323-101355 | Crossref OK; [FT] (CC BY PDF) | Strong as a review. It says the question is open. |
| Ganymede | High-pressure ice under the ocean | Interior structure | Journaux et al. 2020, 10.1007/s11214-019-0633-7 | Crossref OK; [T] | Not used |
| Titan | Convective heat transfer models | Rotating convection | Kvorka and Čadek 2022, 10.1016/j.icarus.2021.114853; Amit et al. 2020, 10.1016/j.icarus.2019.113509 | Crossref OK; [T] | Not used (ADS 405; no Crossref abstract) |

## Per-world verdict

- **Early Earth (open ocean, no lid): rising with depth supported, below the surface mixed layer.**
  - Modern measurements give K rising from ~10⁻⁵ in the pycnocline to ~10⁻⁴–10⁻³ near rough seafloor.
  - The ratio is ≈ 14 for the global upper-1000 m versus the abyss, and ~100 over rough topography. Both sit
    inside ρ = 1–100.
  - Over smooth plains, ρ ≈ 1.
  - Two caveats:
    - The surface mixed layer has much higher D_z, so the full-depth profile is non-monotonic, not exponential.
    - All data are from the modern ocean. Carrying them over to an early ocean with different tides and
      seafloor roughness is an assumption.
- **Europa: unknown.** No observation exists.
  - At low salinity (below S\* = 3.0 under 20 km), models give a stratified top over a convective bottom. That is
    a step rise, and it is already removed by the lid test.
  - At higher salinity, models give a convective, nearly uniform interior (ρ ≈ 1).
  - Tidal or libration dissipation could concentrate at the ice interface (ρ < 1) or at the seafloor (ρ > 1).
  - So rising with depth is neither supported nor contradicted. The evidence leans towards ρ ≈ 1 for the draws
    that pass the lid test.
- **Enceladus: unknown, leaning contradicted for the smooth exponential.**
  - At S < ~20 g/kg the models put a stratified layer under the ice. This is the lid case: 0% with a lid in every
    model.
  - Without a lid, the convective high-salinity case is unstratified (ρ ≈ 1).
  - Strong mixing (κ_z ~10⁻³) at S ≤ 10 can stratify the whole ocean, and Ames et al. explicitly leave the
    profile unknown.
  - The libration and tidal power is large relative to convection, and where it is dissipated is open.
- **Ganymede: unknown.** Models give convection throughout. Basal heat-flux anomalies, from a high-pressure ice
  floor, reach the top. No model gives a D_z depth profile. Rising with depth is not supported by any source found.
- **Titan: unknown.** It is in the strongly rotating convective regime ([AB] Soderlund 2019). The Titan-specific
  model papers could not be read (titles only). No D_z profile was found.

## Is the 1–100 prior reasonable?

**Not for the moons.**

- The models in the table give ρ ≈ 1 (convection), or a step that the lid test already handles, or say plainly
  that the direction is unknown.
- None argues for a smooth rise of up to 100× with depth.
- The lid runs make this sharper. Conditioning on S ≥ S\* keeps exactly the high-salinity regime where Zeng and
  Jansen and Bire et al. give an unstratified interior. So the draws that survive the lid test are those for which
  ρ > 1 has the least support.
- Counting a low-salinity stratified top as ρ > 1 in the no-lid runs would double-count the lid.
- A sharp step is also not the smooth exponential the Eq. 20 onset was tested on. Pinning at a step was not tested.

**For early Earth, 1–100 is defensible** below the mixed layer (ratios ≈ 1 over smooth plains, ≈ 14 globally,
~100 over rough topography). Report it with the mixed-layer caveat.

**A two-sided prior is more honest for the moons.** Use ρ log-uniform on 0.01–100: equal odds of up or down, the
same span, and no claim of a direction. With this prior:

- Half the draws have ρ ≤ 1. They give the step-e result exactly (`sweep3.py:check`).
- The other half are distributed as the present ρ = 1–100 sample.
- So, by arithmetic and without a rerun, each share becomes about 0.5 × step e + 0.5 × current, within sampling
  noise of a few tenths of a point.

| World (P_pved %) | current, ρ = 1–100 | two-sided, ρ = 0.01–100 (arithmetic) |
|---|---|---|
| Early Earth | 67.0–99.9 | 70.9–87.3 |
| Ganymede | 33.5–44.1 | 31.7–37.0 |
| Titan | 25.5–35.8 | 25.5–30.7 |
| Europa, no lid / lid | 29.0–39.6 / 22.0–30.1 | 29.7–35.0 / 22.6–26.6 |
| Enceladus, no lid / lid | 19.9–28.3 / 0.0 | 22.2–26.4 / 0.0 |

Reading:

- The ordering of worlds does not change.
- The upper ends roughly halve their gain over step e.
- Earth's lower bound rises, because the rising-D_z sample lowered it (67.0 < 74.7).
- A narrower prior centred on 1 (for example 0.1–10), which fits "convective, ρ ≈ 1", would move every moon closer
  to its step-e value.

These figures are arithmetic, not runs. A rerun of `sweep3.py`/`multiwall.py` with
`rho = 10 ** rng.uniform(-2, 2)` is needed before any number enters the paper.

## Proposed replacement text (not applied)

**Assumption 7** (replaces `ija_draft.md` line 535):

> 7. *Direction of vertical mixing* (Section 2.6). D_z follows a smooth exponential in depth with bottom-to-top
> ratio ρ. For early Earth's open ocean ρ is drawn log-uniformly from 1 to 100, as measured below the surface mixed
> layer of the modern ocean (Polzin et al., 1997; Ledwell et al., 2000; Waterhouse et al., 2014); the mixed layer
> itself is more strongly mixed, and the modern profile is assumed to hold for the early ocean. No measurement
> constrains the direction in any ice-covered ocean: models give a convective interior with nearly uniform mixing
> at high salinity (Bire et al., 2022; Zeng and Jansen, 2021), a stratified layer under the ice at low salinity
> (treated by the lid test, Section 2.5, not by ρ), and tidal or libration dissipation in boundary layers at the
> ice, the seafloor or both (Jansen et al., 2023; Soderlund et al., 2024). For the moons ρ is therefore drawn
> log-uniformly from 0.01 to 100. Where ρ ≤ 1 no wall reaches the Eq. 20 onset and the result equals the
> uniform-D_z sampling with all other additions on.

**H25** (replaces `ija_draft.md` line 201):

> 25. **D_z direction.** Literature search (`journals/new_analyses_dz_direction.md`): rising with depth is
> supported for Earth below the surface mixed layer (ratio ≈ 14 globally, ~100 over rough topography; Waterhouse
> et al., 2014; Ledwell et al., 2000) and unknown for Europa, Enceladus, Ganymede and Titan, where models give a
> uniform convective interior, a stratified top (the lid case), or boundary-layer dissipation at either end
> (Jansen et al., 2023; Soderlund et al., 2024). Rerun `sweep3.py`/`multiwall.py` for the moons with ρ log-uniform
> on 0.01–100; draws with ρ ≤ 1 equal sweep2 step e. Expected by arithmetic: Ganymede 32–37%, Europa 23–27% (lid)
> and 30–35% (no lid), Enceladus 0% (lid) and 22–26% (no lid). Replace "base model" in assumption 7 by "step e".

Consequential edit (not applied): in the Table 1 row "bottom/top D_z ratio ρ (assumed)", change the moon columns
to "0.01 – 100 (direction unknown)" and keep "1 – 100" for early Earth.

## New references, in the draft's format

Ames et al. 2025, Kang et al. 2022, Ledwell et al. 1993, Munk 1966, Zeng and Jansen 2021 and 2024 are already in
the list. Only papers used as evidence above are listed.

Bire, S., Kang, W., Ramadhan, A., Campin, J.M., Marshall, J. 2022, Exploring Ocean Circulation on Icy Moons Heated From Below, Journal of Geophysical Research: Planets, 127, e2021JE007025. https://doi.org/10.1029/2021JE007025

Ferrari, R., Mashayek, A., McDougall, T.J., Nikurashin, M., Campin, J.M. 2016, Turning Ocean Mixing Upside Down, Journal of Physical Oceanography, 46, 2239–2261. https://doi.org/10.1175/JPO-D-15-0244.1

Jansen, M.F., Kang, W., Kite, E.S., Zeng, Y. 2023, Energetic Constraints on Ocean Circulations of Icy Ocean Worlds, The Planetary Science Journal, 4, 117. https://doi.org/10.3847/PSJ/acda95

Kunze, E., Firing, E., Hummon, J.M., Chereskin, T.K., Thurnherr, A.M. 2006, Global Abyssal Mixing Inferred from Lowered ADCP Shear and CTD Strain Profiles, Journal of Physical Oceanography, 36, 1553–1576. https://doi.org/10.1175/JPO2926.1

Large, W.G., McWilliams, J.C., Doney, S.C. 1994, Oceanic vertical mixing: A review and a model with a nonlocal boundary layer parameterization, Reviews of Geophysics, 32, 363–403. https://doi.org/10.1029/94RG01872

Ledwell, J.R., Montgomery, E.T., Polzin, K.L., St. Laurent, L.C., Schmitt, R.W., Toole, J.M. 2000, Evidence for enhanced mixing over rough topography in the abyssal ocean, Nature, 403, 179–182. https://doi.org/10.1038/35003164

Lobo, A.H., Thompson, A.F., Vance, S.D., Tharimena, S. 2021, A pole-to-equator ocean overturning circulation on Enceladus, Nature Geoscience, 14, 185–189. https://doi.org/10.1038/s41561-021-00706-3

Polzin, K.L., Toole, J.M., Ledwell, J.R., Schmitt, R.W. 1997, Spatial Variability of Turbulent Mixing in the Abyssal Ocean, Science, 276, 93–96. https://doi.org/10.1126/science.276.5309.93

Soderlund, K.M. 2019, Ocean Dynamics of Outer Solar System Satellites, Geophysical Research Letters, 46, 8700–8710. https://doi.org/10.1029/2018GL081880

Soderlund, K.M., Rovira-Navarro, M., Le Bars, M., Schmidt, B.E., Gerkema, T. 2024, The Physical Oceanography of Ice-Covered Moons, Annual Review of Marine Science, 16, 25–53. https://doi.org/10.1146/annurev-marine-040323-101355

Waterhouse, A.F., MacKinnon, J.A., Nash, J.D., Alford, M.H., Kunze, E., Simmons, H.L., Polzin, K.L., St. Laurent, L.C., Sun, O.M., Pinkel, R., et al. 2014, Global Patterns of Diapycnal Mixing from Measurements of the Turbulent Dissipation Rate, Journal of Physical Oceanography, 44, 1854–1872. https://doi.org/10.1175/JPO-D-13-0104.1

Optional, if the table rows that use them are cited in the paper:

Ashkenazy, Y., Tziperman, E. 2021, Dynamic Europa ocean shows transient Taylor columns and convection driven by ice melting and salinity, Nature Communications, 12, 6376. https://doi.org/10.1038/s41467-021-26710-0

Wong, T., Hansen, U., Wiesehöfer, T., McKinnon, W.B. 2022, Layering by Double-Diffusive Convection in the Subsurface Oceans of Europa and Enceladus, Journal of Geophysical Research: Planets, 127, e2022JE007316. https://doi.org/10.1029/2022JE007316

## Limits of this search

- No source gives a measured or modelled D_z depth profile for any ice-covered ocean. Every model found uses a
  constant prescribed κ_z (Jansen et al. 2023 [FT]; Ames et al. 2025 [FT]; Ashkenazy and Tziperman 2021 [FT]).
- Titan- and Ganymede-specific papers (Kvorka and Čadek 2022; Amit et al. 2020; Terra-Nova et al. 2026,
  10.1016/j.icarus.2026.117166; Journaux et al. 2020), and the tidal-dissipation papers (Hay and Matsuyama 2017,
  2019; Rovira-Navarro et al. 2019; Wilson and Kerswell 2018), were verified on Crossref only. Their abstracts
  were not reachable without ADS or ScienceDirect. One of them could hold a dissipation-versus-depth result.
- Kang et al. 2022 (Science Advances) and Zhang et al. 2024 returned 403 and were not re-read.
- The arxiv and openalex MCP servers were down, so no citation walk was done.
