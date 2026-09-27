# Vertical mixing κ_z in icy-moon oceans: what the literature says

Checked 2026-09-27 (routine wave b; OpenAlex, publisher and arXiv full text; builder + separate critic).
Project range (`uncertainty.py`): D_z = 1e-10 to 1e-3 m²/s, log-uniform, log midpoint ≈ 3e-7.

**The bar that matters.** Layers of opposite hands merge only if H·R/(2·D_z) ≤ t (Part 3b, `sphere_layers.py`), so
D_z must be at least H·R/(2t). With H = V/(4πR²) (MEASURED arithmetic, 1 yr = 3.156e7 s):
- Enceladus (H = 39.9 km, R = 232 km): D_z ≳ **1.5e-7** (ocean 1 Gyr old) to **1.5e-4** (1 Myr old).
- Europa (H = 123 km, R = 1460 km): D_z ≳ **7e-7** (4 Gyr) to **3e-5** (100 Myr).

(An earlier draft of this file compared against 2e-9..2e-6 / 4e-9..2e-7. Those are the retired "heal across the
depth" thresholds from README Part 3; the critic caught this before commit.)

| source | DOI | value (m²/s) | kind | where | verified? |
|---|---|---|---|---|---|
| Zeng & Jansen 2021, PSJ 2:151 (Enceladus) | 10.3847/PSJ/ac1114 (arXiv:2101.10530v2) | κ_z ≈ 3e-10 to 3e-3 | scaling estimate from tidal/libration energy 3e-12 to 3e-5 W/m² | Sec. 2.2 | **yes, primary** (arXiv full text re-read 2026-09-27) |
| same paper, same paragraph | same | at the low end "the turbulent diffusivity is weaker than the molecular value (10⁻⁷ for thermal diffusion and 10⁻⁹ for tracer diffusion)" | physical floor | Sec. 2.2 | yes, primary |
| same paper | same | 5e-5 (their low-salinity run) | model input, chosen to resolve the stratified layer | simulations | yes (also quoted by Jansen et al. 2023) |
| Ames, Ferreira, Czaja & Masters 2025, Commun. Earth Environ. | 10.1038/s43247-025-02036-3 | "~10⁻⁷ to 10⁻³" | **not independent:** their ref. 31 is Zeng & Jansen 2021; the 1e-7 low end matches Z&J's thermal floor, not their stated 3e-10 | Results / Fig. 2 caption | yes (full text); runs use 1e-5, 1e-4, 1e-3 |
| Ames et al. 2025 | same | 1 | convective adjustment in unstable regions only, not background mixing | Methods | yes; not counted |
| Jansen, Kang, Kite & Zeng 2023, PSJ | 10.3847/psj/acda95 (arXiv:2206.00732v2) | calls 5e-5 (Z&J) "relatively large"; says Kang et al. use "an even larger eddy diffusivity" (5e-3) | quotes of other papers' inputs | discussion near Figs. 9–10 | yes (full text); tidal/libration energy "highly uncertain" (abstract) |
| Kang, Mittal, Bire, Campin & Marshall 2022, Sci. Adv. 8, eabm4665 | 10.1126/sciadv.abm4665 | 5e-3 | model input (3D run, S = 20 g/kg) | per Jansen et al. 2023 | SECONDARY (paper not opened) |
| Wong, Hansen, Wiesehöfer & McKinnon 2022, JGR Planets | 10.1029/2022je007316 | none given | dimensionless only | abstract | not convertible |
| Rovira-Navarro et al. 2019, Icarus | 10.1016/j.icarus.2018.11.010 | - | - | - | not read |

No Europa-specific κ_z value was found. Europa is treated with the Enceladus range, which is an ASSUMED transfer.

## Plain-language summary

Nobody has measured how fast water mixes up and down inside Enceladus or Europa. The one real estimate (Zeng &
Jansen 2021) works backwards from how much tidal energy might stir the ocean. That energy is uncertain by about
7 powers of ten, so the mixing estimate spans 3e-10 to 3e-3 m²/s. The project's 1e-10 to 1e-3 range matches it.
A 2025 paper quotes the range as 1e-7 to 1e-3, but it cites the same study, so it is not a second opinion.

The lowest values are not physical for our molecules. A dissolved molecule always spreads by itself at about 1e-9 m²/s
(less when cold), and Zeng & Jansen say so in the same paragraph. So the bottom ~14% of the sampled range (below 1e-9)
really means "no stirring, only molecular spreading". This does **not** change any result: all those draws already
fail the merge bar above (they need ≥ 1.5e-7), and they would fail at 1e-9 too. (Checked by reasoning, not by
rerunning `uncertainty.py`: 1e-9 is also below the depth-healing bar ≥ 2e-9 and below the layering bar D_h/333 ≥ 3e-5,
so no condition in the code flips.) It does make the README's
worst-case merge time (1.5e12 yr Enceladus, 2.8e13 yr Europa, quoted at D_z = 1e-10) about 10× too long. At the
1e-9 floor it is ~1.5e11 / ~2.8e12 yr. That is still far longer than the Solar System's age, so the conclusion holds.

Where does the literature sit relative to the bar? The published range straddles it, since the bar itself falls
inside the range. The values modellers actually use in their simulations (5e-5 to 5e-3) clear it for any ocean age.
But those values are chosen partly to make the simulations run, so they are choices, not measurements.
**So the literature does not settle the icy-moon answer. It confirms that vertical mixing is the unknown that decides it.**
