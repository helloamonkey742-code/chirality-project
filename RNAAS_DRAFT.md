# Which Solar System Ocean Could Let the Weak Force Choose Life's Handedness?

**Saptarshi Ghosh**, Cupertino High School, Cupertino, CA, USA

## Abstract

The weak nuclear force makes one mirror-image form of a chiral molecule slightly lower in energy. A large, slow, self-amplifying reaction can turn this small bias into a one-handed result, but only if enough molecules decide together. I apply a stochastic model of this process to fourteen Solar System worlds that have, or may have, liquid-water oceans. I vary every uncertain input. Two factors set the result. First, under a thin ice shell, a mildly salty ocean can form a stable layer below the ice. This layer traps the boundary between regions of opposite hand. Second, without this layer, a thinner ocean merges its layers of opposite hand faster. Of the confirmed oceans, Ganymede (28% of plausible cases) and Europa (22–26%) score highest. Enceladus scores lowest (0–18%). The Uranian moons score highest of all (36–39%), but their oceans are not confirmed.

## Introduction

Life uses left-handed amino acids and right-handed sugars. The parity-violating energy difference (PVED) between mirror-image molecules is about 10⁻¹⁷ kT. Kondepudi & Nelson (1985) proposed that it can select life's hand in a large, slow, self-amplifying reaction. My model (Ghosh, in prep.; code below) finds that about 10³³–10³⁴ molecules must decide together. That needs kilometres of well-mixed water: ponds and pores are too small, but icy-world oceans are not.

Stochastic models show that a small bias selects the hand when it exceeds the noise (Hochberg et al. 2023). Lente (2007) found that the PVED cannot be amplified when self-copying starts from a few chiral molecules. My model instead assumes an open, fed reaction in which the whole population chooses at once. Lunine et al. (1999) discussed abiotic homochirality on Titan. I consider liquid-water oceans only. I ask which ocean world is the most favourable.

## Method

Near the threshold where amplification starts, any such reaction reduces to one equation: dα = (γtα − α³ + g)dt + √ε dW. Here α is the excess of one hand, g = PVED/kT is the bias, ε ≈ 1/N is the noise from counting N molecules, and γ sets how fast conditions sweep through the threshold. The probability that the favoured hand wins is P = Φ(√2 π^(1/4) g (kτ)^(1/4) √N), where k is the reaction rate and τ is the time available. Brute-force simulations confirm it. The weak force "wins" when two conditions are true: P is 97.7% or more, and the reaction finishes in the time available.

Real oceans mix slowly, so the decision occurs in many patches at once. The patches then coarsen: their boundaries move at a rate set by eddy diffusivity times curvature (Allen & Cahn 1979). Stratification slows vertical mixing, so the ocean freezes into stacked layers of opposite hand. These layers merge only as fast as vertical mixing allows.

A stably stratified layer at the top of the ocean (a "lid") stops the boundary completely. A lid forms when the salinity is below a threshold S\*. Below S\*, water at its freezing point becomes lighter as it cools (Zeng & Jansen 2021). I test for a lid with the TEOS-10 equation of state (McDougall & Barker 2011). Under more than about 28 MPa of ice, water at its freezing point sinks as it cools at any salinity, so no lid can form.

For each world, I draw 200,000 input sets. Each input is log-uniform over its plausible range: bias, rate constant, concentration, horizontal mixing, vertical mixing, time and volume. Only five properties change between worlds: body size, ocean depth, ice thickness, age and salinity. All worlds use Europa's chemistry and mixing ranges. This is an assumption, because no measurement tells them apart. Ocean depths come from published interior models (Saur et al. 2015; Vance et al. 2018; Sohl et al. 2014; Hansen et al. 2021; Johnson et al. 2016; Castillo-Rogez et al. 2023; Beuthe et al. 2016; Lainey et al. 2024; Strom et al. 2024). Worlds with no measured salinity get a prior of 1–100 g/kg.

## Results

Figure 1 shows the winning fraction for each world. Worlds without a lid score between 24% and 39%. Ocean depth sets most of this order, because a thinner ocean has fewer layers to merge. Depth uncertainty moves each value by up to about ten points (Ganymede: 22–33% for depths of 50–300 km).

Five worlds form a lid in some draws. Their scores drop (share of draws with a lid in brackets):

- Europa: 29% to 22% (24%)
- Dione: 34% to 21% (38%)
- Miranda: 27% to 9% (67%)
- Mimas: 15% to 5% (68%)
- Enceladus: 22% to 0% (100%)

The measured salinity of Enceladus, 5–20 g/kg (Postberg et al. 2009), is below S\* (18–22.5 g/kg under its thin shell; I use 20 g/kg). Across these choices, its score stays at 0–8%.

I then replaced two inputs with newer measurements (`robust_ranking.py`):

- **Europa.** The lid result assumed 20 km of ice; Juno measured 29 ± 10 km (Levin et al. 2025). No lid can form under more than about 23 km, so with the Juno value Europa has a lid in only 9% of draws and scores 26%. This is within Ganymede's own spread (22–33%), so the two tie.
- **Enceladus.** Higher salinity estimates (Kang et al. 2022) give 4–6% for 7–30 g/kg and 15–18% for 17–40 g/kg.

Of the three confirmed oceans, Ganymede (28%) and Europa (22–26%) are the most favourable. Enceladus (0–18%) is the least favourable under every salinity estimate. Ganymede's ice shell exerts about 200 MPa, which rules out a top lid. Early Earth's open ocean scores 75% for comparison.

## Caveats

These numbers are model results, not detections. None includes the probability that the ocean exists. The Uranian oceans are predictions, and models disagree on their present-day thickness. Titan's ocean is disputed (Petricca et al. 2025).

The model does not include:

- Salt layers deep in large oceans (Vance et al. 2014). If these layers pinned the boundary half of the time, Ganymede would drop to 14%.
- Ammonia, which is likely in outer-Solar-System oceans but missing from TEOS-10.
- Real mixing rates. They are unmeasured for every world.

Magnetized mineral surfaces can also favour one hand (Ozturk et al. 2023). Suppose a planet's magnetic field aligns magnetite across a hemisphere. In the worst case, the weak force keeps its lead only if less than about 2–4 × 10⁻¹⁵ of the deciding chemistry occurs there (Europa, Enceladus).

No measured prebiotic reaction does this yet, although a recent two-reaction network could (Higgs & Blackmond 2025).

Silica nanoparticles suggest that the Enceladus ocean mixes within years (Hsu et al. 2015), which may rule out a stratified top. If so, its salinity is high and it scores near the upper end above.

The testable result is a ranking, not an absolute probability. A deep, unstratified ocean (Ganymede) or a thick-shelled Europa is a better place than Enceladus to look for a hand chosen by the weak force.

Acknowledgments: An AI assistant (Claude, Anthropic) helped write the code and text. I directed the work and checked all results and references.

Code and data: https://github.com/helloamonkey742-code/chirality-project (`ocean_worlds.py`, `robust_ranking.py`, `uncertainty.py`).

![Figure 1. Percentage of 200,000 plausible input draws in which the hand favoured by the weak force wins, for each world. Solid bars include the lid test. Hatched bars show the value if no lid formed. Colour shows the ocean status.](rnaas_fig.png)

## References

- Allen, S. M., & Cahn, J. W. 1979, AcMet, 27, 1085, doi:10.1016/0001-6160(79)90196-2
- Beuthe, M., Rivoldini, A., & Trinh, A. 2016, GeoRL, 43, 10088, doi:10.1002/2016GL070650
- Castillo-Rogez, J. C., et al. 2023, JGRE, 128, e2022JE007432, doi:10.1029/2022JE007432
- Hansen, C. J., et al. 2021, PSJ, 2, 137, doi:10.3847/PSJ/abffd2
- Higgs, P. G., & Blackmond, D. G. 2025, PNAS, 122, e2423683122, doi:10.1073/pnas.2423683122
- Hochberg, D., Buhse, T., Micheau, J.-C., & Ribó, J. M. 2023, PCCP, 25, 31583, doi:10.1039/d3cp03311b
- Hsu, H.-W., et al. 2015, Natur, 519, 207, doi:10.1038/nature14262
- Johnson, B. C., et al. 2016, GeoRL, 43, 10068, doi:10.1002/2016GL070694
- Kang, W., Mittal, T., Bire, S., Campin, J.-M., & Marshall, J. 2022, SciA, 8, eabm4665, doi:10.1126/sciadv.abm4665
- Kondepudi, D. K., & Nelson, G. W. 1985, Natur, 314, 438, doi:10.1038/314438a0
- Lainey, V., et al. 2024, Natur, 626, 280, doi:10.1038/s41586-023-06975-9
- Lente, G. 2007, PCCP, 9, 6134, doi:10.1039/b711546f
- Levin, S. M., et al. 2025, NatAs, 10, 84, doi:10.1038/s41550-025-02718-0
- Lunine, J. I., Beauchamp, J., Smith, M. A., & Nikolaev, E. 1999, in Advances in BioChirality, ed. G. Pályi, C. Zucchi, & L. Caglioti (Amsterdam: Elsevier), 257, doi:10.1016/B978-008043404-9/50017-6
- McDougall, T. J., & Barker, P. M. 2011, Getting Started with TEOS-10 and the GSW Oceanographic Toolbox, SCOR/IAPSO WG127
- Ozturk, S. F., et al. 2023, NatCo, 14, 6351, doi:10.1038/s41467-023-42130-8
- Petricca, F., et al. 2025, Natur, 648, 556, doi:10.1038/s41586-025-09818-x
- Postberg, F., et al. 2009, Natur, 459, 1098, doi:10.1038/nature08046
- Saur, J., et al. 2015, JGRA, 120, 1715, doi:10.1002/2014JA020778
- Sohl, F., et al. 2014, JGRE, 119, 1013, doi:10.1002/2013JE004512
- Strom, C., et al. 2024, PSJ, 5, 226, doi:10.3847/PSJ/ad77d7
- Vance, S., Bouffard, M., Choukroun, M., & Sotin, C. 2014, P&SS, 96, 62, doi:10.1016/j.pss.2014.03.011
- Vance, S. D., et al. 2018, JGRE, 123, 180, doi:10.1002/2017JE005341
- Zeng, Y., & Jansen, M. F. 2021, PSJ, 2, 151, doi:10.3847/PSJ/ac1114
