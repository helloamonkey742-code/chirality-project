# Which Solar System Ocean Could Let the Weak Force Choose Life's Handedness?

**Saptarshi Ghosh**, Cupertino High School, Cupertino, CA, USA

## Abstract

The weak nuclear force makes one mirror-image form of a chiral molecule very slightly lower in energy. A large, slowly driven, self-amplifying reaction could in principle turn that bias into a one-handed outcome, but only if enough molecules decide together. We apply a stochastic model of this process, with every uncertain input varied, to fourteen Solar System worlds known or suspected to hold liquid-water oceans. Two factors set the result. Where the ice shell is thin and the ocean only mildly salty, a stably stratified layer can form beneath the ice; it traps the boundary between regions of opposite hand. Without such a layer, thinner oceans merge their layers of opposite hand faster. Among confirmed oceans, Ganymede (28% of plausible cases) and Europa (22–26%) score highest; Enceladus scores lowest (0–18%). The highest values (36–39%) belong to the Uranian moons, whose oceans are unconfirmed.

## Introduction

Life uses left-handed amino acids and right-handed sugars. Kondepudi & Nelson (1985) proposed that the parity-violating energy difference (PVED) between mirror-image molecules, about 10⁻¹⁷ kT, could select life's hand if a self-amplifying reaction were large and slow enough. Our model (Ghosh, in prep.; code below) finds that about 10³³–10³⁴ molecules must decide together, which means kilometre-scale volumes of well-mixed water. Ponds and pores are too small, but icy-world oceans are not. Stochastic models show a small bias selects the hand when it exceeds the noise (Hochberg et al. 2023). Lente (2007) concluded the PVED cannot be amplified when self-copying starts from a few chiral molecules; the reaction modelled here is open and fed, with the whole population choosing at once (an assumption). Lunine et al. (1999) discussed abiotic homochirality on Titan; we consider liquid-water oceans only. Here we ask which ocean world is most favourable.

## Method

Near the threshold where amplification switches on, any such reaction reduces to dα = (γtα − α³ + g)dt + √ε dW. Here α is the excess of one hand, g = PVED/kT is the bias, ε ≈ 1/N is molecule-counting noise, and γ sets how fast conditions sweep through the threshold. The chance that the favoured hand wins is P = Φ(√2 π^(1/4) g (kτ)^(1/4) √N), checked against brute-force simulations. Physics "wins" when this probability reaches 97.7% and the reaction finishes in the time available.

In a real ocean, mixing is finite. The decision is made in many patches at once, and the patches then coarsen; boundaries move at a rate set by eddy diffusivity times curvature (Allen & Cahn 1979). Stratification makes vertical mixing slow, so the ocean freezes into stacked horizontal layers of opposite hand, which merge only as fast as vertical mixing allows. A stably stratified layer at the top of the ocean pins the boundary entirely. Such a layer forms when the salinity is below a threshold S\*, where water at its freezing point gets lighter as it cools (Zeng & Jansen 2021). We test for such a layer with the TEOS-10 equation of state (McDougall & Barker 2011): under more than about 28 MPa of ice, water at its freezing point sinks as it cools at any salinity, so no layer can form.

For each world we draw 200,000 input sets, log-uniform over plausible ranges of bias, rate constant, concentration, horizontal and vertical mixing, time and volume. Only the body's size, ocean depth, ice thickness, age and salinity change between worlds. All worlds share Europa's chemistry and mixing ranges, an assumption, since no measurement distinguishes them. Ocean depths come from published interior models (Saur et al. 2015; Vance et al. 2018; Sohl et al. 2014; Hansen et al. 2021; Johnson et al. 2016; Castillo-Rogez et al. 2023; Beuthe et al. 2016; Lainey et al. 2024; Strom et al. 2024). Worlds without a measured salinity receive a 1–100 g/kg prior.

## Results

Figure 1 shows the winning fraction per world. Lid-free worlds fall between 24% and 39%, ordered mostly by ocean depth: a thinner ocean has fewer layers to merge. Depth uncertainty moves each value by up to about ten points (Ganymede 22–33% for 50–300 km). Five worlds can form a stratified lid in some draws, and their scores drop: Europa from 29% to 22% (lid in 24% of draws), Dione from 34% to 21% (38%), Miranda from 27% to 9% (67%), Mimas from 15% to 5% (68%), and Enceladus from 22% to 0% (100%). Enceladus forms one in every draw: its measured salinity, 5–20 g/kg (Postberg et al. 2009), lies below S\* (18–22.5 g/kg under its thin shell; 20 g/kg used). Across those choices its score stays at 0–8%.

Two inputs were then replaced by newer measurements (`robust_ranking.py`). Europa's lid result assumed 20 km of ice; Juno measured 29 ± 10 km (Levin et al. 2025). Above about 23 km no lid can form, so drawing the thickness from that measurement pins Europa in only 9% of draws and raises it to 26%. That is within Ganymede's own spread (22–33%), so the two tie. Other salinity estimates for Enceladus are higher (Kang et al. 2022): with 7–30 g/kg it scores 4–6%, and with 17–40 g/kg it scores 15–18%.

Among the three confirmed oceans, Ganymede (28%) and Europa (22–26%) are therefore the most favourable, and Enceladus (0–18%) the least under every salinity estimate. Ganymede's ~200 MPa ice shell rules out a top lid. For comparison, early Earth's open ocean gives 75%.

## Caveats

These numbers are model results, not detections, and none includes the probability that the ocean exists. The Uranian oceans are predictions whose present-day thicknesses disagree between models. Titan's ocean is disputed (Petricca et al. 2025). Not modelled: salt layering deep in large oceans (Vance et al. 2014), which would lower Ganymede to 14% if it pinned the boundary half the time; ammonia, likely present in outer-Solar-System oceans but omitted by TEOS-10; and real mixing rates, which are unmeasured for every world. Magnetized mineral surfaces can also favour one hand (Ozturk et al. 2023); if a planet's field aligns magnetite across a hemisphere, the weak force keeps its lead in the worst case only if under about 2–4 × 10⁻¹⁵ of the deciding chemistry occurs there (Europa, Enceladus). No measured prebiotic reaction yet amplifies handedness in this way, although a recent two-reaction network could (Higgs & Blackmond 2025). The testable result is a ranking, not an absolute probability: a deep, unstratified ocean like Ganymede's or a thick-shelled Europa is a better place to look for a weak-force-chosen hand than Enceladus. Silica nanoparticles suggest Enceladus' ocean mixes within years (Hsu et al. 2015), which may rule out a stratified top; if so, its salinity is high and its score near the upper end above.

Acknowledgments: an AI assistant (Claude, Anthropic) helped write code and text; the author directed the work and checked all results and references.

Code and data: https://github.com/helloamonkey742-code/chirality-project (`ocean_worlds.py`, `robust_ranking.py`, `uncertainty.py`).

![Figure 1. Percentage of 200,000 plausible input draws in which the weak-force-favoured hand wins, by world. Solid bars include the stratified-lid test; hatched bars show the value if no lid formed. Colour gives ocean status.](rnaas_fig.png)

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
