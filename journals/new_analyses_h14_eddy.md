# h14_eddy: in how many sampled input sets is the eddy larger than the wall? (2026-10-21)

Script `h14_eddy.py`, log `h14_eddy.log`, predictions `h14_eddy_predictions.txt` (v1, written before any results run, sha256 `8be9ed63ef5f8742b2af2c97635de5ebf86b6279c969d30bf46a58fa013c04b6`, asserted in the script). Papers in `papers/h14_eddy/` (untracked). No existing file was edited, renamed or deleted; the sweep draws are imported (`lid_melt`, `uncertainty`, `sweep2`, `sweep3`, `multiwall`), N = 200,000, seed 0, the same draws and extra inputs as `coarsen_slow.py`. Nothing is committed. Everything below marked "not applied" is proposed wording only.

Disclosure: the script was run twice. The first run gave sections (1)-(6) of the log and a stress test that read a threshold off the f ladder. I then replaced that stress test by an exact per-draw run and added section (7). Apart from the stress text and the new section (7), the two logs are identical line for line, including every break-even f (checked by diff). The ladder-threshold stress of run 1 gave smaller losses than the exact one (for example Ganymede 3.99 against 4.70 points, Callisto at U = 1 m/s 8.96 against 10.84), so it was a lower bound. The predictions file was not touched. Before the predictions file the only numbers seen were the logged `coarsen_slow.log` values, my own paper-and-pencil estimates, and one timing run of `multiwall.multi` for Ganymede (disclosed in the docstring).

## What I asked, in plain words

Section 2.3 treats eddies as diffusion, which needs the wall width w = sqrt(D_h/k) (k = k2 c) to be larger than the eddy size l. `advect2d.py` tested that at l/w = 0.3, 0.9 and 3.6. The sweep does not draw an eddy size, so the share of sampled input sets with l > w is a function of l: share(l) = mean(w < l) over the sweep's draws. I report it at three eddy sizes taken from a scaling, at published sizes where one exists, and averaged over a log-uniform range, then ask whether E86's x4 / x10 slowdown covers those draws.

## Where the eddy sizes come from

| source | what it gives | where |
|---|---|---|
| Ashkenazy and Tziperman 2021, Nat. Commun. 12:6376, doi 10.1038/s41467-021-26710-0 (Crossref OK 2026-10-21) | Europa. Rhines-type spacing L = sqrt(2U/beta), beta = 2.6e-11 m-1 s-1, U ~ 0.02 m/s, about 40 km. 3D run: columns "about 20-50 km" wide and apart. | "Spacing between the Taylor columns" (`ashkenazy2021.txt:502`); `:160` |
| Kang et al. 2022, Sci. Adv. 8:eabm4665, doi 10.1126/sciadv.abm4665 (OK) | Enceladus. l_GM about 3 km from the Rhines scale (U ~ 1e-3 m/s, beta ~ 1e-10); kappa_GM about 0.3 m2/s. Uses sqrt(U/beta), a factor sqrt(2) below Ashkenazy's convention. | Methods, "We estimate l_GM to be around 3 km" (`kang2022_epmc.xml`) |
| Zhang, Kang and Marshall 2024, Sci. Adv. 10:eadn6857, doi 10.1126/sciadv.adn6857 (OK) | Enceladus. Eddy length "several kilometers", matched to the Rossby deformation radius. | Results, eddy paragraph (`zhang2024_epmc.xml`) |
| Soderlund 2019, GRL 46:8700-8710, doi 10.1029/2018GL081880 (OK) | Rotation rates and radii (Table 1); zonal speeds up to about 1 m/s or more, used only as the upper end U = 1 m/s. | arXiv text, `soderlund2019.txt` |
| Soderlund et al. 2024, Annu. Rev. Mar. Sci. 16:25-53, doi 10.1146/annurev-marine-040323-101355 (OK) | Jet widths outside the equator are set by the Rhines scale. | `soderlund2024.txt:510` |
| Faghmous et al. 2015, Sci. Data 2:150028, doi 10.1038/sdata.2015.28 (OK) | Present-day ocean eddies: radii 25-250 km. Not an early-Earth constraint. | `faghmous2015.txt:54` |
| NASA NSSDC satellite fact sheets (Williams), read 2026-10-21 | Orbital (= synchronous rotation) periods for all moons, today's day for Earth. | Galilean, Saturnian, Uranian, Neptunian and Pluto fact sheets |

For 13 of 15 worlds no eddy-resolving simulation exists, so l is a scaling, not a measurement: l = sqrt(2U/beta), beta = 2 Omega/R at the equator, Omega = 2 pi/P with synchronous rotation assumed, R = L/pi of the sweep, and U = 1e-3 (Kang), 0.02 (Ashkenazy) and 1 m/s (peak zonal jets, upper end only). Early Earth uses today's Omega. Equatorial beta is the largest, so these l are the smallest the scaling gives for each world. Bire et al. 2022 (doi 10.1029/2021JE007025, Crossref OK) was not read: Wiley returned 403 and MIT DSpace 429, and I did not work around either.

## Result 1: the share of draws with l > w (`h14_eddy.log:23-40`, % of all 200,000 draws)

| world | l km (U 1e-3 / 0.02 / 1) | l>w % (U 1e-3 / 0.02 / 1) | log-uniform avg | published l | tied | log line |
|---|---|---|---|---|---|---|
| Enceladus | 2.1 / 9.4 / 66.1 | 8.6 / 17.5 / 34.2 | 20.0 | 10.4 (3 km) | 0.9 | 26 |
| Europa | 8.4 / 37.8 / 267.0 | 16.9 / 28.9 / 48.2 | 31.5 | 23.4-31.4 (20-50 km) | 0.9 | 27 |
| Ganymede | 16.1 / 71.9 / 508.8 | 21.5 / 34.8 / 54.9 | 37.4 | - | 0.9 | 28 |
| Titan | 23.8 / 106.3 / 751.4 | 24.8 / 38.6 / 58.9 | 41.3 | - | 0.9 | 29 |
| Callisto | 23.5 / 105.2 / 743.7 | 24.6 / 38.7 / 58.8 | 41.3 | - | 0.9 | 30 |
| Triton | 10.5 / 46.8 / 330.6 | 18.4 / 31.0 / 50.5 | 33.6 | - | 1.0 | 31 |
| Pluto | 10.2 / 45.7 / 323.0 | 18.2 / 30.5 / 50.0 | 33.2 | - | 0.9 | 32 |
| Titania | 9.7 / 43.5 / 307.3 | 17.9 / 30.3 / 49.8 | 32.9 | - | 0.9 | 33 |
| Oberon | 11.9 / 53.1 / 375.3 | 19.4 / 32.1 / 51.9 | 34.8 | - | 0.9 | 34 |
| Umbriel | 5.8 / 25.8 / 182.6 | 14.5 / 25.8 / 44.5 | 28.4 | - | 1.0 | 35 |
| Ariel | 4.5 / 20.0 / 141.7 | 12.9 / 23.6 / 41.9 | 26.2 | - | 0.9 | 36 |
| Dione | 4.6 / 20.5 / 145.3 | 12.9 / 23.6 / 42.0 | 26.2 | - | 0.9 | 37 |
| Miranda | 2.1 / 9.6 / 67.7 | 8.7 / 17.7 / 34.4 | 20.1 | - | 0.9 | 38 |
| Mimas | 1.6 / 7.2 / 50.7 | 7.4 / 15.6 / 31.5 | 18.0 | - | 0.9 | 39 |
| Early Earth ocean | 9.3 / 41.8 / 295.4 | 5.7 / 18.6 / 44.8 | 22.7 | 13.2-42.4 (25-250 km) | 21.3 | 40 |

(The eddy sizes are `h14_eddy.log:5-19`. "tied" is the cross-check in which l = D_h/U with U log-uniform 1e-3..1 m/s per draw, so D_h and l are no longer independent.)

- At the typical eddy (U = 0.02 m/s) the share is 15.6% (Mimas) to 38.7% (Callisto). At the published sizes it is 10% (Enceladus), 23-31% (Europa) and 13-42% (early Earth). In no world is it a majority at those sizes.
- It passes half only for an eddy at the Rhines scale of a 1 m/s flow, and then only in Callisto, Ganymede, Oberon, Pluto, Titan and Triton (Pluto printed as 50.0%, Titania just under at 49.8%, Europa 48.2%). At U = 0.02 the ordering follows size and slow rotation: Titan and Callisto > Ganymede > Oberon > Triton, Pluto, Titania > Europa > Umbriel > Ariel, Dione > Enceladus, Miranda > Mimas.
- Among the draws that are P_pved at f = 1 the share is higher than among all draws in 14 of 15 worlds and equal in the early Earth (`h14_eddy.log:43-57`; equal for the early Earth, whose P_pved share is almost 100%; for example Callisto 47.8% against 38.7% at U = 0.02, Mimas 27.6% against 15.6%): a draw that can win needs k t >= 10, so a larger k and a smaller wall.
- The tied variant gives about 1% for every moon. The sweep's moon D_h is at most 100 m2/s, so l = D_h/U with U >= 1e-3 m/s is at most 100 km and mostly far below the median wall width of about 310 km (`h14_eddy.log:26-39`). In the early Earth (D_h 1e2-1e4) it is 21.3%. Whether D_h and l are independent (the main table) or tied is a modelling choice I cannot settle from the papers; it changes the moon answer from 10-39% to 1%. The column "kin" in `h14_eddy.log:23-40` gives the share of draws in which D_h/l(0.02) is a speed between 1e-3 and 1 m/s: 0% for Titan and Callisto, 3.6% for Ganymede, 10.5% for Europa, 25.8% for Enceladus, 100% for early Earth. So for the large moons a 0.02 m/s Rhines eddy cannot be the source of any D_h in the sweep's range.

## Result 2: how far past the tested ratio do those draws lie? (`h14_eddy.log:60-74`)

Share of the l > w draws (U = 0.02) with l/w above the tested 3.6: 53.5% (Mimas) to 68.8% (Titan), 37.8% for early Earth. Median l/w 4.1 (Mimas) to 8.9 (Titan), 2.6 for early Earth; 90th percentile 35.5 to 177.6 (10.4 for early Earth). At U = 1 m/s the share beyond 3.6 is 59.9-77.7%. So for the moons more than half of the draws with an eddy larger than the wall lie beyond the one ratio tested.

## Result 3: does E86's x4 / x10 cover them? (`h14_eddy.log:76-173`)

Subset: draws with l > w (U = 0.02) that are P_pved at f = 1; "lost" means P_pved at f = 1 and not at f. Sideways coarsening time x f is the `coarsen_slow.py` slowdown (U.A = 5.5/sqrt(f)); the f = 1 shares reproduce the logged ones within 0.05 points (`h14_eddy.log:21`).

| world | l>w with l/w>3.6 % | subset size (main, no lid) | lost x4 / x10 (no lid; final) | f at which 5% of the subset is lost (no lid; final) | stress p = 1.1, points of N lost (no lid; final) | stress p = 2, points lost | log lines |
|---|---|---|---|---|---|---|---|
| Enceladus | 55.1 | 15273 | 10/24; 0/0 | 200; none (lid) | 0.11; 0.0 | 0.67; 0.0 | 60, 80-81 |
| Europa | 63.7 | 28880 | 16/45; 12/32 | 200; 200 | 0.42; 0.32 | 2.36; 1.79 | 61, 88-89 |
| Ganymede | 66.8 | 37593 | 36/111; 36/111 | 100; 100 | 1.43; 1.43 | 4.70; 4.70 | 62, 96-97 |
| Titan | 68.8 | 33541 | 30/108; 30/108 | 100; 100 | 1.36; 1.36 | 4.33; 4.33 | 63, 102-103 |
| Callisto | 68.6 | 41523 | 32/95; 32/95 | 100; 100 | 1.59; 1.59 | 5.25; 5.25 | 64, 108-109 |
| Triton | 64.6 | 28848 | 13/37; 13/37 | 200; 200 | 0.43; 0.43 | 2.39; 2.39 | 65, 114-115 |
| Pluto | 64.8 | 32416 | 6/21; 6/21 | 500; 500 | 0.37; 0.37 | 2.51; 2.51 | 66, 120-121 |
| Titania | 64.1 | 42567 | 1/8; 1/8 | 1000; 1000 | 0.24; 0.24 | 2.56; 2.56 | 67, 126-127 |
| Oberon | 65.5 | 45237 | 6/12; 6/12 | 1000; 1000 | 0.29; 0.29 | 2.84; 2.84 | 68, 132-133 |
| Umbriel | 61.6 | 39542 | 3/4; 3/4 | 10000; 10000 | 0.09; 0.09 | 1.68; 1.68 | 69, 138-139 |
| Ariel | 60.5 | 33048 | 1/2; 1/2 | 10000; 10000 | 0.06; 0.06 | 1.25; 1.25 | 70, 144-145 |
| Dione | 60.2 | 27703 | 1/3; 1/2 | 10000; 10000 | 0.04; 0.03 | 1.05; 0.65 | 71, 150-151 |
| Miranda | 55.6 | 14906 | 1/1; 1/1 | 10000; 10000 | 0.0; 0.0 | 0.27; 0.08 | 72, 156-157 |
| Mimas | 53.5 | 9550 | 11/24; 5/9 | 200; 200 | 0.08; 0.03 | 0.64; 0.19 | 73, 162-163 |
| Early Earth ocean | 37.8 | 37162 | 0/0; 0/0 | 10000; 10000 | 0.0; 0.0 | 0.0; 0.0 | 74, 168-169 |

(The break-even f is the first rung of the ladder 1, 2, 3, 4, 6, 10, 14, 20, 30, 50, 100, 200, 500, 1000, 1e4, 1e6, 1e10 at which at least 5% of the subset is lost, so the true 5% point lies between that rung and the one below it. For Umbriel, Ariel, Dione, Miranda and the early Earth it lies between x1000 and x10,000. The early Earth loses 2 draws at x1000 and the whole subset at x1e10.)

Reading, in the order of the acceptance test:

1. Within the ratios E86 slowed, the bound holds. In the l > w subset at most 36 draws (0.018% of N) leave P_pved at x4 and 111 (0.056%) at x10, both Ganymede (`h14_eddy.log:96`). Prediction P6 held in all 15 worlds, both rules, with and without the lid.
2. E86 does not cover the draws in the sense that matters. Only the early Earth has at least half of its l > w draws (62.2% inside) within the tested ratio 3.6; the 14 moons have 31-47% inside (`h14_eddy.log:194-209`). The slowdown at the tested ratio, "at least four times", comes from an exploratory analysis made after the run (Section 2.3), so it is a floor, not a ceiling, and nothing is measured at higher ratios. By the decision rule written before the run, 14 of 15 worlds are "covered only if the unmeasured slowdown stays below" the break-even f in the table: x100 (Ganymede, Titan, Callisto), x200 (Enceladus, Europa, Triton, Mimas), x500 (Pluto), x1000 (Titania, Oberon), x10,000 (Umbriel, Ariel, Dione, Miranda). In each case the true 5% point is at or below the number given. Only early Earth is covered outright.
3. The stress test is not a result. If the slowdown of an l > w draw grew like (l/w)^p, run exactly with a per-draw U.A: p = 1.1 is the smallest power that gives 4 at l/w = 3.6; p = 2 is the Damkohler number (13 at 3.6). With p = 1.1 the worlds lose at most 1.59 points of N (Callisto; Ganymede 1.43, Titan 1.36), against shares of 35-44%. With p = 2 they lose up to 5.25 points (Callisto), 4.70 (Ganymede), 4.33 (Titan) and 0.0-2.84 elsewhere, with the early Earth at 0.0. These are P_pved points lost from draws that win at f = 1. At U = 1 m/s the p = 2 losses reach 10.8 points (Callisto; `h14_eddy.log:110-111`). The ordering of the worlds is not tested here.

## Check of the `coarsen_slow.py` margin claim (`h14_eddy.log:175-192`)

`coarsen_slow.py:24-25` says sideways coarsening takes at most 7% of the ocean's age in any draw, "so it must be about 14x slower before it decides a single draw". The 7% is reproduced: the smallest ratio of ocean age to coarsening time over all draws is 14.4 for Ganymede, 15.0 Titan, 17.0 Callisto, 18.4 Enceladus (1/14.4 = 0.069; E86 at `ija_draft.md:187` discloses the same diagnostic as 0.070 and notes it was not in a log, so the paper does not quote it; it is now at `h14_eddy.log:178-192`). The second half is not right. Coarsening adds to the rise time in `multiwall.py:421-422` (`t_sel = max(10/k, (L/U.A)^2/D_h, scale*down)`; `timely = rec & (scale*up + t_sel <= t)`), so a draw that already sits close to the age limit moves at a small slowdown. Over all draws the first loss appears at f = 2 in 10 of 15 worlds (23 draws in Ganymede, 3 in Enceladus), at 3 in Oberon, Ariel and Miranda, at 4 in Dione, and at 200 in the early Earth. The sizes are small: 157 draws (0.08 points) at x10 in Ganymede, matching the logged 0.1 point. The conclusion of E86 stands; the sentence does not.

## Predictions (v1)

7 of 8 held, 1 failed (`h14_eddy.log:212-221`).

| | rule | result |
|---|---|---|
| P1 | published sizes: Enceladus 5-15, Europa 15-30 / 22-38, Earth 6-20 / 30-52 | held: 10.4; 23.4 / 31.4; 13.2 / 42.4 |
| P2 | all worlds 8-45% at U 0.02, Titan and Callisto > Europa > Enceladus | held: 15.6-38.7 |
| P3 | none > 50% at U 0.02; at U 1 the >50% set contains Titan and Callisto and sits inside a named seven | held: Callisto, Ganymede, Oberon, Pluto, Titan, Triton (Pluto at 50.0) |
| P4 | tied: < 8% for every moon, 8-35% for early Earth, below the typical-l share in every world | **FAILED** on the last clause: moons max 1.0, early Earth 21.3, but early Earth's tied share (21.3) is above its typical-l share (18.6); below it in 14 of 15 |
| P5 | share among P_pved draws higher than among all in >= 10 of 15, never by 20 points | held: 15 of 15 as the script counts it (the early Earth ties at 18.6), max difference 12.0 |
| P6 | l > w draws lost <= 200 at x4 and <= 300 at x10 | held: 36 and 111 |
| P7 | <= 10% of the subset lost at x100, <= 50% at x1000 | held, thinly at x100: 10.0% (Ganymede 3743 of 37593, 9.96%); 39.9% at x1000 |
| P8 | 35-75% of the l > w draws beyond l/w = 3.6 in every world | held: 37.8-68.8 |

## Decision rule, applied

P6 holds everywhere, so E86's measured loss is small for the draws it slowed. But 14 of 15 worlds have fewer than half of their l > w draws inside the tested ratio, so by the rule written first they are covered only below a break-even slowdown of x100 to x10,000, and no measurement says what the slowdown is at l/w of 10 or 100. H14 is therefore not closed; what changes is that the question is now narrow: the share is known, and the open part is the slowdown of wall removal at l/w between about 4 and 100 in the 14 moons.

## Proposed wording (not applied)

Section 2.3, "Eddy diffusion", last sentence ("I have not shown how many sampled inputs have eddies larger than the wall width.") replaced by:

> The sweep does not draw an eddy size, so I counted the sampled input sets whose wall width √(D_h/k) is smaller than a chosen eddy size (`h14_eddy.py`, pre-registered). Taking the eddy size as the Rhines scale √(2U/β) of a 0.02 m s⁻¹ flow (Ashkenazy and Tziperman, 2021), with β = 2Ω/R and synchronous rotation, eddies are larger than the wall in 16% of the draws for Mimas, 18% for Enceladus, 29% for Europa, 39% for Callisto and 19% for the early-Earth ocean (`h14_eddy.log:26–40`). Published sizes give 10% for Enceladus (3 km; Kang et al., 2022), 23–31% for Europa (20–50 km) and 13–42% for the early-Earth ocean (present-day eddy radii of 25–250 km; Faghmous et al., 2015). The share passes half only if eddies reach the Rhines scale of a 1 m s⁻¹ flow, and then only for Callisto, Ganymede, Oberon, Pluto, Titan and Triton. In 54–69% of the draws with eddies larger than the wall (38% for the early Earth) the eddy is more than 3.6 wall widths, beyond the largest ratio I tested. These shares come from the sweep's priors and from scalings, not from ocean models, except for the three published sizes.

Section 4, eddies bullet, after "...changes no share by more than 0.1 point (`coarsen_slow.py`, `coarsen_slow.log`)": replace the sentence "in this model coarsening finishes well within the ocean's age, so a slower but finite clearing still finishes in time" with

> sideways coarsening alone takes at most 7% of the ocean's age, but it adds to the time for walls to rise, so a few draws near the age limit change already at a factor of 2 (`h14_eddy.log:175–192`). In the draws with eddies larger than the wall, at most 111 of 200,000 (0.06 points) change at a factor of 10, and 5% of them change only at a factor of about 100 (Ganymede, Titan, Callisto) to between 1,000 and 10,000 (Umbriel, Ariel, Dione, Miranda and the early Earth). The 2D test shows a factor of at least 4 at an eddy of 3.6 wall widths; I have not measured it at the larger ratios that 54–69% of these draws have (38% for the early Earth), so the factor there is not bounded by a measurement. If it grew as the square of the ratio, the large moons would lose about 4–5 points (`h14_eddy.log:97, 103, 109`); this is a stress test, not a result.

E86, last sentence ("H14: the slowdown is now bounded for the shares; the eddy-size question itself stays open.") replaced by:

> H14: the share of draws with an eddy larger than the wall is 16–39% at the Rhines scale of a 0.02 m s⁻¹ flow (`h14_eddy.py`); the x4 and x10 slowdown changes at most 111 of 200,000 of those draws, but 54–69% of them (38% for the early Earth) lie beyond the one ratio tested, so the bound is open there.

H14 (open-item list, item 14), replaced by:

> **14. [Partly resolved by E84, E86 and the eddy-share count.]** **Eddy scale separation.** Eddies exceed the wall width in 16–39% of sampled draws at a typical Rhines-scale eddy (10–42% at published sizes), and 54–69% of those draws (38% for the early Earth) are beyond the largest ratio tested (3.6). The slowdown of wall removal is measured only at that ratio. Open: the slowdown at ratios of 10 to 100, and whether D_h and eddy size are independent (10–39% of moon draws) or tied through D_h = Uℓ (about 1%).

Section 4, "Mixing is the key unknown" bullet (`ija_draft.md:789`), the sentence "Treating eddies as diffusion assumes the wall width exceeds the eddy size (Section 2.3), which is not shown." replaced by:

> Treating eddies as diffusion assumes the wall width exceeds the eddy size (Section 2.3); for a Rhines-scale eddy of a 0.02 m s⁻¹ flow this fails in 16–39% of the sampled draws, and the clearing of walls is measured only up to 3.6 wall widths (`h14_eddy.log:26–40`).

New E-entry (number to be set by the author), new text, no number replaced:

> **E<next>. Eddy size against wall width (2026-10-21; H14; `h14_eddy.py`, `h14_eddy.log`, `journals/new_analyses_h14_eddy.md`; predictions in `h14_eddy_predictions.txt`, 7 of 8 held, P4 failed on one clause).** New sentences in Section 2.3 and the Section 4 eddies and "Mixing is the key unknown" bullets; E86's closing sentence and H14 reworded. Share of the sweep's draws with an eddy larger than the wall width √(D_h/k): 15.6% (Mimas) to 38.7% (Callisto) at the Rhines scale of a 0.02 m s⁻¹ flow, 5.7–24.8% at 10⁻³ m s⁻¹, 31.5–58.9% at 1 m s⁻¹; published sizes give 10.4% (Enceladus, 3 km), 23.4–31.4% (Europa, 20–50 km) and 13.2–42.4% (early Earth, 25–250 km). Of the draws with an eddy larger than the wall, 53.5–68.8% (37.8% for the early Earth) exceed the one tested ratio of 3.6. E86's ×4 and ×10 slowdown changes at most 36 and 111 of 200,000 of those draws; 5% of them change at ×100 (Ganymede, Titan, Callisto) to between ×1,000 and ×10,000. H14 stays open. Not in the model: eddy sizes for 13 of 15 worlds are scalings with synchronous rotation assumed. The `coarsen_slow.py` statement that coarsening "must be about 14x slower before it decides a single draw" is not right (first loss at ×2 in 10 of 15 worlds; `h14_eddy.log:178–192`), with negligible effect in points.

## Caveats

- The shares are fractions of the sweep's priors, not probabilities for the real oceans. Wall width scales as the square root of the D_h/k ratio, whose range is a free input.
- Eddy size is a scaling for 13 worlds, with U uncertain by three orders of magnitude; synchronous rotation is assumed; early Earth uses today's day. Higher latitudes have smaller beta and so larger Rhines scales, so the equatorial values lean low.
- The coarsening slowdown is applied to the sideways term only, as in `coarsen_slow.py`; the stress test applies it per draw and is not a model of eddies.
- One advect2d point (l/w = 3.6, 2D, small box) is the only measurement of wall removal with eddies.

## Found (outside this unit, not fixed)

- `coarsen_slow.py:24-25`: "it must be about 14x slower before it decides a single draw" is wrong as written; losses start at f = 2 in 10 of 15 worlds (`h14_eddy.log:178-192`) because coarsening adds to the rise time (`multiwall.py:421-422`). The 7% figure is right. Negligible in points.
- `uncertainty.py:36-39` lists the Europa D_h bracket as 1e-2..1e2 and says Ashkenazy and Tziperman 2021 "use 30". I did not find 30 in the paper: it sets explicit horizontal tracer diffusion at 5 m2/s (2D) and 2 m2/s (3D) (`ashkenazy2021.txt:402-405`) and diagnoses an effective eddy mixing coefficient of about 300 m2/s with a lower bound of about 200 (`ashkenazy2021.txt:555-559`), above the top of the Europa range (`uncertainty.py:44`). A larger D_h makes the wall wider, so a low bracket overstates the l > w share; I did not rerun with it.
