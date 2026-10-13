# New analyses: base-model results rerun through the current pipeline (2026-10-12)

Script `rerun_base.py` (predictions committed before the run), log `rerun_base.log`. Pipeline: `multiwall.multi` (sweep2 step f, sweep3 recount, deepest-hand cascade) on the seed-0 draws of `sweep2.chains`, N = 200,000, extras as `multiwall.run4`. "main" is the measured cascade time, "x10" the cascade x10 low end. Covers journals/ija_draft.md H19 (line 196), H20 (197) and H27 (204).

Items marked *unchanged* cannot be moved by depth-dependent D_z (reason given); the base-model label can be dropped for them.

## Unchanged

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Europa stratified share | 23.9% (strat_sweep.log:2) | 23.9% | `rerun_base.log:20` | salinity draw only; base label can be dropped |
| Enceladus merge time at D_z = 1e-10 | 1.5e12 yr (uncertainty.log:34) | 1.5e+12 yr; 1.5e+11 yr at the floor 1e-9 | `rerun_base.log:247` | formula; 1e-10 is below the pipeline's D_z floor |
| Europa merge time at D_z = 1e-10 | 2.8e13 yr (uncertainty.log:35) | 2.8e+13 yr; 2.8e+12 yr at the floor 1e-9 | `rerun_base.log:249` | formula; 1e-10 is below the pipeline's D_z floor |
| Bennu Q1 frozen-racemic shares | 96.9 / 73.1 / 64.2 / 43.1% (quench.log:2-3,6-7) | 96.9 / 73.3 / 64.3 / 43.3 | `rerun_base.log:283` | k2 and Bennu draws only; base label can be dropped |

## Sensitivities (3.3, Fig. 5)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Enceladus vertical D_z | 0->43% (uncertainty.log:4) | 5.8->50.7 (x10 0.7->39.2) | `rerun_base.log:25` |  |
| Enceladus concentration | 11->32% (uncertainty.log:5) | 12.3->44.3 (x10 9.0->30.9) | `rerun_base.log:26` |  |
| Enceladus time | 14->29% (uncertainty.log:6) | 20.0->36.6 (x10 12.6->27.2) | `rerun_base.log:28` |  |
| Enceladus bias g | 18->26% (uncertainty.log:7) | 22.2->34.3 (x10 15.9->23.9) | `rerun_base.log:29` |  |
| Enceladus rate k2 | 19->24% (uncertainty.log:8) | 23.5->33.0 (x10 16.7->23.1) | `rerun_base.log:31` |  |
| Enceladus horizontal D_h | 23->21% (uncertainty.log:9) | 29.9->26.6 (x10 21.8->18.0) | `rerun_base.log:33` |  |
| Enceladus volume | 22->21% (uncertainty.log:10) | 29.8->26.8 (x10 21.4->18.4) | `rerun_base.log:34` |  |
| Enceladus chem. factor F | 21->22% (uncertainty.log:11) | 27.9->28.7 (x10 19.7->20.2) | `rerun_base.log:35` |  |
| Europa vertical D_z | 0->57% (uncertainty.log:15) | 6.7->72.5 (x10 0.2->57.7) | `rerun_base.log:49` |  |
| Europa concentration | 23->34% (uncertainty.log:17) | 29.7->49.6 (x10 22.1->35.8) | `rerun_base.log:50` |  |
| Europa time | 23->34% (uncertainty.log:16) | 33.8->45.4 (x10 24.0->33.9) | `rerun_base.log:51` |  |
| Europa rate k2 | 26->31% (uncertainty.log:18) | 35.3->43.9 (x10 25.6->32.3) | `rerun_base.log:52` |  |
| Europa bias g | 26->31% (uncertainty.log:19) | 37.2->42.0 (x10 27.3->30.6) | `rerun_base.log:55` |  |
| Europa volume | 29->28% (uncertainty.log:21) | 41.3->38.0 (x10 30.6->27.3) | `rerun_base.log:56` |  |
| Europa horizontal D_h | 30->28% (uncertainty.log:20) | 40.3->39.0 (x10 30.4->27.5) | `rerun_base.log:58` |  |
| Europa chem. factor F | 28->29% (uncertainty.log:22) | 39.2->40.0 (x10 28.7->29.2) | `rerun_base.log:59` |  |
| Early vertical D_z | 52->97% (uncertainty.log:26) | 99.9->100.0 (x10 99.4->100.0) | `rerun_base.log:74` |  |
| Early time | 59->91% (uncertainty.log:27) | 99.9->100.0 (x10 99.4->99.9) | `rerun_base.log:75` |  |
| Early rate k2 | 75->75% (uncertainty.log:33) | 99.9->100.0 (x10 99.6->99.7) | `rerun_base.log:77` |  |
| Early volume | 78->71% (uncertainty.log:29) | 100.0->99.9 (x10 99.8->99.5) | `rerun_base.log:78` |  |
| Early bias g | 75->75% (uncertainty.log:31) | 99.9->99.9 (x10 99.7->99.7) | `rerun_base.log:79` |  |
| Early concentration | 75->75% (uncertainty.log:32) | 99.9->99.9 (x10 99.7->99.7) | `rerun_base.log:80` |  |
| Early horizontal D_h | 79->70% (uncertainty.log:28) | 99.9->99.9 (x10 99.8->99.5) | `rerun_base.log:82` |  |
| Early chem. factor F | 75->75% (uncertainty.log:30) | 99.9->99.9 (x10 99.7->99.7) | `rerun_base.log:83` |  |

## Breakeven

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Enceladus D_z crossing (m^2/s) | 2.7e-5 (breakeven.log:4) | 1.6e-05 (x10 7.8e-05) | `rerun_base.log:107` |  |
| Europa D_z crossing (m^2/s) | 6.1e-6 (breakeven.log:20) | 1.7e-06 (x10 1e-05) | `rerun_base.log:135` |  |
| Early D_z crossing (m^2/s) | 2.9e-6 (breakeven.log:36) | none (bins 100-100%) (x10 none (bins 98-100%)) | `rerun_base.log:163` |  |
| Enceladus top D_z decile | 66% (breakeven.log:4) | 63% (x10 62%) | `rerun_base.log:107` |  |
| Europa top D_z decile | 86% (breakeven.log:20) | 80% (x10 80%) | `rerun_base.log:135` |  |
| Enceladus highest non-D_z bin | 36% (time) (breakeven.log:6) | 52% (concentration); x10 36% (concentration) | `rerun_base.log:105` |  |
| Europa highest non-D_z bin | 38% (time) (breakeven.log:22) | 50% (time); x10 38% (time) | `rerun_base.log:136` |  |
| Earth time crossing (yr) | 1.4e7 (breakeven.log:37) | none: bins 100-100% | `rerun_base.log:164` |  |
| Earth D_z bin nearest Ledwell 1.1e-5 | 86.9% (breakeven.log:47) | 100.0% (bin 1.1e-05; x10 99.9%) | `rerun_base.log:171` |  |

## Assumed rules and bias (3.3)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| as published | 21.8 / 28.6 / 74.7 (assumed_rules.log:2) | 28.3 / 39.6 / 99.9  (x10 19.9 / 29.0 / 99.7) | `rerun_base.log:189` |  |
| share x0.1 | 19.7 / 27.3 / 74.7 (assumed_rules.log:3) | 25.3 / 37.0 / 99.8  (x10 18.1 / 27.3 / 99.6) | `rerun_base.log:190` |  |
| share x0.28 | 20.7 / 27.9 / 74.7 (assumed_rules.log:4) | 26.7 / 38.3 / 99.9  (x10 18.9 / 28.1 / 99.6) | `rerun_base.log:191` |  |
| share x2.1 | 22.4 / 29.0 / 74.7 (assumed_rules.log:5) | 29.1 / 40.3 / 99.9  (x10 20.4 / 29.4 / 99.7) | `rerun_base.log:192` |  |
| whole ocean | 24.0 / 29.6 / 74.7 (assumed_rules.log:6) | 31.8 / 42.4 / 99.9  (x10 21.9 / 30.5 / 99.7) | `rerun_base.log:193` |  |
| Delta >= 1 | 23.1 / 29.5 / 74.7 (assumed_rules.log:7) | 30.0 / 41.0 / 99.9  (x10 21.0 / 29.9 / 99.7) | `rerun_base.log:194` |  |
| Delta >= 3 | 21.0 / 28.1 / 74.7 (assumed_rules.log:8) | 27.2 / 38.7 / 99.9  (x10 19.2 / 28.4 / 99.7) | `rerun_base.log:195` |  |
| bias x10 | 25.7 / 30.7 / 74.7 (assumed_rules.log:9) | 33.1 / 43.2 / 99.9  (x10 23.1 / 31.2 / 99.7) | `rerun_base.log:196` |  |
| bias x100 | 27.8 / 31.2 / 74.7 (assumed_rules.log:10) | 35.9 / 44.3 / 99.9  (x10 24.9 / 31.8 / 99.7) | `rerun_base.log:197` |  |
| bias unlimited | 28.6 / 31.2 / 74.7 (assumed_rules.log:11) | 37.0 / 44.5 / 99.9  (x10 25.7 / 31.9 / 99.7) | `rerun_base.log:198` |  |

## Drive swings (3.3)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| no wobble | 22 / 29 / 75 (wobble_sweep.log:1) | 28.3 / 39.6 / 99.9  (x10 19.9 / 29.0 / 99.7) | `rerun_base.log:203` |  |
| mild: 0.1-1%, 1e5-1e8 yr | 21 / 28 / 75 (wobble_sweep.log:2) | 27.6 / 38.3 / 99.9  (x10 19.5 / 28.1 / 99.7) | `rerun_base.log:204` |  |
| plausible: 0.1-30%, 1e3-1e8 yr | 20 / 26 / 75 (wobble_sweep.log:3) | 25.6 / 35.8 / 99.8  (x10 18.1 / 26.4 / 99.6) | `rerun_base.log:205` |  |
| harsh: 10-30%, 1e3-1e5 yr | 17 / 24 / 75 (wobble_sweep.log:4) | 21.7 / 31.7 / 99.4  (x10 15.4 / 23.6 / 99.1) | `rerun_base.log:206` |  |

## H20 (Delta -> c Delta)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Enceladus | 28.3 / 0.0 (no lid / final) | c=0.88 27.9 / 0.0; c=0.98 28.2 / 0.0; x10 c=0.88 19.7 / 0.0 | `rerun_base.log:211` |  |
| Europa | 39.6 / 30.1 (no lid / final) | c=0.88 39.3 / 29.9; c=0.98 39.6 / 30.0; x10 c=0.88 28.8 / 21.8 | `rerun_base.log:212` |  |
| Ganymede | 44.1 / 44.1 (no lid / final) | c=0.88 43.9 / 43.9; c=0.98 44.1 / 44.1; x10 c=0.88 33.3 / 33.3 | `rerun_base.log:213` |  |
| Titan | 35.8 / 35.8 (no lid / final) | c=0.88 35.6 / 35.6; c=0.98 35.8 / 35.8; x10 c=0.88 25.4 / 25.4 | `rerun_base.log:214` |  |
| Callisto | 43.5 / 43.5 (no lid / final) | c=0.88 43.2 / 43.2; c=0.98 43.4 / 43.4; x10 c=0.88 32.8 / 32.8 | `rerun_base.log:215` |  |
| Triton | 37.0 / 37.0 (no lid / final) | c=0.88 36.7 / 36.7; c=0.98 37.0 / 37.0; x10 c=0.88 26.3 / 26.3 | `rerun_base.log:216` |  |
| Pluto | 41.1 / 41.1 (no lid / final) | c=0.88 40.7 / 40.7; c=0.98 41.0 / 41.0; x10 c=0.88 30.3 / 30.3 | `rerun_base.log:217` |  |
| Titania | 51.3 / 51.3 (no lid / final) | c=0.88 50.8 / 50.8; c=0.98 51.2 / 51.2; x10 c=0.88 40.5 / 40.5 | `rerun_base.log:218` |  |
| Oberon | 51.4 / 51.4 (no lid / final) | c=0.88 50.9 / 50.9; c=0.98 51.3 / 51.3; x10 c=0.88 40.4 / 40.4 | `rerun_base.log:219` |  |
| Umbriel | 53.7 / 53.7 (no lid / final) | c=0.88 53.1 / 53.1; c=0.98 53.6 / 53.6; x10 c=0.88 43.2 / 43.2 | `rerun_base.log:220` |  |
| Ariel | 49.6 / 49.6 (no lid / final) | c=0.88 49.0 / 49.0; c=0.98 49.5 / 49.5; x10 c=0.88 39.3 / 39.3 | `rerun_base.log:221` |  |
| Dione | 42.0 / 26.2 (no lid / final) | c=0.88 41.5 / 25.9; c=0.98 41.9 / 26.2; x10 c=0.88 31.9 / 19.9 | `rerun_base.log:222` |  |
| Miranda | 27.1 / 9.0 (no lid / final) | c=0.88 26.8 / 8.9; c=0.98 27.0 / 9.0; x10 c=0.88 21.1 / 7.0 | `rerun_base.log:223` |  |
| Mimas | 17.3 / 5.6 (no lid / final) | c=0.88 17.1 / 5.5; c=0.98 17.3 / 5.5; x10 c=0.88 11.4 / 3.7 | `rerun_base.log:224` |  |
| Early | 99.9 / 99.9 (no lid / final) | c=0.88 99.9 / 99.9; c=0.98 99.9 / 99.9; x10 c=0.88 99.7 / 99.7 | `rerun_base.log:225` |  |

## Lids

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Juno Europa (ice 29+/-10 km) | 26.1%, pinned 9.1% (robust_ranking.log:1) | 36.0% (x10 26.3%), pinned 9.1% | `rerun_base.log:230` |  |
| Kang Enceladus S 5-20 g/kg, S* 20.0 | 0.0%, pinned 100% (robust_ranking.log:3) | 0.0% (x10 0.0%), pinned 100.0% | `rerun_base.log:231` |  |
| Kang Enceladus S 5-20 g/kg, S* 22.5 | 0.0%, pinned 100% (robust_ranking.log:4) | 0.0% (x10 0.0%), pinned 100.0% | `rerun_base.log:232` |  |
| Kang Enceladus S 7-30 g/kg, S* 20.0 | 6.0%, pinned 72% (robust_ranking.log:5) | 7.9% (x10 5.6%), pinned 72.2% | `rerun_base.log:233` |  |
| Kang Enceladus S 7-30 g/kg, S* 22.5 | 4.2%, pinned 80% (robust_ranking.log:6) | 5.5% (x10 3.9%), pinned 80.3% | `rerun_base.log:234` |  |
| Kang Enceladus S 17-40 g/kg, S* 20.0 | 17.7%, pinned 19% (robust_ranking.log:7) | 22.9% (x10 16.1%), pinned 18.9% | `rerun_base.log:235` |  |
| Kang Enceladus S 17-40 g/kg, S* 22.5 | 14.7%, pinned 33% (robust_ranking.log:8) | 19.0% (x10 13.4%), pinned 32.5% | `rerun_base.log:236` |  |
| strat alt Enceladus S 5-20, S* 18.0 | 1.6%, p 92.4% (strat_sweep.log:5) | 2.2% (x10 1.5%), pinned 92.4% | `rerun_base.log:237` |  |
| strat alt Enceladus S 5-40, S* 20.0 | 7.4%, p 66.6% (strat_sweep.log:6) | 9.5% (x10 6.7%), pinned 66.5% | `rerun_base.log:238` |  |
| strat alt Enceladus S 1-100, S* 20.0 | 7.7%, p 65.0% (strat_sweep.log:7) | 9.9% (x10 6.9%), pinned 65.1% | `rerun_base.log:239` |  |
| strat alt Europa S 1-100, S* 1.5 | 26.0%, p 8.8% (strat_sweep.log:8) | 36.1% (x10 26.4%), pinned 8.8% | `rerun_base.log:240` |  |
| strat alt Europa S 1-100, S* 19.5 | 10.1%, p 64.6% (strat_sweep.log:9) | 14.0% (x10 10.3%), pinned 64.6% | `rerun_base.log:241` |  |
| strat alt Europa S 3-100, S* 3.0 | 28.5%, p 0.0% (strat_sweep.log:10) | 39.6% (x10 29.0%), pinned 0.0% | `rerun_base.log:242` |  |
| club-sandwich Ganymede (lid p = 0.5) | 14.2% (jupiter_moons.log:7) | 22.1% (x10 16.8%), pinned 50.0% | `rerun_base.log:243` |  |

## Merge (3.2)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Enceladus merge-in-time share | 33% (uncertainty.log:3) | 39.1% (raw D_z: 33.4%) | `rerun_base.log:248` |  |
| Europa merge-in-time share | 34% (uncertainty.log:14) | 39.3% (raw D_z: 33.6%) | `rerun_base.log:250` |  |

## Table 3 depth column

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Ganymede | 50 km 33% / 300 km 22% final (ocean_worlds.log:8) | 50 km 51.6% / 300 km 30.8% final (x10 41.0 / 20.9) | `rerun_base.log:254` |  |
| Titan | 50 km 33% / 400 km 20% final (ocean_worlds.log:10) | 50 km 51.7% / 400 km 26.9% final (x10 41.2 / 17.5) | `rerun_base.log:255` |  |
| Callisto | 20 km 38% / 300 km 22% final (ocean_worlds.log:7) | 20 km 60.9% / 300 km 30.6% final (x10 50.6 / 20.5) | `rerun_base.log:256` |  |
| Triton | 30 km 36% / 200 km 26% final (ocean_worlds.log:9) | 30 km 54.4% / 200 km 33.8% final (x10 43.6 / 23.7) | `rerun_base.log:257` |  |
| Pluto | 50 km 33% / 200 km 27% final (ocean_worlds.log:6) | 50 km 48.5% / 200 km 33.2% final (x10 37.7 / 23.3) | `rerun_base.log:258` |  |
| Titania | 10 km 41% / 50 km 34% final (ocean_worlds.log:5) | 10 km 61.0% / 50 km 46.4% final (x10 51.2 / 36.0) | `rerun_base.log:259` |  |
| Oberon | 10 km 41% / 50 km 34% final (ocean_worlds.log:4) | 10 km 60.7% / 50 km 46.2% final (x10 50.8 / 35.8) | `rerun_base.log:260` |  |
| Umbriel | 10 km 41% / 30 km 37% final (ocean_worlds.log:2) | 10 km 59.0% / 30 km 49.7% final (x10 49.5 / 39.6) | `rerun_base.log:261` |  |
| Ariel | 5 km 44% / 30 km 37% final (ocean_worlds.log:3) | 5 km 62.8% / 30 km 49.6% final (x10 55.2 / 39.6) | `rerun_base.log:262` |  |
| Dione | 30 km 23% / 100 km 20% final (ocean_worlds.log:12) | 30 km 30.9% / 100 km 23.5% final (x10 24.6 / 17.7) | `rerun_base.log:263` |  |
| Miranda | 100 km 9% / 150 km 8% final (ocean_worlds.log:13) | 100 km 8.9% / 150 km 8.3% final (x10 7.1 / 7.0) | `rerun_base.log:264` |  |
| Mimas | 35 km 5% / 55 km 4% final (ocean_worlds.log:14) | 35 km 6.3% / 55 km 5.0% final (x10 4.2 / 3.2) | `rerun_base.log:265` |  |

## k2 extension and Bennu (3.6)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Enceladus k2 1e-12..1: wins, R2 P(win\|Bennu racemic) | 14% wins, P(win\|rac) 2% (review.log:19) | 18.6%, 1.9% (x10 13.1%, 1.4%) | `rerun_base.log:274` |  |
| Enceladus R6 step rule, k2 1e-12..1 (s=1 / s=0.01) | 15 / 17 (review.log:72-73) | 20.8 / 23.3 (x10 14.6 / 16.3); frozen rule 21.2 / 24.3 | `rerun_base.log:274` |  |
| Europa k2 1e-12..1: wins, R2 P(win\|Bennu racemic) | 19%, 4% (review.log:20) | 26.0%, 4.6% (x10 18.8%, 3.0%) | `rerun_base.log:275` |  |
| Europa R6 step rule, k2 1e-12..1 (s=1 / s=0.01) | 21 / 22 (review.log:72-73) | 28.6 / 31.1 (x10 20.7 / 22.7); frozen rule 29.0 / 32.3 | `rerun_base.log:275` |  |
| Early k2 1e-12..1: wins, R2 P(win\|Bennu racemic) | 57%, 15% (review.log:21) | 73.6%, 13.9% (x10 73.3%, 13.8%) | `rerun_base.log:276` |  |
| Early R6 step rule, k2 1e-12..1 (s=1 / s=0.01) | 57 / 55 (review.log:72-73) | 74.3 / 72.9 (x10 74.1 / 72.7); frozen rule 75.1 / 75.6 | `rerun_base.log:276` |  |
| Enceladus R6 step rule, k2 1e-6..1 (s=1 / s=0.01) | 22 / 23 (review.log:70-71) | 28.5 / 31.4 (x10 20.1 / 22.1) | `rerun_base.log:271` |  |
| Europa R6 step rule, k2 1e-6..1 (s=1 / s=0.01) | 29 / 30 (review.log:70-71) | 39.9 / 42.7 (x10 29.2 / 31.4) | `rerun_base.log:272` |  |
| Early R6 step rule, k2 1e-6..1 (s=1 / s=0.01) | 75 / 75 (review.log:70-71) | 99.9 / 99.9 (x10 99.7 / 99.7) | `rerun_base.log:273` |  |

## Magnetite (4)

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Earth weak (L) + field: L / D / pinned | 35 / 35 / 7% (earth.log:4; 34.7 / 34.6 / 7.0 earth_breakeven.log:1) | 46.6 / 46.4 / 6.9% (x10 46.5 / 46.3) | `rerun_base.log:287` |  |
| Earth f_mag where L < 50% | 1.88e-15 (earth_breakeven.log:44) | 5.05e-06 (x10 7.83e-13) | `rerun_base.log:292` | L sits on the coin-flip plateau (within 0.5 point of 50%) from f_mag 3.2e-13 (x10 3.2e-13); the formal crossing is where it finally dips below 50 |
| Enceladus f*: one source / M=1e3 / field | 7.3e-17 / 4.1e-15 / 1.2e-15 (moon_nudges.log:19) | 2.7e-16 / 1.5e-14 / 4.0e-15 (x10 M=1e3 1.4e-14, field 3.7e-15) | `rerun_base.log:296` |  |
| Europa f*: one source / M=1e3 / field | 5.2e-17 / 3.0e-15 / 1.3e-15 (moon_nudges.log:20) | 1.3e-16 / 7.5e-15 / 3.2e-15 (x10 M=1e3 7.3e-15, field 3.3e-15) | `rerun_base.log:297` |  |
| Early f*: one source / M=1e3 / field | 3.9e-17 / 2.3e-15 / 6.0e-16 (moon_nudges.log:21) | 1.1e-16 / 6.0e-15 / 1.6e-15 (x10 M=1e3 6.0e-15, field 1.6e-15) | `rerun_base.log:298` |  |
| Ganymede f*: one source / M=1e3 / field | 4.7e-17 / 2.8e-15 / 1.3e-15 (moon_nudges.log:22) | 1.3e-16 / 7.1e-15 / 3.1e-15 (x10 M=1e3 6.9e-15, field 3.0e-15) | `rerun_base.log:299` |  |

## H27 temperature

| item | old (source) | new | log | note |
|---|---|---|---|---|
| Earth g only (T uniform) | 99.9 (273 K, this log) | 99.9 (x10 99.7) | `rerun_base.log:313` |  |
| Earth Ala k_rac, pipeline cut | 99.9 (273 K, this log) | 97.5 (x10 97.3) | `rerun_base.log:314` |  |
| Earth Ala g + k_rac, cut on every draw | 99.9 (273 K, this log) | 93.9 (x10 93.7) | `rerun_base.log:316` |  |
| Earth Asp g + k_rac, cut on every draw | 99.9 (273 K, this log) | 86.4 (x10 86.2) | `rerun_base.log:322` |  |
| Earth Ala fixed 343 K, cut on every draw | 99.9 (273 K, this log) | 79.9 (x10 79.7) | `rerun_base.log:319` |  |
| Earth Asp fixed 343 K, cut on every draw | 99.9 (273 K, this log) | 61.8 (x10 61.7) | `rerun_base.log:325` |  |
| moons | 273 K | unchanged | `rerun_base.log:313` | oceans under ice sit at the melting point; colder water only slows k_rac and raises g |

## H27 temperature sources (DOIs checked against Crossref)

- Blake, Chang & Lepland 2010, Nature 464, 1029, 10.1038/nature08952: Archean ocean 26-35 C (phosphate oxygen isotopes).
- Robert & Chaussidon 2006, Nature 443, 969, 10.1038/nature05239: ~70 C at 3.5 Ga (cherts, Si and O isotopes).
- Knauth & Lowe 2003, GSA Bulletin 115, 566, 10.1130/0016-7606(2003)115<0566:HACTIF>2.0.CO;2: 55-85 C (Onverwacht cherts).

The sampled range 273-343 K spans the cool and the hot readings.

## Caveats

- Seed collision (rerun_base.log:21): numpy treats default_rng([0, 0]) as default_rng(0), so in sweep2.run / multiwall.run4 Enceladus (world index 0) draws eta, r and rho from the same uniforms as its own bias g, k2 and concentration. Each pair is perfectly rank-correlated, which is why the Enceladus sensitivity rows come in identical pairs. Every Enceladus pipeline number, including the regression targets in multiwall.log, carries this; the other 14 worlds do not.
- Assumed-rule share rows multiply every layer share below 1 (sweep2 top layer, sweep3 deep layer, cascade recursion layer); assumed_rules.py scaled only the top-layer share.
- Old Earth values (74.7, 75) predate the cascade; the pipeline Earth sits at 99.9, so Earth rows mostly show a ceiling.

## Predictions

- P1: held
- P2: FAILED
- P3: held
- P4: held
- P5: FAILED
- P6: held
- P7: held
- P8: FAILED
- P9: held
- P10: held
- P11: held
- P12: FAILED
- P13: FAILED
- P14: FAILED
- P15: held
