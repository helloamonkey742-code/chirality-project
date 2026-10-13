# Two-sided D_z ratio: multiwall Section 4 rerun with ρ on 0.01–100 (2026-10-12)

Script `rho2side.py`, log `rho2side.log`. The predictions were committed before the run (0dbfb4e). Results are
appended to the docstring: 5 of 9 predictions held.

## Plain summary

- The model used to assume vertical mixing D_z always gets stronger with depth: the bottom/top ratio ρ was drawn
  from 1 to 100. `new_analyses_dz_direction.md` found no observation that supports this for the icy moons.
- This run draws ρ from 0.01 to 100 instead, so it is equally likely to get weaker or stronger with depth.
  Everything else is unchanged. `multiwall.py` and `sweep3.py` are imported, not edited.
- When mixing gets weaker with depth (ρ < 1), nothing pins the walls. Both forces push them downward, so the top
  layer's hand fills the ocean. That is exactly the old sweep2 "step e" result, and the script checks this for
  every draw.
- So every world moves about halfway from its old value toward its step-e value. Ganymede drops from 44.1 to
  37.1% (main) and Europa with the lid from 30.1 to 26.6%. Early Earth drops from 99.9 to 87.3%, which is why
  Earth should keep the one-sided prior (its rise with depth is measured).
- The ranking of worlds does not change, and seeds 0–3 agree to within 0.3 points.
- Side finding: Enceladus's extra random numbers were never independent of its own model parameters. This bug
  was already present in sweep2, sweep3 and multiwall. Fixing it moves the Enceladus no-lid numbers by up to
  about 2 points.

## Sign and the choice for ρ < 1

- `G = ln ρ / (2H)`, and ln D(floor) − ln D(top) = ln ρ. So ln ρ < 0 means D_z falls with depth.
- The sharp-interface wall speed (upward) is D(G − 2/r):
  - the gradient term pulls a wall toward low D, which means low tension σ ~ √D;
  - curvature always pulls it in.
- When ρ < 1, G < 0, so both terms push the wall down. There is no onset, no pinning and no reversal.
  - Eq. 20 must not be read with an absolute value. |d ln σ/dr| < 2/r would wrongly allow walls to "pin toward
    the top".
- The choice is that ρ < 1 draws get sweep2 step e: sink time H R/(2 D_z), with D_z the profile mean. This is the
  convention ρ > 1 draws below onset already use, so the treatment is symmetric.
  - It is asserted bit for bit for all 15 worlds, at seeds 0–3, for main and ×10.
  - About 50% of draws per world have ρ < 1.

## Table (P_pved %, seed 0, N = 200000; "no lid / final" where the lid matters)

Mixture = 0.5 × step e + 0.5 × one-sided. The seed range is the two-sided value over seeds 0–3.

| World | step e | one-sided main | **two-sided main** | seeds 0–3 | one-sided ×10 | **two-sided ×10** | seeds 0–3 |
|---|---|---|---|---|---|---|---|
| Early Earth | 74.7 | 99.9 | **87.3** | 87.2–87.4 | 99.7 | **87.2** | 87.1–87.3 |
| Umbriel | 42.3 | 53.7 | **47.9** | 47.7–47.9 | 43.6 | **43.0** | 42.7–43.0 |
| Oberon | 39.6 | 51.4 | **45.5** | 45.5–45.6 | 40.8 | **40.2** | 40.2–40.5 |
| Titania | 39.5 | 51.3 | **45.4** | 45.3–45.5 | 40.9 | **40.2** | 40.1–40.3 |
| Ariel | 40.2 | 49.6 | **44.9** | 44.8–45.0 | 39.7 | **39.9** | 39.8–40.0 |
| Ganymede | 29.9 | 44.1 | **37.1** | 37.1–37.1 | 33.5 | **31.7** | 31.7–31.9 |
| Callisto | 30.0 | 43.5 | **36.7** | 36.7–36.9 | 33.0 | **31.4** | 31.4–31.7 |
| Pluto | 32.3 | 41.1 | **36.7** | 36.6–36.8 | 30.6 | **31.5** | 31.2–31.5 |
| Triton | 29.5 | 37.0 | **33.2** | 33.2–33.3 | 26.5 | **28.0** | 28.0–28.1 |
| Titan | 25.5 | 35.8 | **30.7** | 30.5–30.7 | 25.5 | **25.5** | 25.4–25.5 |
| Europa | 30.4 / 23.1 | 39.6 / 30.1 | **35.0 / 26.6** | 35.0–35.2 / 26.6–26.8 | 29.0 / 22.0 | **29.8 / 22.6** | 29.7–29.8 / 22.6–22.7 |
| Dione | 36.5 / 22.8 | 42.0 / 26.2 | **39.3 / 24.6** | 39.2–39.3 / 24.4–24.7 | 32.2 / 20.1 | **34.4 / 21.5** | 34.3–34.5 / 21.3–21.7 |
| Miranda | 30.1 / 10.0 | 27.1 / 9.0 | **28.6 / 9.5** | 28.6 / 9.3–9.5 | 21.3 / 7.0 | **25.8 / 8.5** | 25.7–25.8 / 8.4–8.6 |
| Mimas | 16.9 / 5.4 | 17.3 / 5.6 | **17.1 / 5.5** | 17.1–17.2 / 5.4–5.6 | 11.6 / 3.7 | **14.2 / 4.5** | 14.2–14.3 / 4.5–4.7 |
| Enceladus | 24.5 / 0.0 | 28.3 / 0.0 | **26.5 / 0.0** | 26.3–26.5 / 0.0 | 19.9 / 0.0 | **20.8 / 0.0** | 20.6–20.8 / 0.0 |

- Every cell is within 0.1 point of the mixture except Enceladus ×10 no lid: 20.8 against 22.2. The seed
  collision below explains it.
- Full three-way splits (P_pved / P_chance / P_none) are in `rho2side.log`.
- The highest P_chance is Umbriel at 10.3, down from 14.9 one-sided.

## Enceladus seed collision (post-hoc, not pre-registered)

- `sweep2.chains` gives each world i its extra draws (eta, r, ρ) from `default_rng([seed, i])`.
- numpy's SeedSequence ignores trailing zeros, so for i = 0 that is `default_rng(seed)`. That is the same stream
  that drew Enceladus's own chain.
- As a result, Enceladus's eta, r and ρ are the same quantiles as its bias g, rate k2 and concentration. This was
  verified True for seeds 0–3.
- The bug predates this work and affects every Enceladus number from sweep2, sweep3 and multiwall.
- Rerun with an independent stream `[seed, 0, 1]` (`python rho2side.py posthoc`), Enceladus no lid, seed 0, P_pved %:

| | step e | step f | one-sided main | two-sided main | one-sided ×10 | two-sided ×10 |
|---|---|---|---|---|---|---|
| as logged (collided) | 24.5 | 3.9 | 28.3 | 26.5 | 19.9 | 20.8 |
| independent | 24.7 | 6.7 | 26.4 | 25.6 | 19.4 | 22.1 |

- With the independent stream, the mixture holds:
  - ×10: (24.7 + 19.4)/2 = 22.05, run 22.1;
  - main: (24.7 + 26.4)/2 = 25.55, run 25.6.
- Seeds 1–3 agree within 0.2.
- Enceladus with the lid stays 0.0 either way.
- Suggested fix: seed the extra draws with `[seed, i, 1]` or `SeedSequence(seed).spawn`. Every world except
  Enceladus would then need a regression check, because their streams would change too unless the fix is limited
  to i = 0.

## Sensitivity: real sink time for ρ < 1 (one-sided caveat)

- The step-e convention takes a wall's sink time as H R/(2 D_z). With D_z falling with depth, the true downward
  speed D(2/r + |G|) is faster.
- Integrating it through the real D(r) gives a sink time no longer than the convention in 100% of ρ < 1 draws,
  except Miranda (96.0%). It never lowers P_pved.
- Changes in main P_pved, no lid:
  - Early Earth +12.7, to 100.0;
  - Ganymede, Callisto, Umbriel, Oberon and Titania about +5;
  - Ariel, Titan, Europa, Pluto and Triton +3.2 to +4.2;
  - Dione, Mimas, Miranda and Enceladus +0.8 to +2.8.
- This overstates the case for the real profile. The same over-conservative convention is used for step e itself
  and for ρ > 1 draws below onset, and those were not recomputed. Report it as an upper side, not a correction.

## Proposed draft sentences (not applied)

**H25** (`ija_draft.md` line 201), replacing the "Expected by arithmetic" sentence:

> Rerun (`rho2side.py`, seeds 0–3, pre-registered): with ρ log-uniform on 0.01–100 for the moons, draws with
> ρ < 1 are never pinned and equal sweep2 step e exactly, and each share is the average of step e and the one-sided
> value to within 0.1 point. P_pved becomes Ganymede 31.7–37.1%, Titan 25.5–30.7%, Europa 29.8–35.0% (no lid) and
> 22.6–26.6% (lid), Enceladus 0% (lid) and 20.8–26.5% (no lid; 22.1–25.6% with independent random streams, see
> below). The ordering of worlds is unchanged. Early Earth keeps ρ on 1–100 (measured); two-sided it would be
> 87.2–87.3%.

**Assumption 7** (line 535, the proposed text in `new_analyses_dz_direction.md`), replacing its last sentence:

> Where ρ < 1 the D_z gradient and the wall curvature both drive walls downward, so no wall pins and the result
> equals the uniform-D_z sampling with all other additions on (sweep2 step e).

**New hole** (Enceladus streams):

> **Enceladus random streams.** The extra draws (η, r, ρ) for world index 0 used `default_rng([seed, 0])`, which
> numpy seeds identically to `default_rng(seed)`, the stream of Enceladus's own chain. With an independent
> stream, Enceladus no-lid P_pved moves by up to 2 points (one-sided main 28.3 → 26.4, step f 3.9 → 6.7); the lid
> result (0%) is unaffected. Rerun sweep2, sweep3 and multiwall for Enceladus with independent streams before
> quoting no-lid Enceladus numbers.
