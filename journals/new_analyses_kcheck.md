# New analysis: is "k > 2 k_rac" too lenient? (kcheck.py, 2026-10-13)

Script `kcheck.py`, output `kcheck.log` (line numbers below refer to it). Predictions P1-P15 were written in the script's docstring and committed before the run (2617e29); 9 held and 6 failed, all listed in the RESULTS block of the docstring. Nothing in `journals/ija_draft.md` was changed; the suggested wording below is not applied.

## Verdict

**Yes, but only in a stated regime, and for the headline shares the damage is small.**

1. As a statement about whether the network amplifies, the test is too lenient. In the draft's own network (Eq. 2, all rate constants equal, ramp ending at 1.6 times the crossing value) a chiral branch exists only if k is about **5.7 times** 2 k_rac, and it holds a purity of 0.9 only if k is about **21.7 times** 2 k_rac. At the draft's factor of 1 the network does nothing: runs with 2 r = k and 2 r = 2k end exactly mixed (log 32-33).
2. The test is close to right only if the reactor runs far past its threshold: with the repo's constants it is adequate once the feed is about six times the crossing value (a0 = 5 needs 1.18; a0 = 7 needs 0.88; log 23-24). The draft's ramp stops at 1.6 times, where the real factor is 5.7.
3. For the headline shares it matters by a few points on the icy worlds and by almost nothing on early Earth. No ranking changes and no claim in the abstract reverses.

The two sets of numbers in `persist.py` (9.2 and 34.7) and here (5.7 and 21.7) are the same thing with different bookkeeping, see "Reconciling". Which pair to quote depends on whether the draw's k2 c is the network's copying rate k at the crossing (then 5.7 and 21.7) or the rate at the end of the ramp (then 9.2 and 34.7).

## Where the test is defined and used

- Definition: `sweep2.py:103` (`K_RAC["Ala"]` = 9.13e-15 per second, so 2 k_rac = 1.826e-14), used at `sweep2.py:133` as `ok &= k2 * c > 2 * krac` under the step-d switch. It is also at `rho2side.py:120`, `sweep3.py:116` (the deep-layer sensitivity only) and in the "cut on every draw" variants of `rerun_base.py` (H27).
- The draft states it in Section 2.6 condition 1 ("the reaction finishes and outpaces racemization: k2 a tau >= 10 and k2 a > 2 k_rac").
- Its origin is `review.py` R7 (`review.log:75-83`): in the normal form the amplifier rate is at most one copying rate k, so the choice stops once 2 r > k.
- The headline pipeline (`multiwall.multi`, `sweep3.outcome3`) does **not** apply it to the reversed-wall (recursion) draws. Applying it to every draw (M = 1 below) is itself a change of 2.4 points for Ganymede.

## Do the test and persist.py answer the same question? No.

With racemization r each way, the difference of the two hands u = L - D obeys du/dt = (k a - kd - 2 r) u + g k a s. So the amplifier grows only if k a > kd + 2 r. The draft's test is this condition with kd set to zero. It is necessary, not sufficient, in two ways:

- kd is not zero. In the draft's units kd = k, and the racemic state eats the feed: its own feed is a_rac = 1.174, not 1.6 (log 3).
- The racemic state loses stability exactly when 2 r = a_rac - 1 (checked: the antisymmetric eigenvalue equals a_rac - 1 - 2 r to 5e-11 for r from 0 to 0.2, log 4-9). This is also where the pure branch ceases to exist. So r_c = (a_rac - 1)/2 = 0.08698 kd, the number `persist.py` found by integration.
- A closed form (k0 neglected) gives 2 a_rac^2 - (2 - ki) a_rac - ki a0 = 0, so r_c = 0.0894, within 3% of the numeric value (log 10).

In units of the draft's test, the ratio k / (2 k_rac) that the network needs is 1 / (2 r_c / kd) = **5.75** for existence, **21.7** for purity of 0.9 and **395** for 0.99 (log 11). The weak bound the draft states, (a0 - 1)/2 = 0.3 kd, is 3.45 times r_c (log 10).

The 5.75 depends on constants the draft does not vary (cross-inhibition ki, spontaneous production k0). Over ki from 0.3 to 3 and k0 from 1e-4 to 1e-2 it runs from 3.1 to 26.1 (log 15-17); it is largest when cross-inhibition is weak. It never falls below 1 at a0 = 1.6.

## Reconciling 9.2 and 34.7 with 5.75 and 21.7

`persist.py` printed "k a0 must exceed 9.2 times 2 r" with a0 = 1.6, so it measured k at the end of the ramp, which is 1.6 times the copying rate at the crossing, where kd = k. In the draft's identification (kd = k = k2 c) the factors are 9.2 / 1.6 = 5.75 and 34.7 / 1.6 = 21.7 (log 11-12 shows both pairs). Both are right; the second pair is the one that matches the draft's own "k2 c". The `new_analyses_persistence.md` replacement bullet that says "about nine times ... 35 times" should say "about six times ... 22 times" if it is used with Table 1's k2 c, or keep "9 and 35" and say that k is measured at the end of the ramp. I would use the first.

## Does a deterministic ramp agree? (P5)

Ramp a0 from 0.8 to 1.6, no noise (log 29-34): at 0.5 r_c it reaches the predicted purity 0.7811 exactly; at 1.05 r_c, r = 0.25 (factor 2 in the draft's test) and r = 0.5 (factor 1) it ends exactly mixed. At 0.8 r_c the run ended mixed, which the pre-registered prediction did not expect (it failed). A post hoc check (scratch, not in the script) shows why: the 1e-8 start excess decays to floating-point zero before the racemic state turns unstable (t = 1691 of 2000), so there was nothing to amplify. Re-seeding 1e-8 at the crossing gives purity 0.529589, equal to the stationary value to six digits. This says the deterministic test needs the noise the real system has; it does not weaken the threshold.

## How much of the paper's sample fails the stricter test? (B1, B5)

The extra mask is k2 c > M 2 k_rac on every draw (win mask and one-hand mask), with M the factor. Share of all 200,000 draws that pass (log 43-58), alanine rate:

- M = 1 (the draft's test): icy worlds 80.9-81.3%, early Earth 100%.
- M = 5.75 (pure branch exists): icy worlds 73.6-74.0%, early Earth 100.00%.
- M = 21.7 (purity 0.9): icy worlds 67.2-67.5%, early Earth 99.2%.
- Same two factors times 1.6 (the 9.2 and 34.7 convention): 71.4-71.7% and 64.8-65.1% for icy worlds, 99.88% and 98.49% for early Earth.

Fraction of the draws that **currently win** (multiwall main, no lid) that fail the stricter test (log 120-135):

- M = 5.75: Europa 11.9%, Ganymede 12.7%, Titan 13.5% (largest), Enceladus 3.8%, Mimas 1.5%, early Earth 0.00%.
- M = 21.7 times 1.6: Europa 20.8%, Ganymede 21.9%, Titan 22.7%, Enceladus 9.2%, Mimas 6.5%, early Earth 1.51%.

## Does any headline share change?

Multiwall main, no lid, P_pved in percent; unmasked first (log 63-77):

| world | unmasked | draft test on every draw | M = 5.75 | M = 21.7 | 21.7 x 1.6 |
|---|---|---|---|---|---|
| Early Earth | 99.9 | 99.9 | 99.9 | 99.1 | 98.4 |
| Ganymede | 44.1 | 41.8 | 38.5 | 35.6 | 34.5 |
| Europa | 39.6 | 37.6 | 34.9 | 32.4 | 31.4 |
| Enceladus | 26.4 | 26.1 | 25.4 | 24.4 | 23.9 |

The abstract's numbers (log 148-156), x10-to-main range, final (with the salt lid):

| quantity | draft now | M = 1 everywhere | M = 5.75 | M = 21.7 x 1.6 |
|---|---|---|---|---|
| Ganymede | 33.5-44.1 | 31.9-41.8 | 29.4-38.5 | 26.2-34.5 |
| Europa | 22.0-30.1 | 20.9-28.5 | 19.3-26.5 | 17.3-23.8 |
| Enceladus, no lid | 19.4-26.4 | 19.3-26.1 | 18.6-25.4 | 17.5-23.9 |
| Enceladus, with lid | 0.0 | 0.0 | 0.0 | 0.0 |
| Early Earth, sweep3 "at least" | 67.0 | 67.0 | 67.0 | 66.0 |

- The largest shift of any abstract number is 2.4 points at M = 1 on every draw, 5.6 at M = 5.75 (Ganymede) and 9.7 at M = 21.7 x 1.6 (log 149-154). The pre-registered "at most 3 at M = 5.75" failed.
- Ranking: early Earth > Ganymede > Europa > Enceladus holds at every level; Enceladus with a lid stays 0.0; the largest "chance" share over the 15 worlds falls from 14.9% to 4.4% (3.7% at the aspartic acid rate, log 113), so "chance at most 15%" is still true and gets weaker.
- Early Earth multiwall main stays 99.9 at M = 5.75 and is 98.4 at the harshest convention (97.2 with the aspartic acid rate, log 158). The 89% figure from the deepest-layer rule (`handoff.py`, `deephand_scale.py`) was not rerun here, so it is not covered by these numbers.
- Regression: with no mask the rerun reproduces the committed multiwall / sweep3 values for all 15 worlds (log 40).

## Suggested wording for the draft (not applied)

Section 2.6 condition 1, replace "k2 a > 2 k_rac" with:

> k2 a exceeds 2 k_rac. In the network of Eq. 2 this is only a necessary condition: with all rate constants equal and the feed raised to 1.6 times its crossing value, a chiral steady state exists only if k2 c is about 5.7 times 2 k_rac (and reaches a purity of 0.9 only above about 22 times), because outflow and feed depletion also count against the amplifier.

Section 4 limitation bullet, add:

> I applied the weaker test k2 c > 2 k_rac to every draw in the main runs except the reversed-wall ones. Replacing it with the network's own threshold (5.7 times) lowers the shares by 0.3-5.6 points for the icy worlds (Ganymede 44.1 to 38.5, Europa 39.6 to 34.9, Enceladus 26.4 to 25.4 without a lid) and by less than 0.1 for early Earth; the order of the worlds and the lid result do not change.

If the author wants the abstract to stay true under the stricter threshold, the three ranges would read Ganymede 29-39%, Europa 19-27%, Enceladus 19-25% (no lid), early Earth "at least 67%". The 21.7 and 395 factors need kd = k = ki = f as chosen; the 5.75 itself moves from 3 to 26 over the unsampled constants (log 15-17), so I would quote one value with that range and not the pure-branch factors alone.

## Caveats

- The map k2 c := kd is the draft's own identification; if the real outflow is slower than the copying rate at the crossing, the factor falls (M_ex = 0.88 at a0 = 7 shows the direction).
- Seed 0 only; I did not repeat other seeds.
- The masked arrays remove draws from both the win mask and the one-hand mask; they do not reclassify them as "chance", so P_none rises (log 114-117).
