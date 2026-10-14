# New analysis: does a pure pool stay pure? (persist.py, 2026-10-13)

Script `persist.py`, output `persist.log` (line numbers below refer to it). Predictions P1-P14 were written in the script's docstring and committed before the run ("persist.py: predictions before run"); they were not changed. 13 held, P11 held only in part (below).

## Plain summary

My main model stops when one hand fills the mixed ocean volume. Avetisov (1999) says that is not the end: a pool made pure would not stay pure once biochemistry that prefers one hand arises, unless the catalysts are about as selective as enzymes. I tested three things with small deterministic rate equations that reuse the draft's constants (k = ki = kd = f = 1, k0 = 1e-3, a0 = 1.6; Bada and Schroeder racemization rates; Table 1 ranges). No new parameter source was added.

1. With no autocatalysis, nothing holds the excess. At 273 K alanine loses half of it in 1.2 million years and 99% of it in 8.0 million years (aspartic acid: 0.40 and 2.7 million years). That is shorter than most of the ocean ages in Table 1. So Avetisov is right about the pool with no restoring process, and the draft's argument for a pure hand already depends on autocatalysis being active.
2. With autocatalysis switched on (open, fed), the pure state is a stable attractor and the mixed state is unstable, but the strength needed is much larger than the draft's test "k > 2 k_rac". At the repo constants the pure state exists only if racemization is below 0.087 of the decay rate kd, which means the copying rate k times the feed must be about 9 times 2 k_rac (35 times for the excess to stay above 90%).
3. A new consumer that uses the right hand more than the wrong hand does not remove the pure state if it works only when the pool is pure (the gated version). Without the gate it can: above a critical selectivity the pool flips to the other hand, it does not just become mixed.

## Numbers

Q1, racemization only (closed pool, excess exp(-2 k_rac t)):
- 273 K, alanine: half in 1.203 Myr, below 1% in 7.992 Myr (log 2). Aspartic acid: 0.402 Myr and 2.673 Myr (5).
- Alanine 308 K: 11.4 kyr and 75.9 kyr; 343 K: 281 yr and 1867 yr (3-4). Aspartic acid 343 K: 20.4 yr and 135 yr (7).
- Against the Table 1 ocean ages (9-14): at 273 K alanine the time to 1% is shorter than the age for 30.1% of the log-uniform Enceladus range and 14.2% for aspartic acid; for the Europa and Early Earth ranges the age is longer in all cases (ages are 12.5-501 times and 1.25-62.6 times the time to 1%). At 308 K and 343 K every age is at least 13 times longer.
- Avetisov's Eq. 1 with only the weak chiral field (16-20): the steady excess is at most 4.3e-6, so the field alone cannot hold purity anywhere in Table 1.

Q2, autocatalysis on (22-65):
- The pure state exists for r below r_c = 0.08698 kd, is stable there (all eigenvalues negative, 24-33), and the mixed state is unstable; above r_c only the mixed state is stable (34-36, 40). The threshold is exactly k a = kd + 2r, the pitchfork where the mixed state loses stability (38).
- Purity at r = 0.01 is 0.957; it stays above 0.9 only for r up to 0.02307 and above 0.99 only for r up to 0.00127 (28, 37). The ODE and the closed formula agree to 6 digits.
- Any starting excess down to 1e-6, of either sign, ends on the pure state of its own sign for r = 0.01 (42-48); an exactly mixed start stays mixed (49).
- Minimum strengths at r = 0.01 (55-57): copying k at least 0.682 for the pure state to exist, 0.790 for purity above 0.9; mutual inhibition ki at least 0.083 and 0.438; with ki = 0 no pure state exists at any r > 0 (54, 58). At r = 0: ki 0.011 and 0.059 (52).
- In the draft's own terms (60): k a0 must exceed 9.2 times 2 r for the pure state to exist, and 34.7 times 2 r for purity above 0.9, against the draft's factor of 1. Share of Table 1 draws that pass: moons 81.0%, 71.5%, 64.9% for factors 1, 9.2, 34.7; early Earth 100%, 99.9%, 98.5% (61-63). These shares use alanine at 273 K, the slowest racemizer, so they are the most favourable case.
- Catalyst selectivity V (Avetisov's scheme 7, a copying step that makes the wrong hand a fraction (1-V)/2 of the time): V must be at least 0.852 for the pure state to exist, 0.956 for purity above 0.9 and 0.9975 for purity above 0.99 at r = 0 (53). At r = 0.01 the numbers are 0.869, 0.975, and purity 0.99 is unreachable (57).

Q3, a consumer with selectivity S (right to wrong ratio), switched on at the pure state:
- Autocatalysis off (closed pool, r = 0.01, c = 0.1) (69-72): the excess is always lost. Time to 1% relative to racemization alone is exactly 1.0000 at S = 1, falls with S, and reaches 0.9947 (gated N = 40) or 0.0665 (always on) at S = infinity. No critical S exists because nothing restores the excess.
- Autocatalysis on, always-on consumer (81-86): critical S_c = 5.86 at c = 0.1 and 1.149 at c = 0.3. For c of 1 or more the excess is lost even at S = 1, as a non-selective drain is just extra loss. Above S_c the pool goes to the other hand (final -0.97 to -1.00), not to a mixed state. For c of 0.03 or less it persists at every S.
- Autocatalysis on, gated consumer (works only in proportion to ((1+eta)/2)^N, the chance an N-link chain is defect-free): N = 40 and N = 130 persist at every S up to infinity and every c from 0.01 to 100, with purity 0.67-0.96 for N = 40 (97-104) and 0.89-0.96 for N = 130 (106-113). N = 10 (a weak gate, 88-95) persists for c up to 1, falls to 0.30-0.63 at c = 3 and flips at c of 10 and 100. This is the part of P11 that failed.
- Enzyme-like selectivity: not sourced (no enzyme selectivity value is in `journals/ija_draft.md`), so no comparison is made. For reference only, the chapter's own chain criterion gives S about 39 for N = 40 and 129 for N = 130 (116).
- Mirror check: the D-pure start in a D-preferring world gives the exact mirror result (64-65, 73, 114).

## Caveats

- Deterministic rate equations only. The draft's noise argument (Kramers flips at small volume) is not included, and the consumer uses a simple form (rate c times a gate times the weight), not a model of real enzymes or polymer synthesis.
- The constants are the repo's order-one values (k = ki = kd = f = 1 in units of 1/k). The thresholds in terms of r/kd depend on a0 and k0; the translation to seconds uses k2 a from Table 1, not a fitted rate. The factors 9.2 and 34.7 are for this parameter set at a0 = 1.6, which is a single point, not a sweep.
- The Table 1 shares assume log-uniform k2 and a and the slowest racemizer (alanine, 273 K); warmer or faster-racemizing amino acids lower them.
- "Persists" means eta of the right hand at least 0.5 at t = 5000; the gated N = 10 results at c = 3 (0.30-0.63) sit near that line and are partly set by it.
- Avetisov's gate N and consumer strength c have no measured values here; they are scanned, not fitted.
- S is defined as the rate ratio of the consumer on the two hands, so gamma = (S-1)/(S+1) as in Avetisov's Eq. 14.

## What this does and does not say about Avetisov's claim

- Supported: a pure pool with no restoring process is lost on a timescale shorter than most ocean ages (Q1). Autocatalysis can restore purity, but only if copying is strong against racemization and flow, by a margin larger than the draft's test (Q2). With V (his scheme 7) the needed catalyst selectivity is 0.96-0.998 for a pool above 90-99% pure, which is enzyme-like in his sense (53, 57).
- Not supported as a general statement: a consumer that prefers the pool's majority hand does not remove a pure state that autocatalysis is keeping alive, unless the consumer works regardless of purity and is strong (c above about 0.1 of the decay rate here), in which case it does not make the pool racemic, it makes the other hand win. A consumer that needs a pure pool (gated, N = 40 or more) leaves the state pure at every S tested.
- Not tested: his conclusion that physical chiral factors cannot matter for the origin of homochirality. This is a statement about the bias g, and these runs set g = 0 or the draft's tiny values; they test the persistence step only. I did not run a case where the field supports the pure state.

## Suggested addition to the Section 4 limitation bullet (not applied)

Current bullet: "The model stops once one hand fills the mixed volume. Whether that hand then lasts is not modelled. Avetisov (1999) argued that ..."

Suggested replacement for the second and following sentences: "A separate small calculation (`persist.py`) shows that without autocatalysis a pure pool loses its excess in about 8 Myr for alanine at 273 K (99% lost), which is shorter than most ocean ages in Table 1, and that autocatalysis must be much stronger than the condition k > 2 k_rac to hold the pure state against racemization and flow (about nine times at the model's constants, and 35 times for a purity above 90%). A consumer that prefers the majority hand does not remove the pure state if it acts only on a pure pool; if it acts regardless of purity and is strong enough, it flips the pool to the other hand. This follows Avetisov's argument that persistence needs strong catalysis, but it is a deterministic test at one parameter set, with no noise and no measured enzyme selectivity."
