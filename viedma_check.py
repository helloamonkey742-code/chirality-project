"""Viedma check: does grinding-driven crystal deracemization behave like the model?

Viedma ripening: a slurry of left (L) and right (D) crystals is ground with glass balls while
the solution recycles material. Grinding makes fragments of the same hand (self-copying);
stirring moves material around (transport). The model says:
  (a) bias x sqrt(N) decides: a small initial excess over N crystals wins once it is far
      above the chance excess ~1/sqrt(N); with no excess the winner is random.
  (b) BOTH a self-copying step and transport are needed for one hand to take over.

This is a POST-HOC reanalysis of numbers stated in open text. Not pre-registered.
Expectations written BEFORE running, after reading the sources:
  V1 Viedma's 5% crystal excess (5.7 g : 6.3 g) in a flask of > 1 million crystals gives
     Delta = 0.05 * sqrt(1e6) = 50, far above 2, so the starting majority must always win.
     The paper says it did (5% L -> 100% L, 5% D -> 100% D). Expect held.
     Caveat: the excess is by mass, not crystal count; crystal sizes differ.
  V2 The "symmetric" 6.0000 g : 6.0000 g mix: weighing to 0.0001 g leaves a mass excess of at
     most 1e-4 / 12 ~ 8e-6, i.e. Delta <= 8e-6 * sqrt(1e6) ~ 0.01, so the winner should be
     chance. The paper says "L or D randomly" but gives NO run counts, so this can only be
     checked qualitatively. The weighing bias would start to steer (Delta = 1) only if more
     than ~1.4e10 crystals decided together; report that number.
  V3 Claim (b): stirring alone (transport, no grinding) should leave the initial state
     unchanged; stirring + glass balls should reach one hand. The paper says stirred flasks
     without balls "maintain ... (several days) the initial enantiomeric excess" and every
     flask with balls reached chiral purity (8 h for 5% excess, 24 h symmetric, 4 g balls,
     600 rpm). Expect held (qualitative).
  V4 A second, independent crystal deracemization (Strecker aminonitrile, grinding-free
     dissolution/recrystallization with racemization in solution): 73 runs, 36 D / 37 L.
     Expect two-sided binomial p > 0.05 vs 50/50, and the split rules out a steering bias
     (Delta >= 1, P = 0.84) with probability < 1e-6.

Sources (DOIs checked against api.crossref.org on 2026-10-07):
  V05: Viedma 2005, Phys. Rev. Lett. 94:065504, doi:10.1103/PhysRevLett.94.065504.
       Numbers read in the open arXiv version cond-mat/0407479 (title differs from the PRL;
       the final PRL text was not checked): "5.7000 g : 6.3000 g", "6.0000 g : 6.0000 g",
       "Every solution contains a population of more than a million of crystals",
       "repeated more than 20 times", "The handedness in this last case is L or D (randomly)".
  M19: Miyagawa, Aiba, Kawamoto, Tokunaga & Kawasaki 2019, Org. Biomol. Chem. 17:1238,
       doi:10.1039/c8ob03092h (abstract): "Among the 73 experiments that were carried out,
       d- and l-enriched isomers occurred 36 and 37 times, respectively."
Not used: Noorduin et al. 2008 JACS 130:1158 (doi:10.1021/ja7106349): open copy could not be
retrieved and the abstract states no run counts or ee values.
"""
from math import comb, erf, sqrt

out = []


def say(s=""):
    print(s)
    out.append(s)


def pmf(k, n, p):
    return comb(n, k) * p**k * (1 - p) ** (n - k)


def two_sided(k, n, p=0.5):
    pk = pmf(k, n, p)
    return min(1.0, sum(pmf(i, n, p) for i in range(n + 1) if pmf(i, n, p) <= pk * (1 + 1e-9)))


phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
P1 = phi(1.0)

n_cryst = 1e6  # lower bound quoted by Viedma

# V1: 5% excess by mass
ee5 = (6.3 - 5.7) / 12.0
d1 = ee5 * sqrt(n_cryst)
say(f"V1  5% excess: ee={ee5:.3f}, Delta >= {d1:.0f} (N >= 1e6). Reported: majority always won.")
v1_reported = True  # "initial 5% L-CEE show 100% L-CEE, ... 5% D-CEE show 100% D-CEE"

# V2: symmetric mix, weighing limit
ee_w = 1e-4 / 12.0
d2 = ee_w * sqrt(n_cryst)
n_steer = 1 / ee_w**2
say(f"\nV2  symmetric: weighing ee <= {ee_w:.1e}, Delta <= {d2:.3f} at N=1e6; "
    f"would steer only if N > {n_steer:.1e}. Reported: random hand (no counts given).")

# V3: transport alone vs transport + grinding (qualitative, from text)
stir_only_changed = False  # "maintain indefinably (several days) the initial ... excess"
balls_reached_purity = True  # "all ... systems where glass balls ... intervene ... chiral purity"
say(f"\nV3  stirring only -> changed: {stir_only_changed}; stirring + balls -> one hand: "
    f"{balls_reached_purity}")

# V4: Strecker run counts
d, l, n = 36, 37, 73
p4 = two_sided(d, n)
p4_bias = max(sum(pmf(i, n, P1) for i in range(0, l + 1)),   # L favoured, saw <= 37 L
              sum(pmf(i, n, P1) for i in range(0, d + 1)))   # D favoured, saw <= 36 D
say(f"\nV4  Strecker: D={d}, L={l} of {n}: p(50/50)={p4:.2f}; "
    f"P(split | Delta=1 bias)={p4_bias:.1e}")

checks = {
    "V1": d1 >= 2 and v1_reported,
    "V2": d2 < 0.1,
    "V3": (not stir_only_changed) and balls_reached_purity,
    "V4a": p4 > 0.05,
    "V4b": p4_bias < 1e-6,
}
say("\nexpectations: " + ", ".join(f"{k} {'held' if v else 'MISSED'}" for k, v in checks.items()))
say("V2 and V3 are qualitative: the paper gives no run-by-run counts.")
with open(__file__.replace(".py", ".log"), "w") as f:
    f.write("\n".join(out) + "\n")

assert checks["V1"], "V1: a 5% excess over >1e6 crystals should always win"
assert checks["V2"], "V2: weighing bias should be far below chance at N=1e6"
assert checks["V3"], "V3: transport alone should not deracemize; with grinding it should"
assert checks["V4a"], "V4a: Strecker runs should fit 50/50"
assert checks["V4b"], "V4b: Strecker split should rule out a steering bias"
