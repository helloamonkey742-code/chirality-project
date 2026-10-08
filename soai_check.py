"""Soai check: do published Soai-reaction run counts behave like the model's claim (a)?

Claim (a) of the model: a bias only controls the outcome once (bias x sqrt(N)) is order 1 or
more, where N = the molecules that decide together. With no bias the outcome is chance (50/50
across runs). In counting form: a racemic set of N molecules has a chance L-R excess of about
sqrt(N) molecules, so an added excess of m molecules wins only when m / sqrt(N) >~ 1.
(The (k tau)^(1/4) sweep-speed factor in the full formula cannot be read from these papers,
so only this sqrt(N) counting form is tested.)

The Soai reaction (pyrimidine-5-carbaldehyde + diisopropylzinc) makes a product that catalyses
its own formation with the same hand, so it is a real self-copying amplifier.

This is a POST-HOC reanalysis of numbers stated in abstracts / open-access text. It is not
pre-registered. Expectations written BEFORE running, after reading the sources:
  S1 Spontaneous runs (no chiral additive) split consistent with 50/50: two-sided exact
     binomial p > 0.05 for each source and pooled.
  S2 The pooled runs rule out a hidden lab-wide bias strong enough to steer by the model's
     own standard (Delta >= 1, i.e. P(S) >= 0.84 or <= 0.16): chance of the observed split
     under such a bias < 1e-6.
  S3 Side-by-side identical pairs (Singleton & Vo): if chance decides, pairs disagree half
     the time, so 12/25 opposite should give p > 0.05 against 0.5. A bias shared within a
     pair (same reagent batch) of Delta >= 1 would make pairs agree (P(opposite) = 2*0.84*0.16
     = 0.27); expect the 12/25 result to reject that at p < 0.05. Uncertain: n = 25 is small.
  S4 "About 60,000 molecules" controlling the outcome means, by the counting form, that the
     hand is decided among N_d <~ 60,000^2 = 3.6e9 molecules (Delta = 1), about 6e-15 mol.
     This is an INFERENCE, not a test: the open text does not give the flask amount, so
     "N_d is far below the molecules in the flask" cannot be checked here. Reported only.
  S5 Sato 2003: ee ~0.00005% (5e-7) controlled the outcome. By the counting form that needs
     the catalyst to hold N >= (1/5e-7)^2 = 4e12 molecules (6.6 pmol, ~1.5 ng at MW 232).
     The catalyst amount is not in the open text, so this is a consistency bound only.
     ASSUMPTION (not from a source): any amount a chemist weighs or doses is >= 1 microgram,
     so the bound is met. Expect: threshold mass < 1 microgram.

Sources (DOIs checked against api.crossref.org on 2026-10-07):
  SV03: Singleton & Vo 2003, Org. Lett. 5:4337, doi:10.1021/ol035605p (abstract):
        "In 54 reactions, R and S enantiomeric products were obtained 27 times each.
         Of 25 pairs of side-by-side identical reactions, 12 afforded opposite enantiomers."
        "... a very few molecules (approximately 60 000) were sufficient to control the
         enantiomeric outcome"
  S03:  Soai et al. 2003, Tetrahedron: Asymmetry 14:185, doi:10.1016/S0957-4166(02)00791-7
        (paywalled). Counts read in the open-access review Soai 2019, Proc. Jpn. Acad. B
        95:89, doi:10.2183/pjab.95.009, sec. on absolute asymmetric synthesis:
        "(S)-1c formed 19 times and (R)-1c formed 18 times" (secondary source, ref. 96 there).
  SA03: Sato, Urabe, Ishiguro, Shibata & Soai 2003, Angew. Chem. Int. Ed. 42:315,
        doi:10.1002/anie.200390105 (abstract): "initial enantiomeric excess of about
        0.00005 % was amplified to >99.5 % ee after three cycles".
        Soai 2019 review adds: amount multiplied "by a factor of ca. 630,000 times".
  MW 232 for alkanol 1c (C14H20N2O) is computed from its structure, not quoted.
Not used: Singleton & Vo 2002 JACS (doi:10.1021/ja027129o) found 48 runs all optically active
but blamed trace chiral impurities; its abstract gives no S/R split. Kawasaki 2006 (silica)
and Suzuki 2010 (amines) report "stochastic" outcomes but their counts were not readable here.
"""
from math import comb, erf, sqrt

out = []


def say(s=""):
    print(s)
    out.append(s)


def pmf(k, n, p):
    return comb(n, k) * p**k * (1 - p) ** (n - k)


def two_sided(k, n, p=0.5):
    # exact binomial test: add up every outcome no more likely than the one seen
    pk = pmf(k, n, p)
    return min(1.0, sum(pmf(i, n, p) for i in range(n + 1) if pmf(i, n, p) <= pk * (1 + 1e-9)))


def tail_le(k, n, p):
    return sum(pmf(i, n, p) for i in range(k + 1))


def tail_ge(k, n, p):
    return sum(pmf(i, n, p) for i in range(k, n + 1))


phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
P1 = phi(1.0)  # chance the biased hand wins when Delta = 1 (0.84)

# S1: spontaneous runs
runs = {"Singleton&Vo 2003": (27, 54), "Soai 2003 (via 2019 review)": (19, 37)}
s_tot = sum(s for s, n in runs.values())
n_tot = sum(n for s, n in runs.values())
runs["pooled"] = (s_tot, n_tot)
say("S1  spontaneous runs, no chiral additive")
p_s1 = {}
for name, (s, n) in runs.items():
    p_s1[name] = two_sided(s, n)
    say(f"    {name:30s} S={s:2d}/{n:2d}  fraction={s/n:.3f}  p(50/50)={p_s1[name]:.2f}")

# S2: could a lab-wide steering bias (Delta >= 1) hide in the pooled split?
p_hiS = tail_le(s_tot, n_tot, P1)  # S favoured at P=0.84, seeing <= 46 S
p_hiR = tail_ge(s_tot, n_tot, 1 - P1)  # R favoured
say(f"\nS2  pooled {s_tot}/{n_tot}: P(split | S favoured, Delta=1)={p_hiS:.1e}, "
    f"P(split | R favoured)={p_hiR:.1e}")

# S3: side-by-side pairs
p_opp_fair = two_sided(12, 25)
q_shared = 2 * P1 * (1 - P1)
p_shared = tail_ge(12, 25, q_shared)
say(f"\nS3  pairs: 12/25 opposite. vs chance (0.5): p={p_opp_fair:.2f}; "
    f"vs shared Delta=1 bias (P_opp={q_shared:.2f}): P(>=12 opposite)={p_shared:.3f}")

# S4: how many molecules decide, if ~60,000 extra molecules are enough
m = 60_000
NA = 6.022e23
for d in (1, 2):
    n_d = (m / d) ** 2
    say(f"\nS4  Delta={d}: N_d <= {n_d:.1e} molecules = {n_d/NA:.1e} mol" if d == 1 else
        f"    Delta={d}: N_d <= {n_d:.1e} molecules = {n_d/NA:.1e} mol  (inference only)")

# S5: smallest controlling ee vs counting noise
ee = 5e-7
n_need = 1 / ee**2
mass_g = n_need / NA * 232
say(f"\nS5  ee=5e-7 beats sqrt(N) noise once N >= {n_need:.1e} molecules "
    f"= {n_need/NA:.1e} mol = {mass_g*1e9:.1f} ng of alkanol")
say("    (catalyst amount not in open text; 1 microgram floor is an assumption)")

checks = {
    "S1": all(p > 0.05 for p in p_s1.values()),
    "S2": max(p_hiS, p_hiR) < 1e-6,
    "S3a": p_opp_fair > 0.05,
    "S3b": p_shared < 0.05,
    "S5": mass_g < 1e-6,
}
say("\nexpectations: " + ", ".join(f"{k} {'held' if v else 'MISSED'}" for k, v in checks.items()))
say("S4 not scored (no flask amount in open text).")
with open(__file__.replace(".py", ".log"), "w") as f:
    f.write("\n".join(out) + "\n")

assert checks["S1"], "S1: spontaneous runs should fit 50/50"
assert checks["S2"], "S2: no steering lab-wide bias should fit the pooled split"
assert checks["S3a"], "S3a: pairs should disagree about half the time"
assert checks["S3b"], "S3b: pairs should rule out a shared Delta=1 bias"
assert checks["S5"], "S5: 5e-7 ee threshold amount should be below 1 microgram"
