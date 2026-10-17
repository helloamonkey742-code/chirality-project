"""Blind check: do unbiased crystal-deracemization runs split by chance, as Section 3.7 claims for lab amplifiers?

WHY. crystal_check.py, soai_check.py and viedma_check.py reanalyse papers the author had already read, so they are not
pre-registered tests. This script is. The predictions below were written and committed BEFORE anyone on the author's side
read the results of the two papers. The run counts are entered into DECK and STEENDAM only in a later commit, so the git
history shows the order.

PROTOCOL.
 - A scout was told to report SETUP ONLY (no outcomes, no titles, no abstracts) and found the sources below. Its report
   contains none of their results.
 - The author is asked not to open either paper until the results commit.
 - After the commit, a separate reader extracts only the run-level hands, and they are entered unchanged.
 - Disclosure: the predictions were drafted by an AI assistant (Claude) that did not open either paper in this session, but
   it may have seen them in training. It recognises Steendam et al. 2014 by name, so B5 is blind for the author only. It
   recalls nothing of Deck et al. 2025.

WHAT THE MODEL SAYS. A bias b steers the hand only when Delta = c * b * sqrt(N) is order 1 or more (Eq. 11; the counting
rule of Section 3.7). With no deliberate bias, a stray bias (trace chiral impurity, PVED) steers a flask-sized pool only if
it is unusually large, so each run's hand should be set by chance. Section 3.7 states this for 91 Soai runs and 73 Strecker
runs read after the fact. Here it is tested blind. Shared stock is the risk: each Deck experiment fills its crystallizers
from one 100 mL stock, so a stray bias in that stock would make the vials of one experiment agree.

SOURCES (DOIs resolved on api.crossref.org on 2026-10-16; titles are withheld until the results are in)
 DK25: Deck, L.T. et al. 2025, J. Am. Chem. Soc., doi:10.1021/jacs.5c02651 (open text PMC12147123). This is NMPA, a
       racemizable conglomerate, with DBU as the racemization catalyst in IPA/ACN 95/5 w/w. Each crystallizer cools a
       racemic solution and then runs temperature cycles. No chiral additive. Each labelled experiment has 8 crystallizers
       of 5 mL from one 100 mL stock, stirred at 1000 rpm. The experiments are:
       - A1, A2: c0 0.035 g/g, cooled at 3.3 C/min, high catalyst;
       - B1, B2: c0 0.048 g/g, 3.3 C/min, high catalyst;
       - C: c0 0.060 g/g, 3.3 C/min, high catalyst;
       - D1, D2: c0 0.048 g/g, 3.3 C/min, low catalyst;
       - G: c0 0.048 g/g, 0.2 C/min, high catalyst;
       - H: c0 0.048 g/g, 0.2 C/min, low catalyst.
       The final hand of each crystallizer is in Table 1.
 ST14: Steendam, R.R.E. et al. 2014, Nat. Commun., doi:10.1038/ncomms6543 (open text PMC4263183). Achiral reactants make
       a conglomerate product, which is ground with glass beads (Viedma-type, DBU in EtOH). 68 identical runs, no bias.

LIMITS. This tests only the model's "no steering" side: any theory in which stray biases are small also predicts 50/50, so
passing is weak support and failing is the informative outcome. A graded-bias test (several bias sizes, run counts at each) was
the goal, but the open-access papers found report one run or an ensemble ee per bias level (Belletti et al. 2018,
doi:10.1021/acs.cgd.8b00856; Hananel et al. 2019, doi:10.1073/pnas.1821923116), so the zone rule could not be scored there.

MODEL ZONE FOR A DECK VIAL (printed by the script). The NMPA in one vial, at most 0.060 g/g x 5 mL x 0.785 g/mL solvent / 252.3
g/mol, is 5.6e20 molecules: an upper bound on the deciding units. Even at the top published PVED (g = 4.8e-16, Section 2.1) and
c = 20, Phi(c g sqrt(N)) = 0.5001, so the PVED sits deep in the "chance" zone. A stray bias steers a vial only if b >= 1/(c sqrt(N))
= 2e-12 to 8e-11 (more if fewer units decide).

PREDICTIONS (written 2026-10-16, before either paper's results were read). "A hand" means a run's final state is reported as
one enantiomer (DK25 Table 1 gives it in parentheses). Runs reported as racemic or undetermined are left out of B1-B3 and counted separately. If more than half of
all Deck runs have no final hand, B1-B3 are reported as NOT SCORABLE, not as held or failed.
 (B1) Deck, pooled over all experiments: the split between the two hands fits 50/50 (two-sided exact binomial p > 0.05).
 (B2) Deck, pooled: the split rules out a lab-wide steering bias at the model's own bar, Delta >= 1 (P = 0.84 for either
      hand). The chance of the split under such a bias must be < 1e-3.
 (B3) Deck, no shared-stock steering: the experiments' splits spread no more than independent coin flips would. The
      statistic is D = sum over experiments of (k - n/2)^2 / (n/4). Under independent 50/50 runs (200,000 simulated sets,
      seed 17), the chance of a D at least as large as observed must be > 0.05. Only experiments with 2 or more runs with a
      hand are included.
 (B4) Deck: no single experiment has all its runs on one hand when it has 8 runs with a hand (chance 2/256 per experiment
      under 50/50). Allowed: at most 1 such experiment over the 9. Under coin flips that limit is exceeded with chance 0.002.
 (B5) Steendam (blind for the author only): the 68 runs fit 50/50 (p > 0.05) and rule out a steering bias at Delta >= 1
      (< 1e-6).
 If B1, B2 or B5 fails, the lab had a steering bias and Section 3.7's "outcomes split as chance predicts" does not extend to
 that system. If B3 or B4 fails, a stray bias in the stock steered 5 mL vials, so b * sqrt(N) >= 1 there. Either outcome is
 reported as found.

RESULTS (blind_check.log, 2026-10-16): 5 of 5 held (B1-B5). Counts were read from the open text by a separate reader after
the predictions commit (ec19c77) and entered unchanged.
 Deck: all 72 crystallizers reached one hand: 42 L, 30 D (p = 0.195 for 50/50; the paper's own text gives the same 30/42 and
 19%). Per experiment L/D: A1 6/2, A2 4/4, B1 6/2, B2 5/3, C 6/2, D1 4/4, D2 3/5, G 5/3, H 3/5. No experiment went 8/0. The spread
 across experiments is what coin flips give (D = 8.00 against an expected 9; P = 0.579), so vials from one stock did not agree.
 A Delta >= 1 steering bias would make a split this even with chance 1.6e-7.
 Steendam: 39 S, 29 R of 68 (p = 0.275; Delta >= 1 bias ruled out, 1.4e-7). Not scored: the same paragraph says runs with freshly
 prepared starting materials again went to S more often, with no counts given.
 Code-only fix after the first run: check() now converts numpy bools, which had made the tally line print "4 of 5" and B3 "True".
 No prediction or threshold was changed.
 Reading: in two systems the model's zone rule puts in the "chance" zone, the hands split as chance predicts, at the level
 of a whole lab and of a shared stock. This is weak support (see LIMITS): it shows no hidden steering bias of Delta >= 1 at the
 5 mL scale, not that the sqrt(N) rule is right.

Run: /opt/anaconda3/bin/python blind_check.py > blind_check.log
"""
from math import comb, erf, sqrt

import numpy as np

# Entered after the predictions were committed (ec19c77), as read from the open text by a separate reader on 2026-10-16.
DECK = {"A1": (6, 2, 0), "A2": (4, 4, 0), "B1": (6, 2, 0), "B2": (5, 3, 0), "C": (6, 2, 0),   # (L, D, no hand), DK25 Table 1
        "D1": (4, 4, 0), "D2": (3, 5, 0), "G": (5, 3, 0), "H": (3, 5, 0)}
STEENDAM = (39, 29)   # (S-1, R-1) enantiopure end states, ST14 Results ("Chiral outcome and rate of the reaction")

RES = {}
phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
P1 = phi(1.0)   # chance the favoured hand wins at Delta = 1 (0.84)


def check(p, ok, detail):
    ok = ok if isinstance(ok, str) else bool(ok)   # numpy bools are not `True`; code-only fix after the first run
    RES[p] = ok
    print(f"  {p} {'HELD' if ok is True else 'FAILED' if ok is False else ok}: {detail}")


def pmf(k, n, p):
    return comb(n, k) * p**k * (1 - p) ** (n - k)


def two_sided(k, n, p=0.5):
    pk = pmf(k, n, p)
    return min(1.0, sum(pmf(i, n, p) for i in range(n + 1) if pmf(i, n, p) <= pk * (1 + 1e-9)))


def steer(k, n):
    """Chance the majority hand wins no more than observed if it were favoured at Delta = 1 (worst case over which hand)."""
    return sum(pmf(i, n, P1) for i in range(max(k, n - k) + 1))


def dispersion(rows, draws=200_000, seed=17):
    n = np.array([a + b for a, b in rows])
    k = np.array([a for a, b in rows])
    d_obs = (((k - n / 2) ** 2) / (n / 4)).sum()
    sim = np.random.default_rng(seed).binomial(n, 0.5, size=(draws, len(n)))
    return d_obs, ((((sim - n / 2) ** 2) / (n / 4)).sum(1) >= d_obs - 1e-9).mean()


def main():
    assert two_sided(36, 72) > 0.99 and steer(36, 72) < 1e-9 and steer(10, 72) == steer(62, 72)
    n_vial = 0.060 * 5 * 0.785 / 252.3 * 6.022e23
    print(f"Deck vial: N <= {n_vial:.1e} molecules; PVED g = 4.8e-16, c = 20: P = {phi(20 * 4.8e-16 * sqrt(n_vial)):.4f} "
          f"(chance zone); steering needs b >= {1 / (20 * sqrt(n_vial)):.0e} to {1 / (0.5 * sqrt(n_vial)):.0e}")
    for lab, n in (("Deck pooled, all 72 with a hand", 72), ("Steendam, 68 runs", 68)):
        ok1 = [k for k in range(n + 1) if two_sided(k, n) > 0.05]
        ok2 = [k for k in range(n + 1) if steer(k, n) < (1e-3 if n == 72 else 1e-6)]
        print(f"{lab}: fits 50/50 if hand 1 = {ok1[0]}-{ok1[-1]}; rules out Delta >= 1 if hand 1 = {ok2[0]}-{ok2[-1]}")
    if DECK is None or STEENDAM is None:
        print("Predictions committed; run counts not yet entered (see the docstring).")
        return
    print("=== Deck et al. 2025, Table 1 ===")
    for lab, (a, b, z) in DECK.items():
        print(f"  {lab:3s} hand 1 {a}, hand 2 {b}, no hand {z}")
    a = sum(v[0] for v in DECK.values())
    b = sum(v[1] for v in DECK.values())
    z = sum(v[2] for v in DECK.values())
    n = a + b
    if z > (n + z) / 2:
        for p in ("B1", "B2", "B3", "B4"):
            check(p, "NOT SCORABLE", f"{z} of {n + z} runs have no final hand")
    else:
        check("B1", two_sided(a, n) > 0.05, f"pooled {a}/{n}, p(50/50) = {two_sided(a, n):.3f}")
        s = steer(a, n)
        check("B2", s < 1e-3, f"P(split | Delta = 1 bias) = {s:.1e}")
        rows = [(x, y) for x, y, _ in DECK.values() if x + y >= 2]
        d, pd = dispersion(rows)
        check("B3", pd > 0.05, f"D = {d:.2f} over {len(rows)} experiments, P(D >= observed | coin flips) = {pd:.3f}")
        full = sum(1 for x, y, _ in DECK.values() if x + y == 8 and min(x, y) == 0)
        check("B4", full <= 1, f"{full} experiment(s) with all 8 runs on one hand")
    sa, sb = STEENDAM
    print("\n=== Steendam et al. 2014 ===")
    check("B5", two_sided(sa, sa + sb) > 0.05 and steer(sa, sa + sb) < 1e-6,
          f"{sa}/{sa + sb}, p(50/50) = {two_sided(sa, sa + sb):.3f}, P(split | Delta = 1) = {steer(sa, sa + sb):.1e}")
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}"
                                                         for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
