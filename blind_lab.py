"""Blind lab test: does the model's P = Phi(Delta), Delta linear in the bias, fit two published ladders of deliberate bias?

WHY.  blind_check.py tested only the "chance" zone (no deliberate bias; hands split 50/50).  The model's other claim is
the shape of the transition: the chance that the biased hand wins is Phi(c * b * sqrt(N)), so at fixed N it rises with
the bias b as a probit that is LINEAR in b (from 0.6 to 0.977 over a factor of about 8 in b), and it never falls below
one half.  Two papers ran ladders of deliberate bias with replicate runs at each step:
 Hawbaker & Blackmond 2019 (Soai autocatalysis, isotopically chiral initiator at nine ee levels, 5-12 runs each, plus
 achiral and racemic controls), and Choi et al. 2024 (NaClO3 Viedma ripening, five initial crystal excesses, 10-17 runs
 each).

PROTOCOL.  The author has not read either paper.  Setup only (conditions, ee levels, run counts, where the outcomes are)
was read by separate scout agents told to withhold every outcome; the predictor saw only their setup reports.  These
predictions are committed before anyone reads an outcome.  A separate reader then enters the run outcomes below, in a
later commit, without changing anything above them.  Disclosed exposure: (i) Section 3.7 of the paper already cites
Hawbaker & Blackmond's qualitative conclusion (lab Soai tests need biases far larger than the PVED); (ii) both titles
were seen, and the Choi title ("Chiral flipping in Viedma deracemization") says some runs end on the initial minority
hand; (iii) the predictor is an AI model (Claude) that may have seen these papers in training; it does not recall their
run counts.  The predictions below were written by the AI from the model; the author reviewed them.

WHAT THE MODEL CAN AND CANNOT SAY.
 Hawbaker: the initiator's isotopic chirality (OCD3 vs OCH3) passes on an unknown fraction s of its ee to the first
 product, so b = s * e and Delta = K * e with K = c * s * sqrt(N), all unknown.  The model therefore cannot predict the
 ee where control sets in; it predicts the SHAPE (one K for all nine levels, P linear-probit in e) and the zero point
 (controls 50/50).  ASSUMPTION: the hand is decided in the first stage, which is the same in every batch (0.1 mmol
 aldehyde and 0.1 mmol initiator 3 in 0.5 mL, SI Fig. 3 caption), so K is one number for all batches; later injections
 (3-7) only amplify.
 Choi: viedma_sim.py's mechanism gives P = Phi(a * eps * sqrt(N)), eps = |CEE_i|, a = 0.565 (95% CI 0.505-0.629), with
 N = the starting crystals.  The scout found two counts: 100 (the level labels 45L/55D ... 25L/75D) and "400 in all
 experiments" (Fig. 3 caption).  The pre-registered range is N = 100-400.  Grinding (10 h) multiplies the particle
 number; if the decision were taken among the ground fragments, N would be far larger and the model would predict no
 flips at all.
 Expected flips (runs ending on the initial minority), a = 0.565, printed by the script:
   CEE_i 0.1 (17 runs): 4.9 at N = 100, 2.2 at N = 400;  0.2 (10): 1.3 / 0.12;  0.3 (10): 0.45 / 0.003;
   0.4 (10): 0.12 / 0.000;  0.5 (10): 0.02 / 0.000.  Total 6.7 / 2.3 of 57.

SOURCES (DOIs checked against api.crossref.org, 2026-10-16).
 HB19: Hawbaker, N.A., Blackmond, D.G. 2019, Nature Chemistry 11, 957-962, doi:10.1038/s41557-019-0321-y; outcomes from
       its Supplementary Information (41557_2019_321_MOESM9_ESM.pdf), Supp. Table 1 (chiral initiator, per run) and Supp.
       Table 2 (controls, per run).  The main text's Table 1 lists slightly different levels (16% R, 1% R); the SI's
       per-run tables govern and any count difference is reported.
 CH24: Choi et al. 2024, J. Phys. Chem. Lett. 15, 4367-4374, doi:10.1021/acs.jpclett.4c00558; outcomes from Fig. 4A.

PREDICTIONS (written 2026-10-16, before any outcome of either paper was read).  "Follows" means a run's final hand is the
hand the bias favours (the initiator's hand; Choi: the initial majority).  Runs with no final hand (HB19 Direction "-";
Choi: not homochiral at the end) are left out of the fits and counted.  A level with more than half its runs without a
hand is dropped; if fewer than 5 HB19 levels or 3 Choi levels remain, that paper's H/C items are NOT SCORABLE.
 (H1) HB19 shape: P(R) = Phi(beta0 + K * e), e = signed initiator ee (R > 0), one K and one lab offset beta0 for all nine
      levels, fits: parametric-bootstrap goodness of fit (deviance, 2,000 sets, seed 17) p > 0.05.
 (H2) HB19 width: the model's P(follow) = Phi(K |e|) (1 parameter) is not worse than a log-threshold curve
      P(follow) = 0.5 + 0.5 Phi(u + v log10|e|) (2 parameters, transition width free) by more than 2.0 log-likelihood
      units.  If the transition is spread over more than a factor of about 8 in ee, this fails.
 (H3) HB19 no anti-steering: no level's follow fraction is below 1/2 at one-sided binomial p < 0.05/9.
 (H4) HB19 controls (achiral + racemic initiator, the model's chance zone), pooled: the split fits 50/50 (two-sided p >
      0.05) and rules out a lab-wide steering bias of Delta >= 1 (chance of the split under it < 1e-3).
 (H5) HB19 controls, no shared-batch steering: the dispersion D (as blind_check.py B3) over control groups with >= 2 runs
      with a hand has P(D >= observed | independent coin flips) > 0.05 (200,000 sets, seed 17).
 (C1) Choi N: with a = 0.565, the maximum-likelihood N_eff in P(follow) = Phi(a eps sqrt(N_eff)) has a 95% profile CI
      that, widened by a's CI (x 0.807 low, x 1.252 high), overlaps N = 100-400.
 (C2) Choi shape: that one-parameter curve fits all five levels: parametric-bootstrap goodness of fit p > 0.05.
 (C3) Choi no anti-steering: no level's follow fraction is below 1/2 at one-sided binomial p < 0.05/5.
 Reported, not predicted: HB19 K, beta0 and e_1 = 1/K (the ee at which the initiator wins 84% of runs); Choi N_eff.
 Readings fixed in advance.  H1 or H2 fails: control does not rise as a probit linear in the bias, the model's core
 form.  H3 or C3 fails: the disfavoured hand wins more than chance, which the model cannot produce.  H4 or H5 fails: the
 lab or a batch had a stray steering bias.  C1 fails high (N_eff >> 400, few flips): the decision is taken among many
 more units than the starting crystals (e.g. the ground fragments).  C1 fails low (N_eff << 100): a few clusters decide,
 or viedma_sim's a does not transfer to this setup.  Every outcome is reported as found.

RESULTS (blind_lab.log, 2026-10-16): 4 of 8 held (H3, H5, C2, C3); 4 FAILED (H1, H2, H4, C1).  Counts were read by a
separate reader after the predictions commit (9f0b6a8) and entered unchanged.
 HB19, 95 initiator runs: the initiator's hand won every run with a hand at 96, 85, 52, 5 and 1.1% ee, lost 1 of 11 at
 39% S, and split 5/5, 6/3 and 6/4 (follow/against) at 0.16% R, 0.10% R and 0.10% S; 6 runs had no hand.  H1 and H2
 FAILED: no single probit linear in ee fits (p = 0.014), and a curve with free width fits better by 6.4 log units
 (fit K = 5.2, e_1 = 19% ee).  Both failures come from the one loss at 39% S: a curve steep enough to give 16 of 16 wins at
 5% and 1.1% makes a loss at 39% all but impossible.  Post hoc (printed, not scored): without that run the linear curve
 fits (p = 0.50, width gap -0.15) with e_1 = 0.48% ee.  The run stays in the score.  Scoring H14, H16 and H23 by the sign
 of their product ee instead of the Direction column (post hoc, not in the script) makes H1 worse (p = 0.002).  H3 held:
 no level went against the initiator.
 HB19 controls, 49 runs: only 21 ended with a hand, 5 R and 16 S.  H4 FAILED (p = 0.027 against 50/50): runs without a
 chiral initiator leaned S, and a steering bias of Delta >= 1 is not ruled out (P = 0.23).  H5 held (P = 0.099): the lean
 is not concentrated in one batch.
 HB19 main text vs SI (reader's note): the main-text Table 1 shows 10 runs at 85% S (SI 11), a 16% R row with no SI
 counterpart, and no 39% S row.  The SI was used, as pre-registered.
 Choi, 57 runs: flips 7/17, 3/10, 2/10, 1/10, 0/10 at CEE_i 0.1-0.5, 13 in all (expected 6.7 at N = 100, 2.3 at N = 400).
 C2 held: one curve Phi(a eps sqrt(N_eff)) fits all five levels (p = 0.90).  C1 FAILED low: N_eff = 34 (95% CI 12-76,
 9-95 widened for a), below the 100-400 starting crystals.  C3 held.  The paper reports flips only, so any runs without a
 final hand are not visible here.
 Code-only addition after the first run: the POST HOC line.  No prediction, threshold or count was changed.
 Reading: the Viedma ladder has the predicted shape, a probit in eps * sqrt(N), but chance matters about 3-12 times more
 than the simulation says: either fewer units decide than there are starting crystals, or a is 0.17-0.33 rather than
 0.565.  The Soai ladder fails the shape test as scored because of one run, and its no-initiator controls leaned S, a stray
 lab bias of the kind the model allows but cannot predict.  Not touched: in both labs, steering needs biases of order
 1e-3 to 1e-1 (Soai: 0.10-0.16% ee is near chance, 1.1% ee wins every run; Viedma: a 10% excess wins 59% of runs), far
 above any PVED.

Run: /opt/anaconda3/bin/python blind_lab.py > blind_lab.log
"""
from math import comb, erf, sqrt

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr, xlogy

RES = {}
phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
A0, A_LO, A_HI = 0.565, 0.505, 0.629   # viedma_sim.log base a and 95% CI

# Setup (from the scouts, before outcomes): HB19 signed ee and run count per level (SI Table 1); Choi |CEE_i| and runs.
HB_LEVELS = {"96R": (0.96, 11), "85S": (-0.85, 11), "52R": (0.52, 11), "39S": (-0.39, 11), "5R": (0.05, 12),
             "1.1R": (0.011, 5), "0.16R": (0.0016, 10), "0.10R": (0.0010, 12), "0.10S": (-0.0010, 12)}
CH_LEVELS = {0.1: 17, 0.2: 10, 0.3: 10, 0.4: 10, 0.5: 10}

# Outcomes: entered after the predictions commit (9f0b6a8), as read by a separate reader on 2026-10-16, unchanged.
# HB19 counts use the SI's Direction column as given (H14, H16, H23 disagree with the sign of their product ee).
HB = {"96R": (11, 0, 0), "85S": (0, 11, 0), "52R": (11, 0, 0), "39S": (1, 10, 0), "5R": (12, 0, 0),   # (n_R, n_S, n_none)
      "1.1R": (4, 0, 1), "0.16R": (5, 5, 0), "0.10R": (6, 3, 3), "0.10S": (4, 6, 2)}                 # SI Table 1
HB_CTRL = {"A achiral": (0, 0, 8), "B achiral": (0, 2, 5), "D achiral": (1, 0, 5), "F racemic": (0, 3, 0),  # SI Table 2
           "F achiral": (0, 3, 0), "G racemic": (1, 2, 0), "G achiral": (0, 3, 0), "H racemic": (1, 1, 6),
           "H achiral": (2, 2, 4)}
CH = {0.1: (10, 7, 0), 0.2: (7, 3, 0), 0.3: (8, 2, 0), 0.4: (9, 1, 0), 0.5: (10, 0, 0)}   # Fig. 4A (n_follow, n_flip,
# n_none); the paper reports flips only, so n_none = 0 means "none reported", not "none confirmed"


def check(p, ok, detail):
    ok = ok if isinstance(ok, str) else bool(ok)
    RES[p] = ok
    print(f"  {p} {'HELD' if ok is True else 'FAILED' if ok is False else ok}: {detail}")


def pmf(k, n, p):
    return comb(n, k) * p**k * (1 - p) ** (n - k)


def two_sided(k, n, p=0.5):
    pk = pmf(k, n, p)
    return min(1.0, sum(pmf(i, n, p) for i in range(n + 1) if pmf(i, n, p) <= pk * (1 + 1e-9)))


def below_half(k, n):
    return sum(pmf(i, n, 0.5) for i in range(k + 1))


def steer(k, n):
    return sum(pmf(i, n, phi(1.0)) for i in range(max(k, n - k) + 1))


def loglik(k, n, P):
    P = np.clip(P, 1e-12, 1 - 1e-12)
    return float(np.sum(xlogy(k, P) + xlogy(n - k, 1 - P)))


def deviance(k, n, P):
    P = np.clip(P, 1e-12, 1 - 1e-12)
    return float(2 * np.sum(xlogy(k, k / (n * P)) + xlogy(n - k, (n - k) / (n * (1 - P)))))


def best(f, starts):
    r = min((minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000)) for x0 in starts),
            key=lambda r: r.fun)
    return r.x, -r.fun


def hb_fit(e, kR, n):
    """P(R) = Phi(beta0 + K e); x = (beta0, log10 K), log10 K clipped to [-1, 7]."""
    f = lambda x: -loglik(kR, n, ndtr(x[0] + 10 ** np.clip(x[1], -1, 7) * e))
    return best(f, [(0.0, 1.0), (0.0, 2.0), (0.0, 3.0)])


def hb_shape(e, kR, n, draws=2000, seed=17):
    x, _ = hb_fit(e, kR, n)
    P = ndtr(x[0] + 10 ** np.clip(x[1], -1, 7) * e)
    d_obs = deviance(kR, n, P)
    sim = np.random.default_rng(seed).binomial(n, P, size=(draws, len(n)))
    d_sim = []
    for s in sim:
        y, _ = hb_fit(e, s, n)
        d_sim.append(deviance(s, n, ndtr(y[0] + 10 ** np.clip(y[1], -1, 7) * e)))
    return x, d_obs, float((np.array(d_sim) >= d_obs - 1e-9).mean())


def hb_width(ae, kf, n):
    _, l_lin = best(lambda x: -loglik(kf, n, ndtr(10 ** np.clip(x[0], -1, 7) * ae)), [(1.0,), (2.0,), (3.0,)])
    _, l_log = best(lambda x: -loglik(kf, n, 0.5 + 0.5 * ndtr(x[0] + x[1] * np.log10(ae))),
                    [(u, v) for u in (-2.0, 2.0, 6.0) for v in (1.0, 3.0)])
    return l_lin, l_log


NGRID = np.logspace(0, 7, 3501)


def ch_curve(eps, a=A0):
    return ndtr(a * np.outer(np.sqrt(NGRID), eps))   # grid x levels


def ch_fit(eps, kf, n):
    P = np.clip(ch_curve(eps), 1e-12, 1 - 1e-12)
    ll = (xlogy(kf, P) + xlogy(n - kf, 1 - P)).sum(1)
    i = int(ll.argmax())
    ok = ll >= ll[i] - 1.92
    return NGRID[i], NGRID[ok].min(), NGRID[ok].max(), P


def ch_shape(eps, kf, n, draws=5000, seed=17):
    Nb, lo, hi, P = ch_fit(eps, kf, n)
    i = int(np.argmin(abs(NGRID - Nb)))
    d_obs = deviance(kf, n, P[i])
    sim = np.random.default_rng(seed).binomial(n, P[i], size=(draws, len(n)))
    ll = xlogy(sim[:, None, :], P[None]) + xlogy(n - sim[:, None, :], 1 - P[None])   # draws x grid x levels
    j = ll.sum(2).argmax(1)
    Pj = P[j]
    d_sim = 2 * (xlogy(sim, sim / (n * Pj)) + xlogy(n - sim, (n - sim) / (n * (1 - Pj)))).sum(1)
    return Nb, lo, hi, d_obs, float((d_sim >= d_obs - 1e-9).mean())


def dispersion(rows, draws=200_000, seed=17):
    n = np.array([a + b for a, b in rows])
    k = np.array([a for a, b in rows])
    d_obs = (((k - n / 2) ** 2) / (n / 4)).sum()
    sim = np.random.default_rng(seed).binomial(n, 0.5, size=(draws, len(n)))
    return d_obs, ((((sim - n / 2) ** 2) / (n / 4)).sum(1) >= d_obs - 1e-9).mean()


def setup_report():
    eps = np.array(list(CH_LEVELS))
    print("Choi, expected flips per level at a = 0.565 (N = 100 / N = 400):")
    tot = [0.0, 0.0]
    for e, n in CH_LEVELS.items():
        f = [n * (1 - phi(A0 * e * sqrt(N))) for N in (100, 400)]
        tot = [t + x for t, x in zip(tot, f)]
        print(f"  CEE_i {e:.1f}, {n:2d} runs: {f[0]:.2f} / {f[1]:.3f}")
    print(f"  total {tot[0]:.1f} / {tot[1]:.1f} of {sum(CH_LEVELS.values())}")
    print(f"HB19: {len(HB_LEVELS)} levels, {sum(v[1] for v in HB_LEVELS.values())} chiral-initiator runs; "
          f"the model's transition P(follow) 0.6 -> 0.977 spans a factor {2.0 / 0.2533:.1f} in ee")
    return eps


def main():
    assert two_sided(24, 48) > 0.99 and steer(24, 48) < 1e-6 and below_half(0, 12) < 0.05 / 9
    eps = setup_report()
    if HB is None or HB_CTRL is None or CH is None:
        print("Predictions committed; run outcomes not yet entered (see the docstring).")
        return

    print("\n=== Hawbaker & Blackmond 2019, SI Tables 1-2 ===")
    keep = []
    for lab, (r, s, z) in HB.items():
        e, n0 = HB_LEVELS[lab]
        flag = "" if r + s + z == n0 else f"  (setup said {n0} runs)"
        drop = z > (r + s + z) / 2
        print(f"  {lab:6s} e = {e:+.4f}: R {r}, S {s}, no hand {z}{' DROPPED' if drop else ''}{flag}")
        if not drop:
            keep.append((e, r, s))
    if len(keep) < 5:
        for p in ("H1", "H2", "H3"):
            check(p, "NOT SCORABLE", f"{len(keep)} levels with a hand")
    else:
        e = np.array([x[0] for x in keep])
        kR = np.array([x[1] for x in keep])
        n = np.array([x[1] + x[2] for x in keep])
        x, d, pg = hb_shape(e, kR, n)
        K = 10 ** np.clip(x[1], -1, 7)
        print(f"  fit: beta0 = {x[0]:+.3f}, K = {K:.3g} (e_1 = 1/K = {100 / K:.3g}% ee)")
        check("H1", pg > 0.05, f"deviance {d:.2f} over {len(e)} levels, bootstrap p = {pg:.3f}")
        kf = np.where(e > 0, kR, n - kR)
        l_lin, l_log = hb_width(abs(e), kf, n)
        check("H2", l_log - l_lin <= 2.0, f"logL linear {l_lin:.2f}, log-threshold {l_log:.2f}, gap {l_log - l_lin:.2f}")
        worst = min(below_half(int(a), int(b)) for a, b in zip(kf, n))
        check("H3", worst >= 0.05 / 9, f"smallest one-sided p(follow < 1/2) = {worst:.3g}")
        # POST HOC, added after the results were seen; not a prediction and not scored.
        k2 = kR.copy()
        k2[e == -0.39] = 0
        y, d2, p2 = hb_shape(e, k2, n)
        l1, l2 = hb_width(abs(e), np.where(e > 0, k2, n - k2), n)
        print(f"  POST HOC, without the one R run at 39% S: K = {10 ** y[1]:.3g} (e_1 = {100 / 10 ** y[1]:.3g}% ee), "
              f"beta0 = {y[0]:+.2f}, GOF p = {p2:.3f}, width gap {l2 - l1:.2f}")
    rows = [(r, s) for r, s, z in HB_CTRL.values()]
    a, b = sum(r for r, _ in rows), sum(s for _, s in rows)
    for lab, (r, s, z) in HB_CTRL.items():
        print(f"  control {lab:12s} R {r}, S {s}, no hand {z}")
    check("H4", two_sided(a, a + b) > 0.05 and steer(a, a + b) < 1e-3,
          f"controls {a} R / {b} S, p(50/50) = {two_sided(a, a + b):.3f}, P(split | Delta = 1) = {steer(a, a + b):.1e}")
    d, pd = dispersion([x for x in rows if sum(x) >= 2])
    check("H5", pd > 0.05, f"D = {d:.2f} over {sum(sum(x) >= 2 for x in rows)} groups, P = {pd:.3f}")

    print("\n=== Choi et al. 2024, Fig. 4A ===")
    keep = []
    for e, (f, g, z) in CH.items():
        drop = z > (f + g + z) / 2
        print(f"  CEE_i {e:.1f}: follow {f}, against {g}, no hand {z}{' DROPPED' if drop else ''}"
              f"{'' if f + g + z == CH_LEVELS[e] else f'  (setup said {CH_LEVELS[e]} runs)'}")
        if not drop:
            keep.append((e, f, g))
    if len(keep) < 3:
        for p in ("C1", "C2", "C3"):
            check(p, "NOT SCORABLE", f"{len(keep)} levels with a hand")
    else:
        ee = np.array([x[0] for x in keep])
        kf = np.array([x[1] for x in keep])
        n = np.array([x[1] + x[2] for x in keep])
        Nb, lo, hi, d, pg = ch_shape(ee, kf, n)
        lo_w, hi_w = lo * (A0 / A_HI) ** 2, hi * (A0 / A_LO) ** 2
        edge = " (upper end at grid edge 1e7: unbounded)" if hi >= NGRID[-1] else ""
        check("C1", lo_w <= 400 and hi_w >= 100,
              f"N_eff = {Nb:.0f}, 95% CI {lo:.0f}-{hi:.0f}, widened {lo_w:.0f}-{hi_w:.0f}{edge}")
        check("C2", pg > 0.05, f"deviance {d:.2f} over {len(ee)} levels, bootstrap p = {pg:.3f}")
        worst = min(below_half(int(a), int(b)) for a, b in zip(kf, n))
        check("C3", worst >= 0.05 / 5, f"smallest one-sided p(follow < 1/2) = {worst:.3g}")
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}"
                                                         for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
