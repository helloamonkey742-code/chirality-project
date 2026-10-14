"""Viedma ripening of NaClO3: stochastic simulation of P(initial majority wins) vs eps and N.

QUESTION. For the planned bench experiment (journals/viedma_plan.md), does the chance that the
initial majority hand wins depend on Delta = eps*sqrt(N) as P = Phi(Delta) (the paper's model),
on neutral drift P = (1+eps)/2, or on P = Phi(a*Delta) with a fitted a? eps = (N_maj-N_min)/N,
N = N_0 = starting crystal count. This script simulates the mechanism and fits the three forms.

MODEL (one run = one slurry). Units: length = initial mean radius, time = ripening time unit
(dr/dt = O(1) for a crystal of size 1). All crystals are NaClO3 cubes-as-spheres, hand h in {0,1}.
 1. Solution is achiral: ONE concentration u (excess over flat-face saturation, in Gibbs-Thomson
    units), shared by both hands. NaClO3 has no handedness once dissolved (Viedma 2005,
    doi:10.1103/PhysRevLett.94.065504), so no solution-phase racemization step is needed.
 2. Growth/dissolution (Gibbs-Thomson / Ostwald ripening; mean-field framework of Lifshitz &
    Slyozov 1961, doi:10.1016/0022-3697(61)90054-3): dr/dt = u - 1/r. A crystal grows if
    u > 1/r, dissolves otherwise. ASSUMPTION: interface-reaction-limited rate (Lifshitz-Slyozov
    use the diffusion-limited form dr/dt ~ (u - 1/r)/r); chosen for simplicity.
    ASSUMPTION: solution reservoir is negligible and relaxes instantly, so u is fixed each step
    by solid-mass conservation: u = sum(r)/sum(r^2). Mass is re-normalised exactly every step.
    Crystals reaching r < 1e-3 are removed (their mass is spread over all crystals by the
    normalisation; hand-blind).
 3. Attrition/grinding (Viedma 2005; Noorduin et al. 2008 Cryst. Growth Des. 8:1675,
    doi:10.1021/cg701211a; Iggland & Mazzotti 2011 Cryst. Growth Des.,
    doi:10.1021/cg2008599 use breakage of this kind): a crystal with r > R_MIN breaks with
    probability 1-exp(-LAMBDA r^2 dt) per step (rate ~ cross-section, ASSUMPTION; the source
    papers do not fix an exponent). A break chips off a same-hand piece of mass fraction
    x ~ U(0.05, 0.5) (ASSUMPTION), mass conserved. LAMBDA = 4 is ASSUMPTION calibrated to ONE
    qualitative fact only: ripening without grinding must be >= ~10x slower than ripening with
    grinding (stirring-only flasks kept their ee for days, grinding flasks finished in hours;
    Viedma 2005, see viedma_check.py V3). Checked at Delta=0, N=1000: no grinding decides at
    t~99, LAMBDA=4 at t~10 (calibration runs, disclosed below).
 4. Same-hand chip retention (the chiral feedback). Ricci, Stillinger & Debenedetti 2013, J. Chem.
    Phys., doi:10.1063/1.4827478 (abstract): size-dependent solubility ALONE is insufficient to
    reproduce Viedma ripening; a solid-phase chiral feedback such as agglomeration is needed.
    Saito & Hyuga 2009, JPSJ 78:104001, doi:10.1143/JPSJ.78.104001 (abstract): with grinding, a
    chance imbalance in cluster numbers is amplified exponentially and the time to homochirality
    does not grow with system size. I implement the minimum feedback consistent with both:
    each chipped piece of hand H is re-attached to a random same-hand crystal (chosen with
    probability ~ r^2) with probability p = S*xA/(1+S*xA), xA = area share of hand H; otherwise
    it stays a free crystal and ripens (usually dissolves, its mass going to all crystals by
    area). The FORM of p and S = 2 are ASSUMPTIONS with no source; they are the most uncertain
    part of the model, so the sensitivity runs vary them (S=0 removes the feedback).
 5. Initial state: N_0 crystals, radii iid lognormal (sigma_ln = 0.25, mean radius 1; a sieved
    fraction, ASSUMPTION), SAME distribution for both hands (size-matched, plan sec. 4). The
    majority hand has round(N(1+eps)/2) crystals; which hand (0/1) is the majority is random
    (mirror arm). Because sizes are iid, subsampling noise (plan sec. 4) is included.
 6. Time stepping: dt = 0.01, one Euler step for ripening then Bernoulli breaks (tau-leaping).
    Stop when one hand holds > 99% of solid mass ("decided", winner = that hand), else at
    t = 60 ("undecided", never dropped; also recorded: which hand holds >50% of mass).
 7. Seeds: np.random.SeedSequence([20261013, N, round(10*Delta), rep, variant]); logged by
    construction. Two worker processes, one BLAS/OMP thread each (2 threads total), niced.

STATISTICS. Cells = (variant, N, Delta). Decided runs only enter the fits (undecided counts are
reported per cell). Fits use Delta_actual = eps_actual*sqrt(N) (rounding of the crystal count).
 H_Phi:   P = Phi(Delta)         (no free parameter)
 H_a:     P = Phi(a*Delta)       (a by maximum likelihood; 95% CI from the profile likelihood)
 H_drift: P = (1+eps)/2          (no free parameter)
 H_aN:    separate a for each N  (collapse test = likelihood ratio H_a vs H_aN, chi2, df = nN-1)
Delta=0 cells carry no information on a but are tested against 50% (Clopper-Pearson).

PREDICTIONS (written 2026-10-13 BEFORE the real grid; base variant, S=2, LAMBDA=4):
 Disclosure of what was seen before writing these: (i) a handful of tiny smoke runs at N=200
 (3 seeds at Delta=0 and 3 at Delta=2; all 3 Delta=2 runs were won by the majority, all runs
 decided in t<26), (ii) calibration runs at N=1000, Delta=0, 8 seeds per LAMBDA in
 {0, 0.5, 2, 8}, which showed decision time and crystal count only, and single runs at
 N=250/16000/64000 to time the code, and (iii) this script's own --smoke mode (N=200, 3 reps
 per Delta, 12 runs, seeds from this script): it ran end to end, all 12 decided in t<16; the
 majority won 1/3, 2/3, 3/3 at Delta = 0.5, 1, 2 (too few to mean anything). I tabulated no
 larger P(win) vs Delta before writing this.
 P1  Delta = 0 control: the (arbitrarily labelled) "majority" hand wins 50% within the 95%
     Clopper-Pearson interval at every N, and hand 0 wins 50% pooled over all N (symmetry).
 P2  Undecided fraction <= 2% in every base cell (amplification from a 1/sqrt(N) imbalance
     finishes well inside t = 60).
 P3  P(win) rises with Delta at every N. Expected pooled values (from the Gaussian-noise
     argument: linear gain + counting noise ~ 1/sqrt(N) gives exactly Phi(a*Delta)):
     Delta=0.5: 0.60-0.75; Delta=1: 0.72-0.90; Delta=2: >= 0.90.
 P4  Collapse: the curves for the five N values (250, 1000, 4000, 16000, 64000) fall on one
     curve of Delta. Likelihood-ratio test H_a vs H_aN: p > 0.05 (my confidence ~50%; the
     largest deviation, if any, is expected at N=250 where eps is up to 0.13 and the
     Gaussian form saturates). No monotonic trend of the per-N a with N.
 P5  Fitted a: point estimate 0.9, range 0.5-1.5. Reasoning: initial size-sampling noise alone
     gives a = 1/CV_mass = 1/0.87 = 1.15; dynamic breakage noise lowers it. Prediction that
     H_Phi (a = 1) is NOT rejected at 5% (likelihood-ratio chi2_1 < 3.84; weak prediction,
     ~55% confidence).
 P6  Neutral drift is rejected decisively: log10 L(H_a) - log10 L(H_drift) > 10 (drift gives
     P <= 0.56 in every cell, the data at Delta = 2 give >= 0.9).
 P7  Decision time is nearly N-independent (Saito & Hyuga): median t_decide at N=64000 is
     between 1x and 3x that at N=250.
 P8  Sensitivity, S = 0 (no chiral feedback), N = 1000: most runs undecided at t = 60
     (> 50%); the majority-mass-lean rate at Delta = 2 is <= 0.70, i.e. no Phi(Delta) curve
     (Ricci et al. 2013: ripening alone is not enough).
 P9  Sensitivity, LAMBDA = 1 (slower grinding), N = 1000: a stays within 0.7x-1.3x of the base
     a (gain and noise both scale with breakage rate); decision time rises >= 2x.
 P10 Caveat stated in advance: every one of these numbers comes from the assumptions above
     (marked ASSUMPTION), not from NaClO3 measurements; the model can only say what the
     mechanism, as written, implies for the bench. It says nothing about oceans.

RESULTS (viedma_sim.log, 2026-10-13): 2440 pre-registered runs (base: N = 250/1000/4000 x 100
reps, 16000 x 80, 64000 x 30, Delta in {0, 0.5, 1, 2}; sensitivity S=0 and LAMBDA=1 at N=1000 x
100 reps), wall time 945 s on 2 workers. Predictions above are unchanged.
 P1  HELD as written. Delta=0: majority-label wins 57/100, 54/100, 60/100, 43/80, 17/30 at
     N = 250..64000, every 95% CI contains 0.5; hand 0 wins 198/410 (CI 0.434-0.532).
     Note (post-hoc, disclosed): the pooled majority-LABEL rate at Delta=0 was 231/410 = 0.563
     (CI excludes 0.5 narrowly). 1800 fresh Delta=0 runs (N = 250, 1000, 4000; unused seeds)
     gave 915/1800 = 0.508, so I read the 0.563 as chance; not a registered check.
 P2  HELD. 0 undecided of 1640 base runs (all decided by t ~ 8-14).
 P3  PARTLY FAILED. Pooled P(win) = 0.661, 0.712, 0.863 at Delta = 0.5, 1, 2. Delta=0.5 inside
     0.60-0.75 (held); Delta=1 is 0.712, just below 0.72 (missed); Delta=2 is 0.863, below the
     >= 0.90 floor (failed; only N=64000 reached 0.90, 27/30). Per-N curves are not strictly
     increasing (N=1000: 0.73, 0.73; N=16000: 0.69, 0.68 at Delta = 0.5, 1).
 P4  HELD. Per-N a (95% CI): 0.57 (0.45-0.70), 0.61 (0.48-0.74), 0.54 (0.42-0.67), 0.51
     (0.38-0.65), 0.61 (0.38-0.86) for N = 250..64000. Collapse LRT chi2 = 1.27, df = 4,
     p = 0.87. No trend of a with N.
 P5  FAILED on the number. a_hat = 0.565, 95% CI 0.505-0.629: inside the stated range 0.5-1.5
     but far from the 0.9 point guess, and H_Phi (a = 1) is rejected (chi2_1 = 157, p << 1e-10).
 P6  HELD. log10 L(H_a)/L(H_drift) = 67.3 (and log10 L(Phi)/L(drift) = 33.2).
 P7  HELD. Median t_decide 7.75 (N=250), 9.5, 10.75, 12.5, 14.0 (N=64000): ratio 1.8, and
     growing slowly with N (about +1.5 per factor 4).
 P8  HELD. S=0: 100/100 undecided in every cell at t = 60 (mass fraction of a hand never beyond
     0.10-0.88); majority-lean rate at Delta = 2 is 0.57.
 P9  FAILED. LAMBDA=1: a = 0.32 (0.21-0.43) vs 0.61 (0.48-0.74) at N=1000 base, ratio 0.52
     (predicted 0.7-1.3); t_decide 17.25 vs 9.5, ratio 1.8 (predicted >= 2).
 POST-HOC extras (--extra mode added after the main run; NOT pre-registered, appended to the
 log): S = 0.5 gives a = 0.41 (0.30-0.53), S = 2 gives 0.61 (0.48-0.74) at the same N=1000,
 S = 8 gives a = 0.545 (0.43-0.67). No feedback variant gives a >= 0.7 and Phi(Delta) (a = 1)
 is rejected in every one.
 SUMMARY: collapse on Delta = eps*sqrt(N) held (the sqrt(N) scaling is robust inside the
 model); neutral drift is rejected; the paper's coefficient Phi(Delta) is rejected in favour of
 Phi(a*Delta), a ~ 0.3-0.6, depending on grinding rate and feedback strength. Bench value:
 90% win needs Delta ~ 2.3.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import beta as beta_dist, chi2, norm

MASTER = 20261013
SIGMA_LN, R_MIN, DT, TMAX, DONE = 0.25, 0.6, 0.01, 60.0, 0.99
BASE = dict(lam=4.0, s=2.0)
VARIANTS = {0: ("base", BASE), 1: ("noretain_S0", dict(lam=4.0, s=0.0)), 2: ("slowgrind_L1", dict(lam=1.0, s=2.0)),
            3: ("retain_S0.5", dict(lam=4.0, s=0.5)), 4: ("retain_S8", dict(lam=4.0, s=8.0))}  # 3,4: post-hoc (--extra)
DELTAS = [0.0, 0.5, 1.0, 2.0]
GRID_BASE = {250: 100, 1000: 100, 4000: 100, 16000: 80, 64000: 30}   # N -> reps
GRID_SENS = {1000: 100}


def run(N, delta, rep, variant):
    seed = np.random.SeedSequence([MASTER, N, int(round(10 * delta)), rep, variant])
    lam, s = VARIANTS[variant][1]["lam"], VARIANTS[variant][1]["s"]
    rng = np.random.default_rng(seed)
    n_maj = int(round(N * (1 + delta / np.sqrt(N)) / 2))
    maj = int(rng.integers(2))
    h = np.full(N, 1 - maj, dtype=np.int8)
    h[:n_maj] = maj
    r = rng.lognormal(-SIGMA_LN**2 / 2, SIGMA_LN, N)
    M0 = np.sum(r**3)
    eps = 2 * n_maj / N - 1
    t, step = 0.0, 0
    while True:
        # ripening: shared achiral u from mass conservation
        u = r.sum() / (r * r).sum()
        r = np.maximum(r + DT * (u - 1.0 / r), 0.0)
        keep = r > 1e-3
        if not keep.all():
            r, h = r[keep], h[keep]
        r *= (M0 / np.sum(r**3)) ** (1 / 3)
        t += DT
        step += 1
        # attrition
        r2 = r * r
        p = np.where(r > R_MIN, 1 - np.exp(-lam * r2 * DT), 0.0)
        ev = np.flatnonzero(rng.random(r.size) < p)
        if ev.size:
            mp = r[ev] ** 3
            chip = rng.uniform(0.05, 0.5, ev.size) * mp
            r[ev] = (mp - chip) ** (1 / 3)
            hc = h[ev]
            att = np.zeros(ev.size, dtype=bool)
            if s > 0:
                A = r2.sum()
                for hand in (0, 1):
                    xa = r2[h == hand].sum() / A
                    att |= (hc == hand) & (rng.random(ev.size) < s * xa / (1 + s * xa))
                m3 = r**3
                for hand in (0, 1):
                    sel = att & (hc == hand)
                    if sel.any():
                        idx = np.flatnonzero(h == hand)
                        cum = np.cumsum(r2[idx])
                        tg = idx[np.minimum(np.searchsorted(cum, rng.random(sel.sum()) * cum[-1]), idx.size - 1)]
                        np.add.at(m3, tg, chip[sel])
                r = m3 ** (1 / 3)
            new = ~att
            r = np.concatenate([r, chip[new] ** (1 / 3)])
            h = np.concatenate([h, hc[new]])
        if step % 25 == 0 or t >= TMAX:
            fm = np.sum(r[h == maj] ** 3) / np.sum(r**3)
            if fm > DONE or fm < 1 - DONE or t >= TMAX:
                dec = fm > DONE or fm < 1 - DONE
                return dict(var=variant, N=N, delta=delta, rep=rep, eps=eps, maj=maj,
                            win=(int(fm > 0.5) if dec else -1), lean=int(fm > 0.5), fm=fm,
                            t=t, ncr=int(r.size), dec=int(dec))


def _job(a):
    return run(*a)


def cp(k, n, alpha=0.05):
    lo = beta_dist.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


def binom_ll(k, n, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return k * np.log(p) + (n - k) * np.log(1 - p)


def fit_a(cells, grid=np.linspace(0.01, 6, 1200)):
    """cells: list of (k, n, Delta>0). Returns a_hat, (lo, hi) from profile likelihood."""
    ll = np.array([sum(binom_ll(k, n, norm.cdf(a * d)) for k, n, d in cells) for a in grid])
    i = ll.argmax()
    ok = grid[ll >= ll[i] - chi2.ppf(0.95, 1) / 2]
    return grid[i], (ok.min(), ok.max()), ll[i]


def analyze(rows, var, out):
    R = [r for r in rows if r["var"] == var]
    if not R:
        return {}
    Ns = sorted({r["N"] for r in R})
    out(f"\n--- variant {VARIANTS[var][0]} {VARIANTS[var][1]} ---")
    out("N      Delta  eps     reps decided undec  wins  P(win)  95%CI(CP)       Phi(D)  drift   lean_rate  med_t")
    cells = []
    for N in Ns:
        for d in DELTAS:
            c = [r for r in R if r["N"] == N and r["delta"] == d]
            if not c:
                continue
            dd = [r for r in c if r["dec"]]
            k, n = sum(r["win"] for r in dd), len(dd)
            eps = c[0]["eps"]
            da = eps * np.sqrt(N)
            lo, hi = cp(k, n) if n else (np.nan, np.nan)
            out(f"{N:<6d} {d:<5.1f} {eps:<7.4f} {len(c):<4d} {n:<7d} {len(c)-n:<5d} {k:<5d} "
                f"{(k/n if n else float('nan')):<7.3f} [{lo:.3f},{hi:.3f}]  {norm.cdf(da):<7.3f} {(1+eps)/2:<7.3f} "
                f"{np.mean([r['lean'] for r in c]):<10.3f} {np.median([r['t'] for r in c]):.2f}")
            cells.append((N, d, k, n, da, eps))
    res = {}
    # Delta=0 control
    z = [(N, k, n) for N, d, k, n, da, e in cells if d == 0 and n]
    for N, k, n in z:
        lo, hi = cp(k, n)
        out(f"P1 Delta=0 N={N}: majority-label wins {k}/{n}, CI [{lo:.3f},{hi:.3f}] -> {'in' if lo <= .5 <= hi else 'OUT'}")
    h0 = sum(r["win"] == 0 for r in R if r["dec"] and r["delta"] == 0 and r["maj"] == 0) + \
        sum(r["win"] == 1 for r in R if r["dec"] and r["delta"] == 0 and r["maj"] == 1)
    n0 = sum(1 for r in R if r["dec"] and r["delta"] == 0)
    if n0:
        lo, hi = cp(h0, n0)
        out(f"P1 Delta=0 pooled: hand 0 wins {h0}/{n0}, CI [{lo:.3f},{hi:.3f}]")
    pos = [(k, n, da) for N, d, k, n, da, e in cells if d > 0 and n]
    if pos:
        ah, (alo, ahi), llh = fit_a(pos)
        ll1 = sum(binom_ll(k, n, norm.cdf(da)) for k, n, da in pos)
        lld = sum(binom_ll(k, n, (1 + e) / 2) for N, d, k, n, da, e in cells if d > 0 and n)
        out(f"FIT {VARIANTS[var][0]}: a_hat={ah:.3f} 95%CI[{alo:.3f},{ahi:.3f}]  logL(a)={llh:.2f} "
            f"logL(Phi(D))={ll1:.2f} logL(drift)={lld:.2f}")
        lr1 = 2 * (llh - ll1)
        out(f"    LRT a=1: chi2={lr1:.2f} p={chi2.sf(max(lr1,0),1):.4f};  log10 L(a)/L(drift)={(llh-lld)/np.log(10):.1f};  "
            f"log10 L(a)/L(Phi)={(llh-ll1)/np.log(10):.2f};  log10 L(Phi)/L(drift)={(ll1-lld)/np.log(10):.1f}")
        res.update(a=ah, alo=alo, ahi=ahi, ll_a=llh, ll_phi=ll1, ll_drift=lld)
        perN, llN = {}, 0.0
        for N in Ns:
            cN = [(k, n, da) for M, d, k, n, da, e in cells if M == N and d > 0 and n]
            if cN:
                a_, (l_, h_), ll_ = fit_a(cN)
                perN[N] = (a_, l_, h_)
                llN += ll_
        if len(perN) > 1:
            lr = 2 * (llN - llh)
            out("    per-N a: " + ", ".join(f"N={N}: {a_:.2f} [{l_:.2f},{h_:.2f}]" for N, (a_, l_, h_) in perN.items()))
            out(f"    collapse LRT (H_a vs H_aN): chi2={lr:.2f} df={len(perN)-1} p={chi2.sf(max(lr,0), len(perN)-1):.4f}")
            res.update(perN=perN, collapse_p=float(chi2.sf(max(lr, 0), len(perN) - 1)))
    for N in Ns:
        tt = [r["t"] for r in R if r["N"] == N and r["dec"]]
        if tt:
            out(f"    N={N}: median t_decide={np.median(tt):.2f} (n={len(tt)}), undecided={sum(1 for r in R if r['N']==N and not r['dec'])}/{sum(1 for r in R if r['N']==N)}")
    return res


def main():
    smoke = "--smoke" in sys.argv
    extra = "--extra" in sys.argv   # POST-HOC, added after the main run; appends to the log
    lines = []

    def out(s=""):
        print(s, flush=True)
        lines.append(s)

    grid_b = {200: 3} if smoke else GRID_BASE
    grid_s = {} if smoke else GRID_SENS
    sets = ((3, GRID_SENS), (4, GRID_SENS)) if extra else ((0, grid_b), (1, grid_s), (2, grid_s))
    jobs = [(N, d, rep, v) for v, g in sets for N, reps in g.items()
            for d in DELTAS for rep in range(reps)]
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    out(f"viedma_sim run: {len(jobs)} runs, master seed {MASTER}, 2 workers x 1 thread, started {time.strftime('%Y-%m-%d %H:%M:%S')}")
    with Pool(2) as pool:
        rows = list(pool.imap_unordered(_job, jobs, chunksize=2))
    out(f"wall time {time.time()-t0:.0f} s")
    for v in ((3, 4) if extra else (0, 1, 2)):
        analyze(rows, v, out)
    out("\nPER-RUN TABLE: var N delta rep eps maj win(-1=undecided) lean fm t ncr")
    for r in sorted(rows, key=lambda r: (r["var"], r["N"], r["delta"], r["rep"])):
        out(f"{r['var']} {r['N']} {r['delta']} {r['rep']} {r['eps']:.5f} {r['maj']} {r['win']} {r['lean']} {r['fm']:.4f} {r['t']:.2f} {r['ncr']}")
    if not smoke:
        with open(__file__.replace(".py", ".log"), "a" if extra else "w") as f:
            f.write(("\nPOST-HOC EXTRA (--extra): S sensitivity, not pre-registered\n" if extra else "") + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
