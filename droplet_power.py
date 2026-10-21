"""Power for the NaClO3 droplet test (journals/droplet_prereg.md).

How many droplets per volume level does the free-width fit P = 1/(1 + (V/Vk)^h) need to tell the counting rule
(every counting form has h <= 1.64) from a sharp critical size (a step), under the same h <= 3 rule as blind_jang.py J2?
Volumes are in units of the pilot's V0, 8 levels log-spaced over 0.2-1.6.  The true V0 may be 0.5-2 times the pilot's
(r), the pipette scatters volumes (CV sigma), and the reader misreads a droplet's outcome with probability eps.
Unscorable runs extend by up to two levels per side (the pre-registered extension rule).

The fit carries a fixed 3% misread term, P_obs = 0.03 + 0.94 P (fit_h).  Without it, 3% misreads make a step near the
top of the range look broad in about 1 run in 6; the first block of output shows this.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python droplet_power.py > droplet_power.log
"""
import numpy as np

from blind_lab import best, loglik
from blind_units import exp_fit

LN2 = np.log(2)
MU_HALF = 1.2564312086261697   # mu/(e^mu - 1) = 1/2: all-nuclei-at-once form has P = 1/2 where exp(-x) does
LAPSE = 0.03
K = 8
REPS = 400


def p_true(kind, x):
    """P_single at x = V / V0_true; every form has P = 1/2 at x = ln 2."""
    if kind == "count":                       # first nucleus + Poisson extras, P = exp(-mu/2); h equiv 1.40
        return np.exp(-x)
    if kind == "all":                         # N ~ Poisson(mu) given N >= 1; h equiv 1.64
        mu = MU_HALF * x / LN2
        return mu / np.expm1(mu)
    if kind == "step":                        # sharp critical size
        return (x < LN2).astype(float)
    return 1 / (1 + (x / LN2) ** kind)        # logistic of given h


def fit_h(v, P, n, lapse=LAPSE):
    """P_obs = lapse + (1 - 2 lapse) / (1 + (V/Vk)^h), binomial likelihood, h in [0.1, 50].  Returns (Vk, h)."""
    lv = np.log10(v)
    curve = lambda x: lapse + (1 - 2 * lapse) / (1 + 10 ** np.clip(np.clip(10 ** x[1], 0.1, 50) * (lv - x[0]), -50, 50))
    x, _ = best(lambda x: -loglik(P * n, n, curve(x)), [(np.median(lv), np.log10(h)) for h in (0.5, 1.5, 5.0, 20.0)])
    return 10 ** x[0], float(np.clip(10 ** x[1], 0.1, 50))


def scorable(k, n):
    P = k / n
    return (n - k).sum() >= 10 and k.sum() >= 10 and P.max() >= 0.6 and P.min() <= 0.4


def trial(rng, kind, n, r, sigma, eps, lapse):
    """One main run.  Returns fitted h, or None if unscorable after the extension rule."""
    v = list(np.geomspace(0.2, 1.6, K))
    q = v[1] / v[0]

    def draw(vi):
        x = vi / r * np.exp(sigma * rng.standard_normal(n))
        return int(((rng.random(n) < p_true(kind, x)) ^ (rng.random(n) < eps)).sum())

    k = [draw(vi) for vi in v]
    for _ in range(2):                        # extension rule: up to two extra levels per side, one at a time
        P = np.array(k) / n
        if P.min() > 0.4:
            v.append(v[-1] * q)
            k.append(draw(v[-1]))
        elif P.max() < 0.6:
            v.insert(0, v[0] / q)
            k.insert(0, draw(v[0]))
    k, nn = np.array(k, float), np.full(len(v), float(n))
    return fit_h(np.array(v), k / nn, nn, lapse)[1] if scorable(k, nn) else None


def correct(kind, n, r, sigma, eps, lapse=LAPSE, seed=1):
    """Share of runs that are scorable and called right: h <= 3 for counting forms, h > 3 otherwise."""
    rng = np.random.default_rng(seed)
    hs = [trial(rng, kind, n, r, sigma, eps, lapse) for _ in range(REPS)]
    want = kind in ("count", "all")
    return np.mean([h is not None and (h <= 3) == want for h in hs]), [h for h in hs if h is not None]


def self_check():
    v = np.geomspace(0.2, 1.6, K)
    nn = np.full(K, 1e5)
    obs = lambda kind: LAPSE + (1 - 2 * LAPSE) * p_true(kind, v)
    assert 1.2 < fit_h(v, obs("count"), nn)[1] < 1.6
    assert 1.4 < fit_h(v, obs("all"), nn)[1] < 1.9
    assert fit_h(v, obs("step"), nn)[1] > 30
    assert abs(p_true("all", np.array([LN2]))[0] - 0.5) < 1e-9


def pilot():
    """Doubling ladder 0.25-16 uL, 20 drops per level; spread of fitted V0 / true V0."""
    v = 0.25 * 2.0 ** np.arange(7)
    n = np.full(7, 20.0)
    for V0 in (1.0, 4.0, 10.0):
        rng = np.random.default_rng(3)
        rat = []
        for _ in range(REPS):
            k = (rng.random((7, 20)) < np.exp(-v[:, None] * np.exp(0.05 * rng.standard_normal((7, 20))) / V0)).sum(1)
            rat.append(exp_fit(v, k / n, n) / V0)
        lo, hi = np.percentile(rat, [2.5, 97.5])
        print(f"  true V0 {V0:>4} uL: fitted / true, 95% range {lo:.2f}-{hi:.2f}")


def main():
    self_check()
    print("Pilot, counting rule, 7 levels 0.25-16 uL doubling, 20 drops each (140 drops):")
    pilot()

    print(f"\nWhy the fit has a misread term: a step at r = 2, pipette CV 0.05, {K} x 20 droplets, share called right")
    for eps in (0.0, 0.03, 0.06):
        print(f"  misread {eps:.2f}: fit without the term {correct('step', 20, 2.0, 0.05, eps, 0.0)[0]:.2f}, "
              f"with it {correct('step', 20, 2.0, 0.05, eps)[0]:.2f}")

    print(f"\nMain run, {K} levels 0.2-1.6 x V0_pilot, {REPS} simulated runs per cell.  Each number is the share of runs")
    print("that are scorable and called right (counting: h <= 3; step and h = 6: h > 3).  r = true V0 / pilot V0.")
    truths = [(k, r) for k in ("count", "all", "step") for r in (0.5, 0.7, 1.0, 1.4, 2.0)] + [(6, 1.0)]
    worst = {}
    for sigma, eps in ((0.05, 0.0), (0.05, 0.03), (0.10, 0.06)):
        print(f"\n pipette CV {sigma:.2f}, misread rate {eps:.2f}")
        for n in (5, 10, 15, 20, 30):
            cells = [(f"{kind}/r{r:g}", correct(kind, n, r, sigma, eps)[0]) for kind, r in truths]
            worst[n] = min(worst.get(n, 1.0), min(c for _, c in cells))
            print(f"  n {n:>2} per level = {K * n:>3}: " + " ".join(f"{name} {c:.2f}" for name, c in cells))
    print("\nWorst cell over every truth and condition: " + ", ".join(f"n {n}: {w:.2f}" for n, w in worst.items()))

    print(f"\nReported, not predicted: central 95% of fitted h under the counting forms, n 20, CV 0.05, misread 0.03, r 0.5-2:")
    for kind in ("count", "all"):
        hs = np.concatenate([correct(kind, 20, r, 0.05, 0.03, seed=11)[1] for r in (0.5, 1.0, 2.0)])
        print(f"  {kind}: h 2.5% {np.percentile(hs, 2.5):.2f}, median {np.median(hs):.2f}, 97.5% {np.percentile(hs, 97.5):.2f}")


if __name__ == "__main__":
    main()
