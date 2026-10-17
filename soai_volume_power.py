"""Power for EXPERIMENT_SOAI.md: how many Soai runs tell m* ∝ sqrt(V) from m* flat in V?

The counting rule with a seed of m net molecules among N_d deciding molecules is
    P(seed hand wins) = Phi(m / sqrt(N_d)),   Delta = m / sqrt(N_d),
so the seed that gives Delta = 1 (P = 0.84) is m* = sqrt(N_d).  Two predictions, fixed before any data:
  H1 one well-mixed pool:      N_d = c* V, so m* ∝ V^0.5  (slope beta = 0.5 on log-log axes)
  H0 local or early events:    N_d set by a few early units, so m* flat (beta = 0)
Anchor (soai_check.py S4): ~60,000 seed molecules controlled Singleton & Vo 2003 runs, so N_d <~ 3.6e9 and
m* ~ 6e4 at their scale.  ASSUMPTION: their scale is the middle volume here (1 mL); their flask volume is not in the
open text.  Stage 1 of the protocol measures m* at 1 mL before the ladders are fixed, so the anchor only sets the
starting bracket; the robustness runs below ask what a residual error of 2x in that calibration costs.

Analysis (pre-registered): one probit fit over all seeded runs, P = Phi(m exp(-theta0 - beta x)), x = ln(V / V_REF),
maximised on a grid; likelihood-ratio tests of beta = 0 and of beta = 0.5, each at alpha = 0.05 (chi2_1 > 3.84).
Mirror seeds (half S, half R at every level) make a stray bias b orthogonal to (theta0, beta) at b = 0, so b is left
out of the power calculation.  Unseeded runs carry no information on beta and are not counted here.
"""
import numpy as np
from scipy.stats import binom, chi2, norm

S_REF = np.sqrt(3.6e9)  # m* at V_REF, soai_check.py S4
V_REF = 1.0  # mL, assumed Singleton & Vo scale (see docstring)
VOLS = np.array([0.1, 1.0, 10.0])  # mL, 100x span
LEVELS, LO, HI = 4, 0.7, 2.8  # seed ladder: LEVELS geometric steps from LO * lower to HI * higher prediction
ALPHA, POWER = 0.05, 0.80
SIMS = 4000
CRIT = chi2.ppf(1 - ALPHA, 1)
X = np.log(VOLS / V_REF)

out = []


def say(s=""):
    print(s)
    out.append(s)


def m_star(beta, anchor=1.0):
    return anchor * S_REF * (VOLS / V_REF) ** beta


def ladder():
    """Seed counts (molecules) per volume, shape (volumes, LEVELS); fixed before data, covers both predictions."""
    lo = np.minimum(m_star(0.0), m_star(0.5)) * LO
    hi = np.maximum(m_star(0.0), m_star(0.5)) * HI
    return np.exp(np.linspace(np.log(lo), np.log(hi), LEVELS).T)


M = ladder()

# Likelihood grid over (theta0, beta); theta0 relative to ln S_REF.  beta = 0 and 0.5 sit exactly on the grid.
T_GRID = np.log(S_REF) + np.arange(-4, 4.0001, 0.05)
B_GRID = np.round(np.arange(-1, 2.0001, 0.01), 2)
B0, B1 = int(np.where(B_GRID == 0.0)[0][0]), int(np.where(B_GRID == 0.5)[0][0])
TT, BB = np.meshgrid(T_GRID, B_GRID, indexing="ij")
Z = M.reshape(1, 1, -1) * np.exp(-TT[..., None] - BB[..., None] * np.repeat(X, LEVELS))  # (T, B, cells)
LP, LQ = norm.logcdf(Z).reshape(-1, Z.shape[-1]), norm.logcdf(-Z).reshape(-1, Z.shape[-1])


def loglik(k, n):
    """k: seed-hand wins per cell, shape (sims, cells); returns log-lik on the grid, shape (sims, T, B)."""
    return (k @ LP.T + (n - k) @ LQ.T).reshape(len(k), len(T_GRID), len(B_GRID))


def lr_tests(k, n):
    """LR statistics for beta = 0 and beta = 0.5, and the grid MLE of beta."""
    ll = loglik(k, n)
    top = ll.max(axis=(1, 2))
    beta_hat = B_GRID[ll.max(axis=1).argmax(axis=1)]
    return 2 * (top - ll[:, :, B0].max(axis=1)), 2 * (top - ll[:, :, B1].max(axis=1)), beta_hat


def simulate(beta, n, anchor=1.0, sims=SIMS, seed=0):
    p = norm.cdf(M / m_star(beta, anchor)[:, None]).ravel()
    k = np.random.default_rng(seed).binomial(n, p, size=(sims, p.size)).astype(float)
    lr0, lr1, _ = np.concatenate([np.array(lr_tests(c, n)) for c in np.array_split(k, 20)], axis=1)
    return lr0, lr1


def power(n, anchor=1.0):
    """(P reject flat | sqrt V true, P reject sqrt V | flat true, P reject flat | flat true)."""
    h1_lr0, _ = simulate(0.5, n, anchor, seed=1)
    h0_lr0, h0_lr1 = simulate(0.0, n, anchor, seed=2)
    return (h1_lr0 > CRIT).mean(), (h0_lr1 > CRIT).mean(), (h0_lr0 > CRIT).mean()


def info(z, n):
    """Fisher information on ln m* from n runs at Delta = z."""
    return n * norm.pdf(z) ** 2 * z**2 / (norm.cdf(z) * norm.sf(z))


def fisher_se(beta, n):
    """Asymptotic SE of beta-hat over all volumes."""
    z = (M / m_star(beta)[:, None]).ravel()
    g = np.stack([np.ones_like(z), np.repeat(X, LEVELS)])
    return np.sqrt(np.linalg.inv((g * info(z, n)) @ g.T)[1, 1])


def se_lnm(j, beta, n):
    """Asymptotic SE of ln m* from one volume's ladder alone (Stage 1b, Arms M and S)."""
    return 1 / np.sqrt(info(M[j] / m_star(beta)[j], n).sum())


say(f"anchor: N_d = {S_REF**2:.1e} at {V_REF:g} mL, m* = {S_REF:.0f} molecules (soai_check.py S4)")
say("predicted m* (molecules) and seed ladder, fixed before data:")
for j, v in enumerate(VOLS):
    say(f"  {v:5.1f} mL  flat {m_star(0)[j]:8.2e}  sqrt V {m_star(0.5)[j]:8.2e}  ladder "
        + " ".join(f"{m:.2e}" for m in M[j]) + f"  ({M[j,0]/6.022e2:.2f}-{M[j,-1]/6.022e2:.1f} zmol)")

target = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)
n_fisher = next(n for n in range(2, 200, 2) if max(fisher_se(0.5, n), fisher_se(0.0, n)) <= 0.5 / target)
say(f"\nasymptotic (Fisher) runs per seed level: {n_fisher}  (SE(beta) <= 0.5/{target:.2f} = {0.5/target:.3f})")

say("\nMonte Carlo, LR tests at alpha = 0.05, %d sims each (runs per seed level, split half S / half R):" % SIMS)
say("   n  runs/vol  P(reject flat | sqrtV)  P(reject sqrtV | flat)  P(reject flat | flat)")
for n in range(2, 61, 2):
    pw1, pw0, size = power(n)
    say(f"  {n:2d}  {n*LEVELS:8d}  {pw1:22.3f}  {pw0:22.3f}  {size:21.3f}")
    if min(pw1, pw0) >= POWER:
        break
N_RUN = n

say(f"\nchosen: {N_RUN} runs per seed level ({N_RUN//2} S-seeded + {N_RUN//2} R-seeded), "
    f"{N_RUN*LEVELS} seeded runs per volume, {N_RUN*LEVELS*len(VOLS)} in all")
rob = {}
for f in (0.5, 2.0):
    rob[f] = power(N_RUN, f)[:2]
    say(f"  if true m* at {V_REF:g} mL is {f:g}x the calibrated value: power {rob[f][0]:.3f} (sqrt V), "
        f"{rob[f][1]:.3f} (flat)")

# One-ladder fits: Stage 1b calibration at 1 mL, and the side arms (M at 10 mL, S at 1 mL) compared with Arm V.
z95, z_det = norm.ppf(0.975), target * np.sqrt(2)  # two independent ladders of equal SE
cal = np.exp(z95 * se_lnm(1, 0.5, N_RUN))
det = {j: np.exp(z_det * max(se_lnm(j, 0.0, N_RUN), se_lnm(j, 0.5, N_RUN))) for j in (1, 2)}
say(f"\none ladder of {N_RUN * LEVELS} runs: Stage 1b fixes m* at {V_REF:g} mL within a factor {cal:.2f} (95%); "
    f"a side arm detects a change in m* of {det[1]:.1f}x at 1 mL, {det[2]:.1f}x at 10 mL (80%, alpha 0.05)")

# Unseeded chance-zone runs: 8 per volume, pooled 24.  Two-sided 50/50 binomial test; power against a Delta = 1 stray bias.
n0 = 8 * len(VOLS)
acc = [k for k in range(n0 + 1) if 2 * min(binom.cdf(k, n0, 0.5), binom.sf(k - 1, n0, 0.5)) > ALPHA]
p_stray = 1 - binom.pmf(acc, n0, norm.cdf(1.0)).sum()
say(f"\nunseeded: {n0} pooled runs, 50/50 accepted for {acc[0]}-{acc[-1]} of one hand; "
    f"P(reject | stray Delta = 1 bias) = {p_stray:.3f}")

with open(__file__.replace(".py", ".log"), "w") as fh:
    fh.write("\n".join(out) + "\n")

# Self-check
assert abs(S_REF - 6.0e4) < 1, "anchor must match soai_check.py S4 (60,000 molecules)"
for j in range(len(VOLS)):
    assert M[j, 0] < min(m_star(0)[j], m_star(0.5)[j]) and M[j, -1] > max(m_star(0)[j], m_star(0.5)[j]), "ladder"
big = 10**6  # expected counts from a huge sample: the grid MLE must return the true slope
for b in (0.0, 0.5):
    k = (big * norm.cdf(M / m_star(b)[:, None]).ravel())[None]
    assert abs(lr_tests(k, big)[2][0] - b) <= 0.01, f"fit does not recover beta = {b}"
pw1, pw0, size = power(N_RUN)
assert min(pw1, pw0) >= POWER, "chosen n must reach 80% power both ways"
assert 0.02 <= size <= 0.09, f"LR test size {size:.3f} should be near 0.05"
assert N_RUN >= n_fisher, "small-sample Monte Carlo should not beat the asymptotic count"
assert cal <= 2.0, "Stage 1b must calibrate within the 2x the robustness runs cover"
assert p_stray > 0.9, "pooled unseeded runs should catch a Delta = 1 stray bias"
print("self-check passed")
