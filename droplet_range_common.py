"""Shared engine for droplet_range_a.py (plateau blind spot) and droplet_range_b.py (volume range).  Simulated drops, not data.

Imports, never edits, droplet_sim2 (curves, FORMS, classify, truths), droplet_power (scorable), exclusion_bak (table, curve).
New here: (1) protocol2, droplet_sim2.protocol with the same random stream when misreads are symmetric, plus a one-way miss,
thinning of mixed drops, an imposed V0 (no pilot) and other level designs; (2) Fitter, the D2 likelihood grids at any slot
vector with the misread rate of each run (policy B of misread_sim.py: e_hat = (k + 0.5)/(n_c + 1)) and, for the amendments
AM1/AM2 only, extra plateau groups; (3) the null cells.  Rules: A = fixed 3% term, c* -1.61; B = measured e_hat, c* -2.21,
NOT SCORABLE above 10% measured misread; AM1 = B + plateau h in the STEP family; AM2 = B + plateau h in both families.

Run (self-check only): UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python droplet_range_common.py
"""
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.special import ndtr

import droplet_sim2 as S
from blind_units import exp_fit
from droplet_power import scorable
from droplet_sim import LADDER
from droplet_sim2 import AG, CG, FORMS, LF, LMU, MU50, Q, SG, V50, XS, classify, curves, protocol

C_A, C_B = -1.61, -2.21                       # droplet_sim2.log, misread_sim.log:5 (c*_B(128)), frozen
NC = 128                                       # control drops
HG = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)
LMUC, CGC = np.arange(-8, 5, 0.1), np.arange(-6, 4, 0.04)     # coarser grids for the plateau groups
NAMES = ("STEP", "SMOOTH", "FLAT", "RISING", "NS")
CONDS8 = [(cv, e) for cv in (0.05, 0.15) for e in (0.0, 0.03, 0.06, 0.10)]
CH = 100
# rule -> (misread used in the fits, fit groups, control cap).  Count groups come first so that ties go to COUNT.
RULES = {"A": ("fixed", ("c0", "s0"), False), "E0": ("zero", ("c0", "s0"), False),
         "B": ("meas", ("c0", "s0"), True), "AM1": ("meas", ("c0", "s0", "sH"), True),
         "AM2": ("meas", ("c0", "cH", "s0", "sH"), True)}


def design(j0, j1, cap):
    """Main levels 0.2 * 8**(j/7) x V0 for j = j0..j1 with cap extension slots each side (design(0, 7, 2) is droplet_sim2's XS)."""
    j = np.arange(j0 - cap, j1 + cap + 1)
    return dict(x=np.log(0.2) + Q * j, main=np.arange(cap, cap + j1 - j0 + 1), cap=cap)


DESIGNS = {"base": design(0, 7, 2), "ext4": design(0, 7, 4), "K10": design(0, 9, 2), "K12": design(-2, 9, 2)}


def fast_exp_fit(V, k, n):
    """blind_units.exp_fit's likelihood (P = exp(-V/V0), k single of n drops per level), solved by Newton on ln(1/V0): it is
    concave, so the maximum is unique and equals exp_fit's up to its tolerance.  No interior maximum or no convergence: exp_fit."""
    if k.sum() == 0 or (n - k).sum() == 0:
        return exp_fit(V, k / n, n)
    t = -np.log(V.mean())
    for _ in range(80):
        e = np.expm1(np.exp(t) * V)
        g = -(k * V).sum() + ((n - k) * V / e).sum()
        d = -((n - k) * V * V * (e + 1) / e ** 2).sum()
        step = np.clip(-g / d / np.exp(t), -1, 1)
        t += step
        if abs(step) < 1e-10:
            return float(np.exp(-t))
    return exp_fit(V, k / n, n)


def protocol2(rng, P, p0, eps, cv, dz=DESIGNS["base"], V0=None, mode="xor", excl=0.0, n=25, exact=False):
    """droplet_sim2.protocol, generalised.  mode 'xor': a read is wrong with probability eps either way (the same random
    stream as droplet_sim2 when excl = 0); 'miss': a mixed drop is read single with probability eps, never the reverse.
    excl: a mixed drop is left out of the count with that probability.  V0 given: no pilot.  Returns (V0, out[2, slots])."""
    def draw(vl, m):
        lv = np.log(vl) + cv * rng.standard_normal(m)
        xtal = rng.random(m) >= np.interp(lv, LF, p0)
        s = rng.random(m) < np.interp(lv, LF, P)
        w = rng.random(m) < eps
        obs = s ^ w if mode == "xor" else s | w
        if excl:
            xtal &= s | (rng.random(m) >= excl)
        return obs[xtal].sum(), xtal.sum()

    if V0 is None:
        k, nn = np.array([draw(v, 20) for v in LADDER], float).T
        V0 = (exp_fit(LADDER[nn > 0], k[nn > 0] / nn[nn > 0], nn[nn > 0]) if exact
              else fast_exp_fit(LADDER[nn > 0], k[nn > 0], nn[nn > 0]))
    kn = {j: draw(V0 * np.exp(dz["x"][j]), n) for j in dz["main"]}
    for _ in range(dz["cap"]):                # extension rule
        k, nn = np.array(list(kn.values()), float).T
        Pl = k[nn > 0] / nn[nn > 0]
        if not Pl.size:
            break
        j = max(kn) + 1 if Pl.min() > 0.4 else min(kn) - 1 if Pl.max() < 0.6 else None
        if j is None:
            break
        kn[j] = draw(V0 * np.exp(dz["x"][j]), n)
    out = np.zeros((2, len(dz["x"])))
    for j, (a, b) in kn.items():
        out[:, j] = a, b
    return V0, out


def controls(rng, reps, eps, mode):
    """Misread count among NC control drops of known outcome, half single-seed and half two-seed."""
    if mode == "xor":
        return (rng.random((reps, NC)) < eps).sum(1)
    return rng.binomial(NC // 2, eps, reps)    # one-way miss: only the two-seed controls can be read wrong


def scorable_v(K, N):
    """droplet_power.scorable for every row."""
    P = K / np.maximum(N, 1)
    hi, lo = np.where(N > 0, P, -np.inf).max(1), np.where(N > 0, P, np.inf).min(1)
    return ((N - K).sum(1) >= 10) & (K.sum(1) >= 10) & (hi >= 0.6) & (lo <= 0.4)


def trend_z(K, N, x):
    xb = (N * x).sum(1, keepdims=True) / N.sum(1, keepdims=True)
    pb = K.sum(1) / N.sum(1)
    return (K * (x - xb)).sum(1) / np.sqrt(pb * (1 - pb) * (N * (x - xb) ** 2).sum(1) + 1e-300)


class Fitter:
    """Log-likelihood grids at ln V/V0 = x.  c0, s0: droplet_sim2's COUNT and STEP grids.  cH, sH: the same with a plateau
    P = h + (1 - h) P0, h in HG, on the coarser grids LMUC, CGC (used by the amendments only)."""

    def __init__(self, x, plateau=False):
        self.x = x
        cur = curves(x)
        self.P = {"c0": cur["count"][0], "s0": cur["step"][0]}
        if plateau:
            Fc = np.concatenate([FORMS[f](np.exp(LMUC[:, None] + a * x)) for f in FORMS for a in AG])
            Fs = np.concatenate([ndtr(-(x - CGC[:, None]) / s) for s in SG])
            self.P["cH"] = np.concatenate([h + (1 - h) * Fc for h in HG])
            self.P["sH"] = np.concatenate([h + (1 - h) * Fs for h in HG])

    def W(self, e, g):
        """(2 slots, columns) so that [K, N] @ W = K log P_obs + (N - K) log(1 - P_obs), P_obs = e + (1 - 2e) P."""
        obs = np.clip(e + (1 - 2 * e) * self.P[g], 1e-12, 1 - 1e-12)
        lp, l1 = np.log(obs).T, np.log1p(-obs).T
        return np.vstack([lp - l1, l1]).astype(np.float64 if g in ("c0", "s0") else np.float32)

    def fit(self, K, N, eh, groups):
        """Per group: max log-likelihood and range of the best-fit curve over the levels with drops; run r uses misread eh[r]."""
        X = np.hstack([K, N])
        out = {g: (np.empty(len(K)), np.empty(len(K))) for g in groups}
        for e in np.unique(eh):
            idx = np.flatnonzero(eh == e)
            for g in groups:
                W = self.W(e, g)
                Xg = X[idx].astype(W.dtype)
                for i in range(0, len(idx), CH):
                    s = idx[i:i + CH]
                    m = Xg[i:i + CH] @ W
                    a = m.argmax(1)
                    Pb, ok = self.P[g][a], N[s] > 0
                    out[g][0][s] = m[np.arange(len(s)), a]
                    out[g][1][s] = np.where(ok, Pb, -np.inf).max(1) - np.where(ok, Pb, np.inf).min(1)
        return out


_FIT = {}


def fitter(name, plateau):
    if (name, plateau) not in _FIT:
        _FIT[(name, plateau)] = Fitter(DESIGNS[name]["x"], plateau)
    return _FIT[(name, plateau)]


def rule_stats(fits, groups):
    """L = best COUNT ll - best STEP ll, and the range of the best-fit curve, among the rule's groups."""
    ll = np.stack([fits[g][0] for g in groups], 1)
    rg = np.stack([fits[g][1] for g in groups], 1)
    c = [i for i, g in enumerate(groups) if g[0] == "c"]
    s = [i for i, g in enumerate(groups) if g[0] == "s"]
    best = ll.argmax(1)
    return dict(L=ll[:, c].max(1) - ll[:, s].max(1), rng=rg[np.arange(len(ll)), best])


def run_batch(task):
    """Simulate and score a list of cells.  cell = (name, P, p0, cv, eps, mode, excl, seed).  task keys: cells, reps, rules,
    design, V0 (None = pilot), cstar (rule -> c*, or None to calibrate), want_v0, exact (pilot fitted by blind_units.exp_fit).  Returns one dict per cell:
    counts[rule] (5 verdict counts) and drops, V0; or, when calibrating, calib[rule] = (runs reaching L, 5th percentile of L)."""
    cells, reps, rules = task["cells"], task["reps"], task["rules"]
    dz = DESIGNS[task["design"]]
    F = fitter(task["design"], any(r in ("AM1", "AM2") for r in rules))
    Ks, Ns, kc, V0s = [], [], [], []
    for name, P, p0, cv, eps, mode, excl, seed in cells:
        rng = np.random.default_rng(seed)
        V0, KN = zip(*(protocol2(rng, P, p0, eps, cv, dz, task["V0"], mode, excl, exact=task.get("exact", False)) for _ in range(reps)))
        KN = np.array(KN)
        Ks.append(KN[:, 0]); Ns.append(KN[:, 1]); V0s.append(np.array(V0))
        kc.append(controls(rng, reps, eps, mode))
    K, N, kc = np.vstack(Ks), np.vstack(Ns), np.concatenate(kc)
    z, sc = trend_z(K, N, F.x), scorable_v(K, N)
    eh = {"meas": (kc + 0.5) / (NC + 1), "fixed": np.full(len(K), S.EPS_FIT), "zero": np.zeros(len(K))}
    st = {}
    for kind in {RULES[r][0] for r in rules}:
        use = [r for r in rules if RULES[r][0] == kind]
        groups = [g for g in ("c0", "cH", "s0", "sH") if any(g in RULES[r][1] for r in use)]
        fits = F.fit(K, N, eh[kind], groups)
        for r in use:
            st[r] = rule_stats(fits, RULES[r][1])
    res = []
    for i, cell in enumerate(cells):
        sl = slice(i * reps, (i + 1) * reps)
        o = {"drops": float(N[sl].sum(1).mean()), "V0": V0s[i] if task.get("want_v0") else None}
        if task["cstar"] is None:
            o["calib"] = {}
            for r in rules:
                ok = sc[sl] & (st[r]["rng"][sl] >= 0.4) & (z[sl] <= 1.645)
                o["calib"][r] = (int(ok.sum()), float(np.percentile(st[r]["L"][sl][ok], 5)) if ok.sum() >= 50 else np.nan)
        else:
            o["counts"] = {}
            for r in rules:
                ns = ~sc[sl] | ((kc[sl] / NC > 0.10) if RULES[r][2] else False)
                lab = classify(dict(L=st[r]["L"][sl], rng=st[r]["rng"][sl], z=z[sl], ns=ns), task["cstar"][r])
                o["counts"][r] = np.array([(lab == v).sum() for v in NAMES])
        res.append(o)
    return res


def pmap(tasks, procs=4, label=""):
    """Run batches on at most `procs` cores, in order, with a progress line per batch on stderr."""
    import sys
    out = [None] * len(tasks)
    with Pool(procs) as pool:
        for n, (i, r) in enumerate(pool.imap_unordered(_tag, list(enumerate(tasks)))):
            out[i] = r
            print(f"  {label} batch {n + 1}/{len(tasks)}", file=sys.stderr, flush=True)
    return out


def _tag(it):
    return it[0], run_batch(it[1])


def batches(cells, size):
    return [cells[i:i + size] for i in range(0, len(cells), size)]


def _cache(name, make):
    d = os.environ.get("DR_CACHE")                # optional scratch cache; results do not depend on it
    f = Path(d) / f"{name}.npz" if d else None
    if f and f.exists():
        z = np.load(f, allow_pickle=False)
        return [(str(n), P, p0) for n, P, p0 in zip(z["names"], z["P"], z["p0"])]
    cells = make()
    if f:
        f.parent.mkdir(parents=True, exist_ok=True)
        np.savez(f, names=np.array([c[0] for c in cells]), P=np.array([c[1] for c in cells]), p0=np.array([c[2] for c in cells]))
    return cells


def null_cells(procs=4):
    """The 87 counting and drying cells of droplet_sim2.truths() (15 counting forms, 72 run-1 drying models)."""
    def make():
        S.Pool = lambda: Pool(procs)           # truths() opens Pool() on every core; keep to `procs`
        T, _ = S.truths()
        return [(n, P, p0) for n, k, P, p0 in T if k == "null"]
    return _cache("null", make)


def floor_cells():
    """Counting forms with a seeded or dominant-site floor: P = h + (1 - h) F(mu), 3 forms x a in 2/3, 1, 5/3 x h in .1, .2, .3."""
    lv, out = LF - np.log(V50), []
    for f in FORMS:
        for a in (2 / 3, 1, 5 / 3):
            mu = MU50[f] * np.exp(a * lv)
            for h in (0.1, 0.2, 0.3):
                out.append((f"{f} a={a:.2f} floor {h}", h + (1 - h) * FORMS[f](mu), np.exp(-mu) if f == "all" else 0 * mu))
    return out


def exclusion_cells(procs=4):
    """exclusion_bak.py's Monte Carlo curves with eta ~ V^e, e = 1/6, 1/3, 1/2 (the physical ones), nu 0/0.3/1, three nucleus sites."""
    def make():
        import exclusion_bak as X
        cfgs = [(i, s) for i in range(len(X.NU)) for s in range(len(X.SITES))]
        with Pool(procs) as pool:
            tabs = pool.map(X.table, cfgs)
        out = []
        for (i, s), (T, _) in zip(cfgs, tabs):
            for e in (1 / 6, 1 / 3, 1 / 2):
                P = X.curve(T, e)
                if P is not None:
                    out.append((f"excl {X.SITES[s]} nu {X.NU[i]} e={e:.2f}", P, 0 * P))
        return out
    return _cache("excl", make)


def probit(h, s, c=1.0):
    """Truth on LF: c (h + (1 - h) Phi(-ln(V/V50)/s))."""
    return c * (h + (1 - h) * ndtr(-(LF - np.log(V50)) / s))


def self_check():
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    assert np.allclose(DESIGNS["base"]["x"], XS) and DESIGNS["base"]["main"].tolist() == list(range(2, 10))
    P = FORMS["first"](MU50["first"] * np.exp(LF - np.log(V50)))
    for seed in range(3):                      # same random stream as droplet_sim2.protocol for symmetric misreads
        a, b = protocol(np.random.default_rng(seed), P, 0 * P, 0.03, 0.05), protocol2(np.random.default_rng(seed), P, 0 * P, 0.03, 0.05, exact=True)
        assert a[0] == b[0] and (a[1] == b[1]).all()
    r = np.random.default_rng(5)               # Newton pilot fit = exp_fit, on 40 random pilots (counting and step shaped)
    for i in range(40):
        n20 = np.full(7, 20.0)
        k = r.binomial(20, np.exp(-LADDER / r.uniform(0.5, 8)) if i % 2 else (LADDER < r.uniform(0.5, 8)) * 0.94 + 0.03).astype(float)
        if 0 < k.sum() < n20.sum():
            assert abs(fast_exp_fit(LADDER, k, n20) / exp_fit(LADDER, k / n20, n20) - 1) < 1e-3
    rng = np.random.default_rng(1)
    zero = 0 * LF
    _, o = protocol2(rng, zero, zero, 0.2, 0.0, V0=3.0, mode="miss", n=2000)      # all mixed: reads single only by a miss, share 0.2
    assert abs(o[0].sum() / o[1].sum() - 0.2) < 0.02 and abs(controls(rng, 4000, 0.2, "miss").mean() - 12.8) < 0.3
    _, o = protocol2(rng, 0.5 + zero, zero, 0.0, 0.0, V0=3.0, excl=0.8, n=4000)    # 0.5 / (0.5 + 0.5 * 0.2) = 5/6 read single
    assert abs(o[0].sum() / o[1].sum() - 5 / 6) < 0.02
    K = np.random.default_rng(0).binomial(25, np.linspace(0.95, 0.05, 8), (6, 8)).astype(float)
    K, N = np.pad(K, ((0, 0), (2, 2))), np.pad(np.full((6, 8), 25.0), ((0, 0), (2, 2)))
    F = Fitter(XS, True)
    ref = S.stats(K, N, XS)                    # the unchanged droplet_sim2 statistics at the fixed 3% term
    f = F.fit(K, N, np.full(6, 0.03), ("c0", "s0"))
    r = rule_stats(f, ("c0", "s0"))
    assert np.allclose(r["L"], ref["L"], atol=1e-8) and np.allclose(r["rng"], ref["rng"], atol=1e-8)
    assert np.allclose(trend_z(K, N, XS), ref["z"]) and (scorable_v(K, N) == [scorable(k[n > 0], n[n > 0]) for k, n in zip(K, N)]).all()
    x = XS[2:10]                               # data from a plateau-0.3 step: the plateau group must win by a wide margin
    Pt = 0.3 + 0.7 * ndtr(-x / 0.1)
    Kp, Np = np.pad(np.round(Pt * 1e4)[None], ((0, 0), (2, 2))), np.pad(np.full((1, 8), 1e4), ((0, 0), (2, 2)))
    f = F.fit(Kp, Np, np.full(1, 0.03), ("c0", "s0", "cH", "sH"))
    assert f["sH"][0][0] > f["s0"][0][0] + 50 and f["sH"][0][0] >= f["cH"][0][0] and F.P["sH"].min() < 0.06 and F.P["cH"].shape[0] == 3 * len(AG) * len(LMUC) * len(HG)


if __name__ == "__main__":
    self_check()
    print("droplet_range_common self-check passed")
