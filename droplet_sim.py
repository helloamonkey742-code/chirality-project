"""Model study for the NaClO3 droplet test (journals/droplet_prereg.md): simulated evaporating drops, not data.

Does the pre-registered protocol (pilot, 8 x 25 main run, extension rule, D2 fit with the 3% misread term, h <= 3)
call the counting rule right when the drops evaporate, as Bak et al. (2023) drops do?  Predictions:
droplet_sim_predictions.txt, sha256 PRED_SHA, written before the first run of the grid.

A drop dries over tau = 0..1 (diffusion-limited: V = V_start (1 - tau)^1.5, drying time ~ V_start^(2/3)), so the
supersaturation is S = (1 - tau)^-1.5.  Primary nuclei form at rate c J(S) x (bulk volume or contact-line length) x
drying time, J = exp(-B / ln^2 S), each of random hand.  Each nucleus leaves a fraction exp(-((tau - tau_k)/w)^p) of
the rate (p 1.5: diffusive depletion, w fixed; p 3: growth front, w ~ V^(-1/3)).  Single = every nucleus shares a
hand, 2^(1 - N) given N >= 1.  Drops with no nucleus by tau = 0.999 have no crystal and are excluded, as in the
prereg ("excl").  Sensitivity, not predicted ("dry"): such drops crystallize once at dryness and count as single.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python droplet_sim.py > droplet_sim.log
"""
import hashlib
from pathlib import Path

import numpy as np

from blind_units import exp_fit
from droplet_power import fit_h, scorable

PRED_SHA = "58a1b610cb577ab94652e1f984805605fc12abf78929f34a473aef1a94cd40a8"
G = 1001
TAU = np.linspace(0, 0.999, G)
LNS = -1.5 * np.log1p(-TAU)                  # ln S
LADDER = 0.25 * 2.0 ** np.arange(7)          # pilot volumes, uL
GRID = np.geomspace(0.03, 300, 65)           # ponytail: curve clamps outside 0.03-300 uL
LG = np.log(GRID)
CAP = 20                                     # 2^-19: a drop with 20 nuclei counts as mixed
REPS = 200
MODELS = [(s, None, None, 0.0) for s in ("line", "bulk")] + [
    (s, p, w, B) for s in ("line", "bulk") for p in (1.5, 3) for w in (0.01, 0.1, 0.5) for B in (0.0, 1.0, 10.0)]


def rate(m, c, v):
    """Nuclei per tau step before quench, drop of v uL."""
    site, p, w1, B = m
    a, g = (5 / 3, (1 - TAU) ** 1.5) if site == "bulk" else (1.0, (1 - TAU) ** 0.5)
    return c * v ** a * g * np.exp(-B / np.maximum(LNS, 1e-12) ** 2) * TAU[1]


def nuclei(m, c, v, n=2000, seed=0):
    """Primary nuclei per drop, drawn one event at a time with the quench of the nuclei before."""
    p, w1 = m[1], m[2]
    lam = rate(m, c, v)
    w = w1 * v ** (-1 / 3) if p == 3 else w1
    rng = np.random.default_rng(seed)
    N, s, ext, act = np.zeros(n, int), np.zeros(n, int), np.zeros((n, G)), np.arange(n)
    while act.size:
        H = np.zeros((act.size, G + 1))
        H[:, 1:] = np.cumsum(lam * np.exp(-ext[act]), 1)
        tgt = H[np.arange(act.size), s[act]] + rng.exponential(size=act.size)
        hit = H[:, -1] >= tgt
        j = (H[hit] >= tgt[hit, None]).argmax(1) - 1      # grid step of the next nucleus
        ev = act[hit]
        N[ev] += 1
        s[ev] = j + 1
        if p is not None:
            ext[ev] += (np.clip(TAU - TAU[j, None], 0, None) / w) ** p
        act = ev[N[ev] < CAP]
    return N


def probs(m, c, v, n=2000, seed=0):
    """(P_single among drops with a crystal, share with no crystal).  No quench: Poisson, exact."""
    if m[1] is None:
        mu = rate(m, c, v).sum()
        return 2 / (np.exp(mu / 2) + 1), np.exp(-mu)
    N = nuclei(m, c, v, n, seed)
    return (np.mean(2.0 ** (1 - N[N > 0])) if (N > 0).any() else 1.0), np.mean(N == 0)


def calibrate(m):
    """c with P_single(3.1 uL) = 1/2 (same random numbers at every c), or None if no c in 1e-6..1e12 gives it."""
    f = lambda lc: probs(m, 10 ** lc, 3.1)[0] - 0.5
    lo, hi = -6.0, 12.0
    if f(hi) > 0 or f(lo) < 0:
        return None
    for _ in range(22):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 10 ** ((lo + hi) / 2)


def protocol(rng, P, p0, eps=0.03, cv=0.05):
    """Prereg Sections 5-6 on simulated drops (P, p0 on GRID).  Returns the D2 fitted h, or None if not scorable."""
    def draw(vl, m):
        lv = np.log(vl) + cv * rng.standard_normal(m)
        xtal = rng.random(m) >= np.interp(lv, LG, p0)
        single = (rng.random(m) < np.interp(lv, LG, P)) ^ (rng.random(m) < eps)
        return single[xtal].sum(), xtal.sum()

    k, n = np.array([draw(v, 20) for v in LADDER], float).T
    V0 = exp_fit(LADDER[n > 0], k[n > 0] / n[n > 0], n[n > 0])
    v = list(V0 * np.geomspace(0.2, 1.6, 8))
    q = v[1] / v[0]
    kn = [draw(x, 25) for x in v]
    for _ in range(2):                        # extension rule
        k, n = np.array(kn, float).T
        Pl = k[n > 0] / n[n > 0]
        if Pl.min() > 0.4:
            v.append(v[-1] * q)
            kn.append(draw(v[-1], 25))
        elif Pl.max() < 0.6:
            v.insert(0, v[0] / q)
            kn.insert(0, draw(v[0], 25))
    k, n = np.array(kn, float).T
    v, k, n = np.array(v)[n > 0], k[n > 0], n[n > 0]
    return fit_h(v, k / n, n)[1] if scorable(k, n) else None


def study(m):
    c = calibrate(m)
    if c is None:
        return None
    P, p0 = np.array([probs(m, c, v, 2000, 1) for v in GRID]).T
    out = {}
    for name, Pc, p0c in (("excl", P, p0), ("dry", p0 + (1 - p0) * P, 0 * p0)):
        at = lambda v: np.interp(np.log(v), LG, Pc)
        V0 = exp_fit(LADDER, at(LADDER))
        lev = V0 * np.geomspace(0.2, 1.6, 8)
        rng = np.random.default_rng(5)
        hs = [protocol(rng, Pc, p0c) for _ in range(REPS)]
        out[name] = dict(V0=V0, h=fit_h(lev, at(lev), np.full(8, 1e5), 0.0)[1], lo=at(lev[0]), hi=at(lev[-1]),
                         p0=np.interp(np.log(lev[0]), LG, p0c), hs=hs)
    return out


def label(m):
    site, p, w1, B = m
    return f"{site:4} " + ("no quench         " if p is None else f"{'diff ' if p == 1.5 else 'front'} w {w1:<4} B {B:<4g}")


def summary(r):
    hs = [h for h in r["hs"] if h is not None]
    held, failed = sum(h <= 3 for h in hs) / REPS, sum(h > 3 for h in hs) / REPS
    med = f"{np.median(hs):5.2f}" if hs else "   - "
    return (f"V0 {r['V0']:6.2f}  h {r['h']:5.2f}  P {r['lo']:.2f}->{r['hi']:.2f}  no xtal@0.2V0 {r['p0']:.2f}  "
            f"HELD {held:.2f} FAILED {failed:.2f} NS {(REPS - len(hs)) / REPS:.2f} med h {med}"), failed


def self_check():
    m = ("line", None, None, 0.0)
    c = 2.0 / rate(m, 1.0, 1.0).sum()        # mu = 2 at 1 uL
    P, p0 = probs(m, c, 1.0)
    assert abs(P - 2 / (np.e + 1)) < 1e-9 and abs(p0 - np.exp(-2)) < 1e-9
    for mm in (m, ("line", 1.5, 1e9, 0.0), ("line", 3, 1e9, 0.0)):   # sampler, and a quench too wide to matter
        N = nuclei(mm, c, 1.0, 20000, 9)
        assert abs(np.mean(N == 0) - p0) < 0.015 and abs(np.mean(2.0 ** (1 - N[N > 0])) - P) < 0.015, mm
    N = nuclei(("line", 1.5, 1e-6, 0.0), c, 1.0, 5000, 9)            # a quench that stops all later nuclei
    assert N.max() == 1


def main():
    assert hashlib.sha256(Path(__file__).with_name("droplet_sim_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore")
    self_check()
    print("Model study, not data.  Each line: true V0 (uL) and true h of the model's curve; true P_single at 0.2 and")
    print(f"1.6 x V0; share of drops with no crystal at 0.2 x V0; then {REPS} runs of the pre-registered protocol.")
    res = {}
    for name in ("excl", "dry"):
        print(f"\n== {name}: " + ("no-crystal drops excluded (prereg)" if name == "excl" else
                                  "no-crystal drops crystallize once at dryness, single (sensitivity, not predicted)"))
        for m in MODELS:
            if m not in res:
                res[m] = study(m)
            print(f"  {label(m)}  " + ("no V0" if res[m] is None else summary(res[m][name])[0]), flush=True)

    ok = {m: r["excl"] for m, r in res.items() if r is not None}
    nq = lambda s: ok[(s, None, None, 0.0)]
    checks = [
        ("M1", nq("line")["h"] <= 1.7, f"line, no quench: h {nq('line')['h']:.2f}"),
        ("M2", nq("bulk")["h"] > 2.0, f"bulk, no quench: h {nq('bulk')['h']:.2f}"),
        ("M3", all(r["h"] <= 3 for r in ok.values()), f"max h {max(r['h'] for r in ok.values()):.2f}"),
        ("M4", all(r["h"] < 2 for m, r in ok.items() if m[1] is not None and m[3] == 10),
         f"max h {max(r['h'] for m, r in ok.items() if m[1] is not None and m[3] == 10):.2f}"),
        ("M5", summary(nq("bulk"))[1] > 0.05, f"bulk, no quench: FAILED {summary(nq('bulk'))[1]:.2f}"),
        ("M6", all(r["lo"] > r["hi"] for r in ok.values()),
         "P rises with V in " + ", ".join(label(m).strip() for m, r in ok.items() if r["lo"] <= r["hi"])),
    ]
    print(f"\nPredictions (excl; {len(res) - len(ok)} model(s) with no V0):")
    for name, held, detail in checks:
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")


if __name__ == "__main__":
    main()
