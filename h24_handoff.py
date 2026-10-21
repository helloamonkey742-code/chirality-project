"""h24_handoff: can a handoff rule close H24, the ~93 % ceiling of the deepest-hand rule on the radial runs? (2026-10-21)

Domains = maximal runs of equal hand; the rule domain holds r*; if it does not survive, the role passes to the domain
above (all 114 plain misses of deephand_scale.log have final hand = hand of the domain above). Variants (frozen on the
600 old runs): V1 threshold 1.25 on d/w; V2, V3 logistic in (d/w, ln drift eps [, ln d/w of the domain above]); chain
exhausted -> the rule domain's hand ('deep'). Tests: 570 fresh runs (seed 20261021, w >= 0.005), 120 thin stacks
(n = 10, w = 0.10, 0.07: the regime where every domain is thinner than a wall width, zero old runs), 354 noise reruns
of old configurations. Descriptive classification of the old misses: `classify`.
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba --with matplotlib \
       python h24_handoff.py run <p> 4 > /dev/null     (p = 0..3, four shards in parallel, writes h24_handoff_runs_<p>.log)
     ... python h24_handoff.py classify > h24_handoff_classify.log ; ... python h24_handoff.py test > h24_handoff_test.log
     (`fit` = the development fit that produced the frozen coefficients, no assertions)
Predictions: h24_handoff_predictions.txt, sha256 090efdb7b1472b28fc1dcdda7bc21aef55e8249452600ca60a849348ceb85afb,
written Wed Oct 21 00:39:24 PDT 2026, before any fresh, thin or noise simulation; asserted below. Never edited after.
Disclosed: `fit` was run before the predictions file existed (old runs only), after scratchpad exploratory fits on
the old runs (logistic fits, CV); the V1-V3 form was chosen with the old runs in view. No radial simulation was run
before the predictions file was hashed. One SMOKE=1 run of h24_handoff_worlds.py (N = 2000, Europa only), also before
the hash, debugged that script: its copy of handoff.hand matched handoff.hand for off / THR 3 / THR 1.25 + fit, and
the bisection shortcut (prediction W2) showed 0 mismatches in 1438 draws for V2 and V3 (so W2 was not blind).
Env note: deephand_scale.one fails under numpy 2 (deephand_scale.py:165, np.bool_ index into a str when hands come
from rng.choice); h24_handoff wraps DH.config to return python ints (types only, no numeric change).
Debug probe after the hash (disclosed): one noise-rerun task was run once in a scratch script to check the numpy-2 patch
(same configuration; Te and final hand equal to the original run); it is re-run in the shards.
"""
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

import handoff as HF            # applies seedfix's default_rng wrapper, imports deephand_scale / multiwall / ...

DH, M = HF.DH, HF.M
_cfg = DH.config                 # numpy 2: deephand_scale.py:165 indexes "+-" with np.bool_ (rng.choice hands) -> TypeError
DH.config = lambda *a, **k: (lambda R, n, lr, rs, u, b, hs: (R, n, lr, rs, u, b, [int(h) for h in hs]))(*_cfg(*a, **k))
HERE = os.path.dirname(os.path.abspath(__file__))
PRED = os.path.join(HERE, "h24_handoff_predictions.txt")
PRED_SHA = "090efdb7b1472b28fc1dcdda7bc21aef55e8249452600ca60a849348ceb85afb"
FROZEN = {"V2": [-2.67, 5.65, 1.5], "V3": [-1.76, 7.33, 1.68, -1.98]}   # development fit, quoted in the predictions file
SEED2 = 20261021                # fresh seed for the out-of-sample runs (old runs: deephand_scale.SEED = 20261012)
CAP, THR1 = 3.0, 1.25           # a domain >= CAP wall widths always holds (365/365 in the old runs); V1's threshold
FRESH_W = DH.WS[:4]             # 0.035, 0.02, 0.01, 0.005 (w = 0.0025 skipped: ~10 min per run)
THIN_W, THIN_RUNS = (0.1, 0.07), 60
LOG = lambda tag: os.path.join(HERE, f"h24_handoff_{tag}.log")


# ---------------------------------------------------------------------------------------------- domains and features
def doms(q):
    """Domains of one ROW run, top -> bottom: (list of dict(hand, dwm, eps, dab), index d0 of the domain containing r*).
    dwm = d/w at the domain middle (handoff.validate's); eps = 2 w_loc(top) (x_top - 1) / r_top, the dimensionless
    drift of the domain's rising top wall (x_top = r_top / r*); dab = dwm of the domain above (nan for the top one)."""
    R, n, lr, w = q["R"], q["n"], q["lr"], q["w"]
    _, _, _, rstar, _, b, _ = DH.config(0, None, lr=lr, u=q["u"], R=R, n=n, alt=True)
    hands = [1 if s == "+" else -1 for s in q["hands"]]
    e = np.concatenate([[R], b, [R - 1.0]])
    cuts = [0] + [m + 1 for m in range(n - 1) if hands[m] != hands[m + 1]] + [n]
    wl = lambda r: w * np.exp(lr * (R - r) / 2)
    out = []
    for d, (a, z) in enumerate(zip(cuts[:-1], cuts[1:])):
        top, bot = e[a], e[z]
        out.append(dict(hand=hands[a], dwm=(top - bot) / wl((top + bot) / 2),
                        eps=2 * wl(top) * (top / rstar - 1) / top if d else np.nan,
                        dab=out[-1]["dwm"] if d else np.nan))
    jl = int(np.sum(b > rstar))
    d0 = next(m for m, (a, z) in enumerate(zip(cuts[:-1], cuts[1:])) if a <= jl < z)
    if "dw" in q:
        assert abs(out[d0]["dwm"] - q["dw"]) <= 1e-3 * max(q["dw"], 1.0) + 1e-3, (q, out[d0]["dwm"])
    return out, d0


def fit_logit(X, y, l2=1e-3):
    X = np.column_stack([np.ones(len(y)), X])
    nll = lambda b: np.sum(np.logaddexp(0, X @ b) - y * (X @ b)) + l2 * np.sum(b[1:] ** 2)
    return minimize(nll, np.zeros(X.shape[1]), method="BFGS").x


COLS = {"V2": ("dwm", "leps"), "V3": ("dwm", "leps", "ldab")}


def feat(dwm, eps, dab, lim=None):
    """Feature dict; ln eps and ln dab clipped to the training range lim = (leps_lo, leps_hi, ldab_lo, ldab_hi)."""
    with np.errstate(all="ignore"):
        le, ld = np.log(eps), np.log(np.maximum(dab, 1e-12))
    if lim is not None:
        le, ld = np.clip(le, lim[0], lim[1]), np.clip(ld, lim[2], lim[3])
    return {"dwm": dwm, "leps": le, "ldab": ld}


def fit_variants(items):
    """Fit V2, V3 on the rule domains with d0 >= 1 and d/w < CAP of items = [(q, doms, d0)]. Returns {V: (b, lim)}."""
    rows = [(q, ds[d0]) for q, ds, d0 in items if d0 >= 1 and ds[d0]["dwm"] < CAP]
    y = np.array([q["ok"] for q, _ in rows], float)
    f0 = feat(np.array([r["dwm"] for _, r in rows]), np.array([r["eps"] for _, r in rows]),
              np.array([r["dab"] for _, r in rows]))
    lim = (f0["leps"].min(), f0["leps"].max(), f0["ldab"].min(), f0["ldab"].max())
    f = feat(np.array([r["dwm"] for _, r in rows]), np.array([r["eps"] for _, r in rows]),
             np.array([r["dab"] for _, r in rows]), lim)
    return {v: (fit_logit(np.column_stack([f[c] for c in cols]), y), lim) for v, cols in COLS.items()}, len(y)


def frozen(fits):
    """The refit on the old runs must still equal the coefficients quoted in the predictions file."""
    assert all(np.allclose(np.round(fits[v][0], 2), c) for v, c in FROZEN.items()), "refit differs from the frozen coefficients"


def holds(var, par, dwm, eps, dab):
    """Does a domain (not the top one) hold? Arrays or scalars. var in {'thr', 'V1', 'V2', 'V3'}."""
    if var == "thr":
        return dwm >= par
    if var == "V1":
        return dwm >= THR1
    b, lim = par
    f = feat(dwm, eps, dab, lim)
    z = b[0] + sum(bi * f[c] for bi, c in zip(b[1:], COLS[var]))
    return (dwm >= CAP) | (z > 0)


def predict(var, par, qd, fb="deep"):
    """Predicted final hand of one ROW run by walking up from the rule domain until a domain holds. The top domain
    holds iff d/w >= THR1 (var 'thr': >= par). If none holds: fb 'deep' -> the rule domain's hand, 'top' -> the top
    domain's hand, 'none' -> 0. Returns (hand, domain used, chain exhausted)."""
    ds, d0 = qd
    for d in range(d0, -1, -1):
        x = ds[d]
        ok = x["dwm"] >= (par if var == "thr" else THR1) if d == 0 else holds(var, par, x["dwm"], x["eps"], x["dab"])
        if ok:
            return x["hand"], d, False
    return {"deep": ds[d0]["hand"], "top": ds[0]["hand"], "none": 0}[fb], d0 if fb == "deep" else -1, True


def agree(rows, cache, var, par=None, fb="deep"):
    """(hit, fired, exhausted) boolean arrays over rows; cache[id] = doms(q). Plain rule: var='thr', par=0."""
    out = []
    for q in rows:
        ds, d0 = cache[id(q)]
        h, d, ex = predict(var, par, (ds, d0), fb)
        out.append((h == q["fin"], d != d0, ex))
    return np.array(out, bool).reshape(-1, 3)


def pct(v):
    return "n=%3d  %5.1f +- %4.1f %%" % (len(v), *(100 * np.array(M.se(v)))) if len(v) else "n=  0"


def table(rows, cache, preds, label):
    """One line per predictor: overall, fired, exhausted; preds = [(name, var, par, fb)]."""
    print(f"  {label}: {len(rows)} runs")
    for name, var, par, fb in preds:
        a = agree(rows, cache, var, par, fb)
        print(f"    {name:22s} overall {pct(a[:, 0])}; fired {pct(a[a[:, 1], 0])}; exhausted chains {a[:, 2].sum()}")


def cv_old(items, rows, cache, reps=5, k=10, seed=0):
    """10-fold x reps CV on the old runs: variants refit on the training folds, predicted on the held-out runs."""
    rng = np.random.default_rng(seed)
    res = {v: [] for v in ("V1", "V2", "V3")}
    for _ in range(reps):
        idx = rng.permutation(len(rows))
        got = {v: np.zeros((len(rows), 3), bool) for v in res}
        for f in np.array_split(idx, k):
            tr = np.setdiff1d(idx, f)
            fits, _ = fit_variants([items[i] for i in tr])
            for v in res:
                par = fits.get(v)
                got[v][f] = agree([rows[i] for i in f], cache, v, par)
        for v in res:
            res[v].append(got[v])
    return {v: np.mean([a[:, 0].mean() for a in r]) for v, r in res.items()}, \
           {v: np.mean([a[a[:, 1], 0].mean() if a[:, 1].any() else np.nan for a in r]) for v, r in res.items()}


def load_old():
    rs = DH.rows()
    cache = {id(q): doms(q) for q in rs}
    return rs, cache, [(q, *cache[id(q)]) for q in rs]


def selfcheck(rs, cache):
    """Hand-checkable: THR = 0 reproduces the ROW ok flag; thr rules reproduce handoff.validate; a 2-domain run."""
    ok = np.array([q["ok"] for q in rs])
    assert (agree(rs, cache, "thr", 0.0)[:, 0] == ok).all()
    for thr in (1.0, 3.0):
        v = HF.validate(rs, thr)
        a = agree(rs, cache, "thr", thr, "none")
        assert (a[:, 0] == v[:, 0]).all() and (a[:, 1] == v[:, 1]).all(), thr
    ds = [dict(hand=1, dwm=5.0, eps=np.nan, dab=np.nan), dict(hand=-1, dwm=0.4, eps=0.1, dab=5.0)]
    assert predict("thr", 1.0, (ds, 1)) == (1, 0, False) and predict("thr", 0.0, (ds, 1)) == (-1, 1, False)
    assert predict("thr", 9.0, (ds, 1), "top")[0] == 1 and predict("thr", 9.0, (ds, 1), "deep")[0] == -1


def main_fit():
    rs, cache, items = load_old()
    selfcheck(rs, cache)
    fits, n = fit_variants(items)
    print(f"training rows (rule domain not the top one, d/w < {CAP:g}): {n}")
    for v, (b, lim) in fits.items():
        print(f"{v}: columns (const, {', '.join(COLS[v])}) coefficients {np.round(b, 2).tolist()}; "
              f"ln eps clip [{lim[0]:.2f}, {lim[1]:.2f}], ln dab clip [{lim[2]:.2f}, {lim[3]:.2f}]")
    preds = [("plain", "thr", 0.0, "none"), ("THR = 1", "thr", 1.0, "none"), ("THR = 3", "thr", 3.0, "none"),
             ("V1 (THR 1.25)", "V1", None, "deep"), ("V2", "V2", fits["V2"], "deep"), ("V3", "V3", fits["V3"], "deep")]
    table(rs, cache, preds, "old 600, in sample")
    cvo, cvf = cv_old(items, rs, cache)
    print("old 600, 10-fold x 5 CV overall agreement: " + "; ".join(f"{v} {100 * a:.1f} % (fired {100 * cvf[v]:.1f})"
                                                                    for v, a in cvo.items()))
    return fits, rs, cache, items



# ----------------------------------------------------------------------------------------------------- runs (sims)
def noise_set(rs, cache):
    """Old runs re-run with fresh noise seeds: every THR = 1 miss, plus every run whose rule domain (not the top one)
    is 0.5 <= d/w < 3 thick at w >= 0.005."""
    a = agree(rs, cache, "thr", 1.0, "none")
    return [q for q, x in zip(rs, a) if not x[0] or (q["w"] >= 0.005 and cache[id(q)][1] >= 1
                                                    and 0.5 <= cache[id(q)][0][cache[id(q)][1]]["dwm"] < 3)]


def tasks(rs, cache):
    """(tag, w, i, rng seed list, sim seed, config kwargs): fresh (F), thin stacks n = 10 at wide walls (T), noise (N)."""
    T = []
    for wi, w in enumerate(FRESH_W):
        T += [("F", w, i, [SEED2, wi, i], SEED2 + 1000 * wi + i, {}) for i in range(DH.RUNS[w])]
    for j, w in enumerate(THIN_W):
        T += [("T", w, i, [SEED2, 10 + j, i], SEED2 + 1000 * (10 + j) + i, {"n": 10}) for i in range(THIN_RUNS)]
    for q in noise_set(rs, cache):
        wi = DH.WS.index(q["w"])
        T += [("N", q["w"], q["i"], [DH.SEED, wi, q["i"]], DH.SEED + 1000 * wi + q["i"] + 7000000 * k, {})
              for k in (1, 2)]
    return T


def main_run(p, P):
    rs, cache, _ = load_old()
    assert hashlib.sha256(open(PRED, "rb").read()).hexdigest() == PRED_SHA, "predictions file changed"
    T = tasks(rs, cache)[p::P]
    orig = {(q["w"], q["i"]): q for q in rs}
    print(f"shard {p}/{P}: {len(T)} sims", flush=True)
    t0 = time.time()
    with open(LOG(f"runs_{p}"), "w") as fh:
        for tag, w, i, rseed, seed, kw in T:
            q = DH.one(w, i, np.random.default_rng(rseed), seed, **kw)
            q.update(tag=tag, seed=seed)
            if tag == "N":                   # same configuration as the old run, different noise seed
                o = orig[(w, i)]
                assert all(q[k] == o[k] for k in ("R", "n", "lr", "u", "hands", "K", "Kd", "dw", "Tc")), (w, i)
            fh.write("ROW " + json.dumps(q) + "\n")
            fh.flush()
    print(f"shard {p}: {time.time() - t0:.0f} s")


def load_runs():
    out = {"F": [], "T": [], "N": []}
    for f in sorted(glob.glob(LOG("runs_*"))):
        for line in open(f):
            if line.startswith("ROW "):
                q = json.loads(line[4:])
                out[q["tag"]].append(q)
    return out


# ---------------------------------------------------------------------------------------------- classify (old runs)
def cv_survival(items, feats, k=10, reps=20, seed=0):
    rows = [(q, ds[d0]) for q, ds, d0 in items if d0 >= 1 and ds[d0]["dwm"] < CAP]
    y = np.array([q["ok"] for q, _ in rows], float)
    f = feat(np.array([r["dwm"] for _, r in rows]), np.array([r["eps"] for _, r in rows]),
             np.array([r["dab"] for _, r in rows]))
    f["lw"] = np.log(f["dwm"])
    X = np.column_stack([f[c] for c in feats])
    rng, acc = np.random.default_rng(seed), []
    for _ in range(reps):
        idx = rng.permutation(len(y))
        pr = np.zeros(len(y))
        for fo in np.array_split(idx, k):
            tr = np.setdiff1d(idx, fo)
            b = fit_logit(X[tr], y[tr])
            pr[fo] = np.column_stack([np.ones(len(fo)), X[fo]]) @ b
        acc.append(np.mean((pr > 0) == y))
    return len(y), np.mean(acc)


def main_classify():
    rs, cache, items = load_old()
    selfcheck(rs, cache)
    ok = np.array([q["ok"] for q in rs])
    print(f"=== h24_handoff classify: the {len(rs)} old runs (deephand_scale.log ROW lines), descriptive ===")
    miss = [(q, *cache[id(q)]) for q in rs if not q["ok"]]
    print(f"plain deepest-hand rule misses {len(miss)}; observed final = hand of the domain above the rule domain in "
          f"{sum(1 for q, ds, d0 in miss if d0 >= 1 and q['fin'] == ds[d0 - 1]['hand'])} of them; "
          f"rule domain is the top one in {sum(1 for q, ds, d0 in miss if d0 == 0)}; unfinished runs (no Te) "
          f"{sum(1 for q in rs if not q['Te'])}")
    a1 = agree(rs, cache, "thr", 1.0, "none")
    ms = [(q, *cache[id(q)], x) for q, x in zip(rs, a1) if not x[0]]
    print(f"\n-- THR = 1 misses: {len(ms)} of {len(rs)} (the 93.3 % ceiling) --")
    print("type U = handed up (rule domain d/w < 1) but the rule domain's own hand won; "
          "type S = kept the rule domain (d/w >= 1) but it lost to the domain above")
    print(f"{'type':4} {'w':>6} {'i':>3} {'R/H':>4} {'n':>2} {'doms':>4} {'K':>2} {'Kd':>2} {'u':>5} {'ln rho':>6} "
          f"{'d/w':>5} {'eps':>7} {'d/w above':>9} {'bottom?':>7} {'Te/Tc':>6} {'Te/tmax':>7}")
    for q, ds, d0, x in sorted(ms, key=lambda t: (t[3][1], -t[0]["w"], t[0]["i"])):
        r = ds[d0]
        print(f"{'U' if x[1] else 'S':4} {q['w']:6g} {q['i']:3d} {q['R']:4g} {q['n']:2d} {len(ds):4d} {q['K']:2d} {q['Kd']:2d} "
              f"{q['u']:5.2f} {q['lr']:6.2f} {r['dwm']:5.2f} {r['eps']:7.4f} {r['dab']:9.2f} {str(d0 == len(ds) - 1):>7} "
              f"{q['Te'] / q['Tc']:6.2f} {q['Te'] / q['tmax']:7.3f}")
    typ = np.array([x[1] for _, _, _, x in ms])
    bott = np.array([d0 == len(ds) - 1 for _, ds, d0, _ in ms])
    dwr = np.array([ds[d0]["dwm"] for _, ds, d0, _ in ms])
    print(f"misses: type U {typ.sum()}, type S {(~typ).sum()}; rule domain is the bottom domain in {bott.sum()}; "
          f"rule d/w in [0,1) {(dwr < 1).sum()}, [1,2) {((dwr >= 1) & (dwr < 2)).sum()}, [2,3) "
          f"{((dwr >= 2) & (dwr < 3)).sum()}, >= 3 {(dwr >= 3).sum()}; by w: "
          + ", ".join(f"{w:g}: {sum(1 for m in ms if m[0]['w'] == w)}" for w in DH.WS)
          + "; by n: " + ", ".join(f"{n}: {sum(1 for m in ms if m[0]['n'] == n)}" for n in DH.NS))
    print("\n-- thickness: P(rule domain holds) by its d/w (runs whose rule domain is not the top one) --")
    sub = [(q, ds[d0]) for q, ds, d0 in items if d0 >= 1]
    for lo, hi in ((0, .5), (.5, 1), (1, 1.5), (1.5, 2), (2, 3), (3, 1e9)):
        g = [q["ok"] for q, r in sub if lo <= r["dwm"] < hi]
        print(f"  d/w [{lo:g}, {hi:g}): {pct(np.array(g))}")
    print("  d/w in [1, 2), by w: " + "; ".join(
        f"{w:g}: {pct(np.array([q['ok'] for q, r in sub if 1 <= r['dwm'] < 2 and q['w'] == w]))}" for w in DH.WS))
    thin = [(q, r, ds, d0) for q, ds, d0 in items if d0 >= 1 and (r := ds[d0])["dwm"] < CAP]
    print(f"\n-- sweep timing, rule domains with 1 <= d/w < 3 and not the top one (held vs lost) --")
    for lab, f in (("held", True), ("lost", False)):
        g = [(q, r) for q, r, ds, d0 in thin if q["ok"] == f and 1 <= r["dwm"] < CAP]
        if not g:
            continue
        te = np.array([q["Te"] / q["Tc"] for q, _ in g])
        tm = np.array([q["Te"] / q["tmax"] for q, _ in g])
        print(f"  {lab}: n={len(g)}  Te/Tc median {np.median(te):.2f} (10-90 %: {np.percentile(te, 10):.2f}-{np.percentile(te, 90):.2f});"
              f"  Te/tmax median {np.median(tm):.3f}, max {tm.max():.3f}")
    print(f"  all {len(rs)} runs: max Te/tmax {max(q['Te'] / q['tmax'] for q in rs):.3f} (>= 1 would mean a cut-off run)")
    print("\n-- layer count: P(rule domain holds), 1 <= d/w < 3, not the top domain --")
    for n in DH.NS:
        print(f"  n = {n}: " + pct(np.array([q["ok"] for q, r, ds, d0 in thin if q["n"] == n and r["dwm"] >= 1])))
    for nd in sorted({len(ds) for _, _, ds, _ in thin}):
        g = [q["ok"] for q, r, ds, d0 in thin if len(ds) == nd and r["dwm"] >= 1]
        if len(g) >= 8:
            print(f"  {nd} domains: {pct(np.array(g))}")
    print("\n-- what predicts holding (training rows: rule domain not the top one, d/w < 3), 10-fold x 20 CV accuracy --")
    for lab, cols in (("d/w", ["dwm"]), ("ln d/w", ["lw"]), ("d/w + ln eps  (V2)", ["dwm", "leps"]),
                      ("d/w + ln d/w above", ["dwm", "ldab"]), ("d/w + ln eps + ln d/w above  (V3)", ["dwm", "leps", "ldab"])):
        n, acc = cv_survival(items, cols)
        print(f"  {lab:36s} n={n}  accuracy {100 * acc:.1f} %")
    print("\n-- d/w-only ceiling: best single threshold on d/w (each domain walked up while d/w < thr; "
          "chain exhausted -> rule domain), tuned on the same runs --")
    grid = np.arange(0.0, 5.001, 0.05)
    sc = [agree(rs, cache, "thr", t, "deep")[:, 0].mean() for t in grid]
    t = int(np.argmax(sc))
    print(f"  old runs: best thr {grid[t]:.2f} -> {100 * sc[t]:.1f} % ({round(sc[t] * len(rs))} / {len(rs)}); "
          f"THR = 1: {100 * sc[20]:.1f} %; plateau (>= best - 0.5 point): thr "
          f"{grid[np.array(sc) >= sc[t] - 0.005].min():.2f}-{grid[np.array(sc) >= sc[t] - 0.005].max():.2f}")
    print()
    main_fit_body(rs, cache, items)


def main_fit_body(rs, cache, items):
    fits, n = fit_variants(items)
    print(f"-- variants V1-V3, fit on the old runs ({n} training rows) --")
    for v, (b, lim) in fits.items():
        print(f"{v}: columns (const, {', '.join(COLS[v])}) coefficients {np.round(b, 2).tolist()}; "
              f"ln eps clip [{lim[0]:.2f}, {lim[1]:.2f}], ln dab clip [{lim[2]:.2f}, {lim[3]:.2f}]")
    preds = [("plain", "thr", 0.0, "none"), ("THR = 1", "thr", 1.0, "none"), ("THR = 3", "thr", 3.0, "none"),
             ("V1 (THR 1.25)", "V1", None, "deep"), ("V2", "V2", fits["V2"], "deep"), ("V3", "V3", fits["V3"], "deep")]
    table(rs, cache, preds, "old runs, in sample")
    cvo, cvf = cv_old(items, rs, cache)
    print("old runs, 10-fold x 5 CV overall agreement: " + "; ".join(f"{v} {100 * a:.1f} % (fired {100 * cvf[v]:.1f})"
                                                                    for v, a in cvo.items()))
    return fits, cvo


# ---------------------------------------------------------------------------------------------------- test (fresh)
def main_test():
    assert hashlib.sha256(open(PRED, "rb").read()).hexdigest() == PRED_SHA, "predictions file changed"
    rs, cache, items = load_old()
    selfcheck(rs, cache)
    R = load_runs()
    fits, cvo, SC = None, None, {}
    print("=== h24_handoff test: variants frozen on the old runs, tested on fresh simulations ===")
    fits, cvo = main_fit_body(rs, cache, items)
    frozen(fits)
    for tag in "FTN":
        print(f"{tag}: {len(R[tag])} runs read")
    # --- noise
    print("\n-- noise: the same configurations, 2 fresh noise seeds each (final hand vs the original run) --")
    orig = {(q["w"], q["i"]): q for q in rs}
    by = {}
    for q in R["N"]:
        by.setdefault((q["w"], q["i"]), []).append(q)
    a1 = agree(rs, cache, "thr", 1.0, "none")
    miss = {(q["w"], q["i"]) for q, x in zip(rs, a1) if not x[0]}
    for lab, keys in (("THR = 1 misses", [k for k in by if k in miss]), ("borderline (not misses)", [k for k in by if k not in miss])):
        same = [r["fin"] == orig[k]["fin"] for k in keys for r in by[k]]
        anyflip = sum(1 for k in keys if any(r["fin"] != orig[k]["fin"] for r in by[k]))
        print(f"  {lab}: {len(keys)} configurations x {len(by[keys[0]]) if keys else 0} reruns; same final hand in "
              f"{sum(same)} / {len(same)} reruns; configurations with any flip {anyflip} / {len(keys)}")
        SC["noise_" + lab[:3]] = sum(same) / max(len(same), 1)
    print("  largest Te change among reruns with the same hand: "
          f"{max((abs(r['Te'] / orig[k]['Te'] - 1) for k in by for r in by[k] if r['fin'] == orig[k]['fin']), default=float('nan')):.3f}")
    # --- fresh
    F = R["F"]
    cF = {id(q): doms(q) for q in F}
    cache.update(cF)
    itF = [(q, *cF[id(q)]) for q in F]
    print(f"\n-- fresh runs (seed {SEED2}): {len(F)} --")
    print(f"fresh runs with rule domain d/w < 3 and not the top one: {sum(1 for _, ds, d0 in itF if d0 >= 1 and ds[d0]['dwm'] < CAP)}; "
          f"plain rule holds in {100 * np.mean([q['ok'] for q in F]):.1f} %; runs without Te: {sum(1 for q in F if not q['Te'])}")
    preds = [("plain", "thr", 0.0, "none"), ("THR = 1", "thr", 1.0, "none"), ("THR = 3", "thr", 3.0, "none"),
             ("V1 (THR 1.25)", "V1", None, "deep"), ("V2", "V2", fits["V2"], "deep"), ("V3", "V3", fits["V3"], "deep")]
    table(F, cache, preds, "fresh, out of sample")
    best = next(v for v in ("V1", "V2", "V3") if cvo[v] >= max(cvo.values()) - 0.002)
    print(f"decision rule (old-run CV overall, ties within 0.2 point -> fewer parameters): best variant = {best}")
    a = agree(F, cache, best, fits.get(best), "deep")
    val = a[:, 0].mean() >= 0.95 and a[a[:, 1], 0].mean() >= 0.85
    SC.update(best=best, ov=a[:, 0].mean(), fi=a[a[:, 1], 0].mean(), val=val,
              A={n: agree(F, cache, v, par, fb)[:, 0].mean() for n, v, par, fb in preds})
    print(f"{best} on fresh runs: {pct(a[:, 0])}; fired {pct(a[a[:, 1], 0])}  -> "
          f"{'VALIDATED' if val else 'NOT validated'} against the pre-registered 95 % overall / 85 % fired")
    ws = np.array([q["w"] for q in F])
    print("  by w: " + "; ".join(f"{w:g}: {pct(a[ws == w, 0])}" for w in FRESH_W))
    ns_ = np.array([q["n"] for q in F])
    print("  by n: " + "; ".join(f"{n}: {pct(a[ns_ == n, 0])}" for n in DH.NS))
    SC["strata"] = min([a[ws == w, 0].mean() for w in FRESH_W] + [a[ns_ == n, 0].mean() for n in DH.NS])
    print("  fallback alternatives (chain exhausted): " + "; ".join(
        f"{fb}: {100 * agree(F, cache, best, fits.get(best), fb)[:, 0].mean():.1f} %" for fb in ("deep", "top", "none")))
    for v in ("V1", "V2", "V3"):
        if v != best:
            b = agree(F, cache, v, fits.get(v), "deep")
            print(f"  {v}: {pct(b[:, 0])}; McNemar vs {best}: {int((a[:, 0] & ~b[:, 0]).sum())} runs only {best} gets, "
                  f"{int((~a[:, 0] & b[:, 0]).sum())} only {v} gets")
    print("  remaining misses of " + best + " on fresh runs:")
    for q, x in zip(F, a):
        if not x[0]:
            ds, d0 = cache[id(q)]
            r = ds[d0]
            print(f"    w {q['w']:g} i {q['i']:3d} n {q['n']:2d} doms {len(ds)} d/w {r['dwm']:.2f} eps {r['eps']:.4f} "
                  f"above {r['dab']:.2f} bottom {d0 == len(ds) - 1} fired {bool(x[1])}")
    print("\n-- d/w-only ceiling on fresh runs (threshold tuned ON the fresh runs = an oracle) --")
    grid = np.arange(0.0, 5.001, 0.05)
    sc = [agree(F, cache, "thr", t, "deep")[:, 0].mean() for t in grid]
    t = int(np.argmax(sc))
    SC["oracle"] = sc[t]
    print(f"  best thr {grid[t]:.2f} -> {100 * sc[t]:.1f} %; thr tuned on the old runs (1.25, V1): "
          f"{100 * agree(F, cache, 'V1', None, 'deep')[:, 0].mean():.1f} %")
    # --- thin stacks and the all-thin regime
    print(f"\n-- thin stacks (w = {THIN_W}, n = 10) and the all-thin regime (every domain down to the rule domain thin) --")
    pool = [(q, "old") for q in rs] + [(q, "fresh") for q in F] + [(q, "thin") for q in R["T"]]
    for q, _ in pool:
        if id(q) not in cache:
            cache[id(q)] = doms(q)
    for lab, lo, hi in (("max d/w < 1", 0, 1), ("max d/w in [1, 1.25)", 1, 1.25), ("max d/w in [1.25, 1.5)", 1.25, 1.5),
                        ("max d/w in [1.5, 2)", 1.5, 2)):
        g = [(q, src) for q, src in pool if lo <= max(x["dwm"] for x in cache[id(q)][0][:cache[id(q)][1] + 1]) < hi]
        if not g:
            print(f"  {lab}: n=0")
            continue
        deep = np.array([q["ok"] for q, _ in g])
        top = np.array([q["fin"] == cache[id(q)][0][0]["hand"] for q, _ in g])
        print(f"  {lab}: n={len(g)} (old {sum(s == 'old' for _, s in g)}, fresh {sum(s == 'fresh' for _, s in g)}, thin "
              f"{sum(s == 'thin' for _, s in g)}); final = rule (deepest) domain's hand {pct(deep)}; "
              f"final = top domain's hand {pct(top)}")
    G = [q for q, _ in pool if max(x["dwm"] for x in cache[id(q)][0][:cache[id(q)][1] + 1]) < 1]
    SC["nG"] = len(G)
    SC["Gdeep"] = np.mean([q["ok"] for q in G]) if G else np.nan
    SC["Gtop"] = np.mean([q["fin"] == cache[id(q)][0][0]["hand"] for q in G]) if G else np.nan
    T = R["T"]
    print(f"  thin stacks: {len(T)} runs, plain rule holds in {pct(np.array([q['ok'] for q in T]))}; "
          f"V3 agreement {pct(agree(T, cache, 'V3', fits['V3'], 'deep')[:, 0])}, with fallback top "
          f"{pct(agree(T, cache, 'V3', fits['V3'], 'top')[:, 0])}")
    scorecard(SC)


def scorecard(SC):
    """Pre-registered predictions P1-P11 (h24_handoff_predictions.txt), pass / fail from the numbers above."""
    A, G = SC["A"], SC["nG"] >= 20
    L = [("P1", SC["val"], f"{SC['best']} overall {100 * SC['ov']:.1f} % (>= 95), fired {100 * SC['fi']:.1f} % (>= 85)"),
         ("P2", SC["best"] == "V3", f"decision rule picked {SC['best']}"),
         ("P3", A["V3"] - A["V1 (THR 1.25)"] >= 0.01, f"V3 - V1 = {100 * (A['V3'] - A['V1 (THR 1.25)']):+.1f} points (>= +1.0)"),
         ("P4", A["V1 (THR 1.25)"] - A["THR = 1"] >= 0.005, f"V1 - THR1 = {100 * (A['V1 (THR 1.25)'] - A['THR = 1']):+.1f} points (>= +0.5)"),
         ("P5", A["THR = 3"] < 0.85, f"THR=3 {100 * A['THR = 3']:.1f} % (< 85)"),
         ("P6", SC["oracle"] <= 0.96 and SC["ov"] - SC["oracle"] >= 0.005,
          f"d/w-only oracle {100 * SC['oracle']:.1f} % (<= 96.0); primary - oracle = {100 * (SC['ov'] - SC['oracle']):+.1f} points (>= +0.5)"),
         ("P7", SC["strata"] >= 0.90, f"worst w or n stratum {100 * SC['strata']:.1f} % (>= 90)"),
         ("P8", SC["noise_THR"] >= 0.90 and SC["noise_bor"] >= 0.85,
          f"reruns with the same hand: THR=1 misses {100 * SC['noise_THR']:.1f} % (>= 90), borderline {100 * SC['noise_bor']:.1f} % (>= 85)"),
         ("P9", G, f"n(G) = {SC['nG']} (>= 20)"),
         ("P10", None if not G else SC["Gtop"] - SC["Gdeep"] >= 0.10,
          f"G: top-hand {100 * SC['Gtop']:.1f} % vs rule-hand {100 * SC['Gdeep']:.1f} % (difference >= 10 points)"),
         ("P11", None if not G else max(SC["Gtop"], SC["Gdeep"]) < 0.85, f"G: best fixed fallback {100 * max(SC['Gtop'], SC['Gdeep']):.1f} % (< 85)")]
    print("\n== scorecard (pass / fail rules in h24_handoff_predictions.txt) ==")
    for k, v, t in L:
        print(f"  {k:3s} {'NOT TESTABLE' if v is None else 'PASS' if v else 'FAIL':12s} {t}")


if __name__ == "__main__":
    P = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    {"fit": main_fit, "classify": main_classify, "test": main_test,
     "run": lambda: main_run(int(sys.argv[2]), P)}[sys.argv[1]]()
