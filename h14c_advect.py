"""h14c_advect: eddy slowdown of wall clearing at l/w = 10 and 30 with D_mol / D_eff well below 1/8, and enough eddy
kicks before freeze-out (open item 14 / H14, journals/ija_draft.md:252; E86 :187; Section 2.3 :479; 2026-10-21).

Question. h14b_advect.py (journals/new_analyses_h14b_advect.md) measured the slowdown S of wall clearing under explicit
eddies at ONE molecular ratio, D_mol / D_eff = 1/8: S50 = 2.0 (l/w 3.6), 6.3 (l/w 10), and S50 > 8 (censored) at l/w 30,
where the first kick came 77 time units AFTER freeze-out (Tk = 277 against a sweep of 200). The real ratio is ~1e-9, and
the molecular reference (no flow) scales as D_eff / D_mol, so the open question is whether S depends on the ratio. This
script raises D_eff (D_mol = 1, h = 1 fixed, so the molecular wall stays as resolved as in h14b) to 32 and 64, measures S
at l/w 3.6, 10 (and 30 at D_eff = 32), sweeps slowly enough at l/w 30 that 5 kicks happen between lam = 0 and the end of the
sweep (gamma = 0.0008), and re-measures D_eff = 8 in the same code as the anchor.

Model, S, S50, S80, t17, Kaplan-Meier, bootstrap: exactly as in h14b_advect.py (imported, not edited); the only changes are
the per-cell D_eff (t17 = (L/5.5)^2 / D_eff, Stat below), the per-cell gamma and the bias g from the ACTUAL gamma (h14b took g
from gamma = 0.01 for every case, which only matters here: g = Delta / delta(1, eps/V, gamma), obtained by passing
Delta (gamma / 0.01)^(1/4) to K.g_for). Cases, bands P1-P14 and the decision rule: h14c_predictions.txt (asserted below).

PREDICTIONS: h14c_predictions.txt, sha256 PRED_SHA (asserted on import).
Run: cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project \
       --with numpy --with scipy --with numba --with matplotlib python h14c_advect.py > h14c_advect.log 2> h14c_advect.err
     ("report" re-reads h14c_advect_raw.pkl; "smoke" = self-check only; 4 worker processes)
"""
import hashlib
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np

import advect2d as AD
import h14b_advect as H
import kz_patches as K
from sim import delta

HERE = os.path.dirname(os.path.abspath(__file__))
PRED_SHA = "522fb6ea19a961bc877d182e6e9f668d96d9f91efc41cc5258d96006b1ce79ea"
assert hashlib.sha256(open(os.path.join(HERE, "h14c_predictions.txt"), "rb").read()).hexdigest() == PRED_SHA

D_MOL, KAPPA2 = 1.0, H.KAPPA2
PKL = os.path.join(HERE, "h14c_advect_raw.pkl")
GUARD_S = 50 * 60
NWORK = 4
G_H14B_G = 6.27                       # h14b_advect.log G: S50 at l/w 10, D_eff 8, gamma 0.01 (anchor for P1b)


def t17(L, Deff):
    return (L / H.T17_A) ** 2 / Deff


def flow_for(Deff, lw):
    """n = 1 box with eddy l = L / (2 pi) = lw * sqrt(Deff / k); Tk = KAPPA2 l^2 / (8 (Deff - D_mol)) rounded to 0.2."""
    L = int(round(2 * np.pi * lw * np.sqrt(Deff)))
    l = L / (2 * np.pi)
    Tk = round(KAPPA2 * l ** 2 / (8 * (Deff - D_MOL)) / 0.2) * 0.2
    return L, dict(n=1, Tk=float(Tk), A=float(np.sqrt(8 * (Deff - D_MOL) * Tk)), l=float(l))


def cases():
    c = {}

    def add(name, Deff, lw, gamma, S, hold17, seed, h=1.0, dt=0.2, deltas=(0.5, 1.0, 2.0), L=None):
        fl = None
        if L is None:
            L, fl = flow_for(Deff, lw)
        if gamma == "auto":                                  # 5 kicks between lam = 0 and the end of the sweep, never faster than 0.01
            gamma = round(min(0.01, 1 / (5 * fl["Tk"])), 4)
        T = t17(L, Deff)
        c[name] = dict(L=float(L), n=1, h=h, dt=dt, S=S, D=D_MOL, flow=fl, noise="x", seed=seed, gamma=gamma, Deff=Deff, lw=lw,
                       deltas=tuple(deltas), hold17=hold17, hold=float(np.ceil(hold17 * T / 2) * 2),
                       cad=max(2.0, round(T / 60 / 0.2) * 0.2), name=name)
    for j, d in enumerate((0.5, 1.0, 2.0, 1.0, 0.5, 1.0, 2.0, 1.0)):
        add(f"V32_30#{j}", 32, 30, "auto", 1, 26, 3110 + j, deltas=(d,))
    for j in range(2):
        add(f"V64_10#{j}", 64, 10, "auto", 7, 48, 3121 + j)
    for j in range(2):
        add(f"V32_10#{j}", 32, 10, "auto", 16, 30, 3107 + j)
    add("K32", 32, None, 0.01, 4, 90, 3132, L=355)
    add("R8_10#0", 8, 10, "auto", 16, 24, 3102)
    add("R8_10#1", 8, 10, "auto", 16, 24, 3103)
    add("J64", 64, None, 0.01, 32, 150, 3131, L=181)
    add("V32_3.6h", 32, 3.6, 0.01, 16, 24, 3105, h=0.5, dt=0.05)
    add("V64_3.6", 64, 3.6, 0.01, 32, 40, 3120)
    add("J32", 32, None, 0.01, 32, 80, 3130, L=128)
    add("V32_3.6", 32, 3.6, 0.01, 64, 24, 3104)
    add("V32_3.6m", 32, 3.6, 0.0556, 32, 24, 3106)
    add("R8_3.6", 8, 3.6, 0.01, 32, 12, 3101)
    return c


ORDER = list(cases())                                          # longest first (as added)


def job(args):
    name, sp, deadline = args
    H.DELTAS = tuple(d * (sp["gamma"] / K.GAMMA) ** 0.25 for d in sp["deltas"])
    r = H.run_case(sp, deadline)
    r["deltas"] = np.repeat(sp["deltas"], sp["S"])
    return name, r


# ----------------------------------------------------------------------------------------------- statistics
def pool(parts):
    """Merge jobs of one case; histories padded with zeros (a job ends early only when all its runs are single-domain)."""
    tg = max((p["tgrid"] for p in parts), key=len)
    nt = len(tg)
    r = dict(spec=parts[0]["spec"], cut=any(p["cut"] for p in parts), elapsed=sum(p["elapsed"] for p in parts), tgrid=tg)
    for k in ("t_single", "tobs", "wall0", "ndom0", "walled_at_t0", "deltas"):
        r[k] = np.concatenate([p[k] for p in parts])
    r["W"] = np.concatenate([np.pad(p["W"], ((0, nt - len(p["tgrid"])), (0, 0))) for p in parts], axis=1)
    r["runs"] = len(r["t_single"])
    return r


class Stat(H.Stat):
    """h14b Stat with the case's own t17 (D_eff differs per case)."""

    def __init__(self, c, B=2000, seed=0):
        self.c, self.t17 = c, t17(c["spec"]["L"], c["spec"]["Deff"])
        self.R, self.nw = c["runs"], int(c["walled_at_t0"].sum())
        self.pt = self._est(np.arange(self.R))
        rng = np.random.default_rng(seed)
        self.bs = np.array([self._est(rng.integers(0, self.R, self.R)) for _ in range(B)])


# ----------------------------------------------------------------------------------------------- self-check and diagnostics
def kick_damping(Tk):
    """-ln(amplitude after one cubic-Lagrange kick) / q^2 / Tk for cos(q x), wavelength 6.4 cells (q = 2 pi 80 / 512): numerical D_num."""
    N, R, m = 512, 2000, 80
    q = 2 * np.pi * m / N
    x = np.arange(N) + 0.5
    a = np.broadcast_to(np.cos(q * x)[None, None, :], (R, 1, N)).astype(np.float32).copy()
    s = np.random.default_rng(3).uniform(-30, 30, (R, 1))
    b = H.kick(a, s, 1.0, 2)[:, 0, :]
    amp = np.abs((b * np.exp(-1j * q * x)[None, :]).mean(1) * 2).mean()
    return -np.log(amp) / q ** 2 / Tk


def self_check():
    t0 = time.time()
    # (a) advect2d.py:246 fix: nk unchanged for every flow it runs and for every Tk < 100; Tk >= 100 gave 0 and now gives 2
    old = lambda Tk: 2 * int(round(100 / Tk / 2))
    new = lambda Tk: max(2, old(Tk))
    nk_tab = {n: (old(f["Tk"]), new(f["Tk"])) for n, f in AD.FLOWS.items()}
    assert nk_tab == {"S": (1000, 1000), "M": (400, 400), "B": (24, 24)}, nk_tab
    assert all(old(Tk) == new(Tk) for Tk in np.arange(0.2, 100, 0.2)) and old(100.0) == 0 and new(100.0) == 2 and old(277.4) == 0
    # (b) the bias is that of the ACTUAL gamma
    for gam in (0.01, 0.0072, 0.0008):
        V = 1066.0 ** 2
        assert abs(K.g_for(1.0 * (gam / K.GAMMA) ** 0.25, V) * delta(1.0, K.EPS / V, gam) - 1.0) < 1e-12
    # (c) every flow: D_eddy = A^2 / (8 Tk) = D_eff - D_mol (construction) and the Lagrangian tracer agrees (fixed tracer_lagrangian)
    diag = dict(nk=nk_tab, tracer={}, damp={})
    for sp in cases().values():
        fl = sp["flow"]
        if fl is None:
            continue
        assert abs(fl["A"] ** 2 / (8 * fl["Tk"]) - (sp["Deff"] - D_MOL)) < 1e-9
        key = (sp["Deff"], round(sp["lw"], 2))
        if key not in diag["tracer"]:
            AD.L = sp["L"]
            d = np.mean([AD.tracer_lagrangian(fl, s) for s in (1, 7)])
            diag["tracer"][key] = (d, sp["Deff"] - D_MOL, fl["Tk"], sp["L"])
            diag["damp"][key] = kick_damping(fl["Tk"])
    # (d) tiny runs through job(): shapes, labels, history, hold
    sp = dict(L=20.0, n=1, h=1.0, dt=0.2, S=2, D=1.0, flow=dict(n=1, Tk=1.0, A=3.0, l=3.0), noise="x", seed=3, hold=10.0, cad=2.0,
              gamma=0.5, Deff=8, lw=None, deltas=(0.5, 1.0, 2.0), name="smoke")
    _, r = job(("smoke", sp, time.time() + 600))
    assert r["runs"] == 6 and (r["deltas"] == [0.5, 0.5, 1, 1, 2, 2]).all() and r["W"].shape == (len(r["tgrid"]), 6) and np.isfinite(r["W"]).all()
    p = pool([r, r])
    assert p["runs"] == 12 and p["W"].shape[1] == 12
    # (e) the case table
    cs = cases()
    for n, s in cs.items():
        N = s["L"] / s["h"]
        assert abs(N - round(N)) < 1e-9 and s["D"] * s["dt"] / s["h"] ** 2 <= 0.25 + 1e-9
        if s["flow"]:
            assert abs(1 / s["gamma"] / s["flow"]["Tk"] - 5) < 0.7 or s["gamma"] == 0.01 or s["gamma"] == 0.0556, (n, 1 / s["gamma"] / s["flow"]["Tk"])
    sib = H.parse_sibling()
    assert len(sib) == 15 and all("lw" in v and "f5" in v for v in sib.values()), sorted(sib)
    print(f"self-check passed ({time.time() - t0:.0f}s): nk fix {nk_tab} and unchanged for all Tk < 100 (0.2 grid); g from actual gamma; "
          f"D_eddy = D_eff - D_mol by construction; tiny job()/pool() consistent; {len(cs)} jobs; 15 worlds parsed from h14_eddy.log", flush=True)
    return diag


# ----------------------------------------------------------------------------------------------- report
def report(res, elapsed, diag):
    parts = {}
    for k, v in res.items():
        parts.setdefault(k.split("#")[0], []).append(v)
    P = {k: pool(v) for k, v in parts.items()}
    st = {k: Stat(P[k], seed=sum(map(ord, k))) for k in P}
    sp_ = {k: P[k]["spec"] for k in P}
    lw = {k: sp_[k]["flow"]["l"] / np.sqrt(sp_[k]["Deff"]) for k in P if sp_[k]["flow"]}
    fmt, line = H.fmt, H.line

    print("\nNK FIX (advect2d.py:246), old -> new: " + ", ".join(f"{n}: {a} -> {b}" for n, (a, b) in diag["nk"].items())
          + "; Tk 31: 4 -> 4; Tk 100: 0 -> 2; Tk 277.4: 0 -> 2; unchanged for every Tk < 100 (0.2 grid, asserted)")
    print("\nREALISED D_eddy (fixed tracer_lagrangian, advect2d.L set to the box; mean of 2 seeds) and cubic-kick numerical diffusion at wavelength 6.4 cells")
    for (Deff, lwk), (d, nom, Tk, L) in diag["tracer"].items():
        print(f"  D_eff {Deff:3d} l/w {lwk:4.1f} L {L:5.0f} Tk {Tk:6.1f}: D_eddy measured {d:.3f} nominal {nom:.0f} (ratio {d / nom:.3f}); "
              f"D_num/D_mol {diag['damp'][(Deff, lwk)]:.4f}")

    print("\nCASES  (box L, cells N, h, dt, gamma, D_eff, D_mol/D_eff, l/w, Tk, A, kappa, t17, kicks between lam=0 and end of sweep, runs, "
          "walled at freeze-out, hold, cens at end, cut, cpu s)")
    for k in sorted(P, key=lambda x: ORDER.index(x) if x in ORDER else ORDER.index(x + "#0")):
        s = sp_[k]
        fl = s["flow"]
        print(f"  {k:9s} L={s['L']:g} N={int(round(s['L'] / s['h']))} h={s['h']:g} dt={s['dt']:g} g={s['gamma']:g} Deff={s['Deff']} 1/{s['Deff']} "
              + (f"l/w={lw[k]:.2f} Tk={fl['Tk']:g} A={fl['A']:.1f} kappa={fl['A'] / fl['l']:.2f} kicks_pre={1 / s['gamma'] / fl['Tk']:.1f} " if fl else "no flow ")
              + f"t17={t17(s['L'], s['Deff']):.1f} runs={P[k]['runs']} walled={st[k].nw} ({st[k].pt[3]:.3f}) hold={s['hold']:g} ({s['hold17']} t17) "
              f"cens_at_end={int(np.isinf(P[k]['t_single'][P[k]['walled_at_t0']]).sum())} cut={P[k]['cut']} cpu={P[k]['elapsed']:.0f}")

    print("\nSLOWDOWN  (time since the end of the sweep / t17(L, D_eff); Kaplan-Meier over runs walled at freeze-out; bootstrap over runs, 2000x)")
    for k in sorted(P, key=lambda x: ORDER.index(x) if x in ORDER else ORDER.index(x + "#0")):
        s = st[k]
        print(f"  {k:9s} {line(s, 0, 'S50')} | {line(s, 1, 'S80')} | {line(s, 2, 'S_half')}")
        print(f"            walled fraction {s.pt[3]:.3f}  wall0 {s.pt[4]:.3f}  mean domains at t0 {P[k]['ndom0'].mean():.2f}  "
              f"cleared by t17 multiples (all walled runs, KM): " + " ".join(f"{m}:{km_cleared(P[k], s, m):.2f}" for m in (4, 8, 12, 16, 24)))

    sv = lambda k, q=0: st[k].pt[q]
    fl_ = lambda k: st[k].cens_floor()

    def band(k, qi, lo, hi):
        v = sv(k, qi)
        return bool(lo <= v <= hi) if np.isfinite(v) else (False if fl_(k) >= hi else "UNDECIDED")

    def atleast(k, qi, x):
        v = sv(k, qi)
        return bool(v >= x) if np.isfinite(v) else (True if fl_(k) >= x else "UNDECIDED")

    def AND(*xs):
        xs = [x if isinstance(x, str) else bool(x) for x in xs]
        if any(x is False for x in xs):
            return False
        return "UNDECIDED" if any(isinstance(x, str) for x in xs) else True

    def slope(keys, xs):
        xs = np.array(xs, float)
        pt = H.fit_p(xs, np.array([sv(k) for k in keys]))
        pb = np.array([H.fit_p(xs, v) for v in np.stack([st[k].bs[:, 0] for k in keys], 1)])
        fin = np.isfinite(pb)
        return pt, (np.percentile(pb[fin], [2.5, 16, 84, 97.5]) if fin.sum() > 100 else None), 1 - fin.mean()

    def sfmt(x):
        return fmt(x[0]) + (f" [68% {x[1][1]:.2f}..{x[1][2]:.2f}; 95% {x[1][0]:.2f}..{x[1][3]:.2f}]" if x[1] is not None else " (bootstrap mostly undefined)") + f" undefined share {x[2]:.3f}"
    # exponents and ratio dependence
    p36 = slope(("V32_3.6", "V32_10"), (lw["V32_3.6"], lw["V32_10"]))
    p36h = slope(("R8_3.6", "R8_10"), (lw["R8_3.6"], lw["R8_10"])) if "R8_10" in st else (np.nan, None, 1)
    p64 = slope(("V64_3.6", "V64_10"), (lw["V64_3.6"], lw["V64_10"]))
    q36 = slope(("R8_3.6", "V32_3.6", "V64_3.6"), (8, 32, 64))
    q10 = slope(("R8_10", "V32_10", "V64_10"), (8, 32, 64))
    f30 = sv("V32_30") if np.isfinite(sv("V32_30")) else fl_("V32_30")
    plo = (np.log(f30 / sv("V32_10")) / np.log(lw["V32_30"] / lw["V32_10"])) if np.isfinite(sv("V32_10")) else np.nan
    print("\nEXPONENT p in S50 ~ (l/w)^p at fixed ratio (S50 point; bootstrap)")
    print(f"  D_eff 8  (3.6 -> 10): {sfmt(p36h)}   [h14b: 1.11, 95% 0.73-1.54, gamma 0.01 for both]")
    print(f"  D_eff 32 (3.6 -> 10): {sfmt(p36)}")
    print(f"  D_eff 64 (3.6 -> 10): {sfmt(p64)}")
    print(f"  D_eff 32 (10 -> 30): lower bound from the longest observation at l/w 30 ({'S50' if np.isfinite(sv('V32_30')) else 'censored floor'} {f30:.2f}): p >= {fmt(plo)}")
    print("RATIO EXPONENT q in S50 ~ (D_eff / D_mol)^q at fixed l/w (least squares over D_eff = 8, 32, 64)")
    print(f"  l/w 3.6: q = {sfmt(q36)}")
    print(f"  l/w 10 : q = {sfmt(q10)}")
    for nm, ks, qq in (("l/w 3.6", ("R8_3.6", "V32_3.6", "V64_3.6"), q36), ("l/w 10", ("R8_10", "V32_10", "V64_10"), q10)):
        if not np.isfinite(qq[0]):
            print(f"  q at {nm} undefined (a S50 is censored); floor-substituted slope (a lower bound if the censored cell has the largest ratio): "
                  f"{fmt(H.fit_p(np.array([8., 32, 64]), np.array([sv(k) if np.isfinite(sv(k)) else fl_(k) for k in ks])))}")
    g_ratio = sv("R8_10") / G_H14B_G
    mol = {k: sv(k) if np.isfinite(sv(k)) else fl_(k) for k in ("J32", "J64", "K32")}
    for kf, km_ in (("V32_3.6", "J32"), ("V32_10", "K32"), ("V64_3.6", "J64")):
        print(f"  flow / molecular at the same box: {kf} {fmt(sv(kf))} vs {km_} {fmt(sv(km_)) if np.isfinite(sv(km_)) else 'censored > ' + fmt(fl_(km_))}"
              f"  ratio {fmt(sv(kf) / mol[km_])}" + ("" if np.isfinite(sv(km_)) else " (upper bound)"))
    if np.isfinite(q10[0]) and q10[0] > 0 and np.isfinite(sv("V64_10")):
        rs = lambda qq: 64 * (100.0 / sv("V64_10")) ** (1 / qq) if qq > 0 else np.inf
        print(f"  ratio D_eff/D_mol at which S50(l/w 10) would reach the smallest world break-even f (100): point q {rs(q10[0]):.3g}"
              + (f", q at 97.5% {rs(q10[1][3]):.3g}, q at 2.5% {rs(q10[1][0]):.3g}" if q10[1] is not None else "") + "  (extrapolation, not a measurement)")

    def lq(ka, kb, ra, rb):
        v = np.log(st[kb].bs[:, 0] / st[ka].bs[:, 0]) / np.log(rb / ra)
        v = v[np.isfinite(v)]
        return np.log(sv(kb) / sv(ka)) / np.log(rb / ra), (np.percentile(v, [2.5, 97.5]) if len(v) > 100 else (np.nan, np.nan))
    print("LOCAL ratio exponents (descriptive, not predictions; does S saturate with the ratio?)")
    for nm, ks in (("l/w 3.6", ("R8_3.6", "V32_3.6", "V64_3.6")), ("l/w 10", ("R8_10", "V32_10", "V64_10"))):
        for (ka, ra), (kb, rb) in ((((ks[0], 8), (ks[1], 32))), (((ks[1], 32), (ks[2], 64)))):
            q_, ci_ = lq(ka, kb, ra, rb)
            print(f"  {nm}: D_eff {ra} -> {rb}: q = {fmt(q_)} [95% {fmt(ci_[0])}..{fmt(ci_[1])}]")
    mm, v = st["V32_3.6m"], st["V32_3.6"]
    rr = mm.bs[:, 0] / v.bs[:, 0]
    rr = rr[np.isfinite(rr)]
    print(f"CONTROL (descriptive) V32_3.6m (gamma 0.0556, pre-time 1.06 t17) vs V32_3.6 (gamma 0.01, 5.9 t17): S50 {fmt(mm.pt[0])} vs {fmt(v.pt[0])}, "
          f"ratio {fmt(mm.pt[0] / v.pt[0])} [95% {fmt(np.percentile(rr, 2.5)) if len(rr) > 100 else 'n/a'}..{fmt(np.percentile(rr, 97.5)) if len(rr) > 100 else 'n/a'}], "
          f"walled {mm.pt[3]:.3f} vs {v.pt[3]:.3f}")

    print("\nPREDICTIONS (bands from h14c_predictions.txt)")
    out = []

    def verdict(n, ok, txt):
        ok = ok if isinstance(ok, str) else bool(ok)
        out.append((n, ok))
        print(f"P{n} {ok if isinstance(ok, str) else ('HELD' if ok else 'FAILED')}  {txt}")
    a = st["R8_3.6"]
    verdict("1a", AND(0.22 <= a.pt[3] <= 0.45, band("R8_3.6", 0, 1.3, 3.0), band("R8_3.6", 1, 2.8, 7.5)),
            f"R8_3.6: walled {a.pt[3]:.3f} (h14b B 0.323), S50 {fmt(a.pt[0])} (B 2.01), S80 {fmt(a.pt[1])} (B 4.37)")
    verdict("1b", AND(0.67 <= g_ratio <= 1.5 if np.isfinite(g_ratio) else "UNDECIDED", st["R8_10"].pt[3] >= 0.9),
            f"R8_10 S50 {fmt(sv('R8_10'))} / 6.27 (h14b G) = {fmt(g_ratio)}; walled {st['R8_10'].pt[3]:.3f}")
    ok2 = [abs(d / nom - 1) <= 0.05 for d, nom, _, _ in diag["tracer"].values()] + [v <= 0.02 for v in diag["damp"].values()]
    verdict(2, all(ok2), f"{len(diag['tracer'])} flows: D_eddy within 5% of nominal {sum(abs(d / nom - 1) <= 0.05 for d, nom, _, _ in diag['tracer'].values())}/{len(diag['tracer'])}; "
                         f"max D_num/D_mol {max(diag['damp'].values()):.4f} (<= 0.02)")
    r3 = sv("V32_3.6") / sv("V32_3.6h") if np.isfinite(sv("V32_3.6")) and np.isfinite(sv("V32_3.6h")) else np.nan
    dw = abs(st["V32_3.6"].pt[3] - st["V32_3.6h"].pt[3])
    verdict(3, AND((0.67 <= r3 <= 1.5) if np.isfinite(r3) else "UNDECIDED", dw <= 0.15),
            f"S50(h=1)/S50(h=0.5) at 1/32, l/w 3.6 = {fmt(r3)}; walled {st['V32_3.6'].pt[3]:.3f} vs {st['V32_3.6h'].pt[3]:.3f} (diff {dw:.3f})")
    verdict(4, AND(band("V32_3.6", 0, 1.6, 9), band("V32_3.6", 1, 3, 22)), f"V32_3.6: S50 {fmt(sv('V32_3.6'))} (guess 3.5), S80 {fmt(sv('V32_3.6', 1))} (guess 8)")
    verdict(5, AND(band("V32_10", 0, 5, 28), atleast("V32_10", 1, 9), st["V32_10"].pt[3] >= 0.9),
            f"V32_10: S50 {fmt(sv('V32_10'))} (guess 12), S80 {fmt(sv('V32_10', 1))} (>= 9), walled {st['V32_10'].pt[3]:.3f}")
    verdict(6, AND(band("V64_10", 0, 5, 45), st["V64_10"].pt[3] >= 0.9), f"V64_10: S50 {fmt(sv('V64_10'))} (guess 18), walled {st['V64_10'].pt[3]:.3f}")
    verdict(7, AND(atleast("V32_30", 0, 8.0 + 1e-9), atleast("V32_30", 0, 25.0), st["V32_30"].pt[3] >= 0.9),
            f"V32_30: S50 {fmt(sv('V32_30')) if np.isfinite(sv('V32_30')) else 'censored > ' + fmt(fl_('V32_30'))} (a: > 8, the old floor; b: > 25 or undefined), "
            f"cleared by 8 / 16 / 24 t17: {km_cleared(P['V32_30'], st['V32_30'], 8):.2f} / {km_cleared(P['V32_30'], st['V32_30'], 16):.2f} / {km_cleared(P['V32_30'], st['V32_30'], 24):.2f}; "
            f"walled {st['V32_30'].pt[3]:.3f}")
    ci36 = p36[1]
    verdict(8, AND((0.5 <= p36[0] <= 1.6) if np.isfinite(p36[0]) else "UNDECIDED", (ci36[3] >= 1.1) if ci36 is not None else "UNDECIDED",
                   (plo >= 0.5) if np.isfinite(plo) else "UNDECIDED"),
            f"p(3.6->10) at 1/32 = {fmt(p36[0])} (97.5% bound {fmt(ci36[3]) if ci36 is not None else 'n/a'}; p = 1.1 {'NOT excluded' if ci36 is not None and ci36[3] >= 1.1 else 'excluded'}); p(10->30) >= {fmt(plo)}")
    ci10 = q10[1]
    verdict(9, AND((-0.1 <= q36[0] <= 0.9) if np.isfinite(q36[0]) else "UNDECIDED", (0.1 <= q10[0] <= 1.0) if np.isfinite(q10[0]) else "UNDECIDED",
                   (ci10[0] > 0) if ci10 is not None else "UNDECIDED"),
            f"q(l/w 3.6) = {fmt(q36[0])}, q(l/w 10) = {fmt(q10[0])} (2.5% bound {fmt(ci10[0]) if ci10 is not None else 'n/a'})")
    rv = lambda kf, km_: (sv(kf) / mol[km_]) if np.isfinite(sv(kf)) else np.nan
    r_a, r_b = rv("V32_3.6", "J32"), rv("V32_10", "K32")
    verdict(10, AND(band("J32", 0, 25, 70), band("J64", 0, 50, 140), band("K32", 0, 35, 85),
                    (r_a <= 0.5) if np.isfinite(r_a) else "UNDECIDED", (r_b <= 0.8) if np.isfinite(r_b) else "UNDECIDED"),
            f"S50_mol J32 {fmt(sv('J32')) if np.isfinite(sv('J32')) else '>' + fmt(fl_('J32'))}, J64 {fmt(sv('J64')) if np.isfinite(sv('J64')) else '>' + fmt(fl_('J64'))}, "
            f"K32 {fmt(sv('K32')) if np.isfinite(sv('K32')) else '>' + fmt(fl_('K32'))}; flow/mol {fmt(r_a)} (<= 0.5), {fmt(r_b)} (<= 0.8)")
    verdict(11, AND(0.3 <= st["V32_3.6"].pt[3] <= 0.95, st["V64_3.6"].pt[3] >= 0.4), f"walled V32_3.6 {st['V32_3.6'].pt[3]:.3f} (0.3-0.95), V64_3.6 {st['V64_3.6'].pt[3]:.3f} (>= 0.4)")
    sh = [(st[k].pt[2] / st[k].pt[0]) if np.isfinite(st[k].pt[0]) and np.isfinite(st[k].pt[2]) else np.nan for k in ("V32_3.6", "V32_10")]
    verdict(12, AND(*[(0.5 <= x <= 2) if np.isfinite(x) else "UNDECIDED" for x in sh]), "S_half/S50: V32_3.6 " + fmt(sh[0]) + ", V32_10 " + fmt(sh[1]))

    # worlds: illustrative extrapolation at the highest measured ratio (NOT a measurement)
    sib = H.parse_sibling()
    s10 = sv("V64_10") if np.isfinite(sv("V64_10")) else fl_("V64_10")
    pp = p36[0] if np.isfinite(p36[0]) else 1.1
    phi = ci36[3] if ci36 is not None else np.inf
    print(f"\nWORLDS (illustrative, not a measurement): S_ext(lw) = S50(V64_10 = {s10:.1f}{'' if np.isfinite(sv('V64_10')) else ', censored floor'}) x (lw / 10)^p at D_eff/D_mol = 64, "
          f"p = {pp:.2f} (point) or {phi:.2f} (97.5%); f_be = smaller non-None main break-even f (h14_eddy.log)")
    cnt = {}
    for wname in sib:
        for U in ("U0.02", "U1"):
            fs, ln = sib[wname]["f5"][U]
            fe = min(f for f in fs if f is not None)
            med, p90 = sib[wname]["lw"][U]
            sp_pt, sp_hi = s10 * (p90 / 10) ** pp, s10 * (p90 / 10) ** phi
            cnt.setdefault(U, [0, 0])
            cnt[U][0] += sp_pt < fe
            cnt[U][1] += sp_hi < fe
            print(f"  {wname:18s} {U:6s} lw med {med:5.1f} p90 {p90:6.1f} (h14_eddy.log:{sib[wname]['lw_line']}) f_be {fe:g} (:{ln}) "
                  f"S_ext(p90) {sp_pt:8.1f} / {sp_hi:9.1f} -> {'below' if sp_pt < fe else 'ABOVE'} / {'below' if sp_hi < fe else 'ABOVE'}")
    print(f"  worlds (of 15) with S_ext(p90 lw) < f_be: U 0.02: {cnt['U0.02'][0]} (p point) / {cnt['U0.02'][1]} (p 97.5%); U 1: {cnt['U1'][0]} / {cnt['U1'][1]}")
    verdict(13, AND(5 <= cnt["U0.02"][0] <= 12, cnt["U1"][0] <= 8), f"U 0.02: {cnt['U0.02'][0]} in [5, 12]; U 1: {cnt['U1'][0]} <= 8")
    cut = any(v["cut"] for v in res.values())
    verdict(14, (elapsed < GUARD_S) and not cut, f"total {elapsed:.0f} s, {NWORK} workers, cut by guard: {cut}")
    held = [n for n, ok in out if ok is True]
    print("\nSUMMARY  held: " + ", ".join(f"P{n}" for n in held) + " | failed: " + ", ".join(f"P{n}" for n, ok in out if ok is False)
          + " | undecided: " + ", ".join(f"P{n}" for n, ok in out if isinstance(ok, str)))
    if ci10 is None or not np.isfinite(q10[0]):
        dec = "C (undetermined: q at l/w 10 not bounded)"
    elif ci10[0] > 0.15:
        dec = "A (ratio-dependent: q 2.5% bound > 0.15) -> item 14 stays OPEN; the slowdown at l/w 10 grows with D_eff/D_mol"
    elif ci10[3] < 0.15:
        dec = "B (ratio-independent to 1/64: q 97.5% bound < 0.15) -> narrow item 14 by l/w alone, world table above"
    else:
        dec = "C (undetermined: q interval straddles 0.15) -> keep open, report measured S and q interval"
    print(f"DECISION (pre-registered): {dec}; q(l/w 10) = {fmt(q10[0])} [95% {fmt(ci10[0]) if ci10 is not None else 'n/a'}..{fmt(ci10[3]) if ci10 is not None else 'n/a'}]")


def km_cleared(c, s, m):
    """KM fraction of walled runs cleared by m t17 (1 - S(m)); the last observation counts as censoring."""
    w = c["walled_at_t0"]
    T, ev = c["tobs"][w] / s.t17, np.isfinite(c["t_single"][w])
    o = np.lexsort((~ev, T))
    T, ev = T[o], ev[o]
    n, S_ = len(T), 1.0
    for i in range(n):
        if T[i] > m:
            break
        if ev[i]:
            S_ *= 1 - 1.0 / (n - i)
    return 1 - S_


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    t0 = time.time()
    if mode == "report":
        d = pickle.load(open(PKL, "rb"))
        report(d["res"], d["elapsed"], d["diag"])
        return
    diag = self_check()
    if mode == "smoke":
        return
    cs = cases()
    print("\njobs (longest first):", flush=True)
    for k in ORDER:
        s = cs[k]
        print(f"  {k}: L={s['L']:g} N={int(round(s['L'] / s['h']))} h={s['h']:g} dt={s['dt']:g} gamma={s['gamma']:g} D_eff={s['Deff']} runs={s['S'] * len(s['deltas'])} "
              f"hold={s['hold']:g} ({s['hold17']} t17) " + (f"l/w={s['flow']['l'] / np.sqrt(s['Deff']):.2f} Tk={s['flow']['Tk']:g} A={s['flow']['A']:.1f}" if s["flow"] else "no flow"), flush=True)
    deadline = t0 + GUARD_S
    res = {}
    with ProcessPoolExecutor(NWORK) as ex:
        futs = [ex.submit(job, (k, cs[k], deadline)) for k in ORDER]
        for fu in as_completed(futs):
            k, r = fu.result()
            res[k] = r
            print(f"  job done: {k} [{r['elapsed']:.0f}s] walled {int(r['walled_at_t0'].sum())}/{r['runs']} cleared {int(np.isfinite(r['t_single'][r['walled_at_t0']]).sum())} "
                  f"cut={r['cut']} (elapsed {time.time() - t0:.0f}s)", flush=True)
    elapsed = time.time() - t0
    pickle.dump(dict(res=res, elapsed=elapsed, diag=diag), open(PKL, "wb"))
    report(res, elapsed, diag)


if __name__ == "__main__":
    main()
