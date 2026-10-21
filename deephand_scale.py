"""deephand_scale: does multiwall's deepest-hand rule hold at ocean-like wall width / depth, and how long does the
cascade really take? (2026-10-12)

Why: multiwall.py found that the layer containing the onset depth r* (the bottom layer if every wall is past onset)
fills the ocean: 336/336 (all past onset) and 94/108 (onset inside), but only for walls 0.035-0.07 of the depth wide;
the 14 failures all had the deciding domain <= 4.4 wall widths thick, and ">= 5 wall widths: 63/63" was post hoc
(ija_draft H24). Its wrapper timed the cascade with a constant-D_z closed form T_c, and the simulated time was
0.006-5.05 x T_c, median 1.2 (H23). Earth's 99.9 % and every main multiwall share rest on both.

Model: multiwall's radial solver, imported unchanged (multiwall.shell / simulate / transit / wrapper_times; numba
Thomas implicit diffusion, explicit reaction, dt = 0.1). Geometry is the WRAPPER's own: constant G = 2/r*, so
x(r) = r/r* (every wall has its own x; walls above r* are past onset and rise, walls below sink) and
D(r) = D_top exp(4 (R - r)/r*), i.e. bottom/top ratio rho = exp(4H/r*), r* = 4H/ln rho, exactly sweep2/sweep3's
r* = 4H/ln(rho). H = 1. Wall width w = sqrt(2 D_top) (the narrowest wall, at the ice); grid dz = w/12 (uniform; the
grid shrinks with w, so cost ~ w^-3).
Configurations ("scale" runs): w = 0.035, 0.02, 0.01, 0.005, 0.0025; R/H cycles 2, 3, 6, 12, 100 (Enceladus ~6, Europa ~12,
Early Earth ~2400); n = 3, 6, 10 layers; ln rho ~ U(1.2 x 4H/R, ln 100) (sweep2's log-uniform rho, kept where the
top wall is clearly past onset, x(R) >= 1.2). Layers are drawn as sweep2 draws them: n - 1 full layers of equal
thickness in stretched units (thickness ~ sqrt(D_z(r)), sweep2's d_w = L sqrt(D_z/D_h) with the local D_z) from the
ice down, plus a bottom remnant u ~ U(0, 1) of a full layer (sweep2/sweep3: the bottom layer is H - K d_w). Hands:
alternating (2 runs in 3, as the wrapper assumes) or random (not all equal). Onset falls inside the stack (a sub-onset
wall sinks) in roughly 1 run in 6 (mostly R/H = 2, 3). Initial state, noise and stopping as multiwall (a run ends when
one hand is gone). tau_max = 3 x the lone-wall transit integral of the slowest deciding wall + 100.
Rule checked: the hand of the layer containing r* (the bottom layer if r* is below the floor) fills the shell.
d/w = thickness of the same-hand domain containing that layer over the local wall width sqrt(2 D) at the domain's
middle (multiwall's post-hoc definition). K = walls above r* (past onset).
Timing: T_c = multiwall.wrapper_times (constant D = depth mean, T_up of the deepest rising wall + T_down of the
shallowest sinking wall). Fit, on finished runs where the rule held:
    log10(T_end / T_c) = c0 + c1 ln(rho) + c2 log10(K + 1) + c3 log10(d/w)   (least squares; h = 10^fit).
Sections (each its own run, `nice -n 10 python deephand_scale.py <s>`):
 (1 <w> [i0 i1]) scale runs at one w (runs i0..i1-1): 150 / 150 / 150 / 120 / 30 runs at w = 0.035 / 0.02 / 0.01 /
     0.005 / 0.0025 (0.0025 in two halves, each < 20 min); one ROW json line per run (read back by section 3).
 (2) convergence: R/H = 6, n = 6, ln rho = 2, u = 0.6, alternating hands, at w = 0.035, 0.01, 0.005, 0.0025, each at
     (dt, dz) and (dt/2, dz/2).
 (3) tables from the ROW lines: P(rule hand fills) +- SE by w, by d/w bin, by K, by n, sub-onset; the timing fit;
     then multiwall section 4 re-timed (import multiwall, not edited): multiwall.multi gets the per-draw h as its
     `scale` (which multiplies T_up and T_down). Per draw, h uses the wrapper's K (clipped to the simulated range),
     ln rho, and the deciding layer's d/w = f H / sqrt(2 D_loc / k) (k = k2 c the reaction rate, D_loc = the drawn
     depth-mean D_z times the exponential profile at the layer's middle; clipped to the simulated range).
     Regression: h = 1 everywhere reproduces multiwall run4(0) main exactly (same arrays) and multiwall.log's main
     table to the printed 0.1. Also an "ocean-scale" variant: each draw the rule decides is weighted by section 1's
     measured P(fill) in its d/w bin (0 below d/w = 1: thinner than a wall, not a layer).
  ASSUMED: d/w and K summarise what matters for the rule (R/H and x enter only through T_c and rho); the wrapper's
  layers are alternating (as sweep3); its K (median 26-120 on the moons) is far beyond the simulated K <= 9, so h is
  extrapolated by clipping.

PREDICTIONS (written 2026-10-12, before any result-producing run):
 (P1) Section 1, pooled over w: when the deciding domain is >= 5 local wall widths thick, the rule hand fills in
      >= 97 % of runs (and in every w level >= 93 %).
 (P2) Section 1: when d/w < 2 the rule hand fills in <= 50 % of runs (pooled).
 (P3) Section 1: at d/w >= 5 the fill fraction at w <= 0.01 and at w >= 0.02 differ by less than 2 combined SE
      (no effect of scale beyond d/w).
 (P4) Section 1: the overall fill fraction (all d/w) at w <= 0.005 (pooled) is >= 85 % and >= the fraction at
      w = 0.035.
 (P5) Section 1: runs with K >= 5 walls past onset and d/w >= 5 fill in >= 95 %.
 (P6) Section 1: runs with a sub-onset wall in the stack (a wall below r*) and d/w >= 5: the rule fails in <= 1 run.
 (P7) Section 1: <= 2 % of runs unfinished at tau_max, at every w.
 (P8) Timing: T_end / T_c over finished rule-held runs has median in [1.0, 3.0] and max <= 30.
 (P9) Timing fit: c1 (ln rho) > 0 (the constant-D form misses the slow top), and the fit's residual RMS log10 is
      <= 0.6 x the RMS of log10(T_end / T_c) (h = 1).
 (P10) Section 2: halving dt and dz changes T_end by < 3 % at every w and keeps the final hand.
 (P11) Section 3 regression: h = 1 reproduces multiwall main exactly (arrays) and multiwall.log to the printed 0.1.
 (P12) Section 3, fitted h: every world's P_pved (no lid and final) lies between multiwall x10 - 0.2 and multiwall
       main + 0.2; Early Earth >= 99.0.
 (P13) Section 3: of the draws the rule decides and counts P_pved, >= 20 % have the deciding layer < 5 local wall
       widths thick in every icy moon; <= 5 % on Early Earth.
 (P14) Section 3, ocean-scale (d/w-weighted, fitted h): Early Earth P_pved >= 95 % (no lid and final); every icy
       moon loses >= 5 points against multiwall main (no lid).

RESULTS (deephand_scale.log, 2026-10-12):
 (P1) HELD. d/w >= 5: 284/284 filled (100.0 +- 0.0 %), every w level 100 % (n = 20, 52, 93, 93, 26). In fact
      d/w >= 3: 365/365; the thickest failing domain was 2.85 wall widths.
 (P2) HELD. d/w < 2: 55/167 = 32.9 +- 3.6 % (d/w < 1: 8.0 +- 2.9 %; 1-2: 60.8 +- 5.5 %).
 (P3) HELD. d/w >= 5: 100.0 +- 0.0 at w <= 0.01 and at w >= 0.02; difference / combined SE = 0.00.
 (P4) HELD. w <= 0.005 pooled: 141/150 = 94.0 +- 1.9 % (>= 85, and > 67.3 +- 3.8 % at w = 0.035). The overall
      fraction rises as w shrinks only because more layers are many wall widths thick at small w.
 (P5) HELD. K >= 5 and d/w >= 5: 81/81.
 (P6) HELD. sub-onset wall in the stack and d/w >= 5: 0 failures in 44 (all d/w: 54/57 = 94.7 +- 3.0 %).
 (P7) HELD. 0 unfinished at every w.
 (P8) HELD. 486 finished rule-held runs: T_end / T_c median 1.731, max 11.1 (min 0.004: a run whose minority hand
      was thin layers that died at t = 16 before any wall moved).
 (P9) FAILED (half). c1 = +0.183 > 0 held, but residual RMS log10 0.392 is 0.79 x the h = 1 RMS 0.494 (needed
      <= 0.6). Fit: log10(T_end/T_c) = -0.570 + 0.183 ln rho - 0.048 log10(K+1) + 0.343 log10(d/w). The scatter
      is mostly at large w (residual RMS 0.60 at w = 0.035, 0.12 at w = 0.0025); median T_end/T_c drifts up with
      shrinking w (1.25 -> 2.21), a trend the fit has no w term for.
 (P10) HELD. dt/2, dz/2 changes T_end by -1.72, -1.93, -1.93, -1.93 % at w = 0.035, 0.01, 0.005, 0.0025; same
       final hand each time.
 (P11) HELD. "EXACT": h = 1 arrays equal multiwall run4(0) main, and all 15 worlds match multiwall.log to 0.1.
 (P12) HELD. Fitted h moves P_pved by 0.0 (Early Earth, 99.9) to -1.4 points (Ganymede, Titan); every world, no lid
       and final, lies in [x10 - 0.2, main + 0.2].
 (P13) FAILED (half). Icy moons: 59-74 % of rule-decided P_pved draws have the deciding layer < 5 wall widths
       (>= 20 % held). Early Earth: 51.7 % (needed <= 5 %); 29.4 % are < 1 wall width. The deciding layer is the
       bottom remnant H - K d_w, drawn anywhere between 0 and one full layer, so about half are thin on every world.
 (P14) FAILED (half). Ocean-scale Early Earth P_pved 89.0 % (no lid and final; needed >= 95). Every icy moon loses
       >= 5 points no lid (Mimas -5.8 to Umbriel -15.3) - held.
 Run notes: section 1 was run in chunks (resumable per (w, i); each ROW is deterministic in (w, i)); jobs at
 w = 0.005 / 0.0025 were killed under machine load 17-18 and resumed from the first missing i, so no run was
 repeated or dropped (600 ROW lines, each (w, i) once). w = 0.0025 has only 30 runs. Section 3 ran in 12 s.
 The ocean-scale variant moves a failed thin-remnant draw to P_none; in the simulations the next layer up then
 fills, and in the wrapper that layer may itself be PVED-set, so 89.0 is a floor, not an estimate.
"""
import json
import os
import re
import sys
import time

import numpy as np

import multiwall as M
import sweep2 as S
import sweep3 as S3
import uncertainty as U

SMOKE = os.environ.get("SMOKE") == "1"     # SMOKE=1 python deephand_scale.py: code paths + cost probe, no results
SEED = 20261012
WS = (0.035, 0.02, 0.01, 0.005, 0.0025)
RUNS = {0.035: 150, 0.02: 150, 0.01: 150, 0.005: 120, 0.0025: 30}
RH, NS, DZW = (2.0, 3.0, 6.0, 12.0, 100.0), (3, 6, 10), 12
BINS = (0.0, 1.0, 2.0, 3.0, 5.0, 10.0, np.inf)
HERE = os.path.dirname(os.path.abspath(__file__))


def config(i, rng, lr=None, u=None, R=None, n=None, alt=None):
    """Wrapper geometry for run i: (R, n, ln rho, r*, u, layer bounds top->bottom, hands top->bottom)."""
    H = 1.0
    R = RH[i % 5] if R is None else R
    n = NS[(i // 5) % 3] if n is None else n
    lr = rng.uniform(1.2 * 4 * H / R, np.log(100)) if lr is None else lr
    u = rng.uniform(0, 1) if u is None else u
    rs, a = 4 * H / lr, lr / (2 * H)                      # a = 2/r*: sqrt(D) ~ exp(a (R - r))
    s = (1 - np.exp(-a * H)) / a * np.arange(1, n) / (n - 1 + u)   # equal stretched thickness, remnant u at bottom
    bounds = R + np.log(1 - a * s) / a
    if (i % 3 != 2) if alt is None else alt:
        hands = [(-1) ** j for j in range(n)]
    else:
        hands = list(rng.choice([-1, 1], n))
        while abs(sum(hands)) == n:
            hands = list(rng.choice([-1, 1], n))
    return R, n, lr, rs, u, bounds, hands


def one(w, i, rng, seed, dt=M.DT, dzw=DZW, tcap=None, **kw):
    R, n, lr, rs, u, bounds, hands = config(i, rng, **kw)
    H, Dt = 1.0, w * w / 2
    xf = lambda r, rs=rs: np.asarray(r) / rs
    _, D, _, _ = M.shell(R, H, Dt, xf, w / dzw, dt)
    walls = [b for b, h1, h2 in zip(bounds, hands[:-1], hands[1:]) if h1 != h2]
    upw, dnw = [b for b in walls if b > rs], [b for b in walls if b < rs]
    Tup = M.transit(min(upw), R, H, D, xf, True) if upw else 0.0
    Tdn = M.transit(max(dnw), R, H, D, xf, False) if dnw else 0.0
    Dm = D(np.linspace(R - H, R, 4001)).mean()
    Tc = sum(M.wrapper_times(R, H, rs, min(upw) if upw else None, max(dnw) if dnw else None, Dm))
    j = int(np.sum(bounds > rs))                         # layer containing r* (bottom if r* below the floor)
    lo, hi = j, j
    while lo > 0 and hands[lo - 1] == hands[j]:
        lo -= 1
    while hi < n - 1 and hands[hi + 1] == hands[j]:
        hi += 1
    e = np.concatenate([[R], bounds, [R - H]])
    dw = (e[lo] - e[hi + 1]) / np.sqrt(2 * D((e[lo] + e[hi + 1]) / 2))
    tmax = 3 * max(Tup, Tdn) + 100 if tcap is None else tcap
    t0 = time.time()
    Te, fin, _, _ = M.simulate(R, H, Dt, xf, bounds, hands, seed, tmax, dz=w / dzw, dt=dt)
    return dict(w=w, i=i, R=R, n=n, lr=round(lr, 4), u=round(u, 4), hands="".join("+-"[int(h < 0)] for h in hands),
                K=len(upw), Kd=len(dnw), dw=round(float(dw), 3), Te=Te, fin=fin, ok=bool(Te and fin == hands[j]),
                Tc=float(Tc), Tm=float(max(Tup, Tdn)), tmax=float(tmax), sec=round(time.time() - t0, 1))


def section1(w, i0=0, i1=None):
    wi, i1 = WS.index(w), RUNS[w] if i1 is None else i1
    print(f"\n=== Section 1: scale runs, w = {w} (w/H at the ice; dz = w/{DZW}, dt = {M.DT}), runs {i0}..{i1 - 1} ===")
    t0 = time.time()
    for i in range(i0, i1):
        q = one(w, i, np.random.default_rng([SEED, wi, i]), SEED + 1000 * wi + i)
        print("ROW " + json.dumps(q), flush=True)
    print(f"section 1 w = {w}: {time.time() - t0:.0f} s")


def section2():
    print("\n=== Section 2: convergence (R/H = 6, n = 6, ln rho = 2, u = 0.6, alternating) ===")
    for w in (0.035, 0.01, 0.005, 0.0025):
        out = []
        for dt, dzw in ((M.DT, DZW), (M.DT / 2, 2 * DZW)):
            q = one(w, 0, np.random.default_rng(0), SEED, dt=dt, dzw=dzw, R=6.0, n=6, lr=2.0, u=0.6, alt=True)
            out.append(q)
            print(f"  w = {w}: dt = {dt}, dz = w/{dzw}: T_end = {q['Te']}, final {q['fin']:+d}, rule hand "
                  f"{'filled' if q['ok'] else 'did NOT fill'}, T_c = {q['Tc']:.0f}, {q['sec']} s", flush=True)
        a, b = out
        print(f"  w = {w}: change in T_end {100 * (b['Te'] / a['Te'] - 1):+.2f} %, same final hand: {a['fin'] == b['fin']}")


def rows():
    return [json.loads(l[4:]) for l in open(os.path.join(HERE, "deephand_scale.log")) if l.startswith("ROW ")]


def fit(rs):
    g = [q for q in rs if q["ok"] and q["Tc"] > 0]
    X = np.array([[1, q["lr"], np.log10(q["K"] + 1), np.log10(q["dw"])] for q in g])
    y = np.array([np.log10(q["Te"] / q["Tc"]) for q in g])
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    return c, X, y


def pbin(rs):
    """Measured P(rule hand fills) per d/w bin (pooled); 0 below d/w = 1."""
    p = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        g = [q["ok"] for q in rs if lo <= q["dw"] < hi]
        p.append(0.0 if lo < 1 or not g else float(np.mean(g)))
    return np.array(p)


def geom(x, L, dz_rng, ex, lid):
    """multiwall.multi's own lines (copied, not edited) for rec, K, ln rho, and the deciding layer's d/w."""
    c, k2, Dh, V = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "volume"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    k, H, lr, t = k2 * c, U.depth(V, L), np.log(ex["rho"]), x["time"]
    R = L / np.pi
    column = H * np.sqrt(Dh / Dz) > L
    d_w = np.minimum(H, L * np.sqrt(Dz / Dh))
    r_w = R - d_w
    rev = column & (lr >= 4 * H / r_w) & (k * t >= 10) & ~lid
    rec = rev & ((H - d_w) * np.sqrt(Dh / Dz) > L) & (lr >= 4 * H / (r_w - d_w))
    with np.errstate(all="ignore"):
        rs = 4 * H / lr
        nw = np.ceil(H / d_w) - 1
        K = np.where(rec, np.maximum(np.minimum(nw, np.floor((r_w - rs) / d_w) + 1), 2), 2)
        f = np.where(rec, np.minimum(d_w, H - K * d_w) / H, 1.0)
        zc = K * d_w + f * H / 2
        Dloc = Dz * np.where(lr > 1e-9, lr / np.expm1(lr), 1.0) * np.exp(lr * zc / H)
        dw = f * H / np.sqrt(2 * Dloc / k)
    return rec, K, lr, dw


def section3(rs=None):
    rs = rows() if rs is None else rs
    assert rs, "no ROW lines in deephand_scale.log: run section 1 first"
    print(f"\n=== Section 3: summary of {len(rs)} scale runs, timing fit, multiwall section 4 re-timed ===")
    fill = lambda g: "n=%3d  %5.1f +- %4.1f %%" % (len(g), *(100 * np.array(M.se([q["ok"] for q in g]))) ) if g else "n=  0"
    print("P(rule hand fills) +- SE; unfinished; median s/run")
    for w in WS:
        g = [q for q in rs if q["w"] == w]
        if g:
            print(f"  w = {w:5}: {fill(g)}   unfinished {sum(q['Te'] is None for q in g)}   "
                  f"{np.median([q['sec'] for q in g]):.0f} s")
    print("by d/w of the deciding domain (rows) x w (columns: all, " + ", ".join(map(str, WS)) + ")")
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        g = [q for q in rs if lo <= q["dw"] < hi]
        print(f"  d/w [{lo:4.0f},{hi:4.0f}): {fill(g)}  |  " +
              "  ".join(fill([q for q in g if q["w"] == w]) for w in WS))
    for name, key, vals in (("K walls past onset", "K", None), ("layers n", "n", None), ("R/H", "R", None)):
        print(f"by {name} (all d/w | d/w >= 5)")
        for v in sorted({q[key] for q in rs}):
            g = [q for q in rs if q[key] == v]
            print(f"  {key} = {v:5}: {fill(g)}  |  {fill([q for q in g if q['dw'] >= 5])}")
    sub = [q for q in rs if q["Kd"] > 0]
    print(f"sub-onset wall in the stack (a wall below r*): {fill(sub)}  |  d/w >= 5: "
          f"{fill([q for q in sub if q['dw'] >= 5])}")
    big = [q for q in rs if q["dw"] >= 5]
    lo_w = [q for q in big if q["w"] <= 0.01]; hi_w = [q for q in big if q["w"] >= 0.02]
    (m1, e1), (m2, e2) = M.se([q["ok"] for q in lo_w]), M.se([q["ok"] for q in hi_w])
    print(f"d/w >= 5: w <= 0.01 {100 * m1:.1f} +- {100 * e1:.1f} vs w >= 0.02 {100 * m2:.1f} +- {100 * e2:.1f}; "
          f"difference / combined SE = {abs(m1 - m2) / max(np.hypot(e1, e2), 1e-12):.2f}")
    print("failures: w R n ln_rho u hands K Kd d/w final")
    for q in rs:
        if not q["ok"]:
            print(f"  {q['w']} {q['R']:5} {q['n']:2} {q['lr']:.2f} {q['u']:.2f} {q['hands']:10} {q['K']} {q['Kd']} "
                  f"{q['dw']:6.2f} {'unfin' if q['Te'] is None else q['fin']}")

    c, X, y = fit(rs)
    res = y - X @ c
    r = 10 ** y
    print(f"\nTiming, {len(y)} finished rule-held runs with T_c > 0: T_end/T_c min {r.min():.3f}, median "
          f"{np.median(r):.3f}, max {r.max():.3f}; RMS log10 (h = 1) {np.sqrt(np.mean(y ** 2)):.3f}")
    print(f"fit log10(T_end/T_c) = {c[0]:+.3f} {c[1]:+.3f} ln(rho) {c[2]:+.3f} log10(K+1) {c[3]:+.3f} log10(d/w); "
          f"residual RMS log10 {np.sqrt(np.mean(res ** 2)):.3f}, max |res| {np.abs(res).max():.3f} "
          f"(x{10 ** np.abs(res).max():.2f}); ratio to h = 1 RMS {np.sqrt(np.mean(res ** 2) / np.mean(y ** 2)):.2f}")
    for w in WS:
        m = np.array([q["w"] == w for q in rs if q["ok"] and q["Tc"] > 0])
        if m.any():
            print(f"  w = {w}: median T_end/T_c {np.median(r[m]):.3f}, residual RMS {np.sqrt(np.mean(res[m] ** 2)):.3f}")
    Kmax = max(q["K"] for q in rs)
    dlo, dhi = min(q["dw"] for q in rs), max(q["dw"] for q in rs)
    p = pbin(rs)
    print(f"h is clipped to the simulated K <= {Kmax} and d/w in [{dlo:.2f}, {dhi:.1f}]; d/w-bin P(fill) for the "
          f"ocean-scale variant: " + ", ".join(f"[{a:g},{b:g}) {v:.3f}" for a, b, v in zip(BINS[:-1], BINS[1:], p)))

    log = {}
    pat = re.compile(r"^(\S.*?)\s+[\d.]+ /\s+[\d.]+ /\s+[\d.]+ \|\s+[\d.]+ /\s+[\d.]+ /\s+[\d.]+\s+\|\|\s+(.*)$")
    for line in open(os.path.join(HERE, "multiwall.log")):
        m = pat.match(line)
        if m and m.group(1).strip() not in log:
            log[m.group(1).strip()] = tuple(re.findall(r"[\d.]+", m.group(2)))
    ref = M.run4(0, full=False)
    _, W = S.chains(0)
    out, bad = {}, []
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([0, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        o = {}
        for lid in (np.zeros(U.N, bool), salt):
            f = lambda **kw: M.multi(x, a[1], a[4], ex, lid, **kw)
            rec, K, lr, dw = geom(x, a[1], a[4], ex, lid)
            Xd = np.column_stack([np.ones(U.N), lr, np.log10(np.clip(K, 0, Kmax) + 1), np.log10(np.clip(dw, dlo, dhi))])
            h = np.where(rec, 10 ** (Xd @ c), 1.0)
            w1, r1, _ = f(scale=np.ones(U.N))
            w0, r0, _ = f()
            if not ((w1 == w0).all() and (r1 == r0).all()):
                bad.append((name, "h=1 arrays differ"))
            w3, r3, _ = f(rule="none")
            wh, rh, _ = f(scale=h)
            w10, r10, _ = f(scale=10.0)
            pb = np.where(dw >= 1, p[np.clip(np.digitize(dw, BINS) - 1, 0, len(p) - 1)], 0.0)
            pw = w3 + (wh & ~w3) * pb
            rw = np.maximum(r3 + (rh & ~r3) * pb, pw)
            dec = wh & ~w3
            for k, v in (("h1", S3.split(w1, r1)), ("x10", S3.split(w10, r10)), ("fit", S3.split(wh, rh)),
                         ("ocean", (100 * pw.mean(), 100 * (rw - pw).mean(), 100 * (1 - rw.mean())))):
                o.setdefault(k, []).append(v)
            o.setdefault("thin", []).append(100 * (dw[dec] < 5).mean() if dec.any() else np.nan)
            o.setdefault("lt1", []).append(100 * (dw[dec] < 1).mean() if dec.any() else np.nan)
            o.setdefault("hmed", []).append(float(np.median(h[rec])) if rec.any() else np.nan)
        if tuple(o["h1"]) != tuple(ref[name]["main"]):
            bad.append((name, "h=1 != multiwall run4 main"))
        if name in log and tuple(f"{v:.1f}" for s in o["h1"] for v in s) != log[name]:
            bad.append((name, "h=1 != multiwall.log", log[name]))
        out[name] = o
    print(f"\nregression h = 1 vs multiwall run4(0) main (arrays) and multiwall.log main (printed 0.1, "
          f"{len(log)} worlds parsed): {'EXACT' if not bad and len(log) == 15 else 'FAILED ' + str(bad)}")
    order = sorted(out, key=lambda n: -out[n]["h1"][1][0])
    fmt = lambda s: " | ".join(f"{a:5.1f} / {b:4.1f} / {c_:5.1f}" for a, b, c_ in s)
    print(f"\nP_pved / P_chance / P_none %, no lid | final, N = {U.N}, seed 0")
    print(f"{'world':18} {'multiwall main':>32} || {'multiwall x10':>32} || {'fitted h':>32} || ocean-scale (d/w-weighted)")
    for n in order:
        print(f"{n:18} " + " || ".join(fmt(out[n][k]) for k in ("h1", "x10", "fit", "ocean")))
    print("\nRule-decided P_pved draws (fitted h): % with deciding layer < 5 | < 1 local wall widths (no lid / final); "
          "median h over recursion draws")
    for n in order:
        q = out[n]
        print(f"{n:18} {q['thin'][0]:5.1f} / {q['thin'][1]:5.1f} | {q['lt1'][0]:5.1f} / {q['lt1'][1]:5.1f} | "
              f"{q['hmed'][0]:7.2f}")
    return out


def smoke():
    """Code paths + a cost probe (seconds per cell-step); nothing result-like is printed."""
    rng = np.random.default_rng(1)
    R, n, lr, rs, u, b, hands = config(0, rng, R=6.0, n=6, lr=2.0, u=0.6, alt=True)
    a = lr / 2   # full stretched layers: thickness x exp(-a (R - r_mid)) roughly equal
    th = -np.diff(np.concatenate([[R], b]))
    assert np.all(np.diff(th) > 0) and abs(R - b[0] - th[0]) < 1e-12, th          # thicker with depth
    s = (1 - np.exp(-a * (R - b))) / a
    assert np.allclose(np.diff(np.concatenate([[0], s])), s[0]), "layers not equal in stretched units"
    q = one(0.035, 0, rng, 1, tcap=50.0, R=6.0, n=6, lr=2.0, u=0.6, alt=True)
    assert q["Tc"] > 0 and q["K"] == 5 and q["Kd"] == 0, q
    q = one(0.035, 0, rng, 1, tcap=50.0, R=2.0, n=6, lr=3.0, u=0.6, alt=True)
    assert q["Kd"] > 0, q                                 # onset inside: sub-onset walls sink
    t0 = time.time()
    q = one(0.005, 0, rng, 1, tcap=2000.0, R=6.0, n=10, lr=2.0, u=0.6, alt=True)
    cost = (time.time() - t0) / (2000 / M.DT * int(round(1 / (0.005 / DZW))))
    U.N = 2000
    x, _, a_, _ = S.chains(0)[1]["Europa"]
    ex = {"eta": np.ones(U.N), "r": np.ones(U.N), "rho": 10 ** np.random.default_rng(0).uniform(0, 2, U.N)}
    rec, K, lr_, dw = geom(x, a_[1], a_[4], ex, np.zeros(U.N, bool))
    assert rec.dtype == bool and np.isfinite(dw[rec]).all()
    import contextlib, io
    fake = [dict(w=w, K=k, Kd=k % 2, n=3, R=6.0, lr=1.0 + k / 3, u=0.5, hands="+-+", dw=d, Te=2.0 * d, fin=1,
                 ok=d > 1, Tc=d * (1 + k / 10), sec=1.0) for w in WS for k in range(4) for d in (0.5, 2.0, 7.0, 20.0)]
    with contextlib.redirect_stdout(io.StringIO()):
        o = section3(fake)
    assert set(o["Europa"]) >= {"h1", "fit", "ocean", "x10"}
    return cost


if __name__ == "__main__":
    if SMOKE:
        c = smoke()
        print(f"smoke ok; cost probe {c * 1e9:.1f} ns per cell-step")
    else:
        {"1": lambda: section1(float(sys.argv[2]), *map(int, sys.argv[3:])), "2": section2, "3": section3}[sys.argv[1]]()
