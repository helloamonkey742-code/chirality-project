"""multiwall: which layer's hand fills a stack of chirality walls past the Eq. 20 onset? (2026-10-12)

Why: sweep3.py counts a reversed draw whose SECOND wall is also past onset ("recursion") as no single hand (P_none);
its bound (iii) instead gives that draw to the next layer down. Neither rule was simulated. pin_onset.log measured a
lone wall: kappa = 0.9999, and above onset (x = G r / 2 > 1) the wall runs UP to the ice, so the deep hand fills.
This script simulates 2-5 stacked walls and measures the rule, then puts the measured rule into sweep3.

Model: pin_onset's radial equation and finite-volume operator, da/dt = a - a^3 + (1/r^2) d/dr(r^2 D(r) da/dr), no
flux at floor and ice, h = 0, implicit diffusion (a numba Thomas solver; checked against pin_onset's splu solve in
the smoke test), explicit reaction, dz = 0.002, dt = 0.1. D_z rises with depth, ln D(r) = ln D_top + 2 int_r^R G,
G(r) = 2 x(r) / r, so the local Eq. 20 strength is x(r) and the sharp-interface wall speed is D (2/r)(x - 1) upward.
Layers are numbered top (under the ice) to bottom; each has its own hand (+-1) and its own x_i; x(r) steps between
layers with a tanh blend of width 0.03 (smooth on the wall-width scale). Initial state: a = the layer's hand with a
tanh wall (width sqrt(2 D)) at each boundary where hands differ, + 1e-3 noise. A run ends when one hand is gone.
tau_max = 5 x the longest predicted time (lone-wall transit integral or the sum over walls, below) + 100; a run still
holding two hands at tau_max is logged "unfinished".
Predicted times (all with the speed law D (2/r)|x - 1|, kappa = 1):
  T_deep = transit of the DEEPEST wall from its start to the ice (one wall width below it) through the real D(r):
           the "max" rule = the layers' crossing times summed bottom-up along ONE wall's path;
  T_sum  = the sum over all walls of each one's lone transit to the ice (walls take turns: the "sum" rule);
  T_conv = sweep3's convention for that wall: H r / (2 D_mean (x - 1)) (whole depth H, depth-mean D, x at the wall).
Sections (each a separate run, `nice -n 10 /opt/anaconda3/bin/python multiwall.py <n> >> multiwall.log`):
 (1) all past onset, R = 16, H = 1, D_top = 6.25e-4 (wall width 0.035): 3, 4, 5 and 6 layers (2-5 walls) x 6 hand
     patterns (alternating; top only differs; bottom only differs; sandwich; 2 random) x 7 x patterns (uniform 1.1,
     1.5, 3, 10; increasing 1.1 -> 10 with depth; decreasing; random from {1.1, 1.5, 3, 10}) x 2 replicates (equal
     layers; layer thickness jittered U(0.8, 1.2)) = 336 runs, seed 20261012. Rules checked on every run: deepest
     layer's hand; largest volume (r^2-weighted); top layer's hand. Volume fraction of each hand tracked every 2 time
     units.
 (2) onset INSIDE the stack, sweep3's own geometry: one constant G = 2 / r* (x = r / r*), R = 2, H = 1 (so x spans
     0.5-1 to 1-2 across the shell; at R = 16 walls near onset move too slowly to simulate), D_top = 2.5e-4. 3-6
     layers, r* at 25 / 50 / 75 % of the way down each layer, hands alternating and random: 108 runs. Walls above r*
     should rise, walls below sink. Times: T_up / T_down = transit integrals of the deepest rising wall to the ice and
     the shallowest sinking wall to the floor (real D(r)); T_c = the wrapper's closed form (constant D = depth mean,
     below); T_conv = sweep3 convention H r_K/(2 D (x_K - 1)) + (r_K - r_floor) r_K / (2 D).
 (3) convergence: section 1, 4 layers, alternating, uniform x = 1.5, equal layers, at (dt, dz) and (dt/2, dz/2).
 (4) sweep3 with the measured rule (sweep3.py, sweep2.py NOT edited; sweep3.outcome3 is called and only its
     "recursion" draws are reclassified). Same draws, seeds and extra draws as sweep3.run. Walls k = 1..n_w sit at
     r_k = r_w - (k - 1) d_w (wall k exists iff H > k d_w); r* = 4H / ln(rho); wall k is past onset iff r_k >= r*;
     K = number past onset (K >= 2 is exactly sweep3's recursion). Rule: the layer CONTAINING r* (layer K + 1; the
     bottom layer if every wall is past onset) fills the ocean; its fraction f = min(d_w, H - K d_w) / H; its hand is
     PVED-set if F Delta(f) >= 2 (sweep3's test on that layer's own N), chance-set otherwise. Time: the deepest
     rising wall's transit T_up = (r*/2D_z)[(R - r_K) + r* ln((R - r*)/(r_K - r*))] (constant D_z, constant G, the
     sharp-interface law integrated) plus t_sel = max(10/k, (L/A)^2/D_h, T_down), T_down = (r*/2D_z)[r* ln((r* -
     r_floor)/(r* - r_{K+1})) - (r_{K+1} - r_floor)] for the walls below r* to sink (0 if K = n_w). Timely iff
     T_up + t_sel <= tau (sequential, conservative, as sweep3). Sensitivities: cascade time x10 (T_up and T_down);
     "sum" rule (T_up summed over walls 1..K). Regression: rule "recursion = P_none" must reproduce every sweep3.log
     main three-way share exactly (to the printed 0.1). Seeds 0-3.
  ASSUMED (as sweep3): D_z drawn = the profile's mean sets every time; depth fractions are volume fractions; hands of
  neighbouring layers are taken to differ (walls everywhere: conservative for time).

PREDICTIONS (written 2026-10-12, before any result-producing run):
 (P1) Section 1: the deepest layer's hand fills the shell in 336/336 runs; 0 unfinished. "Largest volume" and "top
      layer" each fail in >= 25 % of runs.
 (P2) Section 1: the deepest layer's hand never loses more than 0.01 of the volume between checkpoints after t = 20
      (it only grows) in every run.
 (P3) Section 1 timing: T_end / T_deep lies in [0.5, 1.15] in every run, median in [0.85, 1.10]; over runs with
      >= 2 walls the RMS log10 error of the sum rule (T_sum) is more than 2x that of the max rule (T_deep). So the
      walls move together, not in turn.
 (P4) Section 1: T_end / T_conv lies in [0.1, 3] in every run (inside sweep3's x10 sensitivity).
 (P5) Section 2: the hand of the layer containing r* fills the shell in 108/108 runs; 0 unfinished.
 (P6) Section 2: T_end / max(T_up, T_down) lies in [0.6, 1.3] in every run.
 (P7) Section 2: T_end / T_c (the wrapper's closed form, T_up + T_down with D = depth mean) lies in [0.2, 5] in
      every run.
 (P8) Section 2: sweep3's convention overestimates near onset: T_conv > T_end in >= 80 % of runs.
 (P9) Section 3: halving dt and dz changes T_end by < 2 % and keeps the final hand.
 (P10) Section 4 regression: rule "recursion = P_none" reproduces every sweep3.log main three-way share (15 worlds,
      no lid and final) to the printed 0.1, and K >= 2 equals sweep3's recursion mask except for <= 10 draws per world
      (float rounding).
 (P11) Section 4: multiwall P_pved >= sweep3 main (by construction) and <= sweep3 bound (iii) + 0.5 in every world,
      no lid and final.
 (P12) Section 4, no lid: every icy moon gains >= 3 points over sweep3 main and ends <= 6 points below its bound
      (iii); Early Earth 75-86 %.
 (P13) Section 4: cascade time x10 lowers every world by <= 3 points (no lid and final).
 (P14) Section 4: P_chance <= 5 points in every world (no lid and final).
 (P15) Section 4, final: Ganymede > Europa > Enceladus.
 (P16) Section 4: seeds 1-3 reproduce every multiwall share (P_pved / P_chance / P_none, no lid and final) within
      0.5 points of seed 0.
"""
import os
import re
import sys

import numpy as np
from numba import njit
from scipy.special import gammaln

import sweep2 as S
import sweep3 as S3
import uncertainty as U
from ocean import delta_at

SMOKE = os.environ.get("SMOKE") == "1"     # SMOKE=1 python multiwall.py: solver + code paths, prints no results
SEED, DT, DZ, BLEND, EVERY = 20261012, 0.1, 0.002, 0.03, 20
XS = (1.1, 1.5, 3.0, 10.0)


@njit(cache=False)
def evolve(a, lo, cp, iden, dt, nmax, every, w, traj):
    """Implicit-diffusion / explicit-reaction steps (Thomas solve); stores the + hand's volume fraction every
    `every` steps; returns (steps done, checkpoints stored), stopping as soon as one hand is gone."""
    n = a.size
    d = np.empty(n)
    k = 0
    for step in range(nmax):
        for i in range(n):
            d[i] = a[i] + dt * (a[i] - a[i] ** 3)
        d[0] = d[0] * iden[0]
        for i in range(1, n):
            d[i] = (d[i] - lo[i] * d[i - 1]) * iden[i]
        a[n - 1] = d[n - 1]
        for i in range(n - 2, -1, -1):
            a[i] = d[i] - cp[i] * a[i + 1]
        if (step + 1) % every == 0:
            pos, npos = 0.0, 0
            for i in range(n):
                if a[i] > 0:
                    pos += w[i]
                    npos += 1
            traj[k] = pos
            k += 1
            if npos == 0 or npos == n:
                return step + 1, k
    return nmax, k


def shell(R, H, Dt, xfun, dz, dt):
    """Grid, D on a fine grid (ln D = ln Dt + 2 int_r^R 2x/r), and the factored (I - dt A) of pin_onset.setup."""
    n = int(round(H / dz))
    rc = R - H + (np.arange(n) + 0.5) * dz
    rf = R - H + np.arange(1, n) * dz
    fine = np.linspace(R - H, R, 20001)
    G = 2 * xfun(fine) / fine
    lnD = np.log(Dt) + 2 * np.concatenate([np.cumsum(((G[1:] + G[:-1]) / 2 * np.diff(fine))[::-1])[::-1], [0.0]])
    D = lambda r: np.exp(np.interp(r, fine, lnD))
    Df = D(rf)
    up = Df * rf ** 2 / (dz ** 2 * rc[:-1] ** 2)
    dn = Df * rf ** 2 / (dz ** 2 * rc[1:] ** 2)
    main = np.zeros(n); main[:-1] -= up; main[1:] -= dn
    b, c, lo = 1 - dt * main, np.append(-dt * up, 0.0), np.concatenate([[0.0], -dt * dn])
    cp, iden = np.empty(n), np.empty(n)
    iden[0] = 1 / b[0]; cp[0] = c[0] * iden[0]
    for i in range(1, n):
        iden[i] = 1 / (b[i] - lo[i] * cp[i - 1]); cp[i] = c[i] * iden[i]
    return rc, D, (lo, cp, iden), (b, up, dn)


def transit(r0, R, H, D, xfun, up):
    """Lone-wall time from r0 to one wall width below the ice (up) or above the floor, speed D (2/r)|x - 1|."""
    end = R - np.sqrt(2 * D(R)) if up else R - H + np.sqrt(2 * D(R - H))
    r = np.linspace(r0, end, 4001)
    return abs(np.trapezoid(r / (2 * D(r) * np.abs(xfun(r) - 1)), r))


def simulate(R, H, Dt, xfun, bounds, hands, seed, tmax, dz=DZ, dt=DT):
    """(T_end or None, final hand or 0, + hand volume fraction trajectory, D)."""
    rc, D, fac, _ = shell(R, H, Dt, xfun, dz, dt)
    s = np.array(hands, float)
    a = np.full(rc.size, s[-1])
    for b, (hi, lo_) in zip(bounds, zip(s[:-1], s[1:])):
        if hi != lo_:
            a *= -np.tanh((rc - b) / np.sqrt(2 * D(b)))
    a += 1e-3 * np.random.default_rng(seed).standard_normal(rc.size)
    w = rc ** 2 / np.sum(rc ** 2)
    every = max(1, int(round(EVERY * DT / dt)))
    nmax = int(np.ceil(tmax / dt))
    traj = np.empty(nmax // every + 2)
    done, k = evolve(a, *fac, dt, nmax, every, w, traj)
    traj = traj[:k]
    if np.all(a > 0) or np.all(a < 0):
        return done * dt, int(np.sign(a[0])), traj, D
    return None, 0, traj, D


def xstep(xs, bounds):
    """x(r): x_i inside layer i (top to bottom), tanh blends of width BLEND at the boundaries."""
    def f(r):
        out = np.full(np.shape(r), xs[-1], float)
        for b, xa, xb in zip(bounds, xs[:-1], xs[1:]):
            out += (xa - xb) * 0.5 * (1 + np.tanh((np.asarray(r) - b) / BLEND))
        return out
    return f


def wrapper_times(R, H, rs, rK, rK1, D):
    """Section 4's closed-form T_up (deepest rising wall) and T_down (shallowest sinking wall); None = no such wall."""
    up = 0.0 if rK is None else rs / (2 * D) * ((R - rK) + rs * np.log((R - rs) / (rK - rs)))
    down = 0.0 if rK1 is None else rs / (2 * D) * (rs * np.log((rs - (R - H)) / (rs - rK1)) - (rK1 - (R - H)))
    return up, down


def hand_patterns(n, rng):
    p = {"alt": [(-1) ** i for i in range(n)], "top": [1] + [-1] * (n - 1), "bottom": [1] * (n - 1) + [-1],
         "sandwich": [1] + [-1] * (n - 2) + [1]}
    for name in ("rand1", "rand2"):
        h = rng.choice([-1, 1], n)
        while abs(h.sum()) == n:
            h = rng.choice([-1, 1], n)
        p[name] = list(h)
    return p


def layers(n, R, H, rng, jitter):
    th = rng.uniform(0.8, 1.2, n) if jitter else np.ones(n)
    th = H * th / th.sum()
    return R - np.cumsum(th)[:-1], th


def se(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else 0.0


def section1(smoke=False):
    R, H, Dt = 16.0, 1.0, 6.25e-4
    rng = np.random.default_rng(SEED)
    rows, idx = [], 0
    ns = (3,) if smoke else (3, 4, 5, 6)
    for n in ns:
        hp = hand_patterns(n, rng)
        xp = {"u1.1": [1.1] * n, "u1.5": [1.5] * n, "u3": [3.0] * n, "u10": [10.0] * n,
              "inc": list(np.array(XS)[np.round(np.linspace(0, 3, n)).astype(int)]),
              "rand": list(rng.choice(XS, n))}
        xp["dec"] = xp["inc"][::-1]
        for rep in (0, 1):
            bounds, th = layers(n, R, H, rng, rep == 1)
            for hn, hands in hp.items():
                for xn, xs in xp.items():
                    if smoke and (hn != "alt" or xn != "u3" or rep):
                        continue
                    idx += 1
                    xf = xstep(xs, bounds)
                    walls = [b for b, h1, h2 in zip(bounds, hands[:-1], hands[1:]) if h1 != h2]
                    _, D, _, _ = shell(R, H, Dt, xf, DZ, DT)
                    T = [transit(b, R, H, D, xf, True) for b in walls]
                    Tdeep, Tsum = T[-1], sum(T)
                    rr = np.linspace(R - H, R, 4001)
                    Tconv = H * walls[-1] / (2 * D(rr).mean() * (xf(walls[-1]) - 1))
                    tmax = 5 * max(Tdeep, Tsum) + 100
                    if smoke:
                        tmax = 50.0
                    Te, fin, traj, _ = simulate(R, H, Dt, xf, bounds, hands, SEED + idx, tmax)
                    vol = np.diff(np.concatenate([[R - H], bounds[::-1], [R]]) ** 3)[::-1]   # top to bottom
                    vsign = np.sign(np.dot(hands, vol))
                    fb = traj if hands[-1] > 0 else 1 - traj
                    k0 = int(20 / (EVERY * DT))
                    drop = float(np.max(-np.diff(fb[k0:]), initial=0.0))
                    qs = [fb[min(len(fb) - 1, int(q * len(fb)))] for q in (0.0, 0.25, 0.5, 0.75)] if len(fb) else [np.nan] * 4
                    rows.append(dict(n=n, hand=hn, x=xn, rep=rep, nw=len(walls), Te=Te, fin=fin, deep=hands[-1],
                                     top=hands[0], vol=vsign, Tdeep=Tdeep, Tsum=Tsum, Tconv=Tconv, tmax=tmax,
                                     drop=drop, qs=qs))
    if smoke:
        return rows
    print("\n=== Section 1: all walls past onset (R = 16, H = 1, D_top = 6.25e-4, dt = 0.1, dz = 0.002), 336 runs ===")
    print("n hands     x     rep walls  final  T_end    T_deep   T_sum    T_conv   end/deep end/sum end/conv "
          "maxdrop  deep-hand fraction at 0/25/50/75 % of T_end")
    for q in rows:
        Te = q["Te"]
        fin = "unfin" if Te is None else ("DEEP" if q["fin"] == q["deep"] else "other")
        r = lambda T: f"{Te / T:7.3f}" if Te else "    nan"
        print(f"{q['n']} {q['hand']:9} {q['x']:5} {q['rep']}   {q['nw']}    {fin:5} {Te or q['tmax']:8.0f} "
              f"{q['Tdeep']:8.0f} {q['Tsum']:8.0f} {q['Tconv']:8.0f} {r(q['Tdeep'])}  {r(q['Tsum'])} {r(q['Tconv'])} "
              f"{q['drop']:7.4f}  " + " ".join(f"{v:.3f}" for v in q["qs"]))
    summarize1(rows)
    return rows


def summarize1(rows):
    ok = [q for q in rows if q["Te"] is not None]
    print(f"\nSection 1 summary: {len(rows)} runs, {len(rows) - len(ok)} unfinished at tau_max")
    for rule in ("deep", "vol", "top"):
        hit = [q["fin"] == q[rule] for q in rows]
        m, e = se(hit)
        print(f"  rule '{rule}' hand fills: {sum(hit)}/{len(rows)} = {100 * m:.1f} +- {100 * e:.1f} %")
    print("  configuration groups: P(deepest hand fills) +- SE | T_end/T_deep mean +- SE | T_end/T_sum | T_end/T_conv")
    for key in ("n", "hand", "x"):
        for v in sorted({q[key] for q in rows}, key=str):
            g = [q for q in ok if q[key] == v]
            gd = [q["fin"] == q["deep"] for q in rows if q[key] == v]
            pd, pe = se(gd)
            a = [se([q["Te"] / q[T] for q in g]) for T in ("Tdeep", "Tsum", "Tconv")]
            print(f"    {key}={str(v):9} n={len(gd):3}  {100 * pd:5.1f} +- {100 * pe:4.1f} %  " +
                  "  ".join(f"{m:6.3f} +- {e:5.3f}" for m, e in a))
    rd = np.array([q["Te"] / q["Tdeep"] for q in ok])
    print(f"  T_end/T_deep: min {rd.min():.3f}, median {np.median(rd):.3f}, max {rd.max():.3f}")
    rc = np.array([q["Te"] / q["Tconv"] for q in ok])
    print(f"  T_end/T_conv: min {rc.min():.3f}, median {np.median(rc):.3f}, max {rc.max():.3f}")
    multi = [q for q in ok if q["nw"] >= 2]
    for T in ("Tdeep", "Tsum", "Tconv"):
        e = np.sqrt(np.mean([np.log10(q["Te"] / q[T]) ** 2 for q in multi]))
        print(f"  RMS log10(T_end/{T}) over {len(multi)} runs with >= 2 walls: {e:.3f}")
    print(f"  largest loss of the deep hand's volume fraction between checkpoints (t > 20): "
          f"{max(q['drop'] for q in rows):.4f}")


def section2(smoke=False):
    R, H, Dt = 2.0, 1.0, 2.5e-4
    rng = np.random.default_rng(SEED + 1)
    rows, idx = [], 0
    print("\n=== Section 2: onset inside the stack, constant G = 2/r* (x = r/r*), R = 2, H = 1, D_top = 2.5e-4 ===")
    print("n hands j frac  r*     walls final        T_end    T_up   T_down   T_c      T_conv  end/max(up,down) "
          "end/T_c end/conv")
    for n in ((3,) if smoke else (3, 4, 5, 6)):
        bounds, _ = layers(n, R, H, rng, False)
        top = np.concatenate([[R], bounds]); bot = np.concatenate([bounds, [R - H]])
        hp = hand_patterns(n, rng)
        for hn in ("alt", "rand1"):
            hands = hp[hn]
            for j in range(n):
                for fr in (0.25, 0.5, 0.75):
                    if smoke and (j or fr != 0.5):
                        continue
                    idx += 1
                    rs = top[j] - fr * (top[j] - bot[j])
                    xf = lambda r, rs=rs: np.asarray(r) / rs
                    walls = [b for b, h1, h2 in zip(bounds, hands[:-1], hands[1:]) if h1 != h2]
                    upw = [b for b in walls if b > rs]; dnw = [b for b in walls if b < rs]
                    _, D, _, _ = shell(R, H, Dt, xf, DZ, DT)
                    Tup = transit(min(upw), R, H, D, xf, True) if upw else 0.0
                    Tdn = transit(max(dnw), R, H, D, xf, False) if dnw else 0.0
                    rr = np.linspace(R - H, R, 4001)
                    Dm = D(rr).mean()
                    rK = min(upw) if upw else None
                    Tc = sum(wrapper_times(R, H, rs, rK, max(dnw) if dnw else None, Dm))
                    Tconv = (H * rK / (2 * Dm * (rK / rs - 1)) if rK else 0.0) + \
                            ((rK or R) - (R - H)) * (rK or R) / (2 * Dm) * bool(dnw)
                    tmax = 50.0 if smoke else 5 * max(Tup, Tdn) + 100
                    Te, fin, traj, _ = simulate(R, H, Dt, xf, bounds, hands, SEED + 1000 + idx, tmax)
                    want = hands[j]
                    rows.append(dict(n=n, hand=hn, j=j, fr=fr, Te=Te, ok=fin == want, Tm=max(Tup, Tdn), Tc=Tc,
                                     Tconv=Tconv))
                    if not smoke:
                        res = "unfin" if Te is None else ("LAYER(r*)" if fin == want else "other")
                        f = lambda T: f"{Te / T:7.3f}" if (Te and T) else "    nan"
                        print(f"{n} {hn:5} {j} {fr:4.2f} {rs:6.3f}  {len(walls)}    {res:9} {Te or tmax:8.0f} "
                              f"{Tup:8.0f} {Tdn:8.0f} {Tc:8.0f} {Tconv:8.0f}  {f(max(Tup, Tdn))}          "
                              f"{f(Tc)} {f(Tconv)}")
    if smoke:
        return rows
    ok = [q for q in rows if q["Te"] is not None]
    m, e = se([q["ok"] for q in rows])
    print(f"\nSection 2 summary: {len(rows)} runs, {len(rows) - len(ok)} unfinished; layer containing r* fills: "
          f"{sum(q['ok'] for q in rows)}/{len(rows)} = {100 * m:.1f} +- {100 * e:.1f} %")
    for name, key in (("T_end/max(T_up,T_down)", "Tm"), ("T_end/T_c", "Tc"), ("T_end/T_conv", "Tconv")):
        v = np.array([q["Te"] / q[key] for q in ok if q[key] > 0])
        mm, ee = se(v)
        print(f"  {name}: min {v.min():.3f}, median {np.median(v):.3f}, max {v.max():.3f}, mean {mm:.3f} +- {ee:.3f}")
    over = np.mean([q["Tconv"] > q["Te"] for q in ok])
    print(f"  T_conv > T_end in {100 * over:.1f} % of finished runs")
    return rows


def section3():
    R, H, Dt, n = 16.0, 1.0, 6.25e-4, 4
    bounds = R - H * np.arange(1, n) / n
    hands, xf = [1, -1, 1, -1], xstep([1.5] * n, bounds)
    print("\n=== Section 3: convergence, section 1 point (4 layers, alternating, x = 1.5, equal layers) ===")
    out = []
    for dz, dt in ((DZ, DT), (DZ / 2, DT / 2)):
        Te, fin, _, _ = simulate(R, H, Dt, xf, bounds, hands, SEED, 1e6, dz=dz, dt=dt)
        out.append((Te, fin))
        print(f"  dz = {dz:.4f}, dt = {dt:.3f}: T_end = {Te:.1f}, final hand {fin:+d} (deepest layer {hands[-1]:+d})")
    print(f"  change in T_end: {100 * (out[1][0] / out[0][0] - 1):+.2f} %; same final hand: {out[0][1] == out[1][1]}")


# ---------------------------------------------------------------- section 4: sweep3 with the measured rule
def multi(x, L, dz_rng, ex, lid, rule="max", scale=1.0):
    """(P_pved mask, one-hand mask, counts): sweep3.outcome3 with its recursion draws reclassified by the rule."""
    w3, r3, _ = S3.outcome3(x, L, dz_rng, ex, lid)
    if rule == "none":
        return w3, r3, {}
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    k, H, lr = k2 * c, U.depth(V, L), np.log(ex["rho"])
    R = L / np.pi
    column = H * np.sqrt(Dh / Dz) > L
    d_w = np.minimum(H, L * np.sqrt(Dz / Dh))
    r_w = R - d_w
    rev = column & (lr >= 4 * H / r_w) & (k * t >= 10) & ~lid
    rec = rev & ((H - d_w) * np.sqrt(Dh / Dz) > L) & (lr >= 4 * H / (r_w - d_w))   # sweep3's expression
    with np.errstate(all="ignore"):
        rs = 4 * H / lr
        nw = np.ceil(H / d_w) - 1
        K0 = np.minimum(nw, np.floor((r_w - rs) / d_w) + 1)
        mism = int((rec != (rev & (K0 >= 2) & (nw >= 2))).sum())
        K = np.where(rec, np.maximum(K0, 2), 2)            # the <= 10 rounding mismatches (P10) take K = 2
        rK = r_w - (K - 1) * d_w
        rK1 = rK - d_w
        up = rs / (2 * Dz) * ((R - rK) + rs * np.log((R - rs) / (rK - rs)))
        if rule == "sum":
            e = r_w - rs
            slog = K * np.log(d_w) + gammaln(e / d_w + 1) - gammaln(e / d_w - K + 1)   # sum_k ln(r_k - r*)
            up = rs / (2 * Dz) * (K * (R - r_w) + d_w * K * (K - 1) / 2 + rs * (K * np.log(R - rs) - slog))
        down = np.where(K < nw, rs / (2 * Dz) * (rs * np.log((rs - (R - H)) / (rs - rK1)) - (rK1 - (R - H))), 0.0)
        t_sel = np.maximum.reduce([10 / k, (L / U.A) ** 2 / Dh, scale * down])
        timely = rec & (scale * up + t_sel <= t)
        f = np.where(rec, np.minimum(d_w, H - K * d_w) / H, 1.0)
        strong = F * delta_at(c, g, k2, V * ex["r"], t, f) >= 2
    cnt = {"rec": 100 * rec.mean(), "timely": 100 * timely.sum() / max(rec.sum(), 1), "mism": mism,
           "bottom": 100 * (rec & (K == nw)).sum() / max(rec.sum(), 1),
           "Kmed": float(np.median(K[rec])) if rec.any() else np.nan}
    return w3 | (timely & strong), r3 | timely, cnt


def run4(seed, full=True):
    _, W = S.chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        f = lambda lid, **kw: multi(x, a[1], a[4], ex, lid, **kw)
        o = {"main": [S3.split(*f(l)[:2]) for l in (none, salt)]}
        if full:
            o["cnt"] = f(none)[2]
            o["none"] = [S3.split(*f(l, rule="none")[:2]) for l in (none, salt)]
            o["x10"] = [S3.split(*f(l, scale=10.0)[:2]) for l in (none, salt)]
            o["sum"] = [S3.split(*f(l, rule="sum")[:2]) for l in (none, salt)]
            o["bound"] = [S3.split(*S3.outcome3(x, a[1], a[4], ex, l, veto=False)[:2]) for l in (none, salt)]
        out[name] = o
    return out


def sweep3_log():
    pat = re.compile(r"^(\S.*?)\s+([\d.]+) /\s+([\d.]+) /\s+([\d.]+) \|\s+([\d.]+) /\s+([\d.]+) /\s+([\d.]+)\s+\[step f")
    rows = {}
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep3.log")):
        m = pat.match(line)
        if m:
            rows[m.group(1).strip()] = m.groups()[1:]
    assert len(rows) == 15, rows
    return rows


def section4():
    res = run4(0)
    log = sweep3_log()
    bad = [(n, log[n], res[n]["none"]) for n in res
           if tuple(f"{v:.1f}" for s in res[n]["none"] for v in s) != tuple(log[n])]
    mism = {n: res[n]["cnt"]["mism"] for n in res}
    print("\n=== Section 4: sweep3 with the measured rule (sweep3.py, sweep2.py imported, not edited), N = "
          f"{U.N}, seed 0 ===")
    print(f"regression, rule 'recursion = P_none' vs sweep3.log main three-way split (15 worlds, no lid and final, "
          f"printed 0.1): {'EXACT' if not bad else 'FAILED ' + str(bad)}")
    print(f"K >= 2 vs sweep3's recursion mask, mismatched draws per world: {mism}")
    order = sorted(res, key=lambda n: -res[n]["main"][1][0])
    print("\nThree-way split, %: P_pved / P_chance / P_none, no lid | final")
    print(f"{'world':18} {'sweep3 main':>22} | {'multiwall (measured rule)':>36}")
    for n in order:
        s3, m = res[n]["none"], res[n]["main"]
        print(f"{n:18} " + " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in s3) + "  ||  "
              + " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in m))
    print("\nP_pved %, no lid / final: sweep3 main | sweep3 bound (iii) | multiwall | cascade x10 | sum rule")
    for n in order:
        print(f"{n:18} " + " | ".join(f"{res[n][k][0][0]:5.1f} / {res[n][k][1][0]:5.1f}"
                                      for k in ("none", "bound", "main", "x10", "sum")))
    print("\nRecursion draws (no lid): % of all draws | % of them that fill in time | % decided by the bottom layer "
          "(every wall past onset) | median K (walls past onset)")
    for n in order:
        q = res[n]["cnt"]
        print(f"{n:18} {q['rec']:6.1f} | {q['timely']:6.1f} | {q['bottom']:6.1f} | {q['Kmed']:8.0f}")
    print("\nSeed stability, multiwall P_pved / P_chance / P_none (no lid and final), max |seed k - seed 0| over worlds:")
    for s in (1, 2, 3):
        r = run4(s, full=False)
        d = max(abs(r[n]["main"][l][j] - res[n]["main"][l][j]) for n in res for l in (0, 1) for j in range(3))
        print(f"  seed {s}: max diff {d:.2f}   " + ", ".join(f"{n} {r[n]['main'][0][0]:.1f}/{r[n]['main'][1][0]:.1f}"
                                                         for n in r))


def smoke():
    """Numba solver == pin_onset-style splu solve; every section's code path runs (tiny, nothing printed)."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import splu
    R, H, Dt = 16.0, 1.0, 6.25e-4
    bounds = np.array([15.6, 15.3])
    xf = xstep([1.5, 3.0, 10.0], bounds)
    rc, D, fac, (b, up, dn) = shell(R, H, Dt, xf, DZ, DT)
    lu = splu((diags(b) + diags(-DT * up, 1) + diags(-DT * dn, -1)).tocsc())
    a0 = np.tanh((rc - 15.5) / 0.05) + 1e-3 * np.random.default_rng(0).standard_normal(rc.size)
    a1, a2 = a0.copy(), a0.copy()
    for _ in range(200):
        a1 = lu.solve(a1 + DT * (a1 - a1 ** 3))
    evolve(a2, *fac, DT, 200, 10 ** 9, rc ** 2, np.empty(1))
    assert np.max(np.abs(a1 - a2)) < 1e-10, np.max(np.abs(a1 - a2))
    rK, rs = 1.6, 1.5   # closed form == numeric integral of r r*/(2 D (r - r*)) with constant D
    r = np.linspace(rK, 2.0, 200001)
    num = np.trapezoid(r * rs / (2 * 1e-3 * (r - rs)), r)
    assert abs(wrapper_times(2.0, 1.0, rs, rK, None, 1e-3)[0] / num - 1) < 1e-4
    section1(smoke=True)
    section2(smoke=True)
    U.N = 2000
    run4(0)
    sweep3_log()
    print("smoke ok")


if __name__ == "__main__":
    if SMOKE:
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            smoke()
        print("smoke ok")
    else:
        {"1": section1, "2": section2, "3": section3, "4": section4}[sys.argv[1]]()
