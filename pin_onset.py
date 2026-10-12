"""pin_onset: where does a chirality wall get pinned in a SMOOTH exponential D_z(r) profile? (2026-10-12)

Why: sweep2.py step f applies the sharp criterion Eq. 20, |d ln sigma/dr| < 2/r with sigma ~ sqrt(D_z), i.e.
G = (1/2)|d ln D_z/dr| < 2/r, to a smooth exponential profile (bottom/top ratio rho, G = ln(rho)/(2H)): a layered
draw is pinned if ln(rho) >= eq20 * 4H/r_w. two_zone.log measured the onset only for a sharp two-zone step (ratio
1.02-1.05, 2-5x above the sharp criterion at the joint). A smooth profile was never simulated.

Model: two_zone.py's radial equation and grid (R=16, H=1, dz=0.002, dt=0.1, D_top=2.5e-3, same finite-volume
operator, implicit diffusion), da/dt = a - a^3 + h + (1/r^2) d/dr(r^2 D(r) da/dr), no flux at floor and ice, with
    D(r) = D_top * exp(2 G (R - r))      (D grows with depth, the pinning direction of sweep2's rho >= 1),
so |d ln sigma/dr| = G everywhere. Gradient strength x = G / (2/r0), r0 = start radius of the wall; x = 1 is Eq. 20.
+ hand above the wall (under the ice), - hand below; h > 0 favours +, i.e. pushes the wall down.
 kappa = measured critical x / 1 (the sharp criterion), so sweep2's eq20 multiplier = kappa.

Sections:
 (1) fate scan: x = 0.1 .. 20, 6 runs per point (start depth d0 ~ U(0.3, 0.7), initial noise 1e-3, seeds fixed);
     each run ends "over" (- hand gone: wall reached the floor), "reversed" (+ hand gone: wall went up to the ice)
     or "stalled" (neither by t_max = time to move 1.5 depths at the slowest predicted speed). Sign control: D
     DEcreasing with depth at |x| = 5, 20.
 (2) kappa: wall speed measured directly. After 10 time units of relaxation, the wall position (sub-cell zero
     crossing) is fitted linearly over a window of 0.1 r0 / D(r0) (displacement ~0.02-0.08). s = inward speed /
     (D(r_w) 2/r_w); sharp-interface theory (gradient flow, sigma ~ sqrt D, mobility sqrt D) gives s = 1 - x_w,
     x_w = G r_w / 2. kappa = root of a linear fit of s on x_w over x = 0.6 .. 1.4. 8 configurations, mean +- SE.
     Convergence: one configuration at dt/2 and dz/2.
 (3) radius / wall-width dependence: R = 2, 4, 16, 64 (H = 1) x D_top = 6.25e-4, 2.5e-3, 1e-2 (wall width
     sqrt(2D) = 0.035, 0.071, 0.14), 4 configurations each. Bias: h = +-3e-4, +-1e-3 at the base geometry, compared
     with the Eq. 18 shift s -> s + (3/sqrt 2) h r_w / (2 sqrt D(r_w)).
 (4) sweep2 mapping, sweep2.py NOT edited: its chains() and outcome() are imported and step f is evaluated with
     eq20 = kappa, kappa +- SE (and the min/max over section 3's cells), 200,000 draws, seed 0, same extra draws as
     sweep2.run. Regression: eq20 = 1 and 5 reproduce sweep2.log. Per world: no lid / final, three-way split, the
     % of layered (column) draws pinned, and the median x = ln(rho) r_w / (4H) over layered draws (the kappa
     that would unpin half of them).

PREDICTIONS (written 2026-10-12, before any result-producing run):
 (P1) Base geometry (R = 16, D_top = 2.5e-3, h = 0): kappa = 1.00 +- 0.05 (mean of 8 configurations), SE < 0.02.
      Reason: for a profile smooth on the wall-width scale the sharp criterion is exact (gradient flow); the 2-5x of
      two_zone.log is a sharp-joint artifact (the joint, w = 0.02, is narrower than the wall, ~0.07).
 (P2) The speed law holds: a linear fit of s on x_w over all base points has slope -1 +- 0.1 and intercept 1 +- 0.1.
 (P3) No stalls in a smooth profile: every run at x <= 0.9 ends "over"; every run at x >= 1.1 ends "reversed" (the
      wall rises to the ice and the DEEP hand fills the shell); zero "stalled" runs at any x. So above onset the
      wall is not pinned in place; it runs the wrong way.
 (P4) Sign control: with D decreasing with depth, all runs at x = 5 and 20 end "over".
 (P5) No radius / width dependence: kappa within 0.9-1.1 for every cell with R >= 4; within 0.8-1.2 at R = 2.
 (P6) Bias: measured kappa - 1 within 25 % of the Eq. 18 prediction for h = +-3e-4 and +-1e-3 (positive h raises
      kappa). For the real g (<= 3e-16) the shift is negligible (< 1e-3).
 (P7) Convergence: halving dt and dz moves kappa by < 0.01.
 (P8) sweep2 at eq20 = kappa (and kappa +- SE): every world's step-f share, no lid and final, is within 1.0 point of
      sweep2.log step f (Earth 39.9, Europa 9.3 / 7.1, Enceladus 3.9 / 0.0, Ganymede 14.0); the x5 sensitivity is
      not supported. Earth: 100.0 % of layered draws pinned.
"""
import os
import sys
import time

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

import two_zone as TZ

H, DZ, DT, D_TOP, R0 = TZ.H, TZ.DZ, TZ.DT, TZ.D_HI, TZ.R
C18 = 3 / np.sqrt(2)                       # Eq. 18: flat wall speed = C18 * h * sqrt(D)
SMOKE = os.environ.get("SMOKE") == "1"     # tiny smoke test: SMOKE=1 python pin_onset.py
XS_K = (0.6, 0.8, 1.0, 1.2, 1.4)


def setup(x, R, d0, Dt, dz, dt, sign):
    """two_zone.run's radial operator with D(r) = Dt exp(sign 2G (R - r)), G = x 2/r0."""
    n = int(round(H / dz))
    rc = R - H + (np.arange(n) + 0.5) * dz
    rf = R - H + np.arange(1, n) * dz
    r0 = R - d0
    G = x * 2 / r0
    D = lambda r: Dt * np.exp(sign * 2 * G * (R - r))
    Df = D(rf)
    up = Df * rf ** 2 / (dz ** 2 * rc[:-1] ** 2)
    dn = Df * rf ** 2 / (dz ** 2 * rc[1:] ** 2)
    main = np.zeros(n); main[:-1] -= up; main[1:] -= dn
    lu = splu((diags(np.ones(n)) - dt * diags([main, up, dn], [0, 1, -1])).tocsc())
    return rc, r0, G, D, lu


def wall(a, rc):
    """Sub-cell radius of the lowest + crossing (None if one hand is gone)."""
    j = np.flatnonzero(a > 0)
    if j.size == 0 or j[0] == 0:
        return None
    j = j[0]
    return rc[j - 1] + (rc[j] - rc[j - 1]) * (-a[j - 1]) / (a[j] - a[j - 1])


def start(rc, r0, D, seed):
    rng = np.random.default_rng(seed)
    return np.tanh((rc - r0) / np.sqrt(2 * D(r0))) + 1e-3 * rng.standard_normal(rc.size)


def speed(x, R=R0, d0=0.5, Dt=D_TOP, h=0.0, dz=DZ, dt=DT, seed=0):
    """(x_w, s, r_w): normalized inward wall speed s = v / (D(r_w) 2/r_w) at x_w = G r_w / 2."""
    rc, r0, G, D, lu = setup(x, R, d0, Dt, dz, dt, +1)
    a = start(rc, r0, D, seed)
    t_relax, Tw = 10.0, 0.1 * r0 / D(r0)
    n1, n2 = int(round(t_relax / dt)), int(round(Tw / dt))
    every = max(1, n2 // 200)
    ts, rw = [], []
    for k in range(n1 + n2):
        a = lu.solve(a + dt * (a - a ** 3 + h))
        if k >= n1 and (k - n1) % every == 0:
            ts.append((k + 1) * dt); rw.append(wall(a, rc))
    assert None not in rw, f"hand vanished during speed window (x={x}, R={R}, Dt={Dt})"
    v = -np.polyfit(ts, rw, 1)[0]
    rm = float(np.mean(rw))
    return G * rm / 2, v / (D(rm) * 2 / rm), rm, D(rm)


def root(xw, s):
    b1, b0 = np.polyfit(xw, s, 1)
    return -b0 / b1, b1, b0


def kappa(R=R0, Dt=D_TOP, h=0.0, d0=0.5, seed=0, **kw):
    pts = [speed(x, R, d0, Dt, h, seed=seed, **kw) for x in XS_K]
    xw, s, rm, Dw = map(np.array, zip(*pts))
    k, b1, b0 = root(xw, s)
    k_th = root(xw, 1 - xw + C18 * h * rm / (2 * np.sqrt(Dw)))[0]   # sharp theory + Eq. 18, same points
    return k, b1, b0, k_th, xw, s


def configs(n, seed):
    rng = np.random.default_rng(seed)
    return [(float(rng.uniform(0.3, 0.7)), int(rng.integers(1 << 30))) for _ in range(n)]


def fate(x, sign=+1, R=R0, d0=0.5, Dt=D_TOP, seed=0):
    rc, r0, G, D, lu = setup(x, R, d0, Dt, DZ, DT, sign)
    a = start(rc, r0, D, seed)
    t_max = 200 + 1.5 / (Dt * 2 / r0 * abs(1 - sign * x))
    t, k = 0.0, 0
    while t < t_max:
        a = lu.solve(a + DT * (a - a ** 3))
        t += DT; k += 1
        if k % 50 == 0:
            if (a > 0).all():
                return "over", t, R - r0
            if (a < 0).all():
                return "reversed", t, R - r0
    return "stalled", t, R - wall(a, rc)


def section1(n=6):
    print(f"(1) fate scan, R={R0}, H={H}, D_top={D_TOP}, dz={DZ}, dt={DT}; {n} runs per point, d0 ~ U(0.3,0.7)")
    xs = (0.5, 1.5, 5) if SMOKE else (0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.5, 2, 3, 5, 10, 20)
    rows = [(x, +1) for x in xs] + [(x, -1) for x in ((5,) if SMOKE else (5, 20))]
    out = {}
    for i, (x, sg) in enumerate(rows):
        t0 = time.time()
        res = [fate(x, sg, d0=d, seed=s) for d, s in configs(n, 100 + i)]
        c = {o: sum(r[0] == o for r in res) for o in ("over", "reversed", "stalled")}
        p = c["over"] / n
        tm = np.mean([r[1] for r in res])
        print(f"  x={x:5g} D {'grows' if sg > 0 else 'falls'} with depth: over {c['over']}/{n}  reversed {c['reversed']}/{n}  "
              f"stalled {c['stalled']}/{n}  P(over) {p:.2f} +- {np.sqrt(p * (1 - p) / n):.2f}  mean end time {tm:8.0f}  "
              f"[{time.time() - t0:.0f}s]", flush=True)
        out[(x, sg)] = c
    return out


def section2(n=8):
    print(f"\n(2) kappa from wall speed, R={R0}, D_top={D_TOP}, h=0; x = {XS_K}; {n} configurations d0 ~ U(0.3,0.7)")
    ks, allx, alls = [], [], []
    for d, s in configs(2 if SMOKE else n, 7):
        k, b1, b0, k_th, xw, sv = kappa(d0=d, seed=s)
        ks.append(k); allx += list(xw); alls += list(sv)
        print(f"  d0={d:.3f}: kappa {k:.4f} (sharp theory {k_th:.4f})  s at x_w: "
              + " ".join(f"{a:.3f}:{b:+.3f}" for a, b in zip(xw, sv)), flush=True)
    k, se = float(np.mean(ks)), float(np.std(ks, ddof=1) / np.sqrt(len(ks)))
    b1, b0 = np.polyfit(allx, alls, 1)
    print(f"  KAPPA = {k:.4f} +- {se:.4f} (SE, n={len(ks)});  pooled fit s = {b0:.4f} {b1:+.4f} x_w (theory 1 - x_w)")
    kh = kappa(d0=0.5, seed=0, dz=DZ / 2, dt=DT / 2)[0]
    k1 = kappa(d0=0.5, seed=0)[0]
    print(f"  convergence (d0=0.5): dz={DZ}, dt={DT}: {k1:.4f};  dz/2, dt/2: {kh:.4f};  shift {kh - k1:+.4f}", flush=True)
    return k, se, kh - k1


def section3(n=4):
    print(f"\n(3a) radius and wall width: kappa mean +- SE over {n} configurations (h=0)")
    Rs, Ds = ((4, 16), (2.5e-3,)) if SMOKE else ((2, 4, 16, 64), (6.25e-4, 2.5e-3, 1e-2))
    cells = {}
    for R in Rs:
        for Dt in Ds:
            t0 = time.time()
            ks = [kappa(R=R, Dt=Dt, d0=d, seed=s)[0] for d, s in configs(2 if SMOKE else n, 11)]
            cells[(R, Dt)] = (np.mean(ks), np.std(ks, ddof=1) / np.sqrt(len(ks)))
            print(f"  R={R:3g} (r/H~{R - 0.5:.1f})  D_top={Dt:.2e} (width {np.sqrt(2 * Dt):.3f}): kappa {cells[(R, Dt)][0]:.4f} "
                  f"+- {cells[(R, Dt)][1]:.4f}  [{time.time() - t0:.0f}s]", flush=True)
    print(f"\n(3b) bias h at R={R0}, D_top={D_TOP}: measured kappa vs sharp theory + Eq. 18 shift")
    bias = {}
    for h in ((1e-3,) if SMOKE else (-1e-3, -3e-4, 3e-4, 1e-3)):
        r = [kappa(h=h, d0=d, seed=s) for d, s in configs(2 if SMOKE else n, 13)]
        km, kt = np.mean([q[0] for q in r]), np.mean([q[3] for q in r])
        se = np.std([q[0] for q in r], ddof=1) / np.sqrt(len(r))
        bias[h] = (km, kt)
        print(f"  h={h:+.0e}: kappa {km:.4f} +- {se:.4f}; predicted {kt:.4f}; (kappa-1)/(pred-1) = {(km - 1) / (kt - 1):.3f}", flush=True)
    return cells, bias


def section4(k, se, band):
    import sweep2 as S
    import uncertainty as U
    if SMOKE:
        U.N = 2000
    mults = {"kappa-SE": k - se, "kappa": k, "kappa+SE": k + se, "band lo": band[0], "band hi": band[1], "x1 (log)": 1.0, "x5 (log)": 5.0}
    print(f"\n(4) sweep2 step f with eq20 = kappa (sweep2.py imported, not edited); N={U.N}, seed 0")
    print("  multipliers: " + ", ".join(f"{a} {b:.4f}" for a, b in mults.items()))
    _, W = S.chains(0)
    res = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([0, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        L, dz_rng = a[1], a[4]
        Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
        Dh, Hd = x["horizontal D_h"], U.depth(x["volume"], L)
        column = Hd * np.sqrt(Dh / Dz) > L
        r_w = L / np.pi - np.minimum(Hd, L * np.sqrt(Dz / Dh))
        xd = np.log(ex["rho"]) * r_w / (4 * Hd)          # gradient strength in units of Eq. 20
        o = {}
        for lab, m in mults.items():
            wr = [S.outcome(x, L, dz_rng, "f", ex, l, eq20=m) for l in (none, salt)]
            o[lab] = [(100 * w.mean(), 100 * (r & ~w).mean(), 100 * (1 - r.mean())) for w, r in wr]
            o[lab + " pinned"] = 100 * (column & (xd >= m)).sum() / column.sum()
        o["column"], o["xmed"] = 100 * column.mean(), float(np.median(xd[column]))
        res[name] = o
    if U.N == 200_000:
        logged = {"Early Earth ocean": ((39.9, 39.9), (40.0, 40.0)), "Europa": ((9.3, 7.1), (16.3, 12.3)),
                  "Enceladus": ((3.9, 0.0), (16.9, 0.0)), "Ganymede": ((14.0, 14.0), (16.2, 16.2)), "Miranda": ((12.4, 4.1), (30.1, 10.0))}
        for n, (f1, f5) in logged.items():
            for lab, ref in (("x1 (log)", f1), ("x5 (log)", f5)):
                got = (res[n][lab][0][0], res[n][lab][1][0])
                assert max(abs(np.subtract(got, ref))) <= 0.051, f"regression {n} {lab}: got {got}, sweep2.log {ref}"
        print("  regression: eq20 = 1 and 5 reproduce sweep2.log step f (5 worlds, within rounding)")
    order = sorted(res, key=lambda n: -res[n]["kappa"][1][0])
    print(f"\n  step-f share, P_pved %, no lid / final | layered % of draws, % of layered pinned at kappa, median x over layered")
    print(f"  {'world':18}" + "".join(f"{lab:>16}" for lab in mults) + "   layered  pinned@k  median x")
    for n in order:
        print(f"  {n:18}" + "".join(f"{res[n][lab][0][0]:7.1f} /{res[n][lab][1][0]:6.1f} " for lab in mults)
              + f"  {res[n]['column']:6.1f}  {res[n]['kappa pinned']:7.1f}  {res[n]['xmed']:9.3g}")
    print("\n  three-way split at kappa, %: P_pved / P_chance / P_none (no lid) | (final)")
    for n in order:
        (a1, b1, c1), (a2, b2, c2) = res[n]["kappa"]
        print(f"  {n:18} {a1:5.1f} / {b1:5.1f} / {c1:5.1f} | {a2:5.1f} / {b2:5.1f} / {c2:5.1f}")
    return res


def main():
    t0 = time.time()
    sec = sys.argv[1:] or ["1", "2", "3", "4"]
    if "1" in sec:
        section1(2 if SMOKE else 6)
    k = se = None
    if "2" in sec or "4" in sec:
        k, se, _ = section2()
    cells = None
    if "3" in sec:
        cells, _ = section3()
    if "4" in sec:
        band = (min(v[0] for v in cells.values()), max(v[0] for v in cells.values())) if cells else (k, k)
        section4(k, se, band)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
