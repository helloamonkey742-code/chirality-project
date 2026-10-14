"""kcheck.py: is the draft's amplification test "k2 c > 2 k_rac" too lenient?  Evaluate persist.py's stricter factors at
the paper's own parameters and sampling (2026-10-13).

WHERE THE TEST LIVES.  sweep2.py:133 (step d: ok &= k2*c > 2*krac, Ala 9.13e-15 /s), reached by multiwall.py/sweep3.py through
S.outcome for every draw that is not a reversed-wall (recursion) draw; sweep3.outcome3 line 115-116 (deep layer, sensitivity ii
only); rho2side.py:120; rerun_base.py:819-831 (H27 "cut on every draw").  The main multiwall pipeline does NOT apply it to
the reversed-wall draws that are counted late (sweep3 docstring: "the main run does NOT test k > 2 k_rac on the deep layer").
Its origin is review.py R7 (normal form, racemization r in units of k, amplifier rate lam_max = 1 k: "if 2r > lam_max the choice
stops"), i.e. the test is the normal form's own threshold with lam_max taken to be one copying rate k.

THE QUESTION.  The draft's network (Eq. 2, dimensionless k = kd = ki = f = 1, k0 = 1e-3, a0 ramped 0.8 -> 1.6) has, with
racemization r each way, du/dt = (k a - kd - 2r) u + g k a s for u = L - D.  The test is the kd -> 0 limit of the linear
condition k a > kd + 2r: necessary, not sufficient.  In the network kd = k at the crossing (the draft's time unit is 1/k with
k = k2 c, kd = k) and the racemic state's own feed a_rac < a0 (it eats the feed), so the real threshold is
2 r_c = a_rac(a0) - 1 (in kd units), far below the normal form's lam_max = 1.  persist.py found r_c = 0.08698 kd (pure state
exists), r(eta*>=0.9) = 0.02307, r(eta*>=0.99) = 0.00127 and quoted them as 9.2x and 34.7x "k a0 / (2 r)".  That quote used
k a0 = 1.6 k.  If the draw's k2 c is the network's k = kd (the draft's identification), the factors are 1/(2 r) = 5.75x, 21.7x,
394x.  Both conventions are run below (the 1.6 ratio is the ramp end over the crossing).

MODEL FOR THE AMPLIFICATION CHECKS.  persist.py's rhs (reused, imported): same Eq. 2 plus racemization r (-r(L-D), +r(L-D)).
a0 and r are in units of a_c and k.  Draw-level checks use the committed pipeline: multiwall.multi (sweep2 step f -> sweep3
recount -> deepest-hand cascade) on sweep2.chains seed 0, N = 200,000, the seedfix stream for Enceladus ([s, 0, 1], installed
by g_wide.py on import), G_NEW = (4.4e-18, 2.64e-16) as in the draft's Table 3.  A stricter test is an extra mask
k2*c > M * 2 * K_RAC on BOTH the win mask and the one-hand mask of EVERY draw (including the recursion draws the committed
pipeline never tests); P_pved, P_chance, P_none are recomputed from the masked arrays (sweep3.split).  M = 1 on every draw is the
draft's own test applied everywhere.  Levels M:
   1        the draft's test
   1.667    hard necessary bound at a0 = 1.6 ignoring feed depletion: k a0 > kd + 2r  (r < 0.3 kd, M = 1/(2*0.3))
   M_ex     1/(2 r_c) = pure branch exists                       (k2c := kd)
   1.6 M_ex same, ramp-end convention (persist.py's 9.2)
   M_90     1/(2 r90) = eta* >= 0.9                              (k2c := kd)
   1.6 M_90 same, ramp-end convention (persist.py's 34.7)
plus the aspartic acid rate 2.73e-14 at M_ex and M_90.

PREDICTIONS (written 2026-10-13 before any result-producing run; never edited).  Numbers marked "closed form" were worked out by
hand from the equations, and persist.log's r_c, r90, r99 (0.08698, 0.02307, 0.00127) were already known; no draw-level number
(pipeline) had been computed.  A smoke run at G_WIDE_N=2000 with output discarded may be made only to catch crashes.
 Analytic (units kd = k = 1)
 P1  Identity.  At the racemic state (feed a_rac, independent of r) the antisymmetric eigenvalue of the Jacobian equals
     a_rac - 1 - 2r to 1e-5 for r in {0, 0.02, 0.05, 0.087, 0.2}; so the racemic state loses stability at 2r = a_rac - 1 and the
     pure branch exists exactly for r below (a_rac - 1)/2: this equals persist's numeric r_c = 0.08698 to 2e-4.
 P2  Closed form (k0 neglected, f = k = kd = 1): 2 a_rac^2 - (2 - ki) a_rac - ki a0 = 0, so r_c = (a_rac - 1)/2 = 0.0895
     at ki = 1, a0 = 1.6: within 4% of the numeric 0.08698.  The hard necessary bound r_c < (a0 - 1)/2 = 0.3 holds (r_c is 3.4
     times smaller).
 P3  Conversions with k2c := kd: M_ex = 5.75 (+-0.05), M_90 = 21.7 (+-0.3), M_99 = 394 (+-10); ramp-end convention x 1.6:
     9.20, 34.7, 630.  All exceed 1, so the draft's test is below the exists threshold by 5.7x in this network.
 P4  Regime.  At a0 = 1.6, over ki in {0.3, 1, 3} x k0 in {1e-4, 1e-3, 1e-2}, M_ex stays within 2.5-16 (always above 1; the
     weakest cross-inhibition gives the largest).  At the repo ki = 1, k0 = 1e-3, M_ex falls below 1 only for a0 >= about
     5.5: a0 = 5 gives M_ex > 1 and a0 = 7 gives M_ex < 1 (the draft's test is about right only when the reaction runs far
     beyond its instability).
 P5  Direct deterministic ramp (a0: 0.8 -> 1.6 over tau = 2000 /k, then held 3000 /k; start = racemic state with a 1e-8
     relative excess; no noise, g = 0): final eta equals persist's eta*(r) within 0.01 for r = 0.5 r_c and 0.8 r_c; final
     |eta| < 1e-6 for r = 1.05 r_c and also for r = 0.25 and r = 0.5 (that is M = 2 and M = 1: the draft's test is passed, and
     nothing amplifies).
 Draws (committed pipeline)
 P6  Share of each world's actual seed-0 draws with k2 c > M * 2 K_RAC(Ala): the 14 icy worlds (same k2, c ranges): M = 1 80-82%;
     M_ex 72.5-75%; 1.6 M_ex 70-73%; M_90 66-68.5%; 1.6 M_90 64-66.5%.  Early Earth: 100% at M = 1, >= 99.9 at M_ex,
     >= 99.7 at 1.6 M_ex, 99.1 +- 0.4 at M_90, 98.5 +- 0.5 at 1.6 M_90.
 P7  Regression (hard assertion at N = 200,000): the unmasked rerun equals the committed multiwall main / x10 / sweep3 (rule none)
     P_pved, no lid and final, for all 15 worlds (g_wide.committed(): multiwall.log, Enceladus seedfix.log).
 P8  Monotone (hard assertion): for every world, column and lid state, P_pved and P_chance do not rise and P_none does not fall as M rises.
 P9  The draft's own test applied to every draw (M = 1) moves each world's multiwall main no-lid P_pved by <= 0.5 points.
 P10 Share of currently winning draws (multiwall main, no lid) that fail the stricter test, icy worlds: M_ex 3-9%; 1.6 M_90 8-17%.
     Early Earth: <= 0.3% at M_ex, <= 2% at 1.6 M_90.
 P11 Change in multiwall main no-lid P_pved (points, vs the unmasked rerun): Europa and Ganymede -1 to -4 at M_ex, -3 to -8 at
     1.6 M_90; Enceladus -0.8 to -3 at M_ex, -2 to -6 at 1.6 M_90; Early Earth > -0.2 at M_ex and > -2.0 at 1.6 M_90; Early Earth
     sweep3 (rule none, the abstract's "at least 67%") > -1.2 at 1.6 M_90.
 P12 Order unchanged at every M and both columns: final main Earth > Ganymede > Europa > Enceladus; Enceladus final = 0.0;
     the largest P_chance over the 15 worlds (main, no lid) stays <= 15.0.
 P13 Abstract numbers (Ganymede final x10-to-main, Europa final x10-to-main, Enceladus no-lid x10-to-main, Early Earth sweep3
     "at least 67%"): at M_ex and 1.6 M_ex no end of any of them moves by more than 3 points; at 1.6 M_90 at least one end moves by
     3 points or more.
 P14 Aspartic acid (2.73e-14) at M_90: Early Earth multiwall main no-lid P_pved in 96.5-98.7.
 P15 Cost: the whole script finishes in under 20 minutes at 2 threads under nice (expected about 1-2 minutes).

Run: nice -n 10 env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 /opt/anaconda3/bin/python kcheck.py > kcheck.log 2>&1
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

import g_wide as G          # installs the seedfix rng wrapper, imports uncertainty / sweep2 / multiwall
import persist as P
import sweep3 as S3

U, S, M = G.U, G.S, G.M
T0 = time.time()
FULL = U.N == 200_000
KALA, KASP = S.K_RAC["Ala"], S.K_RAC["Asp"]
RES = {}


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


# ------------------------------------------------------------------------------------------------ analytic
def rc_of(**kw):
    """Largest r with a pure branch (persist's eta2 formula), constants as persist.par(**kw)."""
    grid = np.linspace(0.0, 5.0, 5001)
    return P.first_crossing(lambda r: P.eta2_of(P.par(r=r, **kw)), grid, 0.0)


def jac(p, y):
    y = np.array(y, float)
    J = np.empty((3, 3))
    for j in range(3):
        h = 1e-6 * max(abs(y[j]), 1e-3)
        e = np.zeros(3)
        e[j] = h
        J[:, j] = (np.array(P.rhs(0, y + e, p)) - np.array(P.rhs(0, y - e, p))) / (2 * h)
    return J


def closed_arac(ki, a0):
    """k0 neglected, f = k = kd = 1: 2 a^2 - (2 - ki) a - ki a0 = 0."""
    return ((2 - ki) + np.sqrt((2 - ki) ** 2 + 8 * ki * a0)) / 4


def ramp_final(r, tau=2000.0, hold=3000.0, u0=1e-8):
    p = P.par(r=r, a0=0.8)
    y = P.racemic(p)
    y = [y[0], y[1] * (1 + u0), y[2] * (1 - u0)]

    def f(t, yy):
        q = dict(p)
        q["a0"] = 0.8 + 0.8 * min(t / tau, 1.0)
        return P.rhs(t, yy, q)
    sol = solve_ivp(f, [0, tau + hold], y, method="LSODA", rtol=1e-10, atol=1e-20)
    ye = sol.y[:, -1]
    return (ye[1] - ye[2]) / (ye[1] + ye[2])


def analytic():
    print("=== A. Analytic: where does the pure/mixed transition sit, in kd units (k = kd = ki = f = 1, k0 = 1e-3, a0 = 1.6)? ===")
    p0 = P.par()
    arac = P.racemic(p0)[0]
    rc, r90, r99 = rc_of(), None, None
    grid = np.linspace(0.0, 0.2, 20001)
    r90 = P.first_crossing(lambda r: P.eta2_of(P.par(r=r)), grid, 0.81)
    r99 = P.first_crossing(lambda r: P.eta2_of(P.par(r=r)), grid, 0.9801)
    print(f"a_rac (racemic feed, any r) = {arac:.6f};  (a_rac - 1)/2 = {(arac - 1) / 2:.6f};  numeric r_c = {rc:.6f}, "
          f"r(eta*>=0.9) = {r90:.6f}, r(eta*>=0.99) = {r99:.6f}")
    dev = 0.0
    for r in (0.0, 0.02, 0.05, 0.087, 0.2):
        p = P.par(r=r)
        y = P.racemic(p)
        lam = (lambda J: J[1, 1] - J[1, 2])(jac(p, y))
        dev = max(dev, abs(lam - (y[0] - 1 - 2 * r)))
        print(f"  r = {r:<6} antisymmetric eigenvalue {lam:+.6f}   a_rac - 1 - 2r = {y[0] - 1 - 2 * r:+.6f}")
    check("P1", dev < 1e-5 and abs((arac - 1) / 2 - rc) < 2e-4,
          f"max |eigenvalue - (a_rac - 1 - 2r)| = {dev:.2e}; (a_rac - 1)/2 = {(arac - 1) / 2:.6f} vs r_c {rc:.6f}")
    cf = (closed_arac(1.0, 1.6) - 1) / 2
    check("P2", abs(cf / rc - 1) < 0.04 and abs(0.3 / rc - 3.4) < 0.1,
          f"closed form r_c = {cf:.5f} ({100 * (cf / rc - 1):+.1f}% vs numeric); necessary bound 0.3 is {0.3 / rc:.2f} x r_c")
    m_ex, m_90, m_99 = 1 / (2 * rc), 1 / (2 * r90), 1 / (2 * r99)
    print(f"M = 1/(2 r) with k2c := kd:  M_ex = {m_ex:.3f}  M_90 = {m_90:.3f}  M_99 = {m_99:.1f};  x1.6 (ramp end): "
          f"{1.6 * m_ex:.3f}  {1.6 * m_90:.3f}  {1.6 * m_99:.1f}")
    check("P3", abs(m_ex - 5.75) <= 0.05 and abs(m_90 - 21.7) <= 0.3 and abs(m_99 - 394) <= 10 and m_ex > 1,
          f"M_ex {m_ex:.2f}, M_90 {m_90:.2f}, M_99 {m_99:.0f}; x1.6: {1.6 * m_ex:.2f}, {1.6 * m_90:.1f}, {1.6 * m_99:.0f}")

    print("\nM_ex = 1/(2 r_c) over the constants the draft does not sample (a0 = 1.6): rows ki, columns k0")
    grid_m, allm = {}, []
    for ki in (0.3, 1.0, 3.0):
        row = []
        for k0 in (1e-4, 1e-3, 1e-2):
            m = 1 / (2 * rc_of(ki=ki, k0=k0))
            row.append(m)
            allm.append(m)
        grid_m[ki] = row
        print(f"  ki = {ki:<4}  " + "  ".join(f"k0 {k0:g}: {m:6.2f}" for k0, m in zip((1e-4, 1e-3, 1e-2), row))
              + f"   (closed form, k0 = 0: {1 / (closed_arac(ki, 1.6) - 1):.2f})")
    print("M_ex against how far the ramp goes past the crossing (repo ki = 1, k0 = 1e-3):")
    by_a0 = {}
    for a0 in (1.2, 1.6, 2.0, 3.0, 5.0, 7.0, 10.0):
        by_a0[a0] = 1 / (2 * rc_of(a0=a0))
        print(f"  a0 = {a0:<4}  M_ex = {by_a0[a0]:7.3f}   (closed form {1 / (closed_arac(1.0, a0) - 1):7.3f})")
    check("P4", min(allm) >= 2.5 and max(allm) <= 16 and by_a0[5.0] > 1 and by_a0[7.0] < 1,
          f"grid {min(allm):.2f}-{max(allm):.2f}; a0 = 5: {by_a0[5.0]:.2f}, a0 = 7: {by_a0[7.0]:.2f}")

    print("\nDirect deterministic ramp a0: 0.8 -> 1.6 over tau = 2000, hold 3000, 1e-8 relative start excess, g = 0, no noise")
    ok, det = True, []
    for lab, r, kind in (("0.5 r_c", 0.5 * rc, "eta"), ("0.8 r_c", 0.8 * rc, "eta"), ("1.05 r_c", 1.05 * rc, "dead"),
                         ("M = 2 (r = 0.25)", 0.25, "dead"), ("M = 1 (r = 0.5)", 0.5, "dead")):
        e = ramp_final(r)
        tgt = float(np.sqrt(P.asym(P.par(r=r))[2])) if kind == "eta" else 0.0
        good = abs(abs(e) - tgt) < 0.01 if kind == "eta" else abs(e) < 1e-6
        ok &= good
        det.append(f"{lab}: {abs(e):.3g}")
        print(f"  r = {r:.5f} ({lab}): final |eta| = {abs(e):.6g}" + (f"   eta*(r) = {tgt:.6f}" if kind == "eta" else "   (no amplification expected)"))
    check("P5", ok, "; ".join(det))
    return dict(rc=rc, r90=r90, m_ex=m_ex, m_90=m_90)


# ------------------------------------------------------------------------------------------------ draws
def draws(seed=0):
    _, W = S.chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        k = x["rate k2"] * x["concentration"]
        runs = {}
        for li, lid in (("nolid", none), ("final", salt)):
            runs[li, "main"] = M.multi(x, a[1], a[4], ex, lid)[:2]
            runs[li, "x10"] = M.multi(x, a[1], a[4], ex, lid, scale=10.0)[:2]
            runs[li, "none"] = M.multi(x, a[1], a[4], ex, lid, rule="none")[:2]
        out[name] = (k, runs)
    return out


def sp(w, r, kk=None):
    return S3.split(w, r) if kk is None else S3.split(w & kk, r & kk)


def pved(d, name, col, li, K=None, m=None):
    k, runs = d[name]
    w, r = runs[li, col]
    kk = None if m is None else k > m * 2 * K
    return sp(w, r, kk)


def main():
    print(f"kcheck.py, N = {U.N}, seed 0; Enceladus on stream [s, 0, 1] (seedfix wrapper); G_NEW = {S.G_NEW}")
    A = analytic()
    levels = [("draft 1", 1.0, KALA), ("nec 1.667", 1 / 0.6, KALA), ("exists", A["m_ex"], KALA),
              ("exists x1.6", 1.6 * A["m_ex"], KALA), ("pure0.9", A["m_90"], KALA), ("pure0.9 x1.6", 1.6 * A["m_90"], KALA),
              ("Asp exists", A["m_ex"], KASP), ("Asp pure0.9", A["m_90"], KASP)]
    print("\n=== B. Draws: committed pipeline (multiwall.multi) with an extra mask k2*c > M * 2 K_RAC on every draw ===")
    d = draws()
    names = list(d)
    G.lap("draws done")

    print("\n--- B0. Regression of the unmasked rerun against the committed logs ---")
    old, where = G.committed()
    cols = {"none": 0, "main": 4, "x10": 6}
    bad = []
    for n in names:
        for col, off in cols.items():
            got = tuple(f"{pved(d, n, col, li)[0]:.1f}" for li in ("nolid", "final"))
            if FULL and got != tuple(old[n][off:off + 2]):
                bad.append((n, col, where[n], got, old[n][off:off + 2]))
    check("P7", not bad if FULL else True, "15 worlds x (none, main, x10) x (no lid, final) equal the committed values"
          if not bad else f"mismatches {bad}")

    print("\n--- B1. Draw shares with k2 c above M * 2 K_RAC (unconditional, all 200,000 draws), % ---")
    print(f"{'world':18}" + "".join(f"{lab:>14}" for lab, _, _ in levels))
    uncond = {}
    for n in names:
        k = d[n][0]
        uncond[n] = [100 * np.mean(k > m * 2 * K) for _, m, K in levels]
        print(f"{n:18}" + "".join(f"{v:14.2f}" for v in uncond[n]))
    moons = [n for n in names if n != "Early Earth ocean"]
    ok = True
    ranges = [(80, 82), (None, None), (72.5, 75), (70, 73), (66, 68.5), (64, 66.5)]
    for j, (lo, hi) in enumerate(ranges):
        if lo is not None:
            ok &= all(lo <= uncond[n][j] <= hi for n in moons)
    ee = uncond["Early Earth ocean"]
    ok &= ee[0] == 100.0 and ee[2] >= 99.9 and ee[3] >= 99.7 and 98.7 <= ee[4] <= 99.5 and 98.0 <= ee[5] <= 99.0
    check("P6", ok, f"icy worlds M=1 {min(uncond[n][0] for n in moons):.1f}-{max(uncond[n][0] for n in moons):.1f}, exists "
          f"{min(uncond[n][2] for n in moons):.1f}-{max(uncond[n][2] for n in moons):.1f}, x1.6 {min(uncond[n][3] for n in moons):.1f}-"
          f"{max(uncond[n][3] for n in moons):.1f}, pure0.9 {min(uncond[n][4] for n in moons):.1f}-{max(uncond[n][4] for n in moons):.1f}, "
          f"x1.6 {min(uncond[n][5] for n in moons):.1f}-{max(uncond[n][5] for n in moons):.1f}; Earth {', '.join(f'{v:.2f}' for v in ee[:6])}")

    print("\n--- B2. P_pved %, multiwall main, NO LID (rows: world; columns: M levels; first column = committed pipeline, no extra mask) ---")
    base = {n: pved(d, n, "main", "nolid") for n in names}
    hdr = f"{'world':18}{'unmasked':>10}" + "".join(f"{lab:>14}" for lab, _, _ in levels)
    print(hdr)
    tab = {}
    for n in sorted(names, key=lambda n: -base[n][0]):
        row = [pved(d, n, "main", "nolid", K, m) for _, m, K in levels]
        tab[n] = row
        print(f"{n:18}{base[n][0]:10.1f}" + "".join(f"{r[0]:14.1f}" for r in row))

    print("\n--- B3. P_pved %, multiwall main, FINAL (with salt lid); x10 final ---")
    print(hdr)
    fin, fin10 = {}, {}
    for n in sorted(names, key=lambda n: -base[n][0]):
        b, b10 = pved(d, n, "main", "final")[0], pved(d, n, "x10", "final")[0]
        fin[n] = [pved(d, n, "main", "final", K, m)[0] for _, m, K in levels]
        fin10[n] = [pved(d, n, "x10", "final", K, m)[0] for _, m, K in levels]
        print(f"{n:18}{b:10.1f}" + "".join(f"{v:14.1f}" for v in fin[n]) + "   (main)")
        print(f"{'  x10':18}{b10:10.1f}" + "".join(f"{v:14.1f}" for v in fin10[n]))

    print("\n--- B4. Three-way split, multiwall main, no lid: P_chance / P_none (unmasked, then each level) ---")
    chance_max = [0.0] * (len(levels) + 1)
    for n in names:
        b = base[n]
        row = [pved(d, n, "main", "nolid", K, m) for _, m, K in levels]
        chance_max[0] = max(chance_max[0], b[1])
        for j, r in enumerate(row):
            chance_max[j + 1] = max(chance_max[j + 1], r[1])
    print("largest P_chance over the 15 worlds: unmasked " + f"{chance_max[0]:.1f}; " + "; ".join(
        f"{lab} {v:.1f}" for (lab, _, _), v in zip(levels, chance_max[1:])))
    for n in ("Enceladus", "Europa", "Ganymede", "Early Earth ocean"):
        print(f"{n:18} chance / none  unmasked {base[n][1]:4.1f} / {base[n][2]:5.1f} | " + " | ".join(
            f"{r[1]:4.1f} / {r[2]:5.1f}" for r in [pved(d, n, 'main', 'nolid', K, m) for _, m, K in levels]))

    print("\n--- B5. Share of currently winning draws (multiwall main, no lid) that fail the stricter test, % ---")
    print(f"{'world':18}" + "".join(f"{lab:>14}" for lab, _, _ in levels))
    failw = {}
    for n in names:
        k, runs = d[n]
        w = runs["nolid", "main"][0]
        failw[n] = [100 * (1 - (w & (k > m * 2 * K)).sum() / max(w.sum(), 1)) for _, m, K in levels]
        print(f"{n:18}" + "".join(f"{v:14.2f}" for v in failw[n]))

    print("\n--- B6. Early Earth sweep3 rule 'none' (the abstract's 'at least 67%'), no lid, P_pved % ---")
    eb = pved(d, "Early Earth ocean", "none", "nolid")[0]
    e6 = [pved(d, "Early Earth ocean", "none", "nolid", K, m)[0] for _, m, K in levels]
    print(f"  unmasked {eb:.1f} | " + " | ".join(f"{lab} {v:.1f}" for (lab, _, _), v in zip(levels, e6)))

    # ------------------------------------------------------------------ predictions on draws
    print("\n=== PREDICTIONS (draw level) ===")
    viol = []
    for n in names:
        for col in ("main", "x10", "none"):
            for li in ("nolid", "final"):
                prev = pved(d, n, col, li, KALA, 1.0)
                for _, m, K in levels[1:6]:
                    cur = pved(d, n, col, li, K, m)
                    if cur[0] > prev[0] + 1e-9 or cur[1] > prev[1] + 1e-9 or cur[2] < prev[2] - 1e-9:
                        viol.append((n, col, li, m))
                    prev = cur
    assert not viol, f"P8 FAILED (monotone): {viol}"
    check("P8", True, f"{len(names) * 3 * 2} (world, column, lid) cells over M = 1 ... {levels[5][1]:.1f}: no violation")
    d1 = {n: base[n][0] - tab[n][0][0] for n in names}
    check("P9", max(abs(v) for v in d1.values()) <= 0.5, f"M = 1 on every draw moves P_pved by at most "
          f"{max(abs(v) for v in d1.values()):.2f} points ({max(d1, key=lambda n: abs(d1[n]))}); signed: "
          + ", ".join(f"{n.split()[0]} {-d1[n]:+.2f}" for n in ("Enceladus", "Europa", "Ganymede", "Early Earth ocean")))
    j_ex, j_x = 2, 5
    icy = [n for n in names if n != "Early Earth ocean"]
    ok = all(3 <= failw[n][j_ex] <= 9 and 8 <= failw[n][j_x] <= 17 for n in icy)
    ok &= failw["Early Earth ocean"][j_ex] <= 0.3 and failw["Early Earth ocean"][j_x] <= 2.0
    check("P10", ok, f"icy worlds exists {min(failw[n][j_ex] for n in icy):.1f}-{max(failw[n][j_ex] for n in icy):.1f}%, x1.6 pure0.9 "
          f"{min(failw[n][j_x] for n in icy):.1f}-{max(failw[n][j_x] for n in icy):.1f}%; Earth {failw['Early Earth ocean'][j_ex]:.2f}% / "
          f"{failw['Early Earth ocean'][j_x]:.2f}%")
    chg = lambda n, j: tab[n][j][0] - base[n][0]
    ok = (-4 <= chg("Europa", j_ex) <= -1 and -4 <= chg("Ganymede", j_ex) <= -1 and -8 <= chg("Europa", j_x) <= -3
          and -8 <= chg("Ganymede", j_x) <= -3 and -3 <= chg("Enceladus", j_ex) <= -0.8 and -6 <= chg("Enceladus", j_x) <= -2
          and chg("Early Earth ocean", j_ex) > -0.2 and chg("Early Earth ocean", j_x) > -2.0 and e6[j_x] - eb > -1.2)
    check("P11", ok, f"M_ex: Europa {chg('Europa', j_ex):+.2f}, Ganymede {chg('Ganymede', j_ex):+.2f}, Enceladus {chg('Enceladus', j_ex):+.2f}, "
          f"Earth {chg('Early Earth ocean', j_ex):+.2f}; x1.6 M_90: Europa {chg('Europa', j_x):+.2f}, Ganymede {chg('Ganymede', j_x):+.2f}, "
          f"Enceladus {chg('Enceladus', j_x):+.2f}, Earth {chg('Early Earth ocean', j_x):+.2f}, Earth sweep3-none {e6[j_x] - eb:+.2f}")
    ordok = all(fin[names_[0]][j] > fin[names_[1]][j] > fin[names_[2]][j] >= fin[names_[3]][j]
                for names_ in [("Early Earth ocean", "Ganymede", "Europa", "Enceladus")] for j in range(len(levels)))
    ordok &= all(fin["Enceladus"][j] == 0.0 for j in range(len(levels))) and max(chance_max) <= 15.0
    check("P12", ordok, f"final main order Earth > Ganymede > Europa > Enceladus at every level; Enceladus final "
          f"{max(fin['Enceladus']):.1f}; max P_chance {max(chance_max):.1f}")
    # P13: abstract numbers
    def absn(j):
        """j = None: unmasked.  Returns the abstract's numbers (rounded to 1 decimal)."""
        if j is None:
            gf = (pved(d, "Ganymede", "x10", "final")[0], pved(d, "Ganymede", "main", "final")[0])
            ef = (pved(d, "Europa", "x10", "final")[0], pved(d, "Europa", "main", "final")[0])
            en = (pved(d, "Enceladus", "x10", "nolid")[0], pved(d, "Enceladus", "main", "nolid")[0])
            return gf, ef, en, eb
        _, m, K = levels[j]
        gf = (pved(d, "Ganymede", "x10", "final", K, m)[0], pved(d, "Ganymede", "main", "final", K, m)[0])
        ef = (pved(d, "Europa", "x10", "final", K, m)[0], pved(d, "Europa", "main", "final", K, m)[0])
        en = (pved(d, "Enceladus", "x10", "nolid", K, m)[0], pved(d, "Enceladus", "main", "nolid", K, m)[0])
        return gf, ef, en, pved(d, "Early Earth ocean", "none", "nolid", K, m)[0]
    print("\n--- B7. The abstract's numbers (x10-to-main final for Ganymede and Europa; Enceladus no-lid x10-to-main; Earth sweep3 'at least') ---")
    b0 = absn(None)
    print(f"  unmasked        Ganymede {b0[0][0]:.1f}-{b0[0][1]:.1f}  Europa {b0[1][0]:.1f}-{b0[1][1]:.1f}  Enceladus no lid {b0[2][0]:.1f}-{b0[2][1]:.1f}  Earth {b0[3]:.1f}")
    shifts = {}
    for j, (lab, _, _) in enumerate(levels):
        v = absn(j)
        sh = max(abs(a - b) for x, y in zip(v[:3], b0[:3]) for a, b in zip(x, y))
        sh = max(sh, abs(v[3] - b0[3]))
        shifts[lab] = sh
        print(f"  {lab:14}  Ganymede {v[0][0]:.1f}-{v[0][1]:.1f}  Europa {v[1][0]:.1f}-{v[1][1]:.1f}  Enceladus no lid {v[2][0]:.1f}-{v[2][1]:.1f}  "
              f"Earth {v[3]:.1f}   largest shift {sh:.2f}")
    check("P13", shifts["exists"] <= 3 and shifts["exists x1.6"] <= 3 and shifts["pure0.9 x1.6"] >= 3,
          "largest shift of any abstract number: " + ", ".join(f"{k} {v:.2f}" for k, v in shifts.items()))
    e_asp = tab["Early Earth ocean"][7][0]
    check("P14", 96.5 <= e_asp <= 98.7, f"Early Earth multiwall main no lid at Asp M_90: {e_asp:.1f}")
    el = time.time() - T0
    check("P15", el < 1200, f"{el:.0f} s")
    print(f"\n{sum(RES.values())} of {len(RES)} held" + ("" if all(RES.values()) else "; failed: " + ", ".join(k for k, v in RES.items() if not v)))


if __name__ == "__main__":
    main()
