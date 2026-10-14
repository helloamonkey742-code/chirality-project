"""Does a pure (one-handed) monomer pool stay pure?  Tests Avetisov (1999), Ch. 6 of Advances in BioChirality, pp. 69-83.

His claim (pp. 78-82, Eqs. 1, 14-16, 8-9): a pool made pure by a weak chiral field cannot stay pure once
enantioselective biochemistry arises, because the field's advantage factor g is small against the racemization factor
K_r (Eq. 1: -K_r*eta + g*(1-eta^2) = 0) and the new "enantioselective pressure" K_p*gamma opposes it (Eq. 14);
purity made by spontaneous symmetry breaking (autocatalysis) needs catalysts about as selective as enzymes (Eq. 9).
Wording used here: eta = (x_L - x_D)/(x_L + x_D) is the "chiral polarization"; K_r = 1/tau_r the "racemization factor";
gamma = (S-1)/(S+1) the enantioselectivity of the consumer and S = right:wrong rate ratio; K_p its intensity.

MODEL (the draft's Eq. 2 plus racemization, a consumer, and an optional catalyst-selectivity V).  A = achiral feed.
  da/dt = f(a0-a) - 2 k0 a - a (kL L + kD D),          kL = k(1+g), kD = k(1-g)
  dL/dt = k0 a + a[kL (1+V)/2 L + kD (1-V)/2 D] - ki L D - kd L - r (L-D) - c h wL L
  dD/dt = k0 a + a[kD (1+V)/2 D + kL (1-V)/2 L] - ki L D - kd D + r (L-D) - c h wD D
V = 1 is the draft's model (perfect autocatalyst, frank.py); V < 1 is Avetisov's scheme (7), V = (k1-k2)/(k1+k2).
Racemization: L<->D at rate r each way, so eta decays as exp(-2 r t); r = k_rac of the draft (Bada and Schroeder
1975 via Steen et al. 2013; sweep2.py K_RAC; 2 k_rac is the draft's own "k > 2 k_rac" test; Arrhenius as
rerun_base.py k_at).  Constants: k = ki = kd = f = 1, k0 = 1e-3, a0 ramp end 1.6 (frank.py, Eq. 2 text), g = 0 unless
stated.  Consumer: removes L at c h wL, D at c h wD with wL = 1 + gamma, wD = 1 - gamma (so wL/wD = S, mean 1); it
exerts exactly Avetisov's pressure d(eta)/dt = -c h gamma (1-eta^2) (his Eq. 14, K_p = c).  "Right" = the pool's
majority hand at switch-on.  "Works only when pure": gate h = ((1+eta_right)/2)^N, the chance that an N-link chain
is defect-free (Avetisov p. 75: N ~ 40, ~130); N = 0 means always on (his literal K_p*gamma).  Consumer is switched on
suddenly at the pure steady state of the no-consumer system (the hard case).  Deterministic ODEs (no counting noise).
Time unit 1/(k*1).  Q1 uses real seconds/years.
Definitions: "persists" = same-hand eta_right >= 0.5 at t = 5000; <= -0.5 is "flipped" (excess now the other hand).
Analytic reduction used (no consumer, g = 0): du/dt = (k a V - kd - 2r) u, so the pure branch sits at k a V = kd + 2r
(a = c0/(kV), c0 = kd+2r); s = [f(a0-a) - 2 k0 a]/(k a); 1 - eta^2 = [4 k0 a + 2 (k a - kd) s]/(ki s^2).
Pitchfork (racemic state loses stability = Avetisov's p_c = 1) where 1 - eta^2 = 1.

PREDICTIONS (written 2026-10-13 before the run; only toy checks of this model had been run, one-hand-at-a-time
and no result was written down).  Year = 3.15576e7 s.
Q1 (closed pool, autocatalysis off, only racemization, start eta = 1; excess eta(t) = exp(-2 k_rac t)):
 P1  273 K: alanine t_half = 1.20 Myr, t_1% = 7.99 Myr; aspartic acid 0.402 and 2.67 Myr (each within 1%).
 P2  Arrhenius (Ea 93/110 kJ/mol): alanine 308 K t_half 11.4 kyr, t_1% 75.9 kyr; 343 K 281 yr and 1.87 kyr;
     aspartic acid 343 K 20.4 yr and 135 yr (each within 3%).
 P3  vs the paper's tau (Table 1: Enceladus 1e6-1e9, Europa 1e8-4e9, Early Earth 1e7-5e8 yr) at 273 K, alanine:
     t_1% is below tau for the whole Earth and Europa ranges (tau/t_1% = 1.25-63 and 12.5-500) and for 30% of
     Enceladus's log-uniform range (14% for aspartic acid); at 308 K and 343 K every tau of every world exceeds
     t_1% by >= 13x.  So without a restoring process no pure pool survives the ocean age.
 P4  Avetisov Eq. 1 with the field only (g_rate = g k2 a, K_r = 2 k_rac): the steady excess is <= 5e-6 for every
     input in Table 1 (largest at g = 2.64e-16, k2 a = 3e-4 /s: about 4e-6); AF/RF << 1 everywhere.
Q2 (autocatalysis on, open and fed; repo constants; g = 0):
 P5  The pure branch exists, is stable, and is the only stable state for r < r_c = 0.087 kd (+-0.003); above it only
     the racemic state is stable.  The onset is exactly k a V = kd + 2r (pitchfork, p_c = 1).
 P6  Purity eta* (analytic formula): 0.994 at r = 0, 0.957 at r = 0.01; eta* >= 0.9 only for r <= 0.023 and
     eta* >= 0.99 only for r <= 0.0013 (each within 5%).  The full ODE agrees with the formula to 1e-6.
 P7  Every initial excess eta0 > 0, down to 1e-6, ends on the pure branch of its own sign for r < r_c (the racemic
     state is a saddle, so no finite basin boundary); exactly racemic stays racemic.
 P8  Minimum strengths at r = 0.01, other constants fixed: copying k_min = 0.79 (k a0/(kd+2r) = 1.24) for
     eta* >= 0.9, and between 0.64 and 0.70 for the pure state to exist at all; inhibition ki_min = 0.083 (exists)
     and 0.44 (eta* >= 0.9); at r = 0: 0.011 and 0.059.  With ki = 0 no pure state exists for any r > 0.
     In the draft's rate terms the pure state needs k a0 >= 9.2 * (2 r) to exist and >= 34.7 * (2 r) for eta* >= 0.9
     (at r_c and r(0.9) of P5/P6), i.e. 9x and 35x the draft's "k > 2 k_rac" test, at flow = copying rate.
     Share of Table-1 draws (log-uniform k2, a; 200,000, seed 0) with k2 a above 1x/9.2x/34.7x of 2 k_rac:
     moons 81% -> 71% -> 65% (each within 2 points), Early Earth 100% -> 100% (>= 99.7) -> about 98.5%.
 P9  Catalyst selectivity (V, Avetisov Eq. 9 direction): V_min = 0.85 (exists), 0.956 (eta* >= 0.9), 0.9975
     (eta* >= 0.99) at r = 0; at r = 0.01 0.87, 0.975 and none (eta* at V = 1 is only 0.957).  So a pure state of
     practical purity needs V of 0.96-0.998: enzyme-like in Avetisov's sense, not "about 0.5".
Q3 (consumer with selectivity S, switched on at the pure state, r = 0.01):
 P10 Autocatalysis off (closed pool, k = ki = kd = f = k0 = 0, start eta = 1, c = 0.1): the excess is lost for every
     S >= 1, so no critical S exists.  Time to 1%, consumer/no consumer: exactly 1.000 at S = 1 (a non-selective
     drain does not change eta), non-increasing in S; gated N = 40: >= 0.97 even at S = infinity; always on:
     <= 0.2 at S = infinity.  (A consumer could keep eta only if S < 1, i.e. anti-selective.)
 P11 Autocatalysis on, gate N = 40 (and 130, 10): the excess persists for every S up to infinity and every
     c from 0.01 to 100 (the gate throttles the consumer; purity falls to about 0.67-0.96, no critical S).
 P12 Autocatalysis on, always-on consumer (N = 0): a critical S exists and falls as c grows.  c = 0.1:
     S_c between 3 and 9; c = 0.3: between 1 and 1.5; c >= 1: lost even at S = 1 (a non-selective sink is an extra
     flow: it needs k a0 > kd + 2r + c).  Above S_c the pool flips to the other hand (eta_right <= -0.5), not to racemic.
 P13 Enzyme-like selectivities: not sourced in journals/ija_draft.md, so no comparison is made (the chapter's own
     N-link purity criterion gives S about N-1 = 39 and 129, printed for reference only).
Sanity P14: a pure-D start in a D-preferring world gives the exact mirror (eta -> -eta, a unchanged), asserted to
 1e-9 for g = 0 and for g = 1e-17 with the bias flipped; the analytic pure branch equals fsolve of the full ODE
 fixed point to 1e-8 (asserted).

RESULTS (persist.log, 2026-10-13): predictions above are unchanged.  13 of 14 held; P11 held only in part.
 P1 held (log 2, 5).  P2 held (3-4, 7).  P3 held (9-14: 1.25-62.6, 12.5-501, 30.1%, 14.2%, min ratio 13.2).
 P4 held (16-20; largest 4.3e-6).  P5 held (24-38; r_c 0.08698; pitchfork identity to 1e-6).
 P6 held (24-37; 0.994391, 0.957444, r 0.02307, 0.00127; ODE equals formula to 6 digits).
 P7 held (42-49).  P8 held (51-62; k_min 0.7904 and 0.6822, ki_min 0.0833/0.4384 and 0.01119/0.05887, 9.20x and 34.68x,
 shares 81.02/71.48/64.87 and 100/99.89/98.48).  P9 held (53, 57).  P10 held (69-72; 1.0000 at S = 1, 0.9947, 0.0665).
 P11 FAILED for N = 10: at c = 3 the excess falls to 0.30-0.63 (below the 0.5 line at S >= 3, log 93) and at c >= 10 it
 flips (94-95).  It held for N = 40 (97-104; 0.671-0.957) and N = 130 (106-113; 0.89-0.957).
 P12 held (81-86; S_c = 5.86 at c = 0.1, 1.149 at c = 0.3, lost at S = 1 for c >= 1, flips to -0.97).
 P13 held (115-116).  P14 held (64-65, 73, 114).
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

YR = 3.15576e7
K273 = {"Ala": 9.13e-15, "Asp": 2.73e-14}   # 1/s, sweep2.py K_RAC
EA = {"Ala": 93e3, "Asp": 110e3}            # J/mol, rerun_base.py EA
TAU = {"Enceladus": (1e6, 1e9), "Europa": (1e8, 4e9), "Early Earth": (1e7, 5e8)}  # yr, Table 1
BASE = dict(k=1.0, ki=1.0, kd=1.0, f=1.0, k0=1e-3, a0=1.6, g=0.0, V=1.0, r=0.0, c=0.0, gam=0.0, N=0, sg=1)
T_END = 5000.0


def say(*a):
    print(*a, flush=True)


def par(**kw):
    p = dict(BASE)
    p.update(kw)
    return p


def eta_right(y, sg):
    s = y[1] + y[2]
    return sg * (y[1] - y[2]) / max(s, 1e-300)


def rhs(t, y, p):
    a, L, D = y
    L, D = max(L, 0.0), max(D, 0.0)
    kL, kD, V = p["k"] * (1 + p["g"]), p["k"] * (1 - p["g"]), p["V"]
    h = 1.0 if p["N"] == 0 else ((1 + np.clip(eta_right((a, L, D), p["sg"]), -1, 1)) / 2) ** p["N"]
    ch = p["c"] * h
    wL, wD = 1 + p["sg"] * p["gam"], 1 - p["sg"] * p["gam"]
    dL = p["k0"] * a + a * (kL * (1 + V) / 2 * L + kD * (1 - V) / 2 * D) - p["ki"] * L * D - p["kd"] * L \
        - p["r"] * (L - D) - ch * wL * L
    dD = p["k0"] * a + a * (kD * (1 + V) / 2 * D + kL * (1 - V) / 2 * L) - p["ki"] * L * D - p["kd"] * D \
        + p["r"] * (L - D) - ch * wD * D
    da = p["f"] * (p["a0"] - a) - 2 * p["k0"] * a - a * (kL * L + kD * D)
    return [da, dL, dD]


def integrate(p, y0, T=T_END, events=None):
    return solve_ivp(rhs, [0, T], y0, args=(p,), method="LSODA", rtol=1e-10, atol=1e-18, events=events)


def asym(p):
    """Analytic pure branch (no consumer, g = 0). Returns (a, s, eta2) or None."""
    c0 = p["kd"] + 2 * p["r"]
    a = c0 / (p["k"] * p["V"])
    s = (p["f"] * (p["a0"] - a) - 2 * p["k0"] * a) / (p["k"] * a)
    if s <= 0:
        return None
    eta2 = 1 - (4 * p["k0"] * a + 2 * (p["k"] * a - p["kd"]) * s) / (p["ki"] * s * s) if p["ki"] > 0 else -1.0
    return (a, s, eta2) if eta2 > 0 else None


def eta2_of(p):
    q = asym(p)
    return q[2] if q else -1.0


def pure_state(p):
    a, s, e2 = asym(p)
    u = np.sqrt(e2) * s
    return [a, (s + u) / 2, (s - u) / 2]


def racemic(p):
    def g2(a):
        s = (p["f"] * (p["a0"] - a) - 2 * p["k0"] * a) / (p["k"] * a)
        return 2 * p["k0"] * a + (p["k"] * a - p["kd"]) * s - p["ki"] * s * s / 2
    amax = p["f"] * p["a0"] / (p["f"] + 2 * p["k0"])
    a = brentq(g2, 1e-9, amax * (1 - 1e-12))
    s = (p["f"] * (p["a0"] - a) - 2 * p["k0"] * a) / (p["k"] * a)
    return [a, s / 2, s / 2]


def jac_eigs(p, y):
    y = np.array(y, float)
    J = np.empty((3, 3))
    for j in range(3):
        h = 1e-6 * max(abs(y[j]), 1e-3)
        e = np.zeros(3)
        e[j] = h
        J[:, j] = (np.array(rhs(0, y + e, p)) - np.array(rhs(0, y - e, p))) / (2 * h)
    return np.linalg.eigvals(J)


def first_crossing(fn, grid, thr):
    """Smallest-x root of fn(x) = thr along grid, or None."""
    v = [fn(x) - thr for x in grid]
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            return brentq(lambda x: fn(x) - thr, grid[i], grid[i + 1], xtol=1e-12)
    return None


def kr_at(aa, T):
    return K273[aa] * np.exp(EA[aa] / 8.314 * (1 / 273 - 1 / T))


# ------------------------------------------------------------------------------------------------------------- Q1
def q1():
    say("=== Q1. Autocatalysis off, only racemization (closed pool, start eta = 1; eta(t) = exp(-2 k_rac t)) ===")
    ts = {}
    for aa in ("Ala", "Asp"):
        for T in (273, 308, 343):
            k = kr_at(aa, T)
            th, t1 = np.log(2) / (2 * k) / YR, np.log(100) / (2 * k) / YR
            ts[aa, T] = (th, t1)
            # ODE cross-check of the closed form (L<->D at rate k, mirror start D too)
            p = par(r=k, k=0, ki=0, kd=0, f=0, k0=0)
            for y0 in ([0, 1.0, 0.0], [0, 0.0, 1.0]):
                tt = np.log(100) / (2 * k)
                sol = integrate(p, y0, T=tt)
                e = (sol.y[1, -1] - sol.y[2, -1]) / (sol.y[1, -1] + sol.y[2, -1]) * (1 if y0[1] else -1)
                assert abs(e - 0.01) < 1e-6, (aa, T, e)
            say(f"  {aa} {T} K: k_rac {k:.3e}/s  t_half {th:.4g} yr  t_1% {t1:.4g} yr")
    say("  tau (Table 1) vs t_1%: ratio tau/t_1% at the range ends, and log-uniform share of tau below t_1%")
    for aa in ("Ala", "Asp"):
        for T in (273, 308, 343):
            t1 = ts[aa, T][1]
            row = []
            for w, (lo, hi) in TAU.items():
                share = float(np.clip(np.log10(t1 / lo) / np.log10(hi / lo), 0, 1))
                row.append(f"{w} {lo / t1:.3g}-{hi / t1:.3g} below:{100 * share:.1f}%")
            say(f"  {aa} {T} K: " + "; ".join(row))
    say("  Avetisov Eq. 1 with the field only, K_r = 2 k_rac, AF = g*k2*a (rate units); eta_ss solves -K eta + G(1-eta^2)=0")
    worst = 0
    for lbl, g, ka in (("g=1e-17, k2a=3e-4", 1e-17, 3e-4), ("g=2.64e-16, k2a=3e-4", 2.64e-16, 3e-4),
                       ("g=2.64e-16, k2a=1e-4 (moons)", 2.64e-16, 1e-4), ("g=4.4e-18, k2a=1e-13", 4.4e-18, 1e-13)):
        K, G = 2 * K273["Ala"], g * ka
        eta = (-K + np.sqrt(K * K + 4 * G * G)) / (2 * G)
        worst = max(worst, eta)
        say(f"  {lbl}: AF/RF = {G / K:.3e}  eta_ss = {eta:.3e}")
    say(f"  largest eta_ss over these corners of Table 1: {worst:.3e}")
    return ts


# ------------------------------------------------------------------------------------------------------------- Q2
def q2():
    say("\n=== Q2. Autocatalysis on (open, fed), repo constants k=ki=kd=f=1, k0=1e-3, a0=1.6, g=0 ===")
    rows = [0, 0.001, 0.0013, 0.005, 0.01, 0.02, 0.023, 0.05, 0.08, 0.085, 0.09, 0.1, 0.2]
    say("  r    exists  eta*(formula)  eta*(ODE from pure)  max Re(eig) pure   max Re(eig) racemic   fsolve-vs-formula")
    for r in rows:
        p = par(r=r)
        q = asym(p)
        lam_r = max(jac_eigs(p, racemic(p)).real)
        if q:
            y0 = pure_state(p)
            sol = integrate(p, [y0[0], y0[1] * 0.99, y0[2] + y0[1] * 0.01], T=20000)
            ye = sol.y[:, -1]
            eode = (ye[1] - ye[2]) / (ye[1] + ye[2])
            fs = fsolve(lambda y: rhs(0, y, p), y0, xtol=1e-14)
            assert np.allclose(fs, y0, rtol=1e-8, atol=1e-12), (r, fs, y0)
            lam_p = max(jac_eigs(p, y0).real)
            say(f"  {r:<6g} yes   {np.sqrt(q[2]):.6f}   {eode:.6f}   {lam_p:+.3e}   {lam_r:+.3e}   ok")
            assert lam_p < 0 and lam_r > 0 and abs(eode - np.sqrt(q[2])) < 1e-6, (r, lam_p, lam_r, eode)
        else:
            say(f"  {r:<6g} no     -              -                    -                  {lam_r:+.3e}")
            assert lam_r < 0, (r, lam_r)
    rc = first_crossing(lambda r: eta2_of(par(r=r)), np.linspace(1e-6, 0.2, 400), 0)
    r90 = first_crossing(lambda r: eta2_of(par(r=r)), np.linspace(1e-6, 0.2, 400), 0.81)
    r99 = first_crossing(lambda r: eta2_of(par(r=r)), np.linspace(1e-6, 0.2, 4000), 0.9801)
    say(f"  r_c (pure branch exists) = {rc:.5f} kd; r(eta*>=0.9) = {r90:.5f}; r(eta*>=0.99) = {r99:.5f}")
    pc = par(r=rc)
    say(f"  at r_c: k a_rac = {pc['k'] * racemic(pc)[0]:.6f} vs kd + 2r = {pc['kd'] + 2 * rc:.6f} (pitchfork)")
    assert abs(pc["k"] * racemic(pc)[0] - (pc["kd"] + 2 * rc)) < 1e-6
    # survivors on both sides of r_c by full ODE from a 1% excess near the racemic state
    for lbl, r in (("0.8 r_c", 0.8 * rc), ("1.25 r_c", 1.25 * rc)):
        p = par(r=r)
        y = racemic(p)
        sol = integrate(p, [y[0], y[1] * 1.01, y[2] * 0.99], T=2e5)
        ye = sol.y[:, -1]
        say(f"  start at 1% excess, r = {lbl}: final eta = {(ye[1] - ye[2]) / (ye[1] + ye[2]):.6f}")
    # basin
    say("  basin (r = 0.01): initial excess eta0 -> final eta")
    p = par(r=0.01)
    y = racemic(p)
    for e0 in (1e-6, 1e-4, 1e-2, 0.5, 1.0, -1e-6, -0.5, 0.0):
        sol = integrate(p, [y[0], y[1] * (1 + e0), y[2] * (1 - e0)], T=20000)
        ye = sol.y[:, -1]
        say(f"    eta0 {e0:+.0e}: final eta {(ye[1] - ye[2]) / (ye[1] + ye[2]):+.6f}")
    # minimum strengths
    say("  minimum strengths (other constants fixed; thresholds: exists / eta*>=0.9 / eta*>=0.99)")
    thr = ((0.0, "exists"), (0.81, "eta*>=0.9"), (0.9801, "eta*>=0.99"))
    out = {}
    for r in (0.0, 0.01):
        for key, grid in (("k", np.geomspace(0.6, 1e3, 4000)), ("ki", np.geomspace(1e-5, 1e5, 4000)),
                          ("V", np.linspace(0.5, 1.0, 5001))):
            res = []
            for th, nm in thr:
                x = first_crossing(lambda x: eta2_of(par(r=r, **{key: x})), grid, th)
                res.append(x)
                out[r, key, nm] = x
            say(f"    r = {r}: {key}_min  " + "  ".join(f"{nm}: " + ("none" if x is None else f"{x:.4g}")
                                                        for (th, nm), x in zip(thr, res)))
        say(f"    r = {r}, ki = 0: pure state exists? {asym(par(r=r, ki=0.0)) is not None}")
        assert asym(par(r=r, ki=0.0)) is None
    kmin = out[0.01, "k", "eta*>=0.9"]
    say(f"    r = 0.01, k_min(eta*>=0.9) = {kmin:.4g}: k a0/(kd+2r) = {kmin * 1.6 / 1.02:.3f}")
    # draft's rate terms and Table-1 shares
    m_ex, m_90 = 1.6 / (2 * rc), 1.6 / (2 * r90)
    say(f"  pure state needs k a0 >= {m_ex:.2f} x (2r) to exist and >= {m_90:.2f} x (2r) for eta*>=0.9 (draft test: 1 x)")
    rng = np.random.default_rng(0)
    n = 200000
    k2 = 10 ** rng.uniform(-6, 0, n)
    say("  share of Table-1 draws with k2 a above m x 2 k_rac(Ala, 273 K), % (m = 1, exists, eta*>=0.9)")
    for w, (lo, hi) in (("Moons (Enceladus/Europa a)", (1e-12, 1e-4)), ("Early Earth", (1e-7, 3e-4))):
        a = 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)
        say(f"    {w}: " + "  ".join(f"{100 * np.mean(k2 * a > m * 2 * K273['Ala']):.2f}" for m in (1, m_ex, m_90)))
    # mirror
    for g in (0.0, 1e-17):
        p = par(r=0.01, g=g)
        y0 = pure_state(p)
        s1 = integrate(p, [y0[0], y0[1], y0[2]], T=200)
        s2 = integrate(p, [y0[0], y0[2], y0[1]], T=200)
        tol = 1e-12 if g == 0 else 1e-9
        e1 = (s1.y[1, -1] - s1.y[2, -1]) / (s1.y[1, -1] + s1.y[2, -1])
        e2 = (s2.y[1, -1] - s2.y[2, -1]) / (s2.y[1, -1] + s2.y[2, -1])
        assert abs(e1 + e2) < tol and abs(s1.y[0, -1] - s2.y[0, -1]) < tol, (g, e1, e2)
        say(f"  mirror start (g = {g:g}): eta {e1:+.9f} vs {e2:+.9f}  OK")
    return rc, r90


# ------------------------------------------------------------------------------------------------------------- Q3
def S_of(gam):
    return np.inf if gam >= 1 else (1 + gam) / (1 - gam)


def gam_of(S):
    return 1.0 if np.isinf(S) else (S - 1) / (S + 1)


def q3a():
    say("\n=== Q3A. Autocatalysis off: closed pool, racemization r = 0.01 + consumer c = 0.1, start eta = 1 ===")
    say("  time to eta_right < 1% (units 1/k), and ratio to racemization-only (= ln100/(2r) = 230.26)")
    base_t = np.log(100) / 0.02
    ev = [lambda t, y, p: eta_right(y, p["sg"]) - 0.5, lambda t, y, p: eta_right(y, p["sg"]) - 0.01]
    ev[0].terminal, ev[0].direction = False, -1
    ev[1].terminal, ev[1].direction = True, -1
    res = {}
    for N, lbl in ((0, "always on"), (40, "gate N=40"), (130, "gate N=130")):
        row = []
        for S in (1, 2, 10, 100, 1e4, np.inf):
            p = par(r=0.01, k=0, ki=0, kd=0, f=0, k0=0, c=0.1, gam=gam_of(S), N=N)
            sol = integrate(p, [0, 1.0, 0.0], T=3000, events=ev)
            t1 = sol.t_events[1][0] if len(sol.t_events[1]) else np.inf
            res[N, S] = t1 / base_t
            row.append(f"S={S:g}: {t1 / base_t:.4f}")
        say(f"  {lbl}: " + "  ".join(row))
        assert abs(res[N, 1] - 1) < 1e-4
        vals = [res[N, S] for S in (1, 2, 10, 100, 1e4, np.inf)]
        assert all(vals[i + 1] <= vals[i] + 1e-9 for i in range(5)), vals
    say("  no S >= 1 keeps the excess: every run above crossed 1% (finite time); a critical S >= 1 does not exist")
    assert all(np.isfinite(v) for v in res.values())
    # mirror
    p = par(r=0.01, k=0, ki=0, kd=0, f=0, k0=0, c=0.1, gam=gam_of(10), N=40)
    pm = dict(p, sg=-1)
    s1, s2 = integrate(p, [0, 1.0, 0.0], T=100), integrate(pm, [0, 0.0, 1.0], T=100)
    assert np.allclose(s1.y[1, -1], s2.y[2, -1], rtol=1e-9) and np.allclose(s1.y[2, -1], s2.y[1, -1], rtol=1e-9, atol=1e-18)
    say("  mirror start (D pure, D-preferring consumer): identical to the L run with L and D exchanged  OK")


def persists(p0, gam, c, N, sg=1):
    """Switch consumer on at the pure state; return final eta_right and a convergence flag."""
    p = dict(p0, gam=gam, c=c, N=N, sg=sg)
    y0 = pure_state(p0)
    if sg < 0:
        y0 = [y0[0], y0[2], y0[1]]
    sol = integrate(p, y0, T=T_END)
    ye = sol.y[:, -1]
    drift = max(abs(np.array(rhs(0, ye, p))) / (1 + abs(ye)))
    return eta_right(ye, sg), drift < 1e-7


def q3b():
    say("\n=== Q3B. Autocatalysis on: consumer switched on at the pure state (r = 0.01, repo constants) ===")
    say("  final eta_right at t = 5000 (>= 0.5 persists, <= -0.5 flipped to the other hand, else racemic); drift flag = converged")
    p0 = par(r=0.01)
    say(f"  pure state before the consumer: eta* = {np.sqrt(eta2_of(p0)):.6f}")
    Ss = (1, 1.1, 1.25, 1.5, 2, 3, 5, 10, 30, 100, 1e3, 1e4, 1e6, np.inf)
    summary = {}
    for N, lbl in ((0, "always on (N=0)"), (10, "gate N=10"), (40, "gate N=40"), (130, "gate N=130")):
        say(f"  -- {lbl}")
        for c in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 100):
            fin = []
            for S in Ss:
                e, ok = persists(p0, gam_of(S), c, N)
                assert ok, (N, c, S, e)
                fin.append(e)
            keep = [e >= 0.5 for e in fin]
            if all(keep):
                verdict = "persists for every S up to infinity (no critical S)"
            elif not keep[0]:
                verdict = f"lost already at S = 1 (final {fin[0]:+.3f}); at S = inf final {fin[-1]:+.3f}"
            else:
                i = keep.index(False)
                lo, hi = gam_of(Ss[i - 1]), gam_of(Ss[i])
                for _ in range(40):
                    mid = (lo + hi) / 2
                    if persists(p0, mid, c, N)[0] >= 0.5:
                        lo = mid
                    else:
                        hi = mid
                Sc = S_of((lo + hi) / 2)
                mono = all(not x for x in keep[i:])
                verdict = f"critical S_c = {Sc:.4g}; beyond it final at S = inf {fin[-1]:+.3f}; monotone: {mono}"
            summary[N, c] = verdict
            say(f"    c = {c:<5g} eta_right(S=1, 3, 10, 100, inf) = " + ", ".join(f"{fin[Ss.index(s)]:+.3f}" for s in (1, 3, 10, 100, np.inf))
                + f"  -> {verdict}")
    # mirror
    for N in (0, 40):
        e1, _ = persists(p0, gam_of(10), 0.1, N, sg=1)
        e2, _ = persists(p0, gam_of(10), 0.1, N, sg=-1)
        assert abs(e1 - e2) < 1e-9, (N, e1, e2)
    say("  mirror start (D-pure pool with a D-preferring consumer) gives the same eta_right: OK")
    say("  enzyme-like selectivity: not sourced (no value in journals/ija_draft.md); no comparison made.")
    say("  For reference only: the chapter's own chain criterion (Eq. 9, gamma > 1 - 2/N) corresponds to S about N-1 ="
        f" {40 - 1} (N = 40) and {130 - 1} (N = 130).")
    return summary


if __name__ == "__main__":
    q1()
    q2()
    q3a()
    q3b()
    say("\nALL ASSERTS PASSED")
    sys.exit(0)
