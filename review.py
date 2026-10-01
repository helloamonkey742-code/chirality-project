"""Checks added after expert feedback (A. Brandenburg, J. Dworkin, 2026-09-30).

R1 scale    rule of thumb: how many molecules must decide together, by concentration and temperature?
R2 Bennu    Bennu/Ryugu amino acids are racemic, and Bennu is plausibly a piece of a wet parent body.
            If the same amplifying chemistry had run there, would it have finished (and left big excesses
            in gram-sized samples)? Condition the uncertainty sweep on "Bennu stays racemic".
R3 network  connected compartments (Russell-style vent mound) vs one isolated pool: normal form with
            M compartments exchanging at rate q.
R4 long run does the outcome drift or flip if the simulation runs 10-100x longer? Flip rate vs noise.
R5 design   experiment sizing at meteoritic concentrations (<= ~1 mM) and published ee precision.
"""
import numpy as np
from scipy.stats import norm
from ocean import NA, YR, PREF, delta_at
from sim import delta
from design import k2_min
import uncertainty as U

R_GAS = 8.314


def n_needed(g, k_tau, target=2.0):
    """Molecules that must decide together for Delta >= target (from Delta = PREF*g*(k tau)^1/4*sqrt(N))."""
    return (target / (PREF * g * k_tau**0.25)) ** 2


def r1_scale():
    print("R1. How many molecules must decide together? (Delta >= 2, physics wins 97.7%)")
    # the PVED is an energy; the bias is that energy over kT, so colder helps: g ~ 1/T
    dE = 1e-17 * R_GAS * 273          # J/mol, the g = 1e-17 at 0 C used throughout
    print(f"  PVED taken as {dE:.1e} J/mol (g = 1e-17 at 0 C); modern range spans ~1e-19..1e-16 in g")
    for T in (273, 300, 373):
        g = dE / (R_GAS * T)
        n = n_needed(g, 10)
        vols = {c: n / (c * 1000 * NA) for c in (1e-6, 1e-3, 1e-1)}
        print(f"  T={T} K: g={g:.1e}, need N >= {n:.1e} molecules ({n / NA:.1e} mol); "
              + ", ".join(f"at {c:g} M a well-mixed {v:.0e} m^3 (cube {v ** (1 / 3) / 1e3:.2g} km)" for c, v in vols.items()))
    n0 = n_needed(1e-17, 10)
    assert abs(delta_at(1e-6, 1e-17, 1e-3, n0 / (1e-6 * 1000 * NA), 10 / (1e-3 * 1e-6)) - 2) < 1e-6
    print("  -> a warm little pond (1 m^3 at 1 M holds 6e26) is ~7 orders short; only ocean-scale water that")
    print("     is mixed on the decision time can pool enough molecules. Time helps only as (k tau)^1/4.")
    return n0


# Bennu parent body (all ASSUMED ranges; see paper Sec. 3.10):
#   amino acids 10-330 nmol/g rock (Dworkin, pers. comm.: LAP 02342 D+L 330 nmol/g), water/rock 0.1-1 by mass
#   -> 1e-5 .. 3e-3 M in the fluid; aqueous alteration lasting 1-10 Myr.
BENNU_C, BENNU_T = (1e-5, 3e-3), (1e6, 1e7)


def r2_bennu():
    print("\nR2. Bennu/Ryugu are racemic. Would the same amplifier have finished there?")
    ct_lo, ct_hi = BENNU_C[0] * BENNU_T[0] * YR, BENNU_C[1] * BENNU_T[1] * YR
    print(f"  Bennu parent body c*tau = {ct_lo:.0e} .. {ct_hi:.0e} M s -> any amplifier with "
          f"k2 >= {10 / ct_hi:.0e} .. {10 / ct_lo:.0e} /M/s would have finished there")
    # patch size left in Bennu if it had finished: diffusion only (pore water), l = 5.5*sqrt(D t)
    l = 5.5 * np.sqrt(1e-9 * np.array(BENNU_T) * YR)
    print(f"  patches after alteration (pore diffusion only): {l[0]:.0f}-{l[1]:.0f} m across -> a gram sample would sit")
    print("  inside one patch and show a large excess of either hand. Bennu's chiral amino acids (incl. isovaline,")
    print("  which cannot racemize afterwards) are racemic -> no such amplifier finished there.")
    out = {}
    for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):   # paper's range first (Enceladus shares uncertainty.py's draws; later worlds differ), then extended down
        print(f"  rate k2 drawn over {k2_rng[0]:g}..{k2_rng[1]:g} /M/s:")
        for name, args in U.WORLDS.items():
            x, win, _ = U.sweep(*args, k2_rng=k2_rng)
            cB = 10 ** U.RNG.uniform(*np.log10(BENNU_C), U.N)
            tB = 10 ** U.RNG.uniform(*np.log10(BENNU_T), U.N) * YR
            racemic = x["rate k2"] * cB * tB < 10               # Bennu's reaction never finished
            both = (win & racemic).mean()
            cond = both / max(racemic.mean(), 1 / U.N)
            ct = x["concentration"] * x["time"]
            out[name, k2_rng] = (win.mean(), racemic.mean(), both, cond)
            print(f"    {name:18} wins {win.mean():4.0%} | Bennu stays racemic {racemic.mean():4.0%} | both {both:5.1%} | "
                  f"P(win | Bennu racemic) {cond:4.0%} | c*tau above Bennu's max in {(ct > ct_hi).mean():.0%}")
    assert all(v[2] <= v[0] for v in out.values())
    return out


def network(M, q, g, eps_each, gamma, runs=3000, dt=0.02, seed=0, t_hold=200):
    """M compartments, each the normal form with noise eps_each, exchanging at rate q (well-mixed coupling)."""
    rng = np.random.default_rng(seed)
    a = np.zeros((runs, M))
    for t in np.arange(-1 / gamma, 1 / gamma + t_hold, dt):
        lam = min(gamma * t, 1.0)
        a += (lam * a - a**3 + g + q * (a.mean(1, keepdims=True) - a)) * dt + np.sqrt(eps_each * dt) * rng.standard_normal(a.shape)
    m = np.sign(a).mean(1)
    return (m > 0).mean() + 0.5 * (m == 0).mean(), (np.abs(m) == 1).mean()


def r3_network():
    print("\nR3. Connected compartments vs isolated ones (normal form, M = 16, g = 1e-3, gamma = 0.01)")
    M, g, gamma, eps = 16, 1e-3, 0.01, 1e-4 * 16      # each compartment has 1/16 of the molecules
    p_one = norm.cdf(delta(g, eps / M, gamma))
    p_each = norm.cdf(delta(g, eps, gamma))
    print(f"  one pool with all the molecules: P = {p_one:.3f};  a single isolated compartment: P = {p_each:.3f}")
    rows = []
    for q in (0.0, 0.01, 0.1, 1.0):
        p, homo = network(M, q, g, eps, gamma)
        rows.append((q, p, homo))
        print(f"  exchange rate q = {q:<5} (x decision time {q / np.sqrt(gamma):5.2f}): P(majority favoured) = {p:.3f}, "
              f"all 16 one hand = {homo:.0%}")
    assert abs(rows[-1][1] - p_one) < 0.03, "well-connected network should act as one pool"
    assert rows[0][2] < 0.05 and rows[-1][2] > 0.95, "isolation should leave mixed hands, connection one hand"
    # real vent mound: pores ~1 mm linked by diffusion (q = D/l^2 = 1e-3 /s); decision time sqrt(tau/k)
    k2, c, tau, D = 1e-3, 1e-3, 1e5 * YR, 1e-9
    t_dec = np.sqrt(tau / (k2 * c))
    reach = np.sqrt(D * t_dec)
    print(f"  vent mound, k2={k2:g}, c={c:g} M, 1e5 yr: decision time {t_dec / YR:.0f} yr, pore exchange x decision "
          f"time {1e-3 * t_dec:.0e} (fully connected); diffusion reach {reach:.1f} m")
    for label, V in (("single 1 mm^3 pore", 1e-9), (f"diffusion-connected ({reach:.1f} m)^3", reach**3),
                     ("flow-connected 100 m mound", 1e6), ("flow-connected 1 km^3", 1e9)):
        print(f"    {label:32} P(physics) = {norm.cdf(delta_at(c, 1e-17, k2, V, tau)):.3f}")
    return rows


def r4_longrun():
    print("\nR4. Long runs: does the decided hand drift or flip later?")
    g, eps, gamma = 1e-3, 1e-4, 0.01
    rng = np.random.default_rng(0)
    a = np.zeros(4000)
    dt, t, marks = 0.01, -1 / gamma, {}
    for T in (1 / gamma, 1 / gamma + 1e3, 1 / gamma + 1e4):     # hold lam = 1 (reaction saturated) afterwards
        while t < T:
            lam = min(gamma * t, 1.0)
            a += (lam * a - a**3 + g) * dt + np.sqrt(eps * dt) * rng.standard_normal(a.size)
            t += dt
        marks[T] = (a > 0).mean()
    print("  P(favoured) at end of sweep / +1000 / +10000 time units: " + " / ".join(f"{p:.3f}" for p in marks.values()))
    ps = list(marks.values())
    assert max(ps) - min(ps) < 0.01, "outcome drifted after the choice"
    print("  flip rate once decided (lam = 1, no bias) vs noise eps = 1/N:")
    for e in (0.2, 0.15, 0.12, 0.1):
        rng = np.random.default_rng(1)
        a, side, flips, T = np.ones(2000), np.ones(2000), 0, 2000.0
        for _ in range(int(T / 0.01)):
            a += (a - a**3) * 0.01 + np.sqrt(e * 0.01) * rng.standard_normal(a.size)
            new = np.where(a > 0.5, 1, np.where(a < -0.5, -1, side))
            flips += (new != side).sum(); side = new
        kramers = np.sqrt(2) / (2 * np.pi) * np.exp(-1 / (2 * e))
        print(f"    eps={e:<5} measured {flips / (2000 * T):.1e} per unit time, Kramers {kramers:.1e}")
        assert 0.3 < flips / (2000 * T) / kramers < 3
    print("  rate ~ exp(-N/2): 1 um^3 at 1 uM holds ~600 molecules (flips possible); any pond or ocean never flips.")


def r5_design():
    print("\nR5. Experiment at meteoritic concentrations: smallest detectable k2 (/M/s) in 1 year")
    print("  (ee noise per sample: 0.2% optimistic GC-MS; 1% typical; 2.6% Glavin & Dworkin 2009 replicates)")
    print(f"  {'conc':>7} {'seed':>5} " + " ".join(f"{s:>10}" for s in ('0.2%', '1%', '2.6%', '2.6% n=9')))
    for c in (1e-4, 1e-3, 1e-2, 1e-1):
        for e0 in (0.05, 0.2):
            row = [k2_min(c, YR, e0, s) for s in (0.002, 0.01, 0.026)] + [k2_min(c, YR, e0, 0.026 / 3)]
            print(f"  {c:7g} {e0:5.0%} " + " ".join(f"{v:10.0e}" for v in row))
    pond_1mM = 10 / (1e-3 * 1e2 * YR)
    print(f"  nature's floors at realistic concentrations: pond (1 mM, 100 yr) {pond_1mM:.0e}; ocean (1 uM, 1e8 yr) "
          f"{10 / (1e-6 * 1e8 * YR):.0e}")
    assert k2_min(1e-3, YR, 0.2, 0.026 / 3) > pond_1mM, "a 1 mM, 1 yr run cannot reach the pond floor"
    best = k2_min(1e-3, 3 * YR, 0.2, 0.026 / 3)
    print(f"  realistic arm: 1 mM, 20% seed, 3 yr, 9 replicates at 2.6% -> {best:.1e} (pond floor {pond_1mM:.0e})")
    assert best < pond_1mM


if __name__ == "__main__":
    r1_scale()
    r2_bennu()
    r3_network()
    r4_longrun()
    r5_design()
