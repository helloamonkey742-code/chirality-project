"""Checks added after expert feedback (A. Brandenburg, J. Dworkin, 2026-09-30).

R1 scale    rule of thumb: how many molecules must decide together, by concentration and temperature?
R2 Bennu    Bennu/Ryugu amino acids are racemic, and Bennu is plausibly a piece of a wet parent body.
            If the same amplifying chemistry had run there, would it have finished (and left big excesses
            in gram-sized samples)? Condition the uncertainty sweep on "Bennu stays racemic".
R3 network  connected compartments (Russell-style vent mound) vs one isolated pool: normal form with
            M compartments exchanging at rate q.
R4 long run does the outcome drift or flip if the simulation runs 10-100x longer? Flip rate vs noise.
R5 design   experiment sizing at meteoritic concentrations (<= ~1 mM) and published ee precision.
R6 closed   R2 assumed any amplifier that "finished" leaves a lasting excess. That holds only in an open,
            driven system (an ocean fed by vents). A closed rock whose drive runs down must return to 50/50
            (thermodynamics: at equilibrium both hands are equal). Test it, then redo R2 with that rule.
R7 flipping molecules flipping hand at random (racemization, both ways at rate r): does it bias the outcome,
            or only add noise and delay the choice? Normal form gets -2r*alpha drift and noise eps*(1+2r).
R8 beta     beta-decay electrons spin one way and may destroy one hand faster (Vester-Ulbricht; Dreiling &
            Gay 2014: asymmetry A ~ 3e-4, sub-eV electrons, bromocamphor gas). Destruction at rate d (units
            of k) with asymmetry A adds a bias g_beta = d*A/2. How big must d be to rival the PVED?
"""
import numpy as np
from scipy.stats import norm
from scipy.integrate import solve_ivp
from ocean import NA, YR, PREF, delta_at, BODIES
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
#   amino acids 10-330 nmol/g rock, bracketing Bennu (70) and Murchison (253) (Glavin et al. 2025; CR chondrites
#   reach 3,300, Aponte et al. 2020), water/rock 0.1-1 by mass (Lee et al. 2025) -> 1e-5 .. 3e-3 M in the fluid
#   (Bennu alone: 7e-5 .. 7e-4 M); aqueous alteration lasting 1-10 Myr. quench.py Q5 tests other choices.
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
    print(f"  {'conc':>7} {'seed':>5} " + " ".join(f"{s:>10}" for s in ('0.2%', '1%', '2.6%', '2.6% 9 vials')))
    for c in (1e-4, 1e-3, 1e-2, 1e-1):
        for e0 in (0.05, 0.2):
            row = [k2_min(c, YR, e0, s) for s in (0.002, 0.01, 0.026)] + [k2_min(c, YR, e0, 0.026 / 3)]
            print(f"  {c:7g} {e0:5.0%} " + " ".join(f"{v:10.0e}" for v in row))
    pond_1mM = 10 / (1e-3 * 1e2 * YR)
    print(f"  nature's floors at realistic concentrations: pond (1 mM, 100 yr) {pond_1mM:.0e}; ocean (1 uM, 1e8 yr) "
          f"{10 / (1e-6 * 1e8 * YR):.0e}")
    assert k2_min(1e-3, YR, 0.2, 0.026 / 3) > pond_1mM, "a 1 mM, 1 yr run cannot reach the pond floor"
    best = k2_min(1e-3, 3 * YR, 0.2, 0.026 / 3)
    print(f"  realistic arm: 1 mM, 20% seed, 3 yr, 9 independent vials at 2.6% -> {best:.1e} (pond floor {pond_1mM:.1e})")
    assert best < pond_1mM
    # Dworkin et al. 2024 test 5: repeat injections of one vial are not replicates, so noise falls as 1/sqrt(vials)
    three = k2_min(1e-3, 3 * YR, 0.2, 0.026 / np.sqrt(3))
    print(f"  same arm with 3 vials -> {three:.1e}: only {1 - three / pond_1mM:.0%} below the floor, so the plan uses 9 vials")
    assert pond_1mM / three < 1.1, "3-vial margin no longer thin; revisit the 9-vial choice"


def frank_closed(a0, slow=1.0, feed=0.0, t_end=1e8, seed=0.001, kw=0.0, full=False):
    """frank2.py's network with the A<->L and A<->D steps (uncatalysed and autocatalytic) reversible and obeying
    detailed balance: both routes share K = k0/k0r = k/kr = 2/slow (no perpetual A -> L -> A cycle), and L and D
    share the same constants, so equilibrium has no preferred hand. slow scales the reverse rates only. The mutual
    antagonism L + D -> W is one-way unless kw > 0 (W -> L + D), so a closed run drains most of its mass into W
    while the excess rises and fades.
    feed > 0 makes it open (feedstock topped up to a0, everything flows out at rate feed).
    Returns times and ee, starting from seed (default 0.1%); full=True returns times and (a, L, D, W) instead.
    Time is in 1/(k * one concentration unit), with a0 = 200 units (how that maps to k2*c*tau: quench.py Q5)."""
    k, ki, k0, k0r, kr = 1.0, 1.0, 1e-3, 5e-4 * slow, 0.5 * slow

    def rhs(t, y):
        a, L, D, W = y
        nl = k0 * a - k0r * L + k * a * L - kr * L * L          # net A -> L
        nd = k0 * a - k0r * D + k * a * D - kr * D * D
        w = ki * L * D - kw * W                                  # net L + D -> W
        return [-nl - nd + feed * (a0 - a), nl - w - feed * L, nd - w - feed * D, w - feed * W]
    t = np.logspace(-1, np.log10(t_end), 1500)
    y = solve_ivp(rhs, (0, t_end), [a0, 0.01 * (1 + seed), 0.01 * (1 - seed), 0.0], method="LSODA", t_eval=t, rtol=1e-10, atol=1e-14).y
    return (t, y) if full else (t, (y[1] - y[2]) / (y[1] + y[2]))


def r6_closed():
    print("\nR6. Closed rock vs open ocean: does an amplified excess last? (frank2 network, A<->L/D reversible)")
    spans = {}
    for slow, label in ((1.0, "s = 1: K = 2"), (0.01, "s = 0.01: K = 200 (slow reverse)")):
        t, ee = frank_closed(200.0, slow)
        on = t[np.abs(ee) > 0.1 * np.abs(ee).max()]
        spans[slow] = on[0], on[-1]
        print(f"  closed, {label:34}: peak |ee| {np.abs(ee).max():.2f}, above 10% of peak from ~{on[0]:.0e} to ~{on[-1]:.0e} "
              f"(model time ~ k2 c tau, see quench.py Q5); final {ee[-1]:+.3f}")
        assert abs(ee[-1]) < 0.01, "closed system should end racemic"
    t, ee = frank_closed(2.0, 1.0, feed=0.05)
    print(f"  open (fed + outflow)                           : final ee {ee[-1]:+.3f}  (held as long as the drive lasts)")
    assert ee[-1] > 0.99
    print("  -> Bennu's parent body (closed, heat gone after ~10 Myr) stays racemic if its excess never rose")
    print("     (k2 c tau < rise time) OR rose and then faded (k2 c tau > fade time). Redo R2 with that rule:")
    for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):
        for slow, (rise, fade) in spans.items():
            row = []
            for name, args in U.WORLDS.items():
                x, win, _ = U.sweep(*args, k2_rng=k2_rng)
                kct = x["rate k2"] * 10 ** U.RNG.uniform(*np.log10(BENNU_C), U.N) * 10 ** U.RNG.uniform(*np.log10(BENNU_T), U.N) * YR
                rac = (kct < rise) | (kct > fade)
                row.append(f"{name.split()[0]} {(win & rac).mean() / max(rac.mean(), 1 / U.N):4.0%}")
            print(f"    k2 {k2_rng[0]:g}..1, excess {rise:.0e}..{fade:.0e}: Bennu racemic is consistent; P(win | Bennu racemic): "
                  + ", ".join(row))
    return spans


def sweep_sim(g, eps, gamma, r=0.0, runs=4000, dt=0.01, seed=0):
    """sim.simulate's sweep (lam from -1 to +1) with symmetric flipping at rate r per molecule:
    drift -2r*alpha, noise eps*(1+2r) (flips are extra independent events). Returns P(favoured), mean |alpha|."""
    rng = np.random.default_rng(seed)
    a = np.zeros(runs)
    for t in np.arange(-1 / gamma, 1 / gamma, dt):
        a += ((gamma * t - 2 * r) * a - a**3 + g) * dt + np.sqrt(eps * (1 + 2 * r) * dt) * rng.standard_normal(runs)
    return (a > 0).mean(), np.abs(a).mean()


def r7_flipping():
    print("\nR7. Molecules flipping hand at random (racemization at rate r, units of k; g = 1e-3, eps = 1e-4)")
    g, eps, gamma = 1e-3, 1e-4, 0.01
    print(f"  {'r':>5} {'predicted P':>11} {'simulated P':>11} {'final |ee|':>10}")
    for r in (0.0, 0.05, 0.2, 0.45, 0.6):
        p, ee = sweep_sim(g, eps, gamma, r)
        pred = norm.cdf(delta(g, eps * (1 + 2 * r), gamma))
        print(f"  {r:5} {pred if 2 * r < 1 else float('nan'):11.3f} {p:11.3f} {ee:10.3f}" + ("  (no choice made)" if 2 * r >= 1 else ""))
        if 2 * r < 0.9:
            assert abs(p - pred) < 0.03, "flipping should only add noise, not bias"
        else:
            assert ee < 0.2, "flipping faster than amplification should leave the mixture racemic"
    print("  -> flipping is symmetric: no new bias. It adds noise (P drifts toward 50%) and, if faster than the")
    print("     amplifier (2r > lam_max), stops the choice altogether. The weak-force tilt is the only bias.")


def r8_beta():
    print("\nR8. Beta-decay electrons as a second weak-force bias (g_beta = d*A/2, d = destruction rate / k)")
    g_pv, eps, gamma = 1e-3, 1e-4, 0.01           # test case at simulable scale
    # check the mapping: destruction asymmetry in the drift == a shift of g
    rng = np.random.default_rng(3)
    d, A = 0.02, 0.1                               # g_beta = 1e-3, same size as g_pv
    a = np.zeros(4000)
    for t in np.arange(-1 / gamma, 1 / gamma, 0.01):
        a += (gamma * t * a - a**3 + g_pv + d * A / 2 * (1 - a**2)) * 0.01 + np.sqrt(eps * 0.01) * rng.standard_normal(a.size)
    p_sim, p_pred = (a > 0).mean(), norm.cdf(delta(g_pv + d * A / 2, eps, gamma))
    print(f"  check (g_beta = g_pv = 1e-3, aligned): simulated P {p_sim:.3f}, predicted {p_pred:.3f}")
    assert abs(p_sim - p_pred) < 0.03
    A = 3e-4                                       # Dreiling & Gay 2014, sub-eV electrons, bromocamphor
    d_match = 2 * 1e-17 / A
    print(f"  real scale: PVED g = 1e-17. With A = {A:.0e}, beta matches the PVED once d >= {d_match:.0e},")
    print("  i.e. once about 1 molecule in 1e13 is destroyed by polarized electrons per reaction time.")
    print("  P(physics) for a borderline ocean (PVED alone gives Delta = 1, P = 0.84), beta aligned / opposed:")
    for dd in (0.0, 1e-14, d_match, 1e-12):
        gb = dd * A / 2
        ps = [norm.cdf(1.0 * (1e-17 + sgn * gb) / 1e-17) for sgn in (1, -1)]   # Delta scales with total g
        print(f"    d = {dd:7.0e}: g_beta = {gb:.0e}; aligned P = {ps[0]:.3f}, opposed P = {ps[1]:.3f}")
    print("  -> if the lab asymmetry carried over to amino acids in water, even feeble radiolysis would rival the PVED,")
    print("     and the answer would hinge on the unknown sign of A. A for amino acids, at beta (keV-MeV) energies, in")
    print("     water, is unmeasured, so this stays a possible extra bias, not part of the main result.")


if __name__ == "__main__":
    r1_scale()
    r2_bennu()
    r3_network()
    r4_longrun()
    r5_design()
    r6_closed()
    r7_flipping()
    r8_beta()
