"""Quench-time sweep (2026-10-01, prompted by LPSC 2024 abstract #2645).

R6 used a step rule: a closed parent body is racemic if k2*c*tau falls before its excess rises or after it
fades. But Bennu's fluid left partway through percolation (sulfur-rich solvent front frozen in a ~2 mm
grain, doi:10.1111/maps.14335), and solid-state reactions are negligible, so the excess is FROZEN at
whatever the closed network had reached when the water left. Here the frozen value is read off the full
closed-network curve (review.frank_closed) at k2*c*tau, instead of the step rule.

Questions:
 Q1 What share of Bennu-like draws freeze racemic (|ee| < 5.2%, i.e. 2 x the 2.6% measurement error),
    inside the isovaline band seen in meteorites (5-20%; Dworkin et al. 2024: L-Iva 0-20% in 17 meteorites),
    or above it (> 20%)?
 Q2 How wide (in decades of k2*c*tau) is the window that freezes into the 5-20% band? A narrow window means
    reproducing the meteorite spread needs fine-tuned timing.
 Q3 Does P(win | Bennu racemic) change when the step rule is replaced by the frozen-curve rule?
 Q4 Does an excess inherited from before the parent body (20% seed) outlast the closed network?
 Q5 How robust is Q1 to what c means in model time, the racemic cut, the Bennu inputs, and the one-way waste
    step? And how much of the chiral material is still free (not waste) while an excess exists?
Units follow R6: model time (1/(k * one concentration unit), a0 = 200 units) read as k2*c*tau; seed 0.1%.
"""
import numpy as np
from review import frank_closed, BENNU_C, BENNU_T
from ocean import YR
import uncertainty as U

RAC, BAND = 0.052, 0.20          # 2 x 2.6% (Glavin & Dworkin 2009); upper edge of meteoritic Iva ee
failures = []


def assertTrue(c, msg):
    if not c:
        failures.append(msg)


def frozen(kct, t, ee):
    """|ee| frozen at k2*c*tau, read off the closed curve (seed before it starts, racemic after it ends)."""
    return np.interp(np.log10(kct), np.log10(t), np.abs(ee), left=0.001, right=0.0)


def shares(q, rac=RAC):
    """Shares of frozen |ee| that are racemic, inside the meteoritic band, and above it."""
    f = [(q < rac).mean(), ((q >= rac) & (q < BAND)).mean(), (q >= BAND).mean()]
    assertTrue(abs(sum(f) - 1) < 1e-12, "bands do not partition the draws")
    return f


def window(t, ee, lo, hi):
    """Decades of k2*c*tau whose frozen |ee| lies in [lo, hi)."""
    lt = np.log10(t)
    m = (np.abs(ee) >= lo) & (np.abs(ee) < hi)
    return np.sum(np.diff(lt)[m[:-1]])


rng = np.random.default_rng(1)
n = 200_000
curves, draws = {}, {}
for slow, label in ((1.0, "s = 1 (K = 2)"), (0.01, "s = 0.01 (K = 200, reverse 100x slower)")):
    t, ee = frank_closed(200.0, slow)
    curves[slow] = (t, ee)
    a = np.abs(ee)
    assertTrue(abs(frozen(1e-3, t, ee) - 0.001) < 1e-9, f"{label}: frozen value before t[0] is not the seed")
    assertTrue(a[-1] < 0.01, f"{label}: closed curve does not end racemic")
    assertTrue(a.max() > BAND, f"{label}: curve never exceeds the band, Q1 is vacuous")
    w_band, w_hi = window(t, ee, RAC, BAND), window(t, ee, BAND, 2.0)
    print(f"\n{label}: peak |ee| {a.max():.2f}; k2*c*tau window freezing into 5-20%: {w_band:.2f} decades; "
          f"into >20%: {w_hi:.2f} decades (of {np.log10(t[-1] / t[0]):.0f} simulated)")
    for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):
        k2 = 10 ** rng.uniform(*np.log10(k2_rng), n)
        kct = k2 * 10 ** rng.uniform(*np.log10(BENNU_C), n) * 10 ** rng.uniform(*np.log10(BENNU_T), n) * YR
        f = shares(frozen(kct, t, ee))
        draws[slow, k2_rng[0]] = k2, kct, f
        print(f"  k2 {k2_rng[0]:g}..1: frozen racemic {f[0]:6.1%} | in meteoritic band {f[1]:6.1%} | above band {f[2]:6.1%}"
              f"   (k2*c*tau spans {np.log10(kct.min()):.1f}..{np.log10(kct.max()):.1f} decades)")

print("\nQ4. Inherited excess (Cooper & Rios 2016: sugar-acid excesses may predate the parent body)")
for slow, (t0, ee0) in curves.items():
    t, ee = frank_closed(200.0, slow, seed=0.20, t_end=1e10)   # a bigger seed fades later than 1e8
    a, a0 = np.abs(ee), np.abs(ee0)
    ratio = t[a > 0.1 * a.max()][-1] / t0[a0 > 0.1 * a0.max()][-1]
    print(f"  20% seed in the closed network (slow={slow:g}): final {ee[-1]:+.3f}, fades {ratio:.1f}x later than a 0.1% seed")
    assertTrue(abs(ee[-1]) < 0.01, f"closed network should erase an inherited 20% seed too (slow {slow:g})")
    assertTrue(1 < ratio < 3, f"inherited seed changes the fade time by {ratio:.1f}x, outside 1-3x (slow {slow:g})")
print("  so an inherited excess survives only outside a reversible amplifying network, or if the water left first")

R6 = {(1e-6, 1.0): (.22, .29, .75), (1e-6, 0.01): (.23, .30, .75), (1e-12, 1.0): (.15, .21, .57), (1e-12, 0.01): (.17, .22, .55)}  # review.py R6 / PAPER table
print("\nQ3. P(win | Bennu frozen racemic), frozen-curve rule vs R6 step rule (same draws per row)")
for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):
    for slow, (t, ee) in curves.items():
        on = t[np.abs(ee) > 0.1 * np.abs(ee).max()]       # R6: excess above 10% of its peak from on[0] to on[-1]
        row = []
        for name, args in U.WORLDS.items():
            x, win, _ = U.sweep(*args, k2_rng=k2_rng)
            kct = x["rate k2"] * 10 ** rng.uniform(*np.log10(BENNU_C), U.N) * 10 ** rng.uniform(*np.log10(BENNU_T), U.N) * YR
            rac_q = frozen(kct, t, ee) < RAC
            rac_s = (kct < on[0]) | (kct > on[-1])
            pq = (win & rac_q).mean() / max(rac_q.mean(), 1 / U.N)
            ps = (win & rac_s).mean() / max(rac_s.mean(), 1 / U.N)
            row.append(f"{name.split()[0]} {pq:4.0%} (step {ps:4.0%})")
            ref = R6[(k2_rng[0], slow)][len(row) - 1]
            assertTrue(abs(ps - ref) < 0.02, f"step rule {ps:.0%} does not reproduce R6 {ref:.0%} ({name}, k2 {k2_rng[0]:g}, slow {slow:g})")
        print(f"  k2 {k2_rng[0]:g}..1, slow={slow:g}: " + ", ".join(row))

LIT_C, LIT_T = (7e-5, 7e-4), (1e2, 1e7)   # Bennu's 70 nmol/g at water/rock 0.1-1 (Glavin 2025; Lee 2025); water 100 yr..10 Myr
RAC2 = 2 * 0.062                          # 2 x Bennu's own isovaline ee error (Glavin 2025, Extended Data Table 4)
print("\nQ5. How robust is Q1? Frozen racemic | in meteoritic band | above band, same k2 draws as Q1")
for slow, (t, ee) in curves.items():
    _, (a, L, D, W) = frank_closed(200.0, slow, full=True)
    hot = np.abs(ee) >= RAC
    free = (L + D) / (L + D + 2 * W)                  # each W holds one L and one D
    bulk = np.abs(L - D) / (L + D + 2 * W)            # |ee| if the waste were hydrolysed and counted too
    x_pk = (L + D)[np.abs(ee).argmax()]
    tw, eew = frank_closed(200.0, slow, kw=1e-6)
    after = t > t[np.abs(ee) > 0.1 * np.abs(ee).max()][-1]   # once the excess has faded (R6 fade time)
    print(f"\n  slow={slow:g}: chiral material still free while |ee| >= {RAC:.1%}: at most {free[hot].max():.2%} "
          f"(Bennu: 45 of 70 nmol/g), after the excess fades: at most {100 * free[after].max():.1g}%; "
          f"|ee| if the waste were counted: at most {bulk.max():.2%}")
    assertTrue(free[after].max() < 1e-4, f"free chiral material after the fade is not ~0 (slow {slow:g})")
    assertTrue(bulk.max() < 0.01, f"counting the waste should leave the bulk racemic (slow {slow:g})")
    assertTrue(abs(eew[-1]) < 0.01, f"reversible waste: closed curve does not end racemic (slow {slow:g})")
    for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):
        k2, kct, f1 = draws[slow, k2_rng[0]]
        kct_lit = k2 * 10 ** rng.uniform(*np.log10(LIT_C), n) * 10 ** rng.uniform(*np.log10(LIT_T), n) * YR
        rows = {"baseline (= Q1)": shares(frozen(kct, t, ee)),
                "c = whole budget a0 (time /200)": shares(frozen(kct / 200, t, ee)),
                f"c = free L+D at peak (time x{1 / x_pk:.0f})": shares(frozen(kct / x_pk, t, ee)),
                f"racemic below {RAC2:.1%}": shares(frozen(kct, t, ee), RAC2),
                "Bennu-only c, water 1e2..1e7 yr": shares(frozen(kct_lit, t, ee)),
                "waste reversible (kw = 1e-6)": shares(frozen(kct, tw, eew))}
        assertTrue(rows["baseline (= Q1)"] == f1, f"Q5 baseline does not reproduce Q1 (slow {slow:g}, k2 {k2_rng[0]:g})")
        print(f"   k2 {k2_rng[0]:g}..1:")
        for name, f in rows.items():
            print(f"     {name:36} {f[0]:6.1%} | {f[1]:6.1%} | {f[2]:6.1%}")
        moved = max(abs(f[0] - f1[0]) for name, f in rows.items() if name.startswith(("c = ", "Bennu-only")))
        if k2_rng[0] == 1e-12:
            assertTrue(moved < 0.025,f"wide k2 range: racemic share moves {moved:.1%} with the time unit or inputs (slow {slow:g})")
        assertTrue(rows["waste reversible (kw = 1e-6)"][0] >= f1[0], f"reversible waste lowers the racemic share (slow {slow:g}, k2 {k2_rng[0]:g})")

if failures:
    print("FAILED:", *failures, sep="\n  ")
    raise SystemExit(1)
print("\nselfcheck ok")
