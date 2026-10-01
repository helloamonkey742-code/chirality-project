"""Quench-time sweep (2026-10-01, after J. Dworkin's reply and LPSC 2024 abstract #2645).

R6 used a step rule: a closed parent body is racemic if k2*c*tau < 10 (never finished) or > the fade time
(finished, then faded). But Bennu's fluid left partway through percolation (sulfur-rich solvent front frozen in a ~2 mm
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
Units follow R6: time in reaction times (frank_closed's internal unit), seed 0.1%.
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


def window(t, ee, lo, hi):
    """Decades of k2*c*tau whose frozen |ee| lies in [lo, hi)."""
    lt = np.log10(t)
    m = (np.abs(ee) >= lo) & (np.abs(ee) < hi)
    return np.sum(np.diff(lt)[m[:-1]])


rng = np.random.default_rng(1)
n = 200_000
curves = {}
for slow, label in ((1.0, "fast reverse (1/2000)"), (0.01, "slow reverse (100x)")):
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
        q = frozen(kct, t, ee)
        f = [(q < RAC).mean(), ((q >= RAC) & (q < BAND)).mean(), (q >= BAND).mean()]
        assertTrue(abs(sum(f) - 1) < 1e-12, "bands do not partition the draws")
        print(f"  k2 {k2_rng[0]:g}..1: frozen racemic {f[0]:6.1%} | in meteoritic band {f[1]:6.1%} | above band {f[2]:6.1%}"
              f"   (k2*c*tau spans {np.log10(kct.min()):.1f}..{np.log10(kct.max()):.1f} decades)")

R6 = {(1e-6, 1.0): (.22, .29, .75), (1e-6, 0.01): (.25, .31, .75), (1e-12, 1.0): (.15, .21, .57), (1e-12, 0.01): (.12, .17, .42)}  # review.py R6 / PAPER table
print("\nQ3. P(win | Bennu frozen racemic), frozen-curve rule vs R6 step rule (same draws per row)")
for k2_rng in ((1e-6, 1.0), (1e-12, 1.0)):
    for slow, (t, ee) in curves.items():
        fade = t[np.abs(ee) > 0.1 * np.abs(ee).max()][-1]
        row = []
        for name, args in U.WORLDS.items():
            x, win, _ = U.sweep(*args, k2_rng=k2_rng)
            kct = x["rate k2"] * 10 ** rng.uniform(*np.log10(BENNU_C), U.N) * 10 ** rng.uniform(*np.log10(BENNU_T), U.N) * YR
            rac_q = frozen(kct, t, ee) < RAC
            rac_s = (kct < 10) | (kct > fade)
            pq = (win & rac_q).mean() / max(rac_q.mean(), 1 / U.N)
            ps = (win & rac_s).mean() / max(rac_s.mean(), 1 / U.N)
            row.append(f"{name.split()[0]} {pq:4.0%} (step {ps:4.0%})")
            ref = R6[(k2_rng[0], slow)][len(row) - 1]
            assertTrue(abs(ps - ref) < 0.02, f"step rule {ps:.0%} does not reproduce R6 {ref:.0%} ({name}, k2 {k2_rng[0]:g}, slow {slow:g})")
        print(f"  k2 {k2_rng[0]:g}..1, slow={slow:g}: " + ", ".join(row))

if failures:
    print("FAILED:", *failures, sep="\n  ")
    raise SystemExit(1)
print("\nselfcheck ok")
