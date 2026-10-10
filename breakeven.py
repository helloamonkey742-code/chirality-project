"""Break-even values: what would each unknown input have to be for physics to win? (2026-10-09)

uncertainty.py draws every uncertain input at once and reports one share per world (Enceladus 22%, Europa 29%,
Early Earth 75%); its low-half / high-half table says which input matters most. This asks the next question. For
each input, cut its range into thirds and into 10 equal-count bins (20 for vertical D_z), let every other input vary
over its full range exactly as before, and find where the share of winning draws crosses 50%: the value that unknown
would have to reach before physics wins in most of the remaining input space. Same draws, seed and rules as
uncertainty.py (nothing there is changed); its headline numbers are asserted before anything is printed.

PREDICTION (written 2026-10-09, before the first run; other inputs drawn over their ranges):
P1. For both moons, vertical D_z is the only input with any bin at or above 50%; every bin of every other input stays below 50%.
P2. 50% crossing in D_z: Enceladus between 1e-5 and 3e-4 m^2/s; Europa between 3e-6 and 3e-5; Earth drops below 50%
    only for D_z below a value between 1e-6 and 1e-5 (below Ledwell's 1.1e-5).
P3. Ranked by top third minus bottom third, D_z is first in all three worlds; next come concentration and time for the
    moons (either order) and time for Earth. D_z top third: Enceladus 45-65%, Europa 70-85%; Earth bottom third 30-50%, top third >= 95%.
P4. Earth: every third of concentration, k2, g and F lies within 75 +/- 2; Earth's lowest time bin (of 10) lies between 40% and 55%.

RESULTS (breakeven.log, 2026-10-09): all four held.
 P1 held: only vertical D_z has bins >= 50% (Enceladus 61-66%, Europa 74-86%); every other input's bins stay at or
    below 36% (Enceladus) and 38% (Europa).
 P2 held: D_z crossing 2.7e-5 (Enceladus), 6.1e-6 (Europa), 2.9e-6 m^2/s (Earth).
 P3 held: D_z first everywhere (thirds differ by +55, +79, +59 points); then concentration (+27.5, +15) and time
    (+20, +14) for the moons and time (+42) for Earth. D_z top third 55% / 79%; Earth's D_z thirds 41% and 99.7%.
 P4 held: Earth's concentration, k2, g and F thirds all 74.5-74.9%; lowest time bin 49%.
 The predictions leaned on uncertainty.py's half-split table, so this mostly confirms the finer picture matches it.
 New: the break-even values, and a ceiling. Even in the top tenth of D_z the moons level off at 66% (Enceladus) and
 86% (Europa): there the reaction does not finish in 17% / 8% of draws and finishes but still loses in 16% / 6%.
 Earth falls below 50% only for D_z below 2.9e-6 (Ledwell measured 1.1e-5, where its bin is 87%) or for times below
 1.4e7 yr (its range starts at 1e7).
"""
import numpy as np

import uncertainty as U
from ocean import YR

UNITS = {"bias g": "kT", "rate k2": "/M/s", "concentration": "M", "horizontal D_h": "m^2/s",
         "vertical D_z": "m^2/s", "time": "yr", "volume": "m^3", "chem. factor F": ""}
LINEAR = ("volume", "chem. factor F")  # drawn uniform, so interpolate linearly; the rest are log-uniform


def binned(v, win, bins):
    """Win share and median input value in each of `bins` equal-count bins of v."""
    edges = np.quantile(v, np.linspace(0, 1, bins + 1))
    b = np.clip(np.searchsorted(edges, v, side="right") - 1, 0, bins - 1)
    return (np.bincount(b, win.astype(float), bins) / np.bincount(b, minlength=bins),
            np.array([np.median(v[b == i]) for i in range(bins)]))


def cross(mid, rate, log=True):
    """First 50% crossing of a win-share curve, interpolated between points: (value, rising, n_crossings) or None."""
    m = np.log10(mid) if log else np.asarray(mid, float)
    i = np.flatnonzero(np.diff(np.sign(rate - 0.5)))
    if not len(i):
        return None
    j = i[0]
    x = m[j] + (0.5 - rate[j]) / (rate[j + 1] - rate[j]) * (m[j + 1] - m[j])
    return (10 ** x if log else x), rate[j + 1] > rate[j], len(i)


def label(k, rate, mid):
    c = cross(mid, rate, k not in LINEAR)
    if c is None:
        return f"{'always' if rate.min() >= 0.5 else 'never'} >=50% (bins {rate.min():.0%}-{rate.max():.0%})"
    x, rising, n = c
    return f">=50% when {k} {'>=' if rising else '<='} {x:.2g} {UNITS[k]}" + (f"  (crosses {n}x)" if n > 1 else "")


def main():
    U.RNG = np.random.default_rng(0)  # the same stream uncertainty.py's main() uses, worlds in the same order
    for name, args in U.WORLDS.items():
        x, win, _ = U.sweep(*args)
        assert abs(win.mean() - U.README_WIN_PCT[name]) < U.WIN_TOL, f"{name}: {win.mean():.1%} != uncertainty.py"
        rows = []
        for k, v in x.items():
            v = v / YR if k == "time" else v
            thirds = binned(v, win, 3)[0]
            rate, mid = binned(v, win, 10)
            fine = binned(v, win, 20) if k == "vertical D_z" else (rate, mid)
            rows.append((thirds[2] - thirds[0], k, thirds, rate, label(k, *fine), fine))
        rows.sort(key=lambda r: -abs(r[0]))
        print(f"\n{name}: physics wins in {win.mean():.1%} of draws (uncertainty.py / README {U.README_WIN_PCT[name]:.0%})")
        print(f"  {'input':15} {'bottom':>6} {'middle':>6} {'top':>6}  (thirds)   10 bins, low to high (%)            where the share reaches 50%")
        for d, k, t, rate, lab, _ in rows:
            print(f"  {k:15} {t[0]:6.1%} {t[1]:6.1%} {t[2]:6.1%}  {' '.join(f'{r * 100:3.0f}' for r in rate)}   {lab}")
        rate, mid = [r for r in rows if r[1] == "vertical D_z"][0][5]
        print("  vertical D_z, 20 equal-count bins (bin median m^2/s: win share)")
        for i in range(0, 20, 5):
            print("    " + "   ".join(f"{m:.1e}: {r:5.1%}" for m, r in zip(mid[i:i + 5], rate[i:i + 5])))
        if name in U.README_TOP_SPLIT:  # moons: reproduce the README's D_z half split, and the closed bottom third
            two = binned(x["vertical D_z"], win, 2)[0]
            lo, hi = U.README_TOP_SPLIT[name]
            assert abs(two[0] - lo) < U.WIN_TOL and abs(two[1] - hi) < U.WIN_TOL, f"{name}: D_z halves {two}"
            assert [r for r in rows if r[1] == "vertical D_z"][0][2][0] == 0.0, f"{name}: bottom D_z third not 0"


if __name__ == "__main__":
    main()
