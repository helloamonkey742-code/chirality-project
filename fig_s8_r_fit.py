"""r from the +1 share across the 7 D/L competition lanes of Essex SI Fig. S8, with no pure-D control lane.

Model: an L at the primer end stops copying, so L-ended primers sit at +1. For the extended primers
(+1 + +2 + >=+3; the unextended primer band is left out)
    g(ee) = (+1 band) / (+1 + +2 + >=+3) = s + (1 - s) pL(ee; r),     pL = r xL / (r xL + xD)
with xL = (1 - ee)/2, xD = (1 + ee)/2, the same p as arm_g_rates.py (fig9_r: model uses xd / (xd + r xl) = 1 - pL), so r
means "L addition rate / D addition rate on a D primer, racemic-pool definition". s = share of D-ended primers that
stall at +1 anyway. Fit (r, s) by unweighted least squares over the 7 lanes (5 dof); s is linear, so it is solved
exactly for each r on a log grid, clipped to [0, 1].
  - 95% interval: profile of RSS(r) with the F criterion RSS <= RSS_min (1 + F(.95; 1, 5) / 5); lane bootstrap as a check.
  - anchored fit: s(r) set so the 90% ee lane is matched exactly, r from the other 6 lanes (5 dof).
  - robustness: the same on band fractions re-measured by fig_s8_bands.measure for each (margin, baseline) variant.
Caveat: if mixed-hand species also hold D-ended primers at +1 (competitive inhibition), s rises with L, the extra spread
in g is credited to r, and the fitted r is an UPPER bound.

Data: fig_s8_bands.log (competition lanes only; *C_D and *C_L control lanes are not read). Variants need the PNG that
fig_s8_bands.py reads (already rendered; the PDF is not re-rendered).
Run:  /opt/anaconda3/bin/python3 -I fig_s8_r_fit.py | tee fig_s8_r_fit.log   (numpy only; PIL via fig_s8_bands for variants)
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
F95 = 6.607890973703367          # scipy.stats.f.ppf(0.95, 1, 5); dof is 5 for every fit below
RG = np.logspace(-2, 2, 4001)    # r grid, 0.01 - 100 (0.23% steps)
EE = np.array([5, 10, 20, 40, 60, 80, 90]) / 100


def pL(ee, r):
    xl, xd = (1 - ee) / 2, (1 + ee) / 2
    return r * xl / (r * xl + xd)


def g_log():
    """(+1)/(+1 + +2 + >=+3) per competition lane, from the printed band table (1 decimal)."""
    g = {}
    for ln in open(os.path.join(HERE, "fig_s8_bands.log")):
        m = re.match(r"\s+(\d+)% ee\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s", ln)
        if m:
            b = [float(x) for x in m.groups()[1:]]
            g[int(m.group(1))] = b[1] / (b[1] + b[2] + b[3])
    assert sorted(g) == [5, 10, 20, 40, 60, 80, 90], g
    return np.array([g[int(round(e * 100))] for e in EE])


def rss_curve(ee, g, rg, anchor=None):
    """RSS and s for every r on the grid. Free: s least squares. anchor=k: s from lane k, RSS over the others."""
    p = pL(ee[None, :], rg[:, None])
    q = 1 - p
    if anchor is None:
        s = np.clip(((g - p) * q).sum(1) / (q * q).sum(1), 0, 1)
        use = np.ones(len(ee), bool)
    else:
        s = np.clip((g[anchor] - p[:, anchor]) / q[:, anchor], 0, 1)
        use = np.arange(len(ee)) != anchor
    res = (g - p - s[:, None] * q)[:, use]
    return (res ** 2).sum(1), s


def fit(ee, g, rg=RG, anchor=None):
    rss, s = rss_curve(ee, g, rg, anchor)
    i = int(np.argmin(rss))
    ok = np.where(rss <= rss[i] * (1 + F95 / 5))[0]
    r1, _ = rss_curve(ee, g, np.array([1.0]), anchor)
    return dict(r=rg[i], s=s[i], rss=rss[i], lo=rg[ok[0]], hi=rg[ok[-1]], open_hi=ok[-1] == len(rg) - 1,
                open_lo=ok[0] == 0, contiguous=len(ok) == ok[-1] - ok[0] + 1, F1=(r1[0] - rss[i]) / (rss[i] / 5),
                resid=g - (s[i] + (1 - s[i]) * pL(ee, rg[i])))


def bootstrap(ee, g, B=2000, seed=0):
    rg = np.logspace(-2, 2, 801)
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < B:
        ix = rng.integers(0, len(ee), len(ee))
        if len(set(ix)) >= 3:
            out.append(rg[np.argmin(rss_curve(ee[ix], g[ix], rg)[0])])
    return np.array(out)


def fmt(f):
    hi = ">100" if f["open_hi"] else f"{f['hi']:.3f}"
    lo = "<0.01" if f["open_lo"] else f"{f['lo']:.3f}"
    return f"r = {f['r']:.3f} (95% {lo} - {hi}), s = {f['s']:.3f}, rss = {f['rss']:.2e}, F(r=1) = {f['F1']:.2f}"


def variants():
    sys.path.insert(0, HERE)
    import fig_s8_bands as fb                      # isolated mode (-I) does not put the script dir on sys.path
    from PIL import Image
    if not os.path.exists(fb.PNG):
        return None
    im = 255.0 - np.asarray(Image.open(fb.PNG).convert("L"), float)
    out = {}
    for m, bs in ((12, "linear"), (24, "linear"), (12, "const"), (24, "const")):
        res = fb.measure(im, m, bs)
        out[(m, bs)] = np.array([res[int(round(e * 100))][0][1:] for e in EE])      # (7, 3): +1, +2, >=+3
    return out


def selfcheck():
    assert np.allclose(pL(EE, 0.7), 1 - ((1 + EE) / 2) / ((1 + EE) / 2 + 0.7 * (1 - EE) / 2))   # arm_g_rates.fig9_r model
    for r0 in (0.5, 1.2):
        for s0 in (0.1, 0.25):
            g = s0 + (1 - s0) * pL(EE, r0)
            f = fit(EE, g)
            assert abs(f["r"] / r0 - 1) < 0.05 and abs(f["s"] - s0) < 0.02, (r0, s0, f["r"], f["s"])
            fa = fit(EE, g, anchor=6)
            assert abs(fa["r"] / r0 - 1) < 0.05, (r0, s0, fa["r"])
    print("selfcheck OK (synthetic r = 0.5, 1.2 recovered within 5%, free and anchored)")


def report(tag, ee, g, boot=False):
    f, fa = fit(ee, g), fit(ee, g, anchor=6)
    print(f"{tag}")
    print("  free (r, s):      " + fmt(f) + ("" if f["contiguous"] else "  [interval not contiguous]"))
    print("  residuals (pp), lanes 5..90 %ee: " + " ".join(f"{100 * x:+5.1f}" for x in f["resid"]))
    print("  s anchored at 90: " + fmt(fa) + ("" if fa["contiguous"] else "  [interval not contiguous]"))
    print("  residuals (pp):                  " + " ".join(f"{100 * x:+5.1f}" for x in fa["resid"]) + "   (90% lane exact)")
    if boot:
        b = bootstrap(ee, g)
        lo, md, hi = np.percentile(b, [2.5, 50, 97.5])
        print(f"  lane bootstrap (2000): r median {md:.3f}, 95% {lo:.3f} - {hi:.3f}; P(r<1) = {np.mean(b < 1):.3f};"
              f" at grid edge: {np.mean(b >= 99):.3f} high, {np.mean(b <= 0.0101):.3f} low")
    return f, fa


def main():
    g = g_log()
    print("Fig. S8 competition lanes: g = +1 / (+1 + +2 + >=+3), from fig_s8_bands.log (pure-D and pure-L lanes not used)")
    print("  %ee      " + " ".join(f"{int(round(100 * e)):6d}" for e in EE))
    print("  g (%)    " + " ".join(f"{100 * x:6.1f}" for x in g))
    print("  pL r=1   " + " ".join(f"{100 * x:6.1f}" for x in pL(EE, 1.0)))
    print("  pL r=.3  " + " ".join(f"{100 * x:6.1f}" for x in pL(EE, 0.3)))
    print("model g = s + (1-s) pL(ee; r), pL = r xL/(r xL + xD), xL=(1-ee)/2, xD=(1+ee)/2 (same p as arm_g_rates.fig9_r, 1 - xd/(xd + r xl))\n")
    f, fa = report("DEFAULT (margin 12 px, linear baseline), band numbers from the log", EE, g, boot=True)

    print("\nr given s (default bands; s held fixed, r by least squares): the r-s trade-off")
    p = pL(EE[None, :], RG[:, None])
    for s0 in (0.0, 0.1, 0.15, 0.2, 0.25, 0.3):
        rs = ((g - s0 - (1 - s0) * p) ** 2).sum(1)
        print(f"  s = {s0:4.2f}: r = {RG[np.argmin(rs)]:.3f}, rms residual {100 * np.sqrt(rs.min() / len(EE)):.1f} pp")

    print("\nsensitivity: bands re-measured by fig_s8_bands.measure()")
    v = variants()
    rows = []
    if v is None:
        print("  PNG missing; variants skipped")
    else:
        d = v[(12, "linear")]
        gd = d[:, 0] / d.sum(1)
        assert np.abs(gd - g).max() < 0.002, (gd, g)                  # module reproduces the logged numbers
        for (m, bs), b in v.items():
            gv = b[:, 0] / b.sum(1)
            fv, fav = report(f"variant margin {m} px, {bs} baseline: g(%) = " + " ".join(f"{100 * x:.1f}" for x in gv), EE, gv)
            rows.append((m, bs, fv, fav))
        ci = lambda x: f"[{x['lo']:.2f}, {'>100' if x['open_hi'] else format(x['hi'], '.2f')}]"
        print("\nper-variant summary: r [95%], s  (free)  |  r [95%], s  (anchored at 90% ee)")
        for m, bs, fv, fav in rows:
            print(f"  margin {m:2d} {bs:6}: r {fv['r']:.2f} {ci(fv)}, s {fv['s']:.2f}  |  r {fav['r']:.2f} {ci(fav)}, s {fav['s']:.2f}")

    print("\nANSWER")
    hi = np.inf if f["open_hi"] else f["hi"]
    print(f"  free fit: r = {f['r']:.2f}, 95% profile {f['lo']:.2f} - {'unbounded' if f['open_hi'] else format(f['hi'], '.2f')}; "
          f"r = 1 is {'rejected' if f['F1'] > F95 else 'not rejected'} (F = {f['F1']:.2f} vs {F95:.2f}); "
          f"r < 1 is {'supported' if hi < 1 else 'not established'} at 95% by this fit.")
    print("  Upper bound only: competitive inhibition by mixed-hand species would raise s with L and inflate r.")
    selfcheck()


if __name__ == "__main__":
    main()
