"""Fig. S8 band intensities (Essex, Jia, Zhang & Szostak 2026 SI, p. S14): which band count does Fig. 9 use?

Fig. 9 reports "% reduction in primer extension" (TD2 template, 3'-GG-5' at the reaction centre, so +1 and +2 are
templated) relative to a pure-D control. Two readings:
  (a) extension = products past +1  (>= +2)      reduction = 1 - f_a(lane) / f_a(control)
  (b) extension = products past +2  (>= +3)      reduction = 1 - f_b(lane) / f_b(control)
with f = band signal / total lane signal. Fig. 9a draws "L at +1 stops extension", which is reading (a) in the model
(arm_g_rates.py fig9_r(1)); reading (b) is fig9_r(2).

Inputs: the SI page rendered by  pdftoppm -r 300 -png -f 14 -l 14 papers/essex2026/essex2026_si.pdf <scratch>/p14
Run:    /opt/anaconda3/bin/python3 -I fig_s8_bands.py | tee fig_s8_bands.log   (numpy + PIL only; the Homebrew python3 has no numpy)
"""
import numpy as np
from PIL import Image

PNG = "/private/tmp/claude-501/fig_s8/p14-14.png"
# (%ee D, % reduction, SD) from Fig. 9c; copied from arm_g_rates.py FIG9
FIG9 = [(5, 51.9, 0.6), (10, 50.0, 1.0), (20, 42.1, 0.7), (40, 32.0, 3.0), (60, 20.9, 0.5), (80, 11.9, 0.7), (90, 6.1, 0.4)]
Y0, Y1 = 312, 752                       # rows inside the gel frames (frames at ~302 and ~763)
LEFT = {"primer": (322, 445), "D": (462, 598), "L": (614, 748)}   # x ranges, left gel (3 lanes, no dividers)
TOL = 5.0                               # percentage points: one representative gel vs mean of n = 3
NAMES = ["primer", "+1", "+2", ">=+3"]


def smooth(p, k=7):
    return np.convolve(p, np.ones(k) / k, "same")


def right_lanes(im):
    """x ranges of the 7 competition lanes, from the 8 black divider/frame lines."""
    col = im[330:740, 1000:2300].mean(0)
    dark = np.where(col > col.mean() + 50)[0] + 1000
    groups = np.split(dark, np.where(np.diff(dark) > 3)[0] + 1)
    c = [g.mean() for g in groups]
    return c


def profile(im, x0, x1, base="linear"):
    """Row profile of the inverted lane (bands are peaks), background subtracted."""
    p = im[Y0:Y1, x0:x1].mean(1)
    top, bot = p[3:40].mean(), p[-17:].mean()          # empty gel above the +3 band / below the primer
    b = np.interp(np.arange(len(p)), [21, len(p) - 9], [top, bot]) if base == "linear" else bot
    return p - b


def peaks(p):
    """Rows of the +3, +2, +1, primer peaks, searched in windows anchored on the strongest peak (+2)."""
    s = smooth(p)
    a = lambda lo, hi: lo + int(np.argmax(s[lo:hi]))
    i2 = a(500 - Y0, 580 - Y0)
    i1 = a(i2 + 40, i2 + 91)
    ip = a(i1 + 45, i1 + 101)
    i3 = a(i2 - 100, i2 - 39)
    return i3, i2, i1, ip


def valleys(p, pk):
    s = smooth(p)
    return [a + int(np.argmin(s[a:b + 1])) for a, b in zip(pk[:-1], pk[1:])]      # v23, v12, v01


def bands(p, v):
    """Integrals of primer, +1, +2, >=+3 (everything above the +2/+3 valley, smear included)."""
    v23, v12, v01 = v
    return np.array([p[v01:].sum(), p[v12:v01].sum(), p[v23:v12].sum(), p[:v23].sum()])


def fa_fb(b):
    return (b[2] + b[3]) / b.sum(), b[3] / b.sum()      # past +1, past +2


def measure(im, margin=12, base="linear"):
    """Per-lane band fractions for the left gel (primer, D, L) and the 7 competition lanes."""
    out = {}
    pD = profile(im, LEFT["D"][0] + margin - 12, LEFT["D"][1] - margin + 12, base)
    vD = valleys(pD, peaks(pD))                         # left-gel valleys, reused for the primer and L lanes (no +2 peak)
    sL = smooth(profile(im, LEFT["L"][0] + margin - 12, LEFT["L"][1] - margin + 12, base))
    # the D lane has no primer band to bound +1 from below: use the midpoint of the L lane's two bands (primer, +1)
    vD[2] = (int(np.argmax(sL[580 - Y0:625 - Y0])) + 580 - Y0 + int(np.argmax(sL[645 - Y0:700 - Y0])) + 645 - Y0) // 2
    for k, (x0, x1) in LEFT.items():
        if k == "primer":
            continue                                    # lone primer marker, run offset from the other lanes; checked in selfcheck
        p = profile(im, x0 + margin - 12, x1 - margin + 12, base)
        pk = peaks(p) if k == "D" else None
        out[k] = (bands(p, vD), pk, vD)
    c = right_lanes(im)
    assert len(c) == 8, c
    for (ee, _, _), a, b in zip(FIG9, c[:-1], c[1:]):
        p = profile(im, int(a) + 4 + margin, int(b) - 4 - margin, base)
        pk = peaks(p)
        out[ee] = (bands(p, valleys(p, pk)), pk, valleys(p, pk))
    return out


def reductions(res):
    """% reduction per competition lane for each reading against (i) the left-gel pure-D lane and (ii) the competition
    lanes' own linear extrapolation to 100 %ee (the pure-D limit, within one gel)."""
    ee = np.array([e for e, _, _ in FIG9], float)
    f = np.array([fa_fb(res[int(e)][0]) for e in ee])                   # (7, 2): a, b
    ctrl = np.array(fa_fb(res["D"][0]))
    ext = np.array([np.polyval(np.polyfit(ee, f[:, j], 1), 100.0) for j in (0, 1)])
    r2 = [np.corrcoef(ee, f[:, j])[0, 1] ** 2 for j in (0, 1)]
    return 100 * (1 - f / ctrl), 100 * (1 - f / ext), ctrl, ext, r2


def shape(res):
    """Control-free test: Fig. 9 says f(ee) = f(pure D) * (1 - R(ee)); fit the unknown f(pure D) = c per reading."""
    g = 1 - np.array([r for _, r, _ in FIG9]) / 100
    out = []
    for j in (0, 1):
        f = 100 * np.array([fa_fb(res[e][0])[j] for e, _, _ in FIG9])
        c = (f * g).sum() / (g * g).sum()
        out.append((c, np.sqrt((((f - c * g) / (c * g)) ** 2).mean()), (f[0] / f[-1]) / (g[0] / g[-1])))
    return out                                          # per reading: implied pure-D value, rms relative misfit, ratio test


def selfcheck(im, res):
    # integrator on a synthetic lane with known areas (peaks at the rows the gel uses)
    y = np.arange(Y1 - Y0)
    hs = {"3": 30.0, "2": 100.0, "1": 60.0, "p": 8.0}
    g = lambda c, h: h * np.exp(-0.5 * ((y - c) / 10.0) ** 2)
    syn = g(160, hs["3"]) + g(230, hs["2"]) + g(290, hs["1"]) + g(365, hs["p"])
    pk = peaks(syn)
    assert pk == (160, 230, 290, 365), pk
    got = bands(syn, valleys(syn, pk)) / bands(syn, valleys(syn, pk)).sum()
    want = np.array([hs["p"], hs["1"], hs["2"], hs["3"]]) / sum(hs.values())
    assert np.abs(got - want).max() < 0.01, (got, want)
    # lane structure and assignment sanity
    for k, (b, pk, v) in res.items():
        assert abs((b / b.sum()).sum() - 1) < 1e-9 and b.min() > -0.02 * b.sum(), (k, b)
        if pk:
            assert pk[0] < pk[1] < pk[2] < pk[3], (k, pk)        # +3 above +2 above +1 above primer
    fr = {k: v[0] / v[0].sum() for k, v in res.items()}
    ps = smooth(profile(im, *LEFT["primer"]))                    # primer marker lane: one band, in the primer position
    pr = Y0 + int(np.argmax(ps[620 - Y0:700 - Y0])) + 620 - Y0
    assert abs(pr - 670) < 25 and pr - 600 > 40, pr              # near the L lane's lower band (670), well below +1 (600)
    assert fr["D"][0] < 0.06, fr["D"]                            # pure D: primer consumed in 24 h
    assert fr["L"][2:].sum() < 0.2 and fr["L"][1] > 0.4, fr["L"]   # pure L: mostly primer + one nucleotide, little past +1
    assert all(res[e][0][3] / res[e][0].sum() < res["D"][0][3] / res["D"][0].sum() for e, _, _ in FIG9)  # D lane has most >= +3
    print("selfcheck OK")


def main():
    im = 255.0 - np.asarray(Image.open(PNG).convert("L"), float)
    res = measure(im)
    print("Fig. S8 band fractions (% of lane signal, background subtracted, mean over lane width)")
    print(f"  {'lane':>8} " + "".join(f"{n:>8}" for n in NAMES) + f"{'past +1':>9}{'past +2':>9}   peak rows (+3,+2,+1,primer)")
    for k, (b, pk, _) in res.items():
        f = 100 * b / b.sum()
        rows = "" if not pk else tuple(Y0 + i for i in pk)
        label = {"D": "*C_D", "L": "*C_L"}.get(k, f"{k}% ee")
        print(f"  {label:>8} " + "".join(f"{x:8.1f}" for x in f) + f"{f[2] + f[3]:9.1f}{f[3]:9.1f}   {rows}")

    redD, redX, ctrl, ext, r2 = reductions(res)
    print(f"\ncontrol (left gel, *C_D lane): past +1 {100 * ctrl[0]:.1f}%, past +2 {100 * ctrl[1]:.1f}%")
    print(f"control (competition lanes extrapolated to 100 %ee): past +1 {100 * ext[0]:.1f}%, past +2 {100 * ext[1]:.1f}%"
          f"  (linear R2 vs ee: a {r2[0]:.2f}, b {r2[1]:.2f})")
    ref = np.array([r for _, r, _ in FIG9])
    print("\n% reduction in primer extension (Fig. 9c vs this gel)")
    print(f"  {'%ee':>4} {'Fig.9':>7} | {'a: D-lane':>10} {'b: D-lane':>10} | {'a: extrap':>10} {'b: extrap':>10}")
    for i, (e, r, _) in enumerate(FIG9):
        print(f"  {e:4d} {r:7.1f} | {redD[i, 0]:10.1f} {redD[i, 1]:10.1f} | {redX[i, 0]:10.1f} {redX[i, 1]:10.1f}")
    print(f"\nagreement with Fig. 9c over 7 lanes (rms pp, lanes within {TOL:.0f} pp)")
    for name, red in (("vs left-gel *C_D lane", redD), ("vs extrapolated 100 %ee", redX)):
        for j, rd in enumerate("ab"):
            d = red[:, j] - ref
            print(f"  {name:26} reading ({rd}): rms {np.sqrt((d ** 2).mean()):5.1f}, within tol {int((abs(d) <= TOL).sum())}/7,"
                  f" TD2 5% ee: {red[0, j]:5.1f} (Fig. 9: 51.9)")

    print("\ncontrol-free shape test: fit f(ee) = c * (1 - R_Fig9(ee)); c = implied pure-D value, misfit = rms relative residual,")
    print("  ratio = [f(5%)/f(90%)] / [(1-R5)/(1-R90)] (1.00 if the reading reproduces Fig. 9; control cancels)")
    for rd, (c, mis, rat) in zip("ab", shape(res)):
        meas = 100 * ctrl["ab".index(rd)]
        print(f"  reading ({rd}): implied pure-D {c:5.1f}% (left-gel *C_D lane {meas:.1f}%, extrapolated {100 * ext['ab'.index(rd)]:.1f}%),"
              f" misfit {100 * mis:4.1f}%, ratio {rat:.2f}")

    print("\nsensitivity (lane-width margin px, baseline): 5% ee reduction; 90% ee reduction; control past +1/+2; shape misfit a/b")
    for m, bs in ((12, "linear"), (24, "linear"), (12, "const"), (24, "const")):
        rs = measure(im, m, bs)
        rD, rX, ct, ex, _ = reductions(rs)
        sh = shape(rs)
        print(f"  margin {m:2d}, {bs:6}: 5%  D-lane a {rD[0, 0]:5.1f} b {rD[0, 1]:5.1f} extrap a {rX[0, 0]:5.1f} b {rX[0, 1]:5.1f} |"
              f" 90% D-lane a {rD[-1, 0]:5.1f} b {rD[-1, 1]:5.1f} | ctrl {100 * ct[0]:4.1f}/{100 * ct[1]:4.1f} | misfit a {100 * sh[0][1]:4.1f} b {100 * sh[1][1]:4.1f}")
    selfcheck(im, res)


if __name__ == "__main__":
    main()
