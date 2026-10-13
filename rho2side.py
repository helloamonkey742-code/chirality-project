"""rho2side: multiwall Section 4 with the D_z bottom/top ratio rho drawn two-sided, log-uniform 0.01-100 (2026-10-12)

Why: journals/new_analyses_dz_direction.md found no observation or model that fixes the direction of D_z with depth
for any icy moon (convective interior rho ~ 1, a stratified top already handled by the lid test, or boundary-layer
dissipation at the ice, the seafloor or both). Every reversal result drew rho log-uniform on 1-100 (D_z rising with
depth only; sweep2.py step f, sweep3.py:157, multiwall.py:436). This script reruns multiwall Section 4 (main and
cascade time x10) with rho = 10 ** U(-2, 2): equal odds of up or down, the same span.

Model: multiwall.py, sweep3.py, sweep2.py imported, NOT edited. The only change is the rho draw in run4's `ex` dict
(re-implemented here, same RNG and draw order: eta, r, rho from default_rng([seed, world index])); every outcome is
multiwall.multi (rule "max", scale 1 or 10). With the range set back to 10 ** U(0, 2) it is multiwall.run4.

SIGN AND WHAT rho < 1 DOES (checked in the code, and asserted at run time):
  - rho = D_bottom / D_top (sweep2 step f "bottom/top ratio"). multiwall: ln D(r) = ln D_top + 2 int_r^R G, with
    G = ln(rho) / (2H) (sweep3: x = ln(rho) r_w / (4H) = G r_w / 2), so ln D(floor) - ln D_top = 2 G H = ln rho.
    ln rho < 0 means D_z FALLS with depth (strongest under the ice).
  - The wall law used throughout (pin_onset, multiwall) is sharp-interface: upward speed D (2/r)(x - 1), x = G r / 2,
    i.e. D (G - 2/r). The G term moves a wall toward lower D (lower wall tension sigma ~ sqrt(D)); the 2/r term is
    curvature, always inward (down). Eq. 20 (|d ln sigma/dr| < 2/r) is the condition that curvature beats an
    OPPOSING gradient. With rho > 1 the gradient opposes (low D at the top pulls walls up) and above onset walls rise.
    With rho < 1 the low-D side is the floor, so the gradient pulls walls DOWN, the same way as curvature: speed
    D (2/r + |G|) downward. There is no onset, no pinning and no reversal for any |G|; the walls sink faster than in
    uniform D_z and the TOP layer's hand fills (the shallowest hand, which is exactly what step e already gives).
    Applying |d ln sigma/dr| < 2/r with the absolute value would wrongly predict pinning; the signed law does not.
  - CHOICE: rho < 1 draws are given sweep2 step e (the walls sink; the top layer decides; sink time H R/(2 D_z) with
    D_z = the drawn profile mean), i.e. the same convention the rho > 1 draws below onset already use ("a wall's
    sinking speed is not changed by the gradient below onset", "D_z drawn = the profile's mean sets every time").
    That is the symmetric choice. It is what the unmodified code returns (sweep3's pinned mask needs ln rho >= 4H/r_w
    > 0), and this script asserts it draw by draw instead of assuming it.
  - Sensitivity (real profile): for rho < 1 column draws, replace the convention's sink time by the integral of
    dr / (D(r) (2/r + |G|)) from the floor to the ice through the real exponential D(r) (mean = drawn D_z), the same
    span the convention charges. It is expected to be shorter (gradient-assisted; H << R), so the convention is
    conservative.
Seeds 0-3; regression and the step-e identity at seed 0 (and the identity at every seed).

PREDICTIONS (written 2026-10-12, before any result-producing run; "points" are percentage points; "main" = multiwall
rule "max", scale 1; "x10" = cascade time x10; "final" = after the salt lid; "mixture" = 0.5 x sweep2 step e +
0.5 x the one-sided value, both at seed 0):
 (P1) Regression: with rho = 10 ** U(0, 2) this script reproduces multiwall.log Section 4 main three-way split
      (15 worlds, no lid and final) and the x10 P_pved column exactly to the printed 0.1 at seed 0, and equals
      multiwall.run4(0) "main" and "x10" exactly (floats).
 (P2) Fallback identity: in every world, at every seed 0-3, every draw with rho < 1 is unpinned and its main and x10
      masks (P_pved and one-hand, no lid and final) equal sweep2.outcome step e bit for bit.
 (P3) The arithmetic estimates of new_analyses_dz_direction.md hold at seed 0, main and x10 each, rounded to the
      nearest whole percent: Ganymede 32-37 %, Europa (with lid, final) 23-27 %, Enceladus (no lid) 22-26 %.
 (P4) Every world's two-sided P_pved (no lid and final, main and x10) lies within 0.6 points of its mixture.
 (P5) No-lid main P_pved ordering of the five worlds in the dz_direction table is unchanged:
      Early Earth > Ganymede > Europa > Titan > Enceladus; final main: Ganymede > Europa > Enceladus = 0.0.
 (P6) Two-sided main P_chance <= 10 points in every world (no lid and final).
 (P7) Relative to the one-sided run, no-lid main P_pved: Ganymede, Callisto, Umbriel, Oberon and Titania each lose
      5-8 points; Miranda gains (its step e, 30.1, is above its one-sided 27.1); Early Earth loses 10-15 points.
 (P8) Real-profile sink sensitivity: the real sink time is <= the convention's in >= 99 % of rho < 1 column draws in
      every world, and it raises every world's two-sided main P_pved by 0 to 3 points (never lowers it).
 (P9) Seeds 1-3 reproduce every two-sided main share (P_pved / P_chance / P_none, no lid and final) within 0.5 points
      of seed 0.
"""
import os
import sys

import numpy as np

import multiwall as M
import sweep2 as S
import sweep3 as S3
import uncertainty as U
from ocean import delta_at

ONE, TWO = (0.0, 2.0), (-2.0, 2.0)   # log10 rho ranges
GL = np.polynomial.legendre.leggauss(48)
NAMES = ("Early Earth ocean", "Ganymede", "Europa", "Titan", "Enceladus")


def draws(seed, i, lg):
    """multiwall.run4's ex dict, rho range lg (log10); same generator and draw order."""
    rng = np.random.default_rng([seed, i])
    return {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(*lg, U.N)}


def t_sink(R, H, Dz, rho):
    """Sink time of a wall from the ice (R) to the floor (R - H) through D(r) = D_top rho^((R - r)/H), mean = Dz,
    at the signed sharp-interface speed D (2/r - G) downward, G = ln(rho)/(2H) (< 0 for rho < 1). Gauss-Legendre."""
    lr = np.log(rho)
    Dtop = Dz * lr / np.expm1(lr)
    z, w = GL
    r = (R - H / 2)[:, None] + (H / 2)[:, None] * z
    D = Dtop[:, None] * np.exp(lr[:, None] * (R - r) / H[:, None])
    return H / 2 * (w / (D * (2 / r - (lr / (2 * H))[:, None]))).sum(1)


def real_sink(x, L, dz_rng, ex, lid, w_m, r_m):
    """Main masks with the rho < 1 column draws' sink time replaced by t_sink; also the share with t_sink <= conv."""
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    H, R = U.depth(V, L), L / np.pi
    across = (L / U.A) ** 2 / Dh <= t
    column = H * np.sqrt(Dh / Dz) > L
    conv = H * R / (2 * Dz)
    ok = (k2 * c * t >= 10) & ~lid & (k2 * c > 2 * S.K_RAC["Ala"])
    top = np.where(column, np.minimum(1.0, L / (H * np.sqrt(Dh / Dz))), 1.0)
    strong = F * delta_at(c, g, k2, V * ex["r"], t, top) >= 2
    sel = column & (ex["rho"] < 1)
    re = ok & across & (conv <= t)                       # sweep2 step e on these draws (asserted below)
    assert (re[sel] == r_m[sel]).all() and ((re & strong)[sel] == w_m[sel]).all(), "re-implemented step e differs"
    ts = t_sink(R, H[sel], Dz[sel], ex["rho"][sel])
    rr = r_m.copy()
    rr[sel] = (ok & across)[sel] & (ts <= t[sel])
    ww = np.where(sel, rr & strong, w_m)
    return ww, rr, 100 * (ts <= conv[sel]).mean() if sel.any() else np.nan


def run(seed, lg=TWO, full=True):
    _, W = S.chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        ex = draws(seed, i, lg)
        none = np.zeros(U.N, bool)
        lo = ex["rho"] < 1
        o = {}
        for lid_i, l in enumerate((none, salt)):
            we, re = S.outcome(x, a[1], a[4], "e", ex, l)
            for key, kw in (("main", {}), ("x10", {"scale": 10.0})):
                w, r, _ = M.multi(x, a[1], a[4], ex, l, **kw)
                assert (w[lo] == we[lo]).all() and (r[lo] == re[lo]).all(), f"{name}: rho<1 is not step e ({key})"
                o.setdefault(key, []).append(S3.split(w, r))
                if key == "main" and full and lg == TWO:
                    ww, rr, fr = real_sink(x, a[1], a[4], ex, l, w, r)
                    o.setdefault("real", []).append(S3.split(ww, rr))
                    o.setdefault("frac", []).append(fr)
            o.setdefault("e", []).append(S3.split(we, re))
        # rho < 1 never pinned: sweep3's / sweep2's own expression
        Dz = S.remap(x["vertical D_z"], a[4], (S.DZ_FLOOR, a[4][1])) if a[4][0] < S.DZ_FLOOR else x["vertical D_z"]
        H = U.depth(x["volume"], a[1])
        r_w = a[1] / np.pi - np.minimum(H, a[1] * np.sqrt(Dz / x["horizontal D_h"]))
        assert not (lo & (np.log(ex["rho"]) >= 4 * H / r_w)).any() and (r_w > 0).all(), f"{name}: rho<1 pinned"
        o["lo"] = 100 * lo.mean()
        out[name] = o
    return out


def multiwall_log():
    """Section 4 main three-way split (no lid | final) and the x10 P_pved column, as printed in multiwall.log."""
    lines = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "multiwall.log")).read().splitlines()
    s = next(i for i, ln in enumerate(lines) if ln.startswith("=== Section 4"))
    main, x10 = {}, {}
    h1 = next(i for i in range(s, len(lines)) if lines[i].startswith("Three-way split"))
    for ln in lines[h1 + 2: h1 + 17]:
        name, rest = ln[:18].strip(), ln.split("||")[1]
        main[name] = tuple(v.strip() for part in rest.split("|") for v in part.split("/"))
    h2 = next(i for i in range(s, len(lines)) if lines[i].startswith("P_pved %, no lid / final: sweep3 main"))
    for ln in lines[h2 + 1: h2 + 16]:
        name, cols = ln[:18].strip(), ln[18:].split("|")
        x10[name] = tuple(v.strip() for v in cols[3].split("/"))
    assert len(main) == 15 == len(x10), (main, x10)
    return main, x10


def fmt(sp):
    return " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in sp)


def main():
    one = run(0, ONE, full=False)
    if U.N == 200_000:
        main_log, x10_log = multiwall_log()
        bad = [(n, main_log[n], one[n]["main"]) for n in one
               if tuple(f"{v:.1f}" for s in one[n]["main"] for v in s) != main_log[n]]
        bad += [(n, x10_log[n], one[n]["x10"]) for n in one
                if tuple(f"{s[0]:.1f}" for s in one[n]["x10"]) != x10_log[n]]
        ref = M.run4(0)
        bad += [(n, k) for n in one for k in ("main", "x10") if one[n][k] != ref[n][k]]
        print(f"regression, rho = 1-100 vs multiwall.log Section 4 (main three-way split, x10 P_pved; 15 worlds, no lid "
              f"and final, printed 0.1) and vs multiwall.run4(0) floats: {'EXACT' if not bad else 'FAILED ' + str(bad)}")
    two = {0: run(0)}
    for s in (1, 2, 3):
        two[s] = run(s, full=False)
    print("fallback identity: every rho < 1 draw unpinned; its main and x10 masks (P_pved, one hand; no lid and final) "
          "= sweep2 step e bit for bit, all 15 worlds, seeds 0-3: asserted, passed")
    t0 = two[0]
    order = sorted(t0, key=lambda n: -one[n]["main"][1][0])
    print(f"\nN = {U.N}, seed 0. rho < 1 fraction of draws: " + ", ".join(f"{n} {t0[n]['lo']:.1f}" for n in order))
    for key, lab in (("main", "MAIN (rule max)"), ("x10", "CASCADE TIME x10")):
        print(f"\n{lab}. Three-way split, %: P_pved / P_chance / P_none, no lid | final")
        print(f"{'world':18} {'one-sided rho 1-100':>36} || {'two-sided rho 0.01-100':>36}")
        for n in order:
            print(f"{n:18} {fmt(one[n][key])}  ||  {fmt(t0[n][key])}")
    print("\nP_pved %, no lid / final: step e | one-sided main | TWO-SIDED main | mixture main | one-sided x10 | "
          "TWO-SIDED x10 | mixture x10 | two-sided seeds 0-3 min-max main ; x10")
    mix = {}
    for n in order:
        e = [v[0] for v in t0[n]["e"]]
        row = []
        for key in ("main", "x10"):
            o1, o2 = [v[0] for v in one[n][key]], [v[0] for v in t0[n][key]]
            mix[n, key] = [0.5 * (a + b) for a, b in zip(e, o1)]
            row += [o1, o2, mix[n, key]]
        rng_ = {key: [(min(two[s][n][key][l][0] for s in two), max(two[s][n][key][l][0] for s in two)) for l in (0, 1)]
                for key in ("main", "x10")}
        print(f"{n:18} {e[0]:5.1f} / {e[1]:5.1f} | " + " | ".join(f"{p[0]:5.1f} / {p[1]:5.1f}" for p in row)
              + " | " + " ; ".join(" / ".join(f"{a:.1f}-{b:.1f}" for a, b in rng_[k]) for k in ("main", "x10")))
    dev = max(abs(t0[n][k][l][0] - mix[n, k][l]) for n in order for k in ("main", "x10") for l in (0, 1))
    print(f"max |two-sided - mixture| over worlds, no lid and final, main and x10: {dev:.2f} points")
    print("\nSensitivity, rho < 1 column draws sink through the real falling D(r) (gradient-assisted) instead of "
          "H R/(2 D_z):\nworld              two-sided main P_pved no lid / final | real-profile sink | change | "
          "% of rho<1 column draws with t_real <= convention (no lid)")
    for n in order:
        m, r_ = [v[0] for v in t0[n]["main"]], [v[0] for v in t0[n]["real"]]
        print(f"{n:18} {m[0]:5.1f} / {m[1]:5.1f} | {r_[0]:5.1f} / {r_[1]:5.1f} | {r_[0] - m[0]:+5.2f} / "
              f"{r_[1] - m[1]:+5.2f} | {t0[n]['frac'][0]:6.2f}")
    print("\nSeed stability, two-sided main P_pved / P_chance / P_none (no lid and final), max |seed k - seed 0|:")
    for s in (1, 2, 3):
        d = max(abs(two[s][n]["main"][l][j] - t0[n]["main"][l][j]) for n in t0 for l in (0, 1) for j in range(3))
        print(f"  seed {s}: max diff {d:.2f}")
    return one, two


def selfcheck():
    """t_sink: rho -> 1 gives the uniform-D curvature integral (R^2 - (R-H)^2)/(4D); rho = 0.01 matches a dense sum."""
    R, H, Dz = np.array([10.0]), np.array([1.0]), np.array([2e-3])
    assert abs(t_sink(R, H, Dz, np.array([1 - 1e-9]))[0] / ((R ** 2 - (R - H) ** 2) / (4 * Dz))[0] - 1) < 1e-6
    rr = np.linspace(9.0, 10.0, 400001)
    lr = np.log(0.01)
    D = Dz[0] * lr / np.expm1(lr) * np.exp(lr * (10.0 - rr))
    assert abs(t_sink(R, H, Dz, np.array([0.01]))[0] / np.trapezoid(1 / (D * (2 / rr - lr / 2)), rr) - 1) < 1e-6


if __name__ == "__main__":
    U.N = int(sys.argv[1]) if len(sys.argv) > 1 else U.N  # smoke test: python rho2side.py 2000
    selfcheck()
    main()
