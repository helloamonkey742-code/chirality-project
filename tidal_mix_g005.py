"""tidal_mix_g005.py: sensitivity run of tidal_mix.py at mixing efficiency Gamma = 0.05. Pre-registered 2026-10-13.

Why: Munk & Wunsch (1998, p. 2001, Question 8) cite gamma = 0.05 measured in a salt-stabilized fjord (Stigebrandt &
Aure 1989), 4x below tidal_mix.py's mid value 0.2 (range 0.1-0.3). Question: if Gamma = 0.05, how do the interior-band
success shares change, per world, with and without the lid? tidal_mix.py is imported unedited; no existing script, log
or md is edited.

DESIGN (fixed before any result-producing run).
  Band. K_z = Gamma P / (rho A db) is linear in Gamma, so Gamma = 0.05 is applied to the same P range as the original:
  K_low = 0.05 x P_low, K_mid = 0.05 x P_mid, K_high = 0.05 x P_high (original: 0.1 x P_low, 0.2 x P_mid, 0.3 x
  P_high). So K_low is 1/2, K_mid 1/4 and K_high 1/6 of tidal_mix.log's. Interior db = g x 2e-6 (AT21), P rule, g,
  geometry and the nominal-volume ocean are tidal_mix.py's (tidal_mix.power, .geom, .kz). Pluto and Early Earth ocean
  keep the full prior, as there.
  Shares. tidal_mix.run (multiwall.multi, imported unedited, seedfix world-0 stream), seeds 0-3, N = 200000 (same as the
  original full run; the original took 10 s, so no reduction). Band ends below DZ_FLOOR are raised to it; ends above
  1e-3 are not clipped, as in the original.
  Comparisons printed. (a) "G0.2 logged": tidal_mix.py's own band (Gamma 0.1..0.3), recomputed here and checked digit for
  digit against tidal_mix.log's section 4 rows. (b) "G0.05": this run. (c) Exploratory control, not predicted or
  scored: Gamma = 0.2 x the same P_low..P_high (K 4x the G0.05 band), to separate the Gamma effect from the width of
  the original 0.1-0.3 range.
  Regression. The full-prior run must again reproduce multiwall.log / seedfix.log (tidal_mix.targets).

PREDICTIONS (written 2026-10-13 from hand arithmetic on tidal_mix.log, before any result-producing run; never edited).
Delta = this run's interior-band P_pved minus tidal_mix.log's, in percentage points, seed 0. "No lid" / "final" are
the two columns of tidal_mix.log section 4. Hand K values: Gamma 0.05 K_low/K_mid/K_high m^2/s.
  P1  Regression EXACT for the full prior (15 worlds seed 0, seeds 1-3 P_pved), and the recomputed (a) band shares match
      tidal_mix.log section 4 digit for digit.
  P2  K_mid (0.05 x P_mid) is >= 1e-3 for exactly Europa (4.5e-3), Titan (2.7e-3), Callisto (1.05e-3), Triton (0.18) and
      Mimas (2.1e-2); it is < 1e-3 for Enceladus, Ganymede (5.6e-4, joins the lower group), Dione, Miranda (1.2e-4),
      Oberon (3.3e-5), Titania, Umbriel and Ariel.
  P3  K_low < 1e-3 for Europa (7.1e-4) and Ganymede (2.6e-4) and every world below them; Callisto's K_low (1.05e-3)
      stays >= 1e-3 along with Titan (2.4e-3), Triton (0.08) and Mimas (0.021). Ganymede K_high ~1.05e-3, Miranda
      K_high 4.7e-4, Ariel K_low and K_high both ~3e-8.
  P4  Fraction of the band above 1e-3: Europa 90-95 %, Enceladus 20-25 %, Ganymede < 10 %, Miranda 0 %.
  P5  The six worlds whose logged bands were at or above 1e-3 (Europa, Ganymede, Titan, Callisto, Triton, Mimas) change
      by less than 3 points in absolute value, no lid and final (the share model is flat well above ~1e-4; logged
      values 79.5, 80.8, 80.8, 79.8, 78.9, 63.4 no lid; Europa final 60.4, Mimas final 20.4).
  P6  Lower-band worlds fall: Titania by 5-20 points (to 53-68 no lid), Oberon by 1-10, Umbriel by at least 4 (to <= 59.7;
      its band becomes 7e-8..2e-7, Ariel's old range, below the ~1e-7 level where walls start to merge in time).
  P7  Ariel (band a point near 3e-8) stays within -8..+2 of 44.3; Dione (1.8e-7..3.7e-7) within -8..+2 of 51.6 no lid;
      Miranda within -10..0 of 68.8 no lid.
  P8  Enceladus no-lid P_pved is 33-41.3 (logged 41.3) and its final share stays exactly 0.0.
  P9  No sourced world gains more than 1.5 points, no lid or final; the largest fall in the no-lid column belongs to
      Titania or Umbriel and exceeds 5 points.
  P10 Seed stability: max |seed k - seed 0| of the G0.05 band P_pved is < 0.5 for every world, no lid and final.

Run: OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 nice /opt/anaconda3/bin/python tidal_mix_g005.py \
     > tidal_mix_g005.log 2>&1   (TIDAL_N=2000 for a crash-only smoke run; checks are made only at N = 200000.)
"""
import os
import re

import numpy as np

import sweep2 as S
import tidal_mix as T
import uncertainty as U

GAM_NEW = 0.05


def bands_at(W, g):
    """{world: (K_low, K_high)} with K = g P / (rho A db) at the P_low..P_high of tidal_mix.power; plus K_mid."""
    out, mid = {}, {}
    for name, (x, win, a, salt) in W.items():
        p = T.power(name)
        if p is None:
            continue
        A, db = T.geom(a)[2], T.C14[name][5] * T.DRHO_INT
        out[name] = (T.kz(p[0], A, db, g), T.kz(p[2], A, db, g))
        mid[name] = T.kz(p[1], A, db, g)
    return out, mid


def old_bands(W):
    out = {}
    for name, (x, win, a, salt) in W.items():
        p = T.power(name)
        if p is None:
            continue
        A, db = T.geom(a)[2], T.C14[name][5] * T.DRHO_INT
        out[name] = (T.kz(p[0], A, db, T.GAM[0]), T.kz(p[2], A, db, T.GAM[2]))
    return out


def logged_rows():
    """tidal_mix.log section 4: {world: [band no lid P, band final P]} from the 'interior band' columns."""
    rows, on = {}, False
    for l in open(os.path.join(T.HERE, "tidal_mix.log")):
        if l.startswith("Three-way split %"):
            on = True
        elif l.startswith("Seed stability"):
            break
        elif on and "||" in l:
            f = re.findall(r"[-+]?\d+\.\d+", l)
            if len(f) == 14:
                rows[l[:18].strip()] = (float(f[6]), float(f[9]))
    return rows


def main():
    _, W = S.chains(0)
    N = U.N
    print(f"tidal_mix_g005.py, N = {N}. Gamma = {GAM_NEW} sensitivity of the tidal_mix.py interior band\n")
    new, newmid = bands_at(W, GAM_NEW)
    old = old_bands(W)
    ctl, _ = bands_at(W, 0.2)

    print("=== 1. Interior K_z (m^2/s), Gamma 0.05 x (P_low, P_mid, P_high); logged = tidal_mix.log (Gamma 0.1/0.2/0.3) ===")
    print(f"{'world':18} {'K_low':>9} {'K_mid':>9} {'K_high':>9} | {'logged lo':>9} {'logged hi':>9} | {'% band >1e-3':>12} "
          f"{'% eff prior in band':>19} | lo<1e-3 mid<1e-3")
    for n, (lo, hi) in new.items():
        above = 100 * max(0.0, np.log10(hi) - max(np.log10(lo), -3)) / np.log10(hi / lo) if hi > lo else (100.0 if lo > 1e-3 else 0.0)
        pr = (S.DZ_FLOOR, W[n][2][4][1])
        print(f"{n:18} {lo:9.2e} {newmid[n]:9.2e} {hi:9.2e} | {old[n][0]:9.2e} {old[n][1]:9.2e} | {above:12.1f} "
              f"{T.inside((lo, hi), pr):19.1f} | {str(lo < 1e-3):>6} {str(newmid[n] < 1e-3):>8}")
    for n in W:
        if n not in new:
            print(f"{n:18} full prior kept: {T.FLAG[n]}")

    print(f"\n=== 2. Shares, tidal_mix.run (multiwall.multi unedited), N = {N}, seeds 0-3, seedfix world-0 stream ===")
    res = {}
    for s in range(4):
        Ws = W if s == 0 else S.chains(s)[1]
        res[s] = {k: T.run(s, Ws, b) for k, b in (("new", new), ("old", old), ("ctl", ctl))}
    r0 = res[0]
    if N == 200_000:
        t0, ts = T.targets()
        bad = [(n, t0[n]) for n in r0["new"] if tuple(f"{v:.1f}" for s in r0["new"][n]["full"] for v in s) != tuple(t0[n])]
        bad += [(s, n) for s in (1, 2, 3) for n in res[s]["new"]
                if f"{res[s]['new'][n]['full'][0][0]:.1f}/{res[s]['new'][n]['full'][1][0]:.1f}" != ts[s][n]]
        print(f"regression, full prior vs multiwall.log / seedfix.log: {'EXACT' if not bad else 'FAILED ' + str(bad)}")
        lg = logged_rows()
        bad = [(n, lg.get(n), r0["old"][n]["band"]) for n in r0["old"]
               if n not in lg or (f"{r0['old'][n]['band'][0][0]:.1f}", f"{r0['old'][n]['band'][1][0]:.1f}")
               != (f"{lg[n][0]:.1f}", f"{lg[n][1]:.1f}")]
        print(f"recomputed Gamma 0.1-0.3 band vs tidal_mix.log section 4 (P_pved no lid | final, {len(lg)} rows): "
              f"{'EXACT' if not bad else 'FAILED ' + str(bad)}")
    else:
        print(f"regression not checked (N = {N}, smoke)")

    print("\nP_pved %, seed 0, no lid | final.  G0.2 logged = band Gamma 0.1-0.3 (tidal_mix.log); G0.05 = this run; "
          "control = Gamma 0.2 x same P range (exploratory).")
    print(f"{'world':18} {'full prior':>13} || {'G0.2 logged':>13} || {'G0.05':>13} || {'d(G0.05-G0.2)':>15} || "
          f"{'control G0.2':>13} || P_chance / P_none at G0.05 (no lid | final)")
    order = sorted(r0["new"], key=lambda n: -r0["old"][n]["full"][0][0])
    for n in order:
        f, o, w, c = (r0[k][n][b] for k, b in (("old", "full"), ("old", "band"), ("new", "band"), ("ctl", "band")))
        print(f"{n:18} {f[0][0]:5.1f} | {f[1][0]:5.1f} || {o[0][0]:5.1f} | {o[1][0]:5.1f} || {w[0][0]:5.1f} | {w[1][0]:5.1f} || "
              f"{w[0][0] - o[0][0]:+6.1f} | {w[1][0] - o[1][0]:+6.1f} || {c[0][0]:5.1f} | {c[1][0]:5.1f} || "
              f"{w[0][1]:4.1f} / {w[0][2]:5.1f} | {w[1][1]:4.1f} / {w[1][2]:5.1f}")
    print("\nSeed stability, G0.05 band P_pved no lid / final, seeds 1-3, and max |seed k - seed 0|:")
    for n in order:
        d = max(abs(res[s]["new"][n]["band"][l][0] - r0["new"][n]["band"][l][0]) for s in (1, 2, 3) for l in (0, 1))
        print(f"{n:18} " + "  ".join(f"{res[s]['new'][n]['band'][0][0]:5.1f}/{res[s]['new'][n]['band'][1][0]:5.1f}"
                                     for s in (1, 2, 3)) + f"   max diff {d:.2f}")


if __name__ == "__main__":
    main()
