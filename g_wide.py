"""g_wide.py: widen the bias range's upper end from 2.64e-16 kT to 3.8e-16 and 4.8e-16 kT (Zanasi et al. 1999) and
rerun the headline shares (Table 3, Section 3.3) with everything else as in the originals.

WHY. journals/ija_draft.md E72 flagged that Zanasi, Lazzeretti, Ligabue and Soncini (Advances in BioChirality 1999,
Chapter 29, pp. 377-385) tabulate PVED shifts Epv (their Eq. for Delta Epv = Epv - (-Epv) = 2 Epv, p. 377), and that
with PVED = 2 x shift their extended-basis RPA values for the preferred conformations are 7.2-9.0 aeV (g = 3.0-3.8e-16
at 273 K) and the largest tabulated value (valine, 90 deg) is 11.4 aeV (g = 4.8e-16), above the draft's top of 6.2 aeV
(2.64e-16, alanine, Quack et al. 2022). The author decided to read PVED = 2 x the tabulated shift. Section 1 below
recomputes those numbers from the transcribed tables.

WHAT IS SAMPLED. g is drawn: x["bias g"] is log-uniform on G_OLD = (4e-19, 4e-16) and every headline script maps it to
the same quantile of G_NEW = (4.4e-18, 2.64e-16) by sweep2.remap, times eta (log-uniform 0.1-10, unchanged). So the
shares DO depend on the range, and the top is changed here by setting sweep2.G_NEW at run time; no existing script or
log is edited. The draws, the seeds, the extra-input streams (eta, r, rho: the seedfix.py stream fix, re-installed here
by the same four-line wrapper, so Enceladus is on [s, 0, 1] as in the draft's Table 3), N = 200,000 per world, the
lid test and every rule are the originals; the only thing that changes is S.G_NEW[1]. The code called is
multiwall.run4 (which calls sweep3.outcome3 and sweep2.outcome): the same function that made multiwall.log Section 4.
Two facts follow from the code (sweep2.outcome, sweep3.outcome3, ocean.delta_at), and are checked below:
  * the one-hand mask ("reached") and the lid do not depend on g, so P_none cannot change and P_chance + P_pved is
    constant: every point P_pved gains, P_chance loses;
  * the win test is F * delta >= 2 with delta proportional to g (ocean.delta_at), so each draw flips from "chance"
    to "PVED" once and stays there when g goes up: P_pved can only rise.
NOT RERUN (they also call sweep2.G_NEW but are not this table): rho2side (two-sided mixing), tidal_mix and
tidal_mix_g005, handoff and deephand_scale (the 89.0% early-Earth floor), rerun_base (sensitivities), the figure scripts.

PRE-REGISTERED PREDICTIONS (written and committed before any result-producing run; never edited).
Tops: OLD = 2.64e-16, MID = 3.8e-16, WIDE = 4.8e-16. "Rise" = P_pved at the wider top minus P_pved at OLD, same
column, same lid state, percentage points, seed 0, printed to 0.1. Columns: sweep3 main, sweep3 bound (iii),
multiwall main, cascade x10, sum rule. Old values: multiwall.log:543-557 (Enceladus from seedfix.log).
  P1  Reproduction. At OLD, all 15 worlds x 5 columns x (no lid, final) equal the committed printed values
      (multiwall.log:543-557 for 14 worlds, seedfix.log for Enceladus), digit for digit. (Hard assertion.)
  P2  Conversion. With PVED = 2 x shift, the extended-basis RPA preferred-conformation values (largest negative entry
      per amino acid, Tables 3-7) give: glycine, alanine, valine, aspartate all in 7.1-9.1 aeV and 3.0-3.9e-16 kT at
      273 K; serine lower (4.9-5.1 aeV). The largest magnitude over all tabulated amino-acid entries is valine at
      90 degrees, 11.2-11.5 aeV, g between 4.75e-16 and 4.9e-16. Without the factor 2 the same entries are half as big.
  P3  Invariants. At MID and WIDE, P_none is identical to OLD for every world, column and lid state (to the printed
      digit), and P_chance falls by the same amount (within 0.15, rounding) as P_pved rises, for the multiwall main
      column, no lid.
  P4  Monotone. For every world, column and lid state, P_pved(OLD) <= P_pved(MID) <= P_pved(WIDE).
  P5  Size, multiwall main, no lid, WIDE: every world except early Earth rises by 0.3-3.5 points; the mean over the
      14 other worlds is between 0.8 and 2.5.
  P6  Saturation and lids. Early Earth multiwall main stays 99.9 (rise < 0.1) at WIDE; every final (lid) value for
      Enceladus stays 0.0 at both new tops; sweep3 main for early Earth (the "at least 67%" of the abstract) rises by
      0.3-4.5 points.
  P7  Interpolation. For every world with a WIDE rise of at least 1.0 point (multiwall main, no lid), the MID rise
      is 0.45-0.80 of the WIDE rise (the quantile shift in log g at MID is 0.61 of that at WIDE).
  P8  Ranking. At WIDE the confirmed-ocean order by final multiwall main stays Ganymede > Europa > Enceladus, Ganymede
      and Callisto stay above Pluto (no lid, main), and the lid-free world with the lowest main share is still Titan.
  P9  Seeds. For seeds 1, 2, 3 (multiwall main, no lid) the WIDE-minus-OLD rise of every world is within 0.5 points
      of the seed-0 rise.
  P10 Cost. The whole script, 3 tops with the full five-column run at seed 0 and the main column at seeds 1-3,
      finishes in under 20 minutes at 2 threads under nice.

Run: nice -n 10 env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 /opt/anaconda3/bin/python g_wide.py > g_wide.log 2>&1
(A smoke run with G_WIDE_N=2000 is allowed after the predictions commit, only to catch crashes; its numbers are not used.)

RESULTS (g_wide.log, 2026-10-13): N = 200000 per world, 35 s. 8 of 10 held, 2 failed. The predictions were committed in 9aa7353;
fcecf7b fixed one lookup (the seedfix.log row for P1, code only). The N = 2000 smoke runs were made after the predictions commit and
their numbers were seen (they showed the same two failures); they were not used and no prediction was changed.
  P1  HELD.   150 of 150 values reproduce at 2.64e-16 (g_wide.log:16; run 20-34 against multiwall.log:543-557, Enceladus seedfix.log:125).
  P2  HELD.   Gly/Ala/Val/Asp preferred conformations 7.17-8.99 aeV = 3.05-3.82e-16 kT; Ser 4.96 aeV (2.11e-16); largest entry (valine,
              90 deg) 11.35 aeV = 4.825e-16; without the factor 2 half of each (g_wide.log:5-10).
  P3  HELD.   P_none changes by 0.000 points and P_chance falls by exactly the P_pved rise (g_wide.log:173).
  P4  HELD.   150 (world, column, lid) cells, no violation (g_wide.log:174).
  P5  FAILED. Multiwall main, no lid, 4.8e-16: rises 0.29 (Titan) to 1.17 (Umbriel), mean 0.69 over the 14 worlds; registered 0.3-3.5 and
              mean 0.8-2.5, so the effect is smaller than predicted (Titan misses the lower bound by 0.01; the mean misses by 0.11)
              (g_wide.log:122-136, :175).
  P6  FAILED. Early Earth multiwall main stays 99.9 (+0.00) and Enceladus final stays 0.0, as registered, but sweep3 main for early Earth
              (the abstract's "at least 67%") does not move at all, 67.0 -> 67.0, against a registered rise of 0.3-4.5 (g_wide.log:20, :88, :176).
  P7  HELD.   MID/WIDE rise ratio 0.62-0.63 for the four worlds with a WIDE rise of at least 1.0 (g_wide.log:177).
  P8  HELD.   Final Ganymede 44.5 > Europa 30.4 > Enceladus 0.0; Ganymede 44.5 and Callisto 43.9 above Pluto 41.7; lowest lid-free world
              Titan 36.1 (g_wide.log:178).
  P9  HELD.   Seeds 1-3: largest deviation of the rise from seed 0 is 0.08 (Umbriel) (g_wide.log:155-170, :179).
  P10 HELD.   35 s (g_wide.log:180).
Headline numbers, x10-to-main final range, old top -> 4.8e-16 (g_wide.log:139-153): Ganymede 33.5-44.1 -> 33.7-44.5, Europa 22.0-30.1 ->
22.2-30.4 (no lid 29.0-39.6 -> 29.3-40.1), Enceladus final 0.0 (no lid 19.4-26.4 -> 19.9-27.1), Callisto 33.0-43.5 -> 33.2-43.9, Titan
25.5-35.8 -> 25.6-36.1, Uranian moons 39.7-53.7 -> 40.5-54.8, Early Earth 99.7-99.9 unchanged. No world moves by more than 1.2 points.
"""
import os
import re
import sys
import time

import numpy as np

T0 = time.time()
ORIG_RNG = np.random.default_rng


def fixed_rng(seed=None, *a, **kw):
    """seedfix.py's wrapper: the exact two-element seed [s, 0] becomes [s, 0, 1] (world 0 = Enceladus)."""
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng   # before the project modules draw

import uncertainty as U   # noqa: E402
U.N = int(os.environ.get("G_WIDE_N", U.N))   # smoke test only
import sweep2 as S        # noqa: E402
import multiwall as M     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OLD, MID, WIDE = 2.64e-16, 3.8e-16, 4.8e-16
TOPS = (("OLD 2.64e-16", OLD), ("MID 3.8e-16", MID), ("WIDE 4.8e-16", WIDE))
COLS = ("sweep3 main", "sweep3 bound", "multiwall main", "cascade x10", "sum rule")
KEYS = ("none", "bound", "main", "x10", "sum")
FULL = U.N == 200_000
RES = {}


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def lap(lab):
    print(f"[{lab}: {time.time() - T0:.0f} s elapsed]")


# ---------------------------------------------------------------- 1. the conversion (E72 re-derived from the tables)
HARTREE_AEV = 27.211386245988e18    # 1 hartree in attoelectronvolts
KB, EV = 1.380649e-23, 1.602176634e-19
AEV_KT = 1e-18 * EV / (KB * 273.0)  # 1 aeV in units of kT at 273 K
# Extended-basis RPA column (last column), Epv / 1e-20 hartree, Tables 3-7 of Chapter 29 (pp. 381-383), all angles.
ZANASI = {
    "glycine":   [-12.733, -14.179, -8.692, 0.000, 14.179, 12.733, 7.240, 0.000],
    "alanine":   [-12.318, -13.170, -10.977, -5.238, 2.949, 14.873, 14.602, 10.245, 3.476, -3.490],
    "valine":    [-16.525, -9.107, -0.761, 9.364, 20.859, 10.619, -9.064],
    "serine":    [-5.787, -9.105, 1.172, 12.897, 13.484, 4.576],
    "aspartate": [-13.754, -11.028, 4.709, 16.759, 10.383, -4.749],
}


def to_aev(shift):   # PVED = 2 x shift (Chapter 29, p. 377: Delta Epv = Epv - (-Epv) = 2 Epv)
    return 2 * abs(shift) * 1e-20 * HARTREE_AEV


def section1():
    print("=== 1. Zanasi et al. (1999) extended-basis RPA, PVED = 2 x tabulated shift; g = PVED / kT at 273 K ===")
    print(f"1 aeV = {AEV_KT:.4e} kT at 273 K; 6.2 aeV -> g = {6.2 * AEV_KT:.3e} (the draft's top, 2.64e-16)")
    assert abs(6.2 * AEV_KT / 2.64e-16 - 1) < 2e-3
    pref = {}
    for aa, v in ZANASI.items():
        pref[aa] = min(v)                       # largest negative = the conformation preferred in water (p. 382)
        a, b = to_aev(pref[aa]), to_aev(max(v, key=abs))
        print(f"  {aa:10} preferred (shift {pref[aa]:8.3f}e-20 Ha) {a:6.2f} aeV, g = {a * AEV_KT:.3e};"
              f"  largest entry (shift {max(v, key=abs):8.3f}) {b:6.2f} aeV, g = {b * AEV_KT:.3e};"
              f"  without the factor 2: {a / 2:5.2f} aeV, g = {a / 2 * AEV_KT:.3e}")
    big = max(to_aev(max(v, key=abs)) for v in ZANASI.values())
    p4 = [to_aev(pref[a]) for a in ("glycine", "alanine", "valine", "aspartate")]
    ok = (all(7.1 <= e <= 9.1 for e in p4) and all(3.0 <= e * AEV_KT / 1e-16 <= 3.9 for e in p4)
          and 4.9 <= to_aev(pref["serine"]) <= 5.1 and 11.2 <= big <= 11.5 and 4.75e-16 <= big * AEV_KT <= 4.9e-16)
    check("P2", ok, f"Gly/Ala/Val/Asp preferred {min(p4):.2f}-{max(p4):.2f} aeV = {min(p4) * AEV_KT:.2e}-{max(p4) * AEV_KT:.2e}; "
          f"Ser {to_aev(pref['serine']):.2f} aeV; largest {big:.2f} aeV = {big * AEV_KT:.3e}")
    assert abs(big * AEV_KT / 4.8e-16 - 1) < 0.01   # the WIDE top is Zanasi's largest value, rounded
    assert 3.7e-16 < to_aev(pref["valine"]) * AEV_KT < 3.9e-16   # MID 3.8e-16 is the valine preferred value


# ---------------------------------------------------------------- 2. committed values (P1)
def committed():
    """Printed 'P_pved %, no lid / final' rows: 14 worlds from multiwall.log, Enceladus from seedfix.log."""
    row = re.compile(r"^(\S.*?)\s+" + r"\s*\|\s*".join([r"([\d.]+) /\s+([\d.]+)"] * 5) + r"\s*$")
    ml = open(os.path.join(HERE, "multiwall.log")).read().splitlines()
    i0 = next(i for i, ln in enumerate(ml) if ln.startswith("P_pved %, no lid / final:"))
    out, where = {}, {}
    for j in range(i0 + 1, i0 + 16):
        m = row.match(ml[j])
        assert m, ml[j]
        out[m.group(1).strip()], where[m.group(1).strip()] = tuple(m.groups()[1:]), f"multiwall.log:{j + 1}"
    sl = open(os.path.join(HERE, "seedfix.log")).read().splitlines()
    d0 = next(j for j, ln in enumerate(sl) if ln.startswith("=== D. multiwall Section 4"))
    hits = [(j, row.match(ln.split("new:", 1)[1].strip())) for j, ln in enumerate(sl)
            if j > d0 and ln.strip().startswith("new: Enceladus")]
    hits = [(j, m) for j, m in hits if m]
    assert len(hits) == 1, hits   # seedfix.log:125, the multiwall.log:556 row
    j, m = hits[0]
    out["Enceladus"], where["Enceladus"] = tuple(m.groups()[1:]), f"seedfix.log:{j + 1}"
    assert len(out) == 15
    return out, where


def pv(res, n, k):
    return res[n][k][0][0], res[n][k][1][0]    # P_pved no lid, final


def show(res, title):
    order = sorted(res, key=lambda n: -res[n]["main"][1][0])
    print(f"\n--- {title} ---")
    print("P_pved %, no lid / final: sweep3 main | sweep3 bound (iii) | multiwall main | cascade x10 | sum rule")
    for n in order:
        print(f"{n:18} " + " | ".join(f"{pv(res, n, k)[0]:5.1f} / {pv(res, n, k)[1]:5.1f}" for k in KEYS))
    print("Three-way split, multiwall main, %: P_pved / P_chance / P_none, no lid | final")
    for n in order:
        print(f"{n:18} " + " | ".join(f"{a:5.1f} / {b:4.1f} / {c:5.1f}" for a, b, c in res[n]["main"]))
    return order


def main():
    print(f"g_wide.py, N = {U.N}, seed 0 (+ seeds 1-3 for the main column); Enceladus on stream [s, 0, 1] (seedfix wrapper)")
    print(f"G_OLD = {S.G_OLD}; G_NEW restored to {S.G_NEW} after every run")
    section1()
    G0 = S.G_NEW
    assert G0 == (4.4e-18, OLD)
    runs, seeds = {}, {}
    for lab, top in TOPS:
        S.G_NEW = (G0[0], top)
        assert abs(S.remap(S.G_OLD[1], S.G_OLD, S.G_NEW) / top - 1) < 1e-12   # top quantile maps to the top
        runs[lab] = M.run4(0)
        seeds[lab] = {s: M.run4(s, full=False) for s in (1, 2, 3)}
        lap(f"top {lab}")
    S.G_NEW = G0
    ro, rm, rw = runs[TOPS[0][0]], runs[TOPS[1][0]], runs[TOPS[2][0]]
    names = list(ro)
    assert len(names) == 15 and names[0] == "Enceladus"

    print("\n=== 2. Reproduction at the old top (P1) ===")
    old, where = committed()
    if FULL:
        bad = []
        for n in names:
            got = tuple(f"{v:.1f}" for k in KEYS for v in pv(ro, n, k))
            if got != tuple(old[n]):
                bad.append((n, where[n], got, old[n]))
        assert not bad, f"P1 FAILED (world, log line, got, committed): {bad}"
        check("P1", True, f"15 worlds x 5 columns x (no lid, final) = 150 values equal the committed log values "
              f"(rows {where['Early Earth ocean']} to {where['Dione']} and so on; Enceladus {where['Enceladus']})")
    else:
        print("  (smoke run: P1 not evaluated)")

    show(ro, "OLD top 2.64e-16")
    show(rm, "MID top 3.8e-16")
    show(rw, "WIDE top 4.8e-16")

    print("\n=== 3. Rise in P_pved (points), multiwall main, no lid / final, vs OLD ===")
    print(f"{'world':18} {'OLD':>13} | {'MID rise':>14} | {'WIDE rise':>14} | MID/WIDE (no lid) | x10 rise WIDE (no lid / final)")
    rise = {}
    for n in sorted(names, key=lambda n: -ro[n]["main"][1][0]):
        o, m, w = pv(ro, n, "main"), pv(rm, n, "main"), pv(rw, n, "main")
        rise[n] = (m[0] - o[0], w[0] - o[0], m[1] - o[1], w[1] - o[1])
        ratio = rise[n][0] / rise[n][1] if rise[n][1] > 0 else float("nan")
        x = pv(rw, n, "x10"), pv(ro, n, "x10")
        print(f"{n:18} {o[0]:5.1f} / {o[1]:5.1f} | {rise[n][0]:+5.1f} / {rise[n][2]:+5.1f} | {rise[n][1]:+5.1f} / {rise[n][3]:+5.1f}"
              f" | {ratio:5.2f}             | {x[0][0] - x[1][0]:+5.1f} / {x[0][1] - x[1][1]:+5.1f}")

    print("\n=== 4. Table 3 format: final x10-to-main range, (no-lid range); OLD -> MID -> WIDE ===")
    for n in sorted(names, key=lambda n: -ro[n]["main"][1][0]):
        parts = []
        for r in (ro, rm, rw):
            f, nl = (pv(r, n, "x10")[1], pv(r, n, "main")[1]), (pv(r, n, "x10")[0], pv(r, n, "main")[0])
            parts.append(f"{f[0]:.1f}-{f[1]:.1f}" + (f" ({nl[0]:.1f}-{nl[1]:.1f})" if abs(nl[1] - f[1]) > 1e-9 else ""))
        print(f"{n:18} " + "  ->  ".join(parts))

    print("\n=== 5. Seed spread of the WIDE-minus-OLD rise (multiwall main, no lid) ===")
    sp = {}
    for n in names:
        rs = [rise[n][1]] + [seeds[TOPS[2][0]][s][n]["main"][0][0] - seeds[TOPS[0][0]][s][n]["main"][0][0] for s in (1, 2, 3)]
        sp[n] = max(abs(r - rs[0]) for r in rs[1:])
        print(f"{n:18} seed 0 {rs[0]:+5.2f}  seeds 1-3 " + " ".join(f"{r:+5.2f}" for r in rs[1:]) + f"   max |dev| {sp[n]:.2f}")

    print("\n=== PREDICTIONS ===")
    # P3: P_none identical, P_chance falls by the rise (multiwall main, no lid); also P_none for every column/lid
    dn, dc = 0.0, 0.0
    for r in (rm, rw):
        for n in names:
            for li in (0, 1):
                dn = max(dn, abs(r[n]["main"][li][2] - ro[n]["main"][li][2]))
            dc = max(dc, abs((ro[n]["main"][0][1] - r[n]["main"][0][1]) - (r[n]["main"][0][0] - ro[n]["main"][0][0])))
    check("P3", dn < 0.05 and dc < 0.15, f"max change of P_none {dn:.3f} points; max |chance fall - pved rise| {dc:.3f}")
    # P4: monotone, all columns and lid states
    viol = [(n, k, li) for n in names for k in KEYS for li in (0, 1)
            if not (ro[n][k][li][0] <= rm[n][k][li][0] + 1e-9 <= rw[n][k][li][0] + 2e-9)]
    check("P4", not viol, f"{len(names) * len(KEYS) * 2} (world, column, lid) cells; violations {viol}")
    # P5
    others = [n for n in names if n != "Early Earth ocean"]
    r5 = {n: rise[n][1] for n in others}
    mean5 = float(np.mean(list(r5.values())))
    check("P5", all(0.3 <= v <= 3.5 for v in r5.values()) and 0.8 <= mean5 <= 2.5,
          f"WIDE rises min {min(r5.values()):.2f} ({min(r5, key=r5.get)}) max {max(r5.values()):.2f} ({max(r5, key=r5.get)}), mean {mean5:.2f}")
    # P6
    ee3 = pv(rw, "Early Earth ocean", "none")[0] - pv(ro, "Early Earth ocean", "none")[0]
    enc = (rm["Enceladus"]["main"][1][0], rw["Enceladus"]["main"][1][0], rw["Enceladus"]["x10"][1][0])
    check("P6", rise["Early Earth ocean"][1] < 0.1 and pv(rw, "Early Earth ocean", "main")[0] >= 99.85
          and max(enc) == 0.0 and 0.3 <= ee3 <= 4.5,
          f"Early Earth main rise {rise['Early Earth ocean'][1]:+.2f} (value {pv(rw, 'Early Earth ocean', 'main')[0]:.1f}); "
          f"Enceladus final {enc}; Early Earth sweep3 main {pv(ro, 'Early Earth ocean', 'none')[0]:.1f} -> "
          f"{pv(rw, 'Early Earth ocean', 'none')[0]:.1f} ({ee3:+.2f})")
    # P7
    rat = {n: rise[n][0] / rise[n][1] for n in names if rise[n][1] >= 1.0}
    check("P7", all(0.45 <= v <= 0.80 for v in rat.values()),
          f"{len(rat)} worlds with WIDE rise >= 1.0; MID/WIDE min {min(rat.values()):.2f} max {max(rat.values()):.2f}")
    # P8
    fin = lambda r, n: r[n]["main"][1][0]
    nol = lambda r, n: r[n]["main"][0][0]
    lidfree = [n for n in names if abs(nol(rw, n) - fin(rw, n)) < 1e-9 and n != "Early Earth ocean"]
    low = min(lidfree, key=lambda n: nol(rw, n))
    check("P8", fin(rw, "Ganymede") > fin(rw, "Europa") > fin(rw, "Enceladus") and nol(rw, "Ganymede") > nol(rw, "Pluto")
          and nol(rw, "Callisto") > nol(rw, "Pluto") and low == "Titan",
          f"final G/E/En {fin(rw, 'Ganymede'):.1f}/{fin(rw, 'Europa'):.1f}/{fin(rw, 'Enceladus'):.1f}; Ganymede {nol(rw, 'Ganymede'):.1f}, "
          f"Callisto {nol(rw, 'Callisto'):.1f}, Pluto {nol(rw, 'Pluto'):.1f}; lowest lid-free {low} {nol(rw, low):.1f}")
    # P9
    check("P9", max(sp.values()) <= 0.5, f"max seed deviation of the rise {max(sp.values()):.2f} ({max(sp, key=sp.get)})")
    el = time.time() - T0
    check("P10", el < 1200, f"{el:.0f} s")
    print(f"\n{sum(RES.values())} of {len(RES)} held" + ("" if all(RES.values()) else "; failed: " + ", ".join(k for k, v in RES.items() if not v)))


if __name__ == "__main__":
    main()
