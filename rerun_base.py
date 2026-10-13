"""rerun_base.py: every "uniform-D_z base model" result in journals/ija_draft.md H19, plus H20 (Delta -> c Delta) and
H27 (273 K everywhere), recomputed through the current pipeline (2026-10-12).

The pipeline is multiwall.multi (sweep2 step f -> sweep3 recount -> deepest-hand cascade) on the seed-0 draws of
sweep2.chains, with the per-world extras of multiwall.run4 (eta, feed share r, D_z ratio rho). Each result is given
for the main cascade time and for the cascade x10 (the low end). Nothing is edited: model changes are made by
swapping delta_at in sweep2, sweep3 and multiwall for the length of one call (and back), or by passing the k_rac array
into sweep2.outcome. Writes rerun_base.log (line numbers are cited in journals/new_analyses_rerun.md, also written here).

What depth-dependent D_z cannot change (decided from the code, so these keep their value and lose the base label):
 - Europa's stratified share 23.9% (strat_sweep.log:2) is the salinity draw alone; chain B redraws it bit for bit.
 - The merge times H R / (2 D_z) at D_z = 1e-10 are a formula; 1e-10 is now below the D_z floor (1e-9), so the time
   at the floor is printed too.
 - Bennu's frozen-racemic shares (quench.py Q1) use only k2 and the Bennu c, tau draws, no ocean.
 - No lid touches Early Earth, and H27 moves only Earth: the moon oceans sit at the ice's melting point (~273 K).

Model notes:
 - Assumed-rule rows: the layer share multiplier m acts on every share below 1 (sweep2's top layer, sweep3's deep
   layer and the cascade's recursion layer), the one place assumed_rules.py applied it. Bias rows multiply Delta.
 - Lids: a lid only removes draws (w with lid == w without lid & ~lid, asserted), so a lid item's value is the
   no-lid pipeline value with a new mask; the product rule final = no lid x (1 - p) is checked.
 - Magnetite and nudges: earth_breakeven.py's and moon_nudges.py's pushes, put on the pipeline's g (remapped x eta)
   in every Delta call; Earth D_h for the pinning rule is the pipeline Earth draw.
 - H27: Early Earth T ~ uniform 273-343 K (Blake, Chang & Lepland 2010, 26-35 C; Robert & Chaussidon 2006, ~70 C;
   Knauth & Lowe 2003, 55-85 C). k_rac(T) = k_273 exp(Ea/R (1/273 - 1/T)), Ala 93 kJ/mol, Asp 110 kJ/mol;
   g ~ 1/T (a fixed energy over kT).

PREDICTIONS (written before the first result-producing run; never edited):
 P1  Regression: multiwall.log Section 4 main and x10 P_pved (no lid and final) and the main three-way split are
     reproduced exactly (printed 0.1) for all 15 worlds.
 P2  Vertical D_z is the largest |effect| for both moons (main and x10). Main low half: Enceladus <= 10%,
     Europa <= 5%; high half >= 45% for both.
 P3  Early Earth: every input's |effect| <= 1 point (main and x10).
 P4  The horizontal D_h effect is negative for both moons (main and x10).
 P5  Breakeven (main): Enceladus D_z crossing in 1e-6..3e-5 m^2/s, Europa in 3e-7..1e-5. Enceladus: no other input
     has a bin >= 50%. Europa: the highest non-D_z bin <= 60%. Earth: every bin of every input >= 90%.
 P6  Share x0.1 / x0.28 / x2.1 and Delta >= 1 / >= 3 each move a moon <= 3 points, whole ocean <= 6 (main and x10).
     Earth > Europa > Enceladus in every row, Earth >= 99% in every row (main).
 P7  Bias unlimited equals multiwall P_pved + P_chance (Enceladus 37.1, Europa 44.5, +-0.1). Bias x10 recovers
     >= 40% of P_chance, x100 >= 70% (main, both moons).
 P8  The harsh wobble costs the moons <= 8 points and Earth <= 0.5 (main and x10).
 P9  c = 0.88 moves every world <= 1.5 points, c = 0.98 <= 0.3 (main and x10, no lid and final).
 P10 Lid items follow the product rule within 0.6 points, and land within 1 point of (main): Juno Europa 36.0
     (x10 26.4); Kang Enceladus 7-30 g/kg 7.9 (S* 20) / 5.7 (S* 22.5), 17-40 g/kg 22.9 / 19.0; strat alternatives
     Enceladus 2.2, 9.5, 9.9, Europa 36.1, 14.0, 39.6; club-sandwich Ganymede 22.0 (x10 16.8).
 P11 Merge-in-time (creep) share with the D_z floor: Enceladus 37-40%, Europa 38-41%.
 P12 Depth column: the shallow ocean >= the deep one for every moon (main and x10, no lid and final); the largest
     shallow - deep spread <= 20 points (main).
 P13 k2 extension (main): R2 rule P(win | Bennu racemic) Enceladus <= 8%, Europa <= 10%, Earth 10-35%. R6 rule over
     the original k2 range within 3 points of the no-lid main share; R6 extended Enceladus 12-22, Europa 18-32,
     Earth 55-85. The frozen-curve rule within 3 points of the step rule. Q1 frozen-racemic shares reproduce
     96.9 / 73.1 / 64.2 / 43.1 within 2 points.
 P14 Magnetite (main): Earth weak + field L and D each 43-48%, pinned 6-8%; the L < 50% crossing >= 1e-14 (or none);
     Earth field lead-halving f* within x3 of 6.0e-16; many-sources M = 1e3 and field f* within x3 of the old value
     for each of the 4 worlds.
 P15 Temperature, Early Earth, Ala, T uniform 273-343 K (main): the pipeline k_rac cut drops 1-5 points; the k_rac
     cut applied to every draw (g and k_rac both at T) drops 4-10; g alone <= 0.3; the all-draw cut at 273 K changes
     nothing; Asp all-draw drops 8-25. Moons unchanged.
"""
import os
import re
from contextlib import contextmanager

import numpy as np

import sweep2 as S
import sweep3 as S3
import multiwall as MW
import uncertainty as U
import ocean_worlds as OW
from breakeven import binned, cross, LINEAR
from review import frank_closed, BENNU_C, BENNU_T
from ocean import delta_at as ORIG, YR

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "rerun_base.log")
MD = os.path.join(HERE, "journals", "new_analyses_rerun.md")
SEED, N = 0, U.N
SCALES = (("main", 1.0), ("x10", 10.0))
THREE = ("Enceladus", "Europa", "Early Earth ocean")
E = "Early Earth ocean"
MODS = (S, S3, MW)
LINES, MDROWS, CHECKS = [], [], []
FH = None


def say(s=""):
    for part in str(s).split("\n"):
        LINES.append(part)
        FH.write(part + "\n")
    FH.flush()
    print(s)
    return len(LINES)


def md(sec, item, old, new, line, note=""):
    MDROWS.append((sec, item, old, new, f"rerun_base.log:{line}", note))


def check(p, ok, detail):
    CHECKS.append((p, bool(ok)))
    return say(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def pc(m):
    return 100 * np.mean(m)


def logu(rng, lo, hi, n=N):
    return 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)


@contextmanager
def patched(delta=None, krac=None):
    od, oo = [m.delta_at for m in MODS], S.outcome
    try:
        if delta is not None:
            for m in MODS:
                m.delta_at = delta
        if krac is not None:
            S.outcome = lambda *a, **kw: oo(*a, **{**kw, "krac": krac})
        yield
    finally:
        for m, d in zip(MODS, od):
            m.delta_at = d
        S.outcome = oo


def scaled(m=1.0, k=1.0):
    """Layer shares below 1 times m (capped at 1), Delta times k (k may be an array)."""
    def d(c, g, k2, V, tau, f=1.0):
        with np.errstate(all="ignore"):
            f = np.asarray(f, float)
            f = np.where(f < 1, np.minimum(1.0, m * f), f)
            return ORIG(c, g, k2, V, tau, f) * k
    return d


# ---------------------------------------------------------------- draws
A, W = S.chains(SEED)
NAMES = list(W)
EX = {}
for _i, _n in enumerate(NAMES):
    _r = np.random.default_rng([SEED, _i])
    EX[_n] = {"eta": 10 ** _r.uniform(-1, 1, N), "r": 10 ** _r.uniform(-2, 0, N), "rho": 10 ** _r.uniform(0, 2, N)}
Z = np.zeros(N, bool)


def run(n, scale=1.0, lid=None, x=None, **pk):
    xx, _, a, _ = W[n]
    with patched(**pk):
        w, r, _ = MW.multi(xx if x is None else x, a[1], a[4], EX[n], Z if lid is None else lid, scale=scale)
    return w, r


BASE = {}


def short(n):
    return n.split()[0]


# ---------------------------------------------------------------- 0 regression
def mw_log():
    L = open(os.path.join(HERE, "multiwall.log")).read().splitlines()
    i0 = next(i for i, l in enumerate(L) if l.startswith("=== Section 4"))
    three, pp = {}, {}
    for j in range(i0, len(L)):
        l = L[j]
        if l.startswith("Recursion draws"):
            break
        name, nums = l[:18].strip(), re.findall(r"\d+\.\d", l[18:])
        if name not in NAMES:
            continue
        if "  ||  " in l and len(nums) == 12:
            three[name] = (nums[6:9], nums[9:12], j + 1)
        elif len(nums) == 10:
            pp[name] = (nums[4:6], nums[6:8], j + 1)
    assert len(three) == len(pp) == 15, (len(three), len(pp))
    return three, pp


def section0():
    say(f"=== 0. Regression vs multiwall.log Section 4 (seed {SEED}, N = {N}) ===")
    three, pp = mw_log()
    bad = []
    for n in NAMES:
        for lab, sc in SCALES:
            w0, r0 = run(n, sc)
            w1, r1 = run(n, sc, lid=W[n][3])
            assert (w1 == (w0 & ~W[n][3])).all() and (r1 == (r0 & ~W[n][3])).all(), f"{n}: lid is not a pure mask"
            BASE[n, lab] = (w0, r0, w1, r1)
        w0, r0, w1, r1 = BASE[n, "main"]
        got3 = ([f"{v:.1f}" for v in S3.split(w0, r0)], [f"{v:.1f}" for v in S3.split(w1, r1)])
        gotm = [f"{pc(w0):.1f}", f"{pc(w1):.1f}"]
        got10 = [f"{pc(BASE[n, 'x10'][0]):.1f}", f"{pc(BASE[n, 'x10'][2]):.1f}"]
        if got3 != (three[n][0], three[n][1]) or gotm != pp[n][0] or got10 != pp[n][1]:
            bad.append((n, got3, gotm, got10))
    say("world               P_pved main no lid / final (multiwall.log) | x10 no lid / final | three-way main no lid")
    for n in NAMES:
        w0, r0, w1, r1 = BASE[n, "main"]
        s = S3.split(w0, r0)
        ln = say(f"{n:18} {pc(w0):5.1f} / {pc(w1):5.1f} ({'/'.join(pp[n][0])}, multiwall.log:{pp[n][2]}) | "
                 f"{pc(BASE[n, 'x10'][0]):5.1f} / {pc(BASE[n, 'x10'][2]):5.1f} ({'/'.join(pp[n][1])}) | "
                 f"{s[0]:.1f} / {s[1]:.1f} / {s[2]:.1f}")
        BASE[n, "line"] = ln
    check("P1", not bad, "exact for all 15 worlds" if not bad else f"mismatches {bad}")
    salt = W["Europa"][3]
    ln = say(f"Europa stratified (salinity lid) share in chain B: {pc(salt):.1f}% (strat_sweep.log:2: 23.9%); "
             f"Enceladus {pc(W['Enceladus'][3]):.1f}%")
    md("Unchanged", "Europa stratified share", "23.9% (strat_sweep.log:2)", f"{pc(salt):.1f}%", ln,
       "salinity draw only; base label can be dropped")


# ---------------------------------------------------------------- 1 sensitivities (Figure 5)
def inputs(n):
    v = dict(W[n][0])
    v["eta (g factor)"], v["feed share r"], v["D_z ratio rho"] = EX[n]["eta"], EX[n]["r"], EX[n]["rho"]
    return v


def old_sens():
    HEAD = re.compile(r"^(.+): physics wins in (\d+)% of plausible input space$")
    ROW = re.compile(r"^  (.{15}) low half +(\d+)%  high half +(\d+)%   \(effect ([+-]\d+)%\)$")
    out, cur = {}, None
    for i, raw in enumerate(open(os.path.join(HERE, "uncertainty.log")).read().splitlines()):
        m = HEAD.match(raw)
        if m:
            cur = m.group(1)
            out[cur] = {"base": int(m.group(2))}
            continue
        m = ROW.match(raw)
        if m and cur:
            out[cur][m.group(1).strip()] = (int(m.group(2)), int(m.group(3)), i + 1)
    return out


SENS = {}


def section1():
    say("\n=== 1. One-at-a-time sensitivities (Section 3.3, Figure 5): median split of each input, P_pved, no lid ===")
    old = old_sens()
    for n in THREE:
        for lab, _ in SCALES:
            w = BASE[n, lab][0]
            eff = []
            for k, v in inputs(n).items():
                lo = v <= np.median(v)
                plo, phi = pc(w[lo]), pc(w[~lo])
                eff.append((phi - plo, k, plo, phi))
            eff.sort(key=lambda e: -abs(e[0]))
            say(f"{n} ({lab}): PVED-set hand in {pc(w):.1f}% of draws, no lid")
            for d, k, plo, phi in eff:
                ln = say(f"  {k:15} low half {plo:5.1f}%  high half {phi:5.1f}%   (effect {d:+5.1f})")
                SENS[n, lab, k] = (plo, phi, d, ln)
                if lab == "main" and k in old[n]:
                    o = old[n][k]
                    md("Sensitivities (3.3, Fig. 5)", f"{short(n)} {k}", f"{o[0]}->{o[1]}% (uncertainty.log:{o[2]})",
                       f"{plo:.1f}->{phi:.1f} (x10 {SENS.get((n, 'x10', k), (np.nan,) * 2)[0]:.1f}->"
                       f"{SENS.get((n, 'x10', k), (np.nan,) * 2)[1]:.1f})", ln)
            SENS[n, lab, "_top"] = eff[0][1]
            SENS[n, lab, "_eff"] = {k: d for d, k, _, _ in eff}
    # fill x10 into the md rows written before x10 existed
    for i, row in enumerate(MDROWS):
        if row[0].startswith("Sensitivities"):
            n = [m for m in THREE if row[1].startswith(short(m))][0]
            k = row[1][len(short(n)) + 1:]
            plo, phi, _, ln = SENS[n, "main", k]
            x = SENS[n, "x10", k]
            MDROWS[i] = (row[0], row[1], row[2], f"{plo:.1f}->{phi:.1f} (x10 {x[0]:.1f}->{x[1]:.1f})", row[4], row[5])
    ok = all(SENS[n, lab, "_top"] == "vertical D_z" for n in THREE[:2] for lab, _ in SCALES)
    ok &= SENS["Enceladus", "main", "vertical D_z"][0] <= 10 and SENS["Europa", "main", "vertical D_z"][0] <= 5
    ok &= all(SENS[n, "main", "vertical D_z"][1] >= 45 for n in THREE[:2])
    check("P2", ok, "D_z top: " + ", ".join(f"{short(n)} {lab} {SENS[n, lab, '_top']}" for n in THREE[:2] for lab, _ in SCALES)
          + "; main D_z halves " + ", ".join(f"{short(n)} {SENS[n, 'main', 'vertical D_z'][0]:.1f}/"
                                              f"{SENS[n, 'main', 'vertical D_z'][1]:.1f}" for n in THREE[:2]))
    mx = max(abs(d) for lab, _ in SCALES for d in SENS[E, lab, "_eff"].values())
    check("P3", mx <= 1, f"Earth max |effect| {mx:.2f} points")
    dh = [SENS[n, lab, "horizontal D_h"][2] for n in THREE[:2] for lab, _ in SCALES]
    check("P4", all(d < 0 for d in dh), "D_h effects " + ", ".join(f"{d:+.1f}" for d in dh))


# ---------------------------------------------------------------- 2 breakeven
def bvars(n):
    x, a = W[n][0], W[n][2]
    v = dict(x)
    v["bias g"] = S.remap(x["bias g"], S.G_OLD, S.G_NEW)
    if a[4][0] < S.DZ_FLOOR:
        v["vertical D_z"] = S.remap(x["vertical D_z"], a[4], (S.DZ_FLOOR, a[4][1]))
    v["time"] = x["time"] / YR
    v["g x eta"] = v["bias g"] * EX[n]["eta"]
    v["eta"], v["feed share r"], v["D_z ratio rho"] = EX[n]["eta"], EX[n]["r"], EX[n]["rho"]
    return v


BE = {}


def section2():
    say("\n=== 2. Breakeven ('What would the unknowns'): P_pved by equal-count bins of each pipeline input, no lid ===")
    say("    (bias g = remapped draw; D_z = floored; time in yr; 10 bins, D_z crossing from 20 bins)")
    for n in THREE:
        for lab, _ in SCALES:
            w = BASE[n, lab][0]
            say(f"{n} ({lab}): P_pved {pc(w):.1f}%")
            res = {}
            for k, v in bvars(n).items():
                rate, mid = binned(v, w, 10)
                fine = binned(v, w, 20) if k == "vertical D_z" else (rate, mid)
                c = cross(fine[1], fine[0], k not in LINEAR)
                lab_c = (f"{'always' if fine[0].min() >= 0.5 else 'never'} >=50% (bins {100 * fine[0].min():.0f}-"
                         f"{100 * fine[0].max():.0f}%)") if c is None else \
                    f">=50% when {k} {'>=' if c[1] else '<='} {c[0]:.2g}" + (f" (crosses {c[2]}x)" if c[2] > 1 else "")
                ln = say(f"  {k:15} " + " ".join(f"{100 * r:3.0f}" for r in rate) + f"   {lab_c}")
                res[k] = (rate, mid, fine, c, ln)
            rate20, mid20 = res["vertical D_z"][2]
            ln20 = say("  vertical D_z, 20 bins (median: P_pved) " + "  ".join(f"{m:.1e}: {100 * r:4.1f}%"
                                                                            for m, r in zip(mid20, rate20)))
            res["_ln20"] = ln20
            BE[n, lab] = res
    for n, oldc, ln_old in (("Enceladus", "2.7e-5", 4), ("Europa", "6.1e-6", 20), ("Early Earth ocean", "2.9e-6", 36)):
        c = {lab: BE[n, lab]["vertical D_z"][3] for lab, _ in SCALES}
        f = lambda cc, r: "none (" + f"bins {100 * r.min():.0f}-{100 * r.max():.0f}%)" if cc is None else f"{cc[0]:.2g}"
        md("Breakeven", f"{short(n)} D_z crossing (m^2/s)", f"{oldc} (breakeven.log:{ln_old})",
           f"{f(c['main'], BE[n, 'main']['vertical D_z'][2][0])} (x10 {f(c['x10'], BE[n, 'x10']['vertical D_z'][2][0])})",
           BE[n, "main"]["vertical D_z"][4])
    for n, oldv, ln_old in (("Enceladus", "66%", 4), ("Europa", "86%", 20)):
        md("Breakeven", f"{short(n)} top D_z decile", f"{oldv} (breakeven.log:{ln_old})",
           f"{100 * BE[n, 'main']['vertical D_z'][0][-1]:.0f}% (x10 {100 * BE[n, 'x10']['vertical D_z'][0][-1]:.0f}%)",
           BE[n, "main"]["vertical D_z"][4])
    mxo = {}
    for n in THREE:
        for lab, _ in SCALES:
            others = [(100 * BE[n, lab][k][0].max(), k) for k in BE[n, lab] if not k.startswith("_") and k != "vertical D_z"]
            mxo[n, lab] = max(others)
    for n, oldv, ln_old in (("Enceladus", "36% (time)", 6), ("Europa", "38% (time)", 22)):
        md("Breakeven", f"{short(n)} highest non-D_z bin", f"{oldv} (breakeven.log:{ln_old})",
           f"{mxo[n, 'main'][0]:.0f}% ({mxo[n, 'main'][1]}); x10 {mxo[n, 'x10'][0]:.0f}% ({mxo[n, 'x10'][1]})",
           BE[n, "main"][mxo[n, "main"][1]][4])
    ct = BE[E, "main"]["time"][3]
    md("Breakeven", "Earth time crossing (yr)", "1.4e7 (breakeven.log:37)",
       "none: " + f"bins {100 * BE[E, 'main']['time'][0].min():.0f}-{100 * BE[E, 'main']['time'][0].max():.0f}%"
       if ct is None else f"{ct[0]:.2g}", BE[E, "main"]["time"][4])
    rate20, mid20 = BE[E, "main"]["vertical D_z"][2]
    j = int(np.argmin(np.abs(np.log10(mid20 / 1.1e-5))))
    md("Breakeven", "Earth D_z bin nearest Ledwell 1.1e-5", "86.9% (breakeven.log:47)",
       f"{100 * rate20[j]:.1f}% (bin {mid20[j]:.1e}; x10 {100 * BE[E, 'x10']['vertical D_z'][2][0][j]:.1f}%)",
       BE[E, "main"]["_ln20"])
    cE, cU = BE["Enceladus", "main"]["vertical D_z"][3], BE["Europa", "main"]["vertical D_z"][3]
    ok = cE is not None and 1e-6 <= cE[0] <= 3e-5 and cU is not None and 3e-7 <= cU[0] <= 1e-5
    ok &= mxo["Enceladus", "main"][0] < 50 and mxo["Europa", "main"][0] <= 60
    eall = min(100 * BE[E, "main"][k][0].min() for k in BE[E, "main"] if not k.startswith("_"))
    eall = min(eall, 100 * rate20.min())
    ok &= eall >= 90
    check("P5", ok, f"crossings Enceladus {cE and cE[0]:.2g}, Europa {cU and cU[0]:.2g}; max non-D_z bin "
                    f"{mxo['Enceladus', 'main'][0]:.0f}% / {mxo['Europa', 'main'][0]:.0f}%; Earth lowest bin {eall:.1f}%")


# ---------------------------------------------------------------- 3 assumed rules and bias
ROWS = [("as published", 1.0, 1.0, (21.8, 28.6, 74.7), 2), ("share x0.1", 0.1, 1.0, (19.7, 27.3, 74.7), 3),
        ("share x0.28", 0.28, 1.0, (20.7, 27.9, 74.7), 4), ("share x2.1", 2.1, 1.0, (22.4, 29.0, 74.7), 5),
        ("whole ocean", np.inf, 1.0, (24.0, 29.6, 74.7), 6), ("Delta >= 1", 1.0, 2.0, (23.1, 29.5, 74.7), 7),
        ("Delta >= 3", 1.0, 2 / 3, (21.0, 28.1, 74.7), 8), ("bias x10", 1.0, 10.0, (25.7, 30.7, 74.7), 9),
        ("bias x100", 1.0, 100.0, (27.8, 31.2, 74.7), 10), ("bias unlimited", 1.0, np.inf, (28.6, 31.2, 74.7), 11)]


def section3():
    say("\n=== 3. Assumed rules and bias size (assumed_rules.py rows), P_pved %, no lid: main / x10 ===")
    res = {}
    for lab, m, k, old, ol in ROWS:
        for n in THREE:
            for s, sc in SCALES:
                res[lab, n, s] = pc(run(n, sc, delta=scaled(m, k))[0])
        ln = say(f"  {lab:15} " + " | ".join(f"{short(n)} {res[lab, n, 'main']:5.1f} / {res[lab, n, 'x10']:5.1f}"
                                            for n in THREE))
        md("Assumed rules and bias (3.3)", lab, " / ".join(f"{v}" for v in old) + f" (assumed_rules.log:{ol})",
           " / ".join(f"{res[lab, n, 'main']:.1f}" for n in THREE) + "  (x10 "
           + " / ".join(f"{res[lab, n, 'x10']:.1f}" for n in THREE) + ")", ln)
    for n in THREE:
        for s, _ in SCALES:
            assert abs(res["as published", n, s] - pc(BASE[n, s][0])) < 1e-9, "scaled(1, 1) is not the pipeline"
    rule = [r[0] for r in ROWS[1:4] + ROWS[5:7]]
    mv = max(abs(res[l, n, s] - res["as published", n, s]) for l in rule for n in THREE[:2] for s, _ in SCALES)
    mw = max(abs(res["whole ocean", n, s] - res["as published", n, s]) for n in THREE[:2] for s, _ in SCALES)
    order = all(res[l, E, "main"] > res[l, "Europa", "main"] > res[l, "Enceladus", "main"] and res[l, E, "main"] >= 99
                for l, *_ in ROWS)
    check("P6", mv <= 3 and mw <= 6 and order, f"largest moon move {mv:.1f} (rules), {mw:.1f} (whole ocean); "
                                               f"order Earth > Europa > Enceladus and Earth >= 99 in every row: {order}")
    ok, det = True, []
    for n, ref in (("Enceladus", 37.1), ("Europa", 44.5)):
        w, r = BASE[n, "main"][:2]
        unl = res["bias unlimited", n, "main"]
        assert abs(unl - pc(r)) < 1e-9, "unlimited is not the reached share"
        chance = pc(r) - pc(w)
        g10 = (res["bias x10", n, "main"] - pc(w)) / chance
        g100 = (res["bias x100", n, "main"] - pc(w)) / chance
        ok &= abs(unl - ref) <= 0.1 and g10 >= 0.4 and g100 >= 0.7
        det.append(f"{short(n)} unlimited {unl:.1f} (ref {ref}), x10 recovers {100 * g10:.0f}%, x100 {100 * g100:.0f}% of P_chance")
    check("P7", ok, "; ".join(det))


# ---------------------------------------------------------------- 4 wobble
WOB = [("no wobble", (1e-12, 1e-12 * 1.0001), (1, 1.0001), (22, 29, 75), 1),
       ("mild: 0.1-1%, 1e5-1e8 yr", (1e-3, 1e-2), (1e5, 1e8), (21, 28, 75), 2),
       ("plausible: 0.1-30%, 1e3-1e8 yr", (1e-3, 0.3), (1e3, 1e8), (20, 26, 75), 3),
       ("harsh: 10-30%, 1e3-1e5 yr", (0.1, 0.3), (1e3, 1e5), (17, 24, 75), 4)]


def wobble(d_rng, p_rng):
    def d(c, g, k2, V, tau, f=1.0):
        rng = np.random.default_rng(1)  # same draws on every call, as wobble_sweep.py
        dd = 10 ** rng.uniform(*np.log10(d_rng), len(tau))
        P = 10 ** rng.uniform(*np.log10(p_rng), len(tau)) * YR
        return ORIG(c, g, k2, V, tau, f) * (1 + dd * 2 * np.pi * tau / P) ** -0.25
    return d


def section4():
    say("\n=== 4. Drive swings (wobble_sweep.py), P_pved %, no lid: main / x10 ===")
    res = {}
    for lab, d_rng, p_rng, old, ol in WOB:
        for n in THREE:
            for s, sc in SCALES:
                res[lab, n, s] = pc(run(n, sc, delta=wobble(d_rng, p_rng))[0])
        ln = say(f"  {lab:32} " + " | ".join(f"{short(n)} {res[lab, n, 'main']:5.1f} / {res[lab, n, 'x10']:5.1f}"
                                            for n in THREE))
        md("Drive swings (3.3)", lab, " / ".join(f"{v}" for v in old) + f" (wobble_sweep.log:{ol})",
           " / ".join(f"{res[lab, n, 'main']:.1f}" for n in THREE) + "  (x10 "
           + " / ".join(f"{res[lab, n, 'x10']:.1f}" for n in THREE) + ")", ln)
    drop = {(n, s): res[WOB[0][0], n, s] - res[WOB[3][0], n, s] for n in THREE for s, _ in SCALES}
    ok = all(drop[n, s] <= 8 for n in THREE[:2] for s, _ in SCALES) and all(drop[E, s] <= 0.5 for s, _ in SCALES)
    check("P8", ok, "harsh drop " + ", ".join(f"{short(n)} {s} {v:.1f}" for (n, s), v in drop.items()))


# ---------------------------------------------------------------- 5 H20: Delta -> c Delta
def section5():
    say("\n=== 5. H20: Delta -> c Delta (kz_patches.py c = 0.88 1D, 0.98 2D/3D), P_pved %, no lid / final ===")
    say(f"{'world':18} base main | c=0.88 main | c=0.98 main || base x10 | c=0.88 x10 | c=0.98 x10")
    mx = {0.88: 0.0, 0.98: 0.0}
    for n in NAMES:
        cells = {}
        for s, sc in SCALES:
            b = BASE[n, s]
            cells[s, 1.0] = (pc(b[0]), pc(b[2]))
            for c in (0.88, 0.98):
                cells[s, c] = (pc(run(n, sc, delta=scaled(1.0, c))[0]), pc(run(n, sc, lid=W[n][3], delta=scaled(1.0, c))[0]))
                mx[c] = max(mx[c], *np.abs(np.subtract(cells[s, c], cells[s, 1.0])))
        ln = say(f"{n:18} " + " || ".join(" | ".join(f"{cells[s, c][0]:5.1f}/{cells[s, c][1]:5.1f}" for c in (1.0, 0.88, 0.98))
                                          for s, _ in SCALES))
        md("H20 (Delta -> c Delta)", short(n), f"{cells['main', 1.0][0]:.1f} / {cells['main', 1.0][1]:.1f} (no lid / final)",
           f"c=0.88 {cells['main', 0.88][0]:.1f} / {cells['main', 0.88][1]:.1f}; c=0.98 {cells['main', 0.98][0]:.1f} / "
           f"{cells['main', 0.98][1]:.1f}; x10 c=0.88 {cells['x10', 0.88][0]:.1f} / {cells['x10', 0.88][1]:.1f}", ln)
    check("P9", mx[0.88] <= 1.5 and mx[0.98] <= 0.3, f"largest change c=0.88 {mx[0.88]:.2f}, c=0.98 {mx[0.98]:.2f} points")


# ---------------------------------------------------------------- 6 lids
SSTAR_P = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 27, 28, 28.5]   # robust_ranking.py (TEOS-10)
SSTAR_S = [23.85, 22.1, 20.35, 18.6, 16.85, 15.1, 13.36, 11.64, 9.93, 8.24, 6.57, 4.94, 3.34, 1.77, 1.0, 0.24, 0]


def lid_row(n, pinned, label, old, src, expect):
    out = {}
    for s, _ in SCALES:
        w = BASE[n, s][0]
        out[s] = (pc(w), pc(w & ~pinned), pc(w) * (1 - pinned.mean()))
    ln = say(f"  {label:44} pinned {pc(pinned):5.1f}% | main no lid {out['main'][0]:5.1f} final {out['main'][1]:5.1f} "
             f"(product {out['main'][2]:5.1f}) | x10 final {out['x10'][1]:5.1f} (product {out['x10'][2]:5.1f})")
    md("Lids", label, f"{old} ({src})", f"{out['main'][1]:.1f}% (x10 {out['x10'][1]:.1f}%), pinned {pc(pinned):.1f}%", ln)
    return out, expect


def section6():
    say("\n=== 6. Lid items (robust_ranking, strat_sweep, jupiter_moons): no-lid pipeline draws x a fresh lid mask ===")
    rows = []
    rng = np.random.default_rng([SEED, 1, 3])
    h = np.clip(rng.normal(29e3, 10e3, N), 3e3, 60e3)
    s_star = np.interp(920 * 1.315 * h / 1e6, SSTAR_P, SSTAR_S, right=0)
    pinned = logu(rng, 1, 100) < s_star
    say(f"  Juno ice 29+/-10 km: ice < 23 km in {pc(h < 23e3):.1f}% (robust_ranking.log:1: 27%)")
    rows.append(lid_row("Europa", pinned, "Juno Europa (ice 29+/-10 km)", "26.1%, pinned 9.1%", "robust_ranking.log:1",
                        (36.0, 26.4)))
    rng = np.random.default_rng([SEED, 0, 3])
    kang = {((7, 30), 20.0): 7.9, ((7, 30), 22.5): 5.7, ((17, 40), 20.0): 22.9, ((17, 40), 22.5): 19.0}
    old = {((5, 20), 20.0): ("0.0%, pinned 100%", 3), ((5, 20), 22.5): ("0.0%, pinned 100%", 4),
           ((7, 30), 20.0): ("6.0%, pinned 72%", 5), ((7, 30), 22.5): ("4.2%, pinned 80%", 6),
           ((17, 40), 20.0): ("17.7%, pinned 19%", 7), ((17, 40), 22.5): ("14.7%, pinned 33%", 8)}
    for (sr, ss), (o, ol) in old.items():
        p = logu(rng, *sr) < ss
        rows.append(lid_row("Enceladus", p, f"Kang Enceladus S {sr[0]}-{sr[1]} g/kg, S* {ss}", o,
                            f"robust_ranking.log:{ol}", (kang.get((sr, ss)), None)))
    rng = np.random.default_rng([SEED, 0, 4])
    for (sr, ss), o, ol, ex in ((((5, 20), 18.0), "1.6%, p 92.4%", 5, 2.2), (((5, 40), 20.0), "7.4%, p 66.6%", 6, 9.5),
                                (((1, 100), 20.0), "7.7%, p 65.0%", 7, 9.9)):
        rows.append(lid_row("Enceladus", logu(rng, *sr) < ss, f"strat alt Enceladus S {sr[0]}-{sr[1]}, S* {ss}", o,
                            f"strat_sweep.log:{ol}", (ex, None)))
    rng = np.random.default_rng([SEED, 1, 4])
    for (sr, ss), o, ol, ex in ((((1, 100), 1.5), "26.0%, p 8.8%", 8, 36.1), (((1, 100), 19.5), "10.1%, p 64.6%", 9, 14.0),
                                (((3, 100), 3.0), "28.5%, p 0.0%", 10, 39.6)):
        rows.append(lid_row("Europa", logu(rng, *sr) < ss, f"strat alt Europa S {sr[0]}-{sr[1]}, S* {ss}", o,
                            f"strat_sweep.log:{ol}", (ex, None)))
    rng = np.random.default_rng([SEED, 2, 3])
    rows.append(lid_row("Ganymede", rng.random(N) < 0.5, "club-sandwich Ganymede (lid p = 0.5)", "14.2%",
                        "jupiter_moons.log:7", (22.0, 16.8)))
    prod = max(abs(o[s][1] - o[s][2]) for o, _ in rows for s, _ in SCALES)
    far = [(o["main"][1], ex) for o, ex in rows if ex[0] is not None and abs(o["main"][1] - ex[0]) > 1.0]
    far += [(o["x10"][1], ex) for o, ex in rows if ex[1] is not None and abs(o["x10"][1] - ex[1]) > 1.0]
    check("P10", prod <= 0.6 and not far, f"largest |final - product| {prod:.2f}; values > 1 point from expected: {far}")


# ---------------------------------------------------------------- 7 merge
def section7():
    say("\n=== 7. Merge times and merge-in-time share (Section 3.2) ===")
    for n, old_t, ol_t, old_c, ol_c in (("Enceladus", "1.5e12", 34, 33, 3), ("Europa", "2.8e13", 35, 34, 14)):
        V, L = U.WORLDS[n][:2]
        H = U.depth(V, L)
        t10, t9 = (H * (L / np.pi) / (2 * dz) / YR for dz in (1e-10, 1e-9))
        ln = say(f"  {n}: nominal depth {H / 1e3:.1f} km; layers merge within {t10:.1e} yr at D_z = 1e-10, {t9:.1e} yr at "
                 f"the D_z floor 1e-9")
        md("Unchanged", f"{n} merge time at D_z = 1e-10", f"{old_t} yr (uncertainty.log:{ol_t})",
           f"{t10:.1e} yr; {t9:.1e} yr at the floor 1e-9", ln, "formula; 1e-10 is below the pipeline's D_z floor")
        x, a = W[n][0], W[n][2]
        Dh, t, Hx = x["horizontal D_h"], x["time"], U.depth(x["volume"], L)
        across = (L / U.A) ** 2 / Dh <= t
        raw = across & (Hx * (L / np.pi) / (2 * x["vertical D_z"]) <= t)
        fl = across & (Hx * (L / np.pi) / (2 * S.remap(x["vertical D_z"], a[4], (S.DZ_FLOOR, a[4][1]))) <= t)
        ln = say(f"  {n}: merge-in-time share {pc(raw):.1f}% with D_z 1e-10..1e-3 (uncertainty.log:{ol_c}: {old_c}%), "
                 f"{pc(fl):.1f}% with the floor 1e-9..1e-3 (pipeline)")
        md("Merge (3.2)", f"{n} merge-in-time share", f"{old_c}% (uncertainty.log:{ol_c})", f"{pc(fl):.1f}% (raw D_z: {pc(raw):.1f}%)", ln)
        if n == "Enceladus":
            ce = pc(fl)
    check("P11", 37 <= ce <= 40 and 38 <= pc(fl) <= 41, f"floored creep Enceladus {ce:.1f}%, Europa {pc(fl):.1f}%")


# ---------------------------------------------------------------- 8 depth column (Table 3)
def section8():
    say("\n=== 8. Table 3 depth column: chain B's own depth draws (ocean_worlds.py order), P_pved % no lid / final ===")
    olds = {}
    for i, raw in enumerate(open(os.path.join(HERE, "ocean_worlds.log")).read().splitlines()):
        m = re.match(r"^(\S+)\s.*?(\d+) km (\d+)% / (\d+) km (\d+)% \(base", raw)
        if m:
            olds[m.group(1)] = (m.group(3), m.group(5), i + 1)
    U.RNG = np.random.default_rng(SEED)
    salt = lambda s_rng, s_star: (U.logu(*s_rng) < s_star) if s_star else np.zeros(N, bool)
    for name, (s_rng, s_star) in OW.STRAT.items():
        U.sweep(*U.WORLDS[name])
        salt(s_rng, s_star)
    spread, ok = 0.0, True
    for name, (R, H, hs, t, s_rng, s_star, _) in OW.WORLDS.items():
        x0 = U.sweep(*OW.args(R, H, t))[0]
        salt(s_rng, s_star)
        assert all((x0[k] == W[name][0][k]).all() for k in x0), f"{name}: chain B replay drifted"
        cells = {}
        for h in hs:
            ah = OW.args(R, h, t)
            xh = U.sweep(*ah)[0]
            lid_h = salt(s_rng, s_star)
            for s, sc in SCALES:
                w = MW.multi(xh, ah[1], ah[4], EX[name], Z, scale=sc)[0]
                cells[h, s] = (pc(w), pc(w & ~lid_h))
        lo, hi = hs
        for s, _ in SCALES:
            ok &= all(cells[lo, s][j] >= cells[hi, s][j] for j in (0, 1))
        spread = max(spread, cells[lo, "main"][0] - cells[hi, "main"][0], cells[lo, "main"][1] - cells[hi, "main"][1])
        b = BASE[name, "main"]
        ln = say(f"  {name:10} base {H / 1e3:.0f} km {pc(b[0]):5.1f}/{pc(b[2]):5.1f} | " +
                 " | ".join(f"{h / 1e3:.0f} km main {cells[h, 'main'][0]:5.1f}/{cells[h, 'main'][1]:5.1f} "
                            f"x10 {cells[h, 'x10'][0]:5.1f}/{cells[h, 'x10'][1]:5.1f}" for h in hs))
        o = olds[name]
        md("Table 3 depth column", name, f"{lo / 1e3:.0f} km {o[0]}% / {hi / 1e3:.0f} km {o[1]}% final (ocean_worlds.log:{o[2]})",
           f"{lo / 1e3:.0f} km {cells[lo, 'main'][1]:.1f}% / {hi / 1e3:.0f} km {cells[hi, 'main'][1]:.1f}% final "
           f"(x10 {cells[lo, 'x10'][1]:.1f} / {cells[hi, 'x10'][1]:.1f})", ln)
    check("P12", ok and spread <= 20, f"shallow >= deep everywhere: {ok}; largest main spread {spread:.1f} points")


# ---------------------------------------------------------------- 9 k2 extension and Bennu
def section9():
    say("\n=== 9. k2 extension (Section 3.6) and Bennu (review.py R2/R6, quench.py Q1/Q3), P_pved no lid ===")
    curves = {s: frank_closed(200.0, s) for s in (1.0, 0.01)}
    spans = {}
    for s, (t, ee) in curves.items():
        on = t[np.abs(ee) > 0.1 * np.abs(ee).max()]
        spans[s] = (on[0], on[-1])
        say(f"  closed curve slow = {s:g}: excess above 10% of peak for k2 c tau {on[0]:.1g}..{on[-1]:.1g}")
    frozen = lambda kct, t, ee: np.interp(np.log10(kct), np.log10(t), np.abs(ee), left=0.001, right=0.0)
    R = {}
    for n in THREE:
        i = NAMES.index(n)
        x = W[n][0]
        x2 = dict(x)
        x2["rate k2"] = S.remap(x["rate k2"], (1e-6, 1.0), (1e-12, 1.0))
        rng = np.random.default_rng([SEED, i, 2])
        cB, tB = logu(rng, *BENNU_C), logu(rng, *BENNU_T) * YR
        for s, sc in SCALES:
            for kl, xx, w in (("1e-6", x, BASE[n, s][0]), ("1e-12", x2, run(n, sc, x=x2)[0])):
                kct = xx["rate k2"] * cB * tB
                cond = lambda rac: 100 * (w & rac).mean() / max(rac.mean(), 1 / N)
                R[n, s, kl, "win"] = pc(w)
                R[n, s, kl, "R2"] = cond(kct < 10)
                R[n, s, kl, "rac"] = pc(kct < 10)
                for sl, (t, ee) in curves.items():
                    R[n, s, kl, "R6", sl] = cond((kct < spans[sl][0]) | (kct > spans[sl][1]))
                    R[n, s, kl, "Q3", sl] = cond(frozen(kct, t, ee) < 0.052)
                    R[n, kl, "Q1", sl] = pc(frozen(kct, t, ee) < 0.052)
    lines = {}
    for s, _ in SCALES:
        for kl in ("1e-6", "1e-12"):
            for n in THREE:
                lines[n, s, kl] = say(
                    f"  {s:4} k2 {kl:5}..1 {n:18} wins {R[n, s, kl, 'win']:5.1f} | Bennu racemic (R2) {R[n, s, kl, 'rac']:4.1f}% "
                    f"P(win|rac) {R[n, s, kl, 'R2']:5.1f} | R6 step s=1 {R[n, s, kl, 'R6', 1.0]:5.1f} s=0.01 "
                    f"{R[n, s, kl, 'R6', 0.01]:5.1f} | frozen s=1 {R[n, s, kl, 'Q3', 1.0]:5.1f} s=0.01 {R[n, s, kl, 'Q3', 0.01]:5.1f}")
    q1 = {}
    for sl, olds in ((1.0, (96.9, 73.1)), (0.01, (64.2, 43.1))):
        for kl, o in zip(("1e-6", "1e-12"), olds):
            q1[sl, kl] = (R["Enceladus", kl, "Q1", sl], o)
    ln = say("  Q1 frozen-racemic share (Enceladus k2 and Bennu draws; ocean-free): " +
             ", ".join(f"s={sl:g} k2 {kl} {v:.1f}% (quench.log {o})" for (sl, kl), (v, o) in q1.items()))
    md("Unchanged", "Bennu Q1 frozen-racemic shares", "96.9 / 73.1 / 64.2 / 43.1% (quench.log:2-3,6-7)",
       " / ".join(f"{v:.1f}" for v, _ in q1.values()), ln, "k2 and Bennu draws only; base label can be dropped")
    for n, ow, orac in (("Enceladus", "14% wins, P(win|rac) 2%", "15 / 17"), ("Europa", "19%, 4%", "21 / 22"),
                        (E, "57%, 15%", "57 / 55")):
        md("k2 extension and Bennu (3.6)", f"{short(n)} k2 1e-12..1: wins, R2 P(win|Bennu racemic)",
           f"{ow} (review.log:{19 + THREE.index(n)})",
           f"{R[n, 'main', '1e-12', 'win']:.1f}%, {R[n, 'main', '1e-12', 'R2']:.1f}% (x10 {R[n, 'x10', '1e-12', 'win']:.1f}%, "
           f"{R[n, 'x10', '1e-12', 'R2']:.1f}%)", lines[n, "main", "1e-12"])
        md("k2 extension and Bennu (3.6)", f"{short(n)} R6 step rule, k2 1e-12..1 (s=1 / s=0.01)",
           f"{orac} (review.log:72-73)", f"{R[n, 'main', '1e-12', 'R6', 1.0]:.1f} / {R[n, 'main', '1e-12', 'R6', 0.01]:.1f} "
           f"(x10 {R[n, 'x10', '1e-12', 'R6', 1.0]:.1f} / {R[n, 'x10', '1e-12', 'R6', 0.01]:.1f}); frozen rule "
           f"{R[n, 'main', '1e-12', 'Q3', 1.0]:.1f} / {R[n, 'main', '1e-12', 'Q3', 0.01]:.1f}", lines[n, "main", "1e-12"])
    oldr6 = {"Enceladus": "22 / 23", "Europa": "29 / 30", E: "75 / 75"}
    for n in THREE:
        md("k2 extension and Bennu (3.6)", f"{short(n)} R6 step rule, k2 1e-6..1 (s=1 / s=0.01)",
           f"{oldr6[n]} (review.log:70-71)", f"{R[n, 'main', '1e-6', 'R6', 1.0]:.1f} / {R[n, 'main', '1e-6', 'R6', 0.01]:.1f}"
           f" (x10 {R[n, 'x10', '1e-6', 'R6', 1.0]:.1f} / {R[n, 'x10', '1e-6', 'R6', 0.01]:.1f})", lines[n, "main", "1e-6"])
    ok = R["Enceladus", "main", "1e-12", "R2"] <= 8 and R["Europa", "main", "1e-12", "R2"] <= 10 \
        and 10 <= R[E, "main", "1e-12", "R2"] <= 35
    ok &= all(abs(R[n, "main", "1e-6", "R6", sl] - pc(BASE[n, "main"][0])) <= 3 for n in THREE for sl in (1.0, 0.01))
    band = {"Enceladus": (12, 22), "Europa": (18, 32), E: (55, 85)}
    ok &= all(band[n][0] <= R[n, "main", "1e-12", "R6", sl] <= band[n][1] for n in THREE for sl in (1.0, 0.01))
    ok &= all(abs(R[n, "main", kl, "Q3", sl] - R[n, "main", kl, "R6", sl]) <= 3
              for n in THREE for kl in ("1e-6", "1e-12") for sl in (1.0, 0.01))
    ok &= all(abs(v - o) <= 2 for v, o in q1.values())
    check("P13", ok, "R2 ext " + ", ".join(f"{short(n)} {R[n, 'main', '1e-12', 'R2']:.1f}" for n in THREE) +
          "; R6 orig " + ", ".join(f"{R[n, 'main', '1e-6', 'R6', 1.0]:.1f}/{R[n, 'main', '1e-6', 'R6', 0.01]:.1f}" for n in THREE) +
          "; R6 ext " + ", ".join(f"{R[n, 'main', '1e-12', 'R6', 1.0]:.1f}/{R[n, 'main', '1e-12', 'R6', 0.01]:.1f}" for n in THREE) +
          "; |frozen - step| max " + f"{max(abs(R[n, 'main', kl, 'Q3', sl] - R[n, 'main', kl, 'R6', sl]) for n in THREE for kl in ('1e-6', '1e-12') for sl in (1.0, 0.01)):.1f}")


# ---------------------------------------------------------------- 10 magnetite
EGRID = 10 ** np.arange(-20, -1.99, 0.5)     # earth_breakeven.py
NGRID = 10 ** np.arange(-19, -9.99, 0.25)    # moon_nudges.py
KINDS = ("one source", "many sources", "field")
NW = ("Enceladus", "Europa", E, "Ganymede")
OLD_NUDGE = {"Enceladus": (7.3e-17, 1.3e-16, 4.1e-15, 1.2e-15, 19), "Europa": (5.2e-17, 9.4e-17, 3.0e-15, 1.3e-15, 20),
             E: (3.9e-17, 7.4e-17, 2.3e-15, 6.0e-16, 21), "Ganymede": (4.7e-17, 8.7e-17, 2.8e-15, 1.3e-15, 22)}


def earth_field(hand, pinned, f_mag=None):
    xE = W[E][0]
    DH, L = xE["horizontal D_h"], W[E][2][1]

    def biased(c, g, k2, V, tau, f=1.0):
        n = len(tau)
        rng = np.random.default_rng(1)  # earth.py's draw order: P_rev, f_mag, ee, A, s
        P, fm = logu(rng, 1e5, 1e8, n), logu(rng, 1e-8, 1e-2, n)
        if f_mag is not None:
            fm[:] = f_mag
        ee, A_, s = rng.uniform(0.1, 0.6, n), logu(rng, 0.01, 0.5, n), rng.choice([-1, 1], n)
        R = np.minimum(1, 0.53 * np.sqrt(P * YR / np.sqrt(tau / (k2 * c))))
        eps = np.sqrt(DH / (k2 * c)) / (L / np.pi)
        pinned[:] = (eps < 1) & (fm * ee > 0.47 * eps * (A_ + 3 * eps))
        return ORIG(c, 1.0, k2, V, tau, f) * np.where(pinned, 0.0, g + fm * ee * A_ * R * s) * hand
    return biased


def nudge(kind, fv, hand):
    def biased(c, g, k2, V, tau, frac=1.0):
        n = len(tau)
        rng = np.random.default_rng(1)  # identical pushes on every call, as moon_nudges.py
        ee = rng.uniform(0.1, 0.6, n)
        if kind == "one source":
            nu = rng.choice([-1, 1], n)
        elif kind == "many sources":
            nu = rng.standard_normal(n)
        else:
            R = np.minimum(1, 0.53 * np.sqrt(logu(rng, 1e5, 1e8, n) * YR / np.sqrt(tau / (k2 * c))))
            nu = logu(rng, 0.01, 0.5, n) * R * rng.choice([-1, 1], n)
        return ORIG(c, 1.0, k2, V, tau, frac) * (g + fv * ee * nu) * hand
    return biased


def section10():
    say("\n=== 10. Magnetite (Section 4 coin flip): earth_breakeven.py field and moon_nudges.py pushes on the pipeline ===")
    ef = {}
    for s, sc in SCALES:
        for lab, fm in (("as published (f_mag 1e-8..1e-2)", None), ("field share ~0 (f_mag 1e-30)", 1e-30)):
            pinned = np.zeros(N, bool)
            L_, D_ = (pc(run(E, sc, delta=earth_field(h, pinned, fm))[0]) for h in (1, -1))
            ef[s, lab] = (L_, D_, pc(pinned))
            ln = say(f"  Earth {s:4} {lab:32}: L {L_:5.1f}%  D {D_:5.1f}%  no hand {100 - L_ - D_:5.1f}%  pinned {pc(pinned):4.1f}%")
            if fm is None:
                ef[s, "ln"] = ln
    md("Magnetite (4)", "Earth weak (L) + field: L / D / pinned", "35 / 35 / 7% (earth.log:4; 34.7 / 34.6 / 7.0 earth_breakeven.log:1)",
       f"{ef['main', 'as published (f_mag 1e-8..1e-2)'][0]:.1f} / {ef['main', 'as published (f_mag 1e-8..1e-2)'][1]:.1f} / "
       f"{ef['main', 'as published (f_mag 1e-8..1e-2)'][2]:.1f}% (x10 {ef['x10', 'as published (f_mag 1e-8..1e-2)'][0]:.1f} / "
       f"{ef['x10', 'as published (f_mag 1e-8..1e-2)'][1]:.1f})", ef["main", "ln"])
    ecross = {}
    for s, sc in SCALES:
        Ls = []
        for fm in EGRID:
            Ls.append(pc(run(E, sc, delta=earth_field(1, np.zeros(N, bool), fm))[0]) / 100)
        c = cross(EGRID, np.array(Ls))
        ecross[s] = c
        ecross[s, "ln"] = say(f"  Earth {s:4} L share over f_mag 1e-20..1e-2 (half decades): " + " ".join(f"{100 * v:.0f}" for v in Ls)
                              + f" -> L < 50% for f_mag above {'never' if c is None else f'{c[0]:.2e}'}")
    md("Magnetite (4)", "Earth f_mag where L < 50%", "1.88e-15 (earth_breakeven.log:44)",
       f"{'never' if ecross['main'] is None else f'{ecross[chr(109) + chr(97) + chr(105) + chr(110)][0]:.2e}'} "
       f"(x10 {'never' if ecross['x10'] is None else f'{ecross[chr(120) + chr(49) + chr(48)][0]:.2e}'})", ecross["main", "ln"])
    star = {}
    for s, sc in SCALES:
        base = {n: pc(run(n, sc, delta=nudge("one source", 0.0, 1))[0]) for n in NW}
        for kind in KINDS:
            lead = {n: [] for n in NW}
            for fv in NGRID:
                for n in NW:
                    L_, D_ = (pc(run(n, sc, delta=nudge(kind, fv, h))[0]) for h in (1, -1))
                    lead[n].append((L_ - D_) / base[n])
            for n in NW:
                star[s, kind, n] = cross(NGRID, np.array(lead[n]))
        say(f"  {s} lead-halving share f* (lead = (L - D) / L with no push; many sources: f*_G x sqrt(M)):")
        for n in NW:
            st = [star[s, k, n] for k in KINDS]
            fmt = lambda c: "never" if c is None else f"{c[0]:.1e}" + ("*" if c[2] > 1 else "")
            m3 = "never" if st[1] is None else f"{st[1][0] * np.sqrt(1e3):.1e}"
            star[s, "ln", n] = say(f"    {n:18} base L {base[n]:5.1f} | one src {fmt(st[0])} | f*_G {fmt(st[1])} | M=1e3 {m3} | "
                                   f"field {fmt(st[2])}")
    ok_n, det = True, []
    for n in NW:
        o = OLD_NUDGE[n]
        g, fld = star["main", "many sources", n], star["main", "field", n]
        gx, fx = star["x10", "many sources", n], star["x10", "field", n]
        m3 = None if g is None else g[0] * np.sqrt(1e3)
        val = lambda c: "never" if c is None else f"{c[0]:.1e}"
        md("Magnetite (4)", f"{short(n)} f*: one source / M=1e3 / field",
           f"{o[0]:.1e} / {o[2]:.1e} / {o[3]:.1e} (moon_nudges.log:{o[4]})",
           f"{val(star['main', 'one source', n])} / {'never' if m3 is None else f'{m3:.1e}'} / {val(fld)} "
           f"(x10 M=1e3 {'never' if gx is None else f'{gx[0] * np.sqrt(1e3):.1e}'}, field {val(fx)})", star["main", "ln", n])
        ok_n &= m3 is not None and fld is not None and 1 / 3 <= m3 / o[2] <= 3 and 1 / 3 <= fld[0] / o[3] <= 3
        det.append(f"{short(n)} M=1e3 {'never' if m3 is None else f'{m3:.1e}'} field {val(fld)}")
    L_, D_, P_ = ef["main", "as published (f_mag 1e-8..1e-2)"]
    ok = 43 <= L_ <= 48 and 43 <= D_ <= 48 and 6 <= P_ <= 8
    ok &= ecross["main"] is None or ecross["main"][0] >= 1e-14
    fE = star["main", "field", E]
    ok &= fE is not None and 1 / 3 <= fE[0] / 6.0e-16 <= 3
    check("P14", ok and ok_n, f"Earth L {L_:.1f} D {D_:.1f} pinned {P_:.1f}; L<50% crossing "
                              f"{'never' if ecross['main'] is None else f'{ecross[chr(109) + chr(97) + chr(105) + chr(110)][0]:.2e}'}; "
                              + "; ".join(det))


# ---------------------------------------------------------------- 11 H27 temperature
EA = {"Ala": 93e3, "Asp": 110e3}


def k_at(aa, T):
    return S.K_RAC[aa] * np.exp(EA[aa] / 8.314 * (1 / 273 - 1 / np.asarray(T, float)))


def section11():
    say("\n=== 11. H27 temperature: Early Earth T ~ U(273, 343) K; k_rac Arrhenius, g ~ 1/T; moons stay at 273 K ===")
    for aa in EA:
        say(f"  {aa}: k_rac x{k_at(aa, 308) / S.K_RAC[aa]:.0f} at 308 K, x{k_at(aa, 343) / S.K_RAC[aa]:.0f} at 343 K "
            f"(k_273 {S.K_RAC[aa]:.3g}/s, Ea {EA[aa] / 1e3:.0f} kJ/mol); g x{273 / 308:.3f} / x{273 / 343:.3f}")
    T = np.random.default_rng([SEED, 14, 27]).uniform(273, 343, N)
    x = W[E][0]
    kc = x["rate k2"] * x["concentration"]
    say(f"  Earth k2 c spans {kc.min():.1e}..{kc.max():.1e} /s")
    res = {}
    for s, sc in SCALES:
        w0 = BASE[E, s][0]
        res[s, "base 273 K"] = pc(w0)
        res[s, "273 K, cut on every draw"] = pc(w0 & (kc > 2 * S.K_RAC["Ala"]))
        res[s, "g only (T uniform)"] = pc(run(E, sc, delta=scaled(1.0, 273 / T))[0])
        for aa in EA:
            K = k_at(aa, T)
            res[s, f"{aa} k_rac, pipeline cut"] = pc(run(E, sc, krac=K)[0])
            wb = run(E, sc, delta=scaled(1.0, 273 / T), krac=K)[0]
            res[s, f"{aa} g + k_rac, pipeline cut"] = pc(wb)
            res[s, f"{aa} g + k_rac, cut on every draw"] = pc(wb & (kc > 2 * K))
            res[s, f"{aa} share of draws with k2 c <= 2 k_rac(T)"] = pc(kc <= 2 * K)
            for Tf in (308.0, 343.0):
                Kf = k_at(aa, Tf)
                wf = run(E, sc, delta=scaled(1.0, 273 / Tf), krac=Kf)[0]
                res[s, f"{aa} fixed {Tf:.0f} K, cut on every draw"] = pc(wf & (kc > 2 * Kf))
    lns = {}
    for k in [k for s, k in res if s == "main"]:
        lns[k] = say(f"  {k:46} main {res['main', k]:5.1f}  x10 {res['x10', k]:5.1f}")
    b = res["main", "base 273 K"]
    for k in ("g only (T uniform)", "Ala k_rac, pipeline cut", "Ala g + k_rac, cut on every draw",
              "Asp g + k_rac, cut on every draw", "Ala fixed 343 K, cut on every draw", "Asp fixed 343 K, cut on every draw"):
        md("H27 temperature", f"Earth {k}", f"{b:.1f} (273 K, this log)", f"{res['main', k]:.1f} (x10 {res['x10', k]:.1f})", lns[k])
    md("H27 temperature", "moons", "273 K", "unchanged", lns["g only (T uniform)"],
       "oceans under ice sit at the melting point; colder water only slows k_rac and raises g")
    d = {k: b - res["main", k] for k in [k for s, k in res if s == "main"]}
    ok = 1 <= d["Ala k_rac, pipeline cut"] <= 5 and 4 <= d["Ala g + k_rac, cut on every draw"] <= 10
    ok &= d["g only (T uniform)"] <= 0.3 and d["273 K, cut on every draw"] == 0
    ok &= 8 <= d["Asp g + k_rac, cut on every draw"] <= 25
    check("P15", ok, f"drops: Ala pipeline {d['Ala k_rac, pipeline cut']:.2f}, Ala all-draw "
                     f"{d['Ala g + k_rac, cut on every draw']:.2f}, g only {d['g only (T uniform)']:.2f}, 273 K all-draw "
                     f"{d['273 K, cut on every draw']:.2f}, Asp all-draw {d['Asp g + k_rac, cut on every draw']:.2f}; moons at 273 K by construction")


# ---------------------------------------------------------------- md
def write_md():
    out = ["# New analyses: base-model results rerun through the current pipeline (2026-10-12)", "",
           "Script `rerun_base.py` (predictions committed before the run), log `rerun_base.log`. Pipeline: "
           "`multiwall.multi` (sweep2 step f, sweep3 recount, deepest-hand cascade) on the seed-0 draws of "
           "`sweep2.chains`, N = 200,000, extras as `multiwall.run4`. \"main\" is the measured cascade time, "
           "\"x10\" the cascade x10 low end. Covers journals/ija_draft.md H19 (line 196), H20 (197) and H27 (204).", "",
           "Items marked *unchanged* cannot be moved by depth-dependent D_z (reason given); the base-model label can be "
           "dropped for them.", ""]
    for sec in dict.fromkeys(r[0] for r in MDROWS):
        out += [f"## {sec}", "", "| item | old (source) | new | log | note |", "|---|---|---|---|---|"]
        out += [f"| {r[1]} | {r[2]} | {r[3]} | `{r[4]}` | {r[5]} |" for r in MDROWS if r[0] == sec]
        out.append("")
    out += ["## H27 temperature sources (DOIs checked against Crossref)", "",
            "- Blake, Chang & Lepland 2010, Nature 464, 1029, 10.1038/nature08952: Archean ocean 26-35 C (phosphate oxygen isotopes).",
            "- Robert & Chaussidon 2006, Nature 443, 969, 10.1038/nature05239: ~70 C at 3.5 Ga (cherts, Si and O isotopes).",
            "- Knauth & Lowe 2003, GSA Bulletin 115, 566, 10.1130/0016-7606(2003)115<0566:HACTIF>2.0.CO;2: 55-85 C (Onverwacht cherts).",
            "", "The sampled range 273-343 K spans the cool and the hot readings.", "",
            "## Predictions", ""]
    out += [f"- {p}: {'held' if ok else 'FAILED'}" for p, ok in CHECKS]
    out.append("")
    with open(MD, "w") as f:
        f.write("\n".join(out))


def main():
    global FH
    FH = open(LOG, "w")
    say(f"rerun_base.py, seed {SEED}, N = {N}: H19/H20/H27 base-model results through multiwall.multi (main and cascade x10)")
    for sec in (section0, section1, section2, section3, section4, section5, section6, section7, section8, section9,
                section10, section11):
        sec()
    say("\nPredictions: " + ", ".join(f"{p} {'held' if ok else 'FAILED'}" for p, ok in CHECKS))
    write_md()
    say(f"wrote {os.path.relpath(MD, HERE)}")
    FH.close()


if __name__ == "__main__":
    main()
