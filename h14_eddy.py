"""H14: in what share of the sweep's draws is the energetic eddy size l larger than the wall width w = sqrt(D_h / k),
k = k2 * c, and does E86's x4 / x10 slowdown of the sideways coarsening time (coarsen_slow.py) cover those draws?
(2026-10-21; no existing file is edited, everything is imported.)

Draws: sweep2.chains(0) under the lid_melt seedfix, N = 200,000, 15 worlds, the same draws and the same extra inputs
(default_rng([0, world index])) as coarsen_slow.py. l is NOT in the sweep. It is not drawn either: the share of draws
with l > w is mean(w < l), a function of l, reported at U = 1e-3 / 0.02 / 1 m/s of the Rhines-type scale
l = sqrt(2 U / beta) (beta = 2 Omega / R at the equator; Ashkenazy and Tziperman 2021, doi 10.1038/s41467-021-26710-0,
40 km at U 0.02), at the literature size where one exists (Enceladus 3 km: Kang 2022 doi 10.1126/sciadv.abm4665, Zhang
2024 doi 10.1126/sciadv.adn6857; Europa 20-50 km: Ashkenazy 2021; Earth 25-250 km radius: Faghmous 2015
doi 10.1038/sdata.2015.28), as an average over log-uniform l, and with l tied to D_h (l = D_h / U, U log-uniform
1e-3..1 m/s). Omega = 2 pi / P with synchronous rotation ASSUMED, P from the NASA NSSDC satellite fact sheets (read
2026-10-21). Early Earth uses today's Omega.
Coverage: multiwall.multi is rerun at coarsening time x f for a ladder of f up to 1e10 (U.A = 5.5 / sqrt(f), as
coarsen_slow.py), four cases (rule main / x10, no lid / salt lid). Draws P_pved at f = 1 that stop being P_pved at f are
counted inside {l > w}, and compared with E86 (f = 4, 10) and with the break-even f. Stress test (NOT a result): if the
slowdown of each l > w draw were (l/w)^p, run exactly with a per-draw U.A: p = 1.1 is the smallest power that gives >= 4 at the
one measured point l/w = 3.6 (3.6^1.08 = 4), p = 2 is Da (13 at 3.6). Draws with l <= w get no slowdown. Section (7) (the
smallest f that moves any draw, a check of the 14x sentence in the coarsen_slow.py docstring) and the exact per-draw stress
runs were added after the first run of (1)-(6), whose numbers did not change; the first run used a ladder-threshold stress.

Run: UV_CACHE_DIR=~/Library/Caches/uv MPLBACKEND=Agg uv run --no-project --with numpy --with scipy --with matplotlib \
       --with numba python h14_eddy.py > h14_eddy.log

PREDICTIONS: h14_eddy_predictions.txt (v1, 2026-10-21 00:35 PDT), sha256
8be9ed63ef5f8742b2af2c97635de5ebf86b6279c969d30bf46a58fa013c04b6 (asserted below). Disclosure: the l table in the
predictions file came from a scratch arithmetic script (Omega, R, beta only, no draws). After the hash, one timing run of
multiwall.multi for Ganymede (f = 1 / 10 / 1e4, main no lid) was seen: P_pved 44.14 / 44.07 / 21.23%. P6's threshold
(300 draws at x10) was set before that; the predictions were not changed.
"""
import hashlib

import numpy as np

import lid_melt as LM   # applies the seedfix RNG correction before the project modules load
import uncertainty as U
import sweep2 as S
import sweep3 as S3   # noqa: F401  (multiwall imports it)
import multiwall as M

PRED_SHA = "8be9ed63ef5f8742b2af2c97635de5ebf86b6279c969d30bf46a58fa013c04b6"
with open("h14_eddy_predictions.txt", "rb") as fh:
    assert hashlib.sha256(fh.read()).hexdigest() == PRED_SHA, "h14_eddy_predictions.txt changed after the hash was recorded"

# synchronous rotation period, days. NSSDC Saturnian and Uranian satellite fact sheets (orbital period, all marked S);
# NSSDC Small Worlds fact sheet (rotation period, hours / 24); Earth 23.9 h from the same sheet.
PERIOD_D = {"Enceladus": 1.370218, "Mimas": 0.9424218, "Dione": 2.736915, "Titan": 15.945421,
            "Miranda": 1.413479, "Ariel": 2.520379, "Umbriel": 4.144176, "Titania": 8.705867, "Oberon": 13.463234,
            "Europa": 85.2 / 24, "Ganymede": 171.7 / 24, "Callisto": 400.5 / 24, "Pluto": 153.3 / 24,
            "Triton": 141.0 / 24, "Early Earth ocean": 23.9 / 24}
LIT = {"Enceladus": (3e3, 3e3), "Europa": (20e3, 50e3), "Early Earth ocean": (25e3, 250e3)}   # m, see docstring
UU = (1e-3, 0.02, 1.0)   # m/s: Kang 2022 / Ashkenazy 2021 / peak zonal jets (Soderlund 2019)
FL = (1.0, 2.0, 3.0, 4.0, 6.0, 10.0, 14.0, 20.0, 30.0, 50.0, 100.0, 200.0, 500.0, 1e3, 1e4, 1e6, 1e10)
KEYS = (("main", "none"), ("main", "salt"), ("x10", "none"), ("x10", "salt"))
SC = {"main": 1.0, "x10": 10.0}
TESTED = 3.6   # advect2d.py: the only measured l / w
A0 = U.A
N = U.N


def rhines(u, beta):
    return np.sqrt(2 * u / beta)


def ladder_masks(x, a, ex, salt):
    """{(rule, lid): [P_pved mask at each f in FL]}; U.A is restored."""
    none, out = np.zeros(N, bool), {}
    try:
        for f in FL:
            U.A = A0 / np.sqrt(f)
            for rule, lid in KEYS:
                out.setdefault((rule, lid), []).append(M.multi(x, a[1], a[4], ex, none if lid == "none" else salt, scale=SC[rule])[0])
    finally:
        U.A = A0
    return out


def stress_lost(r, l, p):
    """Draws P_pved at f = 1 that stop being P_pved when each draw with l > w has its coarsening slowed by (l/w)^p
    (per-draw U.A), per KEYS. p = 0 must lose nothing."""
    x, a, ex, salt = r["args"]
    none, f = np.zeros(N, bool), np.maximum(l / r["w"], 1.0) ** p
    U.A = A0 / np.sqrt(f)
    try:
        return [int((r["P"][(rule, lid)][0] & ~M.multi(x, a[1], a[4], ex, none if lid == "none" else salt, scale=SC[rule])[0]).sum())
                for rule, lid in KEYS]
    finally:
        U.A = A0


def selfcheck():
    assert abs(rhines(0.02, 2.6e-11) / 1e3 - 39.2) < 0.1                 # Ashkenazy 2021: "approx 40 km"
    assert abs(2 * (2 * np.pi / (PERIOD_D["Enceladus"] * 86400)) - 1.0614e-4) < 2e-7   # Soderlund 2019 Table 1: Omega 5.3e-5
    D, k2, c = np.array([4.0]), np.array([1.0]), np.array([1.0])        # hand case: w = sqrt(4 / 1) = 2 m
    w = np.sqrt(D / (k2 * c))
    assert w[0] == 2.0 and (w < 3.0)[0] and not (w < 1.5)[0]


def main():
    selfcheck()
    _, W = S.chains(0)
    names = list(W)
    print(f"h14_eddy.py, N = {N}, seed 0. l = eddy size, w = sqrt(D_h / (k2 c)); shares in % of all draws")
    R = {n: W[n][2][1] / np.pi for n in names}
    OM = {n: 2 * np.pi / (PERIOD_D[n] * 86400) for n in names}
    BE = {n: 2 * OM[n] / R[n] for n in names}
    LL = {n: [rhines(u, BE[n]) for u in UU] for n in names}
    assert 33e3 < LL["Europa"][1] < 45e3, "Europa typical l should be near Ashkenazy's 40 km"
    print("\n(1) Eddy sizes. P (d) synchronous rotation; R = L/pi of the sweep; beta = 2 Omega / R; l = sqrt(2 U / beta), km")
    print(f"{'world':18} {'P d':>7} {'Omega':>9} {'R km':>7} {'beta':>9} | l(1e-3) l(0.02) l(1) | literature km")
    for n in names:
        lit = f"{LIT[n][0] / 1e3:g}-{LIT[n][1] / 1e3:g}" if n in LIT else "-"
        print(f"{n:18} {PERIOD_D[n]:7.3f} {OM[n]:9.2e} {R[n] / 1e3:7.0f} {BE[n]:9.2e} | "
              + " ".join(f"{v / 1e3:7.1f}" for v in LL[n]) + f" | {lit}")

    # per-world draws, masks and shares
    res, WW = {}, {}
    for i, n in enumerate(names):
        x, win, a, salt = W[n]
        rng = np.random.default_rng([0, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, N), "r": 10 ** rng.uniform(-2, 0, N), "rho": 10 ** rng.uniform(0, 2, N)}
        Dh, k2, c, t = x["horizontal D_h"], x["rate k2"], x["concentration"], x["time"]
        w = np.sqrt(Dh / (k2 * c))
        sh = lambda l, m=None, w=w: 100 * (w < l).mean() if m is None else 100 * (w[m] < l).mean()   # w bound now (loop var)
        tied = 10 ** np.random.default_rng([0, i, 14]).uniform(-3, 0, N)
        P = ladder_masks(x, a, ex, salt)
        r = {"w": w, "P": P, "Dh": Dh, "t": t, "L": a[1], "sh": sh, "args": (x, a, ex, salt)}
        r["typ"] = [sh(l) for l in LL[n]]
        r["avg"] = np.mean([sh(l) for l in np.geomspace(LL[n][0], LL[n][2], 41)])
        r["lit"] = [sh(l) for l in LIT[n]] if n in LIT else None
        r["tied"] = 100 * ((Dh / tied) > w).mean()
        r["speed_ok"] = (Dh / LL[n][1] >= 1e-3) & (Dh / LL[n][1] <= 1.0)   # D_h / l(0.02) is a speed in 1e-3..1 m/s
        r["wmed"] = np.median(w) / 1e3
        res[n] = r
    # f = 1 regression against coarsen_slow.log (logged values, main no lid / final / x10 no lid / final)
    d = max(abs(100 * res[n]["P"][k][0].mean() - v) for n in names for k, v in zip(KEYS, LM.LOGGED[n]))
    assert d <= 0.05, f"f = 1 does not reproduce the logged P_pved (max diff {d:.2f})"
    print(f"\nregression: f = 1 P_pved reproduces coarsen_slow.log / multiwall.log within 0.05 (max diff {d:.3f} points)")

    print("\n(2) Share of ALL draws with l > w. typ = l at U 1e-3 / 0.02 / 1 m/s; avg = over log-uniform l between the ends;"
          "\n    lit = literature l (low, high); tied = l = D_h/U; kin = share of draws where D_h/l(0.02) is a speed in 1e-3..1 m/s")
    print(f"{'world':18} {'median w km':>11} | {'U 1e-3':>6} {'U 0.02':>6} {'U 1':>6} {'avg':>6} | {'lit lo':>6} {'lit hi':>6} | {'tied':>5} | {'kin':>5}")
    for n in names:
        r = res[n]
        lit = f"{r['lit'][0]:6.1f} {r['lit'][1]:6.1f}" if r["lit"] else f"{'-':>6} {'-':>6}"
        print(f"{n:18} {r['wmed']:11.0f} | {r['typ'][0]:6.1f} {r['typ'][1]:6.1f} {r['typ'][2]:6.1f} {r['avg']:6.1f} | {lit} | {r['tied']:5.1f} | {100 * r['speed_ok'].mean():5.1f}")

    print("\n(3) Among P_pved draws at f = 1 (main, no lid): share with l > w at U 0.02 / 1 m/s, against all draws; n = P_pved draws")
    sub_sh = {}
    for n in names:
        r = res[n]
        p0 = r["P"][("main", "none")][0]
        sub_sh[n] = [r["sh"](l, p0) for l in LL[n][1:]]
        print(f"{n:18} n {p0.sum():6d} | U 0.02: {sub_sh[n][0]:5.1f} (all {r['typ'][1]:5.1f}) | U 1: {sub_sh[n][1]:5.1f} (all {r['typ'][2]:5.1f})")

    print("\n(4) l/w among draws with l > w, at U 0.02 and U 1: share beyond the tested 3.6, median, 90th percentile")
    rat = {}
    for n in names:
        r = res[n]
        rat[n] = []
        row = []
        for l in LL[n][1:]:
            q = l / r["w"][r["w"] < l]
            rat[n].append(100 * (q > TESTED).mean())
            row.append(f"n {len(q):6d} >3.6: {rat[n][-1]:5.1f}%  med {np.median(q):7.1f}  p90 {np.percentile(q, 90):8.1f}")
        print(f"{n:18} U 0.02 {row[0]} | U 1 {row[1]}")

    print("\n(5) E86 coverage. Subset = {l > w} and P_pved at f = 1. Lost = in the subset and not P_pved at f. Draws; (% of N).")
    print("    Column order per f: main no lid, main final, x10 no lid, x10 final. Also overall loss (all draws) at f = 10 for the same cases.")
    cov = {}
    for n in names:
        r = res[n]
        P, w = r["P"], r["w"]
        # non-monotone draws (P_pved again at a larger f): reported, not hidden
        gain = max(int(sum((P[k][j + 1] & ~P[k][j]).sum() for j in range(len(FL) - 1))) for k in KEYS)
        flv = {}
        for k in KEYS:
            fl = np.full(N, np.inf)
            for j in range(len(FL) - 1, 0, -1):
                fl[P[k][0] & ~P[k][j]] = FL[j]
            flv[k] = fl
        r["gain"], r["fl"] = gain, flv
        tot10 = [int((P[k][0] & ~P[k][FL.index(10.0)]).sum()) for k in KEYS]
        smin = float(np.min(r["t"] * r["Dh"] / (r["L"] / A0) ** 2))
        print(f"\n{n}: non-monotone draws {gain}; overall lost at x10 {tot10} draws; min over ALL draws of t D_h / (L/A)^2 = {smin:.1f}")
        for lab, l in (("U 0.02", LL[n][1]), ("U 1", LL[n][2])) + ((("lit hi", LIT[n][1]),) if n in LIT else ()):
            sub = w < l
            base = [sub & P[k][0] for k in KEYS]
            nb = [int(b.sum()) for b in base]
            out = {}
            for f in (4.0, 10.0, 100.0, 1000.0, 1e10):
                j = FL.index(f)
                out[f] = [int((base[q] & ~P[k][j]).sum()) for q, k in enumerate(KEYS)]
            # break-even f: the FIRST FL rung where >= 5% of the {l > w} subset is lost, so an UPPER bound on the true 5% point
            #   (every earlier rung stays under 5%). None = not reached by f = 1e10. The verdict below uses min() of these rungs.
            #   The label printed below is regex-parsed by h14b_advect.py:424; keep its wording.
            f5 = []
            for q in (0, 1):
                k = KEYS[q]
                tot = max(nb[q], 1)
                hit = [f for j, f in enumerate(FL) if (base[q] & ~P[k][j]).sum() / tot >= 0.05]
                f5.append(hit[0] if hit else None)
            stress = {p: stress_lost(r, l, p) for p in (1.1, 2.0)}   # draws lost, per KEYS
            cov[(n, lab)] = dict(nb=nb, out=out, f5=f5, stress=stress, sub=sub)
            print(f"  l at {lab:6} subset P_pved n {nb}; lost at x4 {out[4.0]}, x10 {out[10.0]}, x100 {out[100.0][:2]}, "
                  f"x1000 {out[1000.0][:2]}, x1e10 {out[1e10][:2]}")
            print(f"           break-even f (5% of subset lost): main no lid {f5[0]}, main final {f5[1]}; "
                  f"stress (NOT a result) slowdown (l/w)^p, draws lost [main no lid, main final]: p 1.1 {stress[1.1][:2]} = "
                  f"{[round(100 * v / N, 2) for v in stress[1.1][:2]]} points; p 2 {stress[2.0][:2]} = {[round(100 * v / N, 2) for v in stress[2.0][:2]]} points")

    assert stress_lost(res["Ganymede"], LL["Ganymede"][1], 0.0) == [0, 0, 0, 0]   # p = 0 is the f = 1 baseline

    print("\n(7) Check of the coarsen_slow.py docstring (coarsening <= 7% of the ocean's age, 'must be about 14x slower before it decides"
          "\n    a draw'): all draws, main no lid. s* = min over draws of t D_h / (L/A)^2 (coarsening is <= 1/s* of the age at f = 1).")
    print(f"{'world':18} {'min s*':>7} | lost at f = " + " ".join(f"{f:>6g}" for f in (2, 3, 4, 6, 10, 14, 20)) + " | smallest ladder f that moves any draw")
    for n in names:
        r = res[n]
        P0 = r["P"][("main", "none")]
        lost = [int((P0[0] & ~P0[FL.index(float(f))]).sum()) for f in (2, 3, 4, 6, 10, 14, 20)]
        first = next((f for j, f in enumerate(FL) if (P0[0] & ~P0[j]).any()), None)
        r["first"] = first
        print(f"{n:18} {np.min(r['t'] * r['Dh'] / (r['L'] / A0) ** 2):7.1f} | {'':11}" + " ".join(f"{v:6d}" for v in lost) + f" | {first}")

    print("\n(6) Verdict per world at typical l (U = 0.02): E86 covers if P6 holds and >= 50% of the l > w draws have l/w <= 3.6")
    verdict = {}
    for i, n in enumerate(names):
        cv = cov[(n, "U 0.02")]
        p6 = all(cv["out"][10.0][q] <= 300 and cv["out"][4.0][q] <= 200 for q in range(4))
        inside = 100 - rat[n][0]
        f5 = [f for f in cv["f5"] if f is not None]
        if not p6:
            verdict[n] = "NOT COVERED"
        elif inside >= 50:
            verdict[n] = "COVERED"
        else:
            verdict[n] = f"COVERED ONLY IF slowdown beyond l/w = 3.6 stays below x{min(f5) if f5 else '>1e10'}"
        print(f"  {n:18} share l>w {res[n]['typ'][1]:5.1f}%  inside test {inside:5.1f}%  P6 {'yes' if p6 else 'NO'}  -> {verdict[n]}")

    C = []

    def check(p, ok, d):
        C.append(ok)
        print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {d}")

    print("\nPredictions (h14_eddy_predictions.txt v1):")
    sh = lambda n: res[n]["lit"]
    e3, eu, ee = sh("Enceladus")[0], sh("Europa"), sh("Early Earth ocean")
    check("P1", 5 <= e3 <= 15 and 15 <= eu[0] <= 30 and 22 <= eu[1] <= 38 and 6 <= ee[0] <= 20 and 30 <= ee[1] <= 52,
          f"Enceladus 3 km {e3:.1f}; Europa 20/50 km {eu[0]:.1f}/{eu[1]:.1f}; Early Earth 25/250 km {ee[0]:.1f}/{ee[1]:.1f}")
    typ = {n: res[n]["typ"][1] for n in names}
    check("P2", all(8 <= v <= 45 for v in typ.values()) and typ["Titan"] > typ["Europa"] > typ["Enceladus"] and typ["Callisto"] > typ["Europa"],
          f"range {min(typ.values()):.1f}-{max(typ.values()):.1f}; Titan {typ['Titan']:.1f}, Callisto {typ['Callisto']:.1f}, Europa {typ['Europa']:.1f}, Enceladus {typ['Enceladus']:.1f}")
    hi50 = {n for n in names if res[n]["typ"][2] > 50}
    allowed = {"Ganymede", "Titan", "Callisto", "Oberon", "Triton", "Pluto", "Titania"}
    check("P3", max(typ.values()) <= 50 and {"Titan", "Callisto"} <= hi50 <= allowed,
          f"max at U 0.02 {max(typ.values()):.1f}; worlds above 50% at U 1: {sorted(hi50)}")
    moons = [n for n in names if n != "Early Earth ocean"]
    check("P4", all(res[n]["tied"] < 8 for n in moons) and 8 <= res["Early Earth ocean"]["tied"] <= 35 and all(res[n]["tied"] < typ[n] for n in names),
          f"moons max {max(res[n]['tied'] for n in moons):.1f}; Early Earth {res['Early Earth ocean']['tied']:.1f}; tied < typical in {sum(res[n]['tied'] < typ[n] for n in names)} of 15")
    up = sum(sub_sh[n][0] > typ[n] for n in names)
    dmax = max(abs(sub_sh[n][0] - typ[n]) for n in names)
    check("P5", up >= 10 and dmax < 20, f"P_pved share higher than all-draw share in {up} of 15; max |difference| {dmax:.1f} points")
    l10 = max(cov[(n, "U 0.02")]["out"][10.0][q] for n in names for q in range(4))
    l4 = max(cov[(n, "U 0.02")]["out"][4.0][q] for n in names for q in range(4))
    check("P6", l10 <= 300 and l4 <= 200, f"largest lost in the l > w subset: {l4} draws at x4, {l10} at x10 (limits 200 / 300 of N = {N})")
    fr = lambda n, f, q: cov[(n, "U 0.02")]["out"][f][q] / max(cov[(n, "U 0.02")]["nb"][q], 1)
    m100 = max(fr(n, 100.0, q) for n in names for q in (0, 1))
    m1000 = max(fr(n, 1000.0, q) for n in names for q in (0, 1))
    check("P7", m100 <= 0.10 and m1000 <= 0.50, f"largest share of the subset lost: {100 * m100:.1f}% at x100, {100 * m1000:.1f}% at x1000")
    check("P8", all(35 <= rat[n][0] <= 75 for n in names),
          f"share of l > w draws beyond l/w = 3.6: {min(rat[n][0] for n in names):.1f}-{max(rat[n][0] for n in names):.1f}%")
    print(f"\n{sum(C)} of {len(C)} held.")


if __name__ == "__main__":
    main()
