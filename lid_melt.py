"""Two lid questions from K. Soderlund (email, 2026-10-14), run through multiwall's present rule (2026-10-14).

(1) MELT LIDS. Melting of the ice shell can also make a stably stratified layer under the ice-ocean interface, on
    Europa, Ganymede and Enceladus alike (also Ashkenazy & Tziperman 2021 Nat. Commun. 10.1038/s41467-021-26710-0).
    The salinity rule of Section 2.5 does not include this. Here every world gets an extra lid with probability
    p_melt = 0, 0.1, 0.25, 0.5, drawn independently of every other input and of the salt lid (no measurement ties a
    melt layer to D_z, D_h, time, chemistry or salinity). A melt lid is treated exactly like a salt lid (it pins).
    Break-even: how large must Ganymede's p_melt be (with none on Europa) before Ganymede falls below Europa?
(2) COMPOSITION. S* comes from TEOS-10 (Earth-like, mostly NaCl seawater). Sulfate- or carbonate-rich oceans would
    move it. No equation of state for those is used here; S* is scaled by f = 0.5, 0.75, 1.5, 2, 4 for Europa
    (20 km ice, S* 3.0; and Juno's 29 +/- 10 km shell as in robust_ranking.py) and Enceladus (S* 20). The 28.5 MPa
    cutoff is kept, so Ganymede (100+ km of ice) keeps no salt lid at every f.

Base: multiwall.py Section 4, seed 0, rules "main" and cascade x10 (multiwall.log:540-555). Enceladus extras use the
seedfix.py correction (default_rng([s, 0]) -> [s, 0, 1]); its no-lid shares are 26.4 / 19.4 (seedfix.log:9-10).
Melt and new salinity draws use their own streams default_rng([0, i, 7]) and default_rng([0, i, 8]).
Run: nice -n 10 /opt/anaconda3/bin/python lid_melt.py > lid_melt.log

PREDICTIONS (written 2026-10-14, before any run):
 (P1) Regression: p_melt = 0 reproduces multiwall.log main and x10 (no lid and final) within 0.05 points for worlds
      1-14, and seedfix's Enceladus 26.4 / 19.4 (no lid) and 0.0 (final).
 (P2) Independence: each world's final share at p_melt = p equals final(0) x (1 - p) within 0.4 points.
 (P3) With the same p_melt on every world the order Ganymede > Europa > Enceladus never changes.
 (P4) Break-even (Europa p_melt = 0, 20 km ice): Ganymede falls below Europa only above p_melt = 0.30-0.36 (main and
      x10). Against Europa with the Juno shell, 0.15-0.25.
 (P5) At p_melt = 0.5 on all: Ganymede 16-23%, Europa 10-16%, Enceladus 0%.
 (P6) Composition, main rule: Europa (20 km) final ~36 / 24 / 18% at f = 0.5 / 2 / 4 (pinned fraction
      ln(3f)/ln(100)); Enceladus ~13% at f = 0.5, ~5-6% at f = 0.75, 0% at f >= 1. Ganymede stays first at every f.
 (P7) The salt-lid fraction of the fresh salinity draws at f = 1 matches the logged one within 0.5 points
      (Europa 24%, Enceladus 100%).
"""
import os
import sys

import numpy as np

ORIG_RNG = np.random.default_rng


def fixed_rng(seed=None, *a, **kw):   # seedfix.py's correction, copied (importing seedfix would run it)
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng
import uncertainty as U   # noqa: E402
U.N = int(os.environ.get("LIDMELT_N", U.N))   # smoke test only
import sweep2 as S        # noqa: E402
import sweep3 as S3       # noqa: E402
import multiwall as M     # noqa: E402
import ocean_worlds as OW # noqa: E402
import robust_ranking as RR  # noqa: E402

PS = (0.0, 0.1, 0.25, 0.5)
FS = (0.5, 0.75, 1.0, 1.5, 2.0, 4.0)
LOGGED = {  # multiwall.log:541-555 main / x10, no lid | final; Enceladus from seedfix.log:9-10
    "Early Earth ocean": (99.9, 99.9, 99.7, 99.7), "Umbriel": (53.7, 53.7, 43.6, 43.6),
    "Oberon": (51.4, 51.4, 40.8, 40.8), "Titania": (51.3, 51.3, 40.9, 40.9), "Ariel": (49.6, 49.6, 39.7, 39.7),
    "Ganymede": (44.1, 44.1, 33.5, 33.5), "Callisto": (43.5, 43.5, 33.0, 33.0), "Pluto": (41.1, 41.1, 30.6, 30.6),
    "Triton": (37.0, 37.0, 26.5, 26.5), "Titan": (35.8, 35.8, 25.5, 25.5), "Europa": (39.6, 30.1, 29.0, 22.0),
    "Dione": (42.0, 26.2, 32.2, 20.1), "Miranda": (27.1, 9.0, 21.3, 7.0), "Mimas": (17.3, 5.6, 11.6, 3.7),
    "Enceladus": (26.4, 0.0, 19.4, 0.0)}
SALT = {**OW.STRAT, **{n: (w[4], w[5]) for n, w in OW.WORLDS.items() if w[5]}}
CHECK = []


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


def main(seed=0):
    _, W = S.chains(seed)
    names = list(W)
    res, comp = {}, {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        u = np.random.default_rng([seed, i, 7]).random(U.N)
        none = np.zeros(U.N, bool)
        r = {}
        for sc, lab in ((1.0, "main"), (10.0, "x10")):
            r[lab, "nolid"] = 100 * M.multi(x, a[1], a[4], ex, none, scale=sc)[0].mean()
            for p in PS:
                r[lab, p] = 100 * M.multi(x, a[1], a[4], ex, salt | (u < p), scale=sc)[0].mean()
        res[name] = r
        if name in ("Europa", "Enceladus"):
            (s_lo, s_hi), s_star = SALT[name]
            g = np.random.default_rng([seed, i, 8])
            S_new = np.exp(g.uniform(np.log(s_lo), np.log(s_hi), U.N))
            h = np.clip(g.normal(29e3, 10e3, U.N), 3e3, 60e3)
            juno = np.interp(920 * 1.315 * h / 1e6, RR.SSTAR_P, RR.SSTAR_S, right=0)
            c = {"logged_lid": 100 * salt.mean()}
            for f in FS:
                for key, ss in (("base", f * s_star), ("juno", f * juno)):
                    if key == "juno" and name != "Europa":
                        continue
                    lid = S_new < ss
                    c[key, f] = (100 * lid.mean(),) + tuple(
                        100 * M.multi(x, a[1], a[4], ex, lid, scale=sc)[0].mean() for sc in (1.0, 10.0))
            comp[name] = c
    report(names, res, comp)


def report(names, res, comp):
    print(f"lid_melt.py, N = {U.N}, seed 0. P_pved %, rule main | cascade x10.")
    print("\n(1) Melt lid probability p (independent of the salt lid). Columns: no lid, final at p = " +
          ", ".join(str(p) for p in PS))
    for n in names:
        r = res[n]
        print(f"{n:18} main {r['main', 'nolid']:5.1f} | " + " ".join(f"{r['main', p]:5.1f}" for p in PS) +
              f"    x10 {r['x10', 'nolid']:5.1f} | " + " ".join(f"{r['x10', p]:5.1f}" for p in PS))
    print("\nBreak-even p_melt on Ganymede (no melt lid on Europa) before Ganymede < Europa:")
    be = {}
    for lab in ("main", "x10"):
        be[lab] = 1 - res["Europa"][lab, 0.0] / res["Ganymede"][lab, 0.0]
        print(f"  {lab}: Europa 20 km {res['Europa'][lab, 0.0]:.1f} vs Ganymede {res['Ganymede'][lab, 0.0]:.1f} -> p > {be[lab]:.2f}")
    k = {"main": 1, "x10": 2}
    for lab in ("main", "x10"):
        ej = comp["Europa"]["juno", 1.0][k[lab]]
        be[lab, "juno"] = 1 - ej / res["Ganymede"][lab, 0.0]
        print(f"  {lab}: Europa Juno shell {ej:.1f} -> p > {be[lab, 'juno']:.2f}")

    print("\n(2) Composition: S* scaled by f (28.5 MPa cutoff kept). lid % | final main | final x10")
    for n in ("Europa", "Enceladus"):
        c = comp[n]
        print(f"{n} (logged lid fraction {c['logged_lid']:.1f}%):")
        for f in FS:
            line = f"  f {f:4}: " + "{:5.1f}% | {:5.1f} | {:5.1f}".format(*c["base", f])
            if ("juno", f) in c:
                line += "    Juno shell: " + "{:5.1f}% | {:5.1f} | {:5.1f}".format(*c["juno", f])
            print(line)
    print(f"Ganymede: no salt lid at any f (ice > 23 km): {res['Ganymede']['main', 0.0]:.1f} | {res['Ganymede']['x10', 0.0]:.1f}")

    print("\nPredictions:")
    d = max(max(abs(res[n]["main", "nolid"] - LOGGED[n][0]), abs(res[n]["main", 0.0] - LOGGED[n][1]),
                abs(res[n]["x10", "nolid"] - LOGGED[n][2]), abs(res[n]["x10", 0.0] - LOGGED[n][3])) for n in names)
    check("P1", d <= 0.05, f"max |p=0 - logged| = {d:.2f} points over 15 worlds, main and x10")
    d = max(abs(res[n][lab, p] - res[n][lab, 0.0] * (1 - p)) for n in names for lab in ("main", "x10") for p in PS)
    check("P2", d <= 0.4, f"max |final(p) - final(0)(1-p)| = {d:.2f} points")
    ok = all(res["Ganymede"][lab, p] > res["Europa"][lab, p] > res["Enceladus"][lab, p] or
             (res["Enceladus"][lab, p] == 0 and res["Ganymede"][lab, p] > res["Europa"][lab, p])
             for lab in ("main", "x10") for p in PS)
    check("P3", ok, "Ganymede > Europa > Enceladus at every common p (Enceladus 0 counts as last)")
    ok = all(0.30 <= be[l] <= 0.36 for l in ("main", "x10")) and all(0.15 <= be[l, "juno"] <= 0.25 for l in ("main", "x10"))
    check("P4", ok, f"break-even 20 km {be['main']:.2f} / {be['x10']:.2f}; Juno {be['main', 'juno']:.2f} / {be['x10', 'juno']:.2f}")
    g, e, n = (res[w]["main", 0.5] for w in ("Ganymede", "Europa", "Enceladus"))
    g2, e2 = (res[w]["x10", 0.5] for w in ("Ganymede", "Europa"))
    ok = all(16 <= v <= 23 for v in (g, g2)) and all(10 <= v <= 16 for v in (e, e2)) and n == 0 == res["Enceladus"]["x10", 0.5]
    check("P5", ok, f"p = 0.5: Ganymede {g:.1f}/{g2:.1f}, Europa {e:.1f}/{e2:.1f}, Enceladus {n:.1f}")
    ce, cn = comp["Europa"], comp["Enceladus"]
    eu = [ce["base", f][1] for f in (0.5, 2.0, 4.0)]
    en = [cn["base", f][1] for f in FS]
    ok = (abs(eu[0] - 36) <= 1.5 and abs(eu[1] - 24) <= 1.5 and abs(eu[2] - 18) <= 1.5 and 11.5 <= en[0] <= 14.5
          and 4.5 <= en[1] <= 6.5 and max(en[2:]) == 0 and
          all(res["Ganymede"]["main", 0.0] > max(ce["base", f][1], ce["juno", f][1], cn["base", f][1]) for f in FS))
    check("P6", ok, f"Europa f 0.5/2/4: {eu[0]:.1f}/{eu[1]:.1f}/{eu[2]:.1f}; Enceladus by f: "
          + "/".join(f"{v:.1f}" for v in en))
    d = max(abs(c["base", 1.0][0] - c["logged_lid"]) for c in comp.values())
    check("P7", d <= 0.5, f"fresh vs logged lid fraction, max diff {d:.2f} points")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held.")


if __name__ == "__main__":
    main()
