"""sweep2: the planetary success-share sweep rerun with seven additions (2026-10-12).

Same model, draws and seed as uncertainty.py (chain A: Enceladus, Europa, Early Earth) and ocean_worlds.py (chain B:
strat_sweep order, then the 12 moons with their salt-lid draws and burned depth-sensitivity sweeps). Every
uncertainty.py draw is reused (common random numbers); new inputs come from default_rng([seed, world index]).
Steps are cumulative: 0 (as published) -> a -> b -> c -> d -> e -> f; the final step carries split g.

(a) PVED range from the literature (bias g = |dE_pv| / (R*273 K); 1 aeV = 9.6485e-14 J/mol).
    Old prior: log-uniform 4e-19..4e-16. New: log-uniform 4.4e-18..2.64e-16 (same quantile, remapped).
      #  molecule / method                       |dE_pv|                 g(273 K)   sign           source (how checked)
      1  L-Ala aq. conformer, Gly, peptides;     ~1e-14 J/mol (0.10 aeV) 4.4e-18    L lower        Mason & Tranter 1985 Proc R Soc A
         first-generation ab initio                                                                10.1098/rspa.1985.0003 (ABSTRACT ONLY)
      2  Ala, Val, Ser, Asp, first generation    ~1e-14 J/mol            4.4e-18    L in aq. conf. Tranter 1985 Mol Phys 10.1080/00268978500102741
                                                                                    (conformer-     (Crossref + ABSTRACT ONLY)
                                                                                     sensitive)
      3  alanine (Bakasov 1998 / Berger & Quack  6.2 aeV = 6.0e-13 J/mol 2.64e-16   undetermined   Quack, Seyfang & Wichmann 2022 Chem Sci
         2000, tabulated)                                                                          10.1039/d2sc01323a (FULL TEXT, PMC9491092)
      4  Ala gas, 13 conformers                  not in abstract         -          L lower 7/13   Laerdahl et al. 2000 ChemPhysChem 1:60 (ABSTRACT)
      5  Ala                                     not in abstract         -          "L stabilisation unjustified"  Berger & Quack 2000 ChemPhysChem 1:57 (ABSTRACT)
      6  amino acids, second generation          ~10x Mason-Tranter      -          gas L-Ala -, solution sign undetermined
                                                                                                   MacDermott et al. 2009 OLEB 10.1007/s11084-009-9161-x (ABSTRACT)
      7  alpha-methyl amino acids                not stated              -          L negative     MacDermott 2009 10.1007/s11084-009-9162-9 (ABSTRACT)
      8  amino acids, relativistic               not stated              -          no clear L     Wesendrup et al. 2003 JPCA 10.1021/jp022568v (ABSTRACT)
      9  generic light chiral molecule           ~1e-13 J/mol (1.0 aeV)  4.4e-17    -              Lente 2006 JPCA 10.1021/jp0658344 (ABSTRACT ONLY)
     10  hydrated D-glyceraldehyde               not in abstract         -          D lower        Tranter 1986 JCS Chem Commun 10.1039/c39860000060 (ABSTRACT)
     11  deoxyribose, ribose, other pentoses     not in abstract         -          conformer-dependent (D for C2-endo, L for C3-endo)
                                                                                                   Tranter et al. 1992 Proc R Soc A 10.1098/rspa.1992.0037 (ABSTRACT)
     12  Gly, Ala, Ser basis-set study           -                       -          basis-set dependent  Kikuchi & Wang 1990 BCSJ 10.1246/bcsj.63.2751 (ABSTRACT)
    g_min = row 1 (smallest published amino-acid magnitude; not a physical floor: aqueous conformer averaging can be
    smaller and the sign is undetermined). g_max = row 3 (largest verified amino-acid value) = the alanine value.
    Quack 2022 full text: light-atom molecules are "aeV to feV"; its abstract's 1e-11..1e-10 J/mol is for heavier
    molecules, so applying that to amino acids is a misapplication. No sugar magnitude was accessible.
(b) rate bias eta: transition-state PVED biases the rates, k_L/k_D = exp(-ddE_pv^TS / RT); g_eff = eta * g with eta
    log-uniform 0.1-10 (ASSUMED range; no published amino-acid transition-state value found).
(c) D_z floor 1e-9 m^2/s: icy-moon D_z prior 1e-10..1e-3 -> log-uniform 1e-9..1e-3 (same quantile, remapped).
    Earth's 1e-6..1e-4 is untouched.
(d) racemization (review.py R7: no choice unless k > 2 k_rac, k = k2*c). Free alanine in water at 3 C:
    k = 0.046e-5 /yr (Bada & Schroeder 1975 Naturwiss. 10.1007/bf00592179, as tabulated in Steen, Jorgensen &
    Lomstein 2013 PLoS ONE 10.1371/journal.pone.0071648 Table 3, FULL TEXT PMC3741123), taken to 273 K with Steen's
    consensus Ea = 93 kJ/mol (fit to sediment-bound data): k_rac = 2.88e-7 /yr = 9.13e-15 /s. Sensitivity: Asp
    (0.15e-5 /yr, Ea 110 kJ/mol) 2.73e-14 /s. All worlds at 273 K (ASSUMED; no ocean temperatures in the scripts).
(e) feedstock a (the drawn concentration; sets k = k2*a) vs product p (sets N = p V f 1000 N_A), p/a log-uniform
    0.01-1 (ASSUMED): delta_at(a, g, k2, V*p/a, t, f).
(f) depth-dependent D_z: bottom/top ratio rho log-uniform 1-100 (ASSUMED), exponential profile, so sigma ~ sqrt(D)
    gives |d ln sigma / dr| = ln(rho) / (2H). Eq. 20 (two_zone.py P3, |d ln sigma/dr| < 2/r) at the top-layer wall,
    r_w = R - d_w, d_w = min(H, L sqrt(D_z/D_h)): a column draw is pinned if ln(rho) >= 4H/r_w. Non-column draws
    have no stacked walls and are not tested. Sensitivity: threshold x5 (two_zone.log onset lies 2-5x above the
    sharp criterion in ln-contrast). The drawn D_z still sets the merge time (ASSUMED: it is the profile's mean).
(g) outcome split: P_pved = win; P_chance = one hand reached (finishes, not racemized, layers merge, not pinned,
    no salt lid) but F*Delta < 2, so the hand is set by chance; P_none = the rest. Unlimited bias (g -> inf) =
    P_pved + P_chance.
Scale rule (Eq. 13, review.n_needed): N = 2/(sqrt(pi) g^2 sqrt(k tau)) molecules at k tau = 10.

PREDICTIONS (written 2026-10-12, before any full run; points are percentage points of the no-lid share unless
"final" = after the salt lid):
 (P1) (a) PVED range: Enceladus no-lid +1 to +3 points over step 0; Europa +0.3 to +1.5; Earth within 0.3.
 (P2) (b) eta moves no world by more than 1 point relative to step a (no-lid and final).
 (P3) (c) D_z floor: every icy-moon share (no-lid and final, where nonzero) is multiplied by 1.12-1.21 relative to
      step b (exactly 7/6 in expectation: no draw with D_z < 1e-9 wins); Earth within 0.3.
 (P4) (d) racemization (Ala): Earth unchanged (min k2*c = 1e-13 /s > 2 k_rac); the 4-Gyr worlds lose 5-15 % of
      their step-c share (relative); Enceladus loses < 6 % relative; Mimas < 3 % relative.
 (P5) (e) p/a: Enceladus no-lid -1 to -4 points; Europa no-lid -0.3 to -2; Earth moves < 1 point.
 (P6) (f) depth pinning is the largest no-lid change of steps a-f for every world: Earth 36-42 %; Enceladus no-lid
      4-10 % (final 0); Europa no-lid 6-12 %, final 4.5-9 %; Ganymede and Callisto final 10-16 %; the four Uranian
      moons final 8-14 %; Miranda no-lid 9-16 %, final 3-6 %; Mimas final <= 2 %.
 (P7) (g) at the final step P_chance (final, i.e. unlimited minus P_pved) is < 3 points for Europa, < 1 point for
      Earth, 1-4 points for Enceladus no-lid; so the unlimited column for Earth is within 1 point of P_pved.
 (P8) Ganymede stays the highest confirmed ocean at the final step: Ganymede > Europa > Enceladus = 0.
 (P9) Seeds 1-3 reproduce every final-step share (no-lid and final) within 0.5 points of seed 0.
"""
import sys

import numpy as np

import ocean_worlds as OW
import uncertainty as U
from ocean import NA, PREF, delta_at

G_OLD, G_NEW = (4e-19, 4e-16), (4.4e-18, 2.64e-16)
DZ_FLOOR = 1e-9
K_RAC = {"Ala": 9.13e-15, "Asp": 2.73e-14}  # 1/s at 273 K
STEPS = "0abcdef"
G_ALA = 2.64e-16  # Quack 2022 table, alanine 6.2 aeV (= g_max); first-generation alanine = g_min
# logged shares (no lid, final) the all-off step must reproduce, seed 0
REG_B = {"Umbriel": (38.7, 38.7), "Ariel": (37.0, 37.0), "Oberon": (36.5, 36.5), "Titania": (36.4, 36.4),
         "Pluto": (30.2, 30.2), "Callisto": (28.5, 28.5), "Ganymede": (28.3, 28.3), "Triton": (27.7, 27.7),
         "Titan": (24.3, 24.3), "Europa": (28.6, 21.7), "Dione": (33.6, 21.0), "Miranda": (26.9, 8.9),
         "Mimas": (14.7, 4.7), "Enceladus": (21.8, 0.0)}           # ocean_worlds.log, strat_sweep.log
REG_A = {"Enceladus": (21.8, 28.6), "Europa": (28.6, 31.2), "Early Earth ocean": (74.7, 74.7)}  # assumed_rules.log: win, unlimited


def remap(v, old, new):
    """Same quantile of a log-uniform draw, moved to a new log-uniform range."""
    return new[0] * (new[1] / new[0]) ** (np.log(v / old[0]) / np.log(old[1] / old[0]))


def outcome(x, L, dz_rng, step, ex, lid, eq20=1.0, krac=K_RAC["Ala"]):
    """(win, reached one hand) with the options of `step` and all earlier steps on; step '0' is uncertainty.sweep."""
    on = lambda s: STEPS.index(step) >= STEPS.index(s)
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    g = remap(x["bias g"], G_OLD, G_NEW) if on("a") else x["bias g"]
    g = g * ex["eta"] if on("b") else g
    Dz = remap(x["vertical D_z"], dz_rng, (DZ_FLOOR, dz_rng[1])) if on("c") and dz_rng[0] < DZ_FLOOR else x["vertical D_z"]
    H = U.depth(V, L)
    across = (L / U.A) ** 2 / Dh <= t
    column = H * np.sqrt(Dh / Dz) > L
    creep = across & (H * (L / np.pi) / (2 * Dz) <= t)
    top = np.minimum(1.0, L / (H * np.sqrt(Dh / Dz)))
    ok = (k2 * c * t >= 10) & ~lid
    if on("d"):
        ok &= k2 * c > 2 * krac
    if on("f"):
        r_w = L / np.pi - np.minimum(H, L * np.sqrt(Dz / Dh))
        ok &= ~(column & (np.log(ex["rho"]) >= eq20 * 4 * H / r_w))
    reached = ok & np.where(column, creep, across)
    delta = F * delta_at(c, g, k2, V * (ex["r"] if on("e") else 1.0), t, np.where(column, top, 1.0))
    return reached & (delta >= 2), reached


def chains(seed):
    """All worlds' draws in the exact RNG order of the logged runs: {name: (x, published win, sweep args, salt lid)}."""
    U.RNG = np.random.default_rng(seed)
    A = {n: U.sweep(*a)[:2] + (a,) for n, a in U.WORLDS.items()}
    U.RNG = np.random.default_rng(seed)
    salt = lambda s_rng, s_star: (U.logu(*s_rng) < s_star) if s_star else np.zeros(U.N, bool)
    B = {}
    for name, (s_rng, s_star) in OW.STRAT.items():
        a = U.WORLDS[name]
        B[name] = U.sweep(*a)[:2] + (a, salt(s_rng, s_star))
    for name, (R, H, hs, t, s_rng, s_star, _) in OW.WORLDS.items():
        a = OW.args(R, H, t)
        B[name] = U.sweep(*a)[:2] + (a, salt(s_rng, s_star))
        for h in hs:  # burned exactly as ocean_worlds.main draws them
            U.sweep(*OW.args(R, h, t))
            salt(s_rng, s_star)
    return A, {**B, "Early Earth ocean": A["Early Earth ocean"] + (np.zeros(U.N, bool),)}


def run(seed):
    A, W = chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        f = lambda step, lid, **kw: outcome(x, a[1], a[4], step, ex, lid, **kw)
        assert (f("0", none)[0] == win).all(), f"{name}: all-off step is not uncertainty.sweep bit for bit"
        o = {s: tuple(100 * f(s, l)[0].mean() for l in (none, salt)) for s in STEPS}
        o["split"] = [(100 * w.mean(), 100 * (r & ~w).mean(), 100 * (1 - r.mean())) for w, r in (f("f", none), f("f", salt))]
        o["unl0"] = 100 * f("0", none)[1].mean()
        o["eq20x5"] = tuple(100 * f("f", l, eq20=5.0)[0].mean() for l in (none, salt))
        o["asp"] = tuple(100 * f("f", l, krac=K_RAC["Asp"])[0].mean() for l in (none, salt))
        out[name] = o
    for name, (x, win, a) in A.items():  # chain A regression (assumed_rules.log: win and bias unlimited)
        r = outcome(x, a[1], a[4], "0", None, np.zeros(U.N, bool))[1]
        out.setdefault("_A", {})[name] = (100 * win.mean(), 100 * r.mean())
    return out


def scale_table():
    print("\nScale rule (Eq. 13, k tau = 10): N = 2/(sqrt(pi) g^2 sqrt(10))")
    print(f"{'g':>10} {'N needed':>10} | cube side at 1 uM / 1 mM / 0.1 M (m)  | shortfall N/N_avail: 1 m3 @1 M, 1 km3 @1 mM, 10 km3 @1 mM")
    avail = [1 * 1.0 * 1000 * NA, 1e9 * 1e-3 * 1000 * NA, 1e10 * 1e-3 * 1000 * NA]
    for lab, g in (("1e-17", 1e-17), ("g_min", G_NEW[0]), ("g_max", G_NEW[1]), ("Ala", G_ALA)):
        n = (2 / (PREF * g * 10 ** 0.25)) ** 2  # review.n_needed
        assert abs(n / (2 / (np.sqrt(np.pi) * g ** 2 * np.sqrt(10))) - 1) < 1e-12
        side = [(n / (cm * 1000 * NA)) ** (1 / 3) for cm in (1e-6, 1e-3, 0.1)]
        print(f"{lab:>5} {g:.2e} {n:10.2e} | " + " / ".join(f"{s:8.3g}" for s in side) + "  | "
              + " / ".join(f"{n / v:8.3g}" for v in avail))


def main(seed=0):
    res = run(seed)
    if U.N == 200_000 and seed == 0:  # regression: all options off reproduces the logged shares within 0.3 points
        bad = [(n, res[n]["0"], v) for n, v in REG_B.items() if max(abs(np.subtract(res[n]["0"], v))) > 0.3]
        bad += [(n, res["_A"][n], v) for n, v in REG_A.items() if max(abs(np.subtract(res["_A"][n], v))) > 0.3]
        assert not bad, f"REGRESSION FAILED (name, got, logged): {bad}"
        print("regression: all-off step matches ocean_worlds/strat_sweep/assumed_rules logs within 0.3 points; bit-for-bit vs uncertainty.sweep")
    names = [n for n in res if n != "_A"]
    for k, lab in ((0, "NO LID"), (1, "FINAL (salt lid)")):
        print(f"\nLadder, {lab}, % (cumulative: 0 published, a PVED, b eta, c D_z floor, d racemization, e p/a, f depth pinning)")
        print(f"{'world':18}" + "".join(f"{s:>7}" for s in STEPS))
        for n in sorted(names, key=lambda n: -res[n]["f"][1]):
            print(f"{n:18}" + "".join(f"{res[n][s][k]:7.1f}" for s in STEPS))
    print("\nFinal step split (g), %: P_pved / P_chance / P_none [unlimited = P_pved + P_chance]; step-0 unlimited for reference")
    for n in sorted(names, key=lambda n: -res[n]["f"][1]):
        (a1, b1, c1), (a2, b2, c2) = res[n]["split"]
        assert abs(a1 + b1 + c1 - 100) < 1e-9 and abs(a2 + b2 + c2 - 100) < 1e-9
        print(f"{n:18} no lid {a1:5.1f} / {b1:5.1f} / {c1:5.1f} [{a1 + b1:5.1f}]   final {a2:5.1f} / {b2:5.1f} / {c2:5.1f} "
              f"[{a2 + b2:5.1f}]   step-0 unlimited (no lid) {res[n]['unl0']:5.1f}")
    print("\nSensitivity at step f, no lid / final, %: Eq. 20 threshold x5 | Asp k_rac")
    for n in sorted(names, key=lambda n: -res[n]["f"][1]):
        print(f"{n:18} {res[n]['eq20x5'][0]:5.1f} / {res[n]['eq20x5'][1]:5.1f} | {res[n]['asp'][0]:5.1f} / {res[n]['asp'][1]:5.1f}")
    print("\nchain A (uncertainty.py order), step 0: win / unlimited %: "
          + ", ".join(f"{n} {w:.1f} / {u:.1f}" for n, (w, u) in res["_A"].items()))
    return res


if __name__ == "__main__":
    U.N = int(sys.argv[1]) if len(sys.argv) > 1 else U.N  # smoke test: python sweep2.py 2000
    res0 = main(0)
    scale_table()
    print("\nSeed stability, step f (no lid / final), max |seed k - seed 0| over worlds, points:")
    for s in (1, 2, 3):
        r = run(s)
        d = max(abs(r[n]["f"][k] - res0[n]["f"][k]) for n in res0 if n != "_A" for k in (0, 1))
        print(f"  seed {s}: max diff {d:.2f}   " + ", ".join(f"{n} {r[n]['f'][0]:.1f}/{r[n]['f'][1]:.1f}" for n in r if n != "_A"))
