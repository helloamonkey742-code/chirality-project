"""sweep3: sweep2's final step (f) with reversed walls counted, not dropped (2026-10-12).

Why: sweep2.py step f counts a layered draw whose top-layer wall is past Eq. 20 (x = G/(2/r_w) >= 1) as "no single
hand" (P_none). pin_onset.log measured kappa = 1.000 and, in a smooth profile with D_z rising with depth, every run
above onset REVERSED (54/54, 0 stalled): the wall rises to the ice and the deep hand fills the ocean. So those draws
can still end with one hand: the deep layer's.

Model (sweep2.py, pin_onset.py, uncertainty.py, ocean_worlds.py, strat_sweep.py imported, NOT edited). Same draws,
seed and extra draws (eta, r, rho from default_rng([seed, world index])) as sweep2.run; every option of steps a-f on.
Geometry, read from uncertainty.sweep / sweep2.outcome: the ocean (depth H = V / area, R = L/pi) is layered
("column") if H' = H sqrt(D_h/D_z) > L. Layers are L thick in stretched units, so the top layer is
d_w = min(H, L sqrt(D_z/D_h)) deep, the top-layer wall sits at r_w = R - d_w, and the top layer holds a fraction
top = d_w/H of the molecules. G = ln(rho)/(2H) over the whole depth, x = ln(rho) r_w / (4H).
DEEP LAYER (the choice): the bulk below the top-layer wall, depth H_b = H - d_w, outer radius r_w, same L, D_h, D_z
and G. It is judged with the same per-layer rules as a whole world in uncertainty.sweep:
  - one equal-mixing layer if H_b sqrt(D_h/D_z) <= L (i.e. H_b <= d_w: the ocean held only two layers); it heals
    sideways in (L/A)^2/D_h and all of it decides: fraction f_deep = H_b/H;
  - itself layered otherwise; its own top sublayer (thickness d_w, f_deep = d_w/H) decides after its walls sink in
    H_b r_w / (2 D_z) (uncertainty.py's H R/(2 D_z)), and its own top wall, at r_w2 = r_w - d_w, is tested with
    Eq. 20 (same G, kappa = 1). Pinned -> RECURSION -> P_none (counted).
  So f_deep = min(H_b, d_w)/H.
A pinned draw ("reversed") is re-classified only if the pre-pinning checks hold: k2 a tau >= 10 and no salt lid.
  - It fills in time only if t_rise + t_sel <= tau, with t_rise = H r_w / (2 D_z (x - 1)) (pin_onset: wall speed
    D (2/r)(x - 1); H, not d_w, as the distance, so conservative) and t_sel = max(10/(k2 a), (L/A)^2/D_h, bulk
    sink time if the bulk is layered) = the deep layer's own selection time (each of its sweep2 conditions met).
  - Hand: F * Delta_deep >= 2 -> PVED-set (P_pved); one hand but F * Delta_deep < 2 -> chance-set (P_chance).
    Delta_deep = sweep2's delta_at(a, g eta, k2, V p/a, tau, f_deep): its own N, the same k2 a tau.
  - Otherwise (recursion, too late, unfinished, lid) -> P_none.
  - Racemization: the main run does NOT test k > 2 k_rac on the deep layer (as specified); sensitivity (ii) does.
  ASSUMED: D_z drawn = the profile's mean sets every time (as sweep2); a wall's sinking/rising speed is not slowed
  by the gradient below onset (as sweep2); volume fractions are depth fractions (as uncertainty.sweep's `top`).
Sensitivities: (i) t_rise x10; (ii) deep layer must also have k2 a > 2 k_rac (Ala, sweep2's rule);
  (iii, bound only) recursion counted: a draw whose second wall is also pinned is classified by the same deep-layer
  rule instead of P_none (stands in for "the next layer down fills"; NOT the headline).
Regression: recount OFF returns sweep2.outcome step f bit for bit, reproducing sweep2.log step f (and step e) for all
15 worlds within 0.05 points; this script's pinned mask, g, D_z and Delta reproduce sweep2's step e/f arrays exactly.

PREDICTIONS (written 2026-10-12, before any result-producing run; "points" are percentage points of P_pved, no lid
unless "final"):
 (P1) Regression: recount OFF equals sweep2.outcome step f bit for bit for every world (no lid and final), and every
      sweep2.log step e and f share is reproduced within 0.05 points.
 (P2) Early Earth (H' just above L is common because D_h/D_z spans only 4 decades): no lid 60-72 % (step f 39.9),
      P_chance <= 2 points.
 (P3) Every icy moon with no lid gains +1 to +6 points over sweep2 step f, except Enceladus +0.5 to +3; final
      shares gain in proportion to the no-lid share of lid-free draws; Enceladus final stays 0.0.
 (P4) Recursion dominates the moons: in every icy moon >= 80 % of reversed draws (pinned, k2 a tau >= 10, no lid)
      have a pinned second wall (x2/x1 = 1 - d_w/r_w ~ 1); Early Earth 45-65 %.
 (P5) Every world's sweep3 share (no lid and final) stays below its sweep2 step e share.
 (P6) (i) t_rise x10 lowers each world by <= 1.5 points vs the main run; Early Earth by < 0.5.
 (P7) (ii) racemization on the deep layer moves every world by <= 0.3 points.
 (P8) (iii) recursion counted (bound) >= main run everywhere; icy moons <= their step e; Early Earth >= 70 %.
 (P9) Ganymede stays the highest confirmed ocean at the final (lid) share: Ganymede > Europa > Enceladus.
 (P10) P_chance rises by <= 1.5 points over sweep2 step f's P_chance in every world (no lid and final).
 (P11) Seeds 1-3 reproduce every sweep3 share (no lid and final) within 0.5 points of seed 0.
"""
import sys

import numpy as np

import sweep2 as S
import uncertainty as U
from ocean import delta_at

# sweep2.log, seed 0, P_pved % (no lid, final)
LOG_E = {"Early Earth ocean": (74.7, 74.7), "Umbriel": (42.3, 42.3), "Ganymede": (29.9, 29.9), "Titania": (39.5, 39.5),
         "Oberon": (39.6, 39.6), "Callisto": (30.0, 30.0), "Ariel": (40.2, 40.2), "Pluto": (32.3, 32.3),
         "Titan": (25.5, 25.5), "Triton": (29.5, 29.5), "Europa": (30.4, 23.1), "Dione": (36.5, 22.8),
         "Miranda": (30.1, 10.0), "Mimas": (16.9, 5.4), "Enceladus": (24.5, 0.0)}
LOG_F = {"Early Earth ocean": (39.9, 39.9), "Umbriel": (14.4, 14.4), "Ganymede": (14.0, 14.0), "Titania": (13.7, 13.7),
         "Oberon": (13.5, 13.5), "Callisto": (13.4, 13.4), "Ariel": (11.7, 11.7), "Pluto": (9.5, 9.5),
         "Titan": (9.4, 9.4), "Triton": (8.4, 8.4), "Europa": (9.3, 7.1), "Dione": (8.9, 5.5),
         "Miranda": (12.4, 4.1), "Mimas": (5.4, 1.7), "Enceladus": (3.9, 0.0)}
CONFIRMED = ("Ganymede", "Europa", "Enceladus")


def outcome3(x, L, dz_rng, ex, lid, recount=True, rise=1.0, krac_deep=None, veto=True):
    """(P_pved mask, one-hand mask, counts) for sweep2 step f plus the reversed-wall recount."""
    w_f, r_f = S.outcome(x, L, dz_rng, "f", ex, lid)
    if not recount:
        return w_f, r_f, {}
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    k, H, lr = k2 * c, U.depth(V, L), np.log(ex["rho"])
    column = H * np.sqrt(Dh / Dz) > L
    d_w = np.minimum(H, L * np.sqrt(Dz / Dh))
    r_w = L / np.pi - d_w
    pinned = column & (lr >= 4 * H / r_w)              # sweep2's expression, so the mask is identical
    rev = pinned & (k * t >= 10) & ~lid
    if krac_deep is not None:
        rev &= k > 2 * krac_deep
    Hb = H - d_w
    col_b = Hb * np.sqrt(Dh / Dz) > L
    rec = rev & col_b & (lr >= 4 * H / (r_w - d_w))   # the deep layer's own top wall is pinned
    with np.errstate(divide="ignore"):
        t_rise = rise * H * r_w / (2 * Dz * (lr * r_w / (4 * H) - 1))
    t_sel = np.maximum.reduce([10 / k, (L / U.A) ** 2 / Dh, np.where(col_b, Hb * r_w / (2 * Dz), 0.0)])
    timely = t_rise + t_sel <= t
    deep = rev & timely & (~rec if veto else True)
    strong = F * delta_at(c, g, k2, V * ex["r"], t, np.minimum(Hb, d_w) / H) >= 2
    n = rev.sum()
    cnt = {"rev": 100 * rev.mean(), "rec": 100 * rec.mean(), "rec/rev": 100 * rec.sum() / max(n, 1),
           "late": 100 * (rev & ~rec & ~timely).mean(), "chance": 100 * (deep & ~strong).mean()}
    return w_f | (deep & strong), r_f | deep, cnt


def check(x, L, dz_rng, ex, lid):
    """The recount's own pinned mask, g, D_z and Delta reproduce sweep2's step e and f arrays exactly."""
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    H = U.depth(V, L)
    column = H * np.sqrt(Dh / Dz) > L
    pinned = column & (np.log(ex["rho"]) >= 4 * H / (L / np.pi - np.minimum(H, L * np.sqrt(Dz / Dh))))
    w_e, r_e = S.outcome(x, L, dz_rng, "e", ex, lid)
    w_f, r_f = S.outcome(x, L, dz_rng, "f", ex, lid)
    assert (w_f == w_e & ~pinned).all() and (r_f == r_e & ~pinned).all(), "pinned mask differs from sweep2"
    top = np.where(column, np.minimum(1.0, L / (H * np.sqrt(Dh / Dz))), 1.0)
    assert (w_e == r_e & (F * delta_at(c, g, k2, V * ex["r"], t, top) >= 2)).all(), "g / D_z / Delta differ from sweep2"
    return w_e


def split(w, r):
    return 100 * w.mean(), 100 * (r & ~w).mean(), 100 * (1 - r.mean())


def run(seed, full=True):
    _, W = S.chains(seed)
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, i])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        f = lambda lid, **kw: outcome3(x, a[1], a[4], ex, lid, **kw)
        o = {"main": [split(*f(l)[:2]) for l in (none, salt)]}
        if full:
            o["cnt"] = f(none)[2]
            for lid_i, l in enumerate((none, salt)):
                w_off, r_off, _ = f(l, recount=False)
                w_f, r_f = S.outcome(x, a[1], a[4], "f", ex, l)
                assert (w_off == w_f).all() and (r_off == r_f).all(), f"{name}: recount OFF is not sweep2 step f"
                o.setdefault("e", []).append(100 * check(x, a[1], a[4], ex, l).mean())
                o.setdefault("f", []).append(split(w_f, r_f))
            o["rise10"] = [split(*f(l, rise=10.0)[:2]) for l in (none, salt)]
            o["rac"] = [split(*f(l, krac_deep=S.K_RAC["Ala"])[:2]) for l in (none, salt)]
            o["noveto"] = [split(*f(l, veto=False)[:2]) for l in (none, salt)]
        out[name] = o
    return out


def main():
    res = run(0)
    if U.N == 200_000:
        got = {n: {"e": tuple(res[n]["e"]), "f": tuple(v[0] for v in res[n]["f"])} for n in res}
        bad = [(n, k, got[n][k], log[n]) for n in res for k, log in (("e", LOG_E), ("f", LOG_F))
               if max(abs(np.subtract(got[n][k], log[n]))) > 0.05]
        assert not bad, f"REGRESSION FAILED (world, step, got, sweep2.log): {bad}"
        print("regression: recount OFF = sweep2.outcome step f bit for bit (all worlds, no lid and final); "
              "pinned mask, g, D_z, Delta reproduce sweep2 step e/f arrays exactly; sweep2.log steps e and f within 0.05")
    order = sorted(res, key=lambda n: -res[n]["main"][1][0])
    print(f"\nP_pved %, no lid / final: sweep2 step e | sweep2 step f | sweep3 (reversed walls counted)   N={U.N}, seed 0")
    print(f"{'world':18} {'step e':>13} {'step f':>13} {'sweep3':>13} {'gain':>11}")
    for n in order:
        e, fs, m = res[n]["e"], res[n]["f"], res[n]["main"]
        print(f"{n:18} {e[0]:5.1f} / {e[1]:5.1f} {fs[0][0]:5.1f} / {fs[1][0]:5.1f} {m[0][0]:5.1f} / {m[1][0]:5.1f} "
              f"{m[0][0] - fs[0][0]:+5.1f}/{m[1][0] - fs[1][0]:+5.1f}")
    print("\nThree-way split, %: P_pved / P_chance / P_none, no lid | final;  sweep2 step f split for reference")
    for n in order:
        (a1, b1, c1), (a2, b2, c2) = res[n]["main"]
        (d1, e1, g1), (d2, e2, g2) = res[n]["f"]
        assert abs(a1 + b1 + c1 - 100) < 1e-9 and abs(a2 + b2 + c2 - 100) < 1e-9
        print(f"{n:18} {a1:5.1f} / {b1:4.1f} / {c1:5.1f} | {a2:5.1f} / {b2:4.1f} / {c2:5.1f}   "
              f"[step f {d1:5.1f} / {e1:4.1f} / {g1:5.1f} | {d2:5.1f} / {e2:4.1f} / {g2:5.1f}]")
    print("\nReversed draws (no lid), % of all draws: reversed (pinned, k2 a tau >= 10) | recursion (second wall pinned) "
          "| % of reversed that recurse | too late | deep hand chance-set")
    for n in order:
        q = res[n]["cnt"]
        print(f"{n:18} {q['rev']:6.1f} | {q['rec']:6.1f} | {q['rec/rev']:6.1f} | {q['late']:5.1f} | {q['chance']:5.1f}")
    print("\nSensitivities, P_pved %, no lid / final: main | (i) t_rise x10 | (ii) deep layer k > 2 k_rac (Ala) "
          "| (iii) recursion counted (bound) | sweep2 step e")
    for n in order:
        s = [res[n][k] for k in ("main", "rise10", "rac", "noveto")]
        print(f"{n:18} " + " | ".join(f"{v[0][0]:5.1f} / {v[1][0]:5.1f}" for v in s)
              + f" | {res[n]['e'][0]:5.1f} / {res[n]['e'][1]:5.1f}")
    return res


if __name__ == "__main__":
    U.N = int(sys.argv[1]) if len(sys.argv) > 1 else U.N  # smoke test: python sweep3.py 2000
    res0 = main()
    print("\nSeed stability, sweep3 P_pved / P_chance / P_none (no lid and final), max |seed k - seed 0| over worlds, points:")
    for s in (1, 2, 3):
        r = run(s, full=False)
        d = max(abs(r[n]["main"][l][j] - res0[n]["main"][l][j]) for n in res0 for l in (0, 1) for j in range(3))
        print(f"  seed {s}: max diff {d:.2f}   " + ", ".join(f"{n} {r[n]['main'][0][0]:.1f}/{r[n]['main'][1][0]:.1f}" for n in r))
