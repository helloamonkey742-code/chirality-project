"""Extra experiment for the NaClO3 droplet test, POST HOC: model study plus a re-fit of published points, not new data.

1. Can drying drops make a critical size without counting?  The one mechanism we can simulate is exclusion: each
   nucleus depletes a zone ell(dt) = sqrt(dt) + nu dt around it (diffusion, plus a drift nu) and later nuclei inside a
   zone are blocked (Johnson-Mehl thinning).  The nucleation rate e-folds every t_J as the drop concentrates, so
   after the first nucleus (t = 0) candidates arrive with cumulative hazard G1 e^t (units t_J, ell_J).  Sites: the
   drop's volume (hemisphere), its base or its contact ring.  Single = 2^(1 - N); a drop still not covered after 40000
   candidates counts as mixed (flagged).  One control number, eta = drop radius / ell_J, gives P(eta) by Monte Carlo.
   Over volume eta ~ V^e, e = 1/3 - gamma/2 (3-D drops) or 1/2 - gamma (films), drying time ~ V^gamma: e <= 1/2.
2. The frozen droplet_sim2 rule (c* below, from droplet_sim2.log) on those curves, and on Bak et al. (2023) Fig. 2a,
   plus the exclusion e-profile on the Tinc 90 C series.  Bak's data were read before this was written, so every
   result here is POST HOC; expectations E1-E4 are in droplet_sim2_predictions.txt (sha256 in droplet_sim2.PRED_SHA).

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python exclusion_bak.py > exclusion_bak.log
"""
from multiprocessing import Pool

import numba as nb
import numpy as np

from blind_units import BAK_MAIN, BAK_TINC90, exp_fit
from droplet_power import fit_h, scorable
from droplet_sim import LADDER
from droplet_sim2 import EPS_FIT, LF, V50, classify, curves, run_cell, stats

C_STAR = -1.61                                # droplet_sim2.log (calibration seed), frozen
ETA = np.geomspace(0.05, 6.0, 34)
NU = (0.0, 0.3, 1.0)
SITES = ("volume", "base", "ring")
DROPS, CAP, NCAP = 6000, 20, 40000


@nb.njit(cache=True)
def one(eta, nu, site):
    """Accepted nuclei in one drop of radius eta, and 1 if it overflowed."""
    E1 = np.random.exponential(1.0)
    H = E1
    xs = np.zeros((CAP, 3))
    ts = np.zeros(CAP)
    n = 0
    for it in range(NCAP):
        if it > 0:
            H += np.random.exponential(1.0)
        t = np.log(H / E1)
        if it > 0 and np.sqrt(t) + nu * t > 2 * eta:          # the first zone covers the drop
            return n, 0
        if site == 0:                                         # uniform in a hemisphere
            r = eta * np.random.random() ** (1 / 3)
            z = np.random.random()
            ph = 2 * np.pi * np.random.random()
            s = np.sqrt(1 - z * z)
            x0, x1, x2 = r * s * np.cos(ph), r * s * np.sin(ph), r * z
        elif site == 1:                                       # uniform on the base disc
            r = eta * np.sqrt(np.random.random())
            ph = 2 * np.pi * np.random.random()
            x0, x1, x2 = r * np.cos(ph), r * np.sin(ph), 0.0
        else:                                                 # contact ring
            ph = 2 * np.pi * np.random.random()
            x0, x1, x2 = eta * np.cos(ph), eta * np.sin(ph), 0.0
        ok = True
        for k in range(n):
            dt = t - ts[k]
            if (x0 - xs[k, 0]) ** 2 + (x1 - xs[k, 1]) ** 2 + (x2 - xs[k, 2]) ** 2 < (np.sqrt(dt) + nu * dt) ** 2:
                ok = False
                break
        if ok:
            xs[n, 0], xs[n, 1], xs[n, 2] = x0, x1, x2
            ts[n] = t
            n += 1
            if n >= CAP:
                return n, 0
    return n, 1


@nb.njit(cache=True)
def table_row(nu, site, seed, etas, drops):
    np.random.seed(seed)
    P = np.zeros(len(etas))
    ov = np.zeros(len(etas))
    for j in range(len(etas)):
        for _ in range(drops):
            n, o = one(etas[j], nu, site)
            P[j] += 0.0 if o else 2.0 ** (1 - n)
            ov[j] += o
    return P / drops, ov / drops


def table(cfg):
    i, s = cfg
    return table_row(NU[i], s, 1000 + 10 * i + s, ETA, DROPS)


def curve(T, e):
    """P on LF for eta = eta0 (V / V50)^e, eta0 where T crosses 1/2 (first crossing), or None if T never does."""
    j = np.flatnonzero((T[:-1] >= 0.5) & (T[1:] < 0.5))
    if not j.size:
        return None
    j = j[0]
    l0 = np.log(ETA[j]) + (T[j] - 0.5) / (T[j] - T[j + 1]) * np.log(ETA[j + 1] / ETA[j])
    return np.interp(l0 + e * (LF - np.log(V50)), np.log(ETA), T)


def true_h(P):
    at = lambda v: np.interp(np.log(v), LF, P)
    lev = exp_fit(LADDER, at(LADDER)) * np.geomspace(0.2, 1.6, 8)
    return fit_h(lev, at(lev), np.full(8, 1e5), 0.0)[1]


def bak(series, n1):
    """ln V, k, n with n1 drops at the markers whose n is not printed (k = n there)."""
    V = np.array([r[0] for r in series])
    P = np.array([r[1] for r in series])
    n = np.array([n1 if len(r) < 3 and r[1] == 1.0 else r[2] if len(r) == 3 and r[2] else n1 for r in series], float)
    if series is BAK_TINC90:                                  # printed fractions 10/12, 5/6, 2/6, 5/16
        n[7:] = (12, 6, 6, 16)
    return np.log(V), np.round(P * n), n


def ll_masked(K, N, cur, fam, keep):
    _, lp, l1, par = cur[fam]
    return (K @ lp + (N - K) @ l1)[:, keep(par)].max(1)


def self_check():
    a = table_row(0.0, 0, 7, np.array([0.05, 6.0]), 400)
    b = table_row(0.0, 0, 7, np.array([0.05, 6.0]), 400)
    assert (a[0] == b[0]).all(), "seeded Monte Carlo must repeat"
    assert a[0][0] > 0.95 and a[0][1] < 0.05, a
    lv, k, n = bak(BAK_TINC90, 3)
    assert n.tolist() == [3] * 7 + [12, 6, 6, 16] and k.tolist() == [3] * 7 + [10, 5, 2, 5]
    lv, k, n = bak(BAK_MAIN, 6)
    assert n.tolist() == [6, 10, 10] and k.tolist() == [6, 8, 7]


def main():
    assert C_STAR is not None
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    cfgs = [(i, s) for i in range(len(NU)) for s in range(len(SITES))]
    with Pool() as pool:
        tabs = dict(zip(cfgs, pool.map(table, cfgs)))
    print(f"POST HOC.  Model study and re-fit of published points, not new data.  Frozen c* = {C_STAR}.")
    print(f"\n1. Exclusion Monte Carlo, P_single against eta (drop radius / zone scale), {DROPS} drops per point:")
    for (i, s), (T, ov) in tabs.items():
        pick = [0, 8, 13, 17, 21, 25, 29, 33]
        print(f"  {SITES[s]:6} nu {NU[i]:<3} " + " ".join(f"{ETA[j]:.2g}:{T[j]:.2f}" for j in pick)
              + f"  (overflow max {ov.max():.3f})")

    E = (-1 / 6, 0.0, 1 / 6, 1 / 3, 1 / 2)
    print("\n2. Exclusion curves over volume, eta ~ V^e (e = 1/3: 3-D drop, drying time fixed; 1/6: gamma 1/3; 0:")
    print("   gamma 2/3, diffusion-limited; 1/2: film).  True h at 0.2-1.6 x V0, then the frozen rule, CV 0.05,")
    print("   misread 0.03, 500 runs:")
    hmax, stepmax, rows = 0.0, 0.0, []
    with Pool() as pool:
        for (i, s), (T, _) in tabs.items():
            for e in E:
                P = curve(T, e)
                if P is not None:
                    rows.append(((i, s, e), P, pool.apply_async(run_cell, ((P, 0 * P, 0.05, 0.03, 500, [3, i, s, E.index(e)], False),))))
        for (i, s, e), P, job in rows:
            lab = classify(job.get(), C_STAR)
            h = true_h(P)
            if e > 0:
                hmax, stepmax = max(hmax, h), max(stepmax, np.mean(lab == "STEP"))
            print(f"  {SITES[s]:6} nu {NU[i]:<3} e {e:+.2f}: true h {h:5.2f}  "
                  + " ".join(f"{v} {np.mean(lab == v):.2f}" for v in ("STEP", "SMOOTH", "FLAT", "RISING", "NS")))
    missing = [f"{SITES[s]} nu {NU[i]}" for (i, s), (T, _) in tabs.items() if curve(T, 0.0) is None]
    if missing:
        print("  P never crosses 1/2 (no curve): " + ", ".join(missing))

    print("\n3. Bak et al. 2023 Fig. 2a, frozen rule.  n at the P = 1 markers is not printed; tried at 3, 6, 10.")
    print("   L = ll COUNT - ll STEP (s <= 0.6); p = parametric bootstrap share of L <= observed under the best COUNT fit;")
    print("   sharp = 2 (ll best - ll best STEP with s <= 0.15).")
    E3, E4 = [], []
    EG = np.round(np.arange(-0.2, 1.01, 0.05), 2)
    REF = np.linspace(np.log(0.02), np.log(8), 160)
    rng = np.random.default_rng(4)
    for name, series in (("Tinc 90 C", BAK_TINC90), ("main", BAK_MAIN)):
        for n1 in (3, 6, 10):
            x, k, n = bak(series, n1)
            cur = curves(x)
            K, N = k[None], n[None]
            r = stats(K, N, x, cur)
            r["ns"] = np.array([not scorable(k, n)])
            lab = classify(r, C_STAR)[0]
            llc = ll_masked(K, N, cur, "count", lambda p: p > 0)[0]
            lls = ll_masked(K, N, cur, "step", lambda p: p > 0)[0]
            sharp = 2 * (max(llc, lls) - ll_masked(K, N, cur, "step", lambda p: p <= 0.15)[0])
            Pc = cur["count"][0][np.argmax(K @ cur["count"][1] + (N - K) @ cur["count"][2])]
            Ks = rng.binomial(N.astype(int), EPS_FIT + (1 - 2 * EPS_FIT) * Pc, (2000, len(x))).astype(float)
            p = np.mean(stats(Ks, np.repeat(N, 2000, 0), x, cur)["L"] <= r["L"][0])
            line = (f"  {name:9} n={n1:<2}: L {r['L'][0]:+6.2f}  p {p:.3f}  a-hat {r['a'][0]:.2f}  s-hat {r['s'][0]:.3f}  "
                    f"sharp {sharp:5.1f}  z {r['z'][0]:+.2f}  verdict {lab}")
            if series is BAK_TINC90:
                E3.append(lab != "STEP" and sharp > 3.84)
                best = {}
                for (i, s), (T, _) in tabs.items():
                    for e in EG:
                        M = np.interp(REF[:, None] + e * (x - np.log(V50)), np.log(ETA), T)
                        M = np.clip(EPS_FIT + (1 - 2 * EPS_FIT) * M, 1e-12, 1 - 1e-12)
                        best[(i, s, e)] = (k * np.log(M) + (n - k) * np.log1p(-M)).sum(1).max()
                prof = {e: max(v for key, v in best.items() if key[2] == e) for e in EG}
                top = max(prof.values())
                eb = max(prof, key=prof.get)
                E4.append(eb > 0.5 and prof[0.0] - top < -2)
                line += (f"\n      exclusion e-profile (ll - best): "
                         + " ".join(f"{e:g}:{prof[e] - top:.1f}" for e in EG[::2])
                         + f"\n      best e {eb:g}; ll EXCL {top:.2f} vs COUNT {llc:.2f}, STEP {lls:.2f}")
            print(line)

    print("\nExpectations (POST HOC, not blind):")
    for name, held, detail in (
            ("E1", hmax <= 3, f"max true h over physical exclusion curves (0 < e <= 1/2) {hmax:.2f}"),
            ("E2", stepmax <= 0.07, f"max STEP share {stepmax:.2f}"),
            ("E3", all(E3), f"Tinc 90 C not STEP and sharp step rejected at n = 3, 6, 10: {E3}"),
            ("E4", all(E4), f"best e > 1/2 and e = 0 rejected by > 2 at n = 3, 6, 10: {E4}")):
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")


if __name__ == "__main__":
    main()
