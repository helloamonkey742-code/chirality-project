"""h24_handoff_worlds: how the h24_handoff rule moves the world shares, and what the 'no layer holds' draws do. (2026-10-21)

handoff.py's world wrapper (`hand`, copied below as hand2 with one change) where the deciding layer is chosen by the
variant B that `h24_handoff.py test` picked (decision rule in h24_handoff_predictions.txt) instead of d/w >= THR.
The wrapper's layers alternate in hand, so each layer is its own domain; layer 0 holds iff d/w >= 1.25. When no layer
holds the chain is exhausted: 'none' (the draw counts as P_none, as in handoff.py's THR rule), 'deep' (the rule layer
decides) or 'top' (layer 0 decides). Columns: handoff off, THR 3, THR 1 (both handoff.py's), B none / deep / top, V1-V3
(none), B with unclipped features. Pre-registered checks W1-W6 print as a scorecard at the end. Seed 0 tables, seeds
0-3 for the spans, N = 200000, fitted h as in handoff.py.
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba --with matplotlib \
       python h24_handoff_worlds.py > h24_handoff_worlds.log      (after h24_handoff_test.log exists)
Predictions: h24_handoff_predictions.txt, sha256 090efdb7b1472b28fc1dcdda7bc21aef55e8249452600ca60a849348ceb85afb,
written Wed Oct 21 00:39:24 PDT 2026; asserted below (as are the frozen coefficients). Never edited after.
Disclosed: one SMOKE=1 run of this script (N = 2000, Europa only) before the predictions file was hashed: hand2 equalled
handoff.hand for off / THR 3 / THR 1.25 + fit, and the bisection check showed 0 mismatches in 1438 draws (so W2 was
not blind). No other run before the hash. The world shares are a deterministic function of the fitted rule, not a
simulation test of it; the out-of-sample test is `h24_handoff.py test`.
"""
import os
import re
import sys
import time

import numpy as np
from scipy.stats import spearmanr

import h24_handoff as H                 # imports handoff (seedfix wrapper applied first) and deephand_scale

HF, DH, M = H.HF, H.DH, H.M
S, S3, U = HF.S, HF.S3, HF.U
HERE = H.HERE
SMOKE = os.environ.get("SMOKE") == "1"
SEEDS = HF.SEEDS


def hand2(x, L, dz_rng, ex, lid, base, var=None, par=None, fb="none", fit=None, scale=1.0):
    """handoff.hand's lines (copied, not edited) with the deciding layer chosen by H.holds instead of d/w >= thr.
    var None: no handoff (== handoff.hand(thr=None)). var 'thr' (par = threshold) reproduces handoff.hand(thr=par).
    Layers are walked up from layer K (the remnant if K = nw) to layer 0; layer 0 holds iff d/w >= THR1. The first
    layer that holds decides. ASSUMED (checked by chain_check): H.holds is monotone in the layer index, i.e. if layer
    j holds so does every full layer above it, so the deepest holding full layer is found by bisection. If no layer
    holds, fb: 'none' (j = -1, counted P_none), 'deep' (the layer K decides), 'top' (layer 0 decides)."""
    w3, r3 = base
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    k, Hh, lr = k2 * c, U.depth(V, L), np.log(ex["rho"])
    R = L / np.pi
    column = Hh * np.sqrt(Dh / Dz) > L
    d_w = np.minimum(Hh, L * np.sqrt(Dz / Dh))
    r_w = R - d_w
    rev = column & (lr >= 4 * Hh / r_w) & (k * t >= 10) & ~lid
    rec = rev & ((Hh - d_w) * np.sqrt(Dh / Dz) > L) & (lr >= 4 * Hh / (r_w - d_w))
    with np.errstate(all="ignore"):
        rs = 4 * Hh / lr
        nw = np.ceil(Hh / d_w) - 1
        K = np.where(rec, np.maximum(np.minimum(nw, np.floor((r_w - rs) / d_w) + 1), 2), 2)
        rK1 = r_w - (K - 1) * d_w - d_w
        down = np.where(K < nw, rs / (2 * Dz) * (rs * np.log((rs - (R - Hh)) / (rs - rK1)) - (rK1 - (R - Hh))), 0.0)
        Dice = Dz * np.where(lr > 1e-9, lr / np.expm1(lr), 1.0)
        dwf = lambda j, f: f * Hh / np.sqrt(2 * Dice * np.exp(lr * (j * d_w + f * Hh / 2) / Hh) / k)
        f = np.where(rec, np.minimum(d_w, Hh - K * d_w) / Hh, 1.0)
        dw = dwf(K, f)
        j = K
        exh = np.zeros(len(K), bool)
        if var is not None:
            def S_at(ii, jj, ff):
                """Does layer jj (thickness fraction ff of the depth) hold, for the draws ii?"""
                dwm = ff * Hh[ii] / np.sqrt(2 * Dice[ii] * np.exp(lr[ii] * (jj * d_w[ii] + ff * Hh[ii] / 2) / Hh[ii]) / k[ii])
                if var == "thr":
                    return dwm >= par
                wlt = np.sqrt(2 * Dice[ii] * np.exp(lr[ii] * jj * d_w[ii] / Hh[ii]) / k[ii])
                rt = R - jj * d_w[ii]
                eps = 2 * wlt * (rt / rs[ii] - 1) / rt
                dab = d_w[ii] / np.sqrt(2 * Dice[ii] * np.exp(lr[ii] * ((jj - 1) * d_w[ii] + d_w[ii] / 2) / Hh[ii]) / k[ii])
                return np.where(jj == 0, dwm >= H.THR1, H.holds(var, par, dwm, eps, dab))
            idx = np.flatnonzero(rec)
            Ki = K[idx].astype(int)
            fl = (d_w / Hh)[idx]
            jf = np.full(len(idx), -2)
            p0 = S_at(idx, Ki, f[idx])                                     # does the rule layer K hold?
            jf[p0] = Ki[p0]
            td = np.flatnonzero(~p0)                                       # positions (in idx) still to decide
            if td.size:
                a, kk = idx[td], Ki[td]
                s1 = np.zeros(len(td), bool)
                m = kk >= 2
                s1[m] = S_at(a[m], np.ones(m.sum(), int), fl[td][m])        # does layer 1 hold?
                lo, hi = np.ones(len(td), int), kk - 1
                while True:
                    act = s1 & (lo < hi)
                    if not act.any():
                        break
                    mid = (lo + hi + 1) // 2
                    sm = np.zeros(len(td), bool)
                    sm[act] = S_at(a[act], mid[act], fl[td][act])
                    lo = np.where(act & sm, mid, lo)
                    hi = np.where(act & ~sm, mid - 1, hi)
                res = np.where(s1, lo, -2)
                z = ~s1 & (kk >= 1)                                          # layer 0 as the last candidate
                s0 = np.zeros(len(td), bool)
                s0[z] = S_at(a[z], np.zeros(z.sum(), int), fl[td][z])
                res = np.where(z & s0, 0, res)
                jf[td] = res
            ex_ = jf == -2
            jf = np.where(ex_, {"none": -1, "deep": Ki, "top": 0}[fb], jf)
            j = K.copy()
            j[idx] = jf
            exh[idx] = ex_
            f = np.where(rec & (j < K), d_w / Hh, f)
            dw = np.where(j < K, dwf(j, f), dw)
        if fit is not None:
            cf, Kmax, dlo, dhi = fit
            Xd = np.column_stack([np.ones(U.N), lr, np.log10(np.clip(j, 0, Kmax) + 1), np.log10(np.clip(dw, dlo, dhi))])
            scale = np.where(rec, 10 ** (Xd @ cf), 1.0)
        rK = r_w - (j - 1) * d_w
        up = rs / (2 * Dz) * ((R - rK) + rs * np.log((R - rs) / (rK - rs)))
        t_sel = np.maximum.reduce([10 / k, (L / U.A) ** 2 / Dh, scale * down])
        timely = rec & (j >= 0) & (scale * up + t_sel <= t)
        strong = F * HF.delta_at(c, g, k2, V * ex["r"], t, f) >= 2
    nr = max(rec.sum(), 1)
    info = {"fired": 100 * (rec & (j < K)).sum() / nr, "none": 100 * (rec & (j < 0)).sum() / nr,
            "exh": 100 * (rec & exh).sum() / nr, "dw": dw, "rec": rec, "j": j, "K": K}
    return w3 | (timely & strong), r3 | timely, info


def chain_check(x, L, dz_rng, ex, lid, var, par, n=1500):
    """Exact linear walk (layer K, then K-1 ... 0) against hand2's bisection on the first n recursion draws.
    Returns (draws checked, mismatches in the deciding layer)."""
    base = S3.outcome3(x, L, dz_rng, ex, lid)[:2]
    _, _, inf = hand2(x, L, dz_rng, ex, lid, base, var=var, par=par, fb="none")
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    k, Hh, lr = k2 * c, U.depth(V, L), np.log(ex["rho"])
    R = L / np.pi
    d_w = np.minimum(Hh, L * np.sqrt(Dz / Dh))
    with np.errstate(all="ignore"):
        rs = 4 * Hh / lr
        Dice = Dz * np.where(lr > 1e-9, lr / np.expm1(lr), 1.0)
        nw = np.ceil(Hh / d_w) - 1
    idx = np.flatnonzero(inf["rec"])[:n]
    bad = 0
    for i in idx:
        Ki = int(inf["K"][i])
        found = -1
        for jj in range(Ki, -1, -1):
            ff = min(d_w[i], Hh[i] - Ki * d_w[i]) / Hh[i] if jj == Ki else d_w[i] / Hh[i]
            with np.errstate(all="ignore"):
                dwm = ff * Hh[i] / np.sqrt(2 * Dice[i] * np.exp(lr[i] * (jj * d_w[i] + ff * Hh[i] / 2) / Hh[i]) / k[i])
                wlt = np.sqrt(2 * Dice[i] * np.exp(lr[i] * jj * d_w[i] / Hh[i]) / k[i])
                rt = R - jj * d_w[i]
                eps = 2 * wlt * (rt / rs[i] - 1) / rt
                dab = d_w[i] / np.sqrt(2 * Dice[i] * np.exp(lr[i] * ((jj - 1) * d_w[i] + d_w[i] / 2) / Hh[i]) / k[i])
            ok = dwm >= par if var == "thr" else (dwm >= H.THR1 if jj == 0 else bool(H.holds(var, par, dwm, eps, dab)))
            if ok:
                found = jj
                break
        bad += found != int(inf["j"][i])
    return len(idx), int(bad)


def logrows(path):
    """{world: [no-lid (P_pved, P_chance, P_none) per column]} from handoff.log's first table."""
    out, on = {}, False
    for line in open(os.path.join(HERE, path)):
        if line.startswith("P_pved / P_chance / P_none %, no lid"):
            on = True
        elif line.startswith("P_pved / P_chance / P_none %, final"):
            break
        elif on and "||" in line:
            name = line[:18].strip()
            out[name] = [tuple(re.findall(r"[\d.]+", c)) for c in line[18:].split("||")]
    return out


def main():
    T0 = time.time()
    if not SMOKE:
        assert H.hashlib.sha256(open(H.PRED, "rb").read()).hexdigest() == H.PRED_SHA, "predictions file changed"
    rs, cache, items = H.load_old()
    fits, _ = H.fit_variants(items)
    cf, _, _ = DH.fit(rs)
    fit = (cf, max(q["K"] for q in rs), min(q["dw"] for q in rs), max(q["dw"] for q in rs))
    if SMOKE:
        U.N = 2000
        _, W = S.chains(0)
        x, win, a, salt = W["Europa"]
        rng = np.random.default_rng([0, 1])
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        lid = np.zeros(U.N, bool)
        base = S3.outcome3(x, a[1], a[4], ex, lid)[:2]
        for kw in ({}, {"thr": 3.0}, {"thr": 1.25, "fit": fit}):
            w0, r0, _ = HF.hand(x, a[1], a[4], ex, lid, base, **kw)
            v = dict(kw)
            v = {"var": "thr", "par": v.pop("thr")} | v if "thr" in kw else v
            w1, r1, _ = hand2(x, a[1], a[4], ex, lid, base, **v)
            assert (w0 == w1).all() and (r0 == r1).all(), kw
        for var in ("V2", "V3"):
            print(var, "chain check (draws, mismatches):", chain_check(x, a[1], a[4], ex, lid, var, fits[var]))
        print("smoke ok")
        return
    best = re.search(r"best variant = (V\d)", open(H.LOG("test")).read()).group(1)
    H.frozen(fits)
    B = best
    par_of = lambda v, clip=True: None if v == "V1" else (fits[v][0], fits[v][1] if clip else None)
    print(f"h24_handoff_worlds.py, N = {U.N}, seeds {SEEDS}; seed fix applied; B = {B} (read from h24_handoff_test.log); "
          f"coefficients {[(v, np.round(b, 2).tolist()) for v, (b, _) in fits.items()]}")
    ref = logrows("handoff.log")
    cols = [("off", dict(fit=fit)), ("THR3", dict(var="thr", par=3.0, fit=fit)), ("THR1", dict(var="thr", par=1.0, fit=fit)),
            ("V1", dict(var="V1", fit=fit)), ("V2", dict(var="V2", par=par_of("V2"), fit=fit)),
            ("V3", dict(var="V3", par=par_of("V3"), fit=fit)),
            ("B deep", dict(var=B, par=par_of(B), fb="deep", fit=fit)), ("B top", dict(var=B, par=par_of(B), fb="top", fit=fit)),
            ("B unclipped", dict(var=B, par=par_of(B, False), fit=fit))]
    out, bad = {}, []
    for seed in SEEDS:
        _, W = S.chains(seed)
        for i, (name, (x, win, a, salt)) in enumerate(W.items()):
            rng = np.random.default_rng([seed, i])
            ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
            o = out.setdefault((seed, name), {})
            for li, lid in enumerate((np.zeros(U.N, bool), salt)):
                base = S3.outcome3(x, a[1], a[4], ex, lid)[:2]
                for tag, kw in cols:
                    wv, rv, inf = hand2(x, a[1], a[4], ex, lid, base, **kw)
                    o.setdefault(tag, []).append(S3.split(wv, rv))
                    o.setdefault(tag + "_info", []).append({q: inf[q] for q in ("fired", "none", "exh")})
                    if seed == 0 and li == 0 and tag in ("off", "THR3", "THR1"):
                        w0, r0, _ = HF.hand(x, a[1], a[4], ex, lid, base, **({"fit": fit} | ({"thr": kw["par"]} if tag != "off" else {})))
                        if not ((w0 == wv).all() and (r0 == rv).all()):
                            bad.append((name, tag))
                if seed == 0 and li == 0:
                    for var in ("V2", "V3"):
                        o.setdefault("chain", {})[var] = chain_check(x, a[1], a[4], ex, lid, var, fits[var])
            o["B none"], o["B none_info"] = o[B], o[B + "_info"]
            if seed == 0:
                for tag, ci in (("THR3", 4), ("THR1", 5)):
                    if name in ref and tuple(f"{v:.1f}" for v in o[tag][0]) != ref[name][ci][:3]:
                        bad.append((name, tag + " printed", tuple(f"{v:.1f}" for v in o[tag][0]), ref[name][ci][:3]))
        print(f"[seed {seed} done, {time.time() - T0:.0f} s]", flush=True)
    names = sorted({n for _, n in out}, key=lambda n: -out[(0, n)]["off"][0][0])
    print(f"\nregression: hand2 reproduces handoff.hand for handoff off / THR = 3 / THR = 1 (arrays, seed 0, no lid, "
          f"{len(names)} worlds) and handoff.log's printed THR = 3 and THR = 1 columns: {'EXACT' if not bad else 'FAILED ' + str(bad)}")
    tot = [sum(out[(0, n)]["chain"][v][0] for n in names) for v in ("V2", "V3")]
    mis = [sum(out[(0, n)]["chain"][v][1] for n in names) for v in ("V2", "V3")]
    print(f"monotone-chain check (bisection vs exact walk, first 1500 recursion draws per world, no lid, seed 0): "
          f"V2 {mis[0]} mismatches / {tot[0]}, V3 {mis[1]} / {tot[1]}")
    fmt = lambda t: f"{t[0]:5.1f} / {t[1]:4.1f} / {t[2]:5.1f}"
    for cs, lab in ((("off", "THR3", "THR1", "B none", "B deep", "B top"), ""), (("V1", "V2", "V3", "B unclipped"), " (other columns)")):
        for li, ll in ((0, "no lid"), (1, "final")):
            print(f"\nP_pved / P_chance / P_none %, {ll}, N = {U.N}, seed 0, fitted h{lab}; B = {B}; B none = chain exhausted counts as no hand")
            print(f"{'world':18} " + " || ".join(f"{c:>20}" for c in cs))
            for n in names:
                print(f"{n:18} " + " || ".join(fmt(out[(0, n)][c][li]) for c in cs))
    print("\nrecursion draws (%), no lid, seed 0: handoff fires | chain exhausted (no layer holds)")
    cs = ("THR3", "THR1", "B none", "B unclipped")
    print(f"{'world':18} " + " || ".join(f"{c:>13}" for c in cs))
    for n in names:
        print(f"{n:18} " + " || ".join(f"{out[(0, n)][c + '_info'][0]['fired']:5.1f} | {out[(0, n)][c + '_info'][0]['exh']:5.1f}".rjust(13) for c in cs))
    print("\nmovement of the no-layer draws: B deep minus B none (points), no lid | final: P_pved ; P_none")
    print(f"{'world':18} " + " || ".join(f"{c:>26}" for c in ("P_pved", "P_none")))
    M_ = {}
    for n in names:
        o = out[(0, n)]
        d = lambda k: " | ".join(f"{o['B deep'][li][k] - o['B none'][li][k]:+6.1f}" for li in (0, 1))
        M_[n] = [o["B deep"][li][0] - o["B none"][li][0] for li in (0, 1)]
        print(f"{n:18} " + " || ".join(f"{d(k):>26}" for k in (0, 2)))
    print("\nother shifts in P_pved (points), no lid | final, seed 0")
    print(f"{'world':18} " + " || ".join(f"{c:>17}" for c in ("THR1-THR3", "B deep-THR3", "B deep-off", "B top-B none", "unclipped-B none")))
    for n in names:
        o = out[(0, n)]
        d = lambda a, b: " | ".join(f"{o[a][li][0] - o[b][li][0]:+6.1f}" for li in (0, 1))
        print(f"{n:18} " + " || ".join(f"{s:>17}" for s in (d("THR1", "THR3"), d("B deep", "THR3"), d("B deep", "off"),
                                                          d("B top", "B none"), d("B unclipped", "B none"))))
    print(f"\nSeeds 0-3, P_pved % no lid / final for B deep ({B}); span = max - min over seeds")
    span = {}
    for n in names:
        v = np.array([[out[(s, n)]["B deep"][0][0], out[(s, n)]["B deep"][1][0]] for s in SEEDS])
        span[n] = np.ptp(v, axis=0).max()
        print(f"{n:18} " + " ".join(f"{p:5.1f}/{q:5.1f}" for p, q in v) + f"   span {span[n]:.2f}")
    ee = next(n for n in names if n.startswith("Early Earth"))
    rho = spearmanr([out[(0, n)]["B deep"][0][0] for n in names], [out[(0, n)]["THR1"][0][0] for n in names])[0]
    sc = [("W1", not bad, f"regression {'EXACT' if not bad else 'FAILED'}"),
          ("W2", max(mis[0] / tot[0], mis[1] / tot[1]) <= 0.005, f"bisection mismatches V2 {100 * mis[0] / tot[0]:.2f} %, V3 {100 * mis[1] / tot[1]:.2f} % (<= 0.5)"),
          ("W3", 97.0 <= out[(0, ee)]["B deep"][1][0] <= 99.9, f"{ee} final P_pved {out[(0, ee)]['B deep'][1][0]:.1f} % in [97.0, 99.9]"),
          ("W4", all(0 <= m <= 10 for v in M_.values() for m in v) and max(M_[ee]) <= 2.5,
           f"M range {min(m for v in M_.values() for m in v):+.1f} to {max(m for v in M_.values() for m in v):+.1f} (in [0, 10]); {ee} M {M_[ee][0]:+.1f} | {M_[ee][1]:+.1f} (<= 2.5)"),
          ("W5", rho >= 0.98, f"Spearman B deep vs THR = 1 across {len(names)} worlds {rho:.3f} (>= 0.98)"),
          ("W6", max(span.values()) <= 0.5, f"largest seed span {max(span.values()):.2f} points (<= 0.5)")]
    print("\n== scorecard (pass / fail rules in h24_handoff_predictions.txt) ==")
    for k, v, t in sc:
        print(f"  {k:3s} {'PASS' if v else 'FAIL':5s} {t}")
    print(f"\ntotal {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
