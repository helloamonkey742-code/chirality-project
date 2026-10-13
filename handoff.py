"""handoff: replace deephand_scale's "ocean-scale floor" with an estimate. (2026-10-12)

Why: deephand_scale (section 1, 600 simulated stacks) found the rule "the hand of the deciding layer fills the
shell" held 365/365 when that layer was >= 3 local wall widths thick, and when it failed the layer was thin and the
hand of the layer above filled the shell (journals/new_analyses_deephand.md, Section 3). Its ocean-scale variant
counted every failed thin-remnant draw as P_none, so its 89.0 % (Early Earth) is a floor, not an estimate.

Method (a wrapper: nothing is edited, multiwall / sweep2 / sweep3 / deephand_scale are imported).
 * Seed fix: np.random.default_rng is wrapped exactly as seedfix.py does (the seed [s, 0] becomes [s, 0, 1]) before
   any project module is imported, so Enceladus (world 0) gets its own extra draws; worlds 1-14 are unchanged.
 * Handoff rule in multiwall.multi's wrapper (multi's lines copied into hand(), not edited). Layers are d_w thick
   from the ice down (layer j between depths j d_w and (j+1) d_w), the bottom remnant H - nw d_w; K = walls past
   onset, so layer K (the remnant when K = nw) decides. Its d/w = thickness / local wall width sqrt(2 D_loc / k)
   (deephand_scale.geom's definition). If d/w < THR (= 3) the deciding role passes to layer K - 1 (a full layer),
   then K - 2, ..., until a layer is >= THR wall widths thick; if none is (layer 0 too thin), the draw has no hand
   (counted P_none). Full layers keep their thickness and the wall width grows with depth, so the deepest full layer
   that passes has a closed form, j_max = floor(H / (ln rho d_w) ln(d_w^2 k / (2 THR^2 D_ice)) - 1/2).
   The new deciding layer j gets its own PVED-set test, F Delta >= 2 with Delta on its own N (f = d_w / H), and its
   own timing: T_up of the deepest rising wall above it (wall j, multiwall's closed form with K -> j; 0 when j = 0)
   and multiwall's T_down unchanged (ASSUMED: the thin layers below the new deciding layer die faster than any wall
   moves - in the simulations they died before a wall moved; the sink term for onset inside the stack is kept).
   Timing is run twice: h = 1 (multiwall's closed form) and deephand_scale's fitted h (log10 h = c0 + c1 ln rho +
   c2 log10(K + 1) + c3 log10(d/w), coefficients refit from deephand_scale.log's ROW lines by deephand_scale.fit,
   clipped to the simulated K and d/w ranges as there), evaluated with the NEW deciding layer's j and d/w.
 * Secondary (also pre-registered): THR = 1 (hand off only below one wall width), and an exploratory threshold scan
   of the validation (post hoc, labelled so).
 * Validation on deephand_scale's 600 ROW records: each run's layer bounds are rebuilt from its (R, n, ln rho, u)
   with deephand_scale.config, layers of one hand merged into domains (walls exist only between hands; the ROW's d/w
   is a domain's), local wall width sqrt(2 D) with D = (w^2/2) exp(ln rho (R - r)) (multiwall.shell's profile with
   G = 2/r*). The rule domain (containing r*) is handed up the same way; prediction = hand of the first domain >= THR
   wall widths (none -> "no hand", counted a miss). Agreement with the observed final hand +- SE (multiwall.se),
   overall, on runs where the handoff fired, and by w. THR = 0 must reproduce the ROW "ok" flag run by run.
   Pre-registered validation threshold: a handoff rule is VALIDATED if overall agreement >= 95 % and agreement on
   fired runs >= 85 %.
 * Worlds: seeds 0-3, N = 200000, no lid | final. Beside: multiwall main and x10 (seed fix applied; = seedfix's
   values), the ocean-scale floor (recomputed with the seed fix by deephand_scale's own formula).
 * Regressions: (A) handoff off, h = 1: arrays equal multiwall.multi and splits equal seedfix.log's multiwall main
   (worlds 1-14 from multiwall.log section 4, Enceladus from seedfix.log's corrected line) to the printed 0.1, all
   15 worlds; (B) handoff off, fitted h: equals deephand_scale.log's fitted-h column, worlds 1-14; (C) the
   recomputed floor equals deephand_scale.log's ocean-scale column, worlds 1-14.
  ASSUMED: the wrapper's layers alternate in hand (as sweep3/multiwall), so "the layer above" is the other hand and
  a handoff of two layers returns to the first hand; the simulated stacks (w >= 0.0025 of the depth, K <= 9) stand in
  for the moons' far thinner walls and larger K (as deephand_scale assumed).

PREDICTIONS (written 2026-10-12, before any result-producing run):
 (P1) Validation, THR = 3: overall agreement lies in [74, 84] % and is below the plain rule's (no handoff) 486/600
      = 81.0 %; on runs where the handoff fired it lies in [40, 58] %. (In deephand_scale the thin domain itself
      still filled in 97 % of runs at d/w 2-3 and 61 % at 1-2, and a 3-width handoff overrules those.)
 (P2) Validation, THR = 3 is NOT validated (fails the pre-registered >= 95 % overall / >= 85 % fired threshold).
 (P3) Validation, THR = 1: overall agreement in [89, 96] %, fired agreement >= 80 %; it fails the >= 95 % overall
      part of the threshold.
 (P4) Validation: THR = 0 reproduces the ROW "ok" flag in all 600 runs; the rebuilt rule-domain d/w matches the
      ROW's d/w in all 600 (to rounding).
 (P5) Regression (A) EXACT for all 15 worlds (arrays and printed 0.1).
 (P6) Regression (B) EXACT for worlds 1-14; regression (C) EXACT for worlds 1-14.
 (P7) Worlds, THR = 3, seed 0, no lid: the handoff fires in 25-65 % of recursion draws on every world.
 (P8) Worlds, THR = 3: no layer >= 3 wall widths ("no hand") in <= 5 % of recursion draws on every world.
 (P9) Worlds, THR = 3, fitted h, seed 0: every world's P_pved (no lid and final) is >= the ocean-scale floor and
      <= multiwall main (seed fix) + 3.
 (P10) Early Earth, THR = 3: P_pved >= 97 % with fitted h and with h = 1.
 (P11) THR = 3: fitted h and h = 1 differ by <= 2 points in P_pved on every world (no lid and final, seed 0).
 (P12) THR = 3, fitted h: P_pved over seeds 0-3 spans <= 0.5 points on every world (no lid and final).
 (P13) THR = 1 (fitted h) is within 3 points of THR = 3 (fitted h) in P_pved on every world (no lid and final).
"""
import os
import re
import sys
import time

import numpy as np

ORIG_RNG = np.random.default_rng


def fixed_rng(seed=None, *a, **kw):
    """seedfix.py's wrapper (copied): the exact seed [s, 0] becomes [s, 0, 1]."""
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng   # before the project modules import

import deephand_scale as DH   # noqa: E402
import multiwall as M         # noqa: E402
import sweep2 as S            # noqa: E402
import sweep3 as S3           # noqa: E402
import uncertainty as U       # noqa: E402
from ocean import delta_at    # noqa: E402

SMOKE = os.environ.get("SMOKE") == "1"
HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = (0, 1, 2, 3)


def hand(x, L, dz_rng, ex, lid, base, thr=None, fit=None, scale=1.0):
    """multiwall.multi's lines (copied) with the deciding layer handed up while thinner than thr wall widths.
    thr=None: no handoff (== multiwall.multi). fit=(c, Kmax, dlo, dhi): deephand_scale's fitted h replaces scale.
    Returns (P_pved mask, one-hand mask, info)."""
    w3, r3 = base
    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    Dz = S.remap(x["vertical D_z"], dz_rng, (S.DZ_FLOOR, dz_rng[1])) if dz_rng[0] < S.DZ_FLOOR else x["vertical D_z"]
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    k, H, lr = k2 * c, U.depth(V, L), np.log(ex["rho"])
    R = L / np.pi
    column = H * np.sqrt(Dh / Dz) > L
    d_w = np.minimum(H, L * np.sqrt(Dz / Dh))
    r_w = R - d_w
    rev = column & (lr >= 4 * H / r_w) & (k * t >= 10) & ~lid
    rec = rev & ((H - d_w) * np.sqrt(Dh / Dz) > L) & (lr >= 4 * H / (r_w - d_w))
    with np.errstate(all="ignore"):
        rs = 4 * H / lr
        nw = np.ceil(H / d_w) - 1
        K = np.where(rec, np.maximum(np.minimum(nw, np.floor((r_w - rs) / d_w) + 1), 2), 2)
        rK1 = r_w - (K - 1) * d_w - d_w
        down = np.where(K < nw, rs / (2 * Dz) * (rs * np.log((rs - (R - H)) / (rs - rK1)) - (rK1 - (R - H))), 0.0)
        Dice = Dz * np.where(lr > 1e-9, lr / np.expm1(lr), 1.0)               # D_z at the ice (deephand geom)
        dwf = lambda j, f: f * H / np.sqrt(2 * Dice * np.exp(lr * (j * d_w + f * H / 2) / H) / k)
        f = np.where(rec, np.minimum(d_w, H - K * d_w) / H, 1.0)
        dw = dwf(K, f)
        j = K
        if thr is not None:
            jmax = np.floor(H / (lr * d_w) * np.log(d_w ** 2 * k / (2 * thr ** 2 * Dice)) - 0.5)
            j = np.where(dw >= thr, K, np.minimum(K - 1, jmax))
            f = np.where(rec & (j < K), d_w / H, f)
            dw = np.where(j < K, dwf(j, f), dw)
        if fit is not None:
            cf, Kmax, dlo, dhi = fit
            Xd = np.column_stack([np.ones(U.N), lr, np.log10(np.clip(j, 0, Kmax) + 1), np.log10(np.clip(dw, dlo, dhi))])
            scale = np.where(rec, 10 ** (Xd @ cf), 1.0)
        rK = r_w - (j - 1) * d_w
        up = rs / (2 * Dz) * ((R - rK) + rs * np.log((R - rs) / (rK - rs)))
        t_sel = np.maximum.reduce([10 / k, (L / U.A) ** 2 / Dh, scale * down])
        timely = rec & (j >= 0) & (scale * up + t_sel <= t)
        strong = F * delta_at(c, g, k2, V * ex["r"], t, f) >= 2
    nr = max(rec.sum(), 1)
    info = {"fired": 100 * (rec & (j < K)).sum() / nr, "none": 100 * (rec & (j < 0)).sum() / nr,
            "steps": float(np.median((K - j)[rec & (j < K) & (j >= 0)])) if (rec & (j < K) & (j >= 0)).any() else np.nan,
            "dw": dw, "rec": rec, "j": j, "K": K, "dwf": dwf, "d_w": d_w, "H": H}
    return w3 | (timely & strong), r3 | timely, info


def validate(rs, thr):
    """Per ROW run: (predicted hand == observed final, handoff fired, no hand). Domains = runs of one hand."""
    out = []
    for q in rs:
        R, n, lr, w = q["R"], q["n"], q["lr"], q["w"]
        _, _, _, rstar, _, b, _ = DH.config(0, None, lr=lr, u=q["u"], R=R, n=n, alt=True)
        hands = [1 if s == "+" else -1 for s in q["hands"]]
        e = np.concatenate([[R], b, [R - 1.0]])
        cuts = [0] + [m + 1 for m in range(n - 1) if hands[m] != hands[m + 1]] + [n]
        dws = [(e[a] - e[z]) / np.sqrt(w * w * np.exp(lr * (R - (e[a] + e[z]) / 2))) for a, z in zip(cuts[:-1], cuts[1:])]
        jl = int(np.sum(b > rstar))
        d0 = d = next(m for m, (a, z) in enumerate(zip(cuts[:-1], cuts[1:])) if a <= jl < z)
        assert abs(dws[d0] - q["dw"]) <= 1e-3 * max(q["dw"], 1.0) + 1e-3, (q, dws[d0])
        while d >= 0 and dws[d] < thr:
            d -= 1
        pred = hands[cuts[d]] if d >= 0 else 0
        out.append((pred == q["fin"], d != d0, pred == 0))
    return np.array(out, bool)


def pct(v):
    return "n=%3d  %5.1f +- %4.1f %%" % (len(v), *(100 * np.array(M.se(v)))) if len(v) else "n=  0"


def section_validate(rs):
    print(f"\n=== Validation on deephand_scale.log's {len(rs)} ROW runs ===")
    ok = np.array([q["ok"] for q in rs])
    v0 = validate(rs, 0.0)
    print(f"THR = 0 reproduces the ROW 'ok' flag run by run: {bool((v0[:, 0] == ok).all())}; plain rule {pct(ok)}")
    ws, dwr = np.array([q["w"] for q in rs]), np.array([q["dw"] for q in rs])
    res = {}
    for thr in (3.0, 1.0):
        v = validate(rs, thr)
        a, fired, none = v[:, 0], v[:, 1], v[:, 2]
        res[thr] = (M.se(a)[0], M.se(a[fired])[0] if fired.any() else np.nan)
        val = res[thr][0] >= 0.95 and res[thr][1] >= 0.85
        print(f"THR = {thr:g}: agreement {pct(a)}; handoff fired: {pct(a[fired])}; not fired: {pct(a[~fired])}; "
              f"no hand predicted in {none.sum()} runs  -> {'VALIDATED' if val else 'NOT validated'} (>= 95 / >= 85)")
        print("  by w: " + "; ".join(f"{w}: {pct(a[ws == w])}" for w in DH.WS))
        print("  fired, by the rule domain's d/w: " + "; ".join(
            f"[{lo:g},{hi:g}) {pct(a[fired & (dwr >= lo) & (dwr < hi)])}" for lo, hi in ((0, 1), (1, 2), (2, 3))))
    print("exploratory (post hoc) threshold scan, overall agreement | fired:")
    for thr in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 5.0):
        v = validate(rs, thr)
        print(f"  THR = {thr:4g}: {100 * v[:, 0].mean():5.1f} % | {pct(v[v[:, 1], 0])}")
    return res


def logtab(path, start=None):
    """{world: tuple of every printed number after the first '||'} from the first matching table lines."""
    pat = re.compile(r"^(?:\s+new: )?(\S.*?)\s+[\d.]+ /\s+[\d.]+ /\s+[\d.]+ \|\s+[\d.]+ /\s+[\d.]+ /\s+[\d.]+\s+\|\|\s+(.*)$")
    out, on = {}, start is None
    for line in open(os.path.join(HERE, path)):
        on = on or (start in line)
        m = pat.match(line) if on else None
        if m and m.group(1).strip() not in out:
            out[m.group(1).strip()] = [tuple(re.findall(r"[\d.]+", p)) for p in m.group(2).split("||")]
    return out


def p3(s):
    return tuple(f"{v:.1f}" for t in s for v in t)


def worlds(fit, pb):
    mw = logtab("multiwall.log")
    sf = logtab("seedfix.log", "=== D. multiwall Section 4")
    dh = logtab("deephand_scale.log", "P_pved / P_chance / P_none %, no lid | final")
    refA = {n: v[0] for n, v in mw.items()}
    refA["Enceladus"] = sf["Enceladus"][0]
    out, bad = {}, []
    for seed in SEEDS:
        _, W = S.chains(seed)
        for i, (name, (x, win, a, salt)) in enumerate(W.items()):
            rng = np.random.default_rng([seed, i])
            ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
            o = out.setdefault((seed, name), {})
            for lid in (np.zeros(U.N, bool), salt):
                base = S3.outcome3(x, a[1], a[4], ex, lid)[:2]
                hf = lambda **kw: hand(x, a[1], a[4], ex, lid, base, **kw)
                runs = {"main": hf(), "x10": hf(scale=10.0), "fit0": hf(fit=fit), "h1": hf(thr=3.0),
                        "fit": hf(thr=3.0, fit=fit), "thr1": hf(thr=1.0, fit=fit)}
                if seed == 0:
                    wm, rm, _ = M.multi(x, a[1], a[4], ex, lid)
                    if not ((wm == runs["main"][0]).all() and (rm == runs["main"][1]).all()):
                        bad.append((name, "A arrays"))
                    w3, r3 = base
                    wh, rh, inf = runs["fit0"]
                    dw = inf["dw"]
                    p = np.where(dw >= 1, pb[np.clip(np.digitize(dw, DH.BINS) - 1, 0, len(pb) - 1)], 0.0)
                    pw = w3 + (wh & ~w3) * p
                    rw = np.maximum(r3 + (rh & ~r3) * p, pw)
                    o.setdefault("floor", []).append((100 * pw.mean(), 100 * (rw - pw).mean(), 100 * (1 - rw.mean())))
                for k, (wv, rv, inf) in runs.items():
                    o.setdefault(k, []).append(S3.split(wv, rv))
                    o.setdefault(k + "_info", []).append({q: inf[q] for q in ("fired", "none", "steps")})
            if seed == 0:
                if p3(o["main"]) != refA[name]:
                    bad.append((name, "A printed", p3(o["main"]), refA[name]))
                if name != "Enceladus":
                    if p3(o["fit0"]) != dh[name][1]:
                        bad.append((name, "B", p3(o["fit0"]), dh[name][1]))
                    if p3(o["floor"]) != dh[name][2]:
                        bad.append((name, "C", p3(o["floor"]), dh[name][2]))
        print(f"[seed {seed} done, {time.time() - T0:.0f} s]", flush=True)
    ok = lambda tag: not [b for b in bad if b[1].startswith(tag)]
    print(f"\nregression (A) handoff off, h = 1 vs multiwall.multi arrays and seedfix main (printed 0.1, "
          f"{len(refA)} worlds): {'EXACT' if ok('A') and len(refA) == 15 else 'FAILED'}")
    print(f"regression (B) handoff off, fitted h vs deephand_scale.log fitted-h column, worlds 1-14: "
          f"{'EXACT' if ok('B') else 'FAILED'}")
    print(f"regression (C) recomputed floor vs deephand_scale.log ocean-scale column, worlds 1-14: "
          f"{'EXACT' if ok('C') else 'FAILED'}")
    if bad:
        print("  mismatches: " + str(bad))
    return out, mw


def report(out, mw):
    names = sorted({n for _, n in out}, key=lambda n: -out[(0, n)]["main"][0][0])
    fmt = lambda t: f"{t[0]:5.1f} / {t[1]:4.1f} / {t[2]:5.1f}"
    cols = (("main", "multiwall main (seed fix)"), ("x10", "multiwall x10"), ("floor", "ocean-scale floor"),
            ("h1", "handoff, h = 1"), ("fit", "handoff, fitted h"), ("thr1", "THR = 1, fitted h"))
    for li, lab in ((0, "no lid"), (1, "final")):
        print(f"\nP_pved / P_chance / P_none %, {lab}, N = {U.N}, seed 0 (handoff THR = 3 unless stated)")
        print(f"{'world':18} " + " || ".join(f"{c:>22}" for _, c in cols))
        for n in names:
            print(f"{n:18} " + " || ".join(fmt(out[(0, n)][k][li]) for k, _ in cols))
    old = mw["Enceladus"][0]
    print(f"(multiwall.log's collided Enceladus main, for reference: {' / '.join(old[:3])} | {' / '.join(old[3:])})")
    print("\nHandoff (THR = 3, seed 0, no lid): % of recursion draws where it fires | where no layer passes (no hand) "
          "| median layers moved up;  THR = 1: fires | no hand")
    for n in names:
        a, b = out[(0, n)]["fit_info"][0], out[(0, n)]["thr1_info"][0]
        print(f"{n:18} {a['fired']:5.1f} | {a['none']:5.1f} | {a['steps']:4.0f}   ;   {b['fired']:5.1f} | {b['none']:5.1f}")
    print("\nSeeds 0-3, P_pved % no lid / final: handoff fitted h | handoff h = 1; span = max - min over seeds")
    for n in names:
        row, spans = [], []
        for k in ("fit", "h1"):
            v = np.array([[out[(s, n)][k][0][0], out[(s, n)][k][1][0]] for s in SEEDS])
            row.append(" ".join(f"{a:5.1f}/{b:5.1f}" for a, b in v))
            spans.append(np.ptp(v, axis=0).max())
        print(f"{n:18} {row[0]} | {row[1]}   span {spans[0]:.2f} | {spans[1]:.2f}")


def smoke():
    """Self-checks only (closed-form j_max, handoff off == multiwall.multi); nothing result-like is printed."""
    U.N = 2000
    _, W = S.chains(0)
    x, _, a, _ = W["Europa"]
    rng = np.random.default_rng([0, 1])
    ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
    lid = np.zeros(U.N, bool)
    base = S3.outcome3(x, a[1], a[4], ex, lid)[:2]
    wm, rm, _ = M.multi(x, a[1], a[4], ex, lid)
    w0, r0, _ = hand(x, a[1], a[4], ex, lid, base)
    assert (wm == w0).all() and (rm == r0).all()
    for thr in (1.0, 3.0, 30.0):
        _, _, inf = hand(x, a[1], a[4], ex, lid, base, thr=thr)
        j, K, rec = inf["j"], inf["K"], inf["rec"]
        m = rec & (j < K) & (j >= 0)
        assert m.any() and (inf["dw"][m] >= thr * (1 - 1e-9)).all()           # the chosen layer passes
        with np.errstate(all="ignore"):
            nxt = inf["dwf"](j + 1, inf["d_w"] / inf["H"])                    # the full layer just below it fails
        deeper = rec & (j < K - 1)
        assert (nxt[deeper] < thr * (1 + 1e-9)).all()


if __name__ == "__main__":
    T0 = time.time()
    if SMOKE:
        smoke()
        print("smoke ok")
        sys.exit()
    rs = DH.rows()
    cf, _, _ = DH.fit(rs)
    fit = (cf, max(q["K"] for q in rs), min(q["dw"] for q in rs), max(q["dw"] for q in rs))
    print(f"handoff.py, N = {U.N}, seeds {SEEDS}; seed fix applied (default_rng([s, 0]) -> [s, 0, 1]); fitted h "
          f"coefficients {np.round(cf, 3).tolist()}, K clip {fit[1]}, d/w clip [{fit[2]:.2f}, {fit[3]:.1f}]")
    section_validate(rs)
    out, mw = worlds(fit, DH.pbin(rs))
    report(out, mw)
    print(f"\ntotal {time.time() - T0:.0f} s")
