"""Blind units test 2: Jang et al. 2024 (Lab Chip 24, 5055, doi:10.1039/D4LC00658E), NaClO3 antisolvent droplets.

Predictions are in blind_jang_predictions.txt, written before any outcome was read:
  sha256 bd2f04e3f1de0c962099b34e0570e607d19cb3084d5484cebea68ffc4f4ee557, 2026-10-20 17:49 PDT (2026-10-21 00:49 UTC)
  (main() asserts the file still has that hash)
Outcomes are entered below by a separate reader, in a later commit, without changing the predictions file.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python blind_jang.py > blind_jang.log
"""
import hashlib
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, mannwhitneyu

from blind_units import exp_fit, free_width

PRED_SHA = "bd2f04e3f1de0c962099b34e0570e607d19cb3084d5484cebea68ffc4f4ee557"
TOTAL = {35: 379, 40: 618, 45: 276}   # crystal-bearing droplets analysed, per Tsat (text)

# Outcomes: entered after the predictions commit by a separate reader.  None = not yet entered.
# {Tsat: {"L": [d_e um, ...], "D": [...], "LD": [...]}} from the markers of Figs 3c, 5a, 5b (see the data rule).
DROPS = None
RES = {}


def check(p, ok, detail):
    ok = ok if isinstance(ok, str) else bool(ok)
    RES[p] = ok
    print(f"  {p} {'HELD' if ok is True else 'FAILED' if ok is False else ok}: {detail}")


def vol(d):
    return np.pi / 6 * np.asarray(d, float) ** 3 * 1e-6   # um -> nL


def series(s):
    """Per-droplet arrays: V (nL), single (1) / mixed (0), d_e (um)."""
    d = np.array(s["L"] + s["D"] + s["LD"], float)
    y = np.r_[np.ones(len(s["L"]) + len(s["D"])), np.zeros(len(s["LD"]))]
    return vol(d), y, d


def v0(V, y):
    return np.inf if y.all() else exp_fit(V, y, np.ones_like(y))


def self_check():
    rng = np.random.default_rng(5)
    V = vol(rng.uniform(100, 700, 600))
    y = (rng.random(600) < np.exp(-V / 40)).astype(float)
    h_cnt = free_width(V, y, np.ones_like(y))[1]
    h_step = free_width(V, (V < 40).astype(float), np.ones_like(V))[1]
    assert 0.9 < h_cnt < 2.2 and h_step > 10, (h_cnt, h_step)
    assert 25 < v0(V, y) < 60 and v0(V, np.ones_like(V)) == np.inf
    assert abs(vol(1000) - 523.6) < 0.1


def main():
    pred = Path(__file__).with_name("blind_jang_predictions.txt").read_bytes()
    assert hashlib.sha256(pred).hexdigest() == PRED_SHA, "predictions file changed"
    self_check()
    if DROPS is None:
        print("Predictions hashed; outcomes not yet entered (see the docstring).")
        return

    v0s = {}
    for t, s in sorted(DROPS.items()):
        V, y, d = series(s)
        n_mix, n_one = int((y == 0).sum()), int(y.sum())
        print(f"\n=== Tsat {t} C: L {len(s['L'])}, D {len(s['D'])}, mixed {n_mix}; entered {len(y)} of {TOTAL[t]} ===")
        bins = {}
        for b, yi in zip((d // 50).astype(int), y):
            bins.setdefault(b, []).append(yi)
        for b in sorted(bins):
            print(f"  d_e {50 * b}-{50 * b + 50} um: {len(bins[b])} drops, P_single {np.mean(bins[b]):.2f}")
        tag = {35: "J1a", 40: "J1b", 45: "J1c"}[t]
        if n_mix < 5 or n_one < 5:
            check(tag, "NOT SCORABLE", f"{n_mix} mixed, {n_one} single")
        else:
            p = mannwhitneyu(d[y == 0], d[y == 1], alternative="greater").pvalue
            check(tag, p < 0.05, f"mixed median {np.median(d[y == 0]):.0f} um vs single {np.median(d[y == 1]):.0f} um, "
                                 f"one-sided p = {p:.2g}")
        Vk, h = free_width(V, y, np.ones_like(y))
        Ps = [np.mean(v) for v in bins.values() if len(v) >= 10]
        ok = n_mix >= 10 and n_one >= 10 and Ps and max(Ps) >= 0.6 and min(Ps) <= 0.4
        RES.setdefault("_J2", []).append((t, h, Vk, ok))
        v0s[t] = v0(V, y)
        print(f"  reported: mixed share {n_mix / len(y):.3f}, free-width h = {h:.2f} (Vk = {Vk:.3g} nL), "
              f"V0 = {v0s[t]:.3g} nL, L/D {len(s['L'])}/{len(s['D'])}")

    print("\n=== pooled ===")
    j2 = [r for r in RES.pop("_J2") if r[3]]
    if j2:
        check("J2", all(h <= 3 for _, h, _, _ in j2), "; ".join(f"{t} C h = {h:.2f}" for t, h, _, _ in j2) + ", threshold 3")
    else:
        check("J2", "NOT SCORABLE", "no series has >= 10 mixed and >= 10 single with binned P_single spanning 0.4-0.6")
    L = [x for s in DROPS.values() for x in s["L"]]
    D = [x for s in DROPS.values() for x in s["D"]]
    p3 = binomtest(len(L), len(L) + len(D)).pvalue
    check("J3", p3 > 0.05, f"L {len(L)} / D {len(D)}, two-sided p(50/50) = {p3:.3f}")
    p4 = mannwhitneyu(L, D, alternative="two-sided").pvalue
    check("J4", p4 > 0.05, f"L median {np.median(L):.0f} um, D median {np.median(D):.0f} um, two-sided p = {p4:.3f}")
    if sum(v < np.inf for v in v0s.values()) < 2:
        check("J5", "NOT SCORABLE", "fewer than two series with a mixed droplet")
    else:
        check("J5", v0s[35] > v0s[40] > v0s[45], ", ".join(f"V0({t} C) = {v:.3g} nL" for t, v in sorted(v0s.items())))
    held = sum(v is True for v in RES.values())
    print(f"\n{held} of {len(RES)} predictions HELD: " + ", ".join(
        f"{k} {v if isinstance(v, str) else 'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    main()
