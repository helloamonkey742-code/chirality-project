"""h13h14_extreme_forecast: the forecast behind the bands in h13h14_extreme_predictions.txt. Runs NO sweep.

Model (stated, not tested beyond the check below): a draw can only succeed if f * ratio <= 1, ratio = (L/5.5)^2 / (D_h * t),
L = pi R, D_h and t log-uniform over the ranges declared in the sweep (W below, read from sweep2/sweep3/uncertainty by hand).
So p_flip(f) = P(f * ratio > 1) is analytic, and share(f) = base * (1 - p_flip) * lam, with lam the share among draws that do not
flip, relative to the base share (h13h14_combined.log: 0.99-1.00 for main columns at f = 100-1000, Mimas 0.95, Enceladus 1.03).
Check against the logged flip shares: python h13h14_extreme_forecast.py check   (needs h13h14_combined.log)
Run: UV_CACHE_DIR=$HOME/Library/Caches/uv uv run --no-project --with numpy python h13h14_extreme_forecast.py > h13h14_extreme_forecast.out
"""
import re
import sys

import numpy as np

YR = 3.156e7  # ocean.py:26
D, T = (1e-2, 1e2), (1e8, 4e9)
W = {  # name: (R in m, D_h range, t range in yr)
    "Enceladus": (232e3, D, (1e6, 1e9)), "Europa": (1460e3, D, T), "Early Earth ocean": (6.371e6, (1e2, 1e4), (1e7, 5e8)),
    "Ganymede": (2631e3, D, T), "Titan": (2575e3, D, T), "Callisto": (2410e3, D, T), "Triton": (1353e3, D, T),
    "Pluto": (1188e3, D, T), "Titania": (789e3, D, T), "Oberon": (761e3, D, T), "Umbriel": (585e3, D, T),
    "Ariel": (579e3, D, T), "Dione": (561e3, D, T), "Miranda": (236e3, D, (1e8, 5e8)), "Mimas": (198e3, D, (2e6, 2.5e7))}
EE = "Early Earth ocean"
COLS = ["main no lid", "main final", "x10 no lid", "x10 final"]


def pflip(name, f, A=5.5):
    R, (d1, d2), (t1, t2) = W[name]
    a1, a2, b1, b2 = np.log10(d1), np.log10(d2), np.log10(t1 * YR), np.log10(t2 * YR)
    u = np.linspace(a1, a2, 200001)
    return float(np.mean(np.clip((np.log10((np.pi * R / A) ** 2 * f) - u - b1) / (b2 - b1), 0, 1)))


def first_all_flip(name, A=5.5):
    R, (d1, d2), (t1, t2) = W[name]
    return d2 * t2 * YR / (np.pi * R / A) ** 2   # f above this flips every draw


def bases(path="h13h14_combined.log"):
    base, k = {}, None
    for l in open(path):
        m = re.match(r"\[(.+?)\]\s+base share", l)
        if m:
            k = COLS.index(m.group(1))
        elif l.startswith("=== worlds at"):
            k = None
        elif k is not None and (m := re.match(r"  (\S.*?)\s+([\d.]+)\s+[\d.-]+", l)):
            base.setdefault(m.group(1).strip(), [None] * 4)[k] = float(m.group(2))
    return base


def lam_range(n, k):
    if n == EE:
        return (0.995, 1.0007) if k < 2 else (0.99, 1.0033)   # 100/base is the ceiling
    return (0.90, 1.12) if k < 2 else (0.85, 1.12)


def check():
    log = open("h13h14_combined.log").read().splitlines()
    fs = [1, 17, 34, 100, 136, 300, 1000]
    print("logged flip share / analytic, f = 100, 136, 300, 1000, in percent")
    for l in log[2:17]:
        m = re.match(r"\s+(.+?)\s+worst ([\d.]+)\s+1/worst\s+([\d.]+)\s+(.*)", l)
        n, pc = m.group(1), [float(x.strip("%")) for x in m.group(4).split()]
        print(f"  {n:18} " + " ".join(f"{pc[i]:6.3f}/{100 * pflip(n, fs[i]):6.3f}" for i in range(3, 7)))


def main():
    base, rng = bases(), np.random.default_rng(1)
    names, fs = list(W), [3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 1e8]
    print("f above which every draw of the world flips (analytic): " +
          ", ".join(f"{n} {first_all_flip(n):.2e}" for n in sorted(names, key=first_all_flip)[-5:]))
    print("p_flip %, per world and f")
    print(f"{'':18}" + "".join(f"{f:9.0e}" for f in fs))
    for n in names:
        print(f"{n:18}" + "".join(f"{100 * pflip(n, f):9.3f}" for f in fs))
    for f in fs:
        p = np.array([pflip(n, f) for n in names])
        print(f"\nf = {f:.0e}")
        for k in (0, 1):
            b = np.array([base[n][k] for n in names])
            lo, hi = np.array([lam_range(n, k) for n in names]).T
            S = b * (1 - p) * rng.uniform(lo, hi, (4000, len(names)))
            L = b - np.minimum(S, b)
            o = np.argsort(-b)
            rev = np.array([sum(s[o[a]] < s[o[c]] for a in range(15) for c in range(a + 1, 15)) for s in S[:500]])
            ee = names.index(EE)
            rank = (S > S[:, [ee]]).sum(1) + 1
            q = lambda a: f"{int(np.percentile(a, 1))}-{int(np.percentile(a, 99))}"
            ml = L.max(1)
            w, c = np.unique(np.array(names)[L.argmax(1)], return_counts=True)
            print(f"  {COLS[k]:12} counts >=1 {q((L >= 1).sum(1))}, >=5 {q((L >= 5).sum(1))}, >=10 {q((L >= 10).sum(1))}; "
                  f"reversed pairs {q(rev)}; EE rank {q(rank)}; largest loss {np.percentile(ml, 1):.1f}-{np.percentile(ml, 99):.1f} "
                  + ", ".join(f"{x} {100 * y / len(L):.0f}%" for x, y in zip(w, c)))
        e = base[EE]
        print("  EE share (central, lam = 1): main " f"{e[0] * (1 - p[names.index(EE)]):.3f}, x10 {e[2] * (1 - p[names.index(EE)]):.3f}")


if __name__ == "__main__":
    check() if sys.argv[1:] == ["check"] else main()
