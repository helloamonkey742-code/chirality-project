"""h24_handoff_misses: what the misses of V1 / V3 look like on the fresh, thin and old runs (descriptive). (2026-10-21)

No prediction, no pass / fail: it re-reads h24_handoff_runs_*.log (written by `h24_handoff.py run`) after
h24_handoff_test.log was produced and tabulates the misses by type and by the rule domain's d/w. Written after the
test results were seen (so it is exploratory; nothing in the notes rests on it as a test).
  over-fire  : the rule handed up (chosen domain above the rule domain) but the rule domain's hand was right
  under-fire : the rule kept the rule domain, whose hand lost
  even skip  : the rule handed up and the final hand matches neither the chosen nor the rule domain's hand
               (hands alternate, so the winner sits an even number of domains above the rule domain)
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba --with matplotlib \
       python h24_handoff_misses.py > h24_handoff_misses.log
"""
import numpy as np

import h24_handoff as H

BANDS = ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 3), (3, 1e9))


def kinds(rows, cache, var, par):
    """Per row: 0 hit, 1 over-fire, 2 under-fire, 3 even skip."""
    out = []
    for q in rows:
        ds, d0 = cache[id(q)]
        h, d, _ = H.predict(var, par, (ds, d0))
        r = ds[d0]["hand"]
        out.append(0 if h == q["fin"] else 2 if d == d0 else 1 if q["fin"] == r else 3)
    k = np.array(out)
    assert ((k == 0) == H.agree(rows, cache, var, par)[:, 0]).all()
    return k


def band(rows, cache):
    return np.array([next(i for i, (a, b) in enumerate(BANDS) if a <= cache[id(q)][0][cache[id(q)][1]]["dwm"] < b) for q in rows])


def main():
    rs, cache, items = H.load_old()
    fits, _ = H.fit_variants(items)
    H.frozen(fits)
    R = H.load_runs()
    for q in R["F"] + R["T"]:
        cache[id(q)] = H.doms(q)
    print("miss types (over-fire / under-fire / even skip) and agreement by the rule domain's d/w band; V1 = THR 1.25, V3 frozen")
    for lab, rows in (("old (V3 in sample)", rs), ("fresh", R["F"]), ("thin stacks", R["T"])):
        b = band(rows, cache)
        for v, par in (("V1", None), ("V3", fits["V3"])):
            k = kinds(rows, cache, v, par)
            print(f"\n{lab}, {v}: {len(rows)} runs, misses {int((k > 0).sum())} = over {int((k == 1).sum())} + under {int((k == 2).sum())}"
                  f" + even skip {int((k == 3).sum())}")
            print("  rule d/w band      runs  misses  agreement %   (over/under/skip)")
            for i, (lo, hi) in enumerate(BANDS):
                m = b == i
                if m.any():
                    print(f"  [{lo:g}, {hi:g}){'':8}{int(m.sum()):6d}{int((k[m] > 0).sum()):8d}{100 * (k[m] == 0).mean():12.1f}   "
                          f"({int((k[m] == 1).sum())}/{int((k[m] == 2).sum())}/{int((k[m] == 3).sum())})")
    F = R["F"]
    k = kinds(F, cache, "V3", fits["V3"])
    print("\nfresh V3 misses by (w, n): misses / runs")
    for w in H.FRESH_W:
        print(f"  w {w:g}: " + "; ".join(f"n={n}: {int(((k > 0) & (np.array([q['w'] for q in F]) == w) & (np.array([q['n'] for q in F]) == n)).sum())}"
                                         f" / {sum(1 for q in F if q['w'] == w and q['n'] == n)}" for n in (3, 6, 10)))
    print("\nthin stacks V3 by w: " + "; ".join(
        f"w {w:g}: agreement {100 * np.mean([x == 0 for x, q in zip(kinds(R['T'], cache, 'V3', fits['V3']), R['T']) if q['w'] == w]):.1f} % "
        f"of {sum(1 for q in R['T'] if q['w'] == w)}" for w in H.THIN_W))


if __name__ == "__main__":
    main()
