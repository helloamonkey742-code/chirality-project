"""h13h14_extreme: world shares when the sideways-coarsening slowdown f is far beyond 1000, up to f -> infinity. (H13 x H14)

Copy of h13h14_combined.py (which is NOT edited) with the new pre-registration file, strict tie-aware rank reporting,
counts of worlds losing >= 1 / 5 / 10 points, a monotonicity check, the limit check and a `check` mode that compares
the f = 1 and f = 1000 values with h13h14_combined.log. It runs h13_sphere_shares.py UNCHANGED (runpy, same argv) with
coarsen_slow.run wrapped by a cache, so each f is computed once and the script's own lines are reproduced.

f -> infinity: every place the sideways time (L/A)^2/D_h enters is a necessary condition (sweep2.py:127,129,137;
sweep3.py:122-123; multiwall.py:421-422), so above a finite f every draw fails and every share is exactly 0. The limit is
implemented as very large f (1e9 and 1e12; np.inf would make U.A = 0), and the two are compared cell by cell.

Pre-registered in h13h14_extreme_predictions.txt; the hash is asserted here and printed at the top of the log.

Run:   cd <repo> && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy \
         --with numba --with matplotlib python h13h14_extreme.py \
         f=1,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000,100000000,1000000000,1000000000000 \
         mix=h13_sphere_s16_f50.log > h13h14_extreme.log 2> h13h14_extreme.err
Check: python h13h14_extreme.py check          (no sweep; compares f = 1 and f = 1000 with h13h14_combined.log)
Smoke runs before the results run: selfcheck() below on made-up numbers only (no sweep run).
"""
import hashlib
import re
import runpy
import sys

PRED = "h13h14_extreme_predictions.txt"
PRED_SHA = "97a76e595ea94f836c9915d2de0761ed5ea1e673f35cc7fd17c79cff57fa203f"
COLS = ("main no lid", "main final", "x10 no lid", "x10 final")
SCRIPT = "h13_sphere_shares.py"
EE = "Early Earth ocean"
OLD, NEW = "h13h14_combined.log", "h13h14_extreme.log"


def r9(x):
    return round(x, 9)


def losses(R, base):
    """{world: [loss in points per column]} = base - R[f]."""
    return {n: [b - a for a, b in zip(R[n], base[n])] for n in base}


def ranks(R, k):
    """{world: 1 + number of worlds strictly above it}; equal (rounded) shares share a rank."""
    return {n: 1 + sum(r9(R[m][k]) > r9(R[n][k]) for m in R) for n in R}


def ordered_pairs(R, k):
    """Number of pairs of worlds with strictly different shares in column k."""
    ns = list(R)
    return sum(r9(R[a][k]) != r9(R[b][k]) for i, a in enumerate(ns) for b in ns[i + 1:])


def reversed_pairs(R0, R1, k):
    """(reversed, tied) among pairs strictly ordered at f = 1: reversed = now strictly the other way, tied = now equal."""
    ns, rev, tie = list(R0), [], 0
    for i, a in enumerate(ns):
        for b in ns[i + 1:]:
            if r9(R0[a][k]) == r9(R0[b][k]):
                continue
            hi, lo = (a, b) if r9(R0[a][k]) > r9(R0[b][k]) else (b, a)
            if r9(R1[hi][k]) < r9(R1[lo][k]):
                rev.append((hi, lo))
            elif r9(R1[hi][k]) == r9(R1[lo][k]):
                tie += 1
    return rev, tie


def counts(L, k):
    return tuple(sum(v[k] >= T for v in L.values()) for T in (1.0, 5.0, 10.0))


def top(R, k, n=5):
    return sorted(R, key=lambda w: -R[w][k])[:n]


def report(cache, fs):
    base = cache[1.0]
    big = [f for f in fs if f != 1.0]
    print("\n=== per-world loss against f = 1, points (base - f), by column; N = 200,000 so 1 point = 2000 draws ===")
    for k, c in enumerate(COLS):
        print(f"\n[{c}]  base share, then loss at f = " + ", ".join(f"{f:g}" for f in big))
        for n in base:
            print(f"  {n:18} {base[n][k]:6.2f}  " + "  ".join(f"{losses(cache[f], base)[n][k]:7.3f}" for f in big))
    print("\n=== two largest losses per main column (world loss), per f ===")
    for f in big:
        L = losses(cache[f], base)
        s = []
        for k in (0, 1):
            o = sorted(L, key=lambda n: -L[n][k])[:2]
            s.append(f"{COLS[k]}: " + ", ".join(f"{n} {L[n][k]:.2f}" for n in o))
        print(f"f = {f:g}: " + " | ".join(s))
    print("\n=== worlds losing >= 1 / >= 5 / >= 10 points (of 15), per f and column ===")
    for f in big:
        L = losses(cache[f], base)
        print(f"f = {f:g}: " + " | ".join(f"{c}: {'/'.join(str(x) for x in counts(L, k))}" for k, c in enumerate(COLS)))
    print("\n=== strict rank of every world (1 = highest share; ties share a rank), main no lid, then main final ===")
    for k in (0, 1):
        print(f"\n[{COLS[k]}]  f = " + ", ".join(f"{f:g}" for f in fs))
        rk = {f: ranks(cache[f], k) for f in fs}
        for n in base:
            print(f"  {n:18} " + " ".join(f"{rk[f][n]:>3}" for f in fs))
    print("\n=== strictly ordered pairs (of 105) and pairs reversed / tied against f = 1, per f and column ===")
    for f in big:
        s = []
        for k, c in enumerate(COLS):
            rev, tie = reversed_pairs(base, cache[f], k)
            s.append(f"{c}: ordered {ordered_pairs(cache[f], k)}, reversed {len(rev)}, tied {tie}")
        print(f"f = {f:g}: " + " | ".join(s))
    print("\n=== five highest worlds per f, main no lid | main final ===")
    for f in fs:
        print(f"f = {f:g}: " + " | ".join(", ".join(top(cache[f], k)) for k in (0, 1)))
    print("\n=== Early Earth ocean: share, loss and rank in each column ===")
    for f in fs:
        print(f"f = {f:g}: share " + " / ".join(f"{cache[f][EE][k]:.3f}" for k in range(4)) + "  loss " +
              " / ".join(f"{base[EE][k] - cache[f][EE][k]:.4f}" for k in range(4)) + "  rank " +
              " / ".join(str(ranks(cache[f], k)[EE]) for k in range(4)))
    print("\n=== survivors: unrounded non-zero shares (points), f >= 1e7 ===")
    for f in fs:
        if f < 1e7:
            continue
        for k in (0, 2):
            nz = {n: cache[f][n][k] for n in cache[f] if cache[f][n][k] > 0}
            print(f"f = {f:g} {COLS[k]}: " + (", ".join(f"{n} {v:.5f}" for n, v in sorted(nz.items(), key=lambda t: -t[1])) or "none") +
                  f"  | exactly 0.0: {sum(cache[f][n][k] == 0 for n in cache[f])} of 15")
    bad = []
    for n in base:
        for k in range(4):
            for f0, f1 in zip(fs, fs[1:]):
                if cache[f1][n][k] > cache[f0][n][k] + 1e-9:
                    bad.append(f"{n} {COLS[k]} f {f0:g}->{f1:g}: {cache[f0][n][k]:.4f}->{cache[f1][n][k]:.4f}")
    print("\n=== monotone: share non-increasing along the grid, 15 worlds x 4 columns: " + ("yes" if not bad else "NO, " + "; ".join(bad)))
    big_f = [f for f in fs if f >= 1e9]
    if len(big_f) >= 2:
        same = all(cache[big_f[0]][n][k] == cache[f][n][k] for f in big_f[1:] for n in base for k in range(4))
        zero = all(cache[f][n][k] == 0.0 for f in big_f for n in base for k in range(4))
        print(f"=== limit: f = {', '.join(f'{f:g}' for f in big_f)} identical in all 60 cells: {'yes' if same else 'NO'}; all exactly 0.0: {'yes' if zero else 'NO'}")
    print("\nnot a combined result: the 'mix' line above is the sphere factor alone (every draw its own f_t).")


def parse(path):
    """Strings of the values that must agree between the two logs: {key: value-string}."""
    out, col, fl, fpf = {}, None, None, None
    for l in open(path):
        if l.startswith("sideways time / age"):
            fpf = [float(x) for x in l.split("at f = ")[1].split(",")]
            continue
        if fpf and (m := re.match(r"  (\S.*?)\s+worst ([\d.]+)\s+1/worst\s+([\d.]+)\s+(.*)", l)):
            pc = m[4].split()
            out[("flip", m[1], "worst")] = (m[2], m[3])
            for f in (1.0, 1000.0):
                out[("flip", m[1], f)] = pc[fpf.index(f)]
            continue
        if l.startswith("f = 1 against the logged"):
            out[("f1check",)] = l.strip()
        m = re.match(r"(\S.*?)\s+(f\s+[\d.]+:.*)$", l)
        if m and "|" in l and "%" not in l:
            for t in re.finditer(r"f\s+([\d.]+):\s+([\d.]+)/\s*([\d.]+)\s+\|\s+([\d.]+)/\s*([\d.]+)", m[2]):
                if float(t[1]) in (1.0, 1000.0):
                    out[("pved", m[1], float(t[1]))] = t.groups()[1:]
            continue
        if (m := re.match(r"f = (1000\.0|1\.0): largest loss.*", l)):
            out[("largest", float(m[1]))] = l.strip()
        if (m := re.match(r"\[(.+?)\]\s+base share, then loss at f = (.*)", l)):
            col, fl = m[1], [float(x) for x in m[2].split(",")]
            continue
        if l.startswith("==="):
            col = None
        if col and (m := re.match(r"  (\S.*?)\s+(-?\d+\.\d+)\s+((?:-?\d+\.\d+\s*)+)$", l)):
            v = m[3].split()
            out[("base", col, m[1])] = m[2]
            out[("loss", col, m[1], 1000.0)] = v[fl.index(1000.0)]
            continue
        if (m := re.match(r"f = (1|1000): share (.*)$", l)):
            out[("ee", m[1])] = m[2].split("  rank")[0]
    return out


def check():
    a, b = parse(OLD), parse(NEW)
    bad = [k for k in a if k not in b or a[k] != b[k]]
    print(f"REPRO f = 1 and f = 1000, {NEW} against {OLD}: {len(a) - len(bad)} of {len(a)} values identical")
    kinds = {}
    for k in a:
        kinds[k[0]] = kinds.get(k[0], 0) + 1
    print("  compared: " + ", ".join(f"{n} {k}" for k, n in sorted(kinds.items())))
    for k in bad[:20]:
        print(f"  DIFFERENT {k}: old {a[k]} new {b.get(k)}")
    return not bad


def selfcheck():
    a = {"X": [50.0, 50.0, 40.0, 40.0], "Y": [49.0, 49.0, 39.0, 39.0], "Z": [10.0, 10.0, 5.0, 5.0]}
    b = {"X": [47.5, 47.5, 38.0, 38.0], "Y": [49.0, 49.0, 39.0, 39.0], "Z": [10.0, 10.0, 5.0, 5.0]}
    z = {n: [0.0] * 4 for n in a}
    assert abs(losses(b, a)["X"][0] - 2.5) < 1e-12 and losses(b, a)["Y"][0] == 0.0
    assert ranks(a, 0) == {"X": 1, "Y": 2, "Z": 3} and ranks(z, 0) == {"X": 1, "Y": 1, "Z": 1}
    assert ordered_pairs(a, 0) == 3 and ordered_pairs(z, 0) == 0
    assert reversed_pairs(a, b, 0) == ([("X", "Y")], 0) and reversed_pairs(a, a, 0) == ([], 0)
    assert reversed_pairs(a, z, 0) == ([], 3)
    assert counts(losses(b, a), 0) == (1, 0, 0) and counts(losses(z, a), 0) == (3, 3, 3)
    assert top(a, 0, 2) == ["X", "Y"]
    assert hashlib.sha256(b"x").hexdigest().startswith("2d711642")


if __name__ == "__main__":
    selfcheck()
    if sys.argv[1:] == ["check"]:
        sys.exit(0 if check() else 1)
    h = hashlib.sha256(open(PRED, "rb").read()).hexdigest()
    assert h == PRED_SHA, f"{PRED} does not match the recorded hash: {h}"
    print(f"# {PRED} sha256 {h} (asserted)", flush=True)
    import numpy as np
    import h13_sphere  # noqa: F401  (same import order as h13_sphere_shares.py: seedfix is applied by coarsen_slow -> lid_melt)
    import coarsen_slow as C
    cache, _run = {}, C.run

    def run(f, seed=0):
        out = _run(f, seed)
        if np.isscalar(f):
            cache[float(f)] = out
        return out

    C.run = run
    kw = dict(a.split("=") for a in sys.argv[1:])
    fs = sorted(float(v) for v in kw["f"].split(","))
    assert 1.0 in fs and 1000.0 in fs, "f = 1 is the base and f = 1000 the reproduction point"
    sys.argv = [SCRIPT] + sys.argv[1:]
    runpy.run_path(SCRIPT, run_name="__main__")
    report(cache, fs)
