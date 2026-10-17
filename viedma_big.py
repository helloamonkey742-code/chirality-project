"""Bigger Viedma simulation: does the sqrt(N) rule still hold at 256,000 and 1,024,000 crystals?  Imports viedma_sim.py
(its run, fit_a, cp) and edits nothing.

WHY.  viedma_sim.py found that P(initial majority wins) collapses on Delta = eps*sqrt(N) as P = Phi(a*Delta), a = 0.565
(95% CI 0.505-0.629), for N = 250 to 64,000 crystals (a factor 256).  The paper's extrapolation rests on that sqrt(N)
scaling: the excess needed to win shrinks as 1/sqrt(N), not as a fixed fraction.  This script pushes the same mechanism,
unchanged (base variant, S = 2, LAMBDA = 4), 4x and 16x past the largest N run so far, and to a Delta never simulated
(Delta = 3).  At N = 256,000, Delta = 2 is an excess of about 1,000 crystals (0.4%); at N = 1,024,000, Delta = 3 is about
3,000 crystals (0.3%).  If the bias needed were a fixed fraction instead, these excesses would be far too small to steer.

GRID.  N = 256,000: Delta = 1, 2, 3, 30 runs each (90).  N = 1,024,000: Delta = 3, 10 runs.  Seeds are viedma_sim.run's
SeedSequence([20261013, N, round(10*Delta), rep, 0]); these N were never run, so every seed is fresh.  Runs go in 7 chunks
(about 8-9 min each, 2 workers x 1 thread, niced); each finished run is appended to viedma_big_rows.jsonl, and a chunk
skips runs already there.

Run:   for K in 0 1 2 3 4 5 6; do nice -n 10 /opt/anaconda3/bin/python viedma_big.py run $K; done
Check: /opt/anaconda3/bin/python viedma_big.py > viedma_big.log   (reads the rows file and viedma_sim.log; seconds)

LIMITS.  Same as viedma_sim.py: every number comes from the mechanism as written (its ASSUMPTION lines), not from NaClO3
measurements.  This checks that the sqrt(N) scaling holds inside the mechanism over 250 -> 1,024,000 crystals (a factor
4,096).  It is still about 10^24 short of an ocean; no simulation gets there.  Delta = 0 is not rerun (it carries no
information on a).

PREDICTIONS (written 2026-10-16, before any run of the grid).  Seen before writing: one timing run at N = 256,000,
Delta = 0, rep 9999 (a seed outside the grid), printing only its wall time, 53 s.  Reference values are recomputed here
from viedma_sim.log's per-run table (base variant, Delta > 0).
 (V0) Regression: the base a recomputed from viedma_sim.log is 0.565 +- 0.005.
 (V1) Collapse holds at 256,000: a fitted on the N = 256,000 cells alone lies in 0.40-0.75 (base a +- the expected 95%
      sampling width, about 0.18 for 90 runs), and the collapse likelihood-ratio test over N = 250, 1,000, 4,000, 16,000,
      64,000, 256,000 (separate a per N vs one a) gives p > 0.05.
 (V2) New Delta: pooled over N = 256,000 and 1,024,000 (40 runs), the majority wins at Delta = 3 at a rate consistent with
      Phi(3 * base a) = 0.955: exact two-sided binomial p > 0.05.
 (V3) Fixed-fraction rule rejected at the largest size: at N = 1,024,000, Delta = 3 (0.3% excess), the majority wins at
      least 8 of 10 runs (a coin would do that 5.5% of the time; Phi(a*Delta) gives 99%).
 (V4) Decision time keeps growing slowly (viedma_sim: +1.5 time units per factor 4 in N, 14.0 at 64,000): median t_decide
      at 256,000 in 14.0-18.0 and at 1,024,000 in 15.0-21.0.
 (V5) No run undecided by t = 60.
 (V6) Each chunk finishes in under 20 min of wall time (recorded per chunk in the rows file).
"""
import json
import os
import sys
import time
from math import comb
from multiprocessing import Pool

import viedma_sim as V   # first: it sets the BLAS/OMP thread variables to 1 before numpy loads

import numpy as np
from scipy.stats import chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = os.path.join(HERE, "viedma_big_rows.jsonl")
JOBS = ([(1024000, 3.0, rep, 0) for rep in range(10)] +
        [(256000, d, rep, 0) for rep in range(30) for d in (1.0, 2.0, 3.0)])
RES = {}


def chunks():
    """Greedy split of JOBS into chunks of <= 20 cost units (one 256k run = 1, one 1M run = 4.5)."""
    out, cur, c = [], [], 0.0
    for j in JOBS:
        w = 4.5 if j[0] > 500000 else 1.0
        if cur and c + w > 20:
            out.append(cur)
            cur, c = [], 0.0
        cur.append(j)
        c += w
    return out + [cur]


def load():
    if not os.path.exists(ROWS):
        return [], []
    rows, meta = [], []
    for line in open(ROWS):
        d = json.loads(line)
        (meta if "chunk" in d else rows).append(d)
    return rows, meta


def run_chunk(k):
    done = {(r["N"], r["delta"], r["rep"]) for r in load()[0]}
    todo = [j for j in chunks()[k] if (j[0], j[1], j[2]) not in done]
    t0 = time.time()
    print(f"chunk {k}: {len(todo)} runs to do", flush=True)
    with Pool(2) as pool, open(ROWS, "a") as f:
        for r in pool.imap_unordered(V._job, todo):
            f.write(json.dumps({key: (float(v) if isinstance(v, (float, np.floating)) else int(v)) for key, v in r.items()}) + "\n")
            f.flush()
        f.write(json.dumps({"chunk": k, "runs": len(todo), "wall_s": round(time.time() - t0, 1)}) + "\n")
    print(f"chunk {k}: {time.time() - t0:.0f} s", flush=True)


def base_rows():
    rows, f = [], False
    for line in open(os.path.join(HERE, "viedma_sim.log")):
        if line.startswith("PER-RUN TABLE"):
            if f:
                break
            f = True
            continue
        p = line.split()
        if f and len(p) == 11 and p[0] == "0":
            rows.append(dict(N=int(p[1]), delta=float(p[2]), eps=float(p[4]), win=int(p[6]), t=float(p[9])))
    return rows


def cells(rows):
    """(N, Delta_nominal, wins, decided, Delta_actual) for Delta > 0."""
    out = {}
    for r in rows:
        if r["delta"] > 0 and r["win"] >= 0:
            c = out.setdefault((r["N"], r["delta"]), [0, 0, r["eps"] * np.sqrt(r["N"])])
            c[0] += r["win"]
            c[1] += 1
    return [(N, d, k, n, da) for (N, d), (k, n, da) in sorted(out.items())]


def two_sided(k, n, p):
    pmf = [comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(n + 1)]
    return min(1.0, sum(x for x in pmf if x <= pmf[k] * (1 + 1e-9)))


def check(p, ok, detail):
    RES[p] = bool(ok)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def main():
    base, (rows, meta) = base_rows(), load()
    print(f"viedma_big.py: {len(rows)} new runs read from viedma_big_rows.jsonl, {len(base)} base runs from viedma_sim.log")
    cb, cn = cells(base), cells(rows)
    a0, (lo0, hi0), _ = V.fit_a([(k, n, da) for _, _, k, n, da in cb])
    print(f"\nN        Delta  runs  wins  P(win)  95% CI          Phi(a0*Delta)")
    for N, d, k, n, da in cn:
        lo, hi = V.cp(k, n)
        print(f"{N:<8d} {d:<5.1f}  {n:<4d}  {k:<4d}  {k / n:.3f}   [{lo:.3f}, {hi:.3f}]  {norm.cdf(a0 * da):.3f}")
    print(f"\nbase a0 = {a0:.3f} [{lo0:.3f}, {hi0:.3f}] (viedma_sim.log, N = 250-64,000)")
    check("V0", abs(a0 - 0.565) <= 0.005, f"a0 {a0:.3f}")

    c256 = [(k, n, da) for N, d, k, n, da in cn if N == 256000]
    a256, (l256, h256), ll256 = V.fit_a(c256)
    perN = {}
    for N in sorted({c[0] for c in cb} | {256000}):
        perN[N] = V.fit_a([(k, n, da) for M, _, k, n, da in cb + cn if M == N])
    six = [(k, n, da) for N, _, k, n, da in cb + cn if N <= 256000]
    a6, (l6, h6), ll6 = V.fit_a(six)
    lr = 2 * (sum(v[2] for v in perN.values()) - ll6)
    p6 = chi2.sf(max(lr, 0), len(perN) - 1)
    print("per-N a: " + ", ".join(f"{N}: {v[0]:.2f} [{v[1][0]:.2f}, {v[1][1]:.2f}]" for N, v in perN.items()))
    print(f"one a over N = 250-256,000: {a6:.3f} [{l6:.3f}, {h6:.3f}]; collapse LRT chi2 = {lr:.2f}, df = {len(perN) - 1}, p = {p6:.3f}")
    check("V1", 0.40 <= a256 <= 0.75 and p6 > 0.05, f"a(256,000) = {a256:.3f} [{l256:.3f}, {h256:.3f}]; collapse p = {p6:.3f}")

    d3 = [(k, n) for N, d, k, n, _ in cn if d == 3.0]
    k3, n3 = sum(k for k, _ in d3), sum(n for _, n in d3)
    p3 = norm.cdf(3 * a0)
    check("V2", n3 and two_sided(k3, n3, p3) > 0.05, f"Delta = 3: {k3}/{n3} wins vs Phi(3 a0) = {p3:.3f}, p = {two_sided(k3, n3, p3) if n3 else float('nan'):.3f}")

    k1m = [(k, n) for N, d, k, n, _ in cn if N == 1024000]
    check("V3", k1m and k1m[0][0] >= 8, f"N = 1,024,000, Delta = 3: {k1m[0][0] if k1m else '-'} of {k1m[0][1] if k1m else '-'}")

    med = {N: np.median([r["t"] for r in rows if r["N"] == N and r["win"] >= 0]) for N in (256000, 1024000)}
    check("V4", 14.0 <= med[256000] <= 18.0 and 15.0 <= med[1024000] <= 21.0,
          f"median t_decide {med[256000]:.2f} (256,000), {med[1024000]:.2f} (1,024,000); 64,000 was 14.0")
    und = sum(r["win"] < 0 for r in rows)
    check("V5", und == 0 and len(rows) == len(JOBS), f"{und} undecided of {len(rows)} (grid {len(JOBS)})")
    walls = {m["chunk"]: m["wall_s"] for m in meta}
    check("V6", len(walls) == len(chunks()) and max(walls.values()) < 1200, "chunk wall s " + str(walls))
    print(f"\n{sum(RES.values())} of {len(RES)} checked items HELD: " + ", ".join(f"{k} {'HELD' if v else 'FAILED'}" for k, v in sorted(RES.items())))


if __name__ == "__main__":
    if sys.argv[1:2] == ["run"]:
        run_chunk(int(sys.argv[2]))
    else:
        main()
