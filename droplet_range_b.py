"""Volume range of the droplet test (journals/droplet_prereg.md plan: 8 levels about 0.2-1.6 x V0_pilot).  Simulated drops, not data.

Where must the true critical volume V* sit for the frozen D2 rule (rule B: measured misread, c* -2.21) to call a sharp step STEP, and
what do extra levels cost?  u = V*/V0_pilot.  Stages: 1 u-scan with V0_pilot imposed (20 u x s 0.1, 0.15 x misread 0, 0.03, 0.06 x 4 designs
x 300 runs); 2 the known issue (rules B, E0 = both fits without the misread term, A; 2000 runs); 3 full pipeline with the pilot at
absolute V* (13 values x 4 designs x 600 runs, the pilot fitted by blind_units.exp_fit); 4 false STEP of each design on the 87 counting
and drying cells at the frozen c* (2 conditions x 300 runs); 5 P17-P24.  E0 uses the same c* as B (-2.21), not recalibrated.
Designs (droplet_range_common.DESIGNS): base 8 levels 0.2-1.6, extension cap 2; ext4 same with cap 4; K10 10 levels 0.2-2.9; K12 12 levels 0.11-2.9.
Predictions: droplet_range_predictions_v2.txt, sha256 PRED_SHA, `shasum -a 256` and `date` 2026-10-21 01:11:44 PDT, written before this
script's results run.  v2 replaces droplet_range_predictions.txt (sha256 4ddcacbe..., 00:39:04 PDT, never edited) only to lower the run counts
(load average 40-69); predictions and decision rules are unchanged.  Disclosed smoke runs before it (none used for a claim): time1.py, time2.py (see droplet_range_a.py), mkcache.py.
The first launch of this script aborted in self_check (interp tolerance 1e-9; C.LF is not symmetric about 0) before any simulation; the tolerance
was loosened to 1e-5 and the script rerun.  droplet_range_b_pool.py and droplet_range_b_ns.py are POST HOC reads of this script's cache.

Run: DR_CACHE=<scratch dir> UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python droplet_range_b.py > droplet_range_b.log
"""
import hashlib
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import ndtr

import droplet_range_common as C
from droplet_range_common import NAMES

PRED_SHA = "e81b8eb4e8342e99bc0d8c4253c6b1f233290d6005fe1e2f90f19ec10cbdf9dd"
PROCS, REPS, PREPS = 4, 300, 600
U = (0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.35, 0.4, 0.57, 0.8, 1.0, 1.13, 1.39, 1.6, 1.9, 2.26, 2.7, 3.2, 4.5, 6.4)
S1, S2 = (0.35, 1.39), (0.14, 2.8)
V0I = 3.1                                                  # V0_pilot imposed in the u-scan; only u matters
EPS, SS, DS = (0.0, 0.03, 0.06), (0.1, 0.15), ("base", "ext4", "K10", "K12")
VSTAR = (0.1, 0.2, 0.4, 0.7, 1, 2, 3.1, 4, 7, 12, 20, 32, 64)
CSTAR = {"A": C.C_A, "B": C.C_B, "E0": C.C_B}
NOMINAL = {d: 25 * len(C.DESIGNS[d]["main"]) + C.NC for d in DS}
T0 = time.time()


def step(vc, s):
    return ndtr(-(C.LF - np.log(vc)) / s)


def tasks_of(cells, reps, **kw):
    return [dict(cells=b, reps=reps, **kw) for b in C.batches(cells, max(1, 8000 // reps))]


def run_all(jobs, label):
    """jobs: name -> list of tasks.  One pool for all of them; returns name -> flat list of per-cell results."""
    flat = [t for ts in jobs.values() for t in ts]
    out = C.pmap(flat, PROCS, label)
    res, k = {}, 0
    for n, ts in jobs.items():
        res[n] = [o for b in out[k:k + len(ts)] for o in b]
        k += len(ts)
    print(f"  [{time.time() - T0:5.0f}s] {label} done", file=sys.stderr, flush=True)
    return res


def cached(name, fn):
    d = os.environ.get("DR_CACHE")                         # resume aid only; seeds are fixed, so results do not depend on it
    f = Path(d) / f"b_{name}.pkl" if d else None
    if f and f.exists():
        return pickle.loads(f.read_bytes())
    r = fn()
    if f:
        f.write_bytes(pickle.dumps(r))
    return r


def sh(o, rule="B", v="STEP"):
    c = o["counts"][rule]
    return c[NAMES.index(v)] / c.sum()


def edges(share):
    """share: u -> STEP share.  u_lo = smallest grid u with share >= 0.95; u_hi = largest u with share >= 0.95 at every grid u from 0.35 up to it."""
    lo = min((u for u in U if share[u] >= 0.95), default=None)
    hi = None
    for u in U:
        if u < 0.35:
            continue
        if share[u] < 0.95:
            break
        hi = u
    return lo, hi


def self_check():
    C.self_check()
    sh_ = {u: 0.99 if 0.2 <= u <= 2.26 else 0.4 for u in U}
    assert edges(sh_) == (0.2, 2.26)
    sh_[1.0] = 0.9                                         # a hole at u = 1: u_lo unchanged, u_hi falls to the grid point before it
    assert edges(sh_) == (0.2, 0.8)
    assert sh(dict(counts={"B": np.array([8, 1, 1, 0, 0])})) == 0.8 and abs(np.interp(0.0, C.LF, step(1.0, 0.1)) - 0.5) < 1e-5


def main():
    assert hashlib.sha256(Path(__file__).with_name("droplet_range_predictions_v2.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    print("Model study on simulated drops, not data.  Rule B (measured misread from 128 controls, c* -2.21), CV 0.05 unless stated.")
    print(f"u = V* / V0_pilot.  Prior spans: S1 = {S1} (plan's own allowance), S2 = {S2} (0.2-4 x).  Designs, nominal drops (levels x 25 + 128 controls): "
          + ", ".join(f"{d} {NOMINAL[d]}" for d in DS))

    # stage 1: u-scan with V0_pilot imposed
    cells1, jobs = {}, {}
    for d in DS:
        for si, s in enumerate(SS):
            for ei, e in enumerate(EPS):
                cells1[(d, s, e)] = [(f"u={u}", step(u * V0I, s), 0 * C.LF, 0.05, e, "xor", 0.0, [20, DS.index(d), si, ei, ui]) for ui, u in enumerate(U)]
                jobs[("scan", d, s, e)] = tasks_of(cells1[(d, s, e)], REPS, rules=("B", "E0") if d == "base" else ("B",), design=d, V0=V0I, cstar=CSTAR)
    # stage 2: the known issue, 2000 runs
    jobs["issue"] = tasks_of([("u=1.39", step(1.39 * V0I, 0.1), 0 * C.LF, 0.05, 0.03, "xor", 0.0, [21, 0])], 2000, rules=("A", "B", "E0"), design="base", V0=V0I, cstar=CSTAR)
    # stage 3: full pipeline at absolute V*
    for di, d in enumerate(DS):
        jobs[("pipe", d)] = tasks_of([(f"V*={v}", step(v, 0.1), 0 * C.LF, 0.05, 0.03, "xor", 0.0, [22, di, vi]) for vi, v in enumerate(VSTAR)], PREPS,
                                     rules=("B",), design=d, V0=None, cstar=CSTAR, exact=True, want_v0=True)
    # stage 4: false STEP at the frozen c*
    null = C.null_cells(PROCS)
    assert len(null) == 87
    NC2 = [(0.05, 0.03), (0.15, 0.06)]
    for di, d in enumerate(DS):
        jobs[("null", d)] = tasks_of([(n, P, p0, cv, e, "xor", 0.0, [23, di, i, j]) for i, (n, P, p0) in enumerate(null) for j, (cv, e) in enumerate(NC2)],
                                     300, rules=("B",), design=d, V0=None, cstar=CSTAR)
    R = cached("all", lambda: run_all(jobs, "droplet_range_b"))

    # tables
    SC = {k: {u: o for u, o in zip(U, R[("scan",) + k])} for k in cells1}
    print(f"\nu-scan, V0_pilot imposed, STEP share ({REPS} runs).  Columns u = " + " ".join(f"{u:g}" for u in U))
    for d in DS:
        for s in SS:
            for e in EPS:
                print(f"  {d:4} s={s} misread {e:.2f}: " + " ".join(f"{sh(SC[(d, s, e)][u]):.2f}" for u in U))
    print("\nEdges (smallest u with STEP >= 0.95; largest u reached with STEP >= 0.95 at every grid u from 0.35 up)")
    ED = {}
    for d in DS:
        for s in SS:
            for e in EPS:
                ED[(d, s, e)] = edges({u: sh(SC[(d, s, e)][u]) for u in U})
                print(f"  {d:4} s={s} misread {e:.2f}: u_lo {ED[(d, s, e)][0]}, u_hi {ED[(d, s, e)][1]}")
    inside = lambda lo, hi: [u for u in U if lo <= u <= hi]
    lowest = lambda d, s, e, sp: min(sh(SC[(d, s, e)][u]) for u in inside(*sp))
    print("\nLowest STEP share over the grid points of each span")
    for d in DS:
        print(f"  {d:4} " + "; ".join(f"{n} s={s} misread {e:.2f}: {lowest(d, s, e, sp):.2f}" for n, sp in (("S1", S1), ("S2", S2)) for s in SS for e in (0.03, 0.06)))

    I = R["issue"][0]
    print("\nKnown issue: base design, u = 1.39, s = 0.1, CV 0.05, true misread 0.03, 2000 runs, STEP share:")
    print("  rule A (fixed 3% term, c* -1.61) " + f"{sh(I, 'A'):.3f}; rule B (measured) {sh(I, 'B'):.3f}; E0 (no misread term in either fit, c* -2.21) {sh(I, 'E0'):.3f}")
    print("  all verdicts, B: " + " ".join(f"{n} {x:.3f}" for n, x in zip(NAMES, I["counts"]["B"] / 2000)) + "; E0: " + " ".join(f"{n} {x:.3f}" for n, x in zip(NAMES, I["counts"]["E0"] / 2000)))
    print("  E0 and B across the u grid, base, s = 0.1, misread 0.03: " + " ".join(f"u={u:g}: {sh(SC[('base', 0.1, 0.03)][u], 'E0'):.2f}/{sh(SC[('base', 0.1, 0.03)][u]):.2f}" for u in U))

    print(f"\nFull pipeline (pilot ladder 0.25-16 uL, 20 drops; exact pilot fit), s = 0.1, CV 0.05, misread 0.03, {PREPS} runs: STEP share, u = V* / V0_pilot 2.5 / 50 / 97.5 percentiles")
    PI = {d: {v: o for v, o in zip(VSTAR, R[("pipe", d)])} for d in DS}
    PU = {v: np.percentile(v / PI["base"][v]["V0"], [2.5, 50, 97.5]) for v in VSTAR}
    print("  V* (uL):    " + " ".join(f"{v:>6g}" for v in VSTAR))
    for d in DS:
        print(f"  STEP {d:4}  " + " ".join(f"{sh(PI[d][v]):6.3f}" for v in VSTAR))
    print("  u 2.5%      " + " ".join(f"{PU[v][0]:6.2f}" for v in VSTAR) + "\n  u median    " + " ".join(f"{PU[v][1]:6.2f}" for v in VSTAR)
          + "\n  u 97.5%     " + " ".join(f"{PU[v][2]:6.2f}" for v in VSTAR))
    print("  other verdicts, base: " + "; ".join(f"V*={v:g}: " + " ".join(f"{n} {x:.2f}" for n, x in zip(NAMES[1:], PI['base'][v]['counts']['B'][1:] / PREPS)) for v in (0.1, 0.2, 32, 64)))

    nm = {d: [(sh(R[("null", d)][i * 2 + j]), null[i][0], NC2[j]) for i in range(87) for j in range(2)] for d in DS}
    nmax = {d: max(nm[d]) for d in DS}
    print("\nFalse STEP at the frozen c* = -2.21 (87 counting+drying cells, CV 0.05 / misread 0.03 and CV 0.15 / misread 0.06, 300 runs): highest STEP share")
    for d in DS:
        print(f"  {d:4} {nmax[d][0]:.3f} ({nmax[d][1]}, CV {nmax[d][2][0]}, misread {nmax[d][2][1]})")

    meand = {d: np.mean([o["drops"] for u in inside(*S2) for o in [SC[(d, 0.1, 0.03)][u]]]) + C.NC for d in DS}
    print("\nCost (drops per run incl. 128 controls; measured mean over S2 grid points, s = 0.1, misread 0.03): " + ", ".join(f"{d} {meand[d]:.0f} (nominal {NOMINAL[d]})" for d in DS))

    # predictions
    cond = lambda sp, d: all(lowest(d, s, e, sp) >= 0.95 for s in SS for e in (0.03, 0.06))
    b3 = {u: sh(SC[("base", 0.1, 0.03)][u]) for u in U}
    lo, hi = edges(b3)
    klo, khi = edges({u: sh(SC[("K10", 0.1, 0.03)][u]) for u in U})
    ch = [
        ("P17", all(sh(SC[("base", s, e)][u]) >= (0.95 if e < 0.05 else 0.90) for s in SS for e in EPS for u in inside(*S1)),
         "lowest base STEP on S1: " + ", ".join(f"s={s} misread {e}: {lowest('base', s, e, S1):.3f}" for s in SS for e in EPS)),
        ("P18", sh(I, "E0") <= sh(I, "B") - 0.10, f"E0 {sh(I, 'E0'):.3f} vs measured {sh(I, 'B'):.3f}"),
        ("P19", lo is not None and hi is not None and 0.14 <= lo <= 0.28 and 1.9 <= hi <= 3.2, f"base u_lo {lo}, u_hi {hi}"),
        ("P20", all(b3[u] <= 0.5 for u in (0.05, 0.07, 4.5, 6.4)), f"base STEP at u = 0.05 {b3[0.05]:.2f}, 0.07 {b3[0.07]:.2f}, 4.5 {b3[4.5]:.2f}, 6.4 {b3[6.4]:.2f}"),
        ("P21", all(PU[v][0] >= 0.35 and PU[v][2] <= 1.39 for v in (0.7, 2, 3.1, 7)), "central 95% of u: " + "; ".join(f"V*={v:g}: {PU[v][0]:.2f}-{PU[v][2]:.2f}" for v in (0.7, 2, 3.1, 7))),
        ("P22", all(sh(PI["base"][v]) >= 0.95 for v in (1, 2, 4, 7)) and all(sh(PI["base"][v]) < 0.95 for v in (0.1, 64)),
         "base STEP: " + ", ".join(f"V*={v:g} {sh(PI['base'][v]):.3f}" for v in (0.1, 1, 2, 4, 7, 64))),
        ("P23", khi is not None and hi is not None and khi >= hi + 0.5 and nmax["K10"][0] <= 0.07, f"K10 u_hi {khi} vs base {hi}; K10 null max {nmax['K10'][0]:.3f}"),
        ("P24", not cond(S2, "base") and any(cond(S2, d) for d in ("ext4", "K10", "K12")), "S2 coverage (all four conditions >= 0.95): " + ", ".join(f"{d} {cond(S2, d)}" for d in DS)),
    ]
    print("\nPredictions (droplet_range_predictions_v2.txt):")
    for name, held, detail in ch:
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")

    print("\nDecision (b): lowest-cost design (by measured mean drops) whose lowest STEP share over the span is >= 0.95 in (0.05, 0.03) and (0.05, 0.06) at s = 0.1, 0.15,"
          " and whose highest null STEP is <= 0.07")
    for n, sp in (("S1", S1), ("S2", S2)):
        ok = [d for d in sorted(DS, key=meand.get) if cond(sp, d) and nmax[d][0] <= 0.07]
        print(f"  {n} {sp}: " + (f"{ok[0]} ({meand[ok[0]]:.0f} drops per run)" if ok else "no design qualifies"))
    worst = {d: {u: min(sh(SC[(d, s, e)][u]) for s in SS for e in (0.03, 0.06)) for u in U} for d in DS}
    print("  lowest STEP share over s = 0.1, 0.15 and misread 0.03, 0.06, then edges (u_lo, u_hi): " + ", ".join(f"{d} {edges(worst[d])}" for d in DS))
    print(f"\nwall time {time.time() - T0:.0f}s")


if __name__ == "__main__":
    main()
