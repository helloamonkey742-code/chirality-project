"""Plateau blind spot of the frozen D2 droplet rule (journals/droplet_prereg.md Section 2).  Simulated drops, not data.

Rule A = droplet_sim2.log's rule (fixed 3% term, c* -1.61); B = the prereg rule now (measured misread, c* -2.21, NOT SCORABLE above
10% measured); AM1 / AM2 = proposals only: B plus a plateau h in the STEP family only / in both families (droplet_range_common.py).
c*_AM1, c*_AM2 are set on a calibration seed the way c*_B was (misread_sim.py), then every number is on fresh seeds.
Stages: 1 calibration (87 counting cells + 27 floor cells x 8 conditions x 300 runs), 2 false STEP on 87 + exclusion + 27 floor cells
(8 conditions x 400), 3 plateau map h x s x 4 conditions x 600 runs, 4 physical scenarios (ceiling, one-way miss, exclusion), 5 P1-P16.
Predictions: droplet_range_predictions_v2.txt, sha256 PRED_SHA, `shasum -a 256` and `date` 2026-10-21 01:11:44 PDT, written before this
script's results run.  v2 replaces droplet_range_predictions.txt (sha256 4ddcacbe..., 00:39:04 PDT, never edited) only to lower the run counts
(load average 40-69); predictions and decision rules are unchanged.  Disclosed smoke runs before the results run (none used for a claim): scratchpad time1.py (in the predictions
file and v1), time2.py (400 runs per rule set at c* placeholders; its verdict counts showed B STEP 9 of 200 and AM1 198, AM2 173 at h = 0.3,
s = 0.1 with c* = -3.0 for AM), mkcache.py (builds the null / exclusion truth curves; the first launch lacked a __main__ guard and was
killed within a minute, no verdicts computed).

Run: DR_CACHE=<scratch dir> UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python droplet_range_a.py > droplet_range_a.log
"""
import hashlib
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

import droplet_range_common as C
from droplet_range_common import NAMES

PRED_SHA = "e81b8eb4e8342e99bc0d8c4253c6b1f233290d6005fe1e2f90f19ec10cbdf9dd"
PROCS, CALR, EVR, MAPR = 4, 300, 400, 600
M = (0.05, 0.03)                                         # main condition: pipette CV 0.05, true misread 0.03
MAPCONDS = [(0.05, 0.0), M, (0.05, 0.06), (0.15, 0.03)]
H7, S4 = (0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5), (0.05, 0.1, 0.15, 0.3)
R3, R4 = ("B", "AM1", "AM2"), ("A", "B", "AM1", "AM2")
T0 = time.time()


def cells_of(truths, conds, tag, mode="xor", excl=0.0):
    return [(n, P, p0, cv, e, mode, excl, [tag, i, j]) for i, (n, P, p0) in enumerate(truths) for j, (cv, e) in enumerate(conds)]


def run(cells, reps, rules, cstar, label):
    per = max(1, 8000 // reps)
    out = C.pmap([dict(cells=b, reps=reps, rules=rules, design="base", V0=None, cstar=cstar) for b in C.batches(cells, per)], PROCS, label)
    print(f"  [{time.time() - T0:5.0f}s] {label} done", file=sys.stderr, flush=True)
    return [o for b in out for o in b]



def cached(name, fn):
    d = os.environ.get("DR_CACHE")                        # resume aid only; seeds are fixed, so results do not depend on it
    f = Path(d) / f"a_{name}.pkl" if d else None
    if f and f.exists():
        return pickle.loads(f.read_bytes())
    r = fn()
    if f:
        f.write_bytes(pickle.dumps(r))
    return r


def pick(cal, rule, ncell, nc):
    """misread_sim.py's c*: the lowest 5th percentile of L over cells x conditions with at least 50 eligible runs."""
    return min((cal[i * nc + j]["calib"][rule][1], i, j) for i in range(ncell) for j in range(nc) if cal[i * nc + j]["calib"][rule][0] >= 50)


def sh(o, rule, v="STEP"):
    c = o["counts"][rule]
    return c[NAMES.index(v)] / c.sum()


def line(o, rule):
    c = o["counts"][rule] / o["counts"][rule].sum()
    return " ".join(f"{n} {x:.2f}" for n, x in zip(NAMES, c))


def self_check():
    C.self_check()
    z = lambda n, p: dict(calib={"B": (n, p)})
    cal = [z(60, -1.0), z(49, -9.0), z(80, -2.5), z(70, -0.5)]               # cell 1 has too few eligible runs and is skipped
    assert pick(cal, "B", 2, 2) == (-2.5, 1, 0)                              # index = cell * 2 + condition; index 1 is skipped (49 < 50)
    o = dict(counts={"B": np.array([5, 90, 3, 0, 2])})
    assert sh(o, "B") == 0.05 and sh(o, "B", "NS") == 0.02


def main():
    assert hashlib.sha256(Path(__file__).with_name("droplet_range_predictions_v2.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    null, excl, floor = C.null_cells(PROCS), C.exclusion_cells(PROCS), C.floor_cells()
    n0, n1, n2 = len(null), len(excl), len(floor)
    assert n0 == 87 and n2 == 27
    print("Model study on simulated drops, not data.  Rules A (fixed 3%, c* -1.61), B (measured misread, c* -2.21, NS above 10% measured),")
    print("AM1 / AM2 = proposals only (B + plateau h in STEP / in both families).  Protocol droplet_sim2.protocol, 25 drops x 8 levels, 128 controls.")
    print(f"cells: {n0} counting+drying nulls, {n1} exclusion_bak curves, {n2} counting-with-floor.  Runs per cell-condition: calibration {CALR},"
          f" evaluation {EVR}, plateau map {MAPR}.")
    CONDS = C.CONDS8

    # stage 1: calibration on a seed of its own
    cal_cells = null + floor
    cal = cached("cal", lambda: run(cells_of(cal_cells, CONDS, 10), CALR, R3, None, "calibration"))
    nc = len(CONDS)
    cstar = {"A": C.C_A, "B": C.C_B}
    print("\nCalibration: c* = lowest 5th percentile of L over cells x conditions (>= 50 eligible runs: scorable, range >= 0.4, z <= 1.645)")
    for r, nn, lab in (("B", n0, "B, 87 counting+drying cells (c*_B0, to compare with frozen c*_B = -2.21)"),
                       ("AM1", n0, "AM1, 87 cells"), ("AM2", n0, "AM2, 87 cells"),
                       ("AM1", n0 + n2, "AM1, 87 + 27 floor cells (used)"), ("AM2", n0 + n2, "AM2, 87 + 27 floor cells (used)")):
        p5, i, j = pick(cal, r, nn, nc)
        print(f"  c*({lab}) = {p5:.2f}, set by {cal_cells[i][0]} at CV {CONDS[j][0]}, misread {CONDS[j][1]}")
        if nn == n0 + n2:
            cstar[r] = p5
        if lab.startswith("B,"):
            cB0 = p5
    print(f"  frozen A {C.C_A}, B {C.C_B}; amended c*: AM1 {cstar['AM1']:.2f}, AM2 {cstar['AM2']:.2f}")

    # stage 2: false STEP on fresh seeds
    allnull = null + excl + floor
    grp = [("count", range(0, n0)), ("excl", range(n0, n0 + n1)), ("floor", range(n0 + n1, n0 + n1 + n2))]
    ev = cached("ev", lambda: run(cells_of(allnull, CONDS, 11), EVR, R3, cstar, "false STEP"))
    top = {}
    print("\nFalse STEP: highest STEP share over the cells of each group (limit 0.07), by condition")
    print("  CV   misread | " + " | ".join(f"{r:>3} " + " ".join(f"{g:>5}" for g, _ in grp) for r in R3))
    for j, (cv, e) in enumerate(CONDS):
        row = []
        for r in R3:
            vals = []
            for g, rg in grp:
                s, i = max((sh(ev[i * nc + j], r), i) for i in rg)
                vals.append(s)
                if s > top.get((r, g), (-1,))[0]:
                    top[(r, g)] = (s, allnull[i][0], (cv, e))
            row.append("    " + " ".join(f"{v:5.3f}" for v in vals))
        print(f"  {cv:.2f} {e:.2f}    |" + " |".join(row))
    for (r, g), (s, name, c) in top.items():
        print(f"  {r} {g}: max {s:.3f} ({name}, CV {c[0]}, misread {c[1]})")
    fmax = {r: max(top[(r, g)][0] for g, _ in grp) for r in R3}

    # stage 3: plateau map
    truths = [(f"h={h} s={s}", C.probit(h, s), 0 * C.LF) for h in H7 for s in S4]
    pm = cached("map", lambda: run(cells_of(truths, MAPCONDS, 12), MAPR, R4, cstar, "plateau map"))
    MAP = {(r, c, h, s): pm[(hi * len(S4) + si) * len(MAPCONDS) + ci] for hi, h in enumerate(H7) for si, s in enumerate(S4)
           for ci, c in enumerate(MAPCONDS) for r in R4}
    print(f"\nPlateau map: truth P = h + (1 - h) Phi(-ln(V / 3.1 uL) / s), {MAPR} runs per cell.  Verdict shares, main condition (CV 0.05, misread 0.03)")
    for r in R4:
        print(f"rule {r}")
        for h in H7:
            for s in S4:
                print(f"  h={h} s={s} " + line(MAP[(r, M, h, s)], r))
    print("\nSTEP share by plateau h (rows) and s = " + ", ".join(map(str, S4)) + " (columns), every rule and condition")
    for c in MAPCONDS:
        for r in R4:
            print(f"  CV {c[0]} misread {c[1]} rule {r}: " + " | ".join(f"h={h} " + " ".join(f"{sh(MAP[(r, c, h, s)], r):.2f}" for s in S4) for h in H7))
    print("\nNOT SCORABLE share, rule B (rows as above)")
    for c in MAPCONDS:
        print(f"  CV {c[0]} misread {c[1]}: " + " | ".join(f"h={h} " + " ".join(f"{sh(MAP[('B', c, h, s)], 'B', 'NS'):.2f}" for s in S4) for h in H7))

    # stage 4: scenarios, main condition unless stated
    sc = []
    for c in (0.9, 0.8, 0.7):
        for s in (0.1, 0.15):
            sc.append((f"ceiling c={c} s={s}", C.probit(0, s, c), 0 * C.LF, M[0], M[1], "xor", 0.0))
    for m in (0.06, 0.12, 0.20):
        sc.append((f"one-way miss m={m} s=0.1", C.probit(0, 0.1), 0 * C.LF, 0.05, m, "miss", 0.0))
    for u in (0.5, 0.8, 0.9):
        sc.append((f"exclusion u={u} s=0.15", C.probit(0, 0.15), 0 * C.LF, M[0], M[1], "xor", u))
    scc = [(n, P, p0, cv, e, mode, x, [13, i, 0]) for i, (n, P, p0, cv, e, mode, x) in enumerate(sc)]
    so = cached("sc", lambda: run(scc, MAPR, R4, cstar, "scenarios"))
    SC = {c[0]: o for c, o in zip(sc, so)}
    print("\nScenarios (CV 0.05; misread 0.03 except the one-way miss, where controls read wrong at about m / 2), verdict shares")
    for r in R4:
        print(f"rule {r}")
        for c in sc:
            print(f"  {c[0]:26} " + line(SC[c[0]], r))

    # stage 5: predictions
    st = lambda r, c, h, s: sh(MAP[(r, c, h, s)], r)
    am = {r: dict(P12=fmax[r] <= 0.07, P14=all(st(r, M, h, s) >= 0.8 for h in (0.1, 0.2, 0.3) for s in (0.1, 0.15)),
                  P15=all(st(r, M, 0, s) >= 0.95 for s in (0.05, 0.1, 0.15)) and st(r, M, 0, 0.3) >= st("B", M, 0, 0.3) - 0.15) for r in ("AM1", "AM2")}
    flo = top[("B", "floor")]
    ch = [
        ("P1", 0.02 <= st("A", M, 0.3, 0.15) <= 0.09, f"rule A STEP at h=0.3 s=0.15: {st('A', M, 0.3, 0.15):.3f}"),
        ("P2", all(st("B", M, 0, s) >= 0.95 for s in (0.05, 0.1, 0.15)), "B STEP at h=0: " + ", ".join(f"s={s} {st('B', M, 0, s):.3f}" for s in (0.05, 0.1, 0.15))),
        ("P3", all(st("B", c, 0.3, s) <= 0.15 for c in MAPCONDS for s in (0.05, 0.1, 0.15)),
         "B STEP at h=0.3: " + "; ".join(f"CV {c[0]}/{c[1]}: " + " ".join(f"{st('B', c, 0.3, s):.2f}" for s in (0.05, 0.1, 0.15)) for c in MAPCONDS)),
        ("P4", all(max(np.diff([st("B", M, h, s) for h in H7])) <= 0.03 for s in (0.1, 0.15)),
         "B STEP vs h, main: " + "; ".join(f"s={s}: " + " ".join(f"{st('B', M, h, s):.2f}" for h in H7) for s in (0.1, 0.15))),
        ("P5", all(st("B", M, 0.1, s) >= 0.8 and 0.1 <= st("B", M, 0.2, s) <= 0.8 for s in (0.1, 0.15)),
         "B STEP h=0.1 / h=0.2: " + "; ".join(f"s={s}: {st('B', M, 0.1, s):.2f} / {st('B', M, 0.2, s):.2f}" for s in (0.1, 0.15))),
        ("P6", all(st("B", M, 0.5, s) <= 0.15 for s in (0.1, 0.15)), "B STEP at h=0.5: " + ", ".join(f"s={s} {st('B', M, 0.5, s):.3f}" for s in (0.1, 0.15))),
        ("P7", all(sh(SC[f"ceiling c=0.9 s={s}"], "B") >= 0.8 for s in (0.1, 0.15)),
         "B STEP at ceiling 0.9: " + ", ".join(f"s={s} {sh(SC[f'ceiling c=0.9 s={s}'], 'B'):.3f}" for s in (0.1, 0.15))),
        ("P8", all(sh(SC[f"one-way miss m={m} s=0.1"], "B") >= 0.85 for m in (0.06, 0.12)) and sh(SC["one-way miss m=0.2 s=0.1"], "B", "NS") >= 0.3,
         "B STEP m=0.06 " + f"{sh(SC['one-way miss m=0.06 s=0.1'], 'B'):.3f}, m=0.12 {sh(SC['one-way miss m=0.12 s=0.1'], 'B'):.3f}; NS at m=0.20 {sh(SC['one-way miss m=0.2 s=0.1'], 'B', 'NS'):.3f}"),
        ("P9", all(sh(SC[f"exclusion u={u} s=0.15"], "B") >= 0.9 for u in (0.5, 0.8, 0.9)),
         "B STEP: " + ", ".join(f"u={u} {sh(SC[f'exclusion u={u} s=0.15'], 'B'):.3f}" for u in (0.5, 0.8, 0.9))),
        ("P10", -2.9 <= cB0 <= -1.7, f"c*_B0 = {cB0:.2f}"),
        ("P11", -5.0 <= cstar["AM2"] <= -1.5, f"c*_AM2 = {cstar['AM2']:.2f}"),
        ("P12", am["AM2"]["P12"], f"AM2 max false STEP {fmax['AM2']:.3f} (AM1 {fmax['AM1']:.3f}, B {fmax['B']:.3f})"),
        ("P13", flo[0] <= 0.07, f"B max STEP on floor cells {flo[0]:.3f} ({flo[1]}, CV {flo[2][0]}, misread {flo[2][1]})"),
        ("P14", am["AM2"]["P14"], "AM2 STEP main: " + "; ".join(f"h={h}: " + " ".join(f"{st('AM2', M, h, s):.2f}" for s in (0.1, 0.15)) for h in (0.1, 0.2, 0.3))),
        ("P15", am["AM2"]["P15"], "AM2 STEP h=0: " + " ".join(f"s={s} {st('AM2', M, 0, s):.3f}" for s in S4) + f"; B at s=0.3 {st('B', M, 0, 0.3):.3f}"),
        ("P16", st("AM2", M, 0, 0.3) - st("AM1", M, 0, 0.3) >= 0.05, f"AM2 - AM1 STEP at h=0 s=0.3: {st('AM2', M, 0, 0.3):.3f} - {st('AM1', M, 0, 0.3):.3f}"),
    ]
    print("\nPredictions (droplet_range_predictions_v2.txt):")
    for name, held, detail in ch:
        print(f"  {name} {'HELD' if held else 'FAILED'}: {detail}")
    matters = st("B", M, 0.2, 0.1) < 0.8
    held = {n: h for n, h, _ in ch}
    print(f"\nDecision (a): blind spot matters = B STEP(h=0.2, s=0.1) < 0.80 in main: {st('B', M, 0.2, 0.1):.3f} -> {matters}")
    for r in ("AM1", "AM2"):
        print(f"  {r}: P12 {am[r]['P12']}, P14 {am[r]['P14']}, P15 {am[r]['P15']}")
    pick_ = ("no amendment" if not matters else "AM2" if all(am["AM2"].values()) and held["P16"] else
             "AM1 (P16 failed)" if all(am["AM2"].values()) and all(am["AM1"].values()) else
             "AM2 (P16 failed, AM1 fails one of P12, P14, P15)" if all(am["AM2"].values()) else "limit sentence only")
    print(f"  proposal for the notes (not applied): {pick_}")
    print(f"\nwall time {time.time() - T0:.0f}s")


if __name__ == "__main__":
    main()
