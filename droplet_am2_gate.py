"""Three open points of the NaClO3 droplet plan, settled on SIMULATED drops (not data).  Nothing in journals/droplet_prereg.md or any existing file is edited.
(1) Is the failed gate P7 (droplet_ctrl_am2.log:163, 56 of 60) caused by -2.21 being a rounded c*?  misread_sim.py's calibration is reproduced exactly (87 counting+drying
    cells x 12 conditions, 1000 runs, seeds [4, i, j], 192 control draws after the protocol draws, first 128 used), c*_B(128) printed unrounded, and the 60-cell gate rerun with it.
(2) What c* do B and AM2 need at 480 controls?  Same simulated runs; a second block of 480 control draws per run (drawn after the 192-block).  c*_AM2(128) by the same recipe too.
(3) How often does each rule reject a sharp step (sharp > 3.84) when the truth is smooth, on Bak 90 C-like data (n = 3, 6, 10 at the unprinted markers, misread 0-10%)?  The paired
    number is the same rate when the truth is a sharp step.
Predictions: droplet_am2_gate_predictions.txt, sha256 PRED_SHA, written 2026-10-21 PDT before any result of this script.
Disclosed: no timing or cache build was run for this unit.  It reuses droplet_ctrl_am2's copy of the 87 null curves (DR_CACHE, results do not depend on it) and imports
droplet_ctrl_am2 only for its patches (RULES Bn, controls stream), read_log() and the constants.  Per-task results go to AM2_CK (scratchpad) as they finish and a task is skipped when its file
exists, so a run stopped at the BUDGET (default 2400 s of new tasks, inside the 45 min budget) is resumed by running the same command again; every task is a pure function of its seed,
and only a run that finds every task done prints the log.  Run:
cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_am2_gate.py > droplet_am2_gate.log 2> droplet_am2_gate.err
"""
import os

SCR = "/private/tmp/claude-501/-Users-saptarshighosh-Documents-chirality-project/05c76bba-5069-4c19-a20f-58eb08337c55/scratchpad"
for k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):     # another agent shares the CPU: one thread per process
    os.environ.setdefault(k, "1")
os.environ.setdefault("DR_CACHE", SCR + "/ctrl_am2/cache")

import hashlib
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import droplet_ctrl_am2 as P                    # patches RULES (Bn, AM2n) and C.controls at import, in workers too
import droplet_range_common as C

np.seterr(over="ignore", divide="ignore", invalid="ignore")
PRED_SHA = "db59504ecfae5b78f71118b6d306aecd21b788beaf3bc993e01eaf0878251b0e"
HERE = Path(__file__).parent
CK = Path(os.environ.get("AM2_CK", SCR + "/am2_gate"))
BUDGET = float(os.environ.get("AM2_BUDGET", 2400))
PROCS, CAL, REPS, BREPS = 4, 1000, 500, 500
CONDS, SVALS, V, Z = P.CONDS, P.SVALS, P.V, P.Z
NCS, RL = (128, 480), ("B", "AM2")
EV, NB = (0.0, 0.02, 0.04, 0.06, 0.08, 0.10), (3, 6, 10)
C.RULES["Bu"] = ("meas", ("c0", "s0"), False)                          # B with the unrounded c*_B(128), no control cap (the log's policy B)
MISS = {("cnt", 0.05, 0.06), ("cnt", 0.05, 0.08), ("cnt", 0.05, 0.10), (0.15, 0.05, 0.08)}        # droplet_ctrl_am2.log:67, 68, 69, 98
S1 = {0.0: (0.97, 1.0), 0.02: (0.88, 1.0), 0.04: (0.75, 1.0), 0.06: (0.55, 0.95), 0.08: (0.40, 0.85), 0.10: (0.25, 0.75)}
T0 = time.time()
RES = []


def say(s=""):
    print(s)


def held(name, ok, detail):
    RES.append((name, bool(ok)))
    say(f"  {name} {'HELD' if ok else 'FAILED'}: {detail}")


def ckpt(tag, i, deadline, compute):
    """Run compute() and save it as CK/tag_i.npy, unless that file exists or the deadline for starting new tasks has passed."""
    f = CK / f"{tag}_{i}.npy"
    if f.exists() or time.time() > deadline:
        return
    tmp = CK / f"{tag}_{i}.tmp.npy"
    np.save(tmp, compute())
    os.replace(tmp, f)


def stage(fn, tasks, tag):
    """tasks = [(i, ...)].  Runs those without a file; exits (no log) if any is still missing; else returns the arrays in order."""
    todo = [t for t in tasks if not (CK / f"{tag}_{t[0]}.npy").exists()]
    if todo:
        with Pool(PROCS) as pool:
            for n, _ in enumerate(pool.imap_unordered(fn, todo)):
                print(f"  {tag} {n + 1}/{len(todo)}  {time.time() - T0:.0f}s", file=sys.stderr, flush=True)
    left = [t[0] for t in tasks if not (CK / f"{tag}_{t[0]}.npy").exists()]
    if left:
        print(f"[{tag}: {len(left)} of {len(tasks)} tasks still to do at {time.time() - T0:.0f}s; run the same command again]", file=sys.stderr, flush=True)
        sys.exit(0)
    return [np.load(CK / f"{tag}_{t[0]}.npy") for t in tasks]


# ---------------------------------------------------------------- (1), (2): calibration
def cell_runs(P_, p0, cv, e, reps, seed):
    """misread_sim.run_cell's stream: protocol draws, then 192 control draws (first 128 used), then 480 more (new here)."""
    rng = np.random.default_rng(seed)
    KN = np.array([C.protocol2(rng, P_, p0, e, cv, exact=True)[1] for _ in range(reps)])
    w = rng.random((reps, 192)) < e
    w480 = rng.random((reps, 480)) < e
    return KN[:, 0], KN[:, 1], w[:, :128].sum(1), w480.sum(1)


def cal_task(t):
    i, P_, p0, deadline = t

    def compute():
        F = C.fitter("base", True)
        runs = [cell_runs(P_, p0, cv, e, CAL, [4, i, j]) for j, (cv, e) in enumerate(CONDS)]
        K, N = np.vstack([r[0] for r in runs]), np.vstack([r[1] for r in runs])
        z, sc = C.trend_z(K, N, F.x), C.scorable_v(K, N)
        out = np.full((2, 2, len(CONDS), 2), np.nan)                 # [n_c][rule][condition] = (eligible runs, 5th percentile of L)
        for a, (nc, col) in enumerate(zip(NCS, (2, 3))):
            kc = np.concatenate([r[col] for r in runs])
            fits = F.fit(K, N, (kc + 0.5) / (nc + 1), ("c0", "cH", "s0", "sH"))
            for b, r in enumerate(RL):
                st = C.rule_stats(fits, C.RULES[r][1])
                for j in range(len(CONDS)):
                    s = slice(j * CAL, (j + 1) * CAL)
                    ok = sc[s] & (st["rng"][s] >= 0.4) & (z[s] <= 1.645)
                    out[a, b, j] = ok.sum(), np.percentile(st["L"][s][ok], 5) if ok.sum() >= 50 else np.nan
        return out
    ckpt("cal", i, deadline, compute)


def gate_task(t):
    i, task, deadline = t
    ckpt("gate", i, deadline, lambda: np.array([[o["counts"][r] for r in task["rules"]] for o in C.run_batch(task)]))      # (12, 2 rules, 5)


def cstars(cal, names):
    """cal (87, 2, 2, 12, 2) -> {(rule, n_c): (c*, (cell, cv, misread), [3 lowest])}."""
    out = {}
    for a, nc in enumerate(NCS):
        for b, r in enumerate(RL):
            p5 = cal[:, a, b, :, 1]
            low = sorted((p5[i, j], names[i].strip(), *CONDS[j]) for i in range(87) for j in range(12) if not np.isnan(p5[i, j]))[:3]
            out[(r, nc)] = (low[0][0], low[0][1:], low)
    return out


# ---------------------------------------------------------------- (1): the 60-cell gate
def gate_stats(G, r, log):
    """G (91, 12, 2, 5) run counts with rules ('Bn', 'Bu').  The 60 cells of P7: counting maximum + probit s = 0.1..0.3, 12 conditions each."""
    c = G[:, :, ("Bn", "Bu").index(r), :]
    stp = c[:, :, V]
    sh = {"cnt": stp[:87].max(0) / REPS, **{s: stp[87 + k] / REPS for k, s in enumerate(SVALS)}}
    dig, tol, miss, rows = 0, True, [], {}
    for key in ("cnt", *SVALS):
        lv, _, dg = log[key]
        d = np.abs(np.round(sh[key], dg) - lv) < 1e-9
        dig += int(d.sum())
        tol &= bool((np.abs(sh[key] - lv) <= (0.004 if dg == 3 else 0.015) + 1e-9).all())
        miss += [(key, *CONDS[j]) for j in np.flatnonzero(~d)]
        rows[key] = sh[key]
    ns = c[88, :, Z] / REPS
    return dict(dig=dig, tol=tol, miss=miss, cmax=sh["cnt"].max(), ns=ns, ns0=bool((np.round(ns, 2) == 0).all()), rows=rows)


# ---------------------------------------------------------------- (3): Bak 90 C-like smooth and sharp truths
def llmax(F, X, e, g, mask=None):
    W = F.W(e, g)
    W = W if mask is None else W[:, mask]
    X = X.astype(W.dtype)
    return np.concatenate([(X[i:i + 100] @ W).max(1) for i in range(0, len(X), 100)])


def sharp_stats(F, K, N, e, ps0, psH):
    """droplet_ctrl_am2.stage3's 'sharp' per run: 2 (best ll over the rule's groups - best ll among its STEP curves with s <= 0.15)."""
    X = np.hstack([K, N])
    ll = {g: llmax(F, X, e, g) for g in ("c0", "cH", "s0", "sH")}
    s0, sH = llmax(F, X, e, "s0", ps0 <= 0.15), llmax(F, X, e, "sH", psH <= 0.15)
    return {"B": 2 * (np.maximum(ll["c0"], ll["s0"]) - s0), "AM2": 2 * (np.max([ll[g] for g in ll], 0) - np.maximum(s0, sH))}


_BAK = {}


def bak_setup(n1):
    if n1 not in _BAK:
        from blind_units import BAK_TINC90
        from droplet_sim2 import curves
        from exclusion_bak import bak
        x, k, n = bak(BAK_TINC90, n1)
        F = C.Fitter(x, True)
        ps0 = curves(x)["step"][3]
        psH = np.tile(np.repeat(C.SG, len(C.CGC)), len(C.HG))
        assert len(ps0) == F.P["s0"].shape[0] and len(psH) == F.P["sH"].shape[0]
        _BAK[n1] = (F, k, n, ps0, psH)
    return _BAK[n1]


def bak_task(t):
    a, n1, deadline = t

    def compute():
        F, k, n, ps0, psH = bak_setup(n1)
        X = np.hstack([k, n])[None]
        out = np.zeros((2, 2, len(EV)))                              # [truth: smooth, sharp][rule][misread] = share with sharp > 3.84
        for b, e in enumerate(EV):
            pt = (F.P["c0"][np.argmax(X.astype(np.float64) @ F.W(e, "c0"))], F.P["s0"][np.argmax(np.where(ps0 <= 0.15, X @ F.W(e, "s0"), -np.inf))])
            for tr in range(2):
                rng = np.random.default_rng([31, a, b, tr])
                Ks = rng.binomial(n.astype(int), e + (1 - 2 * e) * pt[tr], (BREPS, len(n))).astype(float)
                sh = sharp_stats(F, Ks, np.tile(n, (BREPS, 1)), e, ps0, psH)
                out[tr, :, b] = [(sh[r] > 3.84).mean() for r in RL]
        return out
    ckpt("bak", a, deadline, compute)


# ---------------------------------------------------------------- self-checks (silent; they never print a prediction quantity)
def self_check():
    C.self_check()
    import misread_sim as M
    nulls = C.null_cells(PROCS)
    assert len(nulls) == 87
    F = C.fitter("base", True)
    for i, j in ((3, 5), (60, 4)):                                   # my calibration pipeline against misread_sim.run_cell on the same stream
        _, P_, p0 = nulls[i]
        cv, e = CONDS[j]
        ref = M.run_cell((P_, p0, cv, e, 40, [4, i, j], False))[("B", 128)]
        K, N, kc, _ = cell_runs(P_, p0, cv, e, 40, [4, i, j])
        st = C.rule_stats(F.fit(K, N, (kc + 0.5) / 129, ("c0", "s0")), ("c0", "s0"))
        assert np.allclose(st["L"], ref["L"], atol=1e-8) and ((C.scorable_v(K, N) & (st["rng"] >= 0.4) & (C.trend_z(K, N, F.x) <= 1.645)) == ((~ref["ns"]) & (ref["rng"] >= 0.4) & (ref["z"] <= 1.645))).all()
    rx = r"e (\d\.\d\d): L ([+-]\d+\.\d\d) sharp\s+([\d.]+) (\w+)"       # sharp on the real Bak 90 C series against the two earlier logs (36 values)
    ref = {"B": [re.findall(rx, l) for l in (HERE / "bak_policy_b.log").read_text().splitlines() if l.startswith("Tinc")],
           "AM2": [re.findall(rx, l) for l in (HERE / "droplet_ctrl_am2.log").read_text().splitlines() if l.startswith("Tinc") and " AM2 |" in l]}
    for a, n1 in enumerate(NB):
        F, k, n, ps0, psH = bak_setup(n1)
        for b, e in enumerate(EV):
            sh = sharp_stats(F, k[None], n[None], e, ps0, psH)
            for r in RL:
                assert abs(sh[r][0] - float(ref[r][a][b][2])) <= 0.051, (r, n1, e, sh[r][0], ref[r][a][b])
    print(f"[self-check passed {time.time() - T0:.0f}s]", file=sys.stderr, flush=True)


# ---------------------------------------------------------------- main
def main():
    assert hashlib.sha256((HERE / "droplet_am2_gate_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    CK.mkdir(parents=True, exist_ok=True)
    self_check()
    dl = time.time() + BUDGET
    nulls = C.null_cells(PROCS)
    cal = np.stack(stage(cal_task, [(i, P_, p0, dl) for i, (_, P_, p0) in enumerate(nulls)], "cal"))
    cs = cstars(cal, [n for n, _, _ in nulls])
    cB, cA, cB4, cA4 = cs[("B", 128)][0], cs[("AM2", 128)][0], cs[("B", 480)][0], cs[("AM2", 480)][0]
    steps = [(f"probit s={s}", C.probit(0, s), 0 * C.LF) for s in SVALS]
    cstar = {"Bn": -2.21, "Bu": cB}
    G = np.stack(stage(gate_task, [(i, dict(cells=[(n, P_, p0, cv, e, "xor", 0.0, [5, i, j]) for j, (cv, e) in enumerate(CONDS)], reps=REPS, rules=("Bn", "Bu"), design="base",
                                            V0=None, cstar=cstar, want_v0=False, exact=True), dl) for i, (n, P_, p0) in enumerate(nulls + steps)], "gate"))
    bk = np.stack(stage(bak_task, [(a, n1, dl) for a, n1 in enumerate(NB)], "bak"))      # (3, truth, rule, e)
    log = P.read_log()
    gn, gu = gate_stats(G, "Bn", log), gate_stats(G, "Bu", log)

    say(f"Unit droplet_am2_gate.  Model study on simulated drops, not data.  Predictions sha256 {PRED_SHA}.  Self-checks passed (calibration pipeline = misread_sim.run_cell on 2 cells;")
    say("sharp on the real Bak 90 C series = bak_policy_b.log and droplet_ctrl_am2.log for B and AM2 at all 18 (n, misread) pairs).")
    say(f"\nPART 1 and 2.  Calibration: {len(nulls)} counting+drying cells x {len(CONDS)} conditions, {CAL} runs per cell, seeds [4, i, j]; c* = lowest 5th percentile of L over cells with >= 50 eligible runs.")
    say("  rule  n_c    c* (6 decimals)    set by (cell, CV, misread)                          next two lowest cells")
    for (r, nc), (v, w, low) in cs.items():
        say(f"  {r:<4} {nc:>4}   {v:+.6f}        {w[0]:<22} CV {w[1]:.2f} misread {w[2]:.2f}   " + "; ".join(f"{x[0]:+.4f} {x[1]} {x[2]:.2f}/{x[3]:.2f}" for x in low[1:]))
    say(f"  misread_sim.log:4-6 print c*_B(64/128/192) = -2.61 / -2.21 / -2.03; droplet_ctrl_am2 used -2.21 (B) and -2.40 (AM2, droplet_range_a.log: 300 runs, 8 conditions, 87 + 27 floor cells).")
    say(f"\nGATE.  The 60 cells of P7 ({REPS} runs, seeds [5, i, j], 128 controls, no control cap).  Bn = c* -2.21 (the earlier unit's value); Bu = c* {cB:+.6f} (this unit's calibration), same simulated runs.")
    for key, ttl in (("cnt", "Counting cells: highest STEP share over the 87 cells"), *((s, f"probit s={s}: STEP share") for s in SVALS)):
        lv, line, dg = log[key]
        fm = "{:5.3f}" if dg == 3 else "{:5.2f}"
        say(f"\n{ttl}  (log column: misread_sim.log:{line}-{line + 11}, B 128)")
        say("  CV   misread  log B128        Bn        Bu    Bn-log    Bu-log")
        for j, (cv, e) in enumerate(CONDS):
            a, b = gn["rows"][key][j], gu["rows"][key][j]
            say(f"  {cv:.2f} {e:.2f}    " + "".join(f"{x:>9}" for x in (fm.format(lv[j]), fm.format(a), fm.format(b), f"{a - lv[j]:+.3f}", f"{b - lv[j]:+.3f}")))
    lv, line, _ = log["ns"]
    say(f"\nprobit s=0.15: NOT SCORABLE share (log column misread_sim.log:{line}-{line + 11}): " + " ".join(f"{x:.3f}" for x in lv))
    say("  Bn: " + " ".join(f"{x:.3f}" for x in gn["ns"]) + "\n  Bu: " + " ".join(f"{x:.3f}" for x in gu["ns"]))
    for nm, g in (("Bn", gn), ("Bu", gu)):
        say(f"{nm}: {g['dig']}/60 cells equal the printed digits; all within tolerance {g['tol']}; counting maximum {g['cmax']:.3f}; NS share all 0.00 {g['ns0']}; cells that differ: {g['miss']}")

    for tr, ttl in ((0, "SMOOTH truth (best COUNT curve fitted to the series)"), (1, "SHARP truth (best STEP curve with s <= 0.15 fitted to the series)")):
        say(f"\nPART 3.  {ttl}: share of {BREPS} runs with sharp > 3.84 (SE <= 0.022).  Bak 90 C series, n = 3, 6 or 10 at the seven P = 1 markers, misread e known, seeds [31, n idx, e idx, {tr}].")
        say("  n    rule " + "".join(f"  e={e:.2f}" for e in EV))
        for a, n1 in enumerate(NB):
            for b, r in enumerate(RL):
                say(f"  {n1:<4} {r:<4} " + "".join(f"  {bk[a, tr, b, k]:6.3f}" for k in range(len(EV))))
    sm, sp = bk[:, 0], bk[:, 1]                                       # (3, rule, 6)
    lo, hi = slice(0, 3), slice(3, 6)
    say("\nPaired number (mean over n = 3, 6, 10):  rule   smooth truth, e <= 0.04 | e >= 0.06      sharp truth, e <= 0.04 | e >= 0.06 | max")
    for b, r in enumerate(RL):
        say(f"                                          {r:<4}   {sm[:, b, lo].mean():.3f} | {sm[:, b, hi].mean():.3f}                      {sp[:, b, lo].mean():.3f} | {sp[:, b, hi].mean():.3f} | {sp[:, b].max():.3f}")

    # ------------------------------------------------------------ scoring
    say("\nPredictions (scored as written in droplet_am2_gate_predictions.txt)")
    w = cs[("B", 128)][1]
    cell_ok = w[0] == "dry a=1.67" and (w[1], w[2]) == (0.05, 0.04)
    g1w = -2.2150 <= cB <= -2.2050 and cell_ok
    held("G1", -2.2100 <= cB <= -2.2050 and cell_ok, f"c*_B(128) = {cB:+.6f}, set by {w[0]} at CV {w[1]}, misread {w[2]} (weak form [-2.2150, -2.2050] and same cell: {'HELD' if g1w else 'FAILED'})")
    core = gu["dig"] >= 58 and gu["tol"] and abs(gu["cmax"] - 0.054) <= 0.004 + 1e-9 and gu["ns0"] and not (MISS & set(gu["miss"]))
    ctrl = gn["dig"] == 56 and set(gn["miss"]) == MISS
    held("G2", core and ctrl, f"Bu {gu['dig']}/60 equal digits (misses {gu['miss']}), tolerance {gu['tol']}, counting max {gu['cmax']:.3f}, NS 0.00 {gu['ns0']}; control Bn {gn['dig']}/60 with misses {gn['miss']} (must be 56 and the four of droplet_ctrl_am2.log:67-69, 98)")
    k1 = -2.50 <= cA <= -2.15 and abs(cA + 2.40) <= 0.10 + 1e-9
    held("K1", k1, f"c*_AM2(128) = {cA:+.6f} (band [-2.50, -2.15]: {-2.50 <= cA <= -2.15}; within 0.10 of -2.40: {abs(cA + 2.40) <= 0.10 + 1e-9}; point -2.31)")
    held("K2", -1.84 <= cB4 <= -1.64, f"c*_B(480) = {cB4:+.6f} (band [-1.84, -1.64], point -1.74)")
    held("K3", -1.95 <= cA4 <= -1.65 and -0.20 <= cA4 - cB4 <= 0.05, f"c*_AM2(480) = {cA4:+.6f} (band [-1.95, -1.65], point -1.80); AM2 - B at 480 = {cA4 - cB4:+.4f} (band [-0.20, +0.05])")
    inband = np.array([[S1[e][0] <= sm[a, 0, k] <= S1[e][1] for k, e in enumerate(EV)] for a in range(3)])
    rise = np.diff(sm[:, 0], axis=1).max()
    held("S1", inband.sum() >= 15 and rise <= 0.04 + 1e-9, f"B smooth truth: {inband.sum()}/18 cells in band (cells outside: {[(NB[a], EV[k], round(float(sm[a, 0, k]), 3)) for a, k in zip(*np.where(~inband))]}); largest rise with e {rise:+.3f}")
    held("S2", (sm[:, 1] <= 0.65).all() and (sm[:, 1, hi] <= 0.35).all(), f"AM2 smooth truth: max over 18 cells {sm[:, 1].max():.3f} (<= 0.65), max over the nine cells with e >= 0.06 {sm[:, 1, hi].max():.3f} (<= 0.35)")
    held("S3", (sm[:, 1] <= sm[:, 0] + 0.03 + 1e-9).all() and ((sm[:, 0] - sm[:, 1])[:, lo] >= 0.30 - 1e-9).all(),
         f"max (AM2 - B) {np.max(sm[:, 1] - sm[:, 0]):+.3f} (<= +0.03); min (B - AM2) over the nine cells with e <= 0.04: {np.min((sm[:, 0] - sm[:, 1])[:, lo]):+.3f} (>= 0.30)")
    held("S4", (sp <= 0.15 + 1e-9).all(), f"sharp truth: max share B {sp[:, 0].max():.3f}, AM2 {sp[:, 1].max():.3f} (<= 0.15)")
    say(f"\nHELD {sum(v for _, v in RES)} of {len(RES)}: " + ", ".join(f"{n} {'HELD' if v else 'FAILED'}" for n, v in RES))

    say("\nDecision rules")
    if not g1w:
        say("  D1: G1's weak form FAILED -> my calibration does not reproduce misread_sim; nothing in this unit about c* is usable.")
    elif core:
        say(f"  D1: G1 (weak form) and the gate hold with the unrounded c* ({gu['dig']}/60) -> the rounding of c* to two decimals explains the misses of P7; R1 NOT TRUSTED is lifted for the B column (Bn control reproduced: {ctrl}).  "
            + (f"c*_AM2(128) = {cA:+.4f} is within 0.10 of -2.40, so droplet_ctrl_am2's AM2 columns stand." if abs(cA + 2.40) <= 0.10 + 1e-9 else
               f"c*_AM2(128) = {cA:+.4f} is farther than 0.10 from -2.40 -> droplet_ctrl_am2's AM2 columns must be rerun with the recalibrated c* before being quoted.")
            + "  A prereg that prints power tables to the digit must quote c* to four decimals.")
    else:
        say(f"  D1: G1 holds but the gate with the unrounded c* has {gu['dig']}/60 (needs >= 58 and the conditions of G2) -> rounding is NOT the (whole) cause; the AM2-versus-B comparison stays NOT TRUSTED.")
    d1, d2 = abs(cB4 - cB), abs(cA4 - cB4)
    say(f"  D2: |c*_B(480) - c*_B(128)| = {d1:.4f} -> " + ("> 0.25: W1 (480 controls) cannot be adopted without a c* tied to the control count; the prereg must say which c* goes with which count" if d1 > 0.25
                                                          else "<= 0.25: a c* tied to the control count is not forced by this calibration")
        + f".  |c*_AM2(480) - c*_B(480)| = {d2:.4f} -> " + ("<= 0.05: one threshold could serve both rules at 480" if d2 <= 0.05 else "> 0.05: AM2 needs its own c*")
        + f".  c*_B(480) = {cB4:+.4f}, c*_AM2(480) = {cA4:+.4f}; calibration values only (no false STEP or power evaluated at 480).")
    PA, PB = sm[:, 1], sm[:, 0]
    if (PA <= 0.50).sum() >= 15 and (PB[:, lo] >= 0.75).sum() >= 7:
        say(f"  D3: P_A <= 0.50 in {(PA <= 0.50).sum()}/18 cells and P_B >= 0.75 in {(PB[:, lo] >= 0.75).sum()}/9 low-misread cells -> AM2 gives up most of the ability to reject a sharp step in smooth Bak-like data; adopting AM2 means the prereg "
            "can no longer claim that a SMOOTH label rejects a sharp step (only 'smooth, or a step with a plateau'), to be weighed against B's plateau blind spot.")
    elif (PA >= 0.50).sum() >= 9:
        say(f"  D3: P_A >= 0.50 in {(PA >= 0.50).sum()}/18 cells -> AM2 keeps much of the power; the single-dataset values understated it.")
    else:
        say(f"  D3: mixed (P_A <= 0.50 in {(PA <= 0.50).sum()}/18 cells, P_B >= 0.75 in {(PB[:, lo] >= 0.75).sum()}/9 low-misread cells, P_A >= 0.50 in {(PA >= 0.50).sum()}/18): report the table, propose no sentence.")
    big = [(r, NB[a], EV[k], round(float(sp[a, b, k]), 3)) for b, r in enumerate(RL) for a in range(3) for k in range(6) if sp[a, b, k] > 0.15]
    say("  D4: " + (f"sharp-truth share above 0.15 in {len(big)} cells {big} -> that rule's 'sharp > 3.84 rejects a sharp step' is not a 5%-level statement; any prereg sentence using it must quote the simulated rate." if big
                    else "no sharp-truth share above 0.15 for either rule."))
    say("  D5: nothing is edited in journals/droplet_prereg.md or any existing file; all wording is marked NOT APPLIED.")
    print(f"[done {time.time() - T0:.0f}s]", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
