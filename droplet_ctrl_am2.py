"""Control-drop count and rule AM2 against rule B for the NaClO3 droplet test (journals/droplet_prereg.md Sections 2, 3, 6, 10).
Part 1 is exact binomial arithmetic, Part 2 a model study on simulated drops, not data.  Nothing in the prereg is edited.

Part 1: the 10% control cap (NOT SCORABLE if k/n > 0.10) over n = 128..1024 controls: P(NS) at true misread 6/8/10%, and P(false SCORABLE) = P(k/n <= 0.10)
at true 11-14%; lower caps tried for the pair NS(6%) <= 0.05 and FS(11%) <= 0.10.
Part 2: misread_sim.py's policy-B 128-control tables (87 counting+drying cells, probit s = 0.1/0.15/0.2/0.3, 12 conditions, 500 runs, seeds [5, i, j], controls
drawn as misread_sim.run_cell does) recomputed with droplet_range_common: B (c* -2.21, must reproduce the log first) and AM2 (plateau h in both families,
c* -2.40, droplet_range_a.log:122-149) on the same simulated runs, each without the control cap (what the log and prereg Sections 3 and 10 print) and with it
(what Section 2 says).  Stage 3: Bak et al. (2023) re-fit under B (must reproduce bak_policy_b.log) and AM2 (POST HOC).
Predictions: droplet_ctrl_am2_predictions.txt, sha256 PRED_SHA, written 2026-10-21 01:53:18 PDT before any result of this script.
Disclosed before the results run (no verdict or share printed): scratchpad ctrl_am2/time1.py timed protocol2 (0.73 ms/run exact, 0.13 fast) and Fitter.fit
(0.13 s for 300 runs), and copied the earlier session's null-curve cache (null.npz, the 87 cells; results do not depend on it) to scratchpad ctrl_am2/cache.

Run: cd /Users/saptarshighosh/Documents/chirality-project && DR_CACHE=<scratchpad>/ctrl_am2/cache UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_ctrl_am2.py > droplet_ctrl_am2.log 2> droplet_ctrl_am2.err
"""
import hashlib
import re
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
from scipy.special import ndtr
from scipy.stats import binom

import droplet_range_common as C
from droplet_range_common import NAMES, RULES

PRED_SHA = "e60efb84dfc449661160e4f6becdba3dd277cc7753c9cbc77e6c3c5fd22f1a32"
PROCS, REPS, NCTRL = 4, 500, 128
CONDS = [(cv, e) for cv in (0.05, 0.15) for e in (0.0, 0.02, 0.04, 0.06, 0.08, 0.10)]       # misread_sim.CONDS
SVALS = (0.1, 0.15, 0.2, 0.3)
CSTAR = {"B": -2.21, "AM2": -2.40, "Bn": -2.21, "AM2n": -2.40}                                # n = no control cap
RULES["Bn"], RULES["AM2n"] = ("meas", ("c0", "s0"), False), ("meas", ("c0", "cH", "s0", "sH"), False)
_CONTROLS = C.controls
C.controls = lambda rng, reps, eps, mode: ((rng.random((reps, 192)) < eps)[:, :NCTRL].sum(1) if mode == "xor" else _CONTROLS(rng, reps, eps, mode))   # misread_sim.run_cell's stream (workers import this module too)
NN = np.arange(128, 1025)
HERE = Path(__file__).parent
T0 = time.time()
RES = []
V, Z = NAMES.index("STEP"), NAMES.index("NS")


def say(s=""):
    print(s)


def held(name, ok, detail):
    RES.append((name, bool(ok)))
    say(f"  {name} {'HELD' if ok else 'FAILED'}: {detail}")


# ---------------------------------------------------------------- Part 1
def tails(p, cap=(1, 10)):
    """(NS, FS) over NN: NOT SCORABLE iff k > n * cap, i.e. k > floor(n * cap) in integers; FS = 1 - NS."""
    km = NN * cap[0] // cap[1]
    return binom.sf(km, NN, p), binom.cdf(km, NN, p)


def first_stable(ok, sel=None):
    """(first n where ok, smallest n after which ok holds for every larger n up to 1024), over the n in sel (default all)."""
    n, o = (NN, ok) if sel is None else (NN[sel], ok[sel])
    bad = np.flatnonzero(~o)
    first = int(n[o][0]) if o.any() else None
    stable = int(n[0]) if not bad.size else int(n[bad[-1] + 1]) if bad[-1] + 1 < len(n) else None
    return first, stable


def part1():
    say("PART 1.  Exact binomial arithmetic (scipy.stats.binom), n controls = 128..1024, true misread p the same for every control.")
    say("NS = P(k > n/10), FS (false SCORABLE) = P(k <= n/10).  No cost per control is stated in the prereg (Sections 6-9): counts only.")
    NS = {p: tails(p)[0] for p in (0.06, 0.08, 0.10)}
    FS = {p: tails(p)[1] for p in (0.11, 0.12, 0.13, 0.14)}
    say("\n  n     k_max  " + "".join(f"NS({round(p * 100)}%)  " for p in NS) + "".join(f"FS({round(p * 100)}%)  " for p in FS))
    for n in (128, 129, 130, 139, 140, 160, 192, 240, 256, 320, 400, 480, 500, 509, 510, 512, 560, 640, 768, 800, 1000, 1019, 1020, 1024):
        i = n - 128
        say(f"  {n:<6}{n // 10:<7}" + "".join(f"{NS[p][i]:.3f}    " for p in NS) + "".join(f"{FS[p][i]:.3f}    " for p in FS))
    say("\n  Smallest n meeting the target (first); smallest n after which it holds for every larger n up to 1024 (stable); the same over n % 10 == 0 and")
    say("  n % 80 == 0 (16 controls per level x 8 levels with an integral 10% threshold); min / max of the probability over all n = 128..1024.")
    out = {}
    for lab, d, tgt in [(f"NS({round(p * 100)}%) <= 0.05", NS[p], 0.05) for p in NS] + [(f"FS({round(p * 100)}%) <= 0.10", FS[p], 0.10) for p in FS]:
        ok = d <= tgt
        out[lab] = first_stable(ok)
        a, b = first_stable(ok, NN % 10 == 0), first_stable(ok, NN % 80 == 0)
        say(f"  {lab:<16} first {str(out[lab][0]):>5} stable {str(out[lab][1]):>5} | n%10==0: first {str(a[0]):>5} stable {str(a[1]):>5}"
            f" | n%80==0: first {str(b[0]):>5} stable {str(b[1]):>5} | min {d.min():.3f} at n={int(NN[d.argmin()])}, max {d.max():.3f} at n={int(NN[d.argmax()])}")
    say("\n  Lower caps (NOT SCORABLE iff k > cap * n), pair NS(6%) <= 0.05 and FS(11%) <= 0.10; stable = holds for every larger n up to 1024.")
    pair = {}
    for cp in (6, 7, 8, 9, 10):
        n6, f11 = tails(0.06, (cp, 100))[0] <= 0.05, tails(0.11, (cp, 100))[1] <= 0.10
        a, b, c = first_stable(n6), first_stable(f11), first_stable(n6 & f11)
        pair[cp] = c[1]
        at = ""
        if c[1] is not None:
            i = c[1] - 128
            at = "  | at the stable n: " + ", ".join(f"NS({q}%) {tails(q / 100, (cp, 100))[0][i]:.3f}" for q in (6, 8, 10)) \
                + ", " + ", ".join(f"FS({q}%) {tails(q / 100, (cp, 100))[1][i]:.3f}" for q in (11, 12))
        say(f"  cap {cp:>2}%: NS(6%)<=.05 first {a[0]} stable {a[1]} | FS(11%)<=.10 first {b[0]} stable {b[1]} | pair first {c[0]} stable {c[1]}{at}")
    return NS, FS, out, pair


def score_part1(NS, FS, out, pair):
    say("\nPart 1 predictions")
    f6, s6 = out["NS(6%) <= 0.05"]
    held("P1", f6 == 128 and s6 is not None and 128 <= s6 <= 150 and 0.044 <= NS[0.06].max() <= 0.056,
         f"6%: first {f6}, stable {s6}, max over n {NS[0.06].max():.3f} at n={int(NN[NS[0.06].argmax()])}")
    f8, s8 = out["NS(8%) <= 0.05"]
    saw = f8 is not None and bool((NS[0.08][NN > f8] > 0.05).any())
    held("P2", f8 is not None and 400 <= f8 <= 520 and s8 is not None and 480 <= s8 <= 640 and saw, f"8%: first {f8}, stable {s8}, some n above first has NS > 0.05: {saw}")
    held("P3", not (NS[0.10] <= 0.05).any() and 0.42 <= NS[0.10].min() <= 0.50 and 0.52 <= NS[0.10].max() <= 0.60,
         f"10%: min {NS[0.10].min():.3f} at n={int(NN[NS[0.10].argmin()])}, max {NS[0.10].max():.3f}, any n with NS <= .05: {(NS[0.10] <= 0.05).any()}")
    held("P4", not (FS[0.11] <= 0.10).any() and 0.30 <= FS[0.11][0] <= 0.37 and 0.12 <= FS[0.11].min() <= 0.19,
         f"11%: FS(128) {FS[0.11][0]:.3f}, min {FS[0.11].min():.3f} at n={int(NN[FS[0.11].argmin()])}, any n with FS <= .10: {(FS[0.11] <= 0.10).any()}")
    f12, f13, f14 = (out[f"FS({q}%) <= 0.10"][0] for q in (12, 13, 14))
    held("P5", f12 is not None and 420 <= f12 <= 600 and f13 is not None and 200 <= f13 <= 330 and f14 is not None and f14 <= 170 and 0.07 <= FS[0.14][0] <= 0.13,
         f"first n with FS <= .10: 12% {f12}, 13% {f13}, 14% {f14}; FS(14%, 128) {FS[0.14][0]:.3f}")
    ok89 = [pair[c] for c in (8, 9) if pair[c] is not None]
    held("P6", not (FS[0.11] <= 0.10).any() and bool(ok89) and 330 <= min(ok89) <= 520 and pair[6] is None and pair[7] is None and pair[10] is None,
         "pair stable n by cap: " + ", ".join(f"{c}%: {pair[c]}" for c in pair))
    for c in (8, 9):
        if pair[c] is not None:
            say(f"  P6 price clause (no numeric band in the file, reported not scored): cap {c}%, n = {pair[c]}: NS at a true 8% = {tails(0.08, (c, 100))[0][pair[c] - 128]:.3f}, at a true 10% = {tails(0.10, (c, 100))[0][pair[c] - 128]:.3f}")


# ---------------------------------------------------------------- Part 2
def read_log():
    """misread_sim.log, B128 column of the six tables: {key: (shares[12], first data line, printed digits)}."""
    L = (HERE / "misread_sim.log").read_text().splitlines()
    heads = {"cnt": ("Counting cells: highest STEP share", 3), 0.1: ("probit s=0.1: STEP", 2), 0.15: ("probit s=0.15: STEP", 2), 0.2: ("probit s=0.2: STEP", 2),
             0.3: ("probit s=0.3: STEP", 2), "ns": ("probit s=0.15: NOT SCORABLE", 2)}
    out = {}
    for key, (title, digits) in heads.items():
        k = next(i for i, l in enumerate(L) if l.startswith(title))
        rows = [l.split() for l in L[k + 2:k + 14]]
        assert L[k + 1].split()[:2] == ["CV", "misread"] and [(float(r[0]), float(r[1])) for r in rows] == CONDS
        out[key] = (np.array([float(r[6]) for r in rows]), k + 3, digits)       # r[2:5] = A 64/128/192, r[5:8] = B 64/128/192
    return out


def stage2():
    nulls = C.null_cells(PROCS)
    assert len(nulls) == 87, len(nulls)
    steps = [(f"probit s={s}", C.probit(0, s), 0 * C.LF) for s in SVALS]                  # droplet_sim2.truths: ndtr(-ln(V/V50)/s), p0 = 0
    rules = ("B", "AM2", "Bn", "AM2n")
    tasks = [dict(cells=[(n, P, p0, cv, e, "xor", 0.0, [5, i, j]) for j, (cv, e) in enumerate(CONDS)], reps=REPS, rules=rules, design="base", V0=None,
                  cstar=CSTAR, want_v0=False, exact=True) for i, (n, P, p0) in enumerate(nulls + steps)]
    res = C.pmap(tasks, PROCS, "ctrl")
    return {r: np.array([[res[i][j]["counts"][r] for j in range(len(CONDS))] for i in range(len(tasks))]) for r in rules}   # (91, 12, 5) run counts


def table(title, cols, rows, fmt):
    say(f"\n{title}")
    say("  CV   misread " + " ".join(f"{c:>9}" for c in cols))
    for k, (cv, e) in enumerate(CONDS):
        say(f"  {cv:.2f} {e:.2f}    " + " ".join(fmt.format(r[k]).rjust(9) for r in rows))


def band(share, svals, lo, hi):
    """(min, max) STEP share over the s in svals and both CVs at true misread in [lo, hi]."""
    v = [share[87 + SVALS.index(s)][j] for s in svals for j, c in enumerate(CONDS) if lo - 1e-9 <= c[1] <= hi + 1e-9]
    return min(v), max(v)


def part2(cnt):
    log = read_log()
    stp, nsc = {r: cnt[r][:, :, V] for r in cnt}, {r: cnt[r][:, :, Z] for r in cnt}      # run counts, (91, 12)
    cm = {r: stp[r][:87].max(0) for r in cnt}                                             # counting maximum per condition (runs)
    say("\nPART 2.  Model study, not data.  500 runs per cell, seeds [5, i, j], 128 controls.  'no cap' (Bn, AM2n) = misread_sim.log's policy B and the prereg Sections 3 and 10;")
    say("'cap' (B, AM2) = prereg Section 2 item 1 as written (NOT SCORABLE when more than 10% of the 128 controls are read wrong).")
    say("'Bn' = this run with c* = -2.21 (misread_sim used its unrounded calibrated value, so 0-2 runs per cell can flip).")
    sh = lambda r, k: stp[r][k] / REPS
    ok_digits, ok_tol = [], []
    jobs = [("cnt", "Counting cells: highest STEP share over the 87 cells", {r: cm[r] / REPS for r in cnt})] \
        + [(s, f"probit s={s}: STEP share", {r: sh(r, 87 + k) for r in cnt}) for k, s in enumerate(SVALS)]
    for key, ttl, mine in jobs:
        lv, line, dg = log[key]
        fm = "{:5.3f}" if dg == 3 else "{:5.2f}"
        table(f"{ttl}  (log column: misread_sim.log:{line}-{line + 11}, B 128)", ["log B128", "Bn", "Bn-log", "AM2n", "AM2n-Bn", "B cap", "AM2 cap", "AM2c-Bc"],
              [lv, mine["Bn"], mine["Bn"] - lv, mine["AM2n"], mine["AM2n"] - mine["Bn"], mine["B"], mine["AM2"], mine["AM2"] - mine["B"]], fm)
        ok_digits.append(np.abs(np.round(mine["Bn"], dg) - lv) < 1e-9)
        ok_tol.append(np.abs(mine["Bn"] - lv) <= (0.004 if dg == 3 else 0.015) + 1e-9)
    ok_digits, ok_tol = np.concatenate(ok_digits), np.concatenate(ok_tol)
    lvn, lnl, _ = log["ns"]
    nsn = nsc["Bn"][88] / REPS
    binp = np.array([binom.sf(12, NCTRL, e) for _, e in CONDS])
    table(f"probit s=0.15: NOT SCORABLE share  (log column misread_sim.log:{lnl}-{lnl + 11}, B 128, no cap)", ["log B128", "Bn", "AM2n", "B cap", "AM2 cap", "binomial"],
          [lvn, nsn, nsc["AM2n"][88] / REPS, nsc["B"][88] / REPS, nsc["AM2"][88] / REPS, binp], "{:5.3f}")
    summary = [("s<=0.15, misread <= 0.06", (0.1, 0.15), 0, 0.06), ("s<=0.15, misread 0.08-0.10", (0.1, 0.15), 0.08, 0.10),
               ("s=0.2,  misread <= 0.06", (0.2,), 0, 0.06), ("s=0.2,  misread 0.08-0.10", (0.2,), 0.08, 0.10),
               ("s=0.3,  misread <= 0.06", (0.3,), 0, 0.06), ("s=0.3,  misread 0.08-0.10", (0.3,), 0.08, 0.10)]
    say("\nSummary in the form of prereg Sections 3 and 10 (percent of runs called STEP; prereg quotes: counting max 5.4%, s<=0.15 91-100% and 81-97%, s=0.2 84-100% and 68-91%)")
    say(f"  {'':27}" + "".join(f"{n:>13}" for n in ("log B128", "Bn", "AM2n", "B cap", "AM2 cap")))
    lg = {s: log[s][0] for s in SVALS}
    say(f"  {'counting max, all 12':27}" + "".join(f"{x:>13}" for x in (f"{100 * log['cnt'][0].max():.1f}", *(f"{100 * cm[r].max() / REPS:.1f}" for r in ("Bn", "AM2n", "B", "AM2")))))
    for lab, sv, lo, hi in summary:
        lgb = [lg[s][j] for s in sv for j, c in enumerate(CONDS) if lo - 1e-9 <= c[1] <= hi + 1e-9]
        bd = [band(stp[r] / REPS, sv, lo, hi) for r in ("Bn", "AM2n", "B", "AM2")]
        cells = [f"{100 * min(lgb):.0f}-{100 * max(lgb):.0f}"] + [f"{100 * a:.0f}-{100 * b:.0f}" for a, b in bd]
        say(f"  {lab:27}" + "".join(f"{c:>13}" for c in cells))

    say("\nPredictions, part 2")
    held("P7", ok_digits.sum() >= 58 and ok_tol.all() and abs(cm["Bn"].max() / REPS - 0.054) <= 0.004 + 1e-9 and (np.round(nsn, 2) == 0).all(),
         f"{ok_digits.sum()}/60 cells equal the printed digits, all within tolerance: {ok_tol.all()}, counting max {cm['Bn'].max() / REPS:.3f} (log 0.054), Bn NS share all 0.00: {(np.round(nsn, 2) == 0).all()}")
    gate = dict(RES)["P7"]
    d = stp["AM2n"][87:] - stp["Bn"][87:]                                                   # (4, 12) runs; 25 runs = 0.05
    big = [(SVALS[k], CONDS[j], int(d[k, j]) / REPS) for k in range(4) for j in range(12) if abs(d[k, j]) > 25]
    held("P8", (np.abs(d) <= 25).sum() >= 45 and d.min() >= -30 and d.max() <= 10 and np.abs(d[:2]).max() <= 15 and all(s in (0.2, 0.3) and c[1] >= 0.06 for s, c, _ in big),
         f"{(np.abs(d) <= 25).sum()}/48 cells within 0.05; AM2n-Bn range [{d.min() / REPS:+.3f}, {d.max() / REPS:+.3f}]; s<=0.15 max |d| {np.abs(d[:2]).max() / REPS:.3f}; cells beyond 0.05: {big}")
    dc = cm["AM2n"] - cm["Bn"]
    held("P9", (cm["AM2n"] <= 35).all() and (np.abs(dc) <= 15).all(), f"AM2n counting max {cm['AM2n'].max() / REPS:.3f}, max |AM2n-Bn| per condition {np.abs(dc).max() / REPS:.3f}")
    cen = {0.06: (0.044, 0.03), 0.08: (0.224, 0.05), 0.10: (0.519, 0.06)}
    nsB, nsA = nsc["B"][88], nsc["AM2"][88]
    ok10 = all(abs(nsB[j] / REPS - cen[e][0]) <= cen[e][1] + 1e-9 for j, (_, e) in enumerate(CONDS) if e in cen) \
        and all(nsB[j] <= 5 for j, (_, e) in enumerate(CONDS) if e <= 0.04) and np.abs(nsB - nsA).max() <= 10
    held("P10", ok10, "cap NS share at s=0.15: " + ", ".join(f"CV{c[0]}/{c[1]}: B {nsB[j] / REPS:.3f} AM2 {nsA[j] / REPS:.3f}" for j, c in enumerate(CONDS) if c[1] >= 0.06)
         + f"; at misread <= 0.04 max {max(nsB[j] for j, c in enumerate(CONDS) if c[1] <= 0.04) / REPS:.3f}")
    viol = [(r, SVALS[k], CONDS[j]) for r in ("B", "AM2") for k in range(4) for j in range(12)
            if stp[r][87 + k][j] / REPS > stp[r + "n"][87 + k][j] / REPS * (1 - nsc[r][87 + k][j] / REPS) + 0.03 + 1e-9]
    ex = stp["B"][88][CONDS.index((0.05, 0.10))] / REPS
    held("P11", not viol and 0.30 <= ex <= 0.50, f"violations {len(viol)} {viol[:6]}; B cap, s=0.15, CV 0.05, misread 0.10: {ex:.2f} (no cap {stp['Bn'][88][CONDS.index((0.05, 0.10))] / REPS:.2f})")
    held("P12", all((cm[r] <= cm[r + "n"] + 5).all() for r in ("B", "AM2")),
         f"max (cap - no cap) counting maximum: B {np.max(cm['B'] - cm['Bn']) / REPS:+.3f}, AM2 {np.max(cm['AM2'] - cm['AM2n']) / REPS:+.3f}")

    say("\nDecision rules")
    exceed = [(f"s={SVALS[k]}", CONDS[j], int(d[k, j]) / REPS) for k in range(4) for j in range(12) if abs(d[k, j]) > 25]
    cbad = [(CONDS[j], int(cm["AM2n"][j]) / REPS) for j in range(12) if cm["AM2n"][j] > 35]
    if not gate:
        say("  R1: the gate (P7) FAILED -> the AM2 comparison is NOT TRUSTED; no prereg wording is proposed from it.")
    elif not exceed and not cbad:
        say("  R2: gate held, all 48 step cells within 0.05 and the AM2 counting maximum <= 0.07 in all 12 conditions -> AM2 can replace B without changing the Section 3 / 10 percentages by more than 0.05.")
    else:
        say(f"  R3: gate held, but {len(exceed)} step cells differ by more than 0.05: {exceed}; AM2 counting maximum above 0.07: {cbad} -> restate the Section 3 / 10 ranges from the AM2 columns above.")
    return gate


# ---------------------------------------------------------------- Stage 3
def stage3():
    from blind_units import BAK_MAIN, BAK_TINC90
    from droplet_sim2 import classify, curves
    from exclusion_bak import bak
    ref = [re.findall(r"e (\d\.\d\d): L ([+-]\d+\.\d\d) sharp\s+([\d.]+) (\w+)", l) for l in (HERE / "bak_policy_b.log").read_text().splitlines() if l.startswith(("Tinc", "main"))]
    say("\nSTAGE 3 (POST HOC).  Bak et al. (2023) Fig. 2a re-fit; no control drops, so the misread rate is tried at 0-10% and the unprinted drop counts at P = 1 at n = 3, 6, 10.")
    say("B = c* -2.21 (must reproduce bak_policy_b.log); AM2 = c* -2.40.  sharp = 2 (best ll - best ll among STEP curves with s <= 0.15), > 3.84 rejects a sharp step;")
    say("for AM2 the STEP curves include the plateau family, so a sharp step with a plateau is not rejected.")
    match, tot, row_i = 0, 0, 0
    for name, series in (("Tinc 90 C", BAK_TINC90), ("main", BAK_MAIN)):
        for n1 in (3, 6, 10):
            x, k, n = bak(series, n1)
            K, N = k[None], n[None]
            F = C.Fitter(x, True)
            ps0 = curves(x)["step"][3]
            psH = np.tile(np.repeat(C.SG, len(C.CGC)), len(C.HG))
            assert len(ps0) == F.P["s0"].shape[0] and len(psH) == F.P["sH"].shape[0]
            z, ns = C.trend_z(K, N, x), ~C.scorable_v(K, N)
            X = np.hstack([K, N])
            rows = {"B": [], "AM2": []}
            for e in (0.0, 0.02, 0.04, 0.06, 0.08, 0.10):
                f = F.fit(K, N, np.array([e]), ("c0", "cH", "s0", "sH"))
                m = {g: (X.astype(F.W(e, g).dtype) @ F.W(e, g))[0] for g in ("s0", "sH")}
                for r, grp, cs in (("B", ("c0", "s0"), -2.21), ("AM2", ("c0", "cH", "s0", "sH"), -2.40)):
                    st = C.rule_stats(f, grp)
                    top = max(f[g][0][0] for g in grp)
                    sharp = 2 * (top - max(m[g][p <= 0.15].max() for g, p in (("s0", ps0), ("sH", psH)) if g in grp))
                    lab = classify(dict(L=st["L"], rng=st["rng"], z=z, ns=ns), cs)[0]
                    rows[r].append((e, st["L"][0], sharp, lab))
            for r in ("B", "AM2"):
                say(f"{name:9} n={n1:<2} {r:3} | " + " | ".join(f"e {e:.2f}: L {L:+5.2f} sharp {sh:5.1f} {lab}" for e, L, sh, lab in rows[r]))
            for (e, L, sh, lab), (re_, rL, rs, rlab) in zip(rows["B"], ref[row_i]):
                tot += 1
                match += (f"{e:.2f}" == re_ and f"{L:+.2f}" == rL and f"{sh:.1f}" == rs and lab == rlab)
            row_i += 1
    say(f"\nB reproduces bak_policy_b.log: {match} of {tot} (L, sharp, label) entries equal the printed ones")
    return match, tot


# ---------------------------------------------------------------- self-check and main
def self_check():
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    C.self_check()
    tail = lambda n, e, kmin: sum(comb(n, k) * e ** k * (1 - e) ** (n - k) for k in range(kmin, n + 1))
    for p in (0.06, 0.1):
        assert abs(tails(p)[0][0] - tail(128, p, 13)) < 1e-12 and abs(tails(p)[0][0] + tails(p)[1][0] - 1) < 1e-12
    assert tails(0.06)[0][2] == binom.sf(13, 130, 0.06) and abs(tails(0.06, (8, 100))[0][0] - binom.sf(10, 128, 0.06)) < 1e-15      # n = 130: k > 13; cap 8%: k > 10.24
    assert first_stable(NN >= 500) == (500, 500) and first_stable((NN >= 500) | (NN == 128)) == (128, 500) and first_stable(NN < 500) == (128, None)
    a, b = np.random.default_rng(0), np.random.default_rng(0)
    assert (C.controls(a, 50, 0.05, "xor") == (b.random((50, 192)) < 0.05)[:, :128].sum(1)).all()                             # the monkeypatch is active
    assert np.allclose(C.probit(0, 0.15), ndtr(-(C.LF - np.log(C.V50)) / 0.15))                                               # droplet_sim2.truths' step
    P = C.probit(0, 0.15)
    o = C.run_batch(dict(cells=[("t", P, 0 * P, 0.05, 0.0, "xor", 0.0, [9, 9, 9])], reps=30, rules=("B", "AM2", "Bn", "AM2n"), design="base", V0=None,
                         cstar=CSTAR, want_v0=False, exact=True))[0]["counts"]
    assert all(v.sum() == 30 for v in o.values()) and (o["B"] == o["Bn"]).all() and (o["AM2"] == o["AM2n"]).all()               # eps 0: the cap never bites


def main():
    assert hashlib.sha256((HERE / "droplet_ctrl_am2_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    self_check()
    say(f"Unit droplet_ctrl_am2: control-drop count (exact binomial) and rule AM2 against rule B.  Simulated drops, not data.  Predictions sha256 {PRED_SHA}.")
    NS, FS, out, pair = part1()
    score_part1(NS, FS, out, pair)
    print(f"[part 1 done {time.time() - T0:.0f}s]", file=sys.stderr, flush=True)
    cnt = stage2()
    print(f"[stage 2 simulated {time.time() - T0:.0f}s]", file=sys.stderr, flush=True)
    part2(cnt)
    stage3()
    say(f"\nHELD {sum(v for _, v in RES)} of {len(RES)}: " + ", ".join(f"{n} {'HELD' if v else 'FAILED'}" for n, v in RES))
    print(f"[done {time.time() - T0:.0f}s]", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
