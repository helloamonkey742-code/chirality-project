"""seedfix.py: remove the Enceladus seed collision and regenerate every Enceladus number that depends on it.

THE BUG. sweep2.run (sweep2.py:165), sweep3.run (sweep3.py:156), multiwall.run4 (multiwall.py:435), rho2side.draws
(rho2side.py:96) and rerun_base's EX (rerun_base.py:170) draw each world's extra inputs (eta, r = p/a, rho) from
np.random.default_rng([seed, i]), i = world index. numpy's SeedSequence ignores trailing zeros, so for world 0
(Enceladus) default_rng([seed, 0]) is default_rng(seed): the very generator sweep2.chains uses for Enceladus's own
chain draws. eta, r, rho therefore repeat the quantiles of bias g, rate k2 and concentration draw for draw
(rerun_base.log:21: correlations 1.000000). Worlds 1-14 use [seed, i] with i != 0, which are already independent.

THE FIX (wrapper only; no project script is edited). np.random.default_rng is wrapped so that a call with the exact
two-element seed [s, 0] returns default_rng([s, 0, 1]); every other call (int seeds, [s, i] with i != 0, and the
three-element streams such as [0, 0, 3] in rerun_base) passes through untouched. Only world 0 moves to a new stream.
Why not re-seed every world (e.g. SeedSequence.spawn): worlds 1-14 are already on independent streams, so keeping them
gives an exact regression: every non-Enceladus number must reproduce the logs to the printed digit, which checks that
the wrapper changed nothing else. [s, 0, 1] is the stream rho2side.posthoc used, so its post-hoc Enceladus rows are a
direct prediction. The wrapper records every call site it remaps (printed at the end of seedfix.log).

WHAT IS REGENERATED (seed 0 tables as in the logs; seeds 0-3 for the spread):
  sweep2.main + scale table + seeds 1-3 (vs sweep2.log), sweep3.main + seeds 1-3 (vs sweep3.log, with sweep3's own
  step e/f regression targets for Enceladus taken from the regenerated sweep2), multiwall.section4 (vs multiwall.log
  Section 4, its sweep3 regression read from the regenerated sweep3 text), rho2side.main (vs rho2side.log, its
  multiwall regression read from the regenerated Section 4), and rerun_base sections 0-11 (vs rerun_base.log; the
  log file and journals/new_analyses_rerun.md are NOT rewritten: rerun_base.FH is redirected to memory and neither
  rerun_base.main nor write_md is called). rerun_base section 0 still compares against the real multiwall.log, so its
  Enceladus row and its P1 line are expected to show the old reference beside the new value.
  Every changed line is classified: Enceladus (the line, segment or block names only Enceladus), aggregate (no world
  named: max-over-worlds and check lines), or VIOLATION (a field naming another world changed, or the line count).

PREDICTIONS (written and committed before any result-producing run; never edited). "Old" = the logged value, seed 0,
no lid unless stated.
  P1  At seeds 0-3 the regenerated Enceladus extras are uncorrelated with its base draws: |corr| < 0.05 for
      (log g, log eta), (log k2, log r), (log c, log rho). (Collided streams give 1.000000.)
  P2  Exact regression: in sweep2.log, sweep3.log, multiwall.log Section 4, rho2side.log and rerun_base.log the
      regenerated text has the same number of lines and no VIOLATION: every field naming a non-Enceladus world is
      identical to the printed digit.
  P3  Enceladus with a lid stays 0.0 in every regenerated final column (seed 0 tables and seeds 0-3).
  P4  multiwall one-sided main no lid 28.3 -> 26.4 +- 0.3 (agents' estimate about 26.4, rho2side.log post-hoc).
  P5  sweep2 step f no lid 3.9 -> 6.7 +- 0.3 (agents' estimate about 6.7).
  P6  rho2side cascade x10 gap |two-sided - mixture| for Enceladus 1.4 -> <= 0.15 (agents' estimate about 0.05),
      and the all-world maximum rho2side prints (1.45 points) falls to <= 0.2.
  P7  The other rho2side post-hoc values reproduce within 0.15: step e 24.5 -> 24.7, two-sided main 26.5 -> 25.6,
      one-sided x10 19.9 -> 19.4, two-sided x10 20.8 -> 22.1.
  P8  sweep3 main (reversed walls counted) 6.6 -> between 7.5 and 10.5.
  P9  sweep3 bound (iii) 23.9 -> between 23.0 and 25.5.
  P10 multiwall sum rule 20.1 -> between 19.0 and 21.0.
  P11 multiwall main P_chance 8.8 -> between 9.3 and 10.3; sweep2 step-f split (g) P_chance 2.0 -> <= 1.8.
  P12 sweep2 ladder steps b, c, d, e (no lid) each move by <= 1.5 points (23.3, 27.2, 26.4, 24.5).
  P13 rerun_base Section 1, Enceladus main: the three duplicated pairs (bias g / eta, rate k2 / feed share r,
      concentration / D_z ratio rho) now differ by > 1 point each, and vertical D_z is still the top effect.
  P14 rerun_base Section 6 Kang / strat-alt Enceladus rows: each new main final within 0.6 points of
      old final x (new main no lid / old main no lid).
  P15 rerun_base Section 10 Enceladus main f* values (M=1e3 1.5e-14, field 4.0e-15) stay within a factor 2.
  P16 Seeds 1-3: every regenerated Enceladus no-lid value in the per-seed table is within 0.5 points of seed 0.

Run: OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python seedfix.py > seedfix.log

RESULTS (seedfix.log, 2026-10-12): N = 200000, 137 s. 14 held, 2 failed. (A smoke run at SEEDFIX_N=2000 was made
after the predictions commit, only to catch crashes; its numbers are not used.)
  P1  HELD.   max |corr| 0.0038 over seeds 0-3 and the three pairs (collided: 1.000000 every seed).
  P2  FAILED as registered: 1 VIOLATION, rerun_base.log:305. It is the check line "P14 FAILED: Earth L 46.6 ..."; the
      only change in its Earth-named segment is the verdict word (FAILED -> HELD), flipped by the Enceladus f* values
      later on the line; every number naming another world on that line is identical. Everywhere else 0 violations:
      sweep2.log 85/85, sweep3.log 75/75, multiwall.log Section 4 61/61 and rho2side.log 81/81 lines match except
      Enceladus rows and Enceladus-including aggregates; rerun_base.log 328/328 lines, line counts equal.
  P3  HELD.   all 63 regenerated Enceladus final values are 0.0.
  P4  HELD.   one-sided main 28.3 -> 26.4 (26.36).
  P5  HELD.   step f 3.9 -> 6.7 (6.68).
  P6  HELD.   Enceladus x10 gap 1.45 -> 0.00; all-world max 1.45 -> 0.07.
  P7  HELD.   step e 24.75, two-sided main 25.59, one-sided x10 19.44, two-sided x10 22.10 (post-hoc 24.7/25.6/19.4/22.1).
  P8  HELD.   sweep3 main 6.6 -> 9.0.
  P9  HELD.   bound 23.9 -> 23.6.
  P10 HELD.   sum rule 20.1 -> 19.9.
  P11 HELD.   multiwall P_chance 8.8 -> 9.8; step-f P_chance 2.0 -> 1.5.
  P12 HELD.   ladder 21.8 24.1 23.8 27.8 27.1 24.7 6.7 (b +0.5, c +0.6, d +0.7, e +0.2).
  P13 FAILED: bias g +5.9 vs eta +6.3 differ by only 0.4 (k2 +6.4 vs r +3.0 and concentration +29.7 vs rho +6.0 do
      separate; D_z stays top at +44.0). Independent g and eta now have about equal effects.
  P14 HELD.   max deviation 0.07; finals 7.9->7.4, 5.5->5.2, 22.9->21.3, 19.0->17.7, 2.2->2.0, 9.5->8.8, 9.9->9.3.
  P15 HELD.   M=1e3 1.5e-14 -> 9.0e-15, field 4.0e-15 -> 2.6e-15 (factor 0.6-0.65).
  P16 HELD.   max |seed k - seed 0| 0.28 over the ten Enceladus no-lid columns.
"""
import contextlib
import io
import os
import re
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
T0 = time.time()

# ---------------------------------------------------------------- the fix
ORIG_RNG = np.random.default_rng
REMAPPED = {}


def fixed_rng(seed=None, *a, **kw):
    """default_rng with the world-0 collision removed: the exact seed [s, 0] becomes [s, 0, 1]."""
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        f = sys._getframe(1)
        key = f"{os.path.basename(f.f_code.co_filename)}:{f.f_lineno}"
        REMAPPED[key] = REMAPPED.get(key, 0) + 1
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng   # before the project modules import (rerun_base builds EX at import)

import uncertainty as U   # noqa: E402
U.N = int(os.environ.get("SEEDFIX_N", U.N))   # smoke test only: SEEDFIX_N=2000 (rerun_base reads U.N at import)
import sweep2 as S        # noqa: E402
import sweep3 as S3       # noqa: E402
import multiwall as M     # noqa: E402
import rho2side as R2     # noqa: E402
import uncertainty as U   # noqa: E402
import rerun_base as RB   # noqa: E402

N = U.N
NAMES = list(S.chains(0)[1])
assert NAMES.index("Enceladus") == 0, NAMES
TOK = [n for n in NAMES if n != "Early Earth ocean"] + ["Early", "Earth"]
ENC = "Enceladus"
RES = {}       # prediction -> (held, detail)
TMP = tempfile.mkdtemp(prefix="seedfix_")


def check(p, ok, detail):
    RES[p] = (bool(ok), detail)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def cap(fn, *a, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn(*a, **kw)
    return r, buf.getvalue().splitlines()


def lap(lab):
    print(f"[{lab}: {time.time() - T0:.0f} s elapsed]")


# ---------------------------------------------------------------- per-seed Enceladus table
def ex_from(rng, lg=(0.0, 2.0)):
    return {"eta": 10 ** rng.uniform(-1, 1, N), "r": 10 ** rng.uniform(-2, 0, N), "rho": 10 ** rng.uniform(*lg, N)}


COLS = ("step e", "step f", "sweep3 main", "sweep3 bound", "mw main", "mw x10", "mw x5", "mw sum",
        "2-sided main", "2-sided x10")


def enceladus(seed, rngf):
    _, W = S.chains(seed)
    assert list(W).index(ENC) == 0
    x, _, a, salt = W[ENC]
    ex, ex2 = ex_from(rngf([seed, 0])), ex_from(rngf([seed, 0]), R2.TWO)
    cc = [np.corrcoef(np.log10(x[k]), np.log10(ex[e]))[0, 1]
          for k, e in (("bias g", "eta"), ("rate k2", "r"), ("concentration", "rho"))]
    out = {}
    for li, lid in enumerate((np.zeros(N, bool), salt)):
        L, dz = a[1], a[4]
        out[li] = [S3.split(*S.outcome(x, L, dz, "e", ex, lid)), S3.split(*S.outcome(x, L, dz, "f", ex, lid)),
                   S3.split(*S3.outcome3(x, L, dz, ex, lid)[:2]),
                   S3.split(*S3.outcome3(x, L, dz, ex, lid, veto=False)[:2]),
                   S3.split(*M.multi(x, L, dz, ex, lid)[:2]), S3.split(*M.multi(x, L, dz, ex, lid, scale=10.0)[:2]),
                   S3.split(*M.multi(x, L, dz, ex, lid, scale=5.0)[:2]), S3.split(*M.multi(x, L, dz, ex, lid, rule="sum")[:2]),
                   S3.split(*M.multi(x, L, dz, ex2, lid)[:2]), S3.split(*M.multi(x, L, dz, ex2, lid, scale=10.0)[:2])]
    gap = abs(out[0][9][0] - 0.5 * (out[0][0][0] + out[0][5][0]))
    return out, cc, gap


print(f"seedfix.py, N = {N}: Enceladus (world index {NAMES.index(ENC)}) extra draws moved from default_rng([s, 0]) "
      f"(== default_rng(s)) to default_rng([s, 0, 1]); worlds 1-14 unchanged")
print(f"default_rng([0, 0]) == default_rng(0) in numpy: "
      f"{np.array_equal(ORIG_RNG([0, 0]).uniform(size=8), ORIG_RNG(0).uniform(size=8))}; "
      f"wrapped [0, 0] == [0, 0, 1]: {np.array_equal(fixed_rng([0, 0]).uniform(size=8), ORIG_RNG([0, 0, 1]).uniform(size=8))}")
print("\n=== A. Enceladus, seeds 0-3, collided (as logged) -> fixed. P_pved / P_chance / P_none %, no lid; "
      "final P_pved after '|' ===")
PS = {}
for seed in range(4):
    old, cco, gapo = enceladus(seed, ORIG_RNG)
    new, ccn, gapn = enceladus(seed, fixed_rng)
    PS[seed] = (old, new, ccn, gapo, gapn)
    print(f"seed {seed}: corr (log g, log eta) / (log k2, log r) / (log c, log rho): collided "
          + " / ".join(f"{c:.6f}" for c in cco) + "  fixed " + " / ".join(f"{c:+.4f}" for c in ccn))
    for lab, i in zip(COLS, range(10)):
        o, n_ = old[0][i], new[0][i]
        print(f"  {lab:13} {o[0]:5.1f} / {o[1]:4.1f} / {o[2]:5.1f} | {old[1][i][0]:3.1f}   ->   "
              f"{n_[0]:5.1f} / {n_[1]:4.1f} / {n_[2]:5.1f} | {new[1][i][0]:3.1f}   (P_pved {n_[0] - o[0]:+5.2f})")
    print(f"  x10 |two-sided - mixture (step e + one-sided x10)/2|: {gapo:.2f} -> {gapn:.2f}")
print("\nSpread of the fixed Enceladus no-lid P_pved over seeds 0-3 (min-max):")
for i, lab in enumerate(COLS):
    v = [PS[s][1][0][i][0] for s in range(4)]
    print(f"  {lab:13} {min(v):5.1f}-{max(v):5.1f}  (seed 0 {v[0]:.1f}, max |seed k - seed 0| {max(abs(u - v[0]) for u in v):.2f})")
lap("A")

# ---------------------------------------------------------------- line diff against the logs
SEP = re.compile(r" \|\|? | ; |; |, ")
HDR = re.compile(r"^(\S.*?) \((main|x10)\):")
CHK = re.compile(r"^\s*P\d+ (HELD|FAILED)")


def who(s):
    return {t for t in TOK if t in s}


def compare(name, old, start, new):
    """old: whole log lines; start: index of new[0] in old. Returns {'enc': [...], 'agg': [...], 'bad': [...]}."""
    out = {"enc": [], "agg": [], "bad": []}
    seg = old[start:start + len(new)]
    if len(seg) != len(new):
        out["bad"].append((0, f"line count {len(seg)} in log vs {len(new)} regenerated", ""))
        return out
    ctx = None
    for j, (o, n) in enumerate(zip(seg, new)):
        ln = start + j + 1
        if not o.strip() or o.startswith("==="):
            ctx = None
        m = HDR.match(o)
        if m:
            ctx = m.group(1)
        if o == n:
            continue
        wo, wn = who(o), who(n)
        if CHK.match(o) or (not wo and not wn and ctx is None):
            kind = "enc" if (wo | wn) == {ENC} else "agg" if not (wo | wn) - {ENC} else None
        elif not wo and not wn:
            kind = "enc" if ctx == ENC else "bad"
        elif wo | wn == {ENC}:
            kind = "enc"
        else:
            kind = None
        if kind is None:   # several worlds on the line: every changed segment must name only Enceladus
            so, sn = SEP.split(o), SEP.split(n)
            if len(so) != len(sn):
                kind = "bad"
            else:
                ch = [(a, b) for a, b in zip(so, sn) if a != b]
                names = set().union(*(who(a) | who(b) for a, b in ch))
                kind = "enc" if names == {ENC} else "agg" if not names else "bad"
        out[kind].append((ln, o, n))
    return out


DIFFS = {}


def report(name, old, start, new):
    d = compare(name, old, start, new)
    DIFFS[name] = d
    print(f"\n--- {name}: lines {start + 1}-{start + len(new)} compared ({len(new)} regenerated); "
          f"{len(d['enc'])} Enceladus lines changed, {len(d['agg'])} aggregate lines changed, "
          f"{len(d['bad'])} VIOLATIONS ---")
    for kind, lab in (("enc", "Enceladus"), ("agg", "aggregate (no world named; includes Enceladus)"),
                      ("bad", "VIOLATION")):
        for ln, o, n in d[kind]:
            print(f"  [{lab}] {name}:{ln}\n    old: {o}\n    new: {n}")
    return d


def logl(f):
    return open(os.path.join(HERE, f)).read().splitlines()


# ---------------------------------------------------------------- B. sweep2
print("\n=== B. sweep2 (steps 0-f, step-f split g, sensitivities, seeds 1-3) vs sweep2.log ===")
res2, lines2 = cap(S.main, 0)
_, l2b = cap(S.scale_table)
seeds2 = ["", "Seed stability, step f (no lid / final), max |seed k - seed 0| over worlds, points:"]
for s in (1, 2, 3):
    r = S.run(s)
    d = max(abs(r[n]["f"][k] - res2[n]["f"][k]) for n in res2 if n != "_A" for k in (0, 1))
    seeds2.append(f"  seed {s}: max diff {d:.2f}   " + ", ".join(f"{n} {r[n]['f'][0]:.1f}/{r[n]['f'][1]:.1f}"
                                                              for n in r if n != "_A"))
new2 = lines2 + l2b + seeds2
report("sweep2.log", logl("sweep2.log"), 0, new2)
lap("B")

# ---------------------------------------------------------------- C. sweep3
print("\n=== C. sweep3 (main, three-way split, bound iii, seeds 1-3) vs sweep3.log ===")
e2 = res2[ENC]
S3.LOG_E[ENC] = tuple(round(v, 1) for v in e2["e"])
S3.LOG_F[ENC] = tuple(round(v, 1) for v in e2["f"])
print(f"sweep3's own regression targets for Enceladus set to the regenerated sweep2 steps e {S3.LOG_E[ENC]} and "
      f"f {S3.LOG_F[ENC]} (were 24.5/0.0 and 3.9/0.0); the 14 other worlds keep sweep3.py's targets")
res3, lines3 = cap(S3.main)
seeds3 = ["", "Seed stability, sweep3 P_pved / P_chance / P_none (no lid and final), max |seed k - seed 0| over worlds, points:"]
for s in (1, 2, 3):
    r = S3.run(s, full=False)
    d = max(abs(r[n]["main"][l][j] - res3[n]["main"][l][j]) for n in res3 for l in (0, 1) for j in range(3))
    seeds3.append(f"  seed {s}: max diff {d:.2f}   " + ", ".join(f"{n} {r[n]['main'][0][0]:.1f}/{r[n]['main'][1][0]:.1f}"
                                                              for n in r))
new3 = lines3 + seeds3
report("sweep3.log", logl("sweep3.log"), 0, new3)
open(os.path.join(TMP, "sweep3.log"), "w").write("\n".join(new3) + "\n")
lap("C")

# ---------------------------------------------------------------- D. multiwall Section 4
print("\n=== D. multiwall Section 4 (main, cascade x10, sum rule, three-way split, bound, seeds 1-3) vs multiwall.log ===")
M.__file__ = os.path.join(TMP, "multiwall.py")     # sweep3_log() now reads the regenerated sweep3 text
_, new4 = cap(M.section4)
oldm = logl("multiwall.log")
s4 = next(i for i, l in enumerate(oldm) if l.startswith("=== Section 4")) - 1
report("multiwall.log", oldm, s4, new4)
open(os.path.join(TMP, "multiwall.log"), "w").write("\n".join(new4) + "\n")
lap("D")

# ---------------------------------------------------------------- E. rho2side
print("\n=== E. rho2side main (one- and two-sided, x10, mixture, real-profile sink, seeds 0-3) vs rho2side.log ===")
R2.__file__ = os.path.join(TMP, "rho2side.py")     # multiwall_log() now reads the regenerated Section 4
(one, two), new5 = cap(R2.main)
report("rho2side.log", logl("rho2side.log"), 0, new5)
lap("E")

# ---------------------------------------------------------------- F. rerun_base
print("\n=== F. rerun_base sections 0-11 (H19/H20/H27 tables) vs rerun_base.log; log and md NOT rewritten ===")
RB.FH = io.StringIO()


def rb_all():
    RB.say(f"rerun_base.py, seed {RB.SEED}, N = {RB.N}: H19/H20/H27 base-model results through multiwall.multi "
           f"(main and cascade x10)")
    for sec in (RB.section0, RB.section1, RB.section2, RB.section3, RB.section4, RB.section5, RB.section6,
                RB.section7, RB.section8, RB.section9, RB.section10, RB.section11):
        sec()
        sys.stderr.write(f"rerun_base {sec.__name__} done at {time.time() - T0:.0f} s\n")
    RB.say("\nPredictions: " + ", ".join(f"{p} {'held' if ok else 'FAILED'}" for p, ok in RB.CHECKS))


cap(rb_all)
report("rerun_base.log", logl("rerun_base.log"), 0, RB.LINES)
lap("F")

# ---------------------------------------------------------------- predictions
print("\n=== PREDICTIONS ===")
s0o, s0n = PS[0][0][0], PS[0][1][0]
check("P1", all(abs(c) < 0.05 for s in range(4) for c in PS[s][2]),
      "max |corr| seeds 0-3 " + f"{max(abs(c) for s in range(4) for c in PS[s][2]):.4f}")
bad = {k: len(v["bad"]) for k, v in DIFFS.items()}
check("P2", not any(bad.values()), f"violations per log {bad}")
fin = [PS[s][1][1][i][0] for s in range(4) for i in range(10)]
fin += [res2[ENC][s][1] for s in S.STEPS] + [v[1][0] for k in ("main", "noveto", "rise10", "rac") for v in [res3[ENC][k]]]
fin += [two[s][ENC][k][1][0] for s in two for k in ("main", "x10")] + [one[ENC][k][1][0] for k in ("main", "x10")]
fin += [RB.pc(RB.BASE[ENC, lab][2]) for lab in ("main", "x10")]
check("P3", max(fin) == 0.0, f"max Enceladus final P_pved over {len(fin)} regenerated final values {max(fin):.3g}")
check("P4", abs(s0n[4][0] - 26.4) <= 0.3, f"one-sided main {s0o[4][0]:.1f} -> {s0n[4][0]:.2f}")
check("P5", abs(s0n[1][0] - 6.7) <= 0.3, f"step f {s0o[1][0]:.1f} -> {s0n[1][0]:.2f}")
mx = next((l for l in new5 if l.startswith("max |two-sided - mixture|")), "")
allmax = float(re.findall(r"([\d.]+) points", mx)[0]) if mx else np.nan
check("P6", PS[0][4] <= 0.15 and allmax <= 0.2, f"Enceladus x10 gap {PS[0][3]:.2f} -> {PS[0][4]:.2f}; "
      f"all-world max 1.45 -> {allmax:.2f}")
p7 = [(0, 24.7), (8, 25.6), (5, 19.4), (9, 22.1)]
check("P7", all(abs(s0n[i][0] - v) <= 0.15 for i, v in p7),
      ", ".join(f"{COLS[i]} {s0o[i][0]:.1f} -> {s0n[i][0]:.2f} (post-hoc {v})" for i, v in p7))
check("P8", 7.5 <= s0n[2][0] <= 10.5, f"sweep3 main {s0o[2][0]:.1f} -> {s0n[2][0]:.2f}")
check("P9", 23.0 <= s0n[3][0] <= 25.5, f"sweep3 bound {s0o[3][0]:.1f} -> {s0n[3][0]:.2f}")
check("P10", 19.0 <= s0n[7][0] <= 21.0, f"sum rule {s0o[7][0]:.1f} -> {s0n[7][0]:.2f}")
check("P11", 9.3 <= s0n[4][1] <= 10.3 and s0n[1][1] <= 1.8,
      f"multiwall P_chance {s0o[4][1]:.1f} -> {s0n[4][1]:.2f}; step-f P_chance {s0o[1][1]:.1f} -> {s0n[1][1]:.2f}")
oldlad = dict(zip("0abcdef", (21.8, 24.1, 23.3, 27.2, 26.4, 24.5, 3.9)))   # sweep2.log:18
mv = {s: res2[ENC][s][0] - oldlad[s] for s in "bcde"}
check("P12", all(abs(v) <= 1.5 for v in mv.values()), ", ".join(f"step {s} {oldlad[s]} -> {res2[ENC][s][0]:.1f}"
                                                                 for s in "0abcdef"))
eff = RB.SENS[ENC, "main", "_eff"]
pairs = (("bias g", "eta (g factor)"), ("rate k2", "feed share r"), ("concentration", "D_z ratio rho"))
check("P13", all(abs(eff[a] - eff[b]) > 1 for a, b in pairs) and RB.SENS[ENC, "main", "_top"] == "vertical D_z",
      "; ".join(f"{a} {eff[a]:+.1f} vs {b} {eff[b]:+.1f}" for a, b in pairs) + f"; top {RB.SENS[ENC, 'main', '_top']}")
rbo = logl("rerun_base.log")
LID = re.compile(r"^\s+(Kang|strat alt) Enceladus.*main no lid\s+([\d.]+) final\s+([\d.]+)")
lo_ = [m.groups() for m in map(LID.match, rbo) if m]
ln_ = [m.groups() for m in map(LID.match, RB.LINES) if m]
dev = [abs(float(b[2]) - float(a[2]) * float(b[1]) / float(a[1])) for a, b in zip(lo_, ln_)]
check("P14", len(lo_) == len(ln_) == 9 and max(dev) <= 0.6,
      f"{len(ln_)} rows, max |new final - old final x ratio| {max(dev):.2f}; finals "
      + ", ".join(f"{a[2]}->{b[2]}" for a, b in zip(lo_, ln_)))
FS = re.compile(r"^\s+Enceladus\s+base L\s+([\d.]+) \|.*M=1e3 ([\d.e+-]+) \| field ([\d.e+-]+)")
fo = [m.groups() for m in map(FS.match, rbo) if m][0]
fn = [m.groups() for m in map(FS.match, RB.LINES) if m][0]
rat = [float(b) / float(a) for a, b in zip(fo[1:], fn[1:])]
check("P15", all(0.5 <= r <= 2 for r in rat), f"main M=1e3 {fo[1]} -> {fn[1]}, field {fo[2]} -> {fn[2]}")
sp = max(abs(PS[s][1][0][i][0] - PS[0][1][0][i][0]) for s in (1, 2, 3) for i in range(10))
check("P16", sp <= 0.5, f"max |seed k - seed 0| over the {len(COLS)} Enceladus no-lid columns {sp:.2f}")
print("\nSummary: " + ", ".join(f"{p} {'held' if ok else 'FAILED'}" for p, (ok, _) in RES.items()))
print("Remapped call sites ([s, 0] -> [s, 0, 1]), calls: " + ", ".join(f"{k} x{v}" for k, v in sorted(REMAPPED.items())))
print(f"total {time.time() - T0:.0f} s")
