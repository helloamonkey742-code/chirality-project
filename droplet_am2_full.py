"""Rule AM2 at its full-grid threshold against rule B for the NaClO3 droplet test (journals/droplet_prereg.md Sections 2, 3, 10).
Simulated drops, not data.  New file; it imports droplet_ctrl_am2, droplet_range_a and droplet_range_common and edits none of them.

Question: does AM2 (plateau h in both families) still fix the plateau blind spot of rule B, and what does it cost, once its threshold is the
full-grid calibration c*_AM2(128) = -2.9762901306152343 (droplet_am2_gate.log:6-9) instead of -2.40 (droplet_range_a.log:8)?  Also at the subset
value -2.463367462158203 (true misread {0, .06, .10}, 6 of the 12 conditions).  Rules are added to droplet_range_common.RULES in this process only (same fit
groups and cap as AM2 / B, other c*), so one fit per simulated run feeds every threshold.
Worlds: K = droplet_ctrl_am2.stage2 (87 + 4 probit cells x 12 conditions x 500 runs, seeds [5, i, j], controls drawn as misread_sim does, exact pilot);
R = droplet_range_a stages 2-3 (141 cells x CONDS8 x 400 runs, seeds [11, i, j], 20 cells per batch; plateau map 28 x 4 x 600 runs, seeds [12, i, j], 13 per batch).
Gates (journal D0): the K-world text of droplet_ctrl_am2.part2 / stage3 against droplet_ctrl_am2.log, and the R-world B and AM2 lines against droplet_range_a.log.
Predictions: droplet_am2_full_predictions.txt, sha256 PRED_SHA, written 2026-10-21 03:36 PDT before any result of this script.

Run: cd /Users/saptarshighosh/Documents/chirality-project && UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba python droplet_am2_full.py > droplet_am2_full.log 2> droplet_am2_full.err
Self-check only: the same command with `--check` appended and no redirect.
"""
import hashlib
import io
import os
import pickle
import re
import sys
import time
from contextlib import contextmanager, redirect_stdout
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SCR = Path(os.environ.get("AM2F_SCR", "/private/tmp/claude-501/-Users-saptarshighosh-Documents-chirality-project/05c76bba-5069-4c19-a20f-58eb08337c55/scratchpad"))
WORK = SCR / "am2_full"                                   # checkpoints and the truth-curve cache (resume aids; results do not depend on them)
os.environ.setdefault("DR_CACHE", str(WORK / "cache"))

import droplet_ctrl_am2 as P                              # applies the K-world control patch, adds the Bn / AM2n rules (import level: workers repeat it)
import droplet_range_a as A
import droplet_range_common as C
from droplet_range_common import NAMES, RULES
from droplet_sim2 import classify

PRED_SHA = "4a0900aea24cc052301df5d906824b25759df1b3c8b6c4ca66527880431bd9ce"
PROCS = 3
REPS, CONDS, SVALS = P.REPS, P.CONDS, P.SVALS
C_AMF, C_AMS, C_BSUB = -2.9762901306152343, -2.463367462158203, -2.0587861532873992
C_AM2_K, C_AM2_R = -2.40, -2.403298187255859              # droplet_ctrl_am2 rounds, droplet_range_a does not
T0 = time.time()

for _n, _b in (("AMf", "AM2"), ("AMs", "AM2"), ("Bsub", "B")):    # same groups and cap, other c*
    RULES[_n] = RULES[_b]
    RULES[_n + "n"] = (RULES[_b][0], RULES[_b][1], False)
KC = {"B": C.C_B, "AM2": C_AM2_K, "AMs": C_AMS, "AMf": C_AMF, "Bsub": C_BSUB}
KC.update({r + "n": c for r, c in list(KC.items())})
KRULES = tuple(KC)
RC = {"B": C.C_B, "AM2": C_AM2_R, "AMs": C_AMS, "AMf": C_AMF, "Bsub": C_BSUB}
RRULES = tuple(RC)
FAMS_CAP, FAMS_NC = (["Bsub", "B"], ["AM2", "AMs", "AMf"]), (["Bsubn", "Bn"], ["AM2n", "AMsn", "AMfn"])   # decreasing c*

_K, _R = C.controls, P._CONTROLS                          # K: misread_sim's stream (patched by droplet_ctrl_am2); R: droplet_range_common's


def _controls(rng, reps, eps, mode):
    return (_K if os.environ.get("AM2F_WORLD") == "K" else _R)(rng, reps, eps, mode)


C.controls = _controls                                    # read at call time, so each Pool (created inside `world`) inherits the right world


@contextmanager
def world(w):
    os.environ["AM2F_WORLD"] = w
    try:
        yield
    finally:
        os.environ.pop("AM2F_WORLD", None)


def say(s=""):
    print(s)


def tick(msg):
    print(f"[{time.time() - T0:5.0f}s] {msg}", file=sys.stderr, flush=True)


def stage(name, key, make):
    f = WORK / "ck" / f"{name}.pkl"
    if f.exists():
        k, r = pickle.loads(f.read_bytes())
        if k == key:
            tick(f"{name}: checkpoint")
            return r
    r = make()
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_bytes(pickle.dumps((key, r)))
    tick(f"{name}: simulated")
    return r


def run(cells, reps, rules, cstar, label):
    per = max(1, 8000 // reps)                            # droplet_range_a.run's batch size
    out = C.pmap([dict(cells=b, reps=reps, rules=rules, design="base", V0=None, cstar=cstar) for b in C.batches(cells, per)], PROCS, label)
    return [o for b in out for o in b]


# ---------------------------------------------------------------- simulation stages
def stage_k():
    nulls = C.null_cells(PROCS)
    assert len(nulls) == 87
    steps = [(f"probit s={s}", C.probit(0, s), 0 * C.LF) for s in SVALS]
    tasks = [dict(cells=[(n, Pc, p0, cv, e, "xor", 0.0, [5, i, j]) for j, (cv, e) in enumerate(CONDS)], reps=REPS, rules=KRULES, design="base", V0=None,
                  cstar=KC, want_v0=False, exact=True) for i, (n, Pc, p0) in enumerate(nulls + steps)]
    with world("K"):
        res = C.pmap(tasks, PROCS, "K")
    return {r: np.array([[res[i][j]["counts"][r] for j in range(len(CONDS))] for i in range(len(tasks))]) for r in KRULES}      # (91, 12, 5)


def rcells():
    null, excl, floor = C.null_cells(PROCS), C.exclusion_cells(PROCS), C.floor_cells()
    assert len(null) == 87 and len(floor) == 27
    allnull = null + excl + floor
    n0, n1 = len(null), len(excl)
    grp = [("count", range(0, n0)), ("excl", range(n0, n0 + n1)), ("floor", range(n0 + n1, len(allnull)))]
    return allnull, grp


def stage_ev(allnull):
    with world("R"):
        return run(A.cells_of(allnull, C.CONDS8, 11), A.EVR, RRULES, RC, "R-ev")


def mtruths():
    return [(f"h={h} s={s}", C.probit(h, s), 0 * C.LF) for h in A.H7 for s in A.S4]


def stage_map():
    with world("R"):
        return run(A.cells_of(mtruths(), A.MAPCONDS, 12), A.MAPR, RRULES, RC, "R-map")


def mapo(mp, c, h, s):
    return mp[(A.H7.index(h) * len(A.S4) + A.S4.index(s)) * len(A.MAPCONDS) + A.MAPCONDS.index(c)]


def arr(ev, r):
    return np.array([o["counts"][r] for o in ev])


# ---------------------------------------------------------------- printed forms (same f-strings as droplet_range_a.py)
def false_step(ev, r, allnull, grp):
    nc = len(C.CONDS8)
    rows, top = [], {}
    for j, (cv, e) in enumerate(C.CONDS8):
        vals = []
        for g, rg in grp:
            s, i = max((A.sh(ev[i * nc + j], r), i) for i in rg)
            vals.append(s)
            if s > top.get(g, (-1,))[0]:
                top[g] = (s, allnull[i][0], (cv, e))
        rows.append(vals)
    return rows, top


def fs_lines(rows, top, r):
    t = [f"  {cv:.2f} {e:.2f}    |    " + " ".join(f"{v:5.3f}" for v in rows[j]) for j, (cv, e) in enumerate(C.CONDS8)]
    return t, [f"  {r} {g}: max {s:.3f} ({name}, CV {c[0]}, misread {c[1]})" for g, (s, name, c) in top.items()]


def map_block(mp, r):
    return [f"  h={h} s={s} " + A.line(mapo(mp, A.M, h, s), r) for h in A.H7 for s in A.S4]


def step_rows(mp, r):
    return [f"  CV {c[0]} misread {c[1]} rule {r}: " + " | ".join(f"h={h} " + " ".join(f"{A.sh(mapo(mp, c, h, s), r):.2f}" for s in A.S4) for h in A.H7)
            for c in A.MAPCONDS]


def ns_rows(mp):
    return [f"  CV {c[0]} misread {c[1]}: " + " | ".join(f"h={h} " + " ".join(f"{A.sh(mapo(mp, c, h, s), 'B', 'NS'):.2f}" for s in A.S4) for h in A.H7)
            for c in A.MAPCONDS]


def rd(name):
    return [l.rstrip() for l in (HERE / name).read_text().splitlines()]


# ---------------------------------------------------------------- gates and checks
def trim(L):
    L = list(L)
    while L and not L[0]:
        L.pop(0)
    while L and not L[-1]:
        L.pop()
    return L


def diff(got, ref, first):
    bad = [(first + i, a, b) for i, (a, b) in enumerate(zip(got, ref)) if a != b]
    return len(got) == len(ref) and not bad, bad[:3]


def gate_k(cnt):
    buf = io.StringIO()
    with redirect_stdout(buf):
        P.part2({r: cnt[r] for r in ("B", "AM2", "Bn", "AM2n")})
    got, L = trim(l.rstrip() for l in buf.getvalue().splitlines()), rd("droplet_ctrl_am2.log")
    i0, i1 = next(i for i, l in enumerate(L) if l.startswith("PART 2.")), next(i for i, l in enumerate(L) if l.startswith("STAGE 3"))
    ok, bad = diff(got, trim(L[i0:i1]), i0 + 1)
    return ok, f"{len(got)} printed lines against droplet_ctrl_am2.log:{i0 + 1}-{i1}; first differences {bad}"


def gate_k3():
    buf = io.StringIO()
    with redirect_stdout(buf):
        match, tot = P.stage3()
    got, L = trim(l.rstrip() for l in buf.getvalue().splitlines()), rd("droplet_ctrl_am2.log")
    i0, i1 = next(i for i, l in enumerate(L) if l.startswith("STAGE 3")), next(i for i, l in enumerate(L) if l.startswith("B reproduces bak_policy_b.log"))
    ok, bad = diff(got, L[i0:i1 + 1], i0 + 1)
    return ok and match == tot, f"{len(got)} printed lines against droplet_ctrl_am2.log:{i0 + 1}-{i1 + 1}; B reproduces bak_policy_b.log {match}/{tot}; first differences {bad}", got


def gate_r(ev, mp, allnull, grp):
    L = rd("droplet_range_a.log")
    ip = next(i for i, l in enumerate(L) if l.startswith("Plateau map"))
    ifs = next(i for i, l in enumerate(L) if l.startswith("False STEP"))
    inb = next(i for i, l in enumerate(L) if l.startswith("NOT SCORABLE share"))
    tab = [l for l in L[ifs + 2:ip] if re.match(r"\s+0\.\d\d 0\.\d\d\s+\|", l)]
    assert len(tab) == len(C.CONDS8)
    bad, n = [], 0
    for r, col in (("B", 1), ("AM2", 3)):
        rows, top = false_step(ev, r, allnull, grp)
        t, mx = fs_lines(rows, top, r)
        for j, l in enumerate(t):
            n += 1
            if tab[j].split("|")[col].split() != l.split("|")[1].split():
                bad.append(("false-STEP row", r, tab[j][:20]))
        for l in mx:
            n += 1
            if l not in L[ifs:ip]:
                bad.append(("false-STEP max", l))
        i = next(i for i in range(ip, len(L)) if L[i] == f"rule {r}")
        for k, (a, b) in enumerate(zip(map_block(mp, r), L[i + 1:i + 29])):
            n += 1
            if a != b:
                bad.append(("plateau line", r, i + 2 + k, a, b))
        for l in step_rows(mp, r):
            n += 1
            if l not in L[ip:inb]:
                bad.append(("step row", l[:40]))
    for l in ns_rows(mp):
        n += 1
        if l not in L[inb:]:
            bad.append(("NS row", l[:40]))
    for name, mine, f in (("a_ev.pkl", ev, "a_ev.pkl"), ("a_map.pkl", mp, "a_map.pkl")):
        p = SCR / "droplet_range" / "cache" / f
        if p.exists():
            ref = pickle.loads(p.read_bytes())
            same = len(ref) == len(mine) and all((a["counts"][r] == b["counts"][r]).all() for a, b in zip(mine, ref) for r in ("B", "AM2"))
            n += 1
            if not same:
                bad.append((name, "counts differ"))
        else:
            say(f"  ({name} absent: not compared)")
    return not bad, f"{n} comparisons (log lines {ifs + 1}-{len(L)} of droplet_range_a.log, plus the pickles that exist); differences {bad[:4]}"


def cstar_check():
    """c* literals against their recomputation from the arrays (D0)."""
    ok, out = True, []
    fs = [SCR / "am2_gate" / f"cal_{i}.npy" for i in range(87)]
    if all(f.exists() for f in fs):
        cal = np.stack([np.load(f) for f in fs])           # [cell, n_c (128, 480), rule (B, AM2), condition, (eligible, p5)]
        sub = [j for j, (cv, e) in enumerate(CONDS) if e in (0.0, 0.06, 0.10)]
        for lab, v, lit in (("AM2 full grid", np.nanmin(cal[:, 0, 1, :, 1]), C_AMF), ("AM2 subset {0,.06,.10}", np.nanmin(cal[:, 0, 1, :, 1][:, sub]), C_AMS),
                            ("B subset {0,.06,.10}", np.nanmin(cal[:, 0, 0, :, 1][:, sub]), C_BSUB)):
            good = abs(float(v) - lit) < 1e-12
            ok &= good
            out.append(f"{lab}: recomputed {float(v)!r} literal {lit!r} {'equal' if good else 'DIFFERENT'}")
    else:
        out.append("am2_gate/cal_*.npy absent: AMf, AMs, Bsub literals not recomputed")
    p = SCR / "droplet_range" / "cache" / "a_cal.pkl"
    if p.exists():
        v = A.pick(pickle.loads(p.read_bytes()), "AM2", 114, 8)[0]
        good = abs(v - C_AM2_R) < 1e-12
        ok &= good
        out.append(f"AM2 of droplet_range_a (87 + 27 floor cells): recomputed {v!r} literal {C_AM2_R!r} {'equal' if good else 'DIFFERENT'}")
    else:
        out.append("a_cal.pkl absent: AM2 of droplet_range_a not recomputed")
    return ok, out


def structure(cnt, fams):
    """Exact consequences of classify(): same runs, lower c* -> no more STEP; STEP + SMOOTH, FLAT, RISING, NS unchanged; cap STEP <= no-cap STEP."""
    bad = []
    for fam in fams:
        for a, b in zip(fam, fam[1:]):
            if not (cnt[a][..., 0] >= cnt[b][..., 0]).all():
                bad.append(("STEP order", a, b))
        for r in fam[1:]:
            if not ((cnt[r][..., 2:] == cnt[fam[0]][..., 2:]).all() and (cnt[r][..., :2].sum(-1) == cnt[fam[0]][..., :2].sum(-1)).all()):
                bad.append(("STEP+SMOOTH / FLAT / RISING / NS differ", r, fam[0]))
    for r in cnt:
        if r + "n" in cnt and not (cnt[r][..., 0] <= cnt[r + "n"][..., 0]).all():
            bad.append(("cap STEP > no-cap STEP", r))
    return bad


# ---------------------------------------------------------------- K-world Section 10
SUMMARY = P.part2.__code__ and [("s<=0.15, misread <= 0.06", (0.1, 0.15), 0, 0.06), ("s<=0.15, misread 0.08-0.10", (0.1, 0.15), 0.08, 0.10),
                                ("s=0.2,  misread <= 0.06", (0.2,), 0, 0.06), ("s=0.2,  misread 0.08-0.10", (0.2,), 0.08, 0.10),
                                ("s=0.3,  misread <= 0.06", (0.3,), 0, 0.06), ("s=0.3,  misread 0.08-0.10", (0.3,), 0.08, 0.10)]     # droplet_ctrl_am2.part2's rows
ENDL = ["counting max"] + [f"{lab.strip()} {end}" for lab, *_ in SUMMARY for end in ("low", "high")]


def ends(cnt, r):
    """13 exact shares: the counting maximum over the 87 cells and 12 conditions, then (low, high) of the six summary rows."""
    sh = cnt[r][:, :, 0] / REPS
    out = [sh[:87].max()]
    for _, sv, lo, hi in SUMMARY:
        out += list(P.band(sh, sv, lo, hi))
    return np.array(out)


def pr(x, i):
    return float(f"{100 * x:.1f}" if i == 0 else f"{100 * x:.0f}")      # as printed in the Section 10 summary


def sec10(cnt, rules, title):
    say(f"\n{title}  (percent of runs called STEP; counting max over the 87 cells and 12 conditions, then low-high over the cells of each row)")
    say(f"  {'':27}" + "".join(f"{r:>11}" for r in rules))
    E = {r: ends(cnt, r) for r in rules}
    say(f"  {'counting max, all 12':27}" + "".join(f"{pr(E[r][0], 0):>11}" for r in rules))
    for k, (lab, *_) in enumerate(SUMMARY):
        say(f"  {lab:27}" + "".join(f"{int(pr(E[r][1 + 2 * k], 1))}-{int(pr(E[r][2 + 2 * k], 1))}".rjust(11) for r in rules))


def step_matrix(cnt, a, b, title):
    say(f"\n{title}: STEP share of {a} minus {b}, probit cells (rows: CV, true misread; columns: s), percentage points; {b} share in brackets")
    say("  CV   misread " + "".join(f"{'s=' + str(s):>16}" for s in SVALS))
    for j, (cv, e) in enumerate(CONDS):
        say(f"  {cv:.2f} {e:.2f}    " + "".join(f"{100 * (cnt[a][87 + k, j, 0] - cnt[b][87 + k, j, 0]) / REPS:+7.1f} ({100 * cnt[b][87 + k, j, 0] / REPS:3.0f})" for k in range(4)))


# ---------------------------------------------------------------- Bak
def bak_gate():
    """droplet_am2_gate.log Part 3 (statistic `sharp` has no c*): tables re-averaged, compared with the paired line and with the predictions file."""
    L = rd("droplet_am2_gate.log")
    iS, iH = next(i for i, l in enumerate(L) if l.startswith("PART 3.  SMOOTH")), next(i for i, l in enumerate(L) if l.startswith("PART 3.  SHARP"))
    iP = next(i for i, l in enumerate(L) if l.startswith("Paired number"))
    paired = {l.split()[0]: [float(x) for x in re.findall(r"\d\.\d{3}", l)] for l in L[iP + 1:iP + 3]}
    out, ok = {}, True
    for tag, i0 in (("smooth", iS), ("sharp", iH)):
        T = {}
        for l in L[i0 + 2:i0 + 8]:
            t = l.split()
            T[(int(t[0]), t[1])] = np.array(t[2:], float)
        for r in ("B", "AM2"):
            m = np.mean([T[(n, r)] for n in (3, 6, 10)], 0)
            out[(tag, r)] = (round(float(m[:3].mean()), 3), round(float(m[3:].mean()), 3))
            ok &= paired[r][(0 if tag == "smooth" else 2):(2 if tag == "smooth" else 4)] == list(out[(tag, r)])
    want = {("smooth", "B"): (0.866, 0.374), ("smooth", "AM2"): (0.494, 0.084), ("sharp", "B"): (0.016, 0.022), ("sharp", "AM2"): (0.008, 0.012)}
    return ok and out == want, f"droplet_am2_gate.log:{iS + 1}-{iP + 3}: re-averaged tables {out}, paired line {paired}", iP + 1


def bak_real():
    from blind_units import BAK_MAIN, BAK_TINC90
    from exclusion_bak import bak
    cs = {"AM2 -2.40": C_AM2_K, "AMs -2.4634": C_AMS, "AMf -2.9763": C_AMF}
    rows = []
    for name, series in (("Tinc 90 C", BAK_TINC90), ("main", BAK_MAIN)):
        for n1 in (3, 6, 10):
            x, k, n = bak(series, n1)
            K, N = k[None], n[None]
            F = C.Fitter(x, True)
            z, ns = C.trend_z(K, N, x), ~C.scorable_v(K, N)
            for e in (0.0, 0.02, 0.04, 0.06, 0.08, 0.10):
                f = F.fit(K, N, np.array([e]), ("c0", "cH", "s0", "sH"))
                sa, sb = C.rule_stats(f, RULES["AM2"][1]), C.rule_stats(f, RULES["B"][1])
                lab = lambda s, c: str(classify(dict(L=s["L"], rng=s["rng"], z=z, ns=ns), c)[0])
                rows.append((name, n1, e, float(sb["L"][0]), lab(sb, C.C_B), float(sa["L"][0]), {c: lab(sa, v) for c, v in cs.items()}))
    return rows


# ---------------------------------------------------------------- scoring
SCORES = []


def score(name, ok, detail):
    SCORES.append((name, bool(ok)))
    say(f"  {name} {'HELD' if ok else 'FAILED'}: {detail}")


def bands(label, vals, bs, nd0=1):
    """Endpoint bands (percent): the printed value (counting max one decimal, others integer) against [lo, hi]; returns (hits, text)."""
    hit, txt = 0, []
    for i, (v, (pt, lo, hi)) in enumerate(zip(vals, bs)):
        p = pr(v, i)
        h = lo <= p <= hi
        hit += h
        txt.append(f"    {ENDL[i]:34} {p:6.1f}  band [{lo}, {hi}] point {pt}  {'hit' if h else 'MISS'}  (exact {100 * v:.2f})")
    say(f"  {label}: {hit} of {len(vals)} endpoints in their bands")
    for t in txt:
        say(t)
    return hit == len(vals)


AMF_CAP = [(2.9, 1.7, 4.0), (83, 78, 88), (100, 99, 100), (28, 22, 33), (73, 68, 77), (74, 68, 79), (100, 99, 100), (24, 18, 29), (63, 57, 69), (38, 33, 43), (96, 92, 98), (7, 4, 10), (31, 26, 36)]
AMF_NC = [(2.9, 1.7, 4.0), (87, 82, 91), (100, 99, 100), (76, 70, 81), (93, 89, 96), (79, 74, 84), (100, 99, 100), (62, 56, 67), (86, 80, 90), (41, 36, 45), (97, 94, 98), (33, 28, 37), (50, 45, 53)]


def drop(cnt, a, b):
    """Largest endpoint drop of a against b over the 12 row endpoints, percentage points (exact shares), and the endpoint."""
    d = ends(cnt, b)[1:] - ends(cnt, a)[1:]
    return 100 * max(0.0, d.max()), ENDL[1 + int(d.argmax())], bool((ends(cnt, a)[1:] >= ends(cnt, b)[1:] - 0.05 - 1e-12).all())


def self_check():
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    assert hashlib.sha256((HERE / "droplet_am2_full_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    C.self_check()
    assert all(P.CSTAR[r] == KC[r] for r in P.CSTAR) and set(KRULES) >= {"B", "AM2", "Bn", "AM2n"} and "AM2" in RRULES
    a, b = np.random.default_rng(0), np.random.default_rng(0)
    with world("K"):
        assert (C.controls(a, 50, 0.05, "xor") == (b.random((50, 192)) < 0.05)[:, :128].sum(1)).all()
    a, b = np.random.default_rng(0), np.random.default_rng(0)
    with world("R"):
        assert (C.controls(a, 50, 0.05, "xor") == (b.random((50, 128)) < 0.05).sum(1)).all()
    Pt = C.probit(0, 0.15)
    for w in ("K", "R"):                                  # every threshold on the same runs: the structure holds, the counts add up
        with world(w):
            o = C.run_batch(dict(cells=[("t", Pt, 0 * Pt, 0.05, 0.08, "xor", 0.0, [9, 9, 9])], reps=60, rules=KRULES, design="base", V0=None, cstar=KC,
                                 want_v0=False, exact=True))[0]["counts"]
        cn = {r: v[None] for r, v in o.items()}
        assert all(v.sum() == 60 for v in o.values()) and not structure(cn, FAMS_CAP + FAMS_NC), structure(cn, FAMS_CAP + FAMS_NC)
    bad = {"B": np.array([[5, 90, 3, 0, 2]]), "Bn": np.array([[4, 90, 3, 0, 3]]), "Bsub": np.array([[4, 91, 3, 0, 2]]), "Bsubn": np.array([[4, 91, 3, 0, 2]])}
    assert structure(bad, [["Bsub", "B"]])                # a violated order is caught
    assert ends({r: np.zeros((91, 12, 5)) for r in ("B",)}, "B").shape == (13,) and len(ENDL) == 13 and len(AMF_CAP) == len(AMF_NC) == 13


def main():
    assert hashlib.sha256((HERE / "droplet_am2_full_predictions.txt").read_bytes()).hexdigest() == PRED_SHA
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    if "--check" in sys.argv:
        print("droplet_am2_full self-check passed")
        return
    WORK.mkdir(parents=True, exist_ok=True)
    say(f"Unit droplet_am2_full.  Model study on simulated drops, not data.  Predictions sha256 {PRED_SHA} (droplet_am2_full_predictions.txt, asserted).")
    say("Rule AM2 (plateau h in both families, cap) and rule B at several thresholds on the same simulated runs; nothing existing is edited.  NOT APPLIED: every wording proposal.")
    say(f"Thresholds: B {C.C_B}; Bsub {C_BSUB!r} (B recalibrated on the subset, informational); AM2 {C_AM2_K} (K-world, droplet_ctrl_am2) and {C_AM2_R!r} (R-world, droplet_range_a);")
    say(f"AMs {C_AMS!r} (AM2, subset misread {{0, .06, .10}}); AMf {C_AMF!r} (AM2, full 12-condition grid).  Suffix n = no control cap.  Processes {PROCS}.")

    # ---- simulation
    key = repr((sorted(KC.items()), sorted(RC.items()), REPS, A.EVR, A.MAPR))
    cnt = stage("K", key, stage_k)
    allnull, grp = rcells()
    ev = stage("EV", key, lambda: stage_ev(allnull))
    mp = stage("MAP", key, stage_map)

    # ---- gates and structure
    say("\nGATES AND CHECKS")
    cok, ctext = cstar_check()
    say("  c* literals:" + "".join(f"\n    {t}" for t in ctext))
    gk, tk = gate_k(cnt)
    say(f"  GATE K (droplet_ctrl_am2.part2 text, B/AM2/Bn/AM2n columns): {'PASS' if gk else 'FAIL'}: {tk}")
    g3, t3, bak_text = gate_k3()
    say(f"  GATE K3 (droplet_ctrl_am2.stage3 text): {'PASS' if g3 else 'FAIL'}: {t3}")
    gr, tr = gate_r(ev, mp, allnull, grp)
    say(f"  GATE R (droplet_range_a.log, B and AM2): {'PASS' if gr else 'FAIL'}: {tr}")
    bk, tb, ib = bak_gate()
    say(f"  GATE BAK (droplet_am2_gate.log paired numbers): {'PASS' if bk else 'FAIL'}: {tb}")
    cr = {r: arr(ev, r) for r in RRULES}
    cmap = {r: np.array([[mapo(mp, c, h, s)["counts"][r] for c in A.MAPCONDS for s in A.S4] for h in A.H7]) for r in RRULES}
    sK = structure(cnt, FAMS_CAP + FAMS_NC)
    sR = structure(cr, [["Bsub", "B"], ["AM2", "AMs", "AMf"]])
    sM = structure({r: v for r, v in cmap.items()}, [["Bsub", "B"], ["AM2", "AMs", "AMf"]])
    say(f"  STRUCTURE (STEP order in c*, STEP + SMOOTH / FLAT / RISING / NS equal across thresholds, cap <= no cap): K-world violations {sK}, R false-STEP runs {sR}, R plateau map {sM}")
    D0 = cok and gk and g3 and gr
    say(f"  D0: {'gates hold and c* literals equal their recomputation -> results trusted' if D0 else 'NOT TRUSTED (a gate failed or a c* literal differs)'}")

    # ---- R-world plateau map
    say("\nSECTION 1.  R-world plateau map (droplet_range_a stage 3 seeds [12, i, j], 600 runs per cell): truth P = h + (1 - h) Phi(-ln(V / 3.1 uL) / s)")
    say("Verdict shares, main condition (CV 0.05, misread 0.03)")
    for r in ("AMf", "AMs"):
        say(f"rule {r} (c* {RC[r]!r})")
        for l in map_block(mp, r):
            say(l)
    say("\nSTEP share by plateau h (rows) and s = " + ", ".join(map(str, A.S4)) + " (columns), every condition (B and AM2 lines equal droplet_range_a.log:151-168)")
    for c in A.MAPCONDS:
        for r in ("B", "Bsub", "AM2", "AMs", "AMf"):
            say(step_rows(mp, r)[A.MAPCONDS.index(c)])
    say("\nNOT SCORABLE share, rule B (rows as above; equal to droplet_range_a.log:169-173)")
    for l in ns_rows(mp):
        say(l)
    st = lambda r, c, h, s: A.sh(mapo(mp, c, h, s), r)
    # ---- R-world false STEP
    say("\nSECTION 2.  R-world false STEP (droplet_range_a stage 2, seeds [11, i, j], 400 runs): highest STEP share over the cells of each group, by condition")
    FS = {}
    for r in ("B", "Bsub", "AM2", "AMs", "AMf"):
        rows, top = false_step(ev, r, allnull, grp)
        FS[r] = top
        say(f"rule {r} (c* {RC[r]!r})")
        say("  CV   misread |   count  excl floor")
        t, mx = fs_lines(rows, top, r)
        for l in t:
            say(l)
        for l in mx:
            say(l)
    # ---- K-world
    say("\nSECTION 3.  K-world Section 10 (droplet_ctrl_am2.stage2 seeds [5, i, j], 500 runs, 128 controls)")
    sec10(cnt, ["B", "Bsub", "AM2", "AMs", "AMf"], "With the control cap (prereg Section 2 item 1 as written)")
    sec10(cnt, ["Bn", "Bsubn", "AM2n", "AMsn", "AMfn"], "Without the cap (what the prereg Sections 3 and 10 print)")
    step_matrix(cnt, "AMf", "B", "Cap")
    step_matrix(cnt, "AMfn", "Bn", "No cap")
    step_matrix(cnt, "AMs", "B", "Cap")
    for r in ("B", "AM2", "AMs", "AMf", "Bsub"):
        say(f"\n  K-world counting maximum by condition, {r}: " + " ".join(f"{cnt[r][:87, j, 0].max() / REPS:.3f}" for j in range(len(CONDS))))
    # ---- Bak
    say("\nSECTION 4.  Bak et al. (2023) Fig. 2a re-fit (POST HOC).  `sharp` has no c*, so droplet_am2_gate.log Part 3 (500 runs, seeds [31, a, b, tr]) is not rerun.")
    say(f"  {'OK' if bk else 'MISMATCH'}: {tb}")
    rows = bak_real()
    say("  Real series, rule AM2 (L does not depend on c*) and labels at each threshold; rule B at -2.21 equals bak_policy_b.log (gate K3):")
    gb = [re.findall(r"e (\d\.\d\d): L ([+-]\d+\.\d\d) sharp\s+([\d.]+) (\w+)", l) for l in bak_text if l.startswith(("Tinc", "main"))]
    same_bak = True
    for name in ("Tinc 90 C", "main"):
        sel = [r for r in rows if r[0] == name]
        labs = {c: sorted({r[6][c] for r in sel}) for c in sel[0][6]}
        Ls = [r[5] for r in sel]
        say(f"    {name:9}: AM2 L in [{min(Ls):+.2f}, {max(Ls):+.2f}] over 18 (n, misread) pairs; labels " + "; ".join(f"{c}: {v}" for c, v in labs.items()) + f"; B -2.21: {sorted({r[4] for r in sel})}")
    for i, (name, n1) in enumerate([(a, b) for a in ("Tinc 90 C", "main") for b in (3, 6, 10)]):
        ref_b, ref_a = gb[2 * i], gb[2 * i + 1]
        sel = [r for r in rows if r[0] == name and r[1] == n1]
        same_bak &= all(f"{r[3]:+.2f}" == q[1] and r[4] == q[3] for r, q in zip(sel, ref_b))
        same_bak &= all(f"{r[5]:+.2f}" == q[1] and r[6]["AM2 -2.40"] == q[3] for r, q in zip(sel, ref_a))
    say(f"  my re-computation of L and the labels at B -2.21 / AM2 -2.40 equals the stage-3 text of droplet_ctrl_am2.log for all 36 x 2 rows: {same_bak}")
    tinc_ok = all(r[6][c] == "SMOOTH" for r in rows if r[0] == "Tinc 90 C" for c in r[6])
    main_ok = all(r[6][c] == "NS" for r in rows if r[0] == "main" for c in r[6])

    # ---- scoring
    say("\nSCORING (droplet_am2_full_predictions.txt, as written; bands are for the printed values)")
    bF = {h: st("B", A.M, h, 0.1) for h in (0.2, 0.25, 0.3)}
    f1 = []
    rep = [float(round(st("AM2", A.M, h, 0.1), 2)) for h in (0.2, 0.25, 0.3)]
    rb = [float(round(bF[h], 2)) for h in (0.2, 0.25, 0.3)]
    f1.append(rep == [0.98, 0.96, 0.94] and rb == [0.73, 0.30, 0.05])
    say(f"  F1 reproduction at -2.4034 (droplet_range_a.log:131, 135, 139; B 73, 77, 81): AM2 {rep} (want 0.98/0.96/0.94), B {rb} (want 0.73/0.30/0.05): {f1[-1]}")
    for r, bs in (("AMf", {0.2: (0.92, 0.98), 0.25: (0.87, 0.96), 0.3: (0.80, 0.93)}), ("AMs", {0.2: (0.96, 0.98), 0.25: (0.93, 0.96), 0.3: (0.90, 0.94)})):
        for h, (lo, hi) in bs.items():
            v = st(r, A.M, h, 0.1)
            f1.append(lo - 1e-9 <= v <= hi + 1e-9)
            say(f"    {r} h={h} s=0.1: {v:.4f} band [{lo}, {hi}] {'hit' if f1[-1] else 'MISS'}  (B {bF[h]:.4f}, AM2 -2.4034 {st('AM2', A.M, h, 0.1):.4f})")
    gap = [st("AMf", A.M, 0.2, 0.1) - bF[0.2], st("AMf", A.M, 0.25, 0.1) - bF[0.25]]
    f1.append(gap[0] >= 0.20 - 1e-9 and gap[1] >= 0.55 - 1e-9)
    say(f"    AMf - B at h=0.2 / 0.25 (s=0.1): {gap[0]:.4f} / {gap[1]:.4f} (want >= 0.20 / >= 0.55): {f1[-1]}")
    score("F1", all(f1), f"{sum(f1)} of {len(f1)} items")

    cm = {r: FS[r]["count"][0] for r in FS}
    f2 = [abs(cm["AM2"] - 0.045) < 1e-9, 0.030 - 1e-9 <= cm["AMs"] <= 0.045 + 1e-9, 0.012 - 1e-9 <= cm["AMf"] <= 0.035 + 1e-9]
    kcm = {r: cnt[r][:87, :, 0].max() / REPS for r in ("B", "AM2", "AMs", "AMf", "Bsub", "AMfn")}   # STEP column only
    f2 += [abs(kcm["AM2"] - 0.052) < 1e-9, 0.038 - 1e-9 <= kcm["AMs"] <= 0.052 + 1e-9, 0.017 - 1e-9 <= kcm["AMf"] <= 0.040 + 1e-9]
    score("F2", all(f2), f"R-world counting max AM2 {cm['AM2']:.3f} (want 0.045), AMs {cm['AMs']:.3f} ([0.030, 0.045]), AMf {cm['AMf']:.3f} ([0.012, 0.035]); "
          f"K-world AM2 {kcm['AM2']:.3f} (want 0.052; droplet_ctrl_am2.log:154), AMs {kcm['AMs']:.3f} ([0.038, 0.052]), AMf {kcm['AMf']:.3f} ([0.017, 0.040]); per item {[bool(x) for x in f2]}")
    say(f"    informational: Bsub counting max R {cm['Bsub']:.3f} (predicted [0.045, 0.065]), K {kcm['Bsub']:.3f} ([0.054, 0.075]); B R {cm['B']:.3f}, K {kcm['B']:.3f}; "
        f"exclusion max AM2/AMs/AMf {FS['AM2']['excl'][0]:.3f}/{FS['AMs']['excl'][0]:.3f}/{FS['AMf']['excl'][0]:.3f}, floor max {FS['AM2']['floor'][0]:.3f}/{FS['AMs']['floor'][0]:.3f}/{FS['AMf']['floor'][0]:.3f}")

    say("\n  F3.  Section 10 endpoints (percent as printed; exact in brackets)")
    f3a = bands("AMf, cap", ends(cnt, "AMf"), AMF_CAP)
    f3b = bands("AMfn, no cap", ends(cnt, "AMfn"), AMF_NC)
    f3s = []
    for a, b in (("AMs", "AM2"), ("AMsn", "AM2n")):
        ea, eb = ends(cnt, a), ends(cnt, b)
        okm = bool((ea <= eb + 1e-12).all() and (ea >= eb - 0.03 - 1e-12).all())
        f3s.append(okm)
        say(f"  {a} against {b}: every one of 13 endpoints <= and within 3 points: {okm}; largest gap {100 * (eb - ea).max():.2f} points at {ENDL[int((eb - ea).argmax())]}")
    dr = {}
    for a, b, bd in (("AMf", "B", (4, 11)), ("AMfn", "Bn", (3, 10)), ("AMs", "B", (0, 3)), ("AMsn", "Bn", (0, 3))):
        dr[a] = drop(cnt, a, b)
        say(f"  largest endpoint drop {a} against {b}: {dr[a][0]:.2f} points at {dr[a][1]} (band {bd}); every endpoint AM >= B - 0.05: {dr[a][2]}")
    f3d = [4 <= dr["AMf"][0] <= 11, 3 <= dr["AMfn"][0] <= 10, 0 <= dr["AMs"][0] <= 3, 0 <= dr["AMsn"][0] <= 3]
    score("F3", f3a and f3b and all(f3s) and all(f3d), f"AMf cap bands {f3a}, AMfn bands {f3b}, AMs within 3 points of AM2 {f3s}, largest drops in band {[bool(x) for x in f3d]}")

    score("F4", bk and tinc_ok and main_ok and same_bak, f"paired numbers equal the predictions {bk}; 18 Tinc 90 C AM2 labels SMOOTH at -2.40, -2.4634, -2.9763: {tinc_ok}; 18 main labels NS at every threshold: {main_ok}; "
          f"my L and labels equal the stage-3 text: {same_bak}")
    score("F5", gk and g3 and gr and not (sK or sR or sM), f"K text gate {gk}, K stage-3 gate {g3}, R digit gate {gr}, structure clean {not (sK or sR or sM)}")

    # ---- decision
    say("\nDECISION RULES (at -2.9763; the same tests at -2.4634 are information)")
    d1a = {r: st(r, A.M, 0.2, 0.1) for r in ("AMf", "AMs")}
    d1b = {r: (dr[r][2], dr[r + "n"][2], dr[r][0], dr[r + "n"][0]) for r in ("AMf", "AMs")}
    d2 = {r: (cm[r], kcm[r], cnt[r + "n"][:87, :, 0].max() / REPS) for r in ("AMf", "AMs")}
    D = {}
    for r in ("AMf", "AMs"):
        D[r] = dict(D1a=d1a[r] >= 0.90, D1b=d1b[r][0] and d1b[r][1], D2=max(d2[r][0], d2[r][1], d2[r][2]) <= 0.07)
        say(f"  {r} (c* {RC[r]!r}): D1a STEP(main, s=0.1, h=0.2) = {d1a[r]:.4f} >= 0.90: {D[r]['D1a']}; D1b cap pair {d1b[r][0]} (largest drop {d1b[r][2]:.2f} points), "
            f"no-cap pair {d1b[r][1]} ({d1b[r][3]:.2f} points): {D[r]['D1b']}; D2 counting max R {d2[r][0]:.3f}, K cap {d2[r][1]:.3f}, K no cap {d2[r][2]:.3f} <= 0.07: {D[r]['D2']}")
    f = D["AMf"]
    if not D0:
        rec = "NOT TRUSTED: no recommendation (D0)"
    elif not f["D2"]:
        rec = "reject AM2 (D2 fails: false STEP above 0.07 at the full-grid threshold)"
    elif not f["D1a"]:
        rec = "keep B plus a limit sentence; AM2 does not fix the blind spot at the honest c* (D1a fails)"
    elif f["D1b"]:
        rec = "AM2 at -2.9763 is worth proposing as a replacement for B (wording NOT APPLIED)"
    else:
        rec = ("keep B; add a limit sentence naming the plateau blind spot with B's measured shares (h = 0.2 / 0.25 / 0.3); offer AM2 only as a labelled post hoc sensitivity analysis "
               "(D1a and D2 hold, D1b fails)")
    say(f"  RECOMMENDATION: {rec}")
    say(f"  Expectation check (predictions file): D1a holds {f['D1a']} (95%), D2 holds {f['D2']} (95%), D1b fails at -2.9763 {not f['D1b']} (80%), D1b holds at -2.4634 {D['AMs']['D1b']} (85%)")
    say(f"\nHELD {sum(v for _, v in SCORES)} of {len(SCORES)}: " + ", ".join(f"{n} {'HELD' if v else 'FAILED'}" for n, v in SCORES))
    tick("done")


if __name__ == "__main__":
    main()
