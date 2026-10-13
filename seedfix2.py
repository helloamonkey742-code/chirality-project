"""seedfix2.py: remove the Enceladus seed collision from the four scripts seedfix.py did not regenerate.

THE BUG (as in seedfix.py). np.random.default_rng([s, 0]) == default_rng(s) (SeedSequence drops trailing zeros), so
world 0 (Enceladus) drew its extra inputs (eta, r, rho) from the stream of its own base draws. seedfix.py fixed
sweep2, sweep3, multiwall Section 4, rho2side and rerun_base. Four later scripts repeat the same draw:
  rerun_base.py:170 (EX) and :242, front2d.py:201 (worlds/shares), deephand_scale.py:299 (section3),
  pin_onset.py:254 (section4).

THE FIX: the seedfix.py wrapper, unchanged. np.random.default_rng is wrapped before any project import so that the exact
two-element seed [s, 0] becomes [s, 0, 1]; every other seed passes through. Worlds 1-14 keep their streams, so every
non-Enceladus field must reproduce the old logs to the printed digit (exact regression). No project script is edited.

PER SCRIPT
  rerun_base.py  SKIPPED. seedfix.py Section F already reran all of its result sections (0-11) under this wrapper
                 against rerun_base.log (328/328 lines, the Enceladus rows and Enceladus-including aggregates changed,
                 no other world moved), and the manuscript already quotes those values (E55, E56 cite seedfix.log
                 Section F). Rerunning would repeat Section F byte for byte.
  front2d.py     Part B (worlds()) only. Part A (1D/2D/3D front simulations) draws no world inputs. worlds() asserts
                 its front-off run against multiwall.log Section 4 read from the cwd, so it runs with the cwd set to a
                 temp dir holding a copy of multiwall.log whose six Enceladus lines are replaced by seedfix.log's
                 regenerated ones (each old line checked against the real log first). Compared with front2d.log:106-162.
  deephand_scale Section 3 only (Sections 1-2 are world-agnostic wrapper simulations; their ROW lines are read from
                 the real deephand_scale.log). HERE is pointed at the temp dir only for multiwall.log. Compared with
                 deephand_scale.log:702-900.
  pin_onset.py   Sections 2 and 3 rerun as-is (int seeds; they give the exact kappa, SE and band section 4 needs) and
                 section 4 rerun from its own source with ONE change: its hard-coded sweep2.log regression target for
                 Enceladus ((3.9, 0.0), (16.9, 0.0)), a seed-collided number, becomes seedfix.log's regenerated sweep2
                 step f 6.7/0.0 and x5 19.9/0.0 (seedfix.log:75 and :78). The source is compiled under pin_onset.py's
                 own file name and line numbers, so the wrapper reports the call site as pin_onset.py:254. Section 1 has
                 no world inputs. Compared with pin_onset.log:20-90.
  Timing brackets such as "[42s]" are stripped before comparing. A changed line is classified Enceladus (its changed
  fields name only Enceladus), aggregate (no world named) or VIOLATION (a field naming another world changed, or the
  line count differs).

MANUSCRIPT EXPOSURE (journals/ija_draft.md, read before writing this): no Enceladus level from these three scripts is
quoted. Aggregates that include Enceladus: front2d "no world's P_pved by more than 0.01 points ... at most 0.08 with a
front 5.5 times slower" (Section 4, assumption 8), "2% (Miranda) to 17% (Ganymede)" fronts; deephand "a correction
fitted to these runs changes no share by more than 1.4 points" (E58), "52-74% of decided draws have a deciding layer
under five wall widths" (H24), ocean-scale losses Mimas -5.8 to Umbriel -15.3; pin_onset kappa only (E48).

PREDICTIONS (written and committed before any result-producing run; never edited). "Old" = the logged value.
  P1  The wrapper remaps exactly these call sites and no others: front2d.py:201, deephand_scale.py:299,
      multiwall.py:435 (run4, called by front2d and deephand_scale), pin_onset.py:254.
  P2  Exact regression: front2d.log:106-162, deephand_scale.log:702-900 and pin_onset.log:20-90 have the same line
      counts as the regenerated text and no VIOLATION; pin_onset sections 2-3 (lines 20-51) are identical. front2d
      prints EXACT for both of its regressions, deephand_scale prints EXACT, and the patched pin_onset section 4
      regression passes.
  P3  Enceladus with a lid stays 0.0 in every regenerated final column.
  P4  front2d Enceladus drive row: % u < 3 (all draws) 4.2, the u quartiles 28.1/131/612, % u < 1 1.0 and t_mix 0.0
      are identical (they use only base draws); % u < 3 among draws reaching one hand moves by <= 1.5 from 5.0.
  P5  front2d Enceladus shares, front rule off = 26.4 / 9.8 / 63.9 | 0.0 / 0.0 / 100.0 (seedfix multiwall main),
      and front rule on equal to off at the printed 0.1.
  P6  front2d Enceladus P_pved changes: |main|, |u*=1|, |v x A| <= 0.01 and |v / A| <= 0.08 points, so "0.01" and
      "0.08" in the manuscript stand; the seed-1 maximum stays in [0.30, 0.45] (old 0.36).
  P7  deephand_scale Enceladus multiwall main 26.4 / 9.8 / 63.9 and x10 P_pved 19.4 (seedfix multiwall).
  P8  deephand_scale Enceladus fitted-h P_pved minus main is in [-1.4, -0.5] (old -1.2), so "1.4 points" stands.
  P9  deephand_scale Enceladus ocean-scale P_pved minus main is in [-15.3, -5.8] (old -11.4), so the quoted range of
      ocean-scale losses (Mimas to Umbriel) stands.
  P10 deephand_scale Enceladus thin-layer % (< 5 widths, no lid) is in [51.7, 73.6] (old 68.6) and < 1 width in
      [29.4, 50.8] (old 45.5), so H24's "52-74%" stands; the median h moves by <= 0.2 from 0.79.
  P11 pin_onset Enceladus step f at every kappa-based multiplier and x1 = 6.7 / 0.0, x5 = 19.9 / 0.0, three-way split
      at kappa 6.7 / 1.5 / 91.8 | 0.0 / 0.0 / 100.0 (seedfix sweep2 step f).
  P12 pin_onset Enceladus layered % 95.0 identical (base draws only); pinned@kappa within 0.5 of 84.7 and median x
      within 5% of 3.29 (rho was collided with the concentration draw, which is independent of the geometry that sets
      r_w/H, so the pinned fraction moves only by sampling noise).
  P13 No number quoted in journals/ija_draft.md changes: every aggregate listed under MANUSCRIPT EXPOSURE reproduces.

Run: OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2 nice -n 10 \
     /opt/anaconda3/bin/python seedfix2.py > seedfix2.log

RESULTS (seedfix2.log, 2026-10-13): N = 200000, 70 s. 12 held, 1 failed. (A smoke run at SEEDFIX2_N=2000 was made
before the predictions commit, only to catch crashes; its numbers are not used.)
  P1  HELD.   remapped front2d.py:201 x6, multiwall.py:435 x2, deephand_scale.py:299 x1, pin_onset.py:254 x1; nothing else.
  P2  HELD.   0 violations and 0 aggregate lines changed in all three logs (front2d 57/57, deephand_scale 199/199,
              pin_onset 71/71 lines; sections 2-3 identical); front2d EXACT x2, deephand_scale EXACT, pin_onset
              section 4 regression passed against seedfix's step f 6.7 and x5 19.9.
  P3  HELD.   all 13 regenerated Enceladus final values are 0.0.
  P4  HELD.   drive row 4.2 | 5.0 -> 4.4 | 28.1 131 612 | 1.0 | 0.0 (only the reaching-one-hand column moved).
  P5  HELD.   front rule off = on = 26.4 / 9.8 / 63.9 | 0.0 / 0.0 / 100.0 (was 28.3 / 8.8 / 63.0).
  P6  HELD.   Enceladus changes -0.00 / +0.00 (main, u*=1, v x A), v / A -0.04 -> -0.03; seed-1 max 0.36 unchanged.
  P7  HELD.   deephand main 26.4 / 9.8 / 63.9, x10 19.4 / 6.4 / 74.1.
  P8  FAILED: fitted h 27.1 -> 26.0, shift -1.2 -> -0.4 (registered [-1.4, -0.5]); smaller than predicted, so the
              manuscript's "no share by more than 1.4 points" still stands (Ganymede and Titan set it).
  P9  HELD.   ocean-scale 16.9 / 3.1 / 80.0 -> 17.4 / 2.6 / 80.0; loss vs main -11.4 -> -9.0.
  P10 HELD.   thin 68.6 -> 69.2, < 1 width 45.5 -> 46.6, median h 0.79 -> 0.69.
  P11 HELD.   pin_onset step f 3.9 -> 6.7 / 0.0 in all six kappa/x1 columns, x5 16.9 -> 19.9 / 0.0, split
              3.9 / 2.0 / 94.1 -> 6.7 / 1.5 / 91.8 | 0.0 / 0.0 / 100.0.
  P12 HELD.   layered 95.0 unchanged, pinned@kappa 84.7 -> 84.9 (84.85), median x 3.29 unchanged.
  P13 HELD.   front2d 0.01 / 0.08 / 2-17%, deephand 1.4 points, 52-74% thin, ocean-scale range all reproduce:
              no number in journals/ija_draft.md changes.
"""
import contextlib
import inspect
import io
import os
import re
import shutil
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
T0 = time.time()

# ---------------------------------------------------------------- the fix (copied from seedfix.py)
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


np.random.default_rng = fixed_rng   # before the project modules import

import uncertainty as U   # noqa: E402
SMOKE = "SEEDFIX2_N" in os.environ   # smoke test only: SEEDFIX2_N=2000 (front2d skipped: it asserts the N=200000 log)
U.N = int(os.environ.get("SEEDFIX2_N", U.N))
import sweep2 as S             # noqa: E402
import sweep3 as S3            # noqa: E402
import multiwall as M          # noqa: E402
import front2d as F2           # noqa: E402
import deephand_scale as DH    # noqa: E402
import pin_onset as PO         # noqa: E402

PO.SMOKE = SMOKE
N = U.N
ENC = "Enceladus"
NAMES = list(S.chains(0)[1])
assert NAMES.index(ENC) == 0, NAMES
TOK = [n for n in NAMES if n != "Early Earth ocean"] + ["Early", "Earth"]
TMP = tempfile.mkdtemp(prefix="seedfix2_")
RES = {}
TIMING = re.compile(r"\s*\[\d+s\]")
SEP = re.compile(r" \|\|? | ; |; |, ")


def check(p, ok, detail):
    RES[p] = (bool(ok), detail)
    print(f"  {p} {'HELD' if ok else 'FAILED'}: {detail}")


def cap(fn, *a, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn(*a, **kw)
    return r, buf.getvalue().splitlines()


def logl(f):
    return open(os.path.join(HERE, f)).read().splitlines()


def who(s):
    return {t for t in TOK if t in s}


DIFFS = {}


def report(name, old, start, new):
    """old: whole log; start: 0-based index of new[0] in old. Prints and stores every changed line by class."""
    d = {"enc": [], "agg": [], "bad": []}
    seg = old[start:start + len(new)]
    if len(seg) != len(new) or (start + len(new) < len(old) and old[start + len(new)].strip()
                                and not old[start + len(new)].startswith(("check:", "total"))):
        d["bad"].append((0, f"line count: {len(seg)} log lines vs {len(new)} regenerated", ""))
    for j, (o, n) in enumerate(zip(seg, new)):
        o, n = TIMING.sub("", o), TIMING.sub("", n)
        if o == n:
            continue
        names = who(o) | who(n)
        if ENC in names and names - {ENC}:   # several worlds on one line: changed segments must name only Enceladus
            so, sn = SEP.split(o), SEP.split(n)
            ch = set().union(*(who(a) | who(b) for a, b in zip(so, sn) if a != b)) if len(so) == len(sn) else None
            names = {ENC} if ch == {ENC} else None
        kind = "agg" if names is not None and not names else "enc" if names == {ENC} else "bad"
        d[kind].append((start + j + 1, o, n))
    DIFFS[name] = d
    print(f"\n--- {name}: lines {start + 1}-{start + len(new)} compared ({len(new)} regenerated); "
          f"{len(d['enc'])} Enceladus lines changed, {len(d['agg'])} aggregate lines changed, "
          f"{len(d['bad'])} VIOLATIONS ---")
    for kind, lab in (("enc", "Enceladus"), ("agg", "aggregate (no world named; includes Enceladus)"),
                      ("bad", "VIOLATION")):
        for ln, o, n in d[kind]:
            print(f"  [{lab}] {name}:{ln}\n    old: {o}\n    new: {n}")


def fl(s):
    return [float(v) for v in re.findall(r"[-+]?\d+\.?\d*(?:e[-+]?\d+)?|nan", s)]


def enc_rows(lines):
    return [l for l in lines if l.strip().startswith(ENC)]


print(f"seedfix2.py, N = {N}: Enceladus extra draws moved from default_rng([s, 0]) to default_rng([s, 0, 1]) in "
      f"front2d Part B, deephand_scale Section 3 and pin_onset Section 4; rerun_base skipped (seedfix.py Section F)")

# ---------------------------------------------------------------- patched multiwall.log (seedfix.log Section D)
mw = logl("multiwall.log")
sf = logl("seedfix.log")
pairs = [(int(m.group(1)), sf[i + 1][len("    old: "):], sf[i + 2][len("    new: "):])
         for i, l in enumerate(sf) if (m := re.match(r"\s+\[Enceladus\] multiwall\.log:(\d+)$", l))]
for ln, o, n in pairs:
    assert mw[ln - 1] == o, f"multiwall.log:{ln} does not match seedfix.log's old line"
    mw[ln - 1] = n
open(os.path.join(TMP, "multiwall.log"), "w").write("\n".join(mw) + "\n")
print(f"patched multiwall.log in a temp dir: lines {[p[0] for p in pairs]} replaced by seedfix.log's regenerated "
      f"Enceladus lines (old text checked against the real log)")

# ---------------------------------------------------------------- A. front2d Part B
print("\n=== A. front2d.worlds() (Part B) vs front2d.log:106-162 ===")
saved = {m: m.delta_at for m in (S, S3, M)}
newA = []
if SMOKE:
    print("smoke: skipped (worlds() asserts the N = 200000 multiwall.log)")
else:
    cwd = os.getcwd()
    os.chdir(TMP)
    try:
        _, newA = cap(F2.worlds)
    finally:
        os.chdir(cwd)
        for m, f in saved.items():   # worlds() rebinds delta_at in sweep2/sweep3/multiwall; restore the originals
            m.delta_at = f
    report("front2d.log", logl("front2d.log"), 105, newA)
print(f"[A: {time.time() - T0:.0f} s]")

# ---------------------------------------------------------------- B. deephand_scale Section 3
print("\n=== B. deephand_scale.section3() vs deephand_scale.log:702-900 ===")
rs = DH.rows()                 # ROW lines from the real log
DH.HERE = TMP                  # section3 now reads the patched multiwall.log
outB, newB = cap(DH.section3, rs)
report("deephand_scale.log", logl("deephand_scale.log"), 701, newB)
print(f"[B: {time.time() - T0:.0f} s]")

# ---------------------------------------------------------------- C. pin_onset sections 2-4
print("\n=== C. pin_onset sections 2, 3 and 4 (section 4 with the seedfix Enceladus regression target) "
      "vs pin_onset.log:20-90 ===")
src, first = inspect.getsourcelines(PO.section4)
src = "".join(src)
OLD_T, NEW_T = '"Enceladus": ((3.9, 0.0), (16.9, 0.0))', '"Enceladus": ((6.7, 0.0), (19.9, 0.0))'
assert src.count(OLD_T) == 1
ns = dict(PO.__dict__)
exec(compile("\n" * (first - 1) + src.replace(OLD_T, NEW_T), PO.__file__, "exec"), ns)
(k, se, _), l2 = cap(PO.section2)
(cells, _), l3 = cap(PO.section3)
band = (min(v[0] for v in cells.values()), max(v[0] for v in cells.values()))
outC, l4 = cap(ns["section4"], k, se, band)
newC = l2 + l3 + l4
report("pin_onset.log", logl("pin_onset.log"), 19, newC)
print(f"[C: {time.time() - T0:.0f} s]")

# ---------------------------------------------------------------- predictions
print("\n=== PREDICTIONS ===")
sites = {"front2d.py:201", "deephand_scale.py:299", "multiwall.py:435", "pin_onset.py:254"}
check("P1", set(REMAPPED) == sites, "remapped " + ", ".join(f"{a} x{b}" for a, b in sorted(REMAPPED.items())))
bad = {k_: len(v["bad"]) for k_, v in DIFFS.items()}
agg = {k_: [ln for ln, _, _ in v["agg"]] for k_, v in DIFFS.items()}
same23 = not any(ln <= 51 for kind in ("enc", "agg", "bad") for ln, _, _ in DIFFS["pin_onset.log"][kind])
regs = [l for l in newA + newB if l.startswith("regression")] + [l.strip() for l in l4 if "regression" in l]
check("P2", not any(bad.values()) and same23 and len(regs) == 3 and all("FAILED" not in l for l in regs)
      and sum(l.count("EXACT") for l in regs) == 3,
      f"violations {bad}; aggregate lines {agg}; pin_onset sections 2-3 identical {same23}; " + " | ".join(regs))
eA, eB = enc_rows(newA), enc_rows(newB)
fins = []
if eA:
    sh = fl(eA[1])
    fins += [sh[3], sh[9]]
fins += [outB[ENC][k_][1][0] for k_ in ("h1", "x10", "fit", "ocean")]
fins += [outC[ENC][lab][1][0] for lab in outC[ENC] if isinstance(outC[ENC][lab], list)]
check("P3", max(fins) == 0.0, f"max Enceladus final P_pved over {len(fins)} regenerated finals {max(fins):.3g}")
if eA:
    dr, sh, ch = fl(eA[0]), fl(eA[1]), fl(eA[2])
    check("P4", dr[0] == 4.2 and dr[2:] == [28.1, 131.0, 612.0, 1.0, 0.0] and abs(dr[1] - 5.0) <= 1.5,
          f"drive row {eA[0].split(None, 1)[1]} (old 4.2 | 5.0 | 28.1 131 612 | 1.0 | 0.0)")
    check("P5", sh[:6] == [26.4, 9.8, 63.9, 0.0, 0.0, 100.0] and sh[:6] == sh[6:],
          f"off {sh[:6]}, on {sh[6:]}")
    seed1 = float(re.search(r"= ([\d.]+) points", next(l for l in newA if l.startswith("seed 1"))).group(1))
    check("P6", max(abs(v) for v in ch[:6]) <= 0.01 and max(abs(v) for v in ch[6:]) <= 0.08 and 0.30 <= seed1 <= 0.45,
          f"changes main / u*=1 / vxA / v/A {ch}; seed-1 max {seed1}")
b = outB[ENC]
check("P7", tuple(round(v, 1) for v in b["h1"][0]) == (26.4, 9.8, 63.9) and round(b["x10"][0][0], 1) == 19.4,
      f"main {b['h1'][0][0]:.1f}/{b['h1'][0][1]:.1f}/{b['h1'][0][2]:.1f}, x10 {b['x10'][0][0]:.1f}")
p1 = lambda v: round(v, 1)   # differences of printed values, as the manuscript takes them
dfit, doc = p1(p1(b["fit"][0][0]) - p1(b["h1"][0][0])), p1(p1(b["ocean"][0][0]) - p1(b["h1"][0][0]))
check("P8", -1.4 <= dfit <= -0.5, f"fitted h {b['fit'][0][0]:.1f} - main {b['h1'][0][0]:.1f} = {dfit:+.2f}")
check("P9", -15.3 <= doc <= -5.8, f"ocean-scale {b['ocean'][0][0]:.1f} - main = {doc:+.2f}")
check("P10", 51.7 <= round(b["thin"][0], 1) <= 73.6 and 29.4 <= round(b["lt1"][0], 1) <= 50.8
      and abs(b["hmed"][0] - 0.79) <= 0.2, f"thin {b['thin'][0]:.1f} (old 68.6), < 1 width {b['lt1'][0]:.1f} (old 45.5), "
      f"median h {b['hmed'][0]:.2f} (old 0.79)")
c = outC[ENC]
kl = ("kappa-SE", "kappa", "kappa+SE", "band lo", "band hi", "x1 (log)")
r1 = lambda t: tuple(float(round(v, 1)) for v in t)
check("P11", all(r1((c[l][0][0], c[l][1][0])) == (6.7, 0.0) for l in kl) and r1((c["x5 (log)"][0][0], c["x5 (log)"][1][0]))
      == (19.9, 0.0) and r1(c["kappa"][0]) == (6.7, 1.5, 91.8) and r1(c["kappa"][1]) == (0.0, 0.0, 100.0),
      f"step f {[r1((c[l][0][0], c[l][1][0])) for l in kl]}, x5 {r1((c['x5 (log)'][0][0], c['x5 (log)'][1][0]))}, "
      f"split {r1(c['kappa'][0])} | {r1(c['kappa'][1])}")
check("P12", round(c["column"], 1) == 95.0 and abs(c["kappa pinned"] - 84.7) <= 0.5 and abs(c["xmed"] / 3.29 - 1) <= 0.05,
      f"layered {c['column']:.1f}, pinned@k {c['kappa pinned']:.2f} (old 84.7), median x {c['xmed']:.3g} (old 3.29)")
# P13: the manuscript aggregates, recomputed over all 15 worlds from the regenerated text
q13 = []
if newA:
    i1 = newA.index(next(l for l in newA if l.startswith("P_pved change"))) + 1
    chg = [fl(l)[-8:] for l in newA[i1:i1 + 15]]
    q13 += [max(abs(v) for r in chg for v in r[:2]) <= 0.01 + 1e-9, max(abs(v) for r in chg for v in r[6:8]) == 0.08]
    i0 = newA.index(next(l for l in newA if l.startswith("Drive uniformity"))) + 1
    fr = [fl(l)[0] for l in newA[i0:i0 + 15]]
    q13 += [round(min(fr)) == 2 and round(max(fr)) == 17]
moons = [n for n in outB if n != "Early Earth ocean"]
q13 += [max(abs(p1(p1(outB[n]["fit"][0][0]) - p1(outB[n]["h1"][0][0]))) for n in outB) == 1.4,
        round(min(outB[n]["thin"][0] for n in outB)) == 52 and round(max(outB[n]["thin"][0] for n in outB)) == 74,
        all(-15.3 <= p1(p1(outB[n]["ocean"][0][0]) - p1(outB[n]["h1"][0][0])) <= -5.8 for n in moons)]
check("P13", all(q13), f"front 0.01 / 0.08 / 2-17% , fitted 1.4, thin 52-74%, ocean-scale range: {q13}")
print("\nSummary: " + ", ".join(f"{p} {'held' if ok else 'FAILED'}" for p, (ok, _) in RES.items()))
print("Remapped call sites ([s, 0] -> [s, 0, 1]), calls: " + ", ".join(f"{a} x{b}" for a, b in sorted(REMAPPED.items())))
shutil.rmtree(TMP, ignore_errors=True)
print(f"total {time.time() - T0:.0f} s")
