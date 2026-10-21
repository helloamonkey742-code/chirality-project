"""dh_range.py: per-world headline shares (multiwall P_pved) with the icy-moon D_h upper bound at 1e2 (current), 3e2, 1e3.

Why: uncertainty.py:38-39 says Ashkenazy & Tziperman 2021 "use 30". The paper text has no such value
(papers/h14_eddy/ashkenazy2021.txt:402-405: explicit tracer diffusion 5 (2D) and 2 (3D) m2/s; :553-559: diagnosed about
300 m2/s, rough lower bound about 200), so the bracket 1e-2..1e2 sits below the diagnosed value. This reruns
multiwall.run4 (the Section 4 path that prints the paper's shares) with the D_h range of Enceladus and Europa
(U.WORLDS) and of the 12 other moons (ocean_worlds.DH) overridden. Nothing existing is edited. Same seed, same uniform
stream (common random numbers): only D_h moves. Early Earth keeps 1e2..1e4.
Enceladus is run as logged (seed collision, reproduces multiwall.log) and with seedfix.py's wrapper (the paper's number).

Predictions: dh_range_predictions.txt, sha256 below, hashed 2026-10-21 00:59:38 PDT before any run. (An earlier hash of
the same file, b203cdbe..., taken 00:59:04, was superseded by a log-line-number fix; no run happened in between.)
Smoke (disclosed): DH_SMOKE=1 (N = 2000, no log comparison, no scoring). Two smoke runs on 2026-10-21 before the full run:
the first stopped at "No module named matplotlib" (uncertainty.py imports ocean.py, which imports matplotlib), the second ran
to the end in 4 s. They were used only to debug code paths and the output format; no smoke number is used anywhere.
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba --with matplotlib python dh_range.py > dh_range.log
"""
import hashlib
import itertools
import os
import re
import time

import numpy as np

PRED_SHA256 = "5b75060fbfaca2f1efb64b7e5be7132050dfe5e538409178a44d97d8cfe591e2"
HERE = os.path.dirname(os.path.abspath(__file__))
SMOKE = os.environ.get("DH_SMOKE") == "1"
assert hashlib.sha256(open(os.path.join(HERE, "dh_range_predictions.txt"), "rb").read()).hexdigest() == PRED_SHA256, \
    "dh_range_predictions.txt changed after it was hashed"

def same(a, b):
    """Exact deep equality of run4 outputs (repr would truncate the cnt arrays)."""
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return bool(np.array_equal(a, b))


ORIG_RNG = np.random.default_rng      # seedfix.py's wrapper, switched on only for the 'fix' runs
NFIX = [0, 0]                         # [remapped calls, fixed run4 calls]


def fixed_rng(seed=None, *a, **kw):
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        NFIX[0] += 1
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


import multiwall as M  # noqa: E402  (imports sweep2, sweep3, uncertainty)
import ocean_worlds as OW  # noqa: E402
import uncertainty as U  # noqa: E402

LO, HIS = 1e-2, (1e2, 3e2, 1e3)
BASE, MID, TOP = HIS
EARTH = "Early Earth ocean"
W0, DH0 = dict(U.WORLDS), OW.DH
if SMOKE:
    U.N = 2000


def setrange(hi):
    for n in ("Enceladus", "Europa"):
        v, L, c, _, dz, t = W0[n]
        U.WORLDS[n] = (v, L, c, (LO, hi), dz, t)
    OW.DH = (LO, hi)


def restore():
    U.WORLDS.update(W0)
    OW.DH = DH0


def run(seed, hi, fixed, full):
    setrange(hi)
    assert OW.args(1e6, 1e4, (1e8, 1e9))[3] == (LO, hi) and U.WORLDS["Europa"][3] == (LO, hi)
    np.random.default_rng = fixed_rng if fixed else ORIG_RNG
    NFIX[1] += fixed
    try:
        return M.run4(seed, full=full)
    finally:
        np.random.default_rng = ORIG_RNG
        restore()


def crn_check():
    """Only D_h differs between hi = 1e2 and 1e3 in the raw draws, and log10(D_h)+2 scales by (4+d)/4 exactly."""
    for n in ("Enceladus", "Europa"):
        xs = {}
        for hi in (BASE, TOP):
            setrange(hi)
            U.RNG = ORIG_RNG(0)
            xs[hi] = U.sweep(*U.WORLDS[n])[0]
            d = xs[hi]["horizontal D_h"]
            assert d.min() >= LO * (1 - 1e-9) and d.max() <= hi * (1 + 1e-9), (n, hi, d.min(), d.max())
            restore()
        for k in xs[BASE]:
            if k != "horizontal D_h":
                assert np.array_equal(xs[BASE][k], xs[TOP][k]), (n, k)
        r = (np.log10(xs[TOP]["horizontal D_h"]) + 2) / (np.log10(xs[BASE]["horizontal D_h"]) + 2)
        assert np.abs(r - 5 / 4).max() < 1e-6, (n, np.abs(r - 5 / 4).max())
    return "CRN check ok: Enceladus and Europa raw draws identical except D_h; log10(D_h)+2 scales by 5/4 at 1e3"


t0 = time.time()
print(f"dh_range.py, N = {U.N}, D_h upper bound {HIS} for Enceladus, Europa and the 12 other moons; Early Earth fixed "
      f"at 1e2..1e4; predictions sha256 {PRED_SHA256[:12]}...")
print(crn_check())

R = {}
for hi in HIS:
    R[(hi, "col")] = run(0, hi, False, True)
    R[(hi, "fix")] = run(0, hi, True, True)
SEEDS = {(s, hi): run(s, hi, True, False) for s in (1, 2, 3) for hi in HIS}
assert NFIX[0] == NFIX[1], NFIX                  # the wrapper remapped [s, 0] exactly once per fixed run4 call
ALL = list(R[(BASE, "fix")])
ICY = [n for n in ALL if n != EARTH]
assert len(ALL) == 15 and ALL[0] == "Enceladus"
for hi in HIS:                                   # the wrapper moves Enceladus only
    for n in ALL[1:]:
        assert same(R[(hi, "col")][n], R[(hi, "fix")][n]), (hi, n)

# ---- 0. self-check: hi = 1e2 reproduces the committed logs to the printed 0.1
nums = lambda line: re.findall(r"\d+\.\d", line)
THREE = lambda r: [x for k in ("none", "main") for lid in (0, 1) for x in r[k][lid]]
PVED = lambda r: [r[k][lid][0] for k in ("none", "bound", "main", "x10", "sum") for lid in (0, 1)]


def log_rows(path, first, last, head_line, head_text):
    lines = open(os.path.join(HERE, path)).read().split("\n")
    assert lines[head_line - 1].startswith(head_text), (path, head_line, lines[head_line - 1])
    return {lines[i - 1][:19].strip(): (i, lines[i - 1]) for i in range(first, last + 1)}


print("\n=== 0. Self-check: upper bound 1e2 reproduces the committed logs to the printed 0.1 ===")
if SMOKE:
    print("(skipped in smoke mode)")
else:
    three = log_rows("multiwall.log", 525, 539, 524, "world")
    pved = log_rows("multiwall.log", 542, 556, 541, "P_pved %, no lid / final")
    assert sorted(three) == sorted(ALL) == sorted(pved)
    for n in ALL:
        for rows, fn, what in ((three, THREE, "three-way"), (pved, PVED, "P_pved")):
            got = [f"{x:.1f}" for x in fn(R[(BASE, "col")][n])]
            assert got == nums(rows[n][1]), (n, what, got, nums(rows[n][1]))
    print("Enceladus (collided) and the other 14 worlds: three-way split EXACT vs multiwall.log:525-539; "
          "P_pved no lid / final (sweep3, bound, multiwall, cascade x10, sum rule) EXACT vs multiwall.log:542-556")
    sf = open(os.path.join(HERE, "seedfix.log")).read().split("\n")
    for ln, fn, what in ((122, THREE, "three-way"), (125, PVED, "P_pved")):
        got = [f"{x:.1f}" for x in fn(R[(BASE, "fix")]["Enceladus"])]
        assert "Enceladus" in sf[ln - 1] and got == nums(sf[ln - 1]), (ln, got, sf[ln - 1])
    print("Enceladus (seedfix wrapper): three-way split EXACT vs seedfix.log:122; P_pved row EXACT vs seedfix.log:125")
print("Worlds 1-14 identical between the collided and the fixed run at every bound (exact array equality); "
      f"wrapper remapped {NFIX[0]} calls in {NFIX[1]} fixed runs")

# ---- accessors
G = lambda hi, n, key="main", lid=0, j=0: R[(hi, "fix")][n][key][lid][j]
M0, MF = (lambda hi, n: G(hi, n)), (lambda hi, n: G(hi, n, lid=1))
X0, XF = (lambda hi, n: G(hi, n, "x10")), (lambda hi, n: G(hi, n, "x10", 1))
C0, N0 = (lambda hi, n: G(hi, n, j=1)), (lambda hi, n: G(hi, n, j=2))
CH = lambda f, hi, n: f(hi, n) - f(BASE, n)
ORDER = sorted(ALL, key=lambda n: -MF(BASE, n))


def show(title, f0, f1):
    print(f"\n{title}")
    print(f"{'world':18} " + " | ".join(f"{'bound ' + format(h, '.0e'):>14}" for h in HIS)
          + "  ||  change at 3e2 | change at 1e3     (no lid / final)")
    for n in ORDER:
        lv = " | ".join(f"{f0(h, n):6.1f} /{f1(h, n):6.1f}" for h in HIS)
        ch = " | ".join(f"{CH(f0, h, n):+6.2f} /{CH(f1, h, n):+6.2f}" for h in (MID, TOP))
        print(f"{n:18} {lv}  ||  {ch}")


print("\n=== 1. Headline: multiwall P_pved %, seed 0, Enceladus with the seed fix (the paper's number) ===")
show("1a. main rule", M0, MF)
show("1b. cascade x10 (low-end timing)", X0, XF)
print("\n1c. Enceladus as logged (seed collision), main and x10, no lid / final")
for hi in HIS:
    r = R[(hi, "col")]["Enceladus"]
    print(f"  bound {hi:.0e}: main {r['main'][0][0]:5.1f} / {r['main'][1][0]:5.1f}   x10 {r['x10'][0][0]:5.1f} / "
          f"{r['x10'][1][0]:5.1f}   (fixed: {M0(hi, 'Enceladus'):5.1f} / {MF(hi, 'Enceladus'):5.1f}   "
          f"{X0(hi, 'Enceladus'):5.1f} / {XF(hi, 'Enceladus'):5.1f})")

print("\n=== 2. Where the mass goes (main rule, no lid): P_pved / P_chance / P_none % and their change at 1e3 ===")
print(f"{'world':18} " + " | ".join(f"{'bound ' + format(h, '.0e'):>22}" for h in HIS)
      + "  ||  change at 1e3: P_pved  P_chance  P_none  one-hand")
for n in ORDER:
    lv = " | ".join(f"{M0(h, n):5.1f}/{C0(h, n):4.1f}/{N0(h, n):5.1f}" for h in HIS)
    print(f"{n:18} {lv}  ||  {CH(M0, TOP, n):+14.2f} {CH(C0, TOP, n):+9.2f} {CH(N0, TOP, n):+8.2f} "
          f"{CH(M0, TOP, n) + CH(C0, TOP, n):+9.2f}")

print("\n=== 3. Order of the 15 worlds by final P_pved (main) ===")
for hi in HIS:
    print(f"  bound {hi:.0e}: " + " > ".join(f"{n} {MF(hi, n):.1f}" for n in sorted(ALL, key=lambda n: -MF(hi, n))))
FLIPS = [(a, b, hi, MF(BASE, a) - MF(BASE, b)) for hi in (MID, TOP) for a, b in itertools.permutations(ALL, 2)
         if MF(BASE, a) > MF(BASE, b) and MF(hi, a) < MF(hi, b)]
print("Adjacent gaps at 1e2: " + ", ".join(f"{a}/{b} {MF(BASE, a) - MF(BASE, b):.2f}"
                                         for a, b in zip(ORDER, ORDER[1:]) if MF(BASE, a) - MF(BASE, b) < 2))
print("Pairs whose order flips: " + ("none" if not FLIPS else "; ".join(
    f"{a} falls below {b} at {hi:.0e} (1e2 gap {g:.2f})" for a, b, hi, g in FLIPS)))

print("\n=== 4. Seeds 1-3 (Enceladus fixed): change in main P_pved no lid, 1e2 -> 1e3, points ===")
SCH = {s: {n: SEEDS[(s, TOP)][n]["main"][0][0] - SEEDS[(s, BASE)][n]["main"][0][0] for n in ALL} for s in (1, 2, 3)}
SCH[0] = {n: CH(M0, TOP, n) for n in ALL}
print(f"{'world':18} " + " ".join(f"{'seed ' + str(s):>8}" for s in (0, 1, 2, 3)))
for n in ORDER:
    print(f"{n:18} " + " ".join(f"{SCH[s][n]:+8.2f}" for s in (0, 1, 2, 3)))
print("max |seed k - seed 0| over icy worlds: " + ", ".join(
    f"seed {s} {max(abs(SCH[s][n] - SCH[0][n]) for n in ICY):.2f}" for s in (1, 2, 3)))

# ---- scoring against dh_range_predictions.txt (rules fixed there; no tuning after the run)
print("\n=== 5. Predictions (dh_range_predictions.txt) ===")
SC = []


def verdict(tag, ok, text):
    SC.append(ok)
    print(f"{tag} {'HELD' if ok else 'FAILED'}: {text}")


if SMOKE:
    print("(skipped in smoke mode)")
else:
    verdict("P1", True, "self-check above passed (the script stops with an AssertionError otherwise)")
    verdict("P2", all(same(R[(hi, v)][EARTH], R[(BASE, v)][EARTH]) for hi in HIS for v in ("col", "fix")),
            "Early Earth identical at 3e2 and 1e3 (every field, both variants)")
    down = [n for n in ICY if CH(M0, TOP, n) <= -0.1]
    mono = [n for n in ICY if M0(TOP, n) <= M0(MID, n) <= M0(BASE, n)]
    verdict("P3", len(down) >= 13 and len(mono) >= 13,
            f"M0 falls by >= 0.1 at 1e3 in {len(down)}/14 (need 13); monotone in {len(mono)}/14 (need 13); "
            f"not falling: {[n for n in ICY if n not in down]}")
    bad = [(f.__name__, n, hi, round(CH(f, hi, n), 2)) for f in (M0, X0) for n in ICY
           for hi, lo in ((MID, -1.0), (TOP, -2.0)) if not lo <= CH(f, hi, n) <= 0.0]
    verdict("P4a", not bad, f"all M0 and X0 changes within [-1.0, 0] at 3e2 and [-2.0, 0] at 1e3; outside: {bad}")
    cen = {"Enceladus": (-0.4, -0.9), "Europa": (-0.2, -0.5)}
    bad = [(n, round(CH(M0, MID, n), 2), round(CH(M0, TOP, n), 2)) for n, (a, b) in cen.items()
           if abs(CH(M0, MID, n) - a) > 0.4 or abs(CH(M0, TOP, n) - b) > 0.6]
    verdict("P4b", not bad, "Enceladus (fixed) and Europa M0 changes " + ", ".join(
        f"{n} {CH(M0, MID, n):+.2f}/{CH(M0, TOP, n):+.2f} (central {a:+.1f}/{b:+.1f})" for n, (a, b) in cen.items()))
    bad = [(n, hi) for n in ("Europa", "Dione", "Miranda", "Mimas") for hi in (MID, TOP)
           if not (CH(M0, hi, n) < 0 and CH(MF, hi, n) <= 0 and abs(CH(MF, hi, n)) <= abs(CH(M0, hi, n)))]
    verdict("P5", all(MF(hi, "Enceladus") == 0.0 for hi in HIS) and not bad,
            f"Enceladus final {[MF(hi, 'Enceladus') for hi in HIS]}; lid worlds with MF change not (same sign, no larger): {bad}")
    BLK = [[EARTH], ["Umbriel", "Oberon", "Titania", "Ariel"], ["Ganymede", "Callisto", "Pluto"], ["Triton", "Titan"],
           ["Europa"], ["Dione"], ["Miranda"], ["Mimas"], ["Enceladus"]]
    assert sorted(sum(BLK, [])) == sorted(ALL)
    brk = [(hi, i) for hi in HIS for i in range(len(BLK) - 1)
           if min(MF(hi, n) for n in BLK[i]) <= max(MF(hi, n) for n in BLK[i + 1])]
    verdict("P6", all(g <= 0.7 for *_, g in FLIPS) and not brk,
            "flips: " + ("none" if not FLIPS else "; ".join(f"{a}/{b} at {hi:.0e}, 1e2 gap {g:.2f}"
                                                         for a, b, hi, g in FLIPS))
            + f"; block order broken at (bound, block index): {brk if brk else 'none'}")
    rise = [n for n in ICY if CH(C0, TOP, n) >= 0.05]
    one = {n: CH(M0, TOP, n) + CH(C0, TOP, n) for n in ICY}
    verdict("P7", len(rise) >= 12 and all(abs(v) <= 0.5 for v in one.values()),
            f"P_chance rises >= 0.05 in {len(rise)}/14 (need 12); largest |one-hand change| "
            f"{max(abs(v) for v in one.values()):.2f} ({max(one, key=lambda n: abs(one[n]))}; limit 0.5)")
    q = [n for n in ICY if CH(M0, MID, n) <= -0.3]
    rat = {n: CH(M0, TOP, n) / CH(M0, MID, n) for n in q}
    inb = [n for n in q if 1.7 <= rat[n] <= 2.6]
    verdict("P8", len(q) >= 3 and len(inb) >= 0.8 * len(q),
            f"{len(inb)}/{len(q)} qualifying worlds have a 1e3/3e2 ratio in 1.7-2.6; ratios min {min(rat.values()):.2f} "
            f"max {max(rat.values()):.2f}" if q else "no qualifying world")
    dev = {s: max(abs(SCH[s][n] - SCH[0][n]) for n in ICY) for s in (1, 2, 3)}
    verdict("P9", all(v <= 0.3 for v in dev.values()), "max |seed k - seed 0| of the 1e3 change: "
            + ", ".join(f"seed {s} {v:.2f}" for s, v in dev.items()))
    print(f"Held {sum(SC)} of {len(SC)}.")
print(f"\nwall time {time.time() - t0:.0f} s")
