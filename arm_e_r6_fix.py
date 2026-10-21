"""Arm E R6, POST HOC numerical fix of arm_e_r6.seed_stable (found by the P10 independent check, arm_e_r6_analytic.log).

DEFECT. arm_e_r6.seed_stable writes exp(J tau) w0 with the factors (J11 - l2) and (J22 - l1) etc. The fast eigenvalue l2 is close to
J11, so J11 - l2 cancels to a few digits; the monomer-only seed readout (a-component) is then wrong by 1e-7 .. 4e-4 (lake) and can
change sign or vanish (ocean: seedm_a = None in arm_e_r6_results.json, true value 8e-157, mp 80 digits, scratch check). Only the seed
variants are touched (feed, case (b) and the full-seed b-readout agree with the independent code to 1e-12). arm_e_r6.py, its log and
arm_e_r6_results.json are kept unchanged as the record of what ran.

FIX. (J11 - l1)(J11 - l2) = (J22 - l1)(J22 - l2) = -J12 J21 (characteristic polynomial), so the SMALL factor of each pair is taken from
the other one, which does not cancel. Nothing else changes: this file imports arm_e_r6, replaces seed_stable, redirects the JSON write to
arm_e_r6_fix_results.json (the pre-registered arm_e_r6_results.json is not touched) and runs main() once. No prediction (P1-P12) or
definition is changed; arm_e_r6_predictions.txt is asserted by arm_e_r6.main (sha unchanged). The scored verdicts of the first run stay
the record; this run only replaces the numbers that the defect corrupted (monomer-only seed variant, retention).
Order of events: the first arm_e_r6.py run and the P10 comparison (arm_e_r6_analytic.log) came first; this file was written after seeing
the P10 mismatch; its output is arm_e_r6_fix.log.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with matplotlib python arm_e_r6_fix.py > arm_e_r6_fix.log
"""
import builtins
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arm_e_r6 as m

OLD = m.seed_stable


def seed_stable(b, tau, mono_only=False):
    """exp(J tau) w0 / (x, y) with w0 = (x, y) or (x, 0); small factors from (d1)(d2) = -J12 J21."""
    J11, J12, J21, J22 = b["J"]
    x, y, l1, l2 = b["x"], b["y"], b["l1"], b["l2"]
    w = 0.0 if mono_only else y
    p = J12 * J21

    def pair(Jii):
        d1, d2 = Jii - l1, Jii - l2
        return np.where(np.abs(d1) < np.abs(d2), -p / d2, d1), np.where(np.abs(d2) < np.abs(d1), -p / d1, d2)

    with np.errstate(all="ignore"):
        d11, d12 = pair(J11)
        d21, d22 = pair(J22)
        e1, e2, dl = np.exp(l1 * tau), np.exp(l2 * tau), l1 - l2
        a = (e1 * (d12 * x + J12 * w) - e2 * (d11 * x + J12 * w)) / dl
        bb = (e1 * (J21 * x + d22 * w) - e2 * (J21 * x + d21 * w)) / dl
    return a / x, bb / y


def check():
    """New seed_stable vs a per-point 2x2 eig route at tq = 1.5/|l1| and at tq = 8/|l1| (grid of 5 mM, s = 6); report the old one too."""
    b = m.base(5e-3, 1.0, 6.0)
    J = [np.asarray(v) for v in b["J"]]
    worst = {}
    for mult in (1.5, 8.0):
        tq = mult / np.abs(b["l1"])
        for mono in (False, True):
            ref_a, ref_b = np.zeros(b["x"].shape), np.zeros(b["x"].shape)
            for i in np.ndindex(b["x"].shape):
                M = np.array([[J[0][i], J[1][i]], [J[2][i], J[3][i]]])
                w, V = np.linalg.eig(M)
                w0 = np.array([b["x"][i], 0.0 if mono else b["y"][i]])
                q = (V @ (np.exp(w.real * tq[i]) * np.linalg.solve(V, w0))).real
                ref_a[i], ref_b[i] = q[0] / b["x"][i], q[1] / b["y"][i]
            for nm, f_ in (("new", seed_stable), ("old", OLD)):
                qa, qb = f_(b, tq, mono)
                for ch, q, ref in (("a", qa, ref_a), ("b", qb, ref_b)):
                    ok = np.isfinite(ref) & (np.abs(ref) > 1e-290)
                    dev = np.abs(q - ref)[ok] / np.abs(ref)[ok]
                    worst[(nm, mult, "mono" if mono else "full", ch)] = float(np.nanmax(dev))
    for k, v in worst.items():
        print(f"  seed check {k[0]} t={k[1]}/|l1| {k[2]:4} {k[3]}: max rel dev vs eig {v:.2e}")
    assert max(v for k, v in worst.items() if k[0] == "new") < 1e-6, worst
    assert max(v for k, v in worst.items() if k[0] == "old") > 1e-6, "old routine did not show the defect on this grid (check is not biting)"


def compare():
    """Which numbers of arm_e_r6_results.json (first run) differ from arm_e_r6_fix_results.json (this run)."""
    A = json.load(open(HERE / "arm_e_r6_results.json"))["records"]
    B = json.load(open(HERE / "arm_e_r6_fix_results.json"))["records"]
    print("\nfirst run (arm_e_r6_results.json) vs fixed run (arm_e_r6_fix_results.json), per record/setting: max rel change of e_min and of comp.*")
    for key in A:
        for i, (ra, rb) in enumerate(zip(A[key], B[key])):
            if ra is None or rb is None:
                print(f"  {key:16} {('pond', 'lake', 'ocean')[i]:5} none/none" if ra is rb else f"  {key:16} {i} one side none")
                continue
            ca, cb = {"e_min": ra["e_min"], **ra["comp"]}, {"e_min": rb["e_min"], **rb["comp"]}
            dev, who = 0.0, ""
            for k in ca:
                u, v = ca[k], cb[k]
                if u is None or v is None or not (np.isfinite(u) and np.isfinite(v)):
                    if (u is None) != (v is None) and not (u in (None, float("inf")) and v in (None, float("inf"))):
                        dev, who = float("inf"), f"{k}: {u} -> {v}"
                    continue
                if abs(u - v) / max(abs(u), 1e-300) > dev:
                    dev, who = abs(u - v) / max(abs(u), 1e-300), k
            same_pt = all(ra[k] == rb[k] for k in ("f", "kh", "R", "K"))
            print(f"  {key:16} {('pond', 'lake', 'ocean')[i]:5} same grid point {same_pt!s:5} max rel change {dev:.2e} {who}")


if __name__ == "__main__":
    print("seed_stable check (fix vs per-point eig; old routine shown for comparison)")
    check()
    m.seed_stable = seed_stable
    m.open = lambda p, *a, **k: builtins.open(HERE / "arm_e_r6_fix_results.json" if Path(p).name == "arm_e_r6_results.json" else p, *a, **k)
    m.main()
    compare()
