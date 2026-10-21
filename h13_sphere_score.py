"""h13_sphere_score: reads the h13_sphere_*.log files and applies the pass/fail rules of h13_sphere_predictions.txt (P1-P8, D0-D3).
Only arithmetic on logged lines; a missing log prints "NOT RUN" for the predictions that need it.
Run: python3 h13_sphere_score.py > h13_sphere_score.log
"""
import math
import os
import re

import numpy as np

num = r"([-+0-9.]+|inf|nan)"
tau_R2 = 3.065  # R^2/D in units of tau_sw (h13_sphere_s16_f50.log:4 and the sphere definition tau_sw = (pi R/5.5)^2/D)


def log(name):
    return open(f"h13_sphere_{name}.log").read() if os.path.exists(f"h13_sphere_{name}.log") else None


def summ(name):
    t = log(name)
    m = t and re.search(rf"reached one hand (\d+)/(\d+).*?median {num} p90 {num} p99 {num}", t)
    if not m:
        return None
    w = re.search(r"wins ([0-9.]+) \((\d+)/(\d+)\)", t)
    e = re.search(r"survival e-folds in ([0-9.]+) R", t)
    return dict(ok=int(m[1]), n=int(m[2]), med=float(m[3]), p90=float(m[4]), p99=float(m[5]), win=(int(w[2]), int(w[3])) if w else None,
                tail=float(e[1]) if e else None, l={int(a): float(b) for a, b in re.findall(r"l\((\d+)\) = ([0-9.]+)", t)})


def show(p, ok, text):
    print(f"[{p}] {'HELD' if ok else 'FAILED' if ok is not None else 'NOT RUN'}: {text}")


assert abs(float(re.search(num, "x 1.5 y")[0]) - 1.5) < 1e-12  # parser self-check
S = {k: summ(k) for k in ("s16_f50", "s16_f55", "s32_f50", "s32_f55", "s48_f50", "t32_f50", "t32_f55", "t57_f50")}
sph = [S[k] for k in ("s16_f50", "s16_f55", "s32_f50", "s32_f55", "s48_f50") if S[k]]
bad, tot = sum(s["n"] - s["ok"] for s in sph), sum(s["n"] for s in sph)
show("P1", bad <= 0.01 * tot, f"{bad} of {tot} sphere runs not at one hand at 30 tau_sw (" + ", ".join(f"{k} {S[k]['ok']}/{S[k]['n']}" for k in S if k[0] == "s" and S[k]) + ")")
bands = {"t32_f50": (0.25, 0.55), "t32_f55": (0.10, 0.50), "t57_f50": (0.20, 0.55)}
fr = {k: 1 - S[k]["ok"] / S[k]["n"] for k in bands if S[k]}
show("P2", all(b[0] <= fr[k] <= b[1] for k, b in bands.items()) if len(fr) == 3 else None,
     ", ".join(f"{k} stuck {v:.2f} (band {bands[k]})" for k, v in fr.items()))
t = log("t512")
if t:
    m = re.search(rf"exponent: l ~ t\^{num} .*?fixed 0.5: {num}", t)
    l512 = {int(a): float(b) for a, b in re.findall(r"l\((\d+)\) = ([0-9.]+)", t)}
    e, A = float(m[1]), float(m[2])
    d = {tt: S["s32_f50"]["l"][tt] / l512[tt] - 1 for tt in (20, 40)} if S["s32_f50"] else {}
    show("P3", 0.45 <= e <= 0.55 and 4.5 <= A <= 8.5 and bool(d) and all(abs(v) <= 0.15 for v in d.values()),
         f"exponent {e} (band 0.45-0.55), A {A} (band 4.5-8.5), sphere R=32 vs torus 512 l(20), l(40) differ by " + ", ".join(f"{v:+.1%}" for v in d.values()))
else:
    show("P3", None, "t512 log missing")
a, b = S["s16_f50"], S["s16_f55"]
if a and b:
    show("P4", 2 <= a["med"] <= 6 and 5 <= a["p90"] <= 14 and 9 <= a["p99"] <= 24 and a["med"] > 1.5 and 1.5 <= b["med"] <= 5 and b["med"] < a["med"],
         f"S16a median {a['med']} (2-6) p90 {a['p90']} (5-14) p99 {a['p99']} (9-24); f = 1 needs median <= 1.5: {a['med'] <= 1.5}; S16b median {b['med']} (1.5-5, below S16a)")
    r32 = S["s32_f50"]["med"] / a["med"] if S["s32_f50"] else None
    r48 = S["s48_f50"]["med"] / a["med"] if S["s48_f50"] else None
    show("P5", None if r32 is None else (0.8 <= r32 <= 1.25 and (r48 is None or 0.7 <= r48 <= 1.4) and 0.7 <= a["tail"] <= 1.4),
         f"median ratio S32a/S16a {r32}, S48a/S16a {r48} (S48a part {'run' if r48 else 'NOT RUN'}); S16a tail e-fold {a['tail']} R^2/D (0.7-1.4)")
    if S["s32_f50"] and S["s32_f55"]:
        w0 = [S[k]["win"] for k in ("s16_f50", "s32_f50")]
        w5 = [S[k]["win"] for k in ("s16_f55", "s32_f55")]
        s0, s5 = sum(x[0] for x in w0) / sum(x[1] for x in w0), sum(x[0] for x in w5) / sum(x[1] for x in w5)
        show("P6", 0.4 <= s0 <= 0.6 and 0.8 <= s5 <= 1.0, f"+ hand wins, f0 0.50 pooled {s0:.3f} (0.40-0.60); f0 0.55 pooled {s5:.3f} (0.80-1.00)")
    # decision rule on the S16a numbers (the pre-registered reference run)
    d1 = a["med"] <= 1.5 and a["p90"] <= 2.5
    print(f"D0: share f of Eq. 12 stays 1 if P1 holds. D1: f = 1 stays only if median <= 1.5 and p90 <= 2.5: {d1}. "
          f"D2: f_t = ceil(p90) = {math.ceil(a['p90'])} (median {a['med']}, p99 {a['p99']}); A = 5.5/sqrt(f_t) = {5.5 / math.sqrt(math.ceil(a['p90'])):.2f}. "
          f"D3: p99 {a['p99']} {'<=' if a['p99'] <= 10 else '>'} 10 -> " + ("existing x10 bound covers it" if a["p99"] <= 10 else "rerun coarsen_slow at f = p90 and p99"))
sp = log("special")
if sp:
    ratios = [float(x) for x in re.findall(r"ratio ([0-9.]+)", sp)]
    rates = [float(x) for x in re.findall(r"decay rate ([0-9.]+)", sp)]
    mino = [float(x) for x in re.findall(r"= ([0-9.]+); T1", sp)]
    show("P7", len(ratios) == 4 and all(abs(r - 1) <= 0.10 for r in ratios) and all(2.5 <= r <= 3.5 for r in rates) and all(m >= 0.45 for m in mino),
         f"cap-law ratios {ratios} (within 10%), m=2 decay rates {rates} (2.5-3.5), minority at 30 tau_sw {mino} (>= 0.45)")
sh = log("shell")
if sh:
    T = [float(x) for x in re.findall(r"T1/tau_sw=" + num, sh)]
    n_ok = sum(math.isfinite(x) for x in T)
    med = float(np.median(T)) if n_ok == len(T) else float("nan")
    ref = S["s16_f50"]["med"] if S["s16_f50"] else float("nan")
    show("P8", len(T) == 12 and n_ok == 12 and 0.6 <= med / ref <= 1.7, f"{n_ok}/{len(T)} reached one hand; shell median {med:.2f} tau_sw, ratio to S16a {med / ref:.2f} (0.6-1.7)")
else:
    show("P8", None, "shell log missing")
