"""h14d_advect_posthoc: POST HOC, exploratory, NOT pre-registered. Reads h14d_advect_raw.pkl (and, through h14d_advect, h14c_advect_raw.pkl).
 (a) per-doubling increments of S50 (a log law gives a constant increment per doubling, a power law a growing one) and of kicks to clear;
 (b) S50 = a + b ln(ratio) over the five ratios against the power law: residuals, b with an interval, and the value at ratio 1e9 (NOT a measurement);
 (c) job-level (cluster) bootstrap of q(64 -> 128), as a check on the run-level one (runs in a job share a seed stream; the flow phase is drawn per run).
Run: UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14d_advect_posthoc.py > h14d_advect_posthoc.log
"""
import pickle

import numpy as np

import h14b_advect as H
import h14c_advect as C
import h14d_advect as D          # asserts the predictions-file and imported-file sha256 on import

d = pickle.load(open(D.PKL, "rb"))
res = d["res"]
P, ST = D.stats(res)
KEYS = ("R8_10", "R16_10", "V32_10", "V64_10", "V128_10")
R = np.array([D.RATIO[k] for k in KEYS], float)
LR = np.log(R)
pt = np.array([ST[k].pt[0] for k in KEYS])
bs = np.stack([ST[k].bs[:, 0] for k in KEYS], 1)
ci = lambda x: np.percentile(x[np.isfinite(x)], [2.5, 97.5])
f2 = lambda x: f"{x:.2f}"

print("(a) PER-DOUBLING INCREMENTS (S50 in t17; kicks to clear = S50 t17 / Tk), 95% bootstrap intervals")
kick = np.array([ST[k].t17 / P[k]["spec"]["flow"]["Tk"] for k in KEYS])
for i in range(4):
    inc = bs[:, i + 1] - bs[:, i]
    kc = (bs[:, i + 1] - bs[:, i]) * kick[i + 1]
    print(f"  1/{int(R[i])} -> 1/{int(R[i + 1])}: S50 {f2(pt[i])} -> {f2(pt[i + 1])}, increment {f2(pt[i + 1] - pt[i])} [{f2(ci(inc)[0])}..{f2(ci(inc)[1])}]; "
          f"kicks to clear {pt[i] * kick[i]:.0f} -> {pt[i + 1] * kick[i + 1]:.0f}, increment {(pt[i + 1] * kick[i + 1] - pt[i] * kick[i]):.0f}")
inc_pt = np.diff(pt)
print(f"  increments {', '.join(f2(x) for x in inc_pt)}; last / mean of the first three = {inc_pt[3] / inc_pt[:3].mean():.2f}")
dq = (bs[:, 4] - bs[:, 3]) - (bs[:, 3] - bs[:, 2])
print(f"  increment(64 -> 128) - increment(32 -> 64) = {f2(inc_pt[3] - inc_pt[2])} [{f2(ci(dq)[0])}..{f2(ci(dq)[1])}]  (a power law with q in 0.3-0.4 would make it positive)")

print("\n(b) LOG LAW S50 = a + b ln(ratio) against POWER LAW ln S50 = a + q ln(ratio), both least squares over the five ratios")
def lin(v):
    return np.polyfit(LR, v, 1)
bl, al = lin(pt)
bq, aq = lin(np.log(pt))
res_l = pt - (al + bl * LR)
res_q = np.log(pt) / 1.0 - (aq + bq * LR)
bl_b = np.array([lin(v)[0] for v in bs])
al_b = np.array([lin(v)[1] for v in bs])
print(f"  log law:   b = {f2(bl)} per e-fold [95% {f2(ci(bl_b)[0])}..{f2(ci(bl_b)[1])}], a = {f2(al)}; residuals in S50: " + ", ".join(f"{x:+.2f}" for x in res_l) + f"; rms {np.sqrt((res_l ** 2).mean()):.2f}")
print(f"  power law: q = {f2(bq)}; residuals in S50: " + ", ".join(f"{x:+.2f}" for x in pt - np.exp(aq + bq * LR)) + f"; rms {np.sqrt(((pt - np.exp(aq + bq * LR)) ** 2).mean()):.2f}")
se = np.array([(np.percentile(bs[:, i], 97.5) - np.percentile(bs[:, i], 2.5)) / 3.92 for i in range(5)])
chi_l = (((pt - (al + bl * LR)) / se) ** 2).sum()
chi_q = (((pt - np.exp(aq + bq * LR)) / se) ** 2).sum()
print(f"  chi-square of the five S50 points (SE from the 95% bootstrap width / 3.92), 3 d.o.f.: log law {chi_l:.1f}, power law {chi_q:.1f} (both unweighted fits, so only indicative)")
v9 = al + bl * np.log(1e9)
v9_b = al_b + bl_b * np.log(1e9)
print(f"  log law at ratio 1e9 (NOT a measurement; 4 doublings measured, 23 more extrapolated): S50 = {v9:.0f} [95% {ci(v9_b)[0]:.0f}..{ci(v9_b)[1]:.0f}]; "
      f"power-law fit over all five (q {f2(bq)}): {np.exp(aq + bq * np.log(1e9)):.2g}")
print("  ratio at which the log law reaches f = 100, 200, 500: " + ", ".join(f"{f}: {np.exp((f - al) / bl):.2e}" for f in (100, 200, 500)))

print("\n(c) JOB-LEVEL BOOTSTRAP of q(64 -> 128) (resample whole jobs; 1/64 = 18 jobs of 21 runs, 1/128 = 16 jobs)")
cs = D.cell_sets(res)


def job_ids(parts):
    return np.concatenate([np.full(p["runs"], j) for j, p in enumerate(parts)])


def cluster(key, B=2000, seed=0):
    parts, _ = cs[key]
    c = P[key]
    s = C.Stat(c, B=0)
    jid = job_ids(parts)
    rng = np.random.default_rng(seed)
    nj = len(parts)
    out = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, nj, nj)
        idx = np.concatenate([np.flatnonzero(jid == j) for j in pick])
        out[b] = s._est(idx)[0]
    per = np.array([C.Stat(C.pool([p]), B=0).pt[0] for p in parts])
    return out, per
c64, per64 = cluster("V64_10", seed=64)
c128, per128 = cluster("V128_10", seed=128)
q_c = np.log(c128 / c64) / np.log(2)
q_pt = np.log(ST["V128_10"].pt[0] / ST["V64_10"].pt[0]) / np.log(2)
q_run = np.log(ST["V128_10"].bs[:, 0] / ST["V64_10"].bs[:, 0]) / np.log(2)
print(f"  q(64 -> 128) point {f2(q_pt)}; run-level 95% {f2(ci(q_run)[0])}..{f2(ci(q_run)[1])} (the scored one); job-level 95% {f2(ci(q_c)[0])}..{f2(ci(q_c)[1])}")
print(f"  S50 per job: 1/64 mean {per64.mean():.2f} sd {per64.std(ddof=1):.2f} (18 jobs, 21 runs each, so SE of a job-mean {per64.std(ddof=1) / np.sqrt(len(per64)):.2f}); "
      f"1/128 mean {per128.mean():.2f} sd {per128.std(ddof=1):.2f} (16 jobs, SE {per128.std(ddof=1) / np.sqrt(len(per128)):.2f}); "
      f"run-level bootstrap SE {np.std(ST['V64_10'].bs[:, 0]):.2f} and {np.std(ST['V128_10'].bs[:, 0]):.2f}")
