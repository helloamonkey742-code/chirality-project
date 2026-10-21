"""h14c_advect_posthoc: POST HOC, exploratory, NOT pre-registered. Reads h14c_advect_raw.pkl (written by h14c_advect.py) and asks
what the pre-registered report printed only as a point or not at all:
 (a) the length exponent with the pre-freeze-out time matched: V32_3.6m (gamma 0.0556, 1.06 t17, 5 kicks before freeze-out) against
     V32_10 and V32_30 (same 1.07 t17 and 5 kicks), with a bootstrap, and p(10 -> 30) with an interval (P8 used a point);
 (b) the 8 V32_30 runs one by one;
 (c) the ratio D_eff/D_mol at which S50(l/w 10) reaches each distinct break-even f of h14_eddy.log, and S50(l/w 10) at the real
     ratio ~1e9 IF the power law q held over 7 decades (an extrapolation, not a measurement), against the molecular ceiling.
Run: UV_CACHE_DIR=$HOME/Library/Caches/uv nice -n 10 uv run --no-project --with numpy --with scipy --with numba --with matplotlib python h14c_advect_posthoc.py > h14c_advect_posthoc.log
"""
import pickle

import numpy as np

import h14b_advect as H
import h14c_advect as C          # asserts the predictions-file sha256 on import

d = pickle.load(open(C.PKL, "rb"))
parts = {}
for k, v in d["res"].items():
    parts.setdefault(k.split("#")[0], []).append(v)
P = {k: C.pool(v) for k, v in parts.items()}
ST = {k: C.Stat(P[k], seed=sum(map(ord, k))) for k in P}              # same seeds as the main report
LW = {k: P[k]["spec"]["flow"]["l"] / np.sqrt(P[k]["spec"]["Deff"]) for k in P if P[k]["spec"]["flow"]}
S50 = lambda k: ST[k].pt[0]
BS = lambda k: ST[k].bs[:, 0]


def ci(x, qs=(2.5, 16, 84, 97.5)):
    x = x[np.isfinite(x)]
    return np.percentile(x, qs)


def pexp(ka, kb):
    pt = np.log(S50(kb) / S50(ka)) / np.log(LW[kb] / LW[ka])
    b = np.log(BS(kb) / BS(ka)) / np.log(LW[kb] / LW[ka])
    c = ci(b)
    return pt, c


print("(a) LENGTH EXPONENT p, pre-freeze-out time matched (S50, bootstrap over runs)")
for k in ("V32_3.6", "V32_3.6m", "V32_10", "V32_30"):
    s = P[k]["spec"]
    print(f"  {k:9s} l/w {LW[k]:5.2f}  gamma {s['gamma']:g}  pre-time {1 / s['gamma'] / C.t17(s['L'], s['Deff']):.2f} t17  kicks before freeze-out {1 / s['gamma'] / s['flow']['Tk']:5.1f}  "
          f"S50 {S50(k):.2f} [95% {ci(BS(k))[0]:.2f}-{ci(BS(k))[3]:.2f}]  walled {ST[k].pt[3]:.3f}")
for lab, ka, kb in (("3.6 -> 10, NOT matched (V32_3.6 at gamma 0.01; this is P8)", "V32_3.6", "V32_10"),
                    ("3.6 -> 10, matched (V32_3.6m)", "V32_3.6m", "V32_10"),
                    ("10 -> 30 (both matched; P8 printed only the point)", "V32_10", "V32_30"),
                    ("3.6 -> 30, matched (V32_3.6m, V32_30)", "V32_3.6m", "V32_30")):
    pt, c = pexp(ka, kb)
    print(f"  p({lab}): {pt:.2f}  [68% {c[1]:.2f}..{c[2]:.2f}; 95% {c[0]:.2f}..{c[3]:.2f}]")
m3 = ("V32_3.6m", "V32_10", "V32_30")
x = np.array([LW[k] for k in m3])
pt = H.fit_p(x, np.array([S50(k) for k in m3]))
b = ci(np.array([H.fit_p(x, v) for v in np.stack([BS(k) for k in m3], 1)]))
print(f"  p (least-squares slope over V32_3.6m, V32_10, V32_30): {pt:.2f}  [68% {b[1]:.2f}..{b[2]:.2f}; 95% {b[0]:.2f}..{b[3]:.2f}]")
print(f"  S50(V32_3.6m) / S50(V32_3.6) (same flow, gamma 0.0556 vs 0.01): {S50('V32_3.6m') / S50('V32_3.6'):.2f}  [95% {ci(BS('V32_3.6m') / BS('V32_3.6'))[0]:.2f}..{ci(BS('V32_3.6m') / BS('V32_3.6'))[3]:.2f}]")
print("  the same comparison at D_eff 8 and 64 has no matched l/w 3.6 cell (R8_3.6 and V64_3.6 were run at gamma 0.01 only)")

print("\n(b) V32_30, run by run (time since the end of the sweep / t17 = 1173.9; Delta is the sweep-end bias)")
c = P["V32_30"]
t17_ = C.t17(c["spec"]["L"], c["spec"]["Deff"])
o = np.argsort(c["t_single"])
for i in o:
    print(f"  Delta {c['deltas'][i]:.1f}  walled {bool(c['walled_at_t0'][i])}  domains at freeze-out {int(c['ndom0'][i])}  clear time {c['t_single'][i] / t17_:.2f} t17  (last observation {c['tobs'][i] / t17_:.2f})")
print(f"  median of the 8 clear times {np.median(c['t_single']) / t17_:.2f} t17; KM S50 {S50('V32_30'):.2f}; range {c['t_single'].min() / t17_:.2f}..{c['t_single'].max() / t17_:.2f}")

print("\n(c) RATIO DEPENDENCE at l/w 10 (matched design: 5 kicks and 1.07 t17 before freeze-out at all three ratios)")
r3 = np.array([8.0, 32.0, 64.0])
k3 = ("R8_10", "V32_10", "V64_10")
q_pt = H.fit_p(r3, np.array([S50(k) for k in k3]))
qb = np.array([H.fit_p(r3, v) for v in np.stack([BS(k) for k in k3], 1)])
qc = ci(qb)
print(f"  S50(l/w 10): 1/8 {S50('R8_10'):.2f}, 1/32 {S50('V32_10'):.2f}, 1/64 {S50('V64_10'):.2f}; q = {q_pt:.2f} [95% {qc[0]:.2f}..{qc[3]:.2f}]")
print(f"  molecular reference at l/w 10, 1/32 (K32): {S50('K32'):.2f}; flow / molecular {S50('V32_10') / S50('K32'):.2f}; "
      f"if S_mol is proportional to D_eff/D_mol (h14b J 9.33, K 13.36 at 1/8; here K32 {S50('K32'):.1f} at 1/32), K at 1/64 would be ~{S50('K32') * 2:.0f}")
s64 = S50("V64_10")
sib = H.parse_sibling()
fbe = sorted({min(f for f in sib[w]["f5"][U][0] if f is not None) for w in sib for U in ("U0.02", "U1")})
print("  ratio D_eff/D_mol at which S50(l/w 10) = S50(V64_10) x (ratio/64)^q reaches f:   [point q | q at its 97.5% bound | q at its 2.5% bound]")
for f in fbe:
    print(f"    f = {f:6g}: " + " | ".join(f"{64 * (f / s64) ** (1 / q):.2e}" for q in (q_pt, qc[3], qc[0])))
print("  S50(l/w 10) at ratio 1e9 (the real order of magnitude) IF the power law held over 7 decades (NOT a measurement):")
for lab, q in (("q point", q_pt), ("q 97.5% bound", qc[3]), ("q 2.5% bound", qc[0])):
    print(f"    {lab} {q:.2f}: {s64 * (1e9 / 64) ** q:.3g}")
nw = {f: sum(min(g for g in sib[w]['f5'][U][0] if g is not None) == f for w in sib for U in ('U0.02', 'U1')) for f in fbe}
print("  (world, U) rows by break-even f (of 30): " + ", ".join(f"{f:g}: {n}" for f, n in nw.items()))

print("\n(d) P13 table re-read with the measured plateau: S_ext(lw) = S50(V64_10) x (max(lw, 10) / 10)^p, p = p(10 -> 30) (ratio 64 only)")
pt1030, c1030 = pexp("V32_10", "V32_30")
for lab, pp in (("point", pt1030), ("97.5% bound", c1030[3])):
    cnt = {"U0.02": 0, "U1": 0}
    for w in sib:
        for U in cnt:
            fe = min(f for f in sib[w]["f5"][U][0] if f is not None)
            cnt[U] += s64 * (max(sib[w]["lw"][U][1], 10) / 10) ** pp < fe
    print(f"  p = {pp:.2f} ({lab}): worlds (of 15) with S_ext(p90 l/w) < f_be at ratio 64: U 0.02 {cnt['U0.02']}, U 1 {cnt['U1']}")
print("  p90 l/w above 30 (up to 690, h14_eddy.log:60-74) is outside the measured range; this is a plateau assumption, not a measurement")
