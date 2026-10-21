"""h14b_advect_posthoc: POST HOC, exploratory, NOT pre-registered. Reads h14b_advect_raw.pkl written by h14b_advect.py and
asks the questions the pre-registered report could not answer because S80 and the l/w = 30 S50 are censored:
 (a) clearing fractions (Kaplan-Meier) at fixed multiples of t17 in G, F, H, I against the molecular-only references J, K;
 (b) the local exponent between l/w = 3.6 (B) and 10 (G), and the n = 2 versus n = 1 ratio at l/w = 10 (H / G);
 (c) kicks needed (S50 t17 / Tk);
 (d) an ILLUSTRATIVE extrapolation of S50 with the 3.6 -> 10 exponent against the break-even f of h14_eddy.log (beyond l/w = 30
     and at D_mol/D_eff = 1/8; it is not a measurement).
Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python h14b_advect_posthoc.py > h14b_advect_posthoc.log
"""
import pickle

import numpy as np

import h14b_advect as H          # asserts the predictions-file sha256 on import

res = pickle.load(open(H.PKL, "rb"))["res"]
parts = {}
for k, v in res.items():
    parts.setdefault(k.split("#")[0], []).append(v)
P = {k: H.pool(v) for k, v in parts.items()}
ST = {k: H.Stat(P[k], seed=sum(map(ord, k))) for k in P}        # same seeds as the main report
LW = {k: P[k]["spec"]["L"] / (2 * np.pi * P[k]["spec"]["n"]) / H.W_WALL for k in P if P[k]["spec"]["flow"]}


def km_cleared(c, ts, idx=None):
    """Kaplan-Meier fraction cleared by each t (in t17 units) among runs walled at freeze-out."""
    idx = np.arange(c["runs"]) if idx is None else idx
    w = c["walled_at_t0"][idx]
    T, ev = c["tobs"][idx][w] / H.t17(c["spec"]["L"]), np.isfinite(c["t_single"][idx][w])
    o = np.lexsort((~ev, T))
    T, ev = T[o], ev[o]
    S, out, n = 1.0, [], len(T)
    surv = np.ones(len(ts))
    for i in range(n):
        if ev[i]:
            S *= 1 - 1.0 / (n - i)
            surv[ts >= T[i]] = S
    return 1 - surv


def boot_cleared(k, ts, B=2000):
    c = P[k]
    rng = np.random.default_rng(sum(map(ord, k)) + 7)
    b = np.array([km_cleared(c, ts, rng.integers(0, c["runs"], c["runs"])) for _ in range(B)])
    return km_cleared(c, ts), np.percentile(b, [2.5, 97.5], axis=0)


print("(a) FRACTION CLEARED (KM, runs walled at freeze-out) by t17 multiples; 95% bootstrap in brackets; hold in t17 given")
ts = np.array([1.0, 2.0, 4.0, 6.0, 8.0])
for k in ("B", "G", "F", "H", "I", "J", "K"):
    pt, ci = boot_cleared(k, ts)
    hold = P[k]["spec"]["hold"] / H.t17(P[k]["spec"]["L"])
    print(f"  {k} (L={P[k]['spec']['L']:g}, n={P[k]['spec']['n']}, " + (f"l/w {LW[k]:.1f}" if k in LW else "no flow") + f", hold {hold:.1f}): " +
          "  ".join((f"t={t:g}: {p:.2f} [{lo:.2f}-{hi:.2f}]" if t <= hold + 1e-9 else f"t={t:g}: beyond hold") for t, p, lo, hi in zip(ts, pt, ci[0], ci[1]))
          + f"  (walled at freeze-out {int(P[k]['walled_at_t0'].sum())} of {P[k]['runs']})")
n_i = int(np.isfinite(P["I"]["t_single"][P["I"]["walled_at_t0"]]).sum())
print(f"  I: {n_i} of {P['I']['runs']} cleared within the hold; median over the 14 uncleared runs not defined; "
      f"cleared times (t17): {sorted(np.round(P['I']['t_single'][np.isfinite(P['I']['t_single'])] / H.t17(533.0), 2).tolist())}")

print("\n(b) LOCAL EXPONENT AND n-DEPENDENCE (S50, bootstrap over runs)")
lwB, lwG = LW["B"], LW["G"]
pb = np.log(ST["G"].bs[:, 0] / ST["B"].bs[:, 0]) / np.log(lwG / lwB)
pb = pb[np.isfinite(pb)]
p_pt = np.log(ST["G"].pt[0] / ST["B"].pt[0]) / np.log(lwG / lwB)
p_f = np.log(ST["F"].pt[0] / ST["A"].pt[0]) / np.log(LW["F"] / LW["A"])
pfb = np.log(ST["F"].bs[:, 0] / ST["A"].bs[:, 0]) / np.log(LW["F"] / LW["A"])
pfb = pfb[np.isfinite(pfb)]
print(f"  p(3.6 -> 10) from G/B: {p_pt:.2f}  [68% {np.percentile(pb, 16):.2f}-{np.percentile(pb, 84):.2f}; 95% {np.percentile(pb, 2.5):.2f}-{np.percentile(pb, 97.5):.2f}]")
print(f"  p(3.6 -> 10) from F/A (both h = 0.5): {p_f:.2f}  [95% {np.percentile(pfb, 2.5):.2f}-{np.percentile(pfb, 97.5):.2f}]")
print(f"  p(10 -> 30) lower bound from the censoring floor of I ({ST['I'].cens_floor():.2f}): {np.log(ST['I'].cens_floor() / ST['G'].pt[0]) / np.log(LW['I'] / lwG):.2f}; "
      f"4 of 18 cleared (KM S(8) = {1 - km_cleared(P['I'], np.array([7.99]))[0]:.2f})")
for a, b, lab in (("H", "G", "n=2 / n=1 at l/w = 10"), ("C", "B", "n=2 / n=1 at l/w = 3.6"), ("D", "B", "n=4 / n=1 at l/w = 3.6"), ("E", "B", "n=6 / n=1 at l/w = 3.6")):
    ok = np.isfinite(ST[a].bs[:, 0]) & np.isfinite(ST[b].bs[:, 0])
    r = ST[a].bs[ok, 0] / ST[b].bs[ok, 0]
    print(f"  S50 ratio {lab}: {ST[a].pt[0] / ST[b].pt[0]:.2f}  [95% {np.percentile(r, 2.5):.2f}-{np.percentile(r, 97.5):.2f}]  (share of replicates with both defined {ok.mean():.2f})")

print("\n(c) KICKS NEEDED  N = S50 t17 / Tk  (Tk/t17 is the same 0.237 for every n = 1 case)")
for k in ("A", "B", "C", "D", "E", "F", "G", "H"):
    sp = P[k]["spec"]
    print(f"  {k}: Tk {sp['flow']['Tk']:g}  t17 {H.t17(sp['L']):.1f}  Tk/t17 {sp['flow']['Tk'] / H.t17(sp['L']):.3f}  S50 {ST[k].pt[0]:.2f}  N = {ST[k].pt[0] * H.t17(sp['L']) / sp['flow']['Tk']:.1f}")
print(f"  I: Tk 277.4, first kick 277 t.u. after the start (77 t.u. after freeze-out); {P['I']['spec']['hold'] / 277.4:.0f} kicks in the hold")

print("\n(d) ILLUSTRATIVE EXTRAPOLATION (NOT a measurement): S(lw) = S50(G) x (lw / 10.02)^p with p from (b), compared with f_be")
sib = H.parse_sibling()
for p, nm in ((p_pt, f"p = {p_pt:.2f} (point)"), (float(np.percentile(pb, 97.5)), f"p = {np.percentile(pb, 97.5):.2f} (97.5% bound)"), (1.1, "p = 1.1 (sibling)")):
    tot = {"U0.02": [0, 0], "U1": [0, 0]}
    for w, s_ in sib.items():
        for U in ("U0.02", "U1"):
            fe = min(f for f in s_["f5"][U][0] if f is not None)
            med, p90 = s_["lw"][U]
            tot[U][0] += ST["G"].pt[0] * (max(med, lwG) / lwG) ** p < fe
            tot[U][1] += ST["G"].pt[0] * (max(p90, lwG) / lwG) ** p < fe
    print(f"  {nm}: worlds (of 15) with S < f_be at the median l/w: U 0.02 {tot['U0.02'][0]}, U 1 {tot['U1'][0]}; at the p90 l/w: U 0.02 {tot['U0.02'][1]}, U 1 {tot['U1'][1]}")
print("  rows (p = point estimate): world, U, median S, p90 S, f_be")
for w, s_ in sib.items():
    for U in ("U0.02", "U1"):
        fe = min(f for f in s_["f5"][U][0] if f is not None)
        med, p90 = s_["lw"][U]
        sm, sp90 = (ST["G"].pt[0] * (max(x, lwG) / lwG) ** p_pt for x in (med, p90))
        print(f"    {w:18s} {U:6s} lw med {med:5.1f} -> S {sm:7.1f}   p90 {p90:6.1f} -> S {sp90:8.1f}   f_be {fe:g}   {'p90 below f_be' if sp90 < fe else 'p90 ABOVE f_be'}")

r80 = ST["B"].pt[1] / ST["B"].pt[0]
print(f"\n(d2) SAME WITH AN S80 PROXY: S80 ~ S50(G) x r80 x (lw/10.02)^p, r80 = S80(B)/S50(B) = {r80:.2f} (the pre-registered rule is about S80; G's S80 is censored)")
for p, nm in ((p_pt, f"p = {p_pt:.2f}"), (float(np.percentile(pb, 97.5)), f"p = {np.percentile(pb, 97.5):.2f}")):
    c = {"U0.02": 0, "U1": 0}
    for w, s_ in sib.items():
        for U in ("U0.02", "U1"):
            fe = min(f for f in s_["f5"][U][0] if f is not None)
            c[U] += ST["G"].pt[0] * r80 * (max(s_["lw"][U][1], lwG) / lwG) ** p < fe
    print(f"  {nm}: worlds (of 15) with S80 proxy < f_be at the p90 l/w: U 0.02 {c['U0.02']}, U 1 {c['U1']}")
