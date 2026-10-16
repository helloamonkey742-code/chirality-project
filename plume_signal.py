"""What would an Enceladus plume sample show if the ocean were a chance or weak-force amplifier?  Reading rule of A. Davila
(NASA Ames; email 2026-10-16), run on the model's own Enceladus draws (2026-10-16).

THE RULE (paraphrased).  The ocean's organics are probably dominated by an abiotic, racemic background, so any network signal
is diluted.  An abiotic L-excess (for example from the weak force) could still be detectable in a plume sample.  Unless the
excess is significant AND appears in more than 2-3 amino acids, it would be read as abiotic.  He did not follow our question
whether a plume sample tells an ocean-wide excess from a local one (does the plume carry the hand of the whole ocean, or only
of the layer it comes from?).

(1) LAYERS.  multiwall.py records no layer hands per draw, so they are derived here from the model's own rules.  Every
    Enceladus draw has a lid (salinity 5-20 g/kg below S* = 20), and a lid keeps two hands (Section 3.4): a thin stratified
    layer under the ice (the plume source) and the convecting interior (the bulk).  Assumptions, all stated and none
    measured: (a) the two layers are separate regions with no molecular exchange during the window, so each takes its own
    hand; (b) a layer amplifies when the model's own gates pass (k2 c t >= 10, k2 c > 2 k_rac(Ala), lateral mixing
    (L/A)^2/D_h <= t); the gates do not depend on the layer's share, so they are the same for both layers; (c) a layer with
    volume share f has P(PVED hand) = Phi(F Delta(f)) with the draw's own g, k2, c, V r, t, F (Eqs. 10-12, Delta grows as
    sqrt(f)); (d) the lid layer's share f_lid is unknown: primary 0.1, scanned 0.01, 0.03, 0.1, 0.3; (e) the plume is the lid
    layer (primary); the alternative, that the plume carries interior water through the lid, is also reported.  Without a
    lid (a counterfactual for Enceladus) there is one ocean-wide hand, so plume and bulk agree by construction in every
    one-hand draw.  Layer hands are used as expectations, not sampled.
(2) DILUTION.  ee_plume = f x ee_net x s.  f = share of the sampled amino acid made by the network, scanned log-uniformly
    over 1e-4..1.  ee_net = the network's steady purity eta*(r) with r = k_rac(Ala)/(k2 c) (persist.py's analytic formula,
    k_d := k2 c as in kcheck.py; eta* = 0 where r >= 0.087, i.e. k2 c < 11.5 k_rac).  s = survival against racemization:
    in a pool made at rate P and lost at rate 1/tau_res, the excess is ee_net/(1 + 2 k_rac tau_res); tau_res (how long a
    molecule stays before it is sampled) is unknown: 1e3, 1e5, 1e6, 1e7 yr.  k_rac(Ala) = 9.13e-15 /s, Asp 2.73e-14 /s at
    273 K (sweep2.py K_RAC).  Detection limit: Mora et al. 2022 (abstract only; the papers are not open access) gives
    concentration limits for capillary electrophoresis (0.05-1 uM) and no limit on ee, so the ee limit T is a stated range,
    1, 2, 5, 10% ee.  Steel et al. 2017 (abstract) says the ambient ocean racemizes, with a biological plume ee below 10%.
(3) DAVILA'S RULE.  The number n of amino acids the network makes is unknown (1..6).  Each shares f, a hand and eta*.  The
    reading is "abiotic" if n <= 3 or ee < T, "life" if n > 3 and ee >= T.
    Not used: Davila & McKay 2014 (doi:10.1089/ast.2014.1150).  Its abstract is about "chance and necessity" in
    biochemistry and says nothing on abundance patterns or ee, so it is not cited for a discriminator.

Base: multiwall.py main rule and cascade x10, seed 0, N = 200,000 (U.N), seedfix.py stream ([s, 0] -> [s, 0, 1]), draws
and extras exactly as lid_melt.py (ex streams default_rng([seed, i]), i = Enceladus's index in sweep2.chains).
Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python plume_signal.py > plume_signal.log

PREDICTIONS (written 2026-10-16, before any run).  No number below has been computed from these draws; the ranges come from
the logged Enceladus shares (seedfix.log), kcheck.log's shares of k2 c above M x 2 k_rac (M = 1: 80-82% of icy-world draws;
M_ex = 5.75: 72.5-75%) and the closed forms of the dilution formulas.  Quantities marked A are over draws that pass the gates.
 (P1) Regression: Enceladus no lid, multiwall main 26.4 +- 0.05 % PVED-set and x10 19.4 +- 0.05 %; with the lid (salt mask,
      100% of draws) the final PVED-set share is 0.0 in both.
 (P2) Gates: the share of draws that pass them (A) is 60-85%, and it is the same for the lid layer and the bulk.
 (P3) Layer hands, f_lid = 0.1: the mean P(PVED hand) over A draws is 0.65-0.85 for the bulk and 0.55-0.70 for the lid
      layer, the layer lower than the bulk by at least 0.05, and the layer's value rises with f_lid across 0.01, 0.03, 0.1, 0.3.
 (P4) Disagreement: among A draws the expected share with opposite hands in lid layer and bulk is 0.25-0.50 (f_lid = 0.1).
 (P5) Network purity: of the A draws, 85-95% have eta* > 0 (k2 c above 11.5 k_rac); the median eta* over A draws is >= 0.97.
 (P6) Survival (Ala): s = 0.634 +- 0.01 at tau_res = 1e6 yr; Asp 0.37 +- 0.01; at tau_res = 1e7 yr no draw reaches an ee
      of 2, 5 or 10% at any f (s alone is 0.017).
 (P7) Needed network share (tau_res = 1e5 yr, Ala, over A draws with eta* > 0): median f_req = T / (eta* s) is 0.052-0.058 for
      T = 5% and 0.0104-0.0116 for T = 1%; no draw has f_req below T.
 (P8) Chance of detection, f log-uniform over 1e-4..1, tau_res = 1e5 yr, plume = lid layer, over A draws: T = 1%: 0.40-0.48,
      T = 5%: 0.26-0.32, T = 10%: 0.19-0.24.
 (P9) Davila: no cell of the n <= 3 rows is read as life (by construction); for n >= 4 (T = 5%, tau_res = 1e5 yr) the share of
      ALL lid draws read as life is 0.15-0.25, the same share as those read as abiotic with a detectable excess for n <= 3.
 (P10) Sign: among detectable plumes (T = 5%, tau_res = 1e5 yr) the expected share with the PVED hand is 0.55-0.70 for a
      plume from the lid layer, so five independent plumes all showing the PVED hand has probability 0.05-0.17.
 (P11) If the plume carries interior water, that share is higher than the lid layer's and lies in 0.65-0.85.
 (P12) Without a lid, the share of draws that reach one ocean-wide hand is 30-60% (main rule); plume and bulk agree in all of them.
"""
import os
import sys
import time

import numpy as np

ORIG_RNG = np.random.default_rng


def fixed_rng(seed=None, *a, **kw):   # seedfix.py's correction, copied (importing seedfix would run it)
    if isinstance(seed, (list, tuple)) and len(seed) == 2 and seed[1] == 0:
        seed = [seed[0], 0, 1]
    return ORIG_RNG(seed, *a, **kw)


np.random.default_rng = fixed_rng
from scipy.special import ndtr  # noqa: E402
import uncertainty as U   # noqa: E402
U.N = int(os.environ.get("PLUME_N", U.N))   # smoke test only
import sweep2 as S        # noqa: E402
import multiwall as M     # noqa: E402
import persist as P       # noqa: E402

YR = 3.156e7
FLID = (0.01, 0.03, 0.1, 0.3)
TAUS = (1e3, 1e5, 1e6, 1e7)
TS = (0.01, 0.02, 0.05, 0.10)
NS = (1, 2, 3, 4, 5, 6)
FLO = 1e-4
CHECK = []


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


def eta_table():
    """eta*(r) from persist.py's analytic pure branch (kd = k2 c); zero from the exist limit r_c = 0.08698 upward."""
    R = np.logspace(-6, np.log10(0.0869), 300)
    E = np.array([np.sqrt(max(P.eta2_of(P.par(r=r)), 0.0)) for r in R])
    return lambda r: np.where(r <= R[-1], np.interp(np.log(np.minimum(r, R[-1])), np.log(R), E), 0.0)


def surv(kr, tau):   # steady pool: excess = ee_net / (1 + 2 k_rac tau_res)
    return 1.0 / (1.0 + 2 * kr * tau * YR)


def pdet(e, T):      # P(f >= T/e) for f log-uniform on [1e-4, 1] (e = eta* s)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.clip(np.log(e / T) / np.log(1 / FLO), 0.0, 1.0)


def main(seed=0):
    t0 = time.time()
    _, W = S.chains(seed)
    name = "Enceladus"
    i = list(W).index(name)
    x, win, a, salt = W[name]
    L, dz = a[1], a[4]
    rng = np.random.default_rng([seed, i])
    ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
    none = np.zeros(U.N, bool)
    print(f"plume_signal.py, N = {U.N}, seed 0, Enceladus (sweep2.chains index {i}); lid in {100 * salt.mean():.1f}% of draws.")

    print("\n(0) Whole-ocean outcome, multiwall rule (P_pved % | one-hand %)")
    w0, r0, _ = M.multi(x, L, dz, ex, none)
    w10, r10, _ = M.multi(x, L, dz, ex, none, scale=10.0)
    wl, rl, _ = M.multi(x, L, dz, ex, salt)
    wl10, _, _ = M.multi(x, L, dz, ex, salt, scale=10.0)
    print(f"  no lid main : {100 * w0.mean():5.2f} | {100 * r0.mean():5.2f}      x10: {100 * w10.mean():5.2f} | {100 * r10.mean():5.2f}")
    print(f"  with lid    : {100 * wl.mean():5.2f} | {100 * rl.mean():5.2f}      x10: {100 * wl10.mean():5.2f}")
    print("  No lid: one hand fills the ocean, so plume and bulk carry the same hand in every one-hand draw.")

    c, k2, Dh, t, V, F = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time", "volume", "chem. factor F"))
    g = S.remap(x["bias g"], S.G_OLD, S.G_NEW) * ex["eta"]
    k = k2 * c
    krA, krD = S.K_RAC["Ala"], S.K_RAC["Asp"]
    A = (k * t >= 10) & (k > 2 * krA) & ((L / U.A) ** 2 / Dh <= t)
    phi = lambda f: ndtr(F * M.delta_at(c, g, k2, V * ex["r"], t, f))
    print(f"\n(1) Layers under the lid. Gates pass (A) in {100 * A.mean():.1f}% of draws (the same for both layers).")
    print("  f_lid | mean P(PVED hand | A): lid layer, bulk | opposite hands in layer and bulk, share of A | "
          "layer Delta>=2 share of A, bulk")
    top, bulk, opp = {}, {}, {}
    for fl in FLO_GRID if False else FLID:
        pt, pb = phi(fl), phi(1 - fl)
        top[fl], bulk[fl] = pt[A].mean(), pb[A].mean()
        opp[fl] = (pt * (1 - pb) + (1 - pt) * pb)[A].mean()
        d2t = (F * M.delta_at(c, g, k2, V * ex["r"], t, fl) >= 2)[A].mean()
        d2b = (F * M.delta_at(c, g, k2, V * ex["r"], t, 1 - fl) >= 2)[A].mean()
        print(f"  {fl:5} |  {top[fl]:.3f}  {bulk[fl]:.3f} |  {opp[fl]:.3f} |  {d2t:.3f}  {d2b:.3f}")
    fl = 0.1
    pt, pb = phi(fl), phi(1 - fl)
    print(f"  Over ALL draws at f_lid = 0.1: layer PVED hand {100 * (A * pt).mean():.1f}%, layer other hand "
          f"{100 * (A * (1 - pt)).mean():.1f}%, layer racemic {100 * (~A).mean():.1f}%; bulk PVED {100 * (A * pb).mean():.1f}%, "
          f"other {100 * (A * (1 - pb)).mean():.1f}%; opposite hands {100 * (A * (pt * (1 - pb) + (1 - pt) * pb)).mean():.1f}%, "
          f"same hand {100 * (A * (pt * pb + (1 - pt) * (1 - pb))).mean():.1f}%, both racemic {100 * (~A).mean():.1f}%.")

    eta = np.where(k > 0, eta_table()(krA / k), 0.0)
    print(f"\n(2) Dilution. eta*(k_rac/k2c): >0 in {100 * (eta[A] > 0).mean():.1f}% of A draws; median {np.median(eta[A]):.4f}, "
          f"10th percentile {np.percentile(eta[A], 10):.3f}, share >= 0.9: {100 * (eta[A] >= 0.9).mean():.1f}%.")
    print("  Survival s = 1/(1 + 2 k_rac tau_res): " + "; ".join(
        f"tau {tau:g} yr: Ala {surv(krA, tau):.4f} Asp {surv(krD, tau):.4f}" for tau in TAUS))
    print("  Needed network share f_req = T/(eta* s), Ala, over A draws with eta* > 0 (median [10th-90th pct]; share of A draws "
          "with f_req > 1, i.e. unreachable):")
    freq_med = {}
    mn = np.inf
    for tau in TAUS:
        s = surv(krA, tau)
        row = []
        for T in TS:
            m = A & (eta > 0)
            fr = T / (eta[m] * s)
            mn = min(mn, (fr / T).min())
            freq_med[tau, T] = np.median(fr)
            row.append(f"T {100 * T:g}%: {np.median(fr):.4f} [{np.percentile(fr, 10):.4f}-{np.percentile(fr, 90):.4f}] "
                       f"unreach {100 * (1 - (fr <= 1).mean() * m.sum() / A.sum()):.1f}%")
        print(f"   tau {tau:g} yr: " + " | ".join(row))
    print("  P(plume ee >= T), f log-uniform on [1e-4, 1], plume = lid layer. Over A draws | over all draws:")
    det, pd_ = {}, {}
    for tau in TAUS:
        e = eta * surv(krA, tau)
        for T in TS:
            pd_[tau, T] = pdet(e, T)
            det[tau, T] = (pd_[tau, T][A].mean(), (A * pd_[tau, T]).mean())
        print(f"   tau {tau:g} yr: " + " | ".join(f"T {100 * T:g}%: {det[tau, T][0]:.3f} | {det[tau, T][1]:.3f}" for T in TS))
    print("  (The same numbers hold if the plume carries interior water: eta*, s and the gates do not depend on the layer; "
          "only the sign probability does.)")

    print("\n(3) Davila's rule at tau_res = 1e5 yr (n = amino acids the network makes, all sharing f, hand and eta*). "
          "Share of ALL lid draws read as: abiotic with a detectable excess / life / racemic or below the limit")
    for T in TS:
        w = (A * pd_[1e5, T]).mean()
        print(f"  T {100 * T:g}%: " + "; ".join(
            f"n={n}: {100 * (w if n <= 3 else 0):.1f} / {100 * (w if n > 3 else 0):.1f} / {100 * (1 - w):.1f}" for n in NS))
    w5 = (A * pd_[1e5, 0.05]).mean()
    ok_zero = all(True for n in NS if n <= 3)   # life-read share for n <= 3 is 0 by the rule's definition

    wt = A * pd_[1e5, 0.05]
    q_top = (wt * pt).sum() / wt.sum()
    q_bulk = (wt * pb).sum() / wt.sum()
    print(f"\n(4) Sign (T 5%, tau 1e5 yr, f_lid 0.1): among detectable plumes the PVED-hand share is {q_top:.3f} for a "
          f"lid-layer plume and {q_bulk:.3f} for an interior plume; five independent plumes all PVED: "
          f"{q_top ** 5:.3f} (lid layer) / {q_bulk ** 5:.3f} (interior); all five of one hand if hands were 50/50: {2 * 0.5 ** 5:.3f}.")
    print("  A sign cannot separate the weak force from a chance amplifier with one ocean: a chance amplifier gives each hand "
          "with probability 0.5 and the PVED about q.")

    print("\nPredictions:")
    check("P1", abs(100 * w0.mean() - 26.4) <= 0.05 and abs(100 * w10.mean() - 19.4) <= 0.05 and wl.mean() == 0 and wl10.mean() == 0,
          f"no lid {100 * w0.mean():.2f} / {100 * w10.mean():.2f}; lid {100 * wl.mean():.2f} / {100 * wl10.mean():.2f}")
    check("P2", 0.60 <= A.mean() <= 0.85, f"gates pass in {100 * A.mean():.1f}%")
    ok = 0.65 <= bulk[0.1] <= 0.85 and 0.55 <= top[0.1] <= 0.70 and bulk[0.1] - top[0.1] >= 0.05 and \
        top[0.01] < top[0.03] < top[0.1] < top[0.3]
    check("P3", ok, f"f_lid 0.1: bulk {bulk[0.1]:.3f}, layer {top[0.1]:.3f}; layer by f_lid " +
          "/".join(f"{top[f]:.3f}" for f in FLID))
    check("P4", 0.25 <= opp[0.1] <= 0.50, f"opposite hands {opp[0.1]:.3f} of A draws")
    med = np.median(eta[A])
    check("P5", 0.85 <= (eta[A] > 0).mean() <= 0.95 and med >= 0.97, f"eta*>0 in {100 * (eta[A] > 0).mean():.1f}% of A; median {med:.4f}")
    mx = max(pd_[1e7, T][A].max() for T in (0.02, 0.05, 0.10))
    check("P6", abs(surv(krA, 1e6) - 0.634) <= 0.01 and abs(surv(krD, 1e6) - 0.37) <= 0.01 and mx == 0,
          f"s(1e6) Ala {surv(krA, 1e6):.3f} Asp {surv(krD, 1e6):.3f}; largest detect prob at 1e7 yr, T >= 2%: {mx:.3g}")
    ok = 0.052 <= freq_med[1e5, 0.05] <= 0.058 and 0.0104 <= freq_med[1e5, 0.01] <= 0.0116 and mn >= 1 - 1e-12
    check("P7", ok, f"median f_req {freq_med[1e5, 0.05]:.4f} (T 5%), {freq_med[1e5, 0.01]:.4f} (T 1%); min f_req/T {mn:.3f}")
    d1, d5, d10 = (det[1e5, T][0] for T in (0.01, 0.05, 0.10))
    check("P8", 0.40 <= d1 <= 0.48 and 0.26 <= d5 <= 0.32 and 0.19 <= d10 <= 0.24, f"T 1/5/10%: {d1:.3f} / {d5:.3f} / {d10:.3f}")
    check("P9", ok_zero and 0.15 <= w5 <= 0.25, f"n >= 4 life-read share of all lid draws {w5:.3f} (same share abiotic-read for n <= 3)")
    check("P10", 0.55 <= q_top <= 0.70 and 0.05 <= q_top ** 5 <= 0.17, f"lid-layer PVED share {q_top:.3f}; five all PVED {q_top ** 5:.3f}")
    check("P11", q_bulk > q_top and 0.65 <= q_bulk <= 0.85, f"interior-plume PVED share {q_bulk:.3f} vs layer {q_top:.3f}")
    check("P12", 0.30 <= r0.mean() <= 0.60, f"one-hand share, no lid, main: {100 * r0.mean():.1f}%")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held. {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
