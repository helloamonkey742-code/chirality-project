"""Amplifier vs life vs plain abiotic: a decision rule with error rates.  Reuses plume_signal.py / sweep2.py / uncertainty.py
(imported, never edited).  Written 2026-10-16.

QUESTION.  Given measured ee in n amino acids from a plume, how often is the call (abiotic amplifier / life / plain abiotic)
right, and how many amino acids and what precision are needed for >= 90% correct calls?  Which single feature matters most?

SOURCES (all DOIs resolved through api.crossref.org on 2026-10-16, no email in the URL; abstracts or open pages only, no
downloads).
 - Cohen & Chyba 2000, Icarus 145:272-281, doi 10.1006/icar.1999.6328.  Table II (free amino acids in solution), k in 1/s =
   exp(A - B/T).  This is the ONLY source here with rate laws for seven amino acids; I read it from an open copy.  Values at 273 K:
        aa    A      B       k(273 K) /s   x Ala   Ea kJ/mol   2k tau = 1 at tau (yr)
        Ala  24.91  15863    3.829e-15     1.00    132         4.1e6
        Asp  23.34  14725    5.147e-14    13.44    122         3.1e5
        Glu  22.56  15675    7.270e-16     0.19    130         2.2e7
        Ile  27.24  16913    8.406e-16     0.22    141         1.9e7
        Leu  24.30  15750    3.147e-15     0.82    131         5.0e6
        Phe  19.82  13440    1.687e-13    44.06    112         9.4e4
        Val  23.61  15754    1.556e-15     0.41    131         1.0e7
   The ee decay convention is exp(-2kt), the same as the model's s = 1/(1 + 2 k tau).
 - Model's own k_rac (sweep2.K_RAC, 273 K): Ala 9.13e-15, Asp 2.73e-14 /s, from Bada & Schroeder 1975, Naturwissenschaften 62:71,
   doi 10.1007/bf00592179, as tabulated in Steen, Joergensen, Lomstein 2013, PLoS ONE 8:e71648, doi 10.1371/journal.pone.0071648
   (Table 3, free solution at 3 C: Asx 0.15 and Ala 0.046 in 1e-5 /yr, consistent with the model's values).  The two compilations
   DISAGREE: Asp/Ala is 13.4 in Cohen & Chyba and 3.0 in the model.  So the ordering of k and the ratios are compilation-dependent;
   I run both ("rate set C" = Cohen & Chyba, primary; "rate set M" = Cohen & Chyba with Ala and Asp replaced by the model's).
 - Steen et al. 2013 also give SEDIMENT consensus rates (bound amino acids, not free solution) at 3 C: Asx 3.2, Glx 0.73, Ser 2.5,
   Ala 1.7 (1e-5 /yr), Ea 110/106/101/93 kJ/mol.  They are not free-solution rates and are NOT used.
 - Steel, Davila, McKay 2017, Astrobiology 17:862, doi 10.1089/ast.2017.1673 (abstract only): racemization is fast against
   production in the Enceladus ocean, plume ee for their methanogen biology is below 10%.  Used for the premise that life's amino
   acids racemize with the same s_i once dissolved; their rate constants were NOT obtained.
 - Glavin & Dworkin 2009, PNAS 106:5487, doi 10.1073/pnas.0811618106 (abstract): L-isovaline ee 18.5 +/- 2.6% (Murchison), 15.2 +/-
   4.0% (Orgueil).  Basis for the plain-abiotic "meteorite-like" scenario.
 - Cronin & Pizzarello 1997, Science 275:951, doi 10.1126/science.275.5302.951: DOI verified, abstract NOT read here; cited only
   because Cohen & Chyba cite it for isovaline's resistance to racemization.
 - Named in the request but NOT usable: Bada 1972 JACS 94:1371 (doi 10.1021/ja00759a064; pH dependence, no abstract available to
   me), Bada & Schroeder 1972 EPSL 15:1 (doi 10.1016/0012-821x(72)90022-2; isoleucine in sediments), Smith & Reddy 1989 JOC 54:4529
   (doi 10.1021/jo00280a017; side-chain effect, title only), Smith & Evans 1980 (no open page with numbers found).  No rate from
   them is used.  Serine, norvaline and alpha-aminobutyric acid: no free-solution rate found, so they are LEFT OUT.

ASSUMPTIONS (stated, not measured).
 (a) Plume ee_i = f_i x eta_i x s_i x hand, s_i = 1/(1 + 2 k_i tau), tau = residence time, log-uniform 1e3..1e7 yr.  f_i = share of
     amino acid i that is produced by the process (network or life), log-uniform on 1e-4..1; in 50% of draws one f is shared by all
     amino acids, in 50% each amino acid has its own f.
 (b) Panel: protein amino acids in a fixed order Ala, Asp, Phe, Glu, Leu, Val, Ile (Ile epimerization has Keq = 1.3; treating it
     as a plain ee decay is approximate) plus isovaline ("Iva"), a non-protein chiral alpha-dialkyl amino acid with k = 0 assumed
     (resistant to racemization; qualitative, from the Cohen & Chyba text citing Cronin & Pizzarello).  Panels: protein only
     n = 2, 4, 6, 7, and n-1 protein + Iva for n = 2, 4, 6, 8.
 (c) Scenario A (amplifier): draws from the model's Enceladus draws (seed 0, N = 200,000, sweep2.chains(0)) that pass the Ala
     gates (k t >= 10, k > 2 k_rac(Ala), lateral mixing (L/A)^2/D_h <= t).  Each amino acid i has its own gate k > 2 k_i and purity
     eta_i = eta*(k_i/(k2 c)) from persist.py (PS.eta_table()), zero if the gate fails.  Iva: eta*(0), only the other two gates;
     it carries the network's ee with probability q_np = 0.5, else it is racemic.  Hand: L (sign +) in every draw, the worst case
     for telling A from life (sensitivity: 50/50).
 (d) Scenario B (life): protein amino acids released with purity uniform on 0.97..1 (Steel's premise), then racemizing with the
     same s_i; Iva racemic (ee 0).  Same f and tau priors as A.  Hand L.
 (e) Scenario C (plain abiotic): 50% fully racemic; 50% meteorite-like: Iva ee uniform 2..20% (L), each protein amino acid has
     with probability 0.5 a small L ee uniform 0..5% times s_i.
 (f) Measurement noise: observed ee = true ee + N(0, sigma) per amino acid, sigma = 1, 2, 5% absolute.
 (g) Priors: equal (1/3 each).  Classifier: HistGradientBoosting on the observed ee vector plus count and mean of protein ee
     (a flexible stand-in for the Bayes classifier; it can only under-state the best achievable accuracy).  60,000 train and 60,000
     test samples per class; accuracy standard error about 0.0012.  The model has NO information about real plumes: every number
     below is conditional on A-C above.
 (h) "Visible" = at least one panel amino acid has true |ee| >= 5%.  "acc_V_AB" = fraction of correct calls among visible true-A
     and true-B samples.

Run: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python life_test.py > life_test.log

PREDICTIONS (written 2026-10-16, before any run of this script).  Disclosure: before writing these I read abiotic_vs_life.py/.log,
plume_signal.py/.log, persist.py (headers); I computed the Cohen & Chyba k(273 K) table by a one-line arithmetic check; I ran no
simulation of this script and no probe of its pool.  Ranges below are hand estimates (e.g. visible share ~0.3 per amino acid from
the logged detection chances; density of draws near the gate ~4-6% per e-fold of k from plume_signal.log / abiotic_vs_life.log).
 (P1) Regression, rates and survival: k(273 K) Ala 3.83e-15, Asp 5.15e-14, Glu 7.27e-16, Ile 8.41e-16, Leu 3.15e-15, Phe 1.69e-13,
      Val 1.56e-15, each within 1%; order Phe > Asp > Ala > Leu > Val > Ile > Glu; Asp/Ala 13.4 +/- 0.2; at tau 1e5 yr s_Ala 0.976 +/-
      0.003, s_Asp 0.755 +/- 0.005, s_Phe 0.484 +/- 0.005.
 (P2) Regression, population: with the model's Ala rate the gate-passing share A is 0.786-0.788 and eta*_Ala > 0 in 0.931-0.935 of A;
      with the Cohen & Chyba Ala rate A is larger, 0.79-0.84.
 (P3) Rate-gated purity (A draws, rate set C): share with eta* = 0: Ala 0.04-0.12, Asp 0.10-0.30, Phe 0.18-0.42, Glu 0.00-0.03, with
      the order Phe > Asp > Ala > Glu.
 (P4) Visibility (n = 8 panel with Iva): share of samples with some |ee| >= 5%: A 0.45-0.75, B 0.45-0.75, C 0.35-0.48.
 (P5) Unconditional 3-class accuracy, n = 8 with Iva, sigma 2%: 0.50-0.70.
 (P6) No cell of the 24-cell grid (two panel families x 4 n x 3 sigma) reaches 0.90; the maximum is below 0.75.  So there is no
      minimum n or precision for 90% correct calls: the limit is invisibility (A and B with all ee below noise look racemic).
 (P7) Weak dependence: with Iva, sigma 2%, acc(n = 8) - acc(n = 2) is 0.00-0.12; at n = 8, acc(sigma 1%) - acc(sigma 5%) is 0.02-0.15.
 (P8) A vs B among visible, sigma 2%: protein only n = 7: acc_V_AB 0.50-0.62 (near chance); with Iva n = 8: 0.62-0.82.
 (P9) Decisive feature: removing Iva (n = 8 with Iva -> n = 7 protein only) lowers acc_V_AB (sigma 2%) by >= 0.08, more than replacing the
      protein ee by count and mean only (loses which amino acids, ordering and ratios) and more than taking |ee| (loses sign); the
      last two each lower it by <= 0.04.
 (P10) Sign: with A hand 50/50, signed minus |ee| acc_V_AB is >= 0.05; with A hand L it is <= 0.02 in absolute value.
 (P11) Sensitivities, acc_V_AB (n = 8 with Iva, sigma 2%): q_np = 0: 0.50-0.62; q_np = 1: 0.72-0.90; hand 50/50: 0.70-0.88; rate set M:
       within +/- 0.03 of the main value.
 (P12) Davila's rule (call life if >= 4 of the 7 protein amino acids have observed ee >= 5%, sigma 2%): precision of the life call
       (equal priors) 0.50-0.62; the rest are mostly amplifiers.
 (P13) Confident calls (posterior >= 0.9, n = 8 with Iva, sigma 2%): at most 0.10 of true B samples get a confident life call; 0.08-0.30 of
       true A samples get a confident amplifier call; precision of all confident calls >= 0.90.
 (P14) Best case (hand 50/50, q_np = 1, sigma 1%, n = 8 with Iva): unconditional accuracy still below 0.80.

RESULTS (life_test.log, 2026-10-16): 11 of 14 held.  60 s.  The predictions were committed in d1eafb5 and are unchanged; the only check
before the run was py_compile; no smoke run or probe.  Rate-set M and the main config give the same numbers, so the compilation
disagreement (Asp/Ala 13.4 vs 3.0) does not change the conclusion.
 (P1) HELD   rates, order, Asp/Ala 13.44, s(1e5 yr) Ala 0.976, Asp 0.755, Phe 0.484.
 (P2) HELD   A (model Ala) 0.7872, eta*>0 0.9335; A (Cohen Ala) 0.8036.
 (P3) HELD   eta*=0 shares: Ala 0.049, Asp 0.186, Phe 0.266, Glu 0.002 (Leu 0.042, Val 0.020, Ile 0.005).
 (P4) HELD   visible: A 0.598, B 0.607, C 0.415.
 (P5) HELD   n = 8 with Iva, sigma 2%: 0.597.
 (P6) HELD   best of 24 cells 0.634 (with Iva, n = 8, sigma 1%); no cell near 0.90; acc_V_AB also never reaches 0.90 (max 0.640).
 (P7) HELD   acc(n8) - acc(n2) = 0.088; acc(sigma 1%) - acc(sigma 5%) = 0.086.
 (P8) FAILED protein only n = 7: acc_V_AB 0.526 (inside 0.50-0.62); with Iva n = 8: 0.606 (below the predicted 0.62-0.82, by 0.014).
 (P9) FAILED removing Iva lowers acc_V_AB by 0.0797 (predicted >= 0.08, missed by 0.0003); |ee| 0.013 and count+mean 0.026 held.  The
             direction held: Iva is the largest drop; an Iva-only classifier (0.470) is worse than the full one by 0.136.
 (P10) HELD  hand 50/50: signed minus |ee| 0.194; hand L 0.013.
 (P11) FAILED q_np = 0: 0.546 and hand 50/50: 0.786 and rate set M: 0.606 vs 0.606 held; q_np = 1: 0.673 (predicted 0.72-0.90).
 (P12) HELD  Davila's life call precision 0.523 (47.6% of its life calls are amplifiers; it finds 21.5% of life).
 (P13) HELD  confident life calls for 0.000 of true B; confident amplifier calls for 0.086 of true A; precision of confident calls 0.987.
 (P14) HELD  best case accuracy 0.734 (acc_V_AB 0.841).
 Reading.  (1) In this model a plume cannot be called life with confidence: life is never confirmed (0.000 of B samples reach a 0.9
 posterior), because an amplifier with a racemic non-protein amino acid looks the same.  (2) The amplifier CAN be confirmed, but only
 when a non-protein chiral amino acid (Iva) carries an ee (8.6% of A samples; 0.987 precision), or, if the amplifier's hand is
 random, when the hand is D (A-vs-B accuracy among visible rises from 0.606 to 0.786).  (3) More protein amino acids help only by
 making life or an amplifier visible against racemic, not by separating the two (protein only: acc_V_AB 0.52-0.55 at every n).
 (4) No cell reaches 90%: the largest unconditional accuracy is 0.634, 0.734 in the best case; undetected A and B (about 40% each)
 look like plain abiotic.  (5) Davila's >= 4 amino-acid rule calls amplifiers life about as often as life itself (precision 0.52).
 All of this is conditional on assumptions (a)-(h) and on the amplifier hand being L, q_np = 0.5, f priors shared by life and amplifier.
 CORRECTION (2026-10-16, after review; the lines above are left as written): the P9 line misreads the log.  "Largest drop: Iva
 only (0.136)" is the drop from keeping ONLY Iva, and the Iva-only "visible" set differs from the full panel's; on a fixed visible
 set no-Iva and Iva-only both give 0.524 (life_few.py).  Readings (1) and (5) follow from letting the amplifier act on all seven
 amino acids with life's hand and share prior; with an amplifier limited to 2-3 amino acids, life is confidently called in 25.8%
 and the >= 4 rule has precision 0.992 (life_few.py).  They are superseded by life_few.py.
"""
import os
import sys
import time

import numpy as np

import plume_signal as PS   # installs seedfix's rng correction on import; main() is guarded
import uncertainty as U
import sweep2 as S
from sklearn.ensemble import HistGradientBoostingClassifier

YR = 3.156e7
CHECK = []
CC = {"Ala": (24.91, 15863), "Asp": (23.34, 14725), "Glu": (22.56, 15675), "Ile": (27.24, 16913),
      "Leu": (24.30, 15750), "Phe": (19.82, 13440), "Val": (23.61, 15754)}   # Cohen & Chyba 2000 Table II
PROT = ["Ala", "Asp", "Phe", "Glu", "Leu", "Val", "Ile"]
IVA = 7
M = 60000   # samples per class, train and test


def check(p, ok, detail):
    CHECK.append((p, ok))
    print(f"  [{p}] {'HELD' if ok else 'FAILED'}: {detail}")


RATES_C = {aa: float(np.exp(a - b / 273.0)) for aa, (a, b) in CC.items()}
RATES_M = dict(RATES_C, Ala=S.K_RAC["Ala"], Asp=S.K_RAC["Asp"])


def load(seed=0):
    _, W = S.chains(seed)
    x, win, a, salt = W["Enceladus"]
    L = a[1]
    c, k2, Dh, t = (x[k] for k in ("concentration", "rate k2", "horizontal D_h", "time"))
    return k2 * c, t, (L / U.A) ** 2 / Dh <= t


def pool(raw, rates, tab):
    """Per-draw purity eta_i (N x 8, PROT order then Iva) and the Ala-gate mask (population A)."""
    k, t, mix = raw
    base = (k * t >= 10) & mix
    kk = np.maximum(k, 1e-300)
    E = np.column_stack([tab(rates[aa] / kk) * (base & (k > 2 * rates[aa])) for aa in PROT] + [tab(1e-12 / kk) * base])
    return E, base & (k > 2 * rates["Ala"])


def survival(rates, tau):
    return np.column_stack([PS.surv(rates[aa], tau) for aa in PROT] + [np.ones(len(tau))])


def gen(rng, m, kind, rates, E, Amask, qnp=0.5, hand=1.0, allpop=False):
    tau = 10 ** rng.uniform(3, 7, m)
    s = survival(rates, tau)
    f = 10 ** rng.uniform(-4, 0, (m, 8))
    sh = rng.random(m) < 0.5
    f[sh] = f[sh][:, :1]
    if kind == "A":
        idx = rng.integers(0, len(E), m) if allpop else rng.choice(np.flatnonzero(Amask), m)
        eta = E[idx].copy()
        eta[:, IVA] *= rng.random(m) < qnp
        sign = np.where(rng.random(m) < hand, 1.0, -1.0)[:, None]
        return f * eta * s * sign
    if kind == "B":
        eta = np.repeat(rng.uniform(0.97, 1.0, (m, 1)), 8, axis=1)
        eta[:, IVA] = 0.0
        return f * eta * s
    mete = rng.random(m) < 0.5
    ee = np.zeros((m, 8))
    ee[:, IVA] = rng.uniform(0.02, 0.20, m) * mete
    ee[:, :7] = rng.uniform(0, 0.05, (m, 7)) * (rng.random((m, 7)) < 0.5) * mete[:, None] * s[:, :7]
    return ee


def dataset(seed, rates, E, Amask, **kw):
    rng = np.random.default_rng(seed)
    out = {}
    for split in ("train", "test"):
        ee = np.vstack([gen(rng, M, kd, rates, E, Amask, **(kw if kd == "A" else {})) for kd in "ABC"])
        out[split] = (ee, np.repeat([0, 1, 2], M))
    return out


def Xof(obs, sig, cols, signed=True, summary=False):
    o = obs if signed else np.abs(obs)
    pc = [c for c in cols if c < 7]
    parts = []
    if pc:
        parts += [(o[:, pc] >= 3 * sig).sum(1, keepdims=True), o[:, pc].mean(1, keepdims=True)]
        if not summary:
            parts.append(o[:, pc])
    if IVA in cols:
        parts.append(o[:, [IVA]])
    return np.hstack(parts)


def fit(D, sig, cols, signed=True, summary=False):
    r = {}
    for split in ("train", "test"):
        ee, y = D[split]
        nz = np.random.default_rng(int(1000 * sig) + (1 if split == "train" else 2)).normal(0, sig, ee.shape)
        r[split] = (Xof(ee + nz, sig, cols, signed, summary), y, ee, ee + nz)
    clf = HistGradientBoostingClassifier(max_iter=120, learning_rate=0.1, max_leaf_nodes=15, random_state=0)
    clf.fit(r["train"][0], r["train"][1])
    P = clf.predict_proba(r["test"][0])
    y, ee, obs = r["test"][1], r["test"][2], r["test"][3]
    pred = P.argmax(1)
    vis = (np.abs(ee[:, cols]) >= 0.05).any(1)
    ab = (y < 2) & vis
    return dict(acc=(pred == y).mean(), accV=(pred == y)[ab].mean(), P=P, pred=pred, y=y, vis=vis, obs=obs, ee=ee)


def cols_prot(n):
    return list(range(n))


def cols_iva(n):
    return list(range(n - 1)) + [IVA]


def main(seed=0):
    t0 = time.time()
    tab = PS.eta_table()
    raw = load(seed)
    k, t, mix = raw
    print(f"life_test.py, N = {U.N}, seed 0, Enceladus (sweep2.chains index 0); M = {M} per class and split; "
          f"sigma 1/2/5%; equal priors.")

    print("\n(1) Rates (Cohen & Chyba 2000 Table II, 273 K) and survival")
    print("  " + "; ".join(f"{aa} {RATES_C[aa]:.3e}" for aa in sorted(RATES_C, key=lambda a: -RATES_C[a])) + " /s")
    print(f"  Asp/Ala: Cohen & Chyba {RATES_C['Asp'] / RATES_C['Ala']:.2f}; model {S.K_RAC['Asp'] / S.K_RAC['Ala']:.2f}")
    s5 = {aa: PS.surv(RATES_C[aa], 1e5) for aa in PROT}
    print("  s at tau 1e5 yr: " + "; ".join(f"{aa} {v:.3f}" for aa, v in s5.items()))

    EC, AC = pool(raw, RATES_C, tab)
    EM, AM = pool(raw, RATES_M, tab)
    AmodelAla = (k * t >= 10) & mix & (k > 2 * S.K_RAC["Ala"])
    etapos_M = (tab(S.K_RAC["Ala"] / np.maximum(k, 1e-300))[AmodelAla] > 0).mean()
    print(f"\n(2) Population: gate-passing share with the model's Ala rate {AmodelAla.mean():.4f} "
          f"(eta*_Ala > 0 in {etapos_M:.4f} of it); with Cohen & Chyba Ala {AC.mean():.4f}")
    zero = {aa: (EC[AC][:, i] == 0).mean() for i, aa in enumerate(PROT)}
    print("  Share of A draws (rate set C) with eta* = 0 (network too slow for that amino acid): " +
          "; ".join(f"{aa} {v:.3f}" for aa, v in zero.items()))

    rng = np.random.default_rng(20261016)
    DC = dataset(11, RATES_C, EC, AC)
    vis_share = {}
    ee_te, y_te = DC["test"]
    allcols = cols_iva(8)
    for ci, nm in enumerate("ABC"):
        vis_share[nm] = (np.abs(ee_te[y_te == ci][:, allcols]) >= 0.05).any(1).mean()
    print("\n(3) Visibility: share of samples with some panel |ee| >= 5% (n = 8 with Iva): " +
          "; ".join(f"{nm} {v:.3f}" for nm, v in vis_share.items()))

    print("\n(4) Unconditional 3-class accuracy (rows: panel; columns: sigma 1% / 2% / 5%), equal priors; "
          "and acc_V_AB (A vs B among visible)")
    grid, gridV = {}, {}
    for fam, mk, ns in (("protein only", cols_prot, (2, 4, 6, 7)), ("with Iva", cols_iva, (2, 4, 6, 8))):
        for n in ns:
            for sg in (0.01, 0.02, 0.05):
                r = fit(DC, sg, mk(n))
                grid[(fam, n, sg)], gridV[(fam, n, sg)] = r["acc"], r["accV"]
            print(f"  {fam:12s} n = {n}: acc " + " / ".join(f"{grid[(fam, n, sg)]:.3f}" for sg in (0.01, 0.02, 0.05)) +
                  " | acc_V_AB " + " / ".join(f"{gridV[(fam, n, sg)]:.3f}" for sg in (0.01, 0.02, 0.05)))
    best = max(grid.values())
    ok90 = [key for key, v in grid.items() if v >= 0.90]
    print(f"  Best cell accuracy {best:.3f}; cells >= 0.90: {ok90 if ok90 else 'none'}.  Standard error about 0.0012.")
    okV = [key for key, v in gridV.items() if v >= 0.90]
    print(f"  Cells with acc_V_AB >= 0.90: {okV if okV else 'none'}.")

    print("\n(5) Confusion matrices, n = 8 with Iva, sigma 2% (rows true A/B/C, columns called A/B/C)")
    r = fit(DC, 0.02, cols_iva(8))
    for lab, sel in (("all samples", np.ones(len(r["y"]), bool)), ("visible samples", r["vis"])):
        cm = np.array([[((r["y"] == i) & (r["pred"] == j) & sel).sum() for j in range(3)] for i in range(3)], float)
        cm /= np.maximum(cm.sum(1, keepdims=True), 1)
        print(f"  {lab}:")
        for i, nm in enumerate("ABC"):
            print(f"    true {nm}: " + "  ".join(f"{cm[i, j]:.3f}" for j in range(3)))
    main_r = r
    accV_main = r["accV"]
    P = r["P"]
    conf = P.max(1) >= 0.9
    cov = conf.mean()
    prec = (r["pred"] == r["y"])[conf].mean() if conf.any() else float("nan")
    confB = (conf & (r["pred"] == 1))[r["y"] == 1].mean()
    confA = (conf & (r["pred"] == 0))[r["y"] == 0].mean()
    confC = (conf & (r["pred"] == 2))[r["y"] == 2].mean()
    print(f"\n(6) Confident calls (posterior >= 0.9): coverage {cov:.3f}, precision {prec:.3f}; share of true A/B/C confidently "
          f"called correctly: {confA:.3f} / {confB:.3f} / {confC:.3f}")

    print("\n(7) Feature ablation, n = 8 with Iva, sigma 2%, main config (acc_V_AB / acc)")
    abl = {"full": (accV_main, main_r["acc"])}
    for nm, kw, cols in (("no Iva (n = 7 protein)", {}, cols_prot(7)), ("|ee| (no sign)", dict(signed=False), cols_iva(8)),
                         ("protein = count+mean only", dict(summary=True), cols_iva(8)), ("Iva only", {}, [IVA])):
        rr = fit(DC, 0.02, cols, **kw)
        abl[nm] = (rr["accV"], rr["acc"])
    for nm, (a, b) in abl.items():
        print(f"  {nm:28s} {a:.3f} / {b:.3f}")
    drop = {nm: abl["full"][0] - v[0] for nm, v in abl.items() if nm != "full"}
    top = max(drop, key=drop.get)
    print(f"  Largest drop in acc_V_AB: {top} ({drop[top]:.3f}).")

    print("\n(8) Davila's rule: life if >= 4 of 7 protein amino acids have observed ee >= 5% (sigma 2%, n = 7 protein, equal priors)")
    rp = fit(DC, 0.02, cols_prot(7))
    called = (rp["obs"][:, :7] >= 0.05).sum(1) >= 4
    prec_dav = (rp["y"][called] == 1).mean()
    shareA = (rp["y"][called] == 0).mean()
    rec_dav = called[rp["y"] == 1].mean()
    print(f"  precision of the life call {prec_dav:.3f}; share of life calls that are amplifiers {shareA:.3f}; recall of life {rec_dav:.3f}; "
          f"share of true amplifiers called life {called[rp['y'] == 0].mean():.3f}")

    print("\n(9) Sensitivities (n = 8 with Iva, sigma 2%): acc_V_AB / acc")
    sens = {}
    variants = {"rate set M (model Ala, Asp)": dict(rates=RATES_M, E=EM, A=AM),
                "q_np = 0": dict(rates=RATES_C, E=EC, A=AC, kw=dict(qnp=0.0)),
                "q_np = 1": dict(rates=RATES_C, E=EC, A=AC, kw=dict(qnp=1.0)),
                "A hand 50/50": dict(rates=RATES_C, E=EC, A=AC, kw=dict(hand=0.5)),
                "A over all draws": dict(rates=RATES_C, E=EC, A=AC, kw=dict(allpop=True)),
                "best: hand 50/50, q_np = 1": dict(rates=RATES_C, E=EC, A=AC, kw=dict(hand=0.5, qnp=1.0))}
    Dsens = {}
    for i, (nm, v) in enumerate(variants.items()):
        Dsens[nm] = dataset(21 + i, v["rates"], v["E"], v["A"], **v.get("kw", {}))
        rr = fit(Dsens[nm], 0.02, cols_iva(8))
        sens[nm] = (rr["accV"], rr["acc"])
        print(f"  {nm:28s} {rr['accV']:.3f} / {rr['acc']:.3f}")
    Dh = Dsens["A hand 50/50"]
    sg_full = fit(Dh, 0.02, cols_iva(8))["accV"]
    sg_abs = fit(Dh, 0.02, cols_iva(8), signed=False)["accV"]
    print(f"  Hand 50/50: signed {sg_full:.3f}, |ee| {sg_abs:.3f}")
    rb = fit(Dsens["best: hand 50/50, q_np = 1"], 0.01, cols_iva(8))
    print(f"  Best case at sigma 1%: acc {rb['acc']:.3f}, acc_V_AB {rb['accV']:.3f}")

    print("\n(10) Reading the numbers")
    print("  Unconditional accuracy is capped by invisibility: an amplifier or life with every ee below the noise looks racemic, "
          "which is also what plain abiotic looks like; see (3), (4), (5).")

    print("\nPredictions:")
    ratio = {aa: RATES_C[aa] / RATES_C["Ala"] for aa in RATES_C}
    ord_ok = [aa for aa in sorted(RATES_C, key=lambda a: -RATES_C[a])] == ["Phe", "Asp", "Ala", "Leu", "Val", "Ile", "Glu"]
    ref = dict(Ala=3.83e-15, Asp=5.15e-14, Glu=7.27e-16, Ile=8.41e-16, Leu=3.15e-15, Phe=1.69e-13, Val=1.56e-15)
    check("P1", all(abs(RATES_C[a] / ref[a] - 1) <= 0.01 for a in ref) and ord_ok and abs(ratio["Asp"] - 13.4) <= 0.2 and
          abs(s5["Ala"] - 0.976) <= 0.003 and abs(s5["Asp"] - 0.755) <= 0.005 and abs(s5["Phe"] - 0.484) <= 0.005,
          f"order ok {ord_ok}; Asp/Ala {ratio['Asp']:.2f}; s Ala/Asp/Phe {s5['Ala']:.3f}/{s5['Asp']:.3f}/{s5['Phe']:.3f}")
    check("P2", 0.786 <= AmodelAla.mean() <= 0.788 and 0.931 <= etapos_M <= 0.935 and 0.79 <= AC.mean() <= 0.84,
          f"A(model) {AmodelAla.mean():.4f}; eta>0 {etapos_M:.4f}; A(Cohen) {AC.mean():.4f}")
    check("P3", 0.04 <= zero["Ala"] <= 0.12 and 0.10 <= zero["Asp"] <= 0.30 and 0.18 <= zero["Phe"] <= 0.42 and
          0.0 <= zero["Glu"] <= 0.03 and zero["Phe"] > zero["Asp"] > zero["Ala"] > zero["Glu"],
          "eta*=0 shares " + ", ".join(f"{a} {zero[a]:.3f}" for a in ("Ala", "Asp", "Phe", "Glu")))
    check("P4", 0.45 <= vis_share["A"] <= 0.75 and 0.45 <= vis_share["B"] <= 0.75 and 0.35 <= vis_share["C"] <= 0.48,
          "visible " + ", ".join(f"{a} {v:.3f}" for a, v in vis_share.items()))
    a8 = grid[("with Iva", 8, 0.02)]
    check("P5", 0.50 <= a8 <= 0.70, f"n = 8 with Iva, sigma 2%: {a8:.3f}")
    check("P6", best < 0.75 and not ok90, f"max cell {best:.3f}; cells >= 0.90: {len(ok90)}")
    d_n = grid[("with Iva", 8, 0.02)] - grid[("with Iva", 2, 0.02)]
    d_s = grid[("with Iva", 8, 0.01)] - grid[("with Iva", 8, 0.05)]
    check("P7", 0.0 <= d_n <= 0.12 and 0.02 <= d_s <= 0.15, f"acc(n8)-acc(n2) {d_n:.3f}; acc(1%)-acc(5%) {d_s:.3f}")
    v7, v8 = gridV[("protein only", 7, 0.02)], gridV[("with Iva", 8, 0.02)]
    check("P8", 0.50 <= v7 <= 0.62 and 0.62 <= v8 <= 0.82, f"acc_V_AB protein n=7 {v7:.3f}; with Iva n=8 {v8:.3f}")
    check("P9", drop["no Iva (n = 7 protein)"] >= 0.08 and drop["no Iva (n = 7 protein)"] > drop["protein = count+mean only"] and
          drop["no Iva (n = 7 protein)"] > drop["|ee| (no sign)"] and drop["protein = count+mean only"] <= 0.04 and
          drop["|ee| (no sign)"] <= 0.04,
          "drops " + ", ".join(f"{n}: {d:.3f}" for n, d in drop.items()))
    check("P10", sg_full - sg_abs >= 0.05 and abs(drop["|ee| (no sign)"]) <= 0.02,
          f"hand 50/50 signed-|ee| {sg_full - sg_abs:.3f}; hand L {drop['|ee| (no sign)']:.3f}")
    s = sens
    check("P11", 0.50 <= s["q_np = 0"][0] <= 0.62 and 0.72 <= s["q_np = 1"][0] <= 0.90 and 0.70 <= s["A hand 50/50"][0] <= 0.88 and
          abs(s["rate set M (model Ala, Asp)"][0] - accV_main) <= 0.03,
          f"q0 {s['q_np = 0'][0]:.3f}; q1 {s['q_np = 1'][0]:.3f}; hand {s['A hand 50/50'][0]:.3f}; "
          f"M {s['rate set M (model Ala, Asp)'][0]:.3f} vs main {accV_main:.3f}")
    check("P12", 0.50 <= prec_dav <= 0.62, f"precision of Davila life call {prec_dav:.3f}")
    check("P13", confB <= 0.10 and 0.08 <= confA <= 0.30 and prec >= 0.90,
          f"B confident {confB:.3f}; A confident {confA:.3f}; precision {prec:.3f}")
    check("P14", rb["acc"] < 0.80, f"best case accuracy {rb['acc']:.3f}")
    print(f"\n{sum(ok for _, ok in CHECK)} of {len(CHECK)} held. {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
