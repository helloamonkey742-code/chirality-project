"""Arm G at measured rates: r from Essex et al. 2026, then rna_copying.py at that r.

r = (rate an L nucleotide adds to a D-RNA primer) / (rate a D nucleotide adds), on a D template in a
50/50 pool (rna_copying.py: p = m_R / (m_R + r m_W)). The racemic state can be lost only if r > 1.

Source: Essex, Jia, Zhang & Szostak 2026, ChemRxiv, doi:10.26434/chemrxiv.15005970/v1 (main text + SI).
Every run there: 200 mM HEPES pH 8.0, 100 mM MgCl2, D-RNA primer 1.5 uM, D-RNA template 2.5 uM.

Three routes to r:
  (1) monomer only (Fig. 1, excess free 2-AI suppresses bridged dinucleotides): k_L / k_D.
  (2) bridged pathway, from the kinetic constants. Specificity constants k/K = kobs,max/Km of the four
      A*C species (Fig. 4; +1 nucleotide is added, +2 leaves; a downstream helper oligo is bound),
      weighted by the measured formation rates of hetero- and homochiral bridged dinucleotides from
      a 1:1 *C_D:*C_L pool (Fig. S3). Hydrolysis is the same for both (Fig. S6), so steady-state
      amounts follow formation. One binding site, competitive Michaelis-Menten:
      v_i = (k/K)_i [S_i] / (1 + sum_j [S_j]/K_j), so the L:D addition ratio depends on k/K only,
      at any saturation.
  (3) direct competition (Fig. 9): *C from 5% to 90% ee (D), 20 mM total, primer PD1 on template TD2,
      24 h, % reduction of primer extension relative to pure D. TD2 templates G, G, A, A after the
      primer, so +1 and +2 are matched for C and +3, +4 are not. All primer is used up (Fig. S8 gel),
      so the reduction counts primers stopped by an L. The paper does not say which bands it counted:
        reading A, products past +1 (stopped only by an L at +1):        reduction = 1 - p
        reading B, products past +2 (need D at +1 and at +2):             reduction = 1 - p^2
      with p = x_D / (x_D + r x_L), x = monomer fraction. Both assume no D-ended primer stalls.

Output: arm_g_rates.log.  Run: uv run ... python arm_g_rates.py
"""
import numpy as np
from scipy.optimize import minimize_scalar

import rna_copying as rc

# ------------------------------------------------------------------ CITED  (value, +-)
MONO = {"D": (2.37e-2, 0.08e-2), "L": (1.4e-3, 0.4e-3)}   # p.2, Fig. 1: /h, 20 mM *C, 100 mM 2-AI, SD n=3
CC = {"DD": (19.7, 0.3, 1.07, 0.06, 18.4),                # p.2, Fig. 2: kobs,max /h, +-, Km mM, +-,
      "LL": (0.37, 0.06, 2.1, 0.1, 0.18)}                 #   reported k/K /h/mM (C*C homochiral only)
AC = {"DD": (8.6, 0.2, 2.3, 0.5, 3.69),                   # Fig. 4a,c: key = hand at +1, +2; SE n=2
      "DL": (0.53, 0.09, 11.0, 2.0, 4.8e-2),
      "LD": (1.63, 0.01, 1.34, 0.07, 1.22),
      "LL": (0.26, 0.04, 4.0, 2.0, 6.7e-2)}
FORM = {"all": (2.4e-3, 0.1e-3), "hetero": (1.3e-3, 0.1e-3), "homo": (1.1e-3, 0.1e-3)}  # p.3, Fig. S3:
#   /h/mM, 6 + 6 mM *C; homo = D/D + L/L together (enantiomers share one 31P peak)
HYD = {"D/L": (8.78e-2, 0.54e-2), "D/D": (8.62e-2, 0.64e-2)}  # Fig. S6: A*C hydrolysis /h
TRIMER = {"D": (14.0, 1.0), "L": (7.0, 0.2)}              # p.4, Fig. 5: /h, *C + 0.5 mM D-trimer *GAC
FIG9 = [(5, 51.9, 0.6), (10, 50.0, 1.0), (20, 42.1, 0.7), (40, 32.0, 3.0),  # Fig. 9c: %ee (D),
        (60, 20.9, 0.5), (80, 11.9, 0.7), (90, 6.1, 0.4)]                   #   % reduction, SD n=3


def eff(d):
    return {k: v[0] / v[2] for k, v in d.items()}


def r_bridged(e, het, hom, x=0.5):
    """r from the bridged pathway at D fraction x. Amounts ~ formation: D/D ~ hom (2x)^2 / 2,
    each hetero orientation ~ het 4x(1-x) / 2 (rates measured at x = 1/2)."""
    dd, ll, hx = hom * 2 * x * x, hom * 2 * (1 - x) ** 2, het * 2 * x * (1 - x)
    L = hx * e["LD"] + ll * e["LL"]
    D = dd * e["DD"] + hx * e["DL"]
    return L / D * x / (1 - x)


def fig9_r(n):
    """Per-point r and weighted least-squares r for reading A (n = 1) or B (n = 2)."""
    ee, red, sd = (np.array(c, float) / 100 for c in zip(*FIG9))
    xd, xl = (1 + ee) / 2, (1 - ee) / 2
    per = ((1 - red) ** (-1 / n) - 1) * xd / xl
    model = lambda r: 1 - (xd / (xd + r * xl)) ** n
    fit = minimize_scalar(lambda r: (((model(r) - red) / sd) ** 2).sum(), bounds=(0.01, 10), method="bounded")
    return per, fit.x, fit.fun / (len(FIG9) - 1)


def main():
    rng = np.random.default_rng(0)
    print("Arm G: r from Essex, Jia, Zhang & Szostak 2026 (ChemRxiv), then rna_copying.py at that r")
    print("r = L addition rate / D addition rate on a D primer-template, 50/50 pool; r > 1 needed to lose 50/50")

    print("\n(1) CITED (main text, Figs. 1-5, 9, S3, S6)")
    rm = MONO["L"][0] / MONO["D"][0]
    print(f"  monomer only, 20 mM *C + 100 mM 2-AI:  k_D = {MONO['D'][0]:.2e}, k_L = {MONO['L'][0]:.1e} /h  -> r = {rm:.3f}"
          f" (1/{1 / rm:.0f}; paper ~17-fold)")
    ec, ea = eff(CC), eff(AC)
    print(f"  C*C homochiral: k/K D/D {ec['DD']:.1f}, L/L {ec['LL']:.2f} /h/mM (ratio 1/{ec['DD'] / ec['LL']:.0f}; paper ~102)")
    print("  A*C (+1,+2):    kobs,max /h   Km mM    k/K /h/mM (reported)")
    for k, (km, _, Km, _, rep) in AC.items():
        print(f"     {k[0]}/{k[1]}          {km:6.2f}     {Km:5.2f}    {km / Km:.3f} ({rep})")
    het, hom = FORM["hetero"][0], FORM["homo"][0]
    print(f"  formation from 1:1 *C: hetero {het:.1e}, homo (D/D + L/L) {hom:.1e} /h/mM -> hetero share "
          f"{het / (het + hom):.2f} (0.50 if unbiased)")
    print(f"  hydrolysis A*C: D/L {HYD['D/L'][0]:.3f}, D/D {HYD['D/D'][0]:.3f} /h (equal within error)")
    rt = TRIMER["L"][0] / TRIMER["D"][0]
    print(f"  monomer bridged to a D-trimer: k_D {TRIMER['D'][0]:.0f}, k_L {TRIMER['L'][0]:.1f} /h -> r = {rt:.2f}")

    print("\n(2) bridged pathway from k/K and formation rates")
    rb = r_bridged(ea, het, hom)
    draws = []
    for _ in range(20000):
        e = {k: max(rng.normal(km, s1), 1e-3) / max(rng.normal(Km, s2), 0.05) for k, (km, s1, Km, s2, _) in AC.items()}
        draws.append(r_bridged(e, rng.normal(het, FORM["hetero"][1]), rng.normal(hom, FORM["homo"][1])))
    lo, med, hi = np.percentile(draws, [2.5, 50, 97.5])
    print(f"  r = (het k/K[L,D] + hom k/K[L,L]) / (hom k/K[D,D] + het k/K[D,L]) = {rb:.3f}")
    print(f"  95% interval {lo:.2f}-{hi:.2f} (SE draws, 2e4); P(r > 1) = {np.mean(np.array(draws) > 1):.4f}")
    print(f"  L additions come {100 * het * ea['LD'] / (het * ea['LD'] + hom * ea['LL']):.0f}% through the L(+1)/D(+2) "
          "hetero species: the heterochiral bridge is the main L-adding channel")
    xs = (0.525, 0.6, 0.8, 0.95)
    print("  same r at other pool compositions (D fraction): " + ", ".join(f"{x}: {r_bridged(ea, het, hom, x):.2f}" for x in xs))
    print("  -> r is nearly composition-independent, as Fig. 9's linear plot requires")

    print("\n(3) direct competition, Fig. 9 (PD1/TD2, 24 h)")
    print("  %ee   % red   r (A: past +1)   r (B: past +2)")
    pa, ra, ca = fig9_r(1)
    pb, rb2, cb = fig9_r(2)
    for (ee, red, _), a, b in zip(FIG9, pa, pb):
        print(f"  {ee:4d}   {red:5.1f}        {a:5.2f}            {b:5.2f}")
    print(f"  weighted fit: reading A r = {ra:.2f} (chi2/dof {ca:.1f}), reading B r = {rb2:.2f} (chi2/dof {cb:.1f})")
    print(f"  reading A needs a C_L*C_D bridge (L at +1) to add as fast as C_D*C_D: k/K ratio "
          f"{ra * hom / het:.2f}; measured for A*C: {ea['LD'] / ea['DD']:.2f}")
    print("  reading B agrees with route (2); reading A is the more direct test of 'an L at +1 stops extension'")

    print("\n(4) rna_copying.py at these r (n = 3, c = 0.5, w = 0.02 as LO/HI; feed 0.5-100 mM scanned)")
    rmax = (0.5 * 0.98 / 0.02) ** (1 / 3) - 1
    print(f"  window for losing 50/50: 1 < r < {rmax:.2f}")
    print(f"  {'source':38} {'r':>6} {'lam @20 mM':>11} {'m* (mM)':>8} {'+20% feed, 5000 cycles':>23}")
    cases = [("monomer only (Fig. 1)", rm), ("bridged, k/K (Fig. 4, S3)", rb), ("  95% upper end", hi),
             ("Fig. 9 reading B", rb2), ("monomer + D-trimer (Fig. 5)", rt), ("Fig. 9 reading A", ra),
             ("  lowest per-point (reading A)", pa.min()), ("  highest per-point (reading A)", pa.max())]
    out = {}
    for name, r in cases:
        k = dict(rc.LO, r=r)
        lam = rc.mode(rc.fixed_point(20.0, **k), 20.0, k)[0].real.max()
        ms = rc.threshold(k)[0]
        z0 = np.array([rc.T_IN / 2, rc.T_IN / 2, 20.0 * 1.2 / 2, 20.0 * 0.8 / 2])
        ee = rc.run(z0, 20.0, k, 5000, e=0.2)[1][-1]
        out[name] = ms
        print(f"  {name:38} {r:6.3f} {lam:11.5f} {'none' if ms is None else f'{ms:8.2f}':>8} {ee:23.3f}")

    print("\nANSWER")
    print(f"  Routes (1), (2) and Fig. 9 reading B give r = {rm:.2f}, {rb:.2f} ({lo:.2f}-{hi:.2f}), {rb2:.2f}: below 1,")
    print("  so 50/50 is stable at every feed and arm G follows a seed to a bounded plateau.")
    print(f"  Fig. 9 reading A gives r = {ra:.2f} (per point {pa.min():.2f}-{pa.max():.2f}): inside the window, with")
    print(f"  m* = {out['Fig. 9 reading A']:.1f} mM feed, i.e. an R9 candidate.")
    print("  The paper's own data do not decide between them: it does not say which gel bands Fig. 9 counted.")
    print("  Decisive and cheap: the hand of the +1 nucleotide in Fig. 9's racemic lane (EXPERIMENT.md, Arm G step 0).")

    # transcription checks against the paper's own derived numbers
    assert abs(1 / rm - 17) < 1 and abs(ec["DD"] / ec["LL"] - 102) < 3, "monomer / C*C ratios"
    assert all(abs(km / Km / rep - 1) < 0.04 for km, _, Km, _, rep in AC.values()), "A*C k/K"
    assert all(abs(km / Km / rep - 1) < 0.04 for km, _, Km, _, rep in CC.values()), "C*C k/K"
    assert abs(het + hom - FORM["all"][0]) < 1e-5, "formation partition"
    assert rb < 1 < ra < rmax and out["bridged, k/K (Fig. 4, S3)"] is None and out["Fig. 9 reading A"] is not None
    print("\nself-check passed")


if __name__ == "__main__":
    main()
