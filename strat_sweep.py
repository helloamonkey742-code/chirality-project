"""Feed the two_zone.py pinning result into the uncertainty.py estimate (2026-10-06).

two_zone.py: a slow-mixing (stratified) layer ABOVE faster mixing pins a sinking wall for good once the mixing ratio
exceeds ~1.03, so the ocean keeps two hands. Real stratification jumps D_z by far more than 3%, so any input draw
where the ocean is "stratified over convecting" cannot end with one ocean-wide hand.

When does an icy-moon ocean have a stratified top? (literature, 2026-10-06 lookup)
  - Thermal expansion reverses (cold water lighter) at the freezing point below a salinity S*. Then the ocean stacks
    a stratified upper layer over a convecting interior (Zeng & Jansen 2021 PSJ 10.3847/psj/ac1114; 2024 PSJ
    10.3847/psj/ad0cba). S* from TEOS-10 (gsw 3.6.23, seawater composition):
      Enceladus top of ocean (~2 MPa): 22.5 g/kg (18.0 at the floor);  Europa top (~25 MPa, 20 km ice): 3.0 g/kg.
  - Enceladus salinity 5-20 g/kg (Postberg et al. 2009 Nature 10.1038/nature08046, 0.5-2% salt) -> all below S*.
  - Enceladus poles: no parameter choice gives an unstratified ocean top to bottom under a steady ice shell
    (Ames, Ferreira & Czaja 2025 CEE 10.1038/s43247-025-02036-3), independent of salinity.
  - Europa: weak, salt-dominated stratification but strong convection (Ashkenazy & Tziperman 2021 Nat. Commun.);
    flat ice shell (Kang 2022 ApJ), so the slope-driven stratified layer of Zhang, Kang & Marshall 2026 is weak.
    Europa salinity prior ASSUMED log-uniform 1-100 g/kg (no single published range found in this lookup).

Model: each draw also gets a salinity S; p_strat = P(S < S*). Strict rule: pinned draws lose. Win = old win & ~pinned.
Salinity is drawn independently of every other input (nothing in the lookup ties it to D_z, D_h, time or chemistry).

PREDICTIONS (written before running, 2026-10-06):
 (P1) Independence -> new win = old win x (1 - p_strat), within sampling noise.
 (P2) Enceladus, S log-uniform 5-20 g/kg, S* = 20 (between 18 floor and 22.5 top): p_strat ~1 -> win ~0-2%.
      With Ames 2025 (poles always stratified) the Enceladus answer is "no single hand at the poles" whatever S is.
 (P3) Europa, S log-uniform 1-100, S* = 3: p_strat = log(3)/log(100) = 0.24 -> win ~29% x 0.76 = 22%.
 (P4) uncertainty.py's own numbers (22% / 29%) are unchanged; this script reports beside them, it does not replace them.

RESULTS (strat_sweep.log, 2026-10-06): P1-P4 all held.
 Enceladus 21.8% -> 0.0% (p_strat 100%). Survives only if salinity can exceed S*: 1.6% (S*=18), 7.4-7.7% (S up to 40-100).
 Europa 28.6% -> 21.7% (p_strat 24%). Range 10-29% depending on S* (ice thickness) and the salinity floor.
"""
import numpy as np

import uncertainty as U

# name: (salinity range g/kg, S* g/kg)
SALT = {"Enceladus": ((5, 20), 20.0), "Europa": ((1, 100), 3.0)}


def main():
    out = {}
    for name, ((s_lo, s_hi), s_star) in SALT.items():
        x, win, _ = U.sweep(*U.WORLDS[name])
        S = U.logu(s_lo, s_hi)
        pinned = S < s_star
        new = win & ~pinned
        p = pinned.mean()
        out[name] = (win.mean(), p, new.mean())
        print(f"{name}: old win {win.mean():.1%}  P(stratified top, S<{s_star:g} g/kg) {p:.1%}  "
              f"new win {new.mean():.1%}  (old x (1-p) = {win.mean() * (1 - p):.1%})")
        assert abs(win.mean() - U.README_WIN_PCT[name]) < U.WIN_TOL, "uncertainty.py headline drifted"
        assert abs(new.mean() - win.mean() * (1 - p)) < 0.005, "salinity not independent of the win"
    # sensitivity: how the answer moves with S* (cutoff) and the salinity prior
    print("\nsensitivity (new win):")
    for name, rng, s_star in (("Enceladus", (5, 20), 18.0), ("Enceladus", (5, 40), 20.0), ("Enceladus", (1, 100), 20.0),
                              ("Europa", (1, 100), 1.5), ("Europa", (1, 100), 19.5), ("Europa", (3, 100), 3.0)):
        x, win, _ = U.sweep(*U.WORLDS[name])
        S = U.logu(*rng)
        print(f"  {name:9} S {rng[0]:>3}-{rng[1]:<3} g/kg, S* {s_star:>4}: p_strat {(S < s_star).mean():5.1%}  "
              f"win {(win & (S >= s_star)).mean():5.1%}")
    return out


if __name__ == "__main__":
    main()
