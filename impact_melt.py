"""Titan impact-melt pools: could the weak force pick the hand in a crater's liquid water?

Lunine et al. 1999 proposed Titan's impact-melt pools as a natural lab for this question but computed
nothing (limits_lit.md). This script runs the project's own criterion on them, the same way uncertainty.py
does for the oceans. Draw 200,000 random input sets per pool and count "physics wins" only if ALL hold:
  finishes      k2*c*tau >= 10                    (ocean.py time rule; tau = how long the pool stays liquid)
  one pool      (L/A)^2/D <= tau, A = 5.5         (coarsen3d.py healing rule; L = sqrt(V/H) pool width)
  strong enough F * Delta >= 2                    (ocean.py x real-chemistry factor F)
Layering is skipped: a freezing pool is cooled from the top and should convect (ASSUMED). The pool is
taken as isolated; an inherited excess e0 would compete with the weak force (inherit.py), so the script
also prints the e0 that would override it.

Pools (melt volume, depth, lifetime):
  Small crater (1 km impactor; Lunine 1999: melt lasts thousands of years). Volume SCALED from Selk by
      impactor mass (1/64 of 190-400 km^3 -> 3-6 km^3); lifetime 1e3-1e4 yr. ASSUMED scaling.
  Selk (~80 km crater, 4 km impactor). Wakita et al. 2023 (PSJ, doi:10.3847/psj/acbe40): melt 190-220 km^3
      (complete) to 360-400 km^3 (incipient), torus a few km deep, freezes in up to tens of kyr; estimates
      they cite run from a few hundred yr (thin sheet) to 5,000 and 90,000 yr. Kalousova et al. 2024
      (JGR Planets, doi:10.1029/2023je008107): small pools freeze in <~25 kyr. -> 190-400 km^3, 1-5 km deep,
      3e2-9e4 yr.
  Menrva (~425 km basin). Crosta et al. 2021 (Icarus, doi:10.1016/j.icarus.2021.114679): melt / partial melt
      to ~65 km depth and ~60 km from the centre, and ocean water rises within 10 km of the surface. Melt
      volume 1e3-1e5 km^3, 20-65 km deep, 1e4-1e6 yr are ASSUMED from that geometry. It may join the ocean,
      in which case it is not an isolated pool.
  Lake (1 km^2 x 10 m), reference: the README Part 4 lake, run through the same sweep.
Shared ranges: g 4e-19..4e-16 (uncertainty.py; colder ammonia-water raises g by <=1.5x, ignored);
k2 1e-6..1 /M/s; F U(0.3, 0.7); concentration 1e-8..1e-2 M (ASSUMED: no measurement of a melt pool's
organic load exists; Madan & Pearce 2025, doi:10.3847/psj/ae1c18, argue amino acids form there);
convective D 1e-4..1e1 m^2/s (ASSUMED).

PREDICTIONS (written 2026-10-11, before the first run):
  P1. Selk wins in 35-65% of input space.
  P2. Small crater wins less often than Selk, in 10-40%.
  P3. Menrva wins in 60-95%.
  P4. Ranking: Menrva > Selk > small crater > lake.
  P5. Healing never binds: "heals" >= 95% for every pool.
  P6. At g = 1e-17, k2 = 1e-3, F = 0.5 and mid-range volume/lifetime, Selk needs c between 1e-5 and 1e-3 M.
  P7. Concentration is the most decisive input for Selk (largest low-half/high-half gap).
  P8. The inherited excess that would override the weak force in Selk (c = 1e-4 M, k2 = 1e-3, 1e4 yr)
      is below 1e-12, i.e. far below what meteorites carry (0.185): any seeded excess wins.

RESULTS (impact_melt.log, 2026-10-11):
  P1 PASS  Selk wins in 51%.
  P2 PASS  small crater 33%, below Selk.
  P3 PASS  Menrva 74%.
  P4 PASS  Menrva 74% > Selk 51% > small crater 33% > lake 11%.
  P5 PASS  heals 100% / 98% / 100% / 100%.
  P6 FAIL  Selk needs c >= 7.7e-6 M, just under the predicted 1e-5..1e-3 (bigger volume than guessed helps).
  P7 PASS  concentration decides most for every pool (Selk: 20% low half vs 83% high half).
  P8 PASS  override excess 2e-16 for Selk (1e-15 small crater, 6e-17 Menrva, 2e-14 lake).
  Reading: a Selk-size pool is a coin flip that leans to the weak force (needs ~10 uM organics); the lake
  reference needs ~8 mM, matching README Part 4 ("needs >~ mM"). The whole result hangs on two things
  nobody has measured: the pool's organic concentration, and whether it inherits any excess at all
  (anything above ~1e-16, e.g. from the impactor, overrides the weak force).
"""
import numpy as np
from ocean import delta_at, YR
from inherit import e_need

A = 5.5
N = 200_000
KM3 = 1e9

POOLS = {  # volume m^3, depth m, lifetime yr
    "Small crater (1 km impactor)": ((3 * KM3, 6 * KM3), (5e2, 2e3), (1e3, 1e4)),
    "Selk melt pool": ((190 * KM3, 400 * KM3), (1e3, 5e3), (3e2, 9e4)),
    "Menrva melt zone": ((1e3 * KM3, 1e5 * KM3), (2e4, 6.5e4), (1e4, 1e6)),
    "Lake 1 km^2 x 10 m (reference)": ((1e7, 1e7), (10, 10), (1e3, 1e4)),
}


def sweep(V_rng, H_rng, t_rng, seed=0):
    rng = np.random.default_rng(seed)
    logu = lambda lo, hi: 10 ** rng.uniform(np.log10(lo), np.log10(hi), N)
    x = {"bias g": logu(4e-19, 4e-16), "rate k2": logu(1e-6, 1.0), "concentration": logu(1e-8, 1e-2),
         "mixing D": logu(1e-4, 1e1), "lifetime": logu(*t_rng) * YR, "volume": logu(*V_rng),
         "depth": logu(*H_rng), "chem. factor F": rng.uniform(0.3, 0.7, N)}
    t = x["lifetime"]
    L = np.sqrt(x["volume"] / x["depth"])
    finishes = x["rate k2"] * x["concentration"] * t >= 10
    heals = (L / A) ** 2 / x["mixing D"] <= t
    strong = x["chem. factor F"] * delta_at(x["concentration"], x["bias g"], x["rate k2"], x["volume"], t) >= 2
    win = finishes & heals & strong
    return x, win, {"finishes": finishes.mean(), "heals": heals.mean(), "strong enough": strong.mean()}


def c_min(V, tau, g=1e-17, k2=1e-3, F=0.5):
    """Smallest concentration (M) meeting both the time rule and F*Delta >= 2 (bisection in log c)."""
    lo, hi = -12.0, 1.0
    ok = lambda lc: k2 * 10**lc * tau >= 10 and F * delta_at(10**lc, g, k2, V, tau) >= 2
    assert ok(hi), "no concentration up to 10 M works"
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return 10**hi


def main():
    rates = {}
    for name, (V_rng, H_rng, t_rng) in POOLS.items():
        x, win, parts = sweep(V_rng, H_rng, t_rng)
        rates[name] = win.mean()
        V, tau = np.sqrt(V_rng[0] * V_rng[1]), np.sqrt(t_rng[0] * t_rng[1]) * YR
        print(f"\n{name}: physics wins in {win.mean():.0%} of plausible input space")
        print("  each condition alone: " + ", ".join(f"{k} {v:.0%}" for k, v in parts.items()))
        print(f"  needs c >= {c_min(V, tau):.1e} M at g=1e-17, k2=1e-3, F=0.5 "
              f"(V={V / KM3:.3g} km^3, {tau / YR:.0f} yr)")
        print(f"  inherited excess that overrides the weak force (c=1e-4 M): e0 >= {e_need(1e-4, V, 1e-3, tau):.0e}")
        effects = []
        for k, v in x.items():
            if np.ptp(v) == 0:
                continue
            lo = v <= np.median(v)
            effects.append((win[~lo].mean() - win[lo].mean(), k, win[lo].mean(), win[~lo].mean()))
        effects.sort(key=lambda e: -abs(e[0]))
        for d, k, plo, phi in effects[:4]:
            print(f"  {k:15} low half {plo:5.0%}  high half {phi:5.0%}   (effect {d:+.0%})")
    # sanity: c_min round-trips, a best-case input set wins, a hopeless one loses
    V, tau = 3e11, 1e4 * YR
    c = c_min(V, tau)
    assert abs(0.5 * delta_at(c, 1e-17, 1e-3, V, tau) - 2) < 1e-6 or abs(1e-3 * c * tau - 10) < 1e-6
    assert 0.3 * delta_at(1e-2, 4e-16, 1.0, 4e11, 9e4 * YR) > 2
    assert 0.7 * delta_at(1e-8, 4e-19, 1e-6, 1e7, 1e3 * YR) < 2


if __name__ == "__main__":
    main()
