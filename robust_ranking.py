"""Does the ocean-world ranking survive the newest measurements? (2026-10-06)

Two inputs in ocean_worlds.py were fixed guesses; this replaces each with a published range. The reported numbers
(ocean_worlds.log: Ganymede 28.3%, Europa 21.7%, Enceladus 0.0%) are NOT changed; these sit beside them.

(1) Europa's ice thickness. strat_sweep.py fixed 20 km (ocean-top 24 MPa, S* 3.0 g/kg). Juno's microwave radiometer
    gives a conductive shell of 29 +/- 10 km (Levin et al. 2025 Nat. Astron. 10, 84, 10.1038/s41550-025-02718-0).
    Here: h ~ Normal(29, 10) km clipped to 3-60; p = 920 kg/m^3 * 1.315 m/s^2 * h; S*(p) from the TEOS-10 table in
    verify_sstar.py (gsw 3.6.23), none above 28.5 MPa (~23 km); salinity log-uniform 1-100 g/kg as before.
(2) Enceladus' salinity. strat_sweep.py used 5-20 g/kg (Postberg et al. 2009). Kang et al. 2022 (Sci. Adv.,
    10.1126/sciadv.abm4665) collect other estimates: ice-shell geometry 7-30, geysers >= 17, silica nanoparticles
    < 40 (and a fast-mixing ocean, which Zeng & Jansen 2021 say a stratified top would prevent). Each prior is run with
    S* = 20 (used) and 22.5 (thin-shell top).

PREDICTIONS (written before running):
 (P1) P(h < 23 km) ~ 0.27; Europa pinned ~7-8% of draws; Europa ~25-27%, within Ganymede's own depth spread (22-33%)
      of its 28.3% -> the honest statement is a tie, not "Ganymede first".
 (P2) Enceladus: 5-20 -> 0%; 7-30 -> ~6%; 17-40 -> ~15-18%. Below Ganymede and Europa under every prior.

RESULTS (robust_ranking.log, 2026-10-06): P1 held (ice < 23 km 27%; final 26.1%, pinned 9.1% vs ~7-8% predicted).
 P2 held: 0.0 / 4.2-6.0 / 14.7-17.7%. Europa no longer sits clearly below Ganymede: 26% vs 28%, inside Ganymede's own
 depth spread, so it is a tie. Enceladus stays last under every published salinity, but it is only 0% if its
 salinity is below ~20 g/kg.
"""
import numpy as np

import uncertainty as U

# ocean-top pressure MPa -> S* g/kg (verify_sstar.py output); 0 above 28.5 MPa = no lid possible
SSTAR_P = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 27, 28, 28.5]
SSTAR_S = [23.85, 22.1, 20.35, 18.6, 16.85, 15.1, 13.36, 11.64, 9.93, 8.24, 6.57, 4.94, 3.34, 1.77, 1.0, 0.24, 0]


def sweep(name, seed=0):
    U.RNG = np.random.default_rng(seed)
    return U.sweep(*U.WORLDS[name])[1]


def europa(seed=0):
    win = sweep("Europa", seed)
    h = np.clip(U.RNG.normal(29e3, 10e3, U.N), 3e3, 60e3)
    s_star = np.interp(920 * 1.315 * h / 1e6, SSTAR_P, SSTAR_S, right=0)
    pinned = U.logu(1, 100) < s_star
    return (h < 23e3).mean(), pinned.mean(), win.mean(), (win & ~pinned).mean()


def enceladus(s_rng, s_star, seed=0):
    win = sweep("Enceladus", seed)
    pinned = U.logu(*s_rng) < s_star
    return pinned.mean(), (win & ~pinned).mean()


def main():
    thin, p, raw, eu = europa()
    print(f"Europa, Juno ice 29+/-10 km: ice < 23 km {thin:.0%}, pinned {p:.1%}, no lid {raw:.1%}, final {eu:.1%}"
          f"   (reported, 20 km ice: 21.7%)")
    print("Ganymede (no lid possible): 28.3% reported, 22-33% over 50-300 km ocean")
    enc = {}
    for s_rng in ((5, 20), (7, 30), (17, 40)):
        for s_star in (20.0, 22.5):
            q, enc[s_rng, s_star] = enceladus(s_rng, s_star)
            print(f"Enceladus salinity {s_rng[0]}-{s_rng[1]} g/kg, S* {s_star}: pinned {q:.0%}, final {enc[s_rng, s_star]:.1%}")
    assert max(enc.values()) < eu < 0.33, "ranking claim in RESULTS below drifted"
    return eu, enc


if __name__ == "__main__":
    main()
