"""Every known or suspected ocean world in the uncertainty.py sweep, with strat_sweep.py's stratified-lid check (2026-10-06).

Same model and same chemistry/mixing ranges as Europa (ASSUMED: nothing in the lookup sets these moons apart there).
Per world: radius R, ocean thickness H (with a sensitivity range), ocean age, and the stratified-lid check:
the ocean top is pinned (strat_sweep.py, two_zone.py) when salinity S < S*, the salinity below which cold water at the
freezing point floats. S* from TEOS-10 (gsw 3.6.23) at the ocean-top pressure rho_ice * g * ice thickness:
  Mimas 1.5 MPa -> 22.6 g/kg; Miranda 2.2 -> 21.9; Dione 21 -> 5.6; Umbriel 21 -> 5.6; Ariel 25 -> 2.7.
  Above 28 MPa no S* exists (cold fresh water already sinks): Ganymede, Callisto, Titan, Triton, Pluto, Titania, Oberon.
  TEOS-10 is seawater; ammonia (likely in Titan, Triton, Pluto, Uranian moons) is not modelled.
Salinity prior where unmeasured: log-uniform 1-100 g/kg (ASSUMED, as Europa in strat_sweep.py).
Enceladus and Europa use exactly strat_sweep.py's inputs. The table in the code is the v2 (literature-checked) one below.

Interior inputs (2026-10-06 lookup; ASSUMED where marked):
  Titan     R 2575 km; ocean under ~100 km ice, H 200 km ASSUMED (50-400) (Durante et al. 2019 Icarus
            10.1016/j.icarus.2019.03.003). Ocean confirmed by gravity (tidal Love number).
  Mimas     R 198 km; ocean 20-30 km down, < 25 Myr old (Lainey et al. 2024 Nature 10.1038/s41586-023-06975-9).
            H 30 km ASSUMED (10-50); age 2-25 Myr.
  Dione     R 561 km; ~99 km crust over ~65 km ocean (Beuthe et al. 2016 GRL; Zannoni et al. 2020 Icarus
            10.1016/j.icarus.2020.113713). H 65 km (30-100).
  Titania, Oberon  R 789 / 761 km; thin residual oceans possible under ~190 km ice (Castillo-Rogez et al. 2023 JGR
            Planets). H 30 km ASSUMED (10-50). Unconfirmed.
  Ariel, Umbriel, Miranda  R 579 / 585 / 236 km; thick oceans possible by tidal resonance (Tyler 2025 GRL
            10.1029/2025gl114911). Ariel/Umbriel H 50 km ASSUMED (10-150) under 100 km ice ASSUMED; Miranda H 100 km
            ASSUMED (30-150) under 30 km ice, possibly a past ocean, age 0.1-1 Gyr ASSUMED. Unconfirmed.
  Triton    R 1353 km; likely ocean (Hansen et al. 2021 PSJ 10.3847/psj/abffd2). H 100 km ASSUMED (30-200), ~150 km ice.
  Pluto     (dwarf planet) R 1188 km; ocean > 100 km argued from Sputnik Planitia (Nimmo et al. 2016 Nature), disputed
            by an oceanless impact model (Ballantyne et al. 2024 Nat. Astron. 10.1038/s41550-024-02248-1). H 100 km (50-200).
  Not included (no liquid water ocean, so the model does not apply): Earth's Moon, Phobos, Deimos, Io, Rhea and the
  remaining small moons; Charon (frozen); Ceres (local brine, not a global ocean).

PREDICTIONS (written before running, 2026-10-06):
 (P1) No-lid worlds land near Europa's pre-lid 29%: thinner oceans higher (layers merge faster, merge time H*R/(2 D_z)),
      thicker lower. Titan (200 km) ~20-25%; Titania/Oberon (30 km, small R) ~35-45%; Triton/Pluto ~25-30%.
 (P2) Mimas: its ocean is < 25 Myr old, so layers often cannot merge in time (needs D_z >~ 1e-5): pre-lid ~10-20%,
      and S* 22.6 pins it ~70% of the time -> ~3-6%.
 (P3) Miranda (S* 21.9) -> pinned ~67% -> ~10%. Dione/Umbriel (S* 5.6) pinned ~37% -> ~20%. Ariel (S* 2.7) ~22% -> ~25%.
 (P4) The best confirmed ocean (gravity/magnetic evidence) is Ganymede or Titan, both ~25-29%; the highest numbers
      overall belong to the least certain oceans (Uranian moons).

RESULTS (ocean_worlds.log, 2026-10-06): P1-P4 held (Pluto 30.2%, just over P1's 25-30).
 final: Oberon 36.5, Titania 36.4, Callisto 35.7, Pluto 30.2, Triton 29.7, Ganymede 28.3, Ariel 27.2, Titan 24.3,
 Umbriel 21.9, Europa 21.7, Dione 21.0, Miranda 9.4, Mimas 5.2, Enceladus 0.0 (%).
 Pattern: two things set the number. (1) A salt lid (only small/low-pressure moons can have one) cuts or kills it.
 (2) Without a lid, thinner oceans score higher (layers merge faster); spread 24-37%. Every top-three ocean is
 unconfirmed or only likely; among confirmed oceans Ganymede leads, then Titan.

V2: LITERATURE CHECK (2026-10-06, three parallel searches; key DOIs re-fetched from OpenAlex). Changes to the inputs:
  Titan     ocean DISPUTED: Petricca et al. 2025 Nature 10.1038/s41586-025-09818-x, "Titan's strong tidal dissipation
            precludes a subsurface ocean" (170 km shell over slushy high-pressure ice). Kept in the table, labelled.
            If an ocean exists: >180 km thick under <110 km ice (Sohl et al. 2014 JGR Planets 10.1002/2013JE004512).
  Callisto  ocean 30 -> 100 km (20-300) under ~100 km ice: Cascioli et al. 2026 GRL 10.1029/2025GL120812 favour ~300 km
            under 10-120 km ice; Vance et al. 2018 JGR Planets 10.1002/2017JE005341 give 20-139 km. jupiter_moons.py's
            30 km is superseded.
  Ganymede  R 2634 -> 2631 km (Vance 2018). Ocean 100 km under ~140 km ice confirmed (Vance 2018: 75-117 km).
  Triton    ocean 100 -> 150 km: "~150 km [ice] above an ocean of a similar thickness" (Hansen et al. 2021 PSJ,
            quoting Nimmo & Spencer 2015). Ocean likely, not detected; ~10% ammonia in models.
  Pluto     ice 150-200 -> ~90 km (40-180; Nguyen & McGovern 2024 Icarus, EPSL 2024 10.1016/j.epsl.2024.118974).
            Still 51 MPa at the ocean top -> no lid. Ocean > 100 km (Johnson et al. 2016 GRL) confirmed; still disputed.
  Mimas     ocean 30 -> 45 km (35-55): water+ice layer ~70 km (Denton & Rhoden 2025 EPSL 10.1016/j.epsl.2025.119614)
            minus the 20-30 km shell. Age < 25 Myr confirmed.
  Dione     99 +/- 23 km shell, 65 +/- 30 km ocean confirmed (Beuthe et al. 2016); non-unique (Zannoni et al. 2020).
  Uranus    Castillo-Rogez et al. 2023 JGR Planets 10.1029/2022JE007432 Table 3 (ice/ocean km): Ariel 184/30,
            Umbriel 201/20, Titania 222/30, Oberon 209/30; residual oceans < 30 km (Ariel, Umbriel), < 50 km (Titania,
            Oberon), likely hypersaline (~150 g/kg) with 0.3-1 wt% ammonia. Thicker ice -> Ariel 46 MPa, Umbriel 42 MPa
            at the ocean top -> NO lid possible (v1 had 100 km ice ASSUMED, which wrongly allowed one).
            Ariel ocean 50 -> 30 km (5-30); Umbriel 50 -> 20 km (10-30).
  Miranda   ocean >= 100 km under <= 30 km ice, 100-500 Myr ago (Strom et al. 2024 PSJ 10.3847/psj/ad77d7): H range
            30-150 -> 100-150, age 0.1-1 -> 0.1-0.5 Gyr (ASSUMED as the ocean's lifetime). Past ocean only (likely).
  Enceladus, Europa, Rhea (no ocean, Tortora et al. 2016): confirmed; unchanged.
V2 PREDICTIONS (written before the v2 run): Callisto 36% -> ~29% (v1 sensitivity at 100 km); Triton ~27%; Mimas ~5%;
 Ariel and Umbriel lose their lid -> ~36-38% (thin oceans); Miranda ~8%; Ganymede unchanged ~28%. Confirmed oceans
 ranked: Ganymede > Europa > Enceladus (Titan now disputed).
V2 RESULTS (ocean_worlds.log, 2026-10-06): all V2 predictions held.
 final: Umbriel 38.7, Ariel 37.0, Oberon 36.5, Titania 36.4, Pluto 30.2, Callisto 28.5, Ganymede 28.3, Triton 27.7,
 Titan 24.3 (if it has an ocean), Europa 21.7, Dione 21.0, Miranda 8.9 (past ocean), Mimas 4.7, Enceladus 0.0 (%).
 Confirmed oceans: Ganymede 28% > Europa 22% > Enceladus 0%. Without a lid, the 24-39% spread is set by ocean thickness.
"""
import numpy as np

import uncertainty as U

C, DH, DZ, _ = U.WORLDS["Europa"][2:]
UNK = (1, 100)   # salinity prior, g/kg, where unmeasured (ASSUMED)
GYR = (1e8, 4e9)
# name: (R m, H m, H sensitivity, age yr, salinity range or None, S* or None if no lid possible, ocean status)
WORLDS = {
    "Ganymede": (2631e3, 100e3, (50e3, 300e3), GYR, None, None, "confirmed"),
    "Titan": (2575e3, 200e3, (50e3, 400e3), GYR, None, None, "disputed (2025)"),
    "Callisto": (2410e3, 100e3, (20e3, 300e3), GYR, None, None, "likely"),
    "Triton": (1353e3, 150e3, (30e3, 200e3), GYR, None, None, "likely"),
    "Pluto": (1188e3, 100e3, (50e3, 200e3), GYR, None, None, "disputed"),
    "Titania": (789e3, 30e3, (10e3, 50e3), GYR, None, None, "possible"),
    "Oberon": (761e3, 30e3, (10e3, 50e3), GYR, None, None, "possible"),
    "Umbriel": (585e3, 20e3, (10e3, 30e3), GYR, None, None, "possible"),
    "Ariel": (579e3, 30e3, (5e3, 30e3), GYR, None, None, "possible"),
    "Dione": (561e3, 65e3, (30e3, 100e3), GYR, UNK, 5.6, "possible"),
    "Miranda": (236e3, 100e3, (100e3, 150e3), (1e8, 5e8), UNK, 21.9, "past only (100-500 Myr ago)"),
    "Mimas": (198e3, 45e3, (35e3, 55e3), (2e6, 2.5e7), UNK, 22.6, "likely, < 25 Myr old"),
}
STRAT = {"Enceladus": ((5, 20), 20.0), "Europa": ((1, 100), 3.0)}  # strat_sweep.py, unchanged


def args(R, H, t):
    return (4 * np.pi * R ** 2 * H, np.pi * R, C, DH, DZ, t)


def lid(win, s_rng, s_star):
    """Fraction pinned and the win after the lid; no lid possible -> unchanged."""
    if s_star is None:
        return 0.0, win.mean()
    pinned = U.logu(*s_rng) < s_star
    return pinned.mean(), (win & ~pinned).mean()


def main():
    rows = []
    for name, (s_rng, s_star) in STRAT.items():  # first draws in the module RNG, as in strat_sweep.py
        _, win, _ = U.sweep(*U.WORLDS[name])
        assert abs(win.mean() - U.README_WIN_PCT[name]) < U.WIN_TOL, "uncertainty.py headline drifted"
        rows.append((name, win.mean(), *lid(win, s_rng, s_star), "", "confirmed"))
    for name, (R, H, (h_lo, h_hi), t, s_rng, s_star, status) in WORLDS.items():
        _, win, _ = U.sweep(*args(R, H, t))
        p, new = lid(win, s_rng, s_star)
        sens = [lid(U.sweep(*args(R, h, t))[1], s_rng, s_star)[1] for h in (h_lo, h_hi)]
        rows.append((name, win.mean(), p, new, f"{h_lo / 1e3:.0f} km {sens[0]:.0%} / {h_hi / 1e3:.0f} km {sens[1]:.0%}"
                     f" (base {H / 1e3:.0f} km)", status))
    print(f"{'world':10} {'no lid':>7} {'pinned':>7} {'final':>7}   ocean depth sensitivity          ocean")
    for name, old, p, new, sens, status in sorted(rows, key=lambda r: -r[3]):
        print(f"{name:10} {old:7.1%} {p:7.0%} {new:7.1%}   {sens:32} {status}")
    out = {r[0]: r[3] for r in rows}
    assert out["Ganymede"] > max(out["Titan"], out["Europa"]) > out["Enceladus"] == 0, "RESULTS above drifted"
    assert out["Callisto"] < 0.32 and out["Mimas"] < 0.08, "V2 RESULTS drifted"
    return out


if __name__ == "__main__":
    main()
