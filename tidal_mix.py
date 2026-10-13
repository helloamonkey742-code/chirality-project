"""tidal_mix.py: vertical-mixing time and K_z of each ocean from an energy budget, and what that does to the shares.

Why (2026-10-12): Prof. Francis Nimmo (UCSC) suggested bounding the vertical mixing of each ocean by energy, as Munk &
Wunsch (1998) did for Earth: the potential energy stored in the stratification, divided by the tidal power that
actually reaches the ocean, is a mixing time. His caveat: Earth's stratification is kept up by continuous production of
cold deep water, which is unlikely in the outer solar system. The model's vertical D_z prior for the icy moons is
log-uniform 1e-10..1e-3 m^2/s (uncertainty.py:40, from Zeng & Jansen 2021's 3e-10..3e-3), floored at 1e-9 in sweep2
step c (DZ_FLOOR, the molecular tracer value). This script asks where an energy-budget K_z falls in that prior, and
reruns the per-world shares with D_z restricted to it. No project script, log or journal is edited.

DERIVATION.
  Potential energy. Take a layer of thickness H, z up from its floor, with a linear density profile of buoyancy
  frequency N: rho(z) = rho0 [1 - (N^2/g)(z - H/2)]. Mixing it to uniform rho0 keeps the mass and raises the potential
  energy per unit area by
      PE = int_0^H (rho0 - rho(z)) g z dz = rho0 N^2 int_0^H (z - H/2) z dz = rho0 N^2 (H^3/3 - H^3/4) = rho0 N^2 H^3 / 12.
  A layer with a buoyancy jump db = g drho/rho across it has N^2 = db / H, so PE = rho0 db H^2 / 12.
  Mixing time. Only a fraction Gamma of the turbulent dissipation P (W) goes into potential energy (the mixing
  efficiency; Gamma = 0.2, range 0.1-0.3; Osborn 1980; Munk & Wunsch 1998). So
      tau_mix = PE A / (Gamma P),  A = 4 pi R^2 (the model's ocean volume is 4 pi R^2 H, ocean_worlds.args).
  Diffusivity. Osborn: K_z = Gamma eps / N^2 with eps = P / M, M = rho0 A H (dissipation spread over the ocean). Then
      K_z = Gamma P / (rho0 A db)   (independent of H),   tau_mix = rho0 db H^2 A / (12 Gamma P) = H^2 / (12 K_z).
  So the requested cross-check K_z ~ H^2 / tau holds exactly up to the factor 12 (printed as a check). The model's
  own wall-merge time H R / (2 D_z) (uncertainty.py:64) is printed beside tau_mix.
  Zeng & Jansen 2021 Sec. 2.2 use the same relation for the Enceladus lid: kappa_z = (Gamma E / A) / (rho g |alpha| dT).

WHY THE LID AND THE INTERIOR GET SEPARATE ESTIMATES.
  They are held up by different physics and differ in db by a factor 40.
  Lid (worlds the model gives an S* lid: Enceladus, Europa, Dione, Miranda, Mimas). At low salinity and pressure the
  thermal expansivity alpha < 0 just above freezing, so the water between T_f and T_c (~ +1.4 C) under the ice is
  stably stratified (Zeng & Jansen 2021 Sec. 2.1). Its db = g |alpha| dT with alpha = -4e-5 /K, dT = 2.0 K (their
  8.5 g/kg, 20 bar case). The heat leaking up through it maintains it, so it does not need deep-water formation and
  Nimmo's caveat does not remove it. Its thickness is set by conduction, H_lid = rho c_p kappa dT / Q (their Eq. 1);
  this script scales their own numbers, H_lid = 30 m x (kappa / 1e-7) x (0.03 W/m^2 / Q) x (dT / 2 K), with
  kappa = max(K_lid, 1e-7 molecular), capped at the ocean depth. All of P is assumed to work on the lid (their
  convention, an upper bound). The model's lid is binary (S < S* pins every wall, sweep2.outcome) and has no N^2, so
  the lid K is reported, not fed into the shares.
  Interior (every world). The only published interior density contrast found is Ashkenazy & Tziperman 2021 (Europa,
  50 g/kg, convecting ocean): alpha dT = drho/rho = 2e-6 top to bottom. db_int = g x 2e-6 is used for every world
  (flagged: one world's number). It is 40x weaker than the lid. This is the K that sets the model's D_z (the
  whole-column merge time), so the restricted-prior share run uses the interior band. Nimmo's caveat: with no
  deep-water source the interior relaxes toward neutral (db -> 0), so this K_int is a LOWER bound on interior
  mixing, and tau_int an upper bound.

INPUTS (every DOI checked on Crossref 2026-10-12, curl -s https://api.crossref.org/works/<doi>; title matched).
  [C14]  Chen, Nimmo & Glatzmaier 2014, Icarus 229:11-30, doi:10.1016/j.icarus.2013.10.024. Table 3 (p. 21): ocean
         eccentricity + obliquity tide dissipation (free-surface ocean, no shell, c_D = 0.002), E_rad, solid-body
         tides (k2 = 3/2 column used for the total). Table 2 (p. 15): surface g. Sec. 5.1 (p. 20): a decoupled shell
         (c0 = 0.67) raises obliquity heating "by a factor of roughly 4" -> high = ecc + 4 x obl. Mimas obliquity
         2.5e6 W is bracketed in the table (outside the model's validity), used and flagged.
  [HM19] Hay & Matsuyama 2019, Icarus 319:68-85, doi:10.1016/j.icarus.2018.09.019. Table 3 (accepted manuscript
         p. 20, UA repository hdl 10150/632307): ocean ecc + obl dissipation with an ice shell, c_D = 0.004.
         Triton: thick shell 5.71e10 W, thin shell 2.89e11 W.
  [J23]  Jansen, Kang, Kite, Zeng, Bire, Ferrari 2023, PSJ 4:117, doi:10.3847/psj/acda95. Sec. 4 (tides): Europa ocean
         dissipation 1e9-1e11 W; Enceladus 10-1e4 W, and (citing Rovira-Navarro et al. 2023, doi:10.3847/psj/acae9a)
         "most likely ~10-1e6 W"; Earth deep ocean ~2 TW, eps ~ 1.5e-9 m^2/s^3 (citing Wunsch & Ferrari 2004,
         doi:10.1146/annurev.fluid.36.050802.122121). Sec. 2: Q at the ocean floor 0.02-0.1 (Europa).
  [ZJ21] Zeng & Jansen 2021, PSJ 2:151, doi:10.3847/PSJ/ac1114. Sec. 2.1: lid alpha, dT, H_lid = 30 m (kappa 1e-7) /
         3 km (1e-5) at Q = 0.03 W/m^2 (from ~20 GW, Choblet et al. 2017, doi:10.1038/s41550-017-0289-8). Sec. 2.2:
         Enceladus E ~ 10-1e8 W, upper end = libration O(0.1 GW) (Wilson & Kerswell 2018, doi:10.1016/j.epsl.2018.08.012,
         second-hand), treated by them as an upper bound; kappa_z 3e-10..3e-3.
  [AT21] Ashkenazy & Tziperman 2021, Nat. Commun. 12:6376, doi:10.1038/s41467-021-26710-0. Rotation-scaling paragraph:
         drho/rho = 2e-6, g = 1.314, q = 0.05 W/m^2, rho = 1e3, c_p = 4e3, D = 1e5 m (Europa).
  Gamma: Osborn 1980, doi:10.1175/1520-0485(1980)010<0083:EOTLRO>2.0.CO;2; Munk & Wunsch 1998, doi:10.1016/
         S0967-0637(98)00070-3 (both titles checked; full texts not read; ZJ21 and J23 use Gamma ~ 0.2).
  Read but no usable number: Matsuyama et al. 2018 (doi:10.1016/j.icarus.2018.04.013, values only in figures),
  Beuthe 2016 (doi:10.1016/j.icarus.2016.08.009), Rovira-Navarro et al. 2019 (doi:10.1016/j.icarus.2018.11.010, "several
  orders of magnitude smaller"), Rovira-Navarro et al. 2023 (per-moon dissipation only in figures).
  P_ocean rule: low = smallest published value or range end; mid = HM19 (with a shell) where it exists, else C14
  (Triton: geometric mean of HM19's two shells); high = largest of published values, range ends and C14 ecc + 4 x obl.
  Libration is counted only for Enceladus (ZJ21 gives it). Total heat (C14 rad + solid ecc + solid obl (k2 = 3/2) +
  ocean) is printed as a ceiling, not used in the band. NOT SOURCED: Pluto (in neither C14 nor HM19): full prior kept,
  flagged. Early Earth ocean: J23 gives Earth's eps but no N was read, and its prior (1e-6..1e-4) is not an icy-moon
  prior: full prior kept, flagged.

SHARE RUN. multiwall.multi (sweep3.outcome3 with the measured wall rule), imported unedited, run per world exactly as
multiwall.run4 with seedfix's world-0 stream: extras from default_rng([s, i]) for i >= 1 and default_rng([s, 0, 1]) for
Enceladus (i = 0). Seeds 0-3. Restricted prior: each draw's D_z keeps its quantile and is moved log-uniformly onto the
world's interior band [K_low, K_high] (sweep2.remap from the drawn range a[4]); band ends below DZ_FLOOR are raised to
it (step c); dz_rng = (DZ_FLOOR, 1e-3) is passed so nothing remaps it again. Band ends above 1e-3 are NOT clipped: the
band is the physical estimate; the fraction above the prior is printed. Regression: the full prior must reproduce
multiwall.log's Section 4 three-way main (15 worlds, no lid | final, printed 0.1) with the Enceladus row from
seedfix.log Section D ("new:"), and seedfix.log's corrected seed 1-3 P_pved lines, digit for digit.

PREDICTIONS (written 2026-10-12 from hand estimates before any result-producing run; never edited).
  P1  Regression: the full-prior run reproduces seedfix-corrected multiwall main exactly for all 15 worlds (seed 0,
      three-way, no lid and final) and the seedfix seed 1-3 P_pved lines (no lid / final) exactly.
  P2  Interior K_mid (Gamma = 0.2, P_mid) is above the prior top 1e-3 for Europa, Ganymede, Titan, Callisto, Triton
      and Mimas, and inside [1e-9, 1e-3] for Enceladus, Dione, Miranda, Ariel, Umbriel, Titania and Oberon.
  P3  Europa interior K_low lies between 5e-4 and 5e-3 m^2/s (hand estimate 1.4e-3).
  P4  Enceladus lid at Gamma = 0.2 over E = 10..1e8 W reproduces ZJ21's 3e-10..3e-3 within a factor 1.5 at each end.
  P5  Enceladus lid tau_mix (low..high over P and Gamma) overlaps ZJ21's "hundreds to hundreds of thousands of years".
  P6  Interior tau_mix at mid is shorter than the decision time 1/k_rac (3.5 Myr) for Europa, Ganymede, Titan,
      Callisto, Triton, Mimas, Miranda, Oberon and Titania, and longer for Enceladus, Dione and Ariel.
  P7  Interior tau_mix at mid is shorter than the youngest model ocean age for every sourced world except Enceladus.
  P8  The fraction of the effective prior [1e-9, 1e-3] lying inside the interior band is < 20 % for every sourced
      world except Enceladus, and > 60 % for Enceladus (its band spans 10 W to 1e8 W).
  P9  For the worlds of P2 whose K_mid is above 1e-3 (Europa, Ganymede, Titan, Callisto, Triton, Mimas), the
      restricted-prior P_pved (no lid, seed 0) is higher than the full-prior multiwall P_pved.
  P10 For Ariel, Umbriel and Dione (band predicted ~1e-8..1e-5), the restricted-prior P_pved (no lid, seed 0) is lower
      than the full-prior value.
  P11 Restricted-prior P_pved at seeds 1-3 is within 0.5 points of seed 0 for every world (no lid and final).

Run: OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python tidal_mix.py > tidal_mix.log
     (TIDAL_N=2000 for a crash-only smoke run; the regression is checked only at N = 200000.)
"""
import os
import re

import numpy as np

import multiwall as MW
import sweep2 as S
import sweep3 as S3
import uncertainty as U

HERE = os.path.dirname(os.path.abspath(__file__))
U.N = int(os.environ.get("TIDAL_N", U.N))
RHO, GAM = 1000.0, (0.1, 0.2, 0.3)          # AT21 rho; Gamma low/mid/high
DRHO_INT = 2e-6                              # AT21 drho/rho, Europa interior
ALPHA, DT_LID = 4e-5, 2.0                    # ZJ21 |alpha| (1/K), T_c - T_f (K)
K_MOL = 1e-7                                 # ZJ21 molecular thermal diffusivity
T_DEC = 1 / S.K_RAC["Ala"]                   # decision time, s
YR = 3.156e7
LID_WORLDS = ("Enceladus", "Europa", "Dione", "Miranda", "Mimas")

# C14 Table 3 (W): rad, solid ecc (k2=3/2), solid obl (k2=3/2), ocean ecc, ocean obl; Table 2 g (m/s^2)
C14 = {"Europa": (2.02e11, 7.05e12, 9.70e9, 1.47e7, 3.1e9, 1.31),
       "Ganymede": (4.67e11, 5.46e10, 1.55e9, 1.66e4, 1.3e9, 1.43),
       "Callisto": (3.20e11, 1.65e10, 7.40e8, 8.80e2, 3.8e9, 1.24),
       "Mimas": (4.02e7, 7.59e11, 1.40e8, 2.38e7, 2.5e6, 0.064),
       "Enceladus": (2.73e8, 1.96e10, 8.64e2, 2.32e4, 1.2e1, 0.11),
       "Dione": (2.42e9, 8.44e9, 3.02e5, 4.46e3, 2.2e3, 0.23),
       "Titan": (4.11e11, 4.37e11, 2.35e9, 1.23e5, 1.1e10, 1.35),
       "Miranda": (8.50e7, 1.05e9, 1.21e7, 6.37e2, 2.6e5, 0.079),
       "Ariel": (3.59e9, 4.42e9, 3.40e4, 1.36e3, 5.1e1, 0.27),
       "Umbriel": (2.32e9, 4.08e9, 8.29e4, 9.69e2, 1.8e3, 0.23),
       "Titania": (9.71e9, 3.55e7, 2.46e5, 1.63e1, 1.3e5, 0.38),
       "Oberon": (7.76e9, 5.44e6, 6.98e5, 5.82e3, 3.3e6, 0.35),
       "Triton": (7.12e10, 0.0, 1.60e10, 0.0, 6.6e10, 0.78)}
HM19 = {"Enceladus": 0.959 + 12.0, "Dione": 0.176 + 1.18e4, "Titan": 962 + 1.19e10, "Europa": 1.32e6 + 6.23e9,
        "Ganymede": 1.58e3 + 2.79e9, "Triton": (5.71e10, 2.89e11)}
RANGES = {"Europa": [(1e9, 1e11, "J23")],
          "Enceladus": [(10, 1e4, "J23 tidal"), (10, 1e6, "J23 citing RN23"), (10, 1e8, "ZJ21 incl. libration")]}
Q_LID = {"Enceladus": (0.03, "ZJ21 Sec. 2.1"), "Europa": (0.05, "AT21 q")}   # others: C14 total / A (derived)
FLAG = {"Pluto": "no published ocean tidal power: full prior kept",
        "Early Earth ocean": "Earth calibration only (J23 eps 1.5e-9 m^2/s^3, 2 TW); own prior kept",
        "Mimas": "C14 obliquity value bracketed (outside validity); free-surface ocean, no shell",
        "Callisto": "single source (C14): low = mid", "Titania": "single source (C14): low = mid",
        "Oberon": "single source (C14): low = mid", "Umbriel": "single source (C14): low = mid",
        "Ariel": "single source (C14): low = mid", "Miranda": "single source (C14): low = mid",
        "Enceladus": "high end is libration (ZJ21 upper bound)"}


def power(name):
    """(low, mid, high, total ceiling) ocean power in W, or None where unsourced."""
    if name not in C14:
        return None
    rad, se, so, oe, oo, _ = C14[name]
    vals = [oe + oo, oe + 4 * oo]
    hm = HM19.get(name)
    if isinstance(hm, tuple):
        vals += list(hm)
        mid = float(np.sqrt(hm[0] * hm[1]))
    elif hm is not None:
        vals.append(hm)
        mid = hm
    else:
        mid = oe + oo
    for lo, hi, _ in RANGES.get(name, []):
        vals += [lo, hi]
    return min(vals), mid, max(vals), rad + se + so + oe + oo


def geom(a):
    """(R, H, A) of the model ocean: R = L / pi, H = depth at the nominal volume."""
    R = a[1] / np.pi
    return R, U.depth(a[0], a[1]), 4 * np.pi * R ** 2


def kz(P, A, db, gam):
    return gam * P / (RHO * A * db)


def pct(K, rng):
    lo, hi = np.log10(rng)
    return 100 * np.clip((np.log10(K) - lo) / (hi - lo), 0, 1)


def inside(band, rng):
    """% of the log-uniform prior rng that lies inside band."""
    lo, hi = np.log10(rng)
    a, b = np.log10(band)
    return 100 * max(0.0, min(b, hi) - max(a, lo)) / (hi - lo)


def table(W):
    print("=== 1. Inputs: ocean power P (W), g (m/s^2), model geometry ===")
    print(f"{'world':18} {'P_low':>9} {'P_mid':>9} {'P_high':>9} {'total':>9} {'g':>6} {'R km':>7} {'H km':>6}  flag")
    bands = {}
    for name, (x, win, a, salt) in W.items():
        R, H, A = geom(a)
        p = power(name)
        g = C14[name][5] if name in C14 else np.nan
        if p is None:
            print(f"{name:18} {'-':>9} {'-':>9} {'-':>9} {'-':>9} {'-':>6} {R / 1e3:7.0f} {H / 1e3:6.1f}  {FLAG[name]}")
            continue
        print(f"{name:18} {p[0]:9.2e} {p[1]:9.2e} {p[2]:9.2e} {p[3]:9.2e} {g:6.3f} {R / 1e3:7.0f} {H / 1e3:6.1f}  "
              + FLAG.get(name, ""))
        bands[name] = (p, g, R, H, A)

    print("\n=== 2. Interior: db = g x 2e-6 (AT21). K_z = Gamma P / (rho A db); tau = rho db H^2 A / (12 Gamma P) ===")
    print("K low = Gamma 0.1 x P_low, mid = 0.2 x P_mid, high = 0.3 x P_high. Prior: log-uniform 1e-10..1e-3, "
          "effective 1e-9..1e-3 (step-c floor)")
    print(f"{'world':18} {'K_low':>9} {'K_mid':>9} {'K_high':>9} | {'pct lo/mid/hi (eff. prior)':>27} | "
          f"{'% prior in band':>15} {'% band >1e-3':>12} | {'tau_mid yr':>10} {'tau lo..hi yr':>19} | "
          f"{'H^2/(12 K tau)':>14} | {'merge HR/2K yr':>14} | {'tau/T_dec':>9} {'tau/age_min':>11}")
    out = {}
    for name, (p, g, R, H, A) in bands.items():
        a = W[name][2]
        db = g * DRHO_INT
        K = [kz(P, A, db, G) for P, G in zip(p[:3], GAM)]
        tau = [RHO * db * H ** 2 * A / (12 * G * P) for P, G in zip(p[:3], GAM)]
        chk = H ** 2 / (12 * K[1] * tau[1])
        prior = (S.DZ_FLOOR, a[4][1])
        above = 100 * max(0.0, np.log10(K[2]) - max(np.log10(K[0]), np.log10(a[4][1]))) / np.log10(K[2] / K[0])
        age = a[5][0] * YR
        print(f"{name:18} {K[0]:9.2e} {K[1]:9.2e} {K[2]:9.2e} | "
              + " / ".join(f"{pct(k, prior):6.1f}" for k in K)
              + f" | {inside((K[0], K[2]), prior):15.1f} {above:12.1f} | {tau[1] / YR:10.3g} "
              f"{tau[2] / YR:9.3g}..{tau[0] / YR:9.3g} | {chk:14.6f} | {H * R / (2 * K[1]) / YR:14.3g} | "
              f"{tau[1] / T_DEC:9.3g} {tau[1] / age:11.3g}")
        out[name] = dict(K=K, tau=tau, band=(K[0], K[2]), H=H, age=age, tot=kz(p[3], A, db, GAM[2]))
    print("Ceiling (Gamma 0.3 x C14 total heat, all of it in the ocean): " + ", ".join(
        f"{n} {o['tot']:.1e}" for n, o in out.items()))
    print(f"Decision time 1/k_rac = {T_DEC:.3g} s = {T_DEC / YR / 1e6:.2f} Myr. Nimmo's caveat: no deep-water source "
          "-> interior db smaller -> K larger, tau shorter; these K_int are lower bounds.")

    print("\n=== 3. Lid (worlds with an S* lid): db = g |alpha| dT = g x 4e-5 x 2 K (ZJ21). All of P works on the lid "
          "(ZJ21 convention) ===")
    print(f"{'world':18} {'Q W/m2':>7} {'K_low':>9} {'K_mid':>9} {'K_high':>9} | {'H_lid lo/mid/hi m':>26} | "
          f"{'tau lo/mid/hi yr':>30}  Q source")
    for name in LID_WORLDS:
        p, g, R, H, A = bands[name]
        q, qs = Q_LID.get(name, (p[3] / A, "C14 total / A (derived)"))
        db = g * ALPHA * DT_LID
        K = [kz(P, A, db, G) for P, G in zip(p[:3], GAM)]
        Hl = [min(H, 30.0 * (max(k, K_MOL) / K_MOL) * (0.03 / q) * (DT_LID / 2.0)) for k in K]
        tau = [h ** 2 / (12 * k) for h, k in zip(Hl, K)]
        print(f"{name:18} {q:7.3f} {K[0]:9.2e} {K[1]:9.2e} {K[2]:9.2e} | "
              + " / ".join(f"{h:8.3g}" for h in Hl) + " | " + " / ".join(f"{t / YR:9.3g}" for t in tau) + f"  {qs}")
    p, g, R, H, A = bands["Enceladus"]
    db = g * ALPHA * DT_LID
    print(f"ZJ21 check, Enceladus lid at Gamma 0.2, E = 10 and 1e8 W: {kz(10, A, db, 0.2):.2e} .. "
          f"{kz(1e8, A, db, 0.2):.2e} m^2/s (ZJ21: 3e-10 .. 3e-3)")
    return out


# ---------------------------------------------------------------- shares
def run(seed, W, bands):
    """{name: {'full': [no lid, final], 'band': [...]}} three-way splits, multiwall main rule."""
    out = {}
    for i, (name, (x, win, a, salt)) in enumerate(W.items()):
        rng = np.random.default_rng([seed, 0, 1] if i == 0 else [seed, i])   # seedfix world-0 stream
        ex = {"eta": 10 ** rng.uniform(-1, 1, U.N), "r": 10 ** rng.uniform(-2, 0, U.N), "rho": 10 ** rng.uniform(0, 2, U.N)}
        none = np.zeros(U.N, bool)
        o = {"full": [S3.split(*MW.multi(x, a[1], a[4], ex, l)[:2]) for l in (none, salt)]}
        if name in bands:
            lo, hi = (max(k, S.DZ_FLOOR) for k in bands[name])
            xb = dict(x)
            xb["vertical D_z"] = S.remap(x["vertical D_z"], a[4], (lo, hi))
            o["band"] = [S3.split(*MW.multi(xb, a[1], (S.DZ_FLOOR, 1e-3), ex, l)[:2]) for l in (none, salt)]
        else:
            o["band"] = o["full"]
        out[name] = o
    return out


def targets():
    """Seed-0 three-way main per world and seed 1-3 P_pved strings, multiwall.log with seedfix.log's Enceladus."""
    mw = open(os.path.join(HERE, "multiwall.log")).read().splitlines()
    head = [j for j, l in enumerate(mw) if l.startswith("Three-way split, %: P_pved / P_chance / P_none, no lid | final")]
    assert len(head) == 1, head
    num = r"([\d.]+) /\s+([\d.]+) /\s+([\d.]+)"
    pat = re.compile(rf"^([A-Za-z][A-Za-z ]*?)\s+[\d.].*\|\|\s+{num} \|\s+{num}\s*$")
    t0 = {}
    for l in mw[head[0] + 2: head[0] + 17]:
        m = pat.match(l)
        t0[m.group(1).strip()] = m.groups()[1:]
    sf = open(os.path.join(HERE, "seedfix.log")).read().split("=== D.")[1].split("=== E.")[0].splitlines()
    enc = [l for l in sf if l.strip().startswith("new: Enceladus") and "||" in l]
    assert len(enc) == 1, enc
    t0["Enceladus"] = pat.match(enc[0].split("new: ", 1)[1]).groups()[1:]
    ts = {}
    for l in sf:
        m = re.match(r"\s*new:\s+seed (\d): max diff [\d.]+\s+(.*)$", l)
        if m:
            ts[int(m.group(1))] = {n.strip(): v for n, v in re.findall(r"([A-Za-z][A-Za-z ]*?) ([\d.]+/[\d.]+)", m.group(2))}
    assert len(t0) == 15 and sorted(ts) == [1, 2, 3] and all(len(v) == 15 for v in ts.values()), (t0, ts)
    return t0, ts


def shares(W, bands):
    print(f"\n=== 4. Shares, multiwall main rule (multiwall.multi, imported unedited), N = {U.N}, seeds 0-3, "
          "seedfix world-0 stream ===")
    res = {s: run(s, W if s == 0 else S.chains(s)[1], bands) for s in range(4)}
    r0 = res[0]
    if U.N == 200_000:
        t0, ts = targets()
        bad = [(n, t0[n], r0[n]["full"]) for n in r0 if tuple(f"{v:.1f}" for s in r0[n]["full"] for v in s) != tuple(t0[n])]
        bad += [(s, n, ts[s][n]) for s in (1, 2, 3) for n in res[s]
                if f"{res[s][n]['full'][0][0]:.1f}/{res[s][n]['full'][1][0]:.1f}" != ts[s][n]]
        print("regression, full prior vs multiwall.log Section 4 main with seedfix.log's Enceladus (seed 0, 15 worlds, "
              f"three-way, no lid and final) and seedfix.log seeds 1-3 P_pved: {'EXACT' if not bad else 'FAILED ' + str(bad)}")
    else:
        print(f"regression not checked (N = {U.N}, smoke)")
    print("\nThree-way split %: P_pved / P_chance / P_none, no lid | final. 'band' = D_z log-uniform over the interior "
          "K band (section 2)")
    print(f"{'world':18} {'full prior':>36} || {'interior band':>36} || d P_pved no lid / final")
    order = sorted(r0, key=lambda n: -r0[n]["full"][0][0])
    for n in order:
        f, b = r0[n]["full"], r0[n]["band"]
        tag = "" if n in bands else "   [full prior kept: " + FLAG[n] + "]"
        print(f"{n:18} " + " | ".join(f"{p:5.1f} / {c:4.1f} / {z:5.1f}" for p, c, z in f) + "  ||  "
              + " | ".join(f"{p:5.1f} / {c:4.1f} / {z:5.1f}" for p, c, z in b)
              + f"  || {b[0][0] - f[0][0]:+6.1f} / {b[1][0] - f[1][0]:+6.1f}{tag}")
    print("\nSeed stability, interior band P_pved no lid / final, seeds 1-3, and max |seed k - seed 0|:")
    for n in order:
        d = max(abs(res[s][n]["band"][l][0] - r0[n]["band"][l][0]) for s in (1, 2, 3) for l in (0, 1))
        print(f"{n:18} " + "  ".join(f"{res[s][n]['band'][0][0]:5.1f}/{res[s][n]['band'][1][0]:5.1f}" for s in (1, 2, 3))
              + f"   max diff {d:.2f}")
    return res


def main():
    _, W = S.chains(0)
    print(f"tidal_mix.py, N = {U.N}. Energy-budget vertical mixing (Nimmo's suggestion; Munk & Wunsch 1998)\n")
    out = table(W)
    shares(W, {n: o["band"] for n, o in out.items()})


if __name__ == "__main__":
    main()
