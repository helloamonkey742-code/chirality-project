"""Earth-only test: did anything unique to early Earth pick L for the whole planet?

Earth's one chiral nudge the icy moons lack is a steady magnetic dipole: magnetite magnetized by it
gives ~60% ee (Ozturk 2023), but with OPPOSITE hands in the two hemispheres (README Part 6 open issue).
Patches heal and the majority wins (Part 3), so the whole-planet push from the field is
    g_mag = f_mag * ee_mag * A * R * s
  f_mag  share of the deciding molecules made on magnetized magnetite   1e-8..1e-2  ASSUMED
  ee_mag excess they carry                                              0.1..0.6    Ozturk 2023
  A      hemisphere imbalance in where they are made                    0.01..0.5   ASSUMED
  R      reversal averaging: polarity blocks of length P_rev inside the decision window W = sqrt(tau/k)
         add like a random walk, R = min(1, 0.53*sqrt(P_rev/W)) (0.53 MEASURED in field_flips.py); P_rev 1e5..1e8 yr (Archean rate unknown)
  s      field polarity during the decision: +1 or -1 at random (nothing ties it to L)
Competing nudges, same units (fractional rate bias, like g = dE/kT):
  weak force  g_w = 4e-19..4e-16, sign set by the run (its sign is unsettled, README Part 5)
  meteorites  g_met = f_met * ee_met, L-sign; f_met 1e-4..1 share of amino acids delivered,
              ee_met 1e-3..0.1 excess passed to proteinogenic amino acids (ASSUMED; isovaline is 18.5%)
                    -- same for every Solar System body, so NOT unique to Earth (README Part 4b)
All other inputs and the one-ocean conditions are uncertainty.py's Early Earth draws, unchanged.
Wall pinning (MEASURED, field_walls.py, 6 runs, each within one grid step of the rule): when the wall is
thin (eps = sqrt(D_h/k)/R < 1) each hemisphere picks its own hand, and the field's local push h = f_mag*ee_mag
holds the wall in place, leaving TWO hands, if h > 0.47*eps*(A + 3*eps). Otherwise the bigger side wins.
Pinned draws count as "no one hand".
"""
import numpy as np
import uncertainty as u
from ocean import YR

orig = u.delta_at


def logu(rng, lo, hi, n):
    return 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)


def earth_dh():
    """Re-draw uncertainty.sweep's D_h for Early Earth (same seed and order: g, k2, c, D_h)."""
    r = np.random.default_rng(0)
    for lo, hi in [(4e-19, 4e-16), (1e-6, 1.0), u.WORLDS["Early Earth ocean"][2]]:
        r.uniform(np.log10(lo), np.log10(hi), u.N)
    return 10 ** r.uniform(*np.log10(u.WORLDS["Early Earth ocean"][3]), u.N)


DH = earth_dh()
PINNED = np.zeros(u.N, bool)


def run(weak_sign, mag, met):
    def biased(c, g, k2, V, tau, f=1.0):
        n = len(tau)
        rng = np.random.default_rng(1)  # identical draws on every call within a sweep
        W = np.sqrt(tau / (k2 * c))
        R = np.minimum(1, 0.53 * np.sqrt(logu(rng, 1e5, 1e8, n) * YR / W))  # 0.53 MEASURED, field_flips.py
        g_mag = logu(rng, 1e-8, 1e-2, n) * rng.uniform(0.1, 0.6, n) * logu(rng, 0.01, 0.5, n) * R \
            * rng.choice([-1, 1], n)
        g_met = logu(rng, 1e-4, 1, n) * logu(rng, 1e-3, 0.1, n)
        total = weak_sign * g + mag * g_mag + met * g_met
        if mag:
            r2 = np.random.default_rng(1)  # same order as above, to recover f_mag*ee_mag and A
            r2.uniform(np.log10(1e5), np.log10(1e8), n)
            f_ee = logu(r2, 1e-8, 1e-2, n) * r2.uniform(0.1, 0.6, n)
            A = logu(r2, 0.01, 0.5, n)
            eps = np.sqrt(DH / (k2 * c)) / (u.WORLDS["Early Earth ocean"][1] / np.pi)
            PINNED[:] = (eps < 1) & (f_ee > 0.47 * eps * (A + 3 * eps))
            total = np.where(PINNED, 0.0, total)
        return orig(c, 1.0, k2, V, tau, f) * total * u_sign  # Delta toward the hand being counted
    u.delta_at = biased
    out = []
    global u_sign
    for u_sign in (1, -1):  # count decisive L, then decisive D
        u.RNG = np.random.default_rng(0)
        out.append(u.sweep(*u.WORLDS["Early Earth ocean"])[1].mean())
    u.delta_at = orig
    return out


def main():
    print(f"{'Early Earth ocean':52s} {'L wins':>7} {'D wins':>7} {'no one hand':>11}")
    rows = {}
    for label, ws, mag, met in [
        ("weak force only (favours L)", 1, 0, 0),
        ("weak force only (favours D)", -1, 0, 0),
        ("weak (L) + Earth's magnetic field", 1, 1, 0),
        ("weak (D) + Earth's magnetic field", -1, 1, 0),
        ("weak (D) + magnetic field + meteorites", -1, 1, 1),
    ]:
        L, D = run(ws, mag, met)
        rows[label] = (L, D)
        print(f"{label:52s} {L:7.0%} {D:7.0%} {1 - L - D:11.0%}" + (f"   (field pins two hands in {PINNED.mean():.0%})" if mag else ""))
    assert abs(rows["weak force only (favours L)"][0] - 0.75) < 0.03  # README Part 8
    L, D = rows["weak (L) + Earth's magnetic field"]
    assert abs(L - D) < 0.05, "field should be a coin flip"


if __name__ == "__main__":
    main()
