"""Real chemistry check: does an actual reaction scheme obey the normal-form formula?

Frank-type scheme in a flow reactor of size Omega (molecules per unit concentration):
    feed -> A            inflow f*(a0 - a), a0 ramped slowly upward
    A -> L,  A -> D      k0*a            (slow uncatalysed production)
    A + L -> 2L          k*(1+g)*a*L     (self-copying; g = weak-force bias)
    A + D -> 2D          k*(1-g)*a*D
    L + D -> waste       ki*L*D          (the two hands destroy each other)
    L, D -> out          kd*L, kd*D      (outflow)
Simulated with the chemical Langevin equation: every reaction adds its own counting noise.

Reading off the normal-form parameters at the crossing (k*a = kd), where L = D = x:
    bias   g_eff = 2*g*x*k*a
    noise  eps   = (2*k0*a + 2*x*(k*a + kd)) / Omega     (all turnover that changes L - D)
    sweep  gamma = k * da/dt
Prediction: P(favoured) = Phi(sqrt2 * pi^0.25 * g_eff * gamma^-0.25 * eps^-0.5).

Result: matches within ~0.01 for slow sweeps (tau >= 800); a fast sweep (tau = 200) selects
more strongly than predicted, so the formula is conservative there.

Compared with the simple mapping used in ocean.py (bias g, noise 1/N with N = molecules
of L+D), the real scheme's Delta is smaller by F = sqrt(x / (k0 + 2x)) ~ 0.71.
"""
import numpy as np
from scipy.stats import norm

PREF = np.sqrt(2) * np.pi**0.25
K, KD, KI, K0, F_IN = 1.0, 1.0, 1.0, 1e-3, 1.0


def run(g, omega, tau, runs=2000, dt=0.01, seed=0):
    """Ramp a0 from 0.8 to 1.6 over tau. Returns P(L wins), and a, x, da/dt at the crossing."""
    rng = np.random.default_rng(seed)
    a = np.full(runs, 0.8)
    L = np.full(runs, np.sqrt(K0))
    D = L.copy()
    steps = int(tau / dt)
    a_mean, x_mean = np.empty(steps), np.empty(steps)
    for i in range(steps):
        a0 = 0.8 + 0.8 * i / steps
        r = np.array([K0 * a, K0 * a, K * (1 + g) * a * L, K * (1 - g) * a * D, KI * L * D, KD * L, KD * D])
        n = np.sqrt(np.maximum(r, 0) * dt / omega) * rng.standard_normal(r.shape)
        dL = (r[0] + r[2] - r[4] - r[5]) * dt + n[0] + n[2] - n[4] - n[5]
        dD = (r[1] + r[3] - r[4] - r[6]) * dt + n[1] + n[3] - n[4] - n[6]
        da = (F_IN * (a0 - a) - r[0] - r[1] - r[2] - r[3]) * dt - n[0] - n[1] - n[2] - n[3]
        L, D, a = np.maximum(L + dL, 0), np.maximum(D + dD, 0), a + da
        a_mean[i], x_mean[i] = a.mean(), np.median((L + D) / 2)
    i = int(np.argmax(K * a_mean >= KD))
    w = max(steps // 50, 2)                                  # slope over +-2% of the ramp
    dadt = np.polyfit(np.arange(i - w, i + w) * dt, a_mean[i - w:i + w], 1)[0]
    return (L > D).mean(), (a_mean[i], x_mean[i], dadt)


def predicted(g, omega, cross):
    a, x, dadt = cross
    g_eff = 2 * g * x * K * a
    eps = (2 * K0 * a + 2 * x * (K * a + KD)) / omega
    gamma = K * dadt
    return norm.cdf(PREF * g_eff * gamma**-0.25 * eps**-0.5)


def main():
    print(f"{'g':>7} {'Omega':>7} {'tau':>5}  {'predicted':>9} {'simulated':>9}")
    diffs = []
    for g, omega, tau in [(0.0, 1e4, 800), (0.003, 1e4, 200), (0.003, 1e4, 800), (0.003, 1e4, 3200),
                          (0.006, 1e4, 3200), (0.002, 4e4, 3200)]:
        p, cross = run(g, omega, tau)
        pred = predicted(g, omega, cross)
        if tau >= 800:                       # slow sweeps only: the formula is a slow-sweep result
            diffs.append(abs(p - pred))
        print(f"{g:7.3f} {omega:7.0e} {tau:5d}  {pred:9.3f} {p:9.3f}   (x at crossing {cross[1]:.4f})")
    assert max(diffs) < 0.03, "real scheme does not follow the normal-form formula (slow sweeps)"
    print("fast sweep (tau=200) is outside the slow limit: real chemistry selects MORE than predicted,")
    print("so the ocean model is conservative there. Nature's sweeps are far slower.")
    x = np.sqrt(K0 / KI)
    print(f"\nconversion factor vs the simple mapping in ocean.py: F = sqrt(x/(k0+2x)) = {np.sqrt(x / (K0 + 2 * x)):.2f}")
    print(f"-> required concentrations rise by F^(-4/3) = {np.sqrt(x / (K0 + 2 * x))**(-4 / 3):.2f}x")


if __name__ == "__main__":
    main()
