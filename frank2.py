"""Second reaction scheme: Frank + reverse (wasted) reactions. Does the formula still hold, and
how much does wasted turnover lower the conversion factor F?

    A -> L          k0*a        L -> A          k0r*L          (reversible uncatalysed step)
    A + L -> 2L     k(1+g)*a*L  2L -> A + L     kr*L^2         (reversible self-copying)
    L + D -> waste  ki*L*D      L -> out        kd*L           (+ mirror reactions for D)
    feed -> A       f*(a0 - a), a0 ramped slowly

Normal-form parameters read off at the crossing (L = D = x):
    lam   = k*a - k0r - kd - 2*kr*x                (growth rate of L - D)
    g_eff = 2*g*k*a*x
    eps   = (2*k0*a + 2*k0r*x + 2*k*a*x + 2*kr*x**2 + 2*kd*x) / Omega
    gamma = d(lam)/dt
F versus the simple mapping (bias g per self-copying time, noise 1/N, N = 2*x*Omega):
    F = (g_eff / (g*r)) * sqrt(r / (eps*N)),   r = k*a
"""
import numpy as np
from scipy.stats import norm

PREF = np.sqrt(2) * np.pi**0.25
K, KD, KI, K0, K0R, KR, F_IN = 1.0, 1.0, 1.0, 1e-3, 0.2, 0.5, 1.0


def run(g, omega, tau, runs=2000, dt=0.01, seed=0):
    rng = np.random.default_rng(seed)
    a = np.full(runs, 1.0)
    L = np.full(runs, 0.03)
    D = L.copy()
    steps = int(tau / dt)
    lam, a_m, x_m = np.empty(steps), np.empty(steps), np.empty(steps)
    for i in range(steps):
        a0 = 1.0 + 1.2 * i / steps
        r = np.array([K0 * a, K0 * a, K0R * L, K0R * D, K * (1 + g) * a * L, K * (1 - g) * a * D,
                      KR * L * L, KR * D * D, KI * L * D, KD * L, KD * D])
        n = np.sqrt(np.maximum(r, 0) * dt / omega) * rng.standard_normal(r.shape)
        dL = (r[0] - r[2] + r[4] - r[6] - r[8] - r[9]) * dt + n[0] - n[2] + n[4] - n[6] - n[8] - n[9]
        dD = (r[1] - r[3] + r[5] - r[7] - r[8] - r[10]) * dt + n[1] - n[3] + n[5] - n[7] - n[8] - n[10]
        da = (F_IN * (a0 - a) - r[0] - r[1] + r[2] + r[3] - r[4] - r[5] + r[6] + r[7]) * dt \
            - n[0] - n[1] + n[2] + n[3] - n[4] - n[5] + n[6] + n[7]
        L, D, a = np.maximum(L + dL, 0), np.maximum(D + dD, 0), a + da
        a_m[i], x_m[i] = a.mean(), np.median((L + D) / 2)
        lam[i] = K * a_m[i] - K0R - KD - 2 * KR * x_m[i]
    i = int(np.argmax(lam >= 0))
    w = max(steps // 50, 2)
    gamma = np.polyfit(np.arange(i - w, i + w) * dt, lam[i - w:i + w], 1)[0]
    return (L > D).mean(), (a_m[i], x_m[i], gamma)


def mapped(g, omega, cross):
    a, x, gamma = cross
    g_eff = 2 * g * K * a * x
    eps = (2 * K0 * a + 2 * K0R * x + 2 * K * a * x + 2 * KR * x**2 + 2 * KD * x) / omega
    r, N = K * a, 2 * x * omega
    F = (g_eff / (g * r)) * np.sqrt(r / (eps * N)) if g else np.nan
    return norm.cdf(PREF * g_eff * gamma**-0.25 * eps**-0.5), F


def main():
    print(f"{'g':>7} {'Omega':>7} {'tau':>5}  {'predicted':>9} {'simulated':>9}   F")
    diffs, Fs = [], []
    for g, omega, tau in [(0.0, 1e4, 1600), (0.002, 1e4, 1600), (0.004, 1e4, 1600), (0.002, 4e4, 1600)]:
        p, cross = run(g, omega, tau)
        pred, F = mapped(g, omega, cross)
        diffs.append(abs(p - pred))
        if g:
            Fs.append(F)
        print(f"{g:7.3f} {omega:7.0e} {tau:5d}  {pred:9.3f} {p:9.3f}   {F:.2f}")
    assert max(diffs) < 0.03, "second scheme does not follow the formula"
    F = float(np.mean(Fs))
    print(f"\nF = {F:.2f} (scheme 1 in frank.py: 0.70) -> concentrations x{F**(-4/3):.1f}")


if __name__ == "__main__":
    main()
