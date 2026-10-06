"""Simulation check for earth.py's reversal factor R = min(1, sqrt(P_rev / W)).

The field flips polarity every ~P_rev, so its push on the whole planet changes sign in blocks. Only the
push during the decision window W (= 1/sqrt(gamma) in sim.py's units, sqrt(tau/k) in real time) counts.
Normal form as in sim.py, with a bias that flips sign in blocks of length P (random phase):
    d(alpha) = (gamma*t*alpha - alpha^3 + G*s(t)) dt + sqrt(eps) dW
Measured: the effective bias G_eff = Phi^-1(P(alpha ends with the sign of the window-weighted bias))
/ (Delta per unit bias), against the predicted typical size G*R.
"""
import numpy as np
from scipy.stats import norm
from sim import delta

EPS, GAMMA, RUNS = 1e-4, 0.01, 4000
W = 1 / np.sqrt(GAMMA)
D1 = delta(1.0, EPS, GAMMA)          # Delta per unit (steady) bias
G = 3.0 / D1                          # steady bias -> Delta = 3


def simulate(P, seed=0):
    rng = np.random.default_rng(seed)
    dt, t0 = 0.01, 6 * W
    ts = np.arange(-t0, t0, dt)
    phase = rng.uniform(0, P, RUNS)
    nblocks = int(2 * t0 / P) + 2
    signs = rng.choice([-1.0, 1.0], (RUNS, nblocks))
    a = np.zeros(RUNS)
    weighted = np.zeros(RUNS)
    for t in ts:
        s = signs[np.arange(RUNS), ((t + t0 + phase) // P).astype(int)]
        weighted += s * np.exp(-GAMMA * t**2 / 2) * dt
        a += (GAMMA * t * a - a**3 + G * s) * dt + np.sqrt(EPS * dt) * rng.standard_normal(RUNS)
    b = weighted / np.sqrt(2 * np.pi / GAMMA)  # window-weighted bias, in units of G (1 = steady)
    return a, b


def main():
    print(f"{'P/W':>6} {'rms weighted bias':>18} {'R=min(1,sqrt(P/W))':>19} {'follows its sign':>16} {'predicted':>9}")
    for P_over_W in [0.03, 0.1, 0.3, 1, 3, 10]:
        a, b = simulate(P_over_W * W)
        rms = np.sqrt((b**2).mean())
        follows = (np.sign(a) == np.sign(b)).mean()
        pred = norm.cdf(3 * np.abs(b)).mean()  # Phi(Delta * |b|) averaged over the realized biases
        print(f"{P_over_W:6.2f} {rms:18.3f} {min(1, np.sqrt(P_over_W)):19.3f} {follows:16.3f} {pred:9.3f}", flush=True)


if __name__ == "__main__":
    main()
