"""Part 8 sweep with a wobbling environment (wobble.py). Worst-case rule, MEASURED as a lower bound in
wobble.py (every simulated row sat at or above it): Delta *= (1 + r)^(-1/4), r = delta * 2*pi*tau / P,
where delta = relative swing of the control parameter and P = the wobble period."""
import numpy as np
import uncertainty as u
from ocean import YR

orig = u.delta_at

def run(d_rng, p_rng):
    def wobbly(c, g, k2, V, tau, f=1.0):
        rng = np.random.default_rng(1)  # same draws on every call within a sweep
        d = 10 ** rng.uniform(*np.log10(d_rng), len(tau))
        P = 10 ** rng.uniform(*np.log10(p_rng), len(tau)) * YR
        return orig(c, g, k2, V, tau, f) * (1 + d * 2 * np.pi * tau / P) ** -0.25
    u.delta_at = wobbly
    out = {}
    for name, w in u.WORLDS.items():
        u.RNG = np.random.default_rng(0)
        out[name] = u.sweep(*w)[1].mean()
    u.delta_at = orig
    return out

for label, d_rng, p_rng in [("no wobble", (1e-12, 1e-12 * 1.0001), (1, 1.0001)),
                            ("mild: 0.1-1% swing, 1e5-1e8 yr", (1e-3, 1e-2), (1e5, 1e8)),
                            ("plausible: 0.1-30% swing, 1e3-1e8 yr", (1e-3, 0.3), (1e3, 1e8)),
                            ("harsh: 10-30% swing, 1e3-1e5 yr", (0.1, 0.3), (1e3, 1e5))]:
    r = run(d_rng, p_rng)
    print(f"{label:40s} " + "  ".join(f"{k} {v:.0%}" for k, v in r.items()))
