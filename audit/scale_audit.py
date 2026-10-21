"""Blind reimplementation audit of the scale rule (journals/ija_draft.md Eqs 10-13).

Written from the paper text alone (model statement, Secs 2.1-2.2, Eqs 6-13). No other script
was opened before the predictions below were saved.
Prediction block written and hashed BEFORE any run:
  date:   Mon Oct 19 23:36:30 PDT 2026
  sha256: f635e8d62ef6b22ca0fea4ed4966bb7e9c4555a7e0e63c9d7938cec9f506543a
  (shasum -a 256 of the exact text between the BEGIN and END lines; main() re-hashes and asserts)
-----BEGIN PREDICTIONS-----
Rebuilt from journals/ija_draft.md model statement, Secs 2.1-2.2, Eqs 6-13 only.
RULE: Delta = sqrt(2) pi^(1/4) g (k tau)^(1/4) sqrt(N)  [Eq 10 with gamma=1/(k tau), eps=1/N]
      N(Delta=2) = 2 / (sqrt(pi) g^2 sqrt(k tau))
P1 N(g=1e-17,   k tau=10) = 3.568e33
P2 N(g=4.4e-18, k tau=10) = 1.843e34
P3 N(g=2.64e-16,k tau=10) = 5.120e30
P4 side of cube holding P1's N (N = c V 1000 N_A, N_A=6.02214076e23, f=1): 1 uM 18.10 km; 1 mM 1.810 km; 0.1 M 0.3899 km
P5 pond (1 m3 at 1 M = 6.022e26 molecules) shortfall, k tau=10: g=4.4e-18 3.06e7; g=1e-17 5.93e6; g=2.64e-16 8.5e3
MC of Eq 6, dalpha=(gamma t alpha - alpha^3 + g)dt + sqrt(eps) dW, alpha(t0)=0, t0=-5/sqrt(gamma) to +5/sqrt(gamma), 200000 paths, Euler, 4000 steps;
 P_MC = fraction ending alpha>0. Prediction P_MC = Phi(Delta):
 A: k tau=1e2 N=1e8 g=3.359e-5  Delta=2 -> P=0.9772
 B: k tau=1e4 N=1e8 g=1.062e-5  Delta=2 -> P=0.9772
 C: k tau=1e2 N=1e10 g=1.680e-6 Delta=1 -> P=0.8413
Pass criterion: |P_MC - Phi(Delta)| < 4 standard errors + 0.002 at every point. A miss beyond that means the prefactor or an exponent of Eq 13 is wrong.
Expected disagreement: none beyond finite-step Euler error (< 0.002).
-----END PREDICTIONS-----
Run from the repo dir: /opt/anaconda3/bin/python3 audit/scale_audit.py | tee audit/scale_audit.log
"""
import hashlib, re, sys
from math import pi, sqrt
from pathlib import Path
from statistics import NormalDist
import numpy as np
from scipy.integrate import quad

NA = 6.02214076e23
HERE = Path(__file__).resolve().parent
PHI = NormalDist().cdf


def delta(g, kt, n):          # Eq 11
    return sqrt(2) * pi ** 0.25 * g * kt ** 0.25 * sqrt(n)


def n_needed(g, kt):          # Eq 13
    return 2 / (sqrt(pi) * g * g * sqrt(kt))


def side_km(n, c):            # Eq 12 with f = 1: N = c V 1000 N_A, cube side = V^(1/3)
    return (n / (c * 1000 * NA)) ** (1 / 3) / 1e3


def g_for(d, kt, n):          # invert Eq 11
    return d / (sqrt(2) * pi ** 0.25 * kt ** 0.25 * sqrt(n))


def mc(g, eps, gamma, npaths=200000, nsteps=4000, seed=1):
    """Euler-Maruyama of Eq 6 with lambda = gamma t, alpha(t0)=0; returns fraction ending alpha>0."""
    rng = np.random.default_rng(seed)
    T = 5 / sqrt(gamma)
    dt = 2 * T / nsteps
    a = np.zeros(npaths)
    z = np.empty(npaths)
    sq = sqrt(eps * dt)
    for i in range(nsteps):
        t = -T + i * dt
        rng.standard_normal(out=z)
        a += (gamma * t * a - a ** 3 + g) * dt + sq * z
    return float((a > 0).mean())


def main():
    out = []
    p = lambda s: (out.append(s), print(s))
    # 0. the docstring block is the one hashed before the run
    blk = re.search(r"-----BEGIN PREDICTIONS-----\n(.*?)-----END PREDICTIONS-----", __doc__, re.S).group(1)
    sha = hashlib.sha256(blk.encode()).hexdigest()
    assert sha == "f635e8d62ef6b22ca0fea4ed4966bb7e9c4555a7e0e63c9d7938cec9f506543a", sha
    p(f"prediction block sha256 verified: {sha}")
    p(f"script sha256 at run: {hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}")

    # 1. Eq 10 prefactor from Eq 9 by numerical integration (gamma=0.37, g=1, eps=1)
    gam = 0.37
    mean = quad(lambda s: np.exp(-gam * s * s / 2), -np.inf, np.inf)[0]
    var = quad(lambda s: np.exp(-gam * s * s), -np.inf, np.inf)[0]
    d_num = mean / sqrt(var)
    d_cls = sqrt(2) * pi ** 0.25 * gam ** -0.25
    assert abs(d_num / d_cls - 1) < 1e-9
    p(f"Eq 9 -> Eq 10 prefactor by quadrature: {d_num:.8f} vs closed form {d_cls:.8f}")

    # 2. N at the three g, k tau = 10
    N = {g: n_needed(g, 10) for g in (1e-17, 4.4e-18, 2.64e-16)}
    for g, n in N.items():
        assert abs(delta(g, 10, n) - 2) < 1e-12
        p(f"N(g={g:g}, k tau=10) = {n:.4e}  [Delta check = {delta(g, 10, n):.12f}]")
    assert f"{N[1e-17]:.3e}" == "3.568e+33" and f"{N[4.4e-18]:.3e}" == "1.843e+34" and f"{N[2.64e-16]:.3e}" == "5.120e+30"

    # 3. cube sides for g=1e-17
    sides = {c: side_km(N[1e-17], c) for c in (1e-6, 1e-3, 0.1)}
    for c, s in sides.items():
        p(f"cube side, c={c:g} M, N={N[1e-17]:.4e}: {s:.4f} km")
    assert [f"{s:.4g}" for s in sides.values()] == ["18.1", "1.81", "0.3899"]

    # 4. pond (1 m3 at 1 M) shortfall at k tau = 10
    pond = 1.0 * 1000 * NA
    sf = {g: N[g] / pond for g in (4.4e-18, 1e-17, 2.64e-16)}
    p(f"pond 1 m3 at 1 M = {pond:.4e} molecules; shortfall N/pond: " + ", ".join(f"g={g:g}: {v:.3g}" for g, v in sf.items()))
    assert [f"{v:.3g}" for v in sf.values()] == ["3.06e+07", "5.93e+06", "8.5e+03"]

    # 5. Monte Carlo of Eq 6
    p("MC of Eq 6 (200000 paths, 4000 Euler steps, t from -5/sqrt(gamma) to +5/sqrt(gamma)):")
    ok = True
    for i, (nm, kt, n, d) in enumerate((("A", 1e2, 1e8, 2), ("B", 1e4, 1e8, 2), ("C", 1e2, 1e10, 1))):
        g = g_for(d, kt, n)
        pm = mc(g, 1 / n, 1 / kt, seed=20261019 + i)
        pt = PHI(delta(g, kt, n))
        se = sqrt(pt * (1 - pt) / 200000)
        dm = NormalDist().inv_cdf(pm)
        good = abs(pm - pt) < 4 * se + 0.002
        ok &= good
        p(f"  {nm}: k tau={kt:g} N={n:g} g={g:.4e}  Delta_formula={delta(g, kt, n):.4f}  P_formula={pt:.5f}  "
          f"P_MC={pm:.5f} (SE {se:.5f}, diff {pm - pt:+.5f} = {(pm - pt) / se:+.2f} SE)  Delta_MC=Phi^-1(P_MC)={dm:.4f}  "
          f"{'PASS' if good else 'FAIL'}")
    assert ok, "MC disagrees with Phi(Delta) beyond the pre-registered tolerance"
    p("self-check: all asserts passed")
    (HERE / "scale_audit_result.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
