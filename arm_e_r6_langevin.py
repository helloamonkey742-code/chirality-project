"""Arm E R6, P11: chemical-Langevin (CLE) simulation of the full 5-species model against the linear-noise formulas used for e_min.

What is tested. The e_min of arm_e_r6.py / arm_e_r6_analytic.py rests on three linear-noise statements: the stationary readout sd (Lyapunov),
P(correct sign) = Phi(Delta) with Delta = mean/sd, and, for the unstable mode, Delta_b with accumulated noise. Here the 14 channels are
integrated directly (Euler-Maruyama, M trajectories, N = Omega*C0 = 3e6 molecules of feed species) and compared. Channels, rates and the
linear-noise values come from arm_e_r6_analytic.py (props, NU, racemic, jac, analyse, lyap); the sampler itself shares no code with it.

FIXED BEFORE THE RUN (operational reading of P11 in arm_e_r6_predictions.txt, whose sha is asserted):
  Points   V1/V2a: pond record a|feed_a|E of arm_e_r6_results.json (monomer readout, engaged subset): f = 1.778e-6 /s, kh = 3.162e-4,
           R = 1, K = 1e-6, F = 1, s = 6, C0 = 1e-2 M. V2b: pond record b|feed_strict|U: F = 553.03, f = 1e-6 /s, kh = 3.162e-3, R = 1, K = 1e-4
           (lambda_max = 1.038e-8 /s). Omega is reset to N/C0 with N = 3e6; the linear-noise values are recomputed at that Omega (they scale 1/sqrt(Omega)).
  V1       g = 0, start at the deterministic racemic state, T = 4/f, sd of (L - D)/2 / x over trajectories vs analyse()['sa']. Pass: ratio within 5%.
  V2a      g = Delta*sa/|Ga| for Delta = 1 and 2, same T; P(sign of (L - D) equals sign of Ga*g). Pass: 0.841 +/- 0.04 and 0.977 +/- 0.03.
           The finite-T linear-noise value (mean response and variance at T from the 2x2 odd block, start at the racemic state) is also printed.
  V2b      T = 3/lambda, g = 1/Delta_b_feed(T) per unit g (so Delta_b(T) = 1; also 2/..), PRIMARY statistic = sign of the unstable-mode projection
           w = lv.(a, b) (the object Delta_b describes); sign of (L - D) is printed as well and does not decide the verdict. Pass: |P - Phi(Delta_b(T))| <= 0.05.
           A g = 0 population gives sd(w) vs the linear-noise value (printed, not scored). |a|/x at T is printed as a linearity check.
  M        20000 per population in V1/V2a; 4000 per population in V2b and V2b' (P11 asks for >= 4000). dt * (fastest 5x5 rate) = 0.01 (V1/V2a), 0.05 (V2b, V2b').
  Cores    3 (multiprocessing.Pool(3)). Seeds: numpy SeedSequence(20261021) children.
  V2b'     ADDED AFTER THE SMOKE RUN, NOT PART OF P11: same point and T as V2b but N = 1e12, where the linear regime holds. Rules Q1-Q3 are in
           arm_e_r6_langevin_predictions.txt (sha256 4f4b5e3051eaa12d1b64a7959ad2bb083fe337145e705bf99b21442ebdaafaa7, written Wed Oct 21 01:37:55 PDT 2026,
           asserted below). Reason: the smoke showed N = 3e6 saturates by T (see below), so V2b alone does not test the linear accumulated-noise formula.
DEBUG RUNS BEFORE THE FINAL RUN (disclosed): `python arm_e_r6_langevin.py --smoke`, M = 200 per population, no result read as a verdict, output left in the scratchpad (not a deliverable):
    V1 ratio 1.03; V2a P(sign) 0.785 (Delta 1), 0.990 (Delta 2); V2b P(sign of w) 0.845, 0.970 but median |L - D|/(2x) 0.107 and sd(w) 0.33 of the
    linear value, i.e. N = 3e6 is saturated at T = 3/lambda (ee_sat 0.129, arm_e_r6.log:104). Led to V2b' and to the 3-core pool. The first smoke also had a
    comment claiming ee_sat ~ 0.03; that was a guess and is replaced by the logged 0.129.
    A second smoke (same flags, same seeds, so V1/V2a/V2b numbers identical) exercised the V2b' code path: Q1-Q3 all passed at M = 200; not a verdict either.
    Both smokes were run after arm_e_r6_predictions.txt (P11) was frozen; the addendum file was frozen between them and the second smoke.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python arm_e_r6_langevin.py > arm_e_r6_langevin.log
"""
import hashlib
import json
import sys
import time
from math import erf, sqrt
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.linalg import eig, expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arm_e_r6_analytic as A

PRED_SHA = "235e93aa8a0e4e93bd8d68cecd97befa58eae407b22873ddb90d0f616ca3b9ab"   # arm_e_r6_predictions.txt (P11 is in it)
ADD_SHA = "4f4b5e3051eaa12d1b64a7959ad2bb083fe337145e705bf99b21442ebdaafaa7"    # arm_e_r6_langevin_predictions.txt (V2b')
N_MOL = 3e6
Phi = lambda d: 0.5 * (1 + erf(d / sqrt(2)))
SMOKE = "--smoke" in sys.argv


def props_v(n, p, g):
    """Vectorised channel rates (M, 14) in M/s, same list as arm_e_r6_analytic.props."""
    L, D, A_, LL, DD = n.T
    rL, rD = p["k1"] * LL + p["k2"] * DD + p["ku"], p["k1"] * DD + p["k2"] * LL + p["ku"]
    f, h, kh = p["f"], p["f"] * p["C0"] / 2, p["kh"]
    return np.stack([h * (1 + g), h * (1 - g), f * L, f * D, f * A_, f * LL, f * DD, kh * L * L, kh * D * D, 2 * p["R"] * kh * L * D,
                     rL * L, rL * p["K"] * A_, rD * D, rD * p["K"] * A_], 1)


def cle(p, Om, g, T, dt_rho, seed, n0):
    """Euler-Maruyama of dn = nu^T a dt + nu^T sqrt(a dt / Om) xi; g = per-trajectory feed bias. Returns n(T) (M, 5)."""
    rho = np.abs(np.linalg.eigvals(A.jac(n0, p))).max()
    steps = int(np.ceil(T * rho / dt_rho))
    dt = T / steps
    rng = np.random.default_rng(seed)
    n = np.tile(n0, (len(g), 1))
    for _ in range(steps):
        a = np.maximum(props_v(n, p, g), 0.0)
        n += (a * dt + np.sqrt(a * dt / Om) * rng.standard_normal(a.shape)) @ A.NU
    return n, steps


def odd(p, n0):
    """Odd 2x2 block Jo, diffusion Do (M/s), feed response ho (per unit g), from the 5-species objects."""
    J5 = A.jac(n0, p)
    D5 = A.NU.T @ (A.props(n0, p)[:, None] * A.NU)
    h5 = np.zeros(5)
    h5[:2] = p["f"] * p["C0"] / 2 * np.array([1, -1])
    return A.BP @ J5 @ A.B, A.BP @ D5 @ A.BP.T, A.BP @ h5


def part_a(seed):
    R = json.load(open(HERE / "arm_e_r6_results.json"))["records"]["a|feed_a|E"][0]
    p, Om = A.params(R), N_MOL / R["C0"]
    n0 = A.racemic(p)
    r = A.analyse(p, 4 / p["f"], Om)
    x, sa, Ga = r["x"], r["sa"], r["comp"]["Ga"]
    Jo, Do, ho = odd(p, n0)
    T, M = 4 / p["f"], (200 if SMOKE else 20000)
    g = np.repeat([0.0, 1 * sa / abs(Ga), 2 * sa / abs(Ga)], M)
    t0 = time.time()
    n, steps = cle(p, Om, g, T, 0.01, seed, n0)
    a_rel = (n[:, 0] - n[:, 1]) / 2 / x
    out = [f"V1/V2a pond point f = {p['f']:.4g}, kh = {p['kh']:.4g}, R = {p['R']:.0f}, K = {p['K']:.0e}; Omega = {Om:.4g} (N = {Om * p['C0']:.3g}), T = 4/f = {T:.4g} s, "
           f"{steps} steps, M = {M} per population, {time.time() - t0:.0f} s",
           f"  linear noise: x = {x:.4e} M, sa = {sa:.4e}, Ga = {Ga:.5f}, lambda1 = {r['l1']:.4e}, m = {r['m']:.4f}"]
    sd = a_rel[:M].std(ddof=1)
    # finite-T linear noise, start at the racemic state: mean = (e^{Jt} - I) J^{-1} ho g, Sigma(T) = Sinf - e^{JT} Sinf e^{JT}'
    E = expm(Jo * T)
    mean1 = (E - np.eye(2)) @ np.linalg.solve(Jo, ho)
    Sinf = A.lyap(Jo, Do / Om)
    ST = Sinf - E @ Sinf @ E.T
    sdT = np.sqrt(ST[0, 0]) / x
    out.append(f"  V1: simulated sd of (L - D)/(2x) = {sd:.4e}  vs stationary linear noise {sa:.4e}: ratio {sd / sa:.4f}  (finite-T linear noise {sdT:.4e}: ratio {sd / sdT:.4f}); "
               f"mean {a_rel[:M].mean():.2e} (se {sd / np.sqrt(M):.1e})")
    res = {"V1": abs(sd / sa - 1) <= 0.05}
    for k, (d, lo_hi) in enumerate(((1.0, 0.04), (2.0, 0.03)), 1):
        y = a_rel[k * M:(k + 1) * M]
        gk = g[k * M]
        pc = float(np.mean(np.sign(y) == np.sign(Ga * gk)))
        dT = abs(mean1[0] * gk / x) / sdT
        res[f"V2a{k}"] = abs(pc - Phi(d)) <= lo_hi
        out.append(f"  V2a Delta = {d:.0f}: g = {gk:.4e}; P(sign correct) = {pc:.4f} (se {np.sqrt(pc * (1 - pc) / M):.4f}) vs Phi({d:.0f}) = {Phi(d):.4f}, tolerance +/- {lo_hi}; "
                   f"finite-T Delta = {dT:.4f} -> Phi = {Phi(dT):.4f}; mean response {y.mean() / (Ga * gk):.4f} of the stationary value")
    return out, res


def part_b(seed, N, lin):
    """Unstable pond (b) point at N molecules. lin=False: P11 V2b (scored on the sign only). lin=True: V2b' (Q1-Q3 of the addendum)."""
    R = json.load(open(HERE / "arm_e_r6_results.json"))["records"]["b|feed_strict|U"][0]
    p, Om = A.params(R), N / R["C0"]
    n0 = A.racemic(p)
    lam = A.analyse(p, 1.0, Om)["l1"]
    T = 3 / lam
    r = A.analyse(p, T, Om)
    x, db = r["x"], r["Delta_b_feed"]
    Jo, Do, ho = odd(p, n0)
    w_, VL = eig(Jo, left=True, right=False)
    k = int(np.argmax(w_.real))
    lv = VL[:, k].real
    assert np.isclose(w_[k].real, lam, rtol=1e-6)
    sd_w = np.sqrt((lv @ Do @ lv) * np.expm1(2 * lam * T) / (2 * lam * Om))
    mean_w1 = (lv @ ho) * np.expm1(lam * T) / lam                    # per unit g
    mean_a1 = ((expm(Jo * T) - np.eye(2)) @ np.linalg.solve(Jo, ho))[0]   # linear mean of a per unit g
    assert np.isclose(abs(mean_w1) / sd_w, db, rtol=1e-6), (mean_w1 / sd_w, db)
    M = 200 if SMOKE else 4000
    g = np.repeat([0.0, 1 / db, 2 / db], M)
    t0 = time.time()
    n, steps = cle(p, Om, g, T, 0.05, seed, n0)
    a = (n[:, 0] - n[:, 1]) / 2
    w = lv[0] * a + lv[1] * (n[:, 3] - n[:, 4]) / 2
    tag = "V2b'" if lin else "V2b"
    out = [f"{tag} pond point F = {R['F']:.2f}, f = {p['f']:.4g}, kh = {p['kh']:.4g}, R = {p['R']:.0f}, K = {p['K']:.0e}; Omega = {Om:.4g} (N = {Om * p['C0']:.3g}), "
           f"lambda = {lam:.4e} /s, T = 3/lambda = {T:.4g} s, {steps} steps, M = {M} per population, {time.time() - t0:.0f} s",
           f"  linear noise: Delta_b(T) per unit g = {db:.4e}; sd(w) = {sd_w:.4e}"]
    s0, p95 = w[:M].std(ddof=1), np.percentile(np.abs(a[:M]), 95) / x
    out.append(f"  g = 0: simulated sd(w) = {s0:.4e} vs linear {sd_w:.4e}: ratio {s0 / sd_w:.4f}; median |L - D|/(2x) = {np.median(np.abs(a[:M])) / x:.3e}, "
               f"95th percentile {p95:.3e} (linearity check; ee_sat of this point is 0.129, arm_e_r6.log:104)")
    res = {}
    if lin:
        res["Q1"], res["Q3"] = abs(s0 / sd_w - 1) <= 0.05, p95 <= 0.02
    for kk, d, tol in ((1, 1.0, 0.04), (2, 2.0, 0.03)):
        sl = slice(kk * M, (kk + 1) * M)
        gk = g[kk * M]
        pw = float(np.mean(np.sign(w[sl]) == np.sign((lv @ ho) * gk)))
        pa = float(np.mean(np.sign(a[sl]) == np.sign(mean_a1 * gk)))
        out.append(f"  {tag} Delta_b = {d:.0f}: g = {gk:.4e}; P(sign of w) = {pw:.4f} (se {np.sqrt(pw * (1 - pw) / M):.4f}) vs Phi = {Phi(d):.4f}; "
                   f"mean w / linear mean = {w[sl].mean() / (mean_w1 * gk):.3f}; (P(sign of L - D) = {pa:.4f}, not scored)")
        if lin:
            res[f"Q2{'ab'[kk - 1]}"] = abs(pw - Phi(d)) <= tol
        elif kk == 1:
            res["V2b"] = abs(pw - Phi(d)) <= 0.05
    return out, res


def self_check():
    """Vectorised rates equal arm_e_r6_analytic.props, the racemic state is a fixed point, the sampler keeps it in the mean (no-noise limit)."""
    R = json.load(open(HERE / "arm_e_r6_results.json"))["records"]["a|feed_a|E"][0]
    p = A.params(R)
    n0 = A.racemic(p)
    assert np.allclose(props_v(n0[None, :], p, np.array([0.0]))[0], A.props(n0, p), rtol=1e-13)
    assert np.abs(A.drift(n0, p)).max() < 1e-12 * p["C0"] * p["f"]
    n, _ = cle(p, 1e30, np.zeros(3), 1 / p["f"], 0.05, 1, n0)       # Omega -> infinity: deterministic, stays at the fixed point
    assert np.allclose(n, n0, rtol=1e-6)


def main():
    for f_, want in (("arm_e_r6_predictions.txt", PRED_SHA), ("arm_e_r6_langevin_predictions.txt", ADD_SHA)):
        sha = hashlib.sha256((HERE / f_).read_bytes()).hexdigest()
        assert sha == want, (f_, sha)
    self_check()
    s = np.random.SeedSequence(20261021).spawn(3)
    print(f"pre-registration files: arm_e_r6_predictions.txt sha256 {PRED_SHA}; arm_e_r6_langevin_predictions.txt sha256 {ADD_SHA}\n"
          "self-check passed: vectorised channels = arm_e_r6_analytic.props, racemic fixed point, deterministic limit\n" + ("SMOKE RUN (M = 200): not a result\n" if SMOKE else ""), flush=True)
    t0 = time.time()
    with Pool(3) as pool:
        jobs = [pool.apply_async(part_a, (s[0],)), pool.apply_async(part_b, (s[1], N_MOL, False)), pool.apply_async(part_b, (s[2], 1e12, True))]
        (oa, resa), (ob, resb), (oc, resc) = [j.get() for j in jobs]
    for line in oa + ob + oc:
        print(line)
    res = {**resa, **resb}
    print(f"\nP11 {'held' if all(res.values()) else 'FAILED'}: " + ", ".join(f"{k} {'pass' if v else 'FAIL'}" for k, v in res.items()))
    print(f"Q (addendum, V2b') {'held' if all(resc.values()) else 'FAILED'}: " + ", ".join(f"{k} {'pass' if v else 'FAIL'}" for k, v in resc.items()) + f"   (wall {time.time() - t0:.0f} s, 3 cores)")


if __name__ == "__main__":
    main()
