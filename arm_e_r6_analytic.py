"""Arm E R6, independent check (P10) of arm_e_r6.py: the full 5-species chemical-Langevin calculation from the reaction channels.

Species L, D, A, LL, DD (M). Channels (rate in M/s; stoichiometry NU), Higgs & Blackmond 2025 destructive mechanism:
  feed L f C0/2 (1+g); feed D f C0/2 (1-g); outflow f n for every species; L+L->LL kh L^2; D+D->DD kh D^2; L+D->LD (sink) 2 R kh L D;
  L->A r_L L; A->L r_L K A; D->A r_D D; A->D r_D K A;   r_L = k1 LL + k2 DD + ku, r_D = k1 DD + k2 LL + ku.
Nothing is imported from arm_e_rates.py or arm_e_r6.py and no formula of theirs is reused: racemic state = brentq on the channel drift;
Jacobian = complex-step derivative of the drift; odd block (a, b) = ((L-D)/2, (LL-DD)/2) = projection of the 5x5 Jacobian (closure checked);
noise = 5x5 Lyapunov equation J S + S J' + D5/Omega = 0 (scipy Bartels-Stewart); seed propagation and the unstable mode = LAPACK eig of
the 2x2 odd block (left and right eigenvectors). Variants and definitions are those of arm_e_r6.py's docstring and arm_e_r6_predictions.txt.
Constants KB_SUM, KU are re-typed from Yu 2024 SI Table 5 (reverse rates 0.91 and 4.04 mM/h; the same arithmetic as arm_e_rates.py:75-76),
so P (the racemic breakdown rate) is also compared.

PRE-REGISTRATION (P10). Inputs read from arm_e_r6_results.json in --predict mode: ONLY the grid point of each optimum record
(setting index, C0, F, s, f, kh, R, K, tau, Omega). Disclosure: arm_e_r6.log had been read before this file was written, so the values
are not blind; they are computed by code that shares no formula with arm_e_r6.py, which is what the check needs. The predictions file
arm_e_r6_analytic_predictions.txt (every number of this calculation at every unique optimum point) was written by `--predict` BEFORE
the comparison run; the comparison run asserts the file hash, re-computes the numbers (they must reproduce the file to 1e-9) and only
then reads the other columns of arm_e_r6_results.json. Pass rule (P10): max relative deviation <= 1e-5 over J11..J22, D, x, y, A, P,
margin, lambda1,2 and every e_min component (this covers G, Sigma_aa and Sigma_bb through 2 s/G, and Delta_b through u_feed = 2/Delta_b).
    file sha256 efd13dc79a9342445fadcc98ccfac6a7cbd8058a0873a696e7952d23e301edad, written Wed Oct 21 01:21:39 PDT 2026 (shasum at 01:21:40)
    (one failed --predict attempt at 01:21:24 stopped at the input check before any number was computed: ocean.py uses NA = 6.022e23, I had typed the CODATA value)
After the first comparison run (which printed the table, then crashed formatting None sa/sb of the unstable points in the "extras" lines) I
changed ONLY: the extras print (None-safe), and added `--json NAME` (compare against another results file; POST HOC, not P10). No number
this file computes changed (the 1465 predicted numbers still reproduce to 1e-9). The table showed seedm_a deviating (see arm_e_r6_analytic.log):
P10 is scored on that first table. The cause is a cancellation in arm_e_r6.seed_stable, fixed in arm_e_r6_fix.py (post hoc); the second run of
this file against arm_e_r6_fix_results.json is reported separately (arm_e_r6_analytic_fix.log) and does not rescore P10.

Run: python3 arm_e_r6_analytic.py --predict > arm_e_r6_analytic_predictions.txt   (needs numpy, scipy)
     UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python arm_e_r6_analytic.py > arm_e_r6_analytic.log
     ... python arm_e_r6_analytic.py --json arm_e_r6_fix_results.json > arm_e_r6_analytic_fix.log      (post hoc)
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig, solve_continuous_lyapunov
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
PRED_SHA = "efd13dc79a9342445fadcc98ccfac6a7cbd8058a0873a696e7952d23e301edad"
KB_SUM = (4.04 - 0.91) * 1e-3 / 3600 / (0.100 * 0.050)   # k1b + k2b, /M/s
KU = 0.91e-3 / 3600 / 0.100                               # uncatalysed breakdown, /s
NA, VOL = 6.022e23, (np.pi, 1e7, 1.332e18)               # ocean.py:25 NA; scorecard.SETTINGS volumes (m^3)
NU = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [-1, 0, 0, 0, 0], [0, -1, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, -1],
               [-2, 0, 0, 1, 0], [0, -2, 0, 0, 1], [-1, -1, 0, 0, 0], [-1, 0, 1, 0, 0], [1, 0, -1, 0, 0], [0, -1, 1, 0, 0], [0, 1, -1, 0, 0]], float)
B = np.array([[1, 0], [-1, 0], [0, 0], [0, 1], [0, -1]], float)    # odd subspace: delta = B (a, b)
E = np.array([[1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]], float)  # even subspace
BP, EP = np.linalg.pinv(B), np.linalg.pinv(E)
SCAL = ("x", "y", "A", "P", "m", "l1", "l2", "J11", "J12", "J21", "J22", "Daa", "Dab", "Dbb", "cat")


def props(n, p, g=0.0):
    L, D, A, LL, DD = n
    rL, rD = p["k1"] * LL + p["k2"] * DD + p["ku"], p["k1"] * DD + p["k2"] * LL + p["ku"]
    f, h, kh = p["f"], p["f"] * p["C0"] / 2, p["kh"]
    return np.array([h * (1 + g), h * (1 - g), f * L, f * D, f * A, f * LL, f * DD, kh * L * L, kh * D * D, 2 * p["R"] * kh * L * D,
                     rL * L, rL * p["K"] * A, rD * D, rD * p["K"] * A])


def drift(n, p, g=0.0):
    return NU.T @ props(n, p, g)


def racemic(p):
    def fill(L):
        LL = drift([L, L, 0, 0, 0], p)[3] / p["f"]                                   # d(LL)/dt = kh L^2 - f LL
        d0 = drift([L, L, 0, LL, LL], p)[2]
        d1 = drift([L, L, 1, LL, LL], p)[2] - d0                                    # d(A)/dt is linear in A
        return np.array([L, L, -d0 / d1, LL, LL])
    L = brentq(lambda L: drift(fill(L), p)[0], 0.0, p["C0"] / 2, xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=500)
    return fill(L)


def jac(n, p):
    J = np.empty((5, 5))
    for j in range(5):
        h = 1e-20 * max(n[j], 1e-6 * n.max())
        z = n.astype(complex)
        z[j] += 1j * h
        J[:, j] = drift(z, p).imag / h
    return J


def lyap(J, Dm):
    st, sd = 1 / np.abs(J).max(), 1 / np.abs(Dm).max()
    return solve_continuous_lyapunov(J * st, -Dm * st * sd) / sd


def analyse(p, tau, Om, guard=5.0):
    """Everything arm_e_r6 reports at one grid point, from the 5-species channels."""
    n = racemic(p)
    J5 = jac(n, p)
    a = props(n, p)
    D5 = (NU.T * a) @ NU                                   # sum_j a_j nu_j nu_j' (per Omega)
    Jo, Do, Je = BP @ J5 @ B, BP @ D5 @ BP.T, EP @ J5 @ E
    scale = np.abs(J5).max()
    clos = max(np.abs(J5 @ B - B @ Jo).max(), np.abs(J5 @ E - E @ Je).max()) / scale
    assert clos < 1e-9, clos                               # odd and even sectors decouple
    x, A, y = n[0], n[2], n[3]
    lo = np.sort(np.linalg.eigvals(Jo).real)[::-1]
    l1, l2 = lo
    ev_max = np.linalg.eigvals(Je).real.max()
    h5 = np.zeros(5)
    h5[:2] = p["f"] * p["C0"] / 2 * np.array([1, -1])      # d drift / d g
    ho = BP @ h5
    CA, CB = BP[0], BP[1]
    out = dict(x=x, y=y, A=A, P=(p["k1"] + p["k2"]) * y + p["ku"], l1=l1, l2=l2, J11=Jo[0, 0], J12=Jo[0, 1], J21=Jo[1, 0], J22=Jo[1, 1],
               Daa=Do[0, 0], Dab=Do[0, 1], Dbb=Do[1, 1], even_max=ev_max, closure=clos)
    inf, comp = float("inf"), {}
    with np.errstate(all="ignore"):
        out["m"] = Jo[0, 1] * Jo[1, 0] / (Jo[0, 0] * Jo[1, 1])
        out["cat"] = (p["k1"] + p["k2"]) * y / out["P"]
        st = l1 < 0
        ok = st and p["f"] * tau >= guard and -l1 * tau >= guard
        Ga = Gb = sa = sb = float("nan")
        if st:
            dl = -np.linalg.solve(J5, h5)
            Ga, Gb = CA @ dl / x, CB @ dl / y
            S = lyap(J5, D5 / Om)
            assert np.abs(J5 @ S + S @ J5.T + D5 / Om).max() <= 1e-6 * np.abs(D5 / Om).max()
            sa, sb = np.sqrt(CA @ S @ CA) / x, np.sqrt(CB @ S @ CB) / y
            w, V = np.linalg.eig(Jo)
            q = lambda w0: (V @ (np.exp(w.real * tau) * np.linalg.solve(V, w0))).real
            qs = {"feed": (Ga, Gb), "seed": (q(np.array([x, y]))[0] / x, q(np.array([x, y]))[1] / y),
                  "seedm": (q(np.array([x, 0.0]))[0] / x, q(np.array([x, 0.0]))[1] / y)}
            for k, (qa, qb) in qs.items():
                for nm, qq, s_ in (("a", qa, sa), ("b", qb, sb)):
                    lit = 2 * s_ / abs(qq)
                    comp[f"{k}_{nm}"], comp[f"{k}S_{nm}"] = (lit, max(lit, 0.01 / abs(qq))) if ok else (inf, inf)
        else:
            for k in ("feed", "seed", "seedm"):
                for nm in "ab":
                    comp[f"{k}_{nm}"] = comp[f"{k}S_{nm}"] = inf
        out.update(sa=sa, sb=sb, Saa_Om=sa * sa * x * x * Om, Sbb_Om=sb * sb * y * y * Om)    # Sigma/Omega convention of arm_e_r6.lyap
        un = {k: inf for k in ("u_feed", "u_feedS", "u_seed", "u_seedS", "u_seedm", "u_seedmS")}
        if l1 > 0:                                         # unstable odd mode: accumulated noise vs mean, growth gate to 1% output
            w, VL, VR = eig(Jo, left=True, right=True)
            k = int(np.argmax(w.real))
            lam, lv, rv = w[k].real, VL[:, k].real, VR[:, k].real
            rv = rv / (lv @ rv)
            lt = lam * tau
            sd = np.sqrt((lv @ Do @ lv) * (-np.expm1(-2 * lt)) / (2 * lam * Om))
            gain = max(abs(rv[0]) / x, abs(rv[1]) / y)
            for nm, c in (("feed", (lv @ ho) * (-np.expm1(-lt)) / lam), ("seed", lv @ np.array([x, y])), ("seedm", lv @ np.array([x, 0.0]))):
                d = 2 * sd / abs(c)
                gate = np.exp(min(np.log(0.01) - lt - np.log(abs(c) * gain), 700.0))   # 700: exp overflow cap
                un[f"u_{nm}"], un[f"u_{nm}S"] = d, max(d, gate)
            out["Delta_b_feed"] = abs((lv @ ho) * (-np.expm1(-lt)) / lam) / sd
        comp.update(un)
        for k in ("feed", "seed", "seedm"):
            comp[f"{k}_lit"] = np.fmin(np.fmin(comp[f"{k}_a"], comp[f"{k}_b"]), comp[f"u_{k}"])
            comp[f"{k}_strict"] = np.fmin(np.fmin(comp[f"{k}S_a"], comp[f"{k}S_b"]), comp[f"u_{k}S"])
    out["comp"] = {k: (float(v) if np.isfinite(v) else None) for k, v in comp.items()}
    out["comp"]["Ga"], out["comp"]["Gb"] = (float(Ga) if st and np.isfinite(Ga) else None), (float(Gb) if st and np.isfinite(Gb) else None)
    return {k: (v if isinstance(v, dict) else (float(v) if np.isfinite(v) else None)) for k, v in out.items()}


# ------------------------------------------------------------------ inputs
def params(rec):
    F, s = rec["F"], rec["s"]
    return dict(f=rec["f"], C0=rec["C0"], kh=rec["kh"], R=rec["R"], K=rec["K"], k1=F * KB_SUM / (1 + s), k2=F * KB_SUM * s / (1 + s), ku=KU)


def points(name="arm_e_r6_results.json"):
    """Unique optimum grid points of a results json: id -> (record used for parameters, [record labels]). Parameters only."""
    R = json.load(open(HERE / name))["records"]
    pts, idx = {}, {}
    for key, lst in R.items():
        for i, rec in enumerate(lst):
            if rec is None:
                continue
            assert np.isclose(rec["Om"], VOL[i] * 1000 * NA, rtol=1e-12)      # i = 0, 1, 2 = pond, lake, ocean
            pid = f"S{i}_F{rec['F']:.10g}_s{rec['s']:.6g}_f{rec['f']:.10g}_kh{rec['kh']:.10g}_R{rec['R']:.0f}_K{rec['K']:.10g}"
            pts.setdefault(pid, (rec, []))[1].append(key)
            idx[(key, i)] = pid
    return pts, idx


def flat(r):
    d = {k: v for k, v in r.items() if k != "comp"}
    d.update({f"comp.{k}": v for k, v in r["comp"].items()})
    return d


def predict():
    lines = ["# arm_e_r6_analytic.py pre-registration (--predict): independent 5-species values at every unique optimum grid point.",
             "# Only the parameters (setting, F, s, f, kh, R, K, tau, Omega) of arm_e_r6_results.json were read. none = not finite / not applicable.",
             f"# written: {time.strftime('%a %b %d %H:%M:%S %Z %Y')}"]
    for pid, (rec, labels) in points()[0].items():
        r = analyse(params(rec), rec["tau"], rec["Om"])
        lines.append(f"# {pid} used by: {', '.join(labels)}")
        lines += [f"{pid}.{k} = {'none' if v is None else repr(float(v))}" for k, v in flat(r).items()]
    return "\n".join(lines)


def self_check():
    """No reaction (kh = k1 = k2 = ku = 0): L and D are independent birth-death pools; G = 1, sd = 1/sqrt(N), seed decays as e^(-f tau)."""
    p = dict(f=1e-6, C0=5e-3, kh=0.0, R=1.0, K=1e-3, k1=0.0, k2=0.0, ku=0.0)
    Om, tau = 1e20, 1e7
    N = Om * p["C0"]
    r = analyse(p, tau, Om)
    c = r["comp"]
    assert np.isclose(r["x"], 2.5e-3, rtol=1e-12) and np.isclose(c["Ga"], 1.0, rtol=1e-10) and np.isclose(r["sa"], 1 / np.sqrt(N), rtol=1e-9)
    assert np.isclose(c["feed_a"], 2 / np.sqrt(N), rtol=1e-9) and np.isclose(c["seed_a"], 2 / np.sqrt(N) * np.exp(p["f"] * tau), rtol=1e-8)
    assert np.isclose(r["l1"], -p["f"], rtol=1e-9) and r["even_max"] < 0


def main():
    sha = hashlib.sha256((HERE / "arm_e_r6_analytic_predictions.txt").read_bytes()).hexdigest()
    assert sha == PRED_SHA, f"analytic predictions file changed: {sha}"
    self_check()
    log = lambda *a: print(" ".join(str(x) for x in a), flush=True)
    log(f"pre-registration file arm_e_r6_analytic_predictions.txt: sha256 {sha}")
    log("self-check passed: no-reaction case (G = 1, sd = 1/sqrt(N), e_min = 2/sqrt(N), seed decays as e^(-f tau), lambda = -f)\n")
    pred = {}
    for line in (HERE / "arm_e_r6_analytic_predictions.txt").read_text().splitlines():
        if line and not line.startswith("#"):
            k, v = line.split(" = ")
            pred[k] = None if v == "none" else float(v)
    jname = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else "arm_e_r6_results.json"
    pts, idx = points(jname)
    live, nrep, nnew = {}, 0, 0
    for pid, (rec, _) in pts.items():
        live[pid] = analyse(params(rec), rec["tau"], rec["Om"])
        if f"{pid}.x" not in pred:                         # grid point not in the pre-registered file (only with --json)
            nnew += 1
            continue
        for k, v in flat(live[pid]).items():
            a = pred[f"{pid}.{k}"]
            assert (a is None) == (v is None) and (a is None or np.isclose(a, v, rtol=1e-9, atol=0)), (pid, k, a, v)   # reproduces the file
            nrep += 1
    log(f"{len(pts)} unique optimum points ({nnew} not in the pre-registered file); {nrep} numbers recomputed and reproduced from the predictions file to rel 1e-9")
    log(f"comparison target: {jname}" + ("" if jname == "arm_e_r6_results.json" else "   (POST HOC, not the pre-registered P10 test)") + "\n")
    R = json.load(open(HERE / jname))["records"]
    log(f"{'record':26} {'setting':6} {'e_min arm_e_r6':>14} {'e_min analytic':>14} {'max rel dev':>11}  worst quantity")
    worst_all, ncmp = 0.0, 0
    for key, lst in R.items():
        for i, rec in enumerate(lst):
            if rec is None:
                log(f"{key:26} {('pond', 'lake', 'ocean')[i]:6} none (F_c = inf)")
                continue
            lv = live[idx[(key, i)]]
            pairs = [(k, lv[k], rec[k]) for k in SCAL] + [(f"comp.{k}", lv["comp"][k], v) for k, v in rec["comp"].items()]
            dev, who = 0.0, ""
            for k, a, b in pairs:
                ncmp += 1
                if (a is None) != (b is None):
                    dev, who = float("inf"), f"{k}: finite/none mismatch ({a} vs {b})"
                    break
                if a is not None and abs(a - b) / max(abs(b), 1e-300) > dev:
                    dev, who = abs(a - b) / max(abs(b), 1e-300), k
            worst_all = max(worst_all, dev)
            ml = lv["comp"][rec["key"]]
            log(f"{key:26} {('pond', 'lake', 'ocean')[i]:6} {rec['e_min']:14.6e} {(ml if ml is not None else float('nan')):14.6e} {dev:11.2e}  {who}")
    log(f"\ncompared {ncmp} numbers against arm_e_r6_results.json; max relative deviation {worst_all:.3e}")
    log("independent extras (not in the JSON), per unique point: even-sector max Re(lambda), odd/even closure, readout sd, Sigma/Omega, Delta_b")
    fm = lambda v, f: "none" if v is None else format(v, f)
    for pid, lv in live.items():
        log(f"  {pid}: even_max {fm(lv['even_max'], '.2e')} /s, closure {fm(lv['closure'], '.1e')}, sa {fm(lv['sa'], '.3e')} sb {fm(lv['sb'], '.3e')}, "
            f"Sigma_aa/Omega {fm(lv['Saa_Om'], '.3e')} Sigma_bb/Omega {fm(lv['Sbb_Om'], '.3e')}" + (f", Delta_b(feed, per unit g) {lv['Delta_b_feed']:.3e}" if "Delta_b_feed" in lv else ""))
    ok = worst_all <= 1e-5
    tag = "P10" if jname == "arm_e_r6_results.json" else f"post hoc check vs {jname} (not P10)"
    log(f"\n{tag} {'held' if ok else 'FAILED'}: max relative deviation {worst_all:.3e} (rule <= 1e-5) over {ncmp} numbers")
    log(f"even sector stable at every point: {all(lv['even_max'] < 0 for lv in live.values())}")


if __name__ == "__main__":
    if "--predict" in sys.argv:
        self_check()
        print(predict())
    else:
        main()
