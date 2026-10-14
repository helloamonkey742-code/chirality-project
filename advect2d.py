"""advect2d: does Eq. 14 survive when eddies are explicit advection instead of an eddy diffusivity? (H29;
journals/ija_draft.md Section 2.3, "Eddy diffusion" paragraph, 2026-10-13.)

THE UNTESTED ASSUMPTION. Section 2.3 treats ocean eddies as an eddy diffusivity D_h in Eq. 14 and says: "This needs the
wall width sqrt(D/k) to exceed the eddy size. I have not shown that it does for every sampled input; larger eddies could
fold and stretch walls in ways Eq. 14 does not include." This script puts the folding and stretching in: the same
stochastic field equation in 2D with an explicit incompressible chaotic flow u,

    d(alpha) = (lam alpha - alpha^3 + g + D_mol lap(alpha) - u . grad(alpha)) dt + sqrt(eps / h^2) dW,
    lam = clip(gamma t, -1, 1), gamma = 0.01, eps = 1e-4 per unit area (same conventions as kz_patches.py / front2d.py),

and compares it with the draft's assumption: the same grid, no flow, D = the measured D_eff of the flow.

MODEL. Periodic box L = 64 (units: the reaction rate k = 1 at lam = 1, cell size 1 in the no-flow arms; h = 0.5 in the
flow arms), D_eff target 8, so the draft's coherence length is xi = sqrt(D_eff) gamma^-1/4 = 8.9 (L = 7.2 xi) and the
wall width is w = sqrt(D_eff / k) = 2.83. Flow: the random-phase alternating sine shear (a kicked, discrete-time version of
Pierrehumbert's alternating flow). Every Tk, alternately along x and along y, the field is displaced exactly along
streamlines, alpha(x, y) -> alpha(x - A sin(2 pi n y / L + phi), y), with a new uniform random phase phi per kick per
run (cubic Lagrange interpolation along the shear direction; it is circulant, so it conserves the mean of alpha exactly),
and the reaction-diffusion-noise step (explicit Euler, dt = 0.05) runs between kicks (operator splitting). Tracers move
as independent random walks of variance A^2 / 2 per kick, so the eddy diffusivity is D_eddy = A^2 / (8 Tk) exactly;
D_mol = 1 is the molecular part, and A is set so D_mol + D_eddy = 8:  A = sqrt(56 Tk). Eddy size = 1 / k_flow =
L / (2 pi n); the shear strength per kick is kappa = A k_flow (stretching; kappa >~ 1 is chaotic).

    regime            n    eddy l   l / w   Tk     A      kappa   U = A/Tk   Da = k l / U
    S small eddies    12   0.85     0.30    0.10   2.37   2.8     23.7       0.04
    M matched         4    2.55     0.90    0.25   3.74   1.5     15.0       0.17
    B big eddies      1    10.2     3.6     4.0    14.97  1.5     3.7        2.7

(Eddy size and Da are not free choices at fixed D_eff = U l: l / w = sqrt(Da) with w = sqrt(D_eff/k). So "eddy size vs
wall width" is the same axis as the Damkohler number Da = k l / U, the ratio of the eddy turnover time to the reaction
time.) The three regimes keep D_eff fixed and move l / w across 0.3, 0.9, 3.6. Flow arms use h = 0.5 (dt = 0.05,
D_mol dt / h^2 = 0.2, the Euler limit is 0.25) and 128 x 128 cells; all metrics are computed after 2 x 2 block averaging
onto a 64 x 64 grid so that they are on the same lattice as the no-flow arm (h = 1, dt = 0.025, D = 8). Further reference
arm: no flow, D = D_mol = 1 (h = 1, dt = 0.05), which shows what the flow arms would do if the eddies did nothing.

D_eff MEASUREMENT (same kick code). (i) Lagrangian: 400 flow realisations x 50 tracers, <dx^2 + dy^2> / (4 t), eddy part
(exact map). (ii) Eulerian: a long-wave mode cos(K x) of a passive tracer (advection + D_mol only, same kick and diffusion
code and h) in a box of 4 L (B) or L (S, M); its amplitude decays as exp(-K^2 D_eff t), which is the characteristic
function of the tracer displacement and so is tracer variance growth; D_eff(Eulerian) = -ln(A(T) / A(0)) / (K^2 T),
T = 1 / (K^2 8), 4 realisations. The Eulerian number includes D_mol and any numerical diffusion from the interpolation.
The no-flow comparison arm uses D = 8.0 (the design value); P1 says how close the measurement is to that.

OUTCOMES. Delta_total = 0.5, 1, 2 (g set from the whole volume L^2 as kz_patches.g_for, i.e. the eddy-diffusion
prediction), 64 seeds per Delta per arm. Score = majority hand of sign(alpha) on the 64 x 64 grid (ties count 1/2) at the
end of the sweep (lam = 1) = "end", at lam = 0.4 = "early", and at 80 time units after the end of the sweep (flow and
diffusion running, lam = 1) = "final". Fit P = Phi(a Delta) by maximum likelihood (kz_patches.fit_probit): a = 1 is the
draft's whole-volume Eq. 11; the no-flow D_eff arm gives the reference a (kz_patches' 2D pooled c = 0.98 +- 0.05 for
D = 1; theory 0.96). Walls: wall length / L (sign changes of nearest-neighbour pairs on the 64 x 64 grid, divided by L^2
and multiplied by L; a single straight wall across a periodic box counts 2), number of domains with area >= 8
(periodic, 4-connected), at t_after = 0, 5, 10, 20, 40, 80; the first time (every 2 t.u.) the whole grid has one sign
("time to single domain", censored at 80), pooled over the three Deltas (g is ~5e-5, too small to move walls).
Precision: 192 runs per arm gives SE(a) ~ 0.1, so P3 can only exclude differences larger than ~0.2 in a.

SMOKE RUNS (disclosed). Before writing these predictions I ran nothing. After the first commit I ran only tiny
debugging smokes (L = 16, 4 seeds, gamma = 0.5) that print step times, shapes and NaN/finite checks; no outcome
statistic (P, a, D_eff) was inspected from them; the smoke printed two wall-length and two domain-count values per
arm from 4 runs at L = 16 and gamma = 0.5 (values of order 5-12 and 2-5), which say nothing about L = 64. Timing runs
(120 steps at full size, gamma = 0.5, output discarded except elapsed seconds) set only the seed count and order of
jobs. Design numbers were chosen from the algebra above and the cost estimates only.

PREDICTIONS (written 2026-10-13, before any result-producing run):
 P1. Flow calibration: the Lagrangian eddy part is within 5% of 7.0 in S, M, B; the Eulerian total D_eff is within 15%
     of 8.0 in S and M and within 25% in B.
 P2. No-flow control: kz_patches.check() (uniform limit) passes within 3 SE, and a(no flow, D = 8, end) lies in
     [0.85, 1.15].
 P3. The outcome follows the eddy-diffusion prediction: for each of S, M, B, a(flow, end) is within 2 combined SE of
     a(no flow, D = 8, end).
 P4. Probit shape: in every flow arm P(end) increases with Delta and lies within 2.5 SE of Phi(a Delta) at each Delta.
 P5. N is the whole volume: a(flow, end) lies in [0.8, 1.2] in each of S, M, B.
 P6. Walls at freeze-out (t_after = 0): wall length relative to the no-flow D = 8 arm is in [0.75, 1.25] for S and
     >= 1.2 for B (eddies larger than the wall leave scales below l to molecular D only).
 P7. Walls are not torn away faster than D_eff diffusion carries them: wall length at t_after = 20 relative to the
     no-flow D = 8 arm lies in [0.5, 2] for S, is >= 0.5 for M, and is >= 1.3 for B.
 P8. Time to one domain: the fraction of runs with a single domain by t_after = 80 differs from the no-flow D = 8 arm by
     less than 0.15 for S, and is lower by at least 0.15 for B.
 P9. Monotone in eddy size: wall length at t_after = 20 satisfies S <= M <= B (each step not below the previous by more
     than 1 combined SE).
 P10. In B the pattern at freeze-out is set by molecular D alone: wall length at t_after = 0 for B lies within
     [0.7, 1.4] x the no-flow D = 1 arm.
 P11. Wall time under 20 minutes (nice -n 10, 2 threads).

Usage: OMP_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python advect2d.py > advect2d.log     (smoke: python advect2d.py smoke)

RESULTS (advect2d.log, 2026-10-13):
Run: OMP_NUM_THREADS=2 nice -n 10, 775 s. Disclosure: the first full run crashed in the wall report (KeyError: 5; probes 5
and 10 are not multiples of the 2 t.u. cadence because probe recording sat inside the cadence block). Its outcome table
(the a values) had already printed and I saw it; the results were not saved. I moved the probe block out of the cadence
block (commit "advect2d: fix probe-time bug", predictions unchanged) and reran with the same seeds, which reproduced
identical outcome numbers; advect2d.log is the rerun. P1-P11 above are untouched.
 P1 HELD. Lagrangian eddy part 6.95 / 7.10 / 7.04 (S/M/B, target 7.0); Eulerian total 8.05 / 8.35 / 8.34 (target 8.0;
    +4.4% in M and B).
 P2 HELD. Uniform-limit check passed (+1.13 SE); a(no flow, D = 8) = 0.929 +- 0.112.
 P3 HELD. a(end): S 1.019 +- 0.121, M 0.855 +- 0.107, B 1.103 +- 0.134; differences from the no-flow arm are
    +0.55, -0.47, +1.00 combined SE. (With SE(a) ~ 0.11 this only excludes differences of ~0.2 or more.)
 P4 FAILED, on one point. 8 of 9 points lie within 1.8 SE of Phi(a Delta); M at Delta = 2 gives P = 1.000 (64 of 64) and
    is +5.58 SE only because the SE floor for P = 1 with n = 64 is 0.008. That is an artifact of the floor, not a
    misfit of the probit, but the pre-registered rule says FAILED. P increases with Delta in every arm.
 P5 HELD. a(end): S 1.02, M 0.855, B 1.10, all inside [0.8, 1.2].
 P6 FAILED. Wall(t0) / no-flow D = 8: S 1.39, M 4.66, B 28.9 (no-flow 0.03). The B half (>= 1.2) held by a wide margin;
    the S half (in [0.75, 1.25]) failed on counts of 2 vs 3 non-single runs out of 192 (0.04 +- 0.03 vs 0.03 +- 0.02),
    which is noise at that level.
 P7 HELD. Wall(t20) / no-flow D = 8: S 1.37, M 3.74, B 21.2.
 P8 FAILED. Single domain by t_after = 80: S 0.990, M 0.984, B 0.943, no flow D = 8 0.984, no flow D = 1 0.448. B is
    lower than the D = 8 arm by 0.04, not 0.15, because 71% of B runs were already single at freeze-out.
 P9 HELD. Wall(t20): S 0.043 <= M 0.117 <= B 0.663.
 P10 FAILED. Wall(t0) B / no-flow D = 1 = 0.56 (0.92 vs 1.63), outside [0.7, 1.4]. The big eddies do smooth the
    pattern, about half as much as the molecular-D pattern, but far less than D_eff = 8 does (which leaves 0.03).
 P11 HELD. 775 s.
Post-hoc diagnostics (scratch script on the saved raw results, NOT in this file; labelled exploratory): among runs that
are not single-domain at freeze-out (S 2/192, M 4/192, B 56/192, no flow D = 8 3/192, no flow D = 1 151/192), the
fraction single by t_after = 10 / 20 / 40 / 80 is B 0.14 / 0.25 / 0.50 / 0.80, no flow D = 1 0.08 / 0.11 / 0.17 / 0.30,
M 0 / 0 / 0 / 0.25, S and no flow D = 8 0 throughout (too few runs). The majority hand at the end of the sweep equals
the one 80 t.u. later in 99.5% of B and M runs. For comparison the draft's Eq. 17 with D_eff = 8 puts the single-domain
time at ~17 t.u.
Reading: the OUTCOME (which hand wins) follows the eddy-diffusion prediction in all three regimes. The WALLS do not in
the big-eddy regime (eddies about 3.6 wall widths): walls survive at freeze-out in 29% of runs (1.6% at D_eff) and are
removed at least 4x slower than D_eff diffusion predicts, though much faster than molecular D alone.
Limits: SE(a) ~ 0.11, so differences below ~0.2 in a are invisible; box is 7.2 xi; kicks are operator-split; flow arms
h = 0.5 vs no-flow h = 1; no convergence test in h, dt or Tk; Da = 2.7 is the only slow-eddy case.
"""
import os
import pickle
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.ndimage import label
from scipy.stats import norm

import kz_patches as K

GAMMA, EPS = K.GAMMA, K.EPS
L, D_MOL, D_EFF = 64.0, 1.0, 8.0
DELTAS, SEEDS = (0.5, 1.0, 2.0), 64
PROBES, CADENCE, HOLD, AMIN = (0, 5, 10, 20, 40, 80), 2.0, 80.0, 8.0
FLOWS = {"S": dict(n=12, Tk=0.10), "M": dict(n=4, Tk=0.25), "B": dict(n=1, Tk=4.0)}
for _f in FLOWS.values():
    _f["A"] = float(np.sqrt(8 * (D_EFF - D_MOL) * _f["Tk"]))


def kick(a, s, h, axis):
    """Shear displacement: along the last axis (axis=2) or the middle one (axis=1), cell row r moves by s[run, r]
    (physical length). Cubic Lagrange interpolation; circulant with weights summing to 1, so the mean is conserved."""
    N = a.shape[-1]
    b = a if axis == 2 else a.swapaxes(1, 2)
    sc = s / h
    m = np.floor(sc)
    u = (1 - (sc - m)).astype(np.float32)[..., None]
    base = np.arange(N)[None, None, :] - m.astype(np.intp)[:, :, None]
    ws = ((-2, -u * (u - 1) * (u - 2) / 6), (-1, (u + 1) * (u - 1) * (u - 2) / 2),
          (0, -(u + 1) * u * (u - 2) / 2), (1, (u + 1) * u * (u - 1) / 6))
    out = sum(w * np.take_along_axis(b, (base + off) % N, axis=2) for off, w in ws)
    return np.ascontiguousarray(out if axis == 2 else out.swapaxes(1, 2))


def shear(rng, flow, runs, N, h, axis):
    phi = rng.uniform(0, 2 * np.pi, runs)[:, None]
    coord = (np.arange(N) + 0.5) * h
    return flow["A"] * np.sin(2 * np.pi * flow["n"] * coord[None, :] / L + phi), axis


def coarse(a, f):
    R, N = a.shape[0], a.shape[1]
    return a if f == 1 else a.reshape(R, N // f, f, N // f, f).mean((2, 4))


def n_domains(s, amin=AMIN):
    """Periodic 4-connected domains of both signs with area >= amin (cells)."""
    tot = 0
    for val in (True, False):
        lab, nl = label(s == val)
        if nl == 0:
            continue
        par = list(range(nl + 1))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for p, q in list(zip(lab[0], lab[-1])) + list(zip(lab[:, 0], lab[:, -1])):
            if p and q:
                par[find(p)] = find(q)
        area = np.bincount(lab.ravel(), minlength=nl + 1)
        roots = {}
        for i in range(1, nl + 1):
            r = find(i)
            roots[r] = roots.get(r, 0) + area[i]
        tot += sum(v >= amin for v in roots.values())
    return tot


def wall_len(s):
    return ((s != np.roll(s, 1, 1)).sum((1, 2)) + (s != np.roll(s, 1, 2)).sum((1, 2))) / L


def run_arm(S, seed, Lbox=L, gamma=GAMMA, hold=HOLD, flow=None, D=D_EFF, h=1.0, dt=None, deltas=DELTAS):
    """Sweep lam -1 -> 1 then hold at 1. flow=None: no advection. Returns per-run arrays (sorted by Delta block)."""
    N = int(round(Lbox / h))
    dt = dt or min(0.05, 0.2 * h * h / D)
    assert D * dt / h ** 2 <= 0.25 + 1e-9, "Euler stability"
    rng = np.random.default_rng(seed)
    runs = S * len(deltas)
    g = np.repeat([K.g_for(d, Lbox ** 2) for d in deltas], S)
    a = np.zeros((runs, N, N), np.float32)
    gb = g.astype(np.float32).reshape(runs, 1, 1)
    noise, dtf, Dl = np.float32(np.sqrt(EPS / h ** 2 * dt)), np.float32(dt), np.float32(D / h ** 2)
    nsweep, nhold = int(round(2 / gamma / dt)), int(round(hold / dt))
    ke = int(round(flow["Tk"] / dt)) if flow else 0
    f = int(round(1 / h))
    probe = {nsweep - 1 + int(round(t / dt)): t for t in PROBES if t <= hold}
    cad = {nsweep - 1 + j * int(round(CADENCE / dt)) for j in range(int(hold / CADENCE) + 1)}
    out = dict(wall={}, ndom={}, t_single=np.full(runs, np.inf), deltas=np.repeat(deltas, S))
    kn = 0
    for i in range(nsweep + nhold):
        if flow and i % ke == 0:
            s, ax = shear(rng, flow, runs, N, h, 2 if kn % 2 == 0 else 1)
            a = kick(a, s, h, ax)
            kn += 1
        lam = np.float32(min(1.0, -1.0 + gamma * dt * i))
        a += (lam * a - a ** 3 + gb + Dl * K.lap2(a)) * dtf + noise * rng.standard_normal(a.shape, dtype=np.float32)
        if "early" not in out and lam >= 0.4:
            out["early"] = K.majority(coarse(a, f))[0]
        if i == nsweep - 1:
            out["end"] = K.majority(coarse(a, f))[0]
        if i in cad:
            c = coarse(a, f) > 0
            one = c.all((1, 2)) | (~c).all((1, 2))
            t_after = (i - (nsweep - 1)) * dt
            out["t_single"][one & np.isinf(out["t_single"])] = t_after
        if i in probe:
            c = coarse(a, f) > 0
            out["wall"][probe[i]] = wall_len(c)
            out["ndom"][probe[i]] = np.array([n_domains(x) for x in c])
    out["final"] = K.majority(coarse(a, f))[0]
    assert np.isfinite(a).all(), "NaN/inf in field"
    return out


def tracer_lagrangian(flow, seed):
    """Exact-map tracers: eddy diffusivity <dx^2+dy^2>/(4 t), 400 realisations x 50 tracers."""
    rng = np.random.default_rng(seed)
    R, P = 400, 50
    nk = 2 * int(round(100 / flow["Tk"] / 2))
    pos = rng.uniform(0, L, (R, P, 2))
    disp = np.zeros((R, P, 2))
    for j in range(nk):
        ax = j % 2                                    # kick along x uses y, along y uses x
        phi = rng.uniform(0, 2 * np.pi, (R, 1))
        d = flow["A"] * np.sin(2 * np.pi * flow["n"] * pos[:, :, 1 - ax] / L + phi)
        pos[:, :, ax] += d
        disp[:, :, ax] += d
    return (disp ** 2).sum(-1).mean() / (4 * nk * flow["Tk"])


def tracer_eulerian(flow, seed, h=0.5, runs=4):
    Lb = 4 * L if flow["n"] == 1 else L
    N, dt = int(round(Lb / h)), 0.05
    Kw = 2 * np.pi / Lb
    T = 1 / (Kw ** 2 * D_EFF)
    steps = int(round(T / dt))
    ke = int(round(flow["Tk"] / dt))
    rng = np.random.default_rng(seed)
    x = (np.arange(N) + 0.5) * h
    c = np.broadcast_to(np.cos(Kw * x)[None, None, :], (runs, N, N)).astype(np.float32).copy()
    ex = np.exp(1j * Kw * x)

    def amp(c):
        return np.abs((c.mean(1) * ex[None, :]).mean(-1)).mean() * 2
    a0, kn = amp(c), 0
    for i in range(steps):
        if i % ke == 0:
            s, ax = shear(rng, dict(flow, n=flow["n"]), runs, N, h, 2 if kn % 2 == 0 else 1)
            c = kick(c, s, h, ax)
            kn += 1
        c += np.float32(D_MOL * dt / h ** 2) * K.lap2(c)
    return -np.log(amp(c) / a0) / (Kw ** 2 * steps * dt)


def job(spec):
    t0 = time.time()
    kind, name = spec[:2]
    if kind == "arm":
        r = run_arm(**spec[2])
    elif kind == "lag":
        r = tracer_lagrangian(FLOWS[name], 1)
    else:
        r = tracer_eulerian(FLOWS[name], 2)
    return (kind, name), r, time.time() - t0


def sem(x):
    x = np.asarray(x, float)
    return x.std(ddof=1) / np.sqrt(len(x))


def fit(out, key):
    ks, ns = [], []
    for d in DELTAS:
        s = out[key][out["deltas"] == d]
        ks.append(s.sum()), ns.append(len(s))
    return K.fit_probit(DELTAS, ks, ns)


def wmean(o, t):
    return o["wall"][t].mean(), sem(o["wall"][t])


def verdict(n, ok, txt):
    print(f"P{n} {'HELD' if ok else 'FAILED'}  {txt}")


def main(smoke=False):
    global L
    t0 = time.time()
    if smoke:
        L = 16.0
        kw = dict(S=4, Lbox=16.0, gamma=0.5, hold=8.0)
    else:
        K.check()
        kw = dict(S=SEEDS)
    arms = [("arm", nm, dict(kw, seed=700 + j, flow=FLOWS[nm], D=D_MOL, h=0.5, dt=0.05)) for j, nm in enumerate("SMB")]
    arms += [("arm", "noflow8", dict(kw, seed=710, D=D_EFF, h=1.0, dt=0.025))]
    if not smoke:
        arms += [("arm", "noflow1", dict(kw, seed=711, D=D_MOL, h=1.0, dt=0.05))]
    jobs = arms + ([] if smoke else [(k, nm) for k in ("lag", "eul") for nm in "SMB"])
    print(f"\nL={L:g}, D_mol={D_MOL}, D_eff target={D_EFF}, xi={np.sqrt(D_EFF) * GAMMA ** -.25:.2f}, w={np.sqrt(D_EFF):.2f}; "
          f"{kw['S']} seeds x Deltas {DELTAS} per arm", flush=True)
    for k, v in FLOWS.items():
        ell = L / (2 * np.pi * v["n"])
        print(f"  flow {k}: n={v['n']} l={ell:.2f} (l/w={ell / np.sqrt(D_EFF):.2f}) Tk={v['Tk']} A={v['A']:.2f} "
              f"kappa={v['A'] * 2 * np.pi * v['n'] / L:.2f} U={v['A'] / v['Tk']:.1f} Da={ell / (v['A'] / v['Tk']):.2f}", flush=True)
    res = {}
    with ThreadPoolExecutor(2) as ex:
        for key, r, dt_ in ex.map(job, jobs):
            res[key] = r
            print(f"  job done: {key} [{dt_:.0f}s]", flush=True)
    if os.environ.get("ADVECT_PKL"):
        pickle.dump(res, open(os.environ["ADVECT_PKL"], "wb"))
    arm = {k[1]: v for k, v in res.items() if k[0] == "arm"}
    if smoke:
        for nm, o in arm.items():
            print(nm, {k: bool(np.isfinite(v).all()) for k, v in o.items() if hasattr(v, "shape")}, o["wall"][0][:2], o["ndom"][0][:2])
        print(f"smoke total {time.time() - t0:.0f}s")
        return

    print("\nD_eff (flow arms: design 8.0; eddy part design 7.0)")
    for nm in "SMB":
        print(f"  {nm}: Lagrangian eddy part {res[('lag', nm)]:.3f} (+D_mol = {res[('lag', nm)] + D_MOL:.3f}); "
              f"Eulerian total {res[('eul', nm)]:.3f}")

    print("\nOUTCOME  P(favoured hand) at Delta = " + ", ".join(map(str, DELTAS)) + "; fit P = Phi(a Delta)")
    a = {}
    for nm, o in arm.items():
        for key in ("early", "end", "final"):
            c, se = fit(o, key)
            a[(nm, key)] = (c, se)
            cells = "  ".join(f"D={d:g} P={K.pse(o[key][o['deltas'] == d])[0]:.3f}+-{K.pse(o[key][o['deltas'] == d])[1]:.3f}"
                              for d in DELTAS)
            print(f"  {nm:8s} {key:5s} {cells}  ->  a = {c:.3f} +- {se:.3f}", flush=True)

    print("\nWALLS  wall length / L (mean +- SE), number of domains (area >= 8), pooled over Deltas")
    for nm, o in arm.items():
        print(f"  {nm:8s} wall " + "  ".join(f"t{t}:{wmean(o, t)[0]:.2f}+-{wmean(o, t)[1]:.2f}" for t in PROBES))
        print(f"  {'':8s} ndom " + "  ".join(f"t{t}:{o['ndom'][t].mean():.2f}" for t in PROBES))
        ts = o["t_single"]
        print(f"  {'':8s} single-domain fraction by t_after = " +
              "  ".join(f"{t}:{(ts <= t).mean():.3f}" for t in (0, 10, 20, 40, 80)) +
              f"  median time {np.median(ts):.1f}" + ("" if np.isfinite(np.median(ts)) else " (>80, censored)"))

    print("\nPREDICTIONS")
    ok = all(abs(res[("lag", n)] / 7.0 - 1) < 0.05 for n in "SMB") and all(
        abs(res[("eul", n)] / 8.0 - 1) < tol for n, tol in (("S", .15), ("M", .15), ("B", .25)))
    verdict(1, ok, " ".join(f"{n}: lag {res[('lag', n)]:.2f} eul {res[('eul', n)]:.2f}" for n in "SMB"))
    c8, s8 = a[("noflow8", "end")]
    verdict(2, 0.85 <= c8 <= 1.15, f"check passed (asserted); a(noflow D=8) = {c8:.3f} +- {s8:.3f}")
    z = {n: (a[(n, "end")][0] - c8) / np.hypot(a[(n, "end")][1], s8) for n in "SMB"}
    verdict(3, all(abs(v) < 2 for v in z.values()), "a(flow) - a(noflow8) in combined SE: " + " ".join(f"{n}:{v:+.2f}" for n, v in z.items()))
    bad, zs = [], []
    for n in "SMB":
        c = a[(n, "end")][0]
        ps = [K.pse(arm[n]["end"][arm[n]["deltas"] == d]) for d in DELTAS]
        zs += [(n, d, (p - norm.cdf(c * d)) / s) for d, (p, s) in zip(DELTAS, ps)]
        if not (ps[0][0] < ps[1][0] < ps[2][0]):
            bad.append(f"{n} not monotone")
    verdict(4, not bad and all(abs(v) < 2.5 for _, _, v in zs), "P - Phi(a D) in SE: " + " ".join(f"{n}@{d:g}:{v:+.2f}" for n, d, v in zs) + " " + "; ".join(bad))
    verdict(5, all(0.8 <= a[(n, "end")][0] <= 1.2 for n in "SMB"), "a(end) " + " ".join(f"{n}:{a[(n, 'end')][0]:.3f}" for n in "SMB"))
    w0 = {n: wmean(arm[n], 0) for n in arm}
    w20 = {n: wmean(arm[n], 20) for n in arm}
    r0 = {n: w0[n][0] / w0["noflow8"][0] for n in "SMB"}
    r20 = {n: w20[n][0] / w20["noflow8"][0] for n in "SMB"}
    verdict(6, 0.75 <= r0["S"] <= 1.25 and r0["B"] >= 1.2, "wall(t0)/noflow8: " + " ".join(f"{n}:{r0[n]:.2f}" for n in "SMB") + f" (noflow8 {w0['noflow8'][0]:.2f})")
    verdict(7, 0.5 <= r20["S"] <= 2 and r20["M"] >= 0.5 and r20["B"] >= 1.3, "wall(t20)/noflow8: " + " ".join(f"{n}:{r20[n]:.2f}" for n in "SMB") + f" (noflow8 {w20['noflow8'][0]:.3f})")
    f80 = {n: (arm[n]["t_single"] <= 80).mean() for n in arm}
    verdict(8, abs(f80["S"] - f80["noflow8"]) < 0.15 and f80["noflow8"] - f80["B"] >= 0.15, "single by t80: " + " ".join(f"{n}:{v:.3f}" for n, v in f80.items()))
    verdict(9, w20["M"][0] >= w20["S"][0] - np.hypot(w20["M"][1], w20["S"][1]) and w20["B"][0] >= w20["M"][0] - np.hypot(w20["B"][1], w20["M"][1]),
            "wall(t20) S, M, B: " + " ".join(f"{n}:{w20[n][0]:.3f}+-{w20[n][1]:.3f}" for n in "SMB"))
    rb = w0["B"][0] / w0["noflow1"][0]
    verdict(10, 0.7 <= rb <= 1.4, f"wall(t0) B / noflow1 = {rb:.2f} (noflow1 {w0['noflow1'][0]:.2f})")
    el = time.time() - t0
    verdict(11, el < 1200, f"total {el:.0f}s")


if __name__ == "__main__":
    main("smoke" in sys.argv)
