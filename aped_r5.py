"""APED scorecard R5 (fast enough): published NCA-formation lower bounds run through aped.py's time
rescaling (sigma = s* a), scored against aped_r5_predictions.txt.

Predictions hashed before any rate from this unit was passed through the copied model:
  sha256 2fb0f2759ac8c188fc5319dd42eaecb515d3d7ad01210b6446886aa7a7b1c715, 2026-10-21 00:31 PDT (07:31 UTC)
Disclosure: (1) a few arithmetic-only python3 -I one-liners (the a_lo values and ratios of the P2-P9 data) were run
before the hash, to write correct hand values into the file; nothing went through the model. (2) The predictions were
NOT blind for P2-P7 (see the file's preamble and aped.py's POST-HOC paragraphs). (3) P9 was registered with the dose
factor 9.2e-4 = 2.3e-5/0.025; the exact 0.09 mol% x 25 mM = 2.25e-5 M gives 9.0e-4; both are printed, no verdict differs.
Sources (read by me): Leman, Orgel & Ghadiri, Science 306:283 (2004), doi:10.1126/science.1102722, Table 1 (p. 284),
notes 5, 8, 14; Danger et al., JACS 128:7412 (2006), doi:10.1021/ja061339+, Table 1 entries 2 and 3. Both PDFs are
copies sent by C. McKay (2026-10-20), not open-access downloads. Metal entries formed sulfide precipitates (Table 1 note).
Not found: any printed NCA formation rate constant or time course (searches in journals/new_analyses_aped_r5.md).
The model functions below are copied from aped.py (which runs a 86 s scan on import) and unchanged.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy python aped_r5.py > aped_r5.log
"""
import hashlib
import math

import numpy as np
from scipy.integrate import solve_ivp

PRED_SHA = "2fb0f2759ac8c188fc5319dd42eaecb515d3d7ad01210b6446886aa7a7b1c715"
YR, H, D = 3.156e7, 3600.0, 86400.0
SETTINGS = {"Warm little pond": (1e-2, 1e2 * YR), "Lake / lagoon": (1e-4, 1e4 * YR), "Open ocean": (1e-6, 1e8 * YR)}
S_STAR, X_STAR, SIG_BASE = 0.0063166, 0.3852, 0.0024332  # aped.log L89
NEED = {n: 10 / (tau * S_STAR) for n, (_, tau) in SETTINGS.items()}  # a needed, /s (aped.log L99-L101)
K_LOSS = math.log(2) / (10 * H)  # thiocarbamate loss, half-life 10 h (Leman note 8): upper bound on k_cyc
DOSE = 0.09e-2  # volcanic-gas COS mole fraction (Leman note 14) = dissolved COS as a fraction of the pure-COS saturation (note 5), Henry's law

# (id, route, t [s], f_lo, f_alt)  f_lo = residues-in-product basis (y2/2 + 2 y3/3), f_alt = molecules basis
def leman(y2, y3=0.0):
    return y2 / 2 + 2 * y3 / 3, y2 + 2 * y3


ROWS = [("L1", "Phe+COS, no additive, 25C, pH 8.7", 56 * H, *leman(.068)),
        ("L5", "Phe+COS+PbCl2, pH 9.1", 11 * H, *leman(.481, .016)),
        ("L6", "Phe+COS+PbCl2, pH 8.2", 3.5 * H, *leman(.342, .012)),
        ("L7", "Phe+COS+PbCl2, Pacific water", 3 * H, *leman(.203, .016)),
        ("L8", "Phe+COS+CdCl2", 16 * H, *leman(.400, .006)),
        ("L9", "Phe+COS+K3Fe(CN)6 100 mM", 20 * 60.0, *leman(.363)),
        ("D2", "cyanate CAA 20 mM alone, 80C", 28 * D, .017, None),
        ("D3", "cyanate CAA 20 mM + Ala 100 mM, 80C", 53 * D, (.127 + 2 * .019) * 20 / 120, None)]


def a_lo(f, t):
    return math.log(1 / (1 - f)) / t


def letters(lo, hi=math.inf):
    return {n: "Y" if lo >= NEED[n] else "N" if hi < NEED[n] else "?" for n in NEED}


def row_letter(L):
    v = set(L.values())
    return "?" if "?" in v else v.pop() if len(v) == 1 else "P"


# ---- copied from aped.py (unchanged) ----
N = ["l", "d", "ls", "ds", "ll", "ld", "dl", "dd"]
IX = {n: i for i, n in enumerate(N)}
W4 = np.array([2, 2, 4, 4.0])
SL, SD = [0, 2, 4, 5], [1, 3, 7, 6]
BASE = dict(ctot=2.0, kA=1.0, kH=1.0, kP=1.0, kD=1.0, kE=1.0, kR=1e-4, al=0.1, be=0.1, ga=0.1,
            Kact=1e-6, Kdim=1e-3, Kfx=1e-9, Vy=1.0, act="B")


def net(x, p):
    Vx, Vy, Ka, Kd, al, be, ga = x / p["Kfx"], p["Vy"], p["Kact"], p["Kdim"], p["al"], p["be"], p["ga"]
    kA, kH, kP, kD, kE, kR = (p[k] for k in ("kA", "kH", "kP", "kD", "kE", "kR"))
    rx = []
    for h, o in ("ld", "dl"):
        hh, ho, oh = h + h, h + o, o + h
        fA, rA = (kA, kA / (Ka * Vx)) if p["act"] == "A" else (kA * x, kA * p["Kfx"] / Ka)
        rx += [([h], [h + "s"], fA, rA),
               ([h + "s"], [h], kH, kH * Ka * Vy),
               ([h + "s", h], [hh], kP, kP * Ka * Vy / Kd),
               ([h + "s", o], [ho], al * kP, al * kP * Ka * Vy / (ga * Kd)),
               ([hh], [h, h], kD, kD * Kd),
               ([ho], [h, o], be * kD, be * ga * kD * Kd),
               ([oh], [hh], kE, ga * kE)]
    rx.append((["l"], ["d"], kR, kR))
    Rm, Pm = np.zeros((len(rx), 8)), np.zeros((len(rx), 8))
    for j, (a, b, _, _) in enumerate(rx):
        for s in a:
            Rm[j, IX[s]] += 1
        for s in b:
            Pm[j, IX[s]] += 1
    return (Pm - Rm).T, Rm, Pm, np.array([r[2] for r in rx]), np.array([r[3] for r in rx])


def prod(z, m):
    return np.prod(z ** m, axis=1)


def dprod(z, m):
    Dm = np.empty_like(m)
    for i in range(8):
        t = z ** m
        t[:, i] = m[:, i] * z[i] ** np.maximum(m[:, i] - 1, 0)
        Dm[:, i] = t.prod(1)
    return Dm


def rhs(z, n):
    S, Rm, Pm, kf, kr = n
    return S @ (kf * prod(z, Rm) - kr * prod(z, Pm))


def jac(z, n):
    S, Rm, Pm, kf, kr = n
    return S @ (kf[:, None] * dprod(z, Rm) - kr[:, None] * dprod(z, Pm))


def anti(J):
    return J[np.ix_(SL, SL)] - J[np.ix_(SL, SD)]


def sym(J):
    return J[np.ix_(SL, SL)] + J[np.ix_(SL, SD)]


def full(y):
    z = np.empty(8)
    z[SL], z[SD] = y, y
    return z


def newton(y, n, p):
    for _ in range(60):
        z = full(y)
        F, J = rhs(z, n)[SL], sym(jac(z, n))
        F[0], J[0] = W4 @ y - p["ctot"], W4
        dy = np.linalg.solve(J, -F)
        t = 1.0
        while np.any(y + t * dy <= 0) and t > 1e-8:
            t /= 2
        y = y + t * dy
        if np.max(np.abs(t * dy) / y) < 1e-10 and t == 1.0:
            return y
    return None


def racemic(x, p, y0=None):
    n = net(x, p)
    if y0 is not None:
        y = newton(np.array(y0, float), n, p)
        if y is not None:
            return y
    y = np.array([p["ctot"] / 2, 1e-30, 1e-30, 1e-30])
    for T in (1e4, 1e7, 1e10):
        sol = solve_ivp(lambda t, v: rhs(full(v), n)[SL], [0, T], y, method="Radau",
                        jac=lambda t, v: sym(jac(full(v), n)), rtol=1e-9, atol=1e-22)
        y = np.maximum(sol.y[:, -1], 1e-30)
        yn = newton(y, n, p)
        if yn is not None:
            return yn
    raise RuntimeError(f"no symmetric fixed point at x={x}")


def lam(x, p, y):
    ev = np.linalg.eigvals(anti(jac(full(y), net(x, p))))
    return ev[np.argmax(ev.real)]
# ---- end copy ----


def scan_ratio(r, kR):
    """Leading difference-mode eigenvalue at x*, with activation a = r*kH slowed against the other base rates."""
    p = dict(BASE, kA=r * BASE["kH"] / X_STAR, kR=kR)
    try:
        y = racemic(X_STAR, p)
    except RuntimeError:
        return None
    L = lam(X_STAR, p, y)
    resid = np.abs(rhs(full(y), net(X_STAR, p))).max() / np.abs(net(X_STAR, p)[3]).max()
    return L, resid


def main():
    with open("aped_r5_predictions.txt", "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == PRED_SHA, "predictions file changed"
    assert abs(a_lo(0.034, 56 * H) - 1.716e-7) < 1e-10  # hand check: ln(1/0.966) / 201600 s
    y0 = racemic(X_STAR, BASE)
    assert abs(lam(X_STAR, BASE, y0).real - SIG_BASE) < 2e-7, "copied model does not reproduce aped.log L89"
    print(f"a needed (aped.log L99-L101): " + ", ".join(f"{n} {v:.3e}" for n, v in NEED.items()) + " /s;  s* = %.7g" % S_STAR)
    print(f"copied model at x* = {X_STAR} M, kA = 1, kR = 1e-4: sigma = {lam(X_STAR, BASE, y0).real:.5g} /s (aped.log L89: {SIG_BASE})\n")

    res = {}
    print("Lower bounds a_lo = ln(1/(1-f))/t; ratio = a_lo / a_needed; letters (pond, lake, ocean)")
    for id_, route, t, f, fa in ROWS:
        lo = a_lo(f, t)
        hi = K_LOSS if id_ == "L1" else math.inf  # thiocarbamate-loss bound holds for the plain route only
        L = letters(lo, hi)
        res[id_] = (lo, hi, L)
        rat = "  ".join(f"{lo / NEED[n]:9.3g}" for n in NEED)
        sg = lo * S_STAR
        print(f"{id_}  {route:37} f {f:.4f}  t {t / H:7.1f} h  a_lo {lo:.3e} /s  sigma {sg:.2e} /s  ratio {rat}  "
              f"letters {''.join(L.values())}")
    lo1, hi1, L1 = res["L1"]
    alt = a_lo(ROWS[0][4], ROWS[0][2])
    print(f"\nL1 molecule-basis sensitivity: f {ROWS[0][4]:.3f}  a_lo {alt:.3e} /s  letters {''.join(letters(alt, hi1).values())}")
    print(f"L1 upper bound a_hi = ln2/(10 h) = {hi1:.3e} /s = {hi1 / NEED['Warm little pond']:.1f} x pond need")
    dose = {}
    for d in (9.2e-4, DOSE):
        Ld = letters(lo1 * d, hi1 * d)
        dose[d] = {n: lo1 * d / NEED[n] for n in NEED}
        print(f"L1 dose-scaled by {d:.1e} (a proportional to dissolved COS, unmeasured): a_lo {lo1 * d:.2e} /s  "
              f"ratio " + "  ".join(f"{v:.3g}" for v in dose[d].values()) + f"  letters {''.join(Ld.values())}")

    ratios = (0.385, 3.85e-2, 3.85e-3, 3.85e-4, 3.85e-5, 3.85e-6)
    print("\nModel check at x* (copied aped.py functions, other base rates kept): a/kH, sigma(kR=0), sigma/a, "
          "sigma(kR=1e-4), residual")
    s0, s1 = [], []
    for r in ratios:
        a0, a1 = scan_ratio(r, 0.0), scan_ratio(r, 1e-4)
        s0.append(a0 and a0[0])
        s1.append(a1 and a1[0])
        f0 = f"{a0[0].real:+.4e} (Im {a0[0].imag:.1e}) sigma/a {a0[0].real / r:.4e}" if a0 else "no fixed point"
        f1 = f"{a1[0].real:+.4e}" if a1 else "no fixed point"
        rs = f"{max(a0[1], a1[1]):.1e}" if a0 and a1 else "n/a"
        print(f"  a/kH {r:.3g}: kR=0 {f0};  kR=1e-4 {f1};  resid {rs}")

    conv = all(s is not None for s in s0 + s1)
    p10 = conv and all(S_STAR / 3 <= s.real / r <= 3 * S_STAR for s, r in zip(s0, ratios))
    p11 = conv and all(s.real < 0 for s, r in zip(s1, ratios) if r <= 3.85e-3) and s1[0].real > 0
    rowl = row_letter(L1)
    res_ = [
        ("P1", True, "no direct NCA-formation rate printed in Leman 2004 or Danger 2006 (end points only); OA searches found none"),
        ("P2", 1.63e-7 <= lo1 <= 1.80e-7 and all(abs(lo1 / NEED[n] / v - 1) < 0.02 for n, v in zip(NEED, (0.34, 34, 3.4e5))),
         f"a_lo {lo1:.3e}; ratios " + ", ".join(f"{lo1 / NEED[n]:.3g}" for n in NEED)),
        ("P3", 1.90e-5 <= hi1 <= 1.95e-5 and abs(hi1 / NEED["Warm little pond"] - 38) < 1, f"a_hi {hi1:.3e}, {hi1 / NEED['Warm little pond']:.1f}x pond need"),
        ("P4", ''.join(L1.values()) == "?YY" and ''.join(letters(alt, hi1).values()) == "?YY", f"letters {''.join(L1.values())}; molecule basis {''.join(letters(alt, hi1).values())}"),
        ("P5", 1.60e-4 <= res["L9"][0] <= 1.75e-4 and ''.join(res["L9"][2].values()) == "YYY", f"L9 a_lo {res['L9'][0]:.3e}, letters {''.join(res['L9'][2].values())}"),
        ("P6", all(''.join(res[k][2].values()) == "YYY" for k in ("L5", "L6", "L7", "L8")) and 3.7e-6 <= res["L8"][0] <= 4.2e-6,
         "L5-L8 pond ratios " + ", ".join(f"{res[k][0] / NEED['Warm little pond']:.1f}" for k in ("L5", "L6", "L7", "L8")) + f"; L8 a_lo {res['L8'][0]:.3e}"),
        ("P7", 5.5e-9 <= res["D3"][0] <= 6.7e-9 and ''.join(res["D3"][2].values()) == "?YY",
         f"D3 a_lo {res['D3'][0]:.3e}; ratios " + ", ".join(f"{res['D3'][0] / NEED[n]:.3g}" for n in NEED) + f"; letters {''.join(res['D3'][2].values())}"),
        ("P8", rowl == "?", f"row letter {rowl} (pond/lake/ocean {'/'.join(L1.values())})"),
        ("P9", all(dose[d]["Lake / lagoon"] < 1 and dose[d]["Open ocean"] > 1 and dose[d]["Warm little pond"] < 1e-3 for d in dose),
         "dose-scaled ratios (lake, ocean) " + "; ".join(f"{d:.1e}: {dose[d]['Lake / lagoon']:.3g}, {dose[d]['Open ocean']:.3g}" for d in dose)),
        ("P10", p10, "sigma/a (kR=0) " + ", ".join(f"{s.real / r:.3e}" if s is not None else "n/a" for s, r in zip(s0, ratios)) + f" vs [{S_STAR / 3:.3e}, {3 * S_STAR:.3e}]"),
        ("P11", p11, "sigma(kR=1e-4) " + ", ".join(f"{s.real:+.2e}" if s is not None else "n/a" for s in s1)),
    ]
    print()
    for name, ok, why in res_:
        print(f"{name} {'HELD  ' if ok else 'FAILED'} {why}")
    print(f"{sum(ok for _, ok, _ in res_)} of {len(res_)} held.  R5 = {rowl}"
          + ("" if p10 else "  (P10 failed: letters hold under the aped.py time rescaling only)") + " (scorecard.py not changed).")


if __name__ == "__main__":
    main()
