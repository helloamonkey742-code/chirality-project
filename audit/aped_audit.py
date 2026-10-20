"""Independent re-implementation audit of the Plasson, Bersini & Commeyras 2004 APED model.

PREDICTIONS_SHA256 = 350d96c44d8b55e8cb2d2a23b290a8d9ad5ec1279035fe0cffa16c1c51a8c41b
  (shasum -a 256 of the text between the markers below, i.e. the 44-line file P.txt)
PREDICTIONS_WRITTEN = Mon Oct 19 23:38:52 PDT 2026   (output of `date`, before any run)
Session: audit-aped (blinded to aped.py / aped.log / scorecard / README / draft until after the run).

Run: /opt/anaconda3/bin/python3 audit/aped_audit.py | tee audit/aped_audit.log

<<<PREDICTIONS
SOURCE: Plasson, Bersini & Commeyras, PNAS 101:16733 (2004), open-access full text
  PMC534711 (HTML + Fig. 1 + Methods matrix images M1/M2). No aped.py / aped.log /
  scorecard / README / draft consulted. Saved under papers/audit-aped/.

READINGS (R)
R1 Mass action. 16 elementary steps of Fig. 1, rate = k * product of reactant conc.
   Activation a*[L], a*[D] (FIRST ORDER in the deactivated monomer; the paper's
   Methods matrix M1 has d[L*]/d[L] = a - x_l and d[L]/d[L] = -a - x_l - alpha*x_d).
   Deactivation b*[L*]. Polymerisation p*[X*][Y] (heterochiral alpha*p, activated
   residue written on the left). Hydrolysis h*[XX], beta*h*[XY]. Epimerisation of the
   left (N-terminal) residue e*[XY->XX], gamma*e*[XX->XY].
   c = l + d + l* + d* + 2(ll + dd + ld + dl) conserved.
   The paper has no explicit fuel species: drive == activation rate a (and a/b).
R2 Parameter sets. S0 = Fig. 3 (b=0, beta=gamma=0, a=h=e=p=1). S1 = Methods/Fig. 4
   (b=0, a=h=e=p=1, alpha=0.1, beta=gamma=0.01; c = 2 M per Methods, c = 1 M per the
   Fig. 4 caption: the paper disagrees with itself). S2 = Fig. 5 (a=1e-8, b=5e-4,
   p=2e-2, alpha=0.35, h=e=1e-7, beta=0.2, gamma=0.3).

PREDICTIONS (P)
P1 For b = 0 the racemic fixed point exists iff c > 2a/(p(1+alpha)), whatever beta,
   gamma, h, e are (algebra from Fig. 1). alpha = 0.1 -> c_exist = 1.8182 M.
   So c = 1 M (Fig. 4 caption) is a dead system for S1; Fig. 4 must have used 2 M.
P2 S0, alpha = 0.1: racemic state is unstable (one real positive eigenvalue in the
   antisymmetric block) from its birth, for all c up to at least 5 M:
   c_min = 1.8182 M. S1 (beta = gamma = 0.01): same c_min to 1e-6 relative.
P3 Drive window, S0, c = 2, alpha = 0.1, b = 0: unstable for all 1e-6 <= a < 1.1
   (= p c (1+alpha)/2); no racemic state (dead) for a >= 1.1. Upper edge 1.1 exactly.
P4 Drive window, S2 at c = 1 M, a scanned with b = 5e-4 fixed: a = 1e-8 is inside the
   unstable window. Edge values: no a-priori number; they will simply be reported.
P5 S2 c_min lies between 20 mM and 100 mM (Fig. 5: 20 mM flat, 100 mM still rising).
   At c = 1 M, integration from ee0 = 0.01 ends at |ee| in [0.65, 0.75]
   (text: about 70%; Fig. 5 plateau about 0.695).
P6 Sign test, S0, c = 2, beta = gamma = 0: alpha = 0.1 unstable real; alpha = 10
   stable (damped oscillation); alpha = 100 unstable OSCILLATORY (complex pair,
   Re > 0). alpha_c(c = 2) about 17 (Fig. 3F, by eye, +-2). So racemic instability
   alone does not need alpha < 1, but a stable asymmetric end state does.
P7 Sign test, full model, c = 2, b = 0, a=h=e=p=1, (alpha,beta,gamma) in {0.1,1,10}^3:
   unstable-real iff alpha < 1 and at least one of beta, gamma < 1.
   Specifically (0.1,0.1,0.1) U; (0.1,0.1,10) U; (0.1,10,0.1) U; (0.1,1,1) S;
   (0.1,10,10) S; no alpha = 10 entry is unstable-real.
P8 Time integration from ee0 = 0.01 agrees with the Jacobian: |ee| grows iff the
   antisymmetric block has Re(lambda) > 0.
P9 Self-check: my Jacobian equals the paper's matrix M1 entry for entry at random
   states when b = beta = gamma = 0 (max abs diff < 1e-12).
PREDICTIONS>>>
"""
import hashlib, re, sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

# ---- integrity: predictions block above must hash to the value recorded at the top
src = open(__file__).read()
blk = src.split("<<<PREDICTIONS\n", 1)[1].split("PREDICTIONS>>>", 1)[0]
h = hashlib.sha256(blk.encode()).hexdigest()
rec = re.search(r"PREDICTIONS_SHA256 = (\w+)", src).group(1)
print("predictions sha256 recomputed:", h, "| recorded:", rec, "| MATCH" if h == rec else "| MISMATCH")
assert h == rec

# ---- model: state order follows the paper's matrix M1: l d l* d* ll dd dl ld
IDX = dict(l=0, d=1, ls=2, ds=3, ll=4, dd=5, dl=6, ld=7)
W = np.array([1, 1, 1, 1, 2, 2, 2, 2.])
N = np.zeros((8, 16))
def _st(col, **k):
    for n, v in k.items(): N[IDX[n], col] += v
_st(0, l=-1, ls=1); _st(1, d=-1, ds=1); _st(2, ls=-1, l=1); _st(3, ds=-1, d=1)
_st(4, ls=-1, l=-1, ll=1); _st(5, ds=-1, l=-1, dl=1); _st(6, ls=-1, d=-1, ld=1); _st(7, ds=-1, d=-1, dd=1)
_st(8, ll=-1, l=2); _st(9, dl=-1, l=1, d=1); _st(10, ld=-1, l=1, d=1); _st(11, dd=-1, d=2)
_st(12, ld=-1, dd=1); _st(13, dd=-1, ld=1); _st(14, ll=-1, dl=1); _st(15, dl=-1, ll=1)
assert np.allclose(W @ N, 0), "residue conservation violated by stoichiometry"

def rates(x, P):
    l, d, ls, ds, ll, dd, dl, ld = x
    a, b, p, al, h, be, e, ga = (P[k] for k in "a b p al h be e ga".split())
    return np.array([a*l, a*d, b*ls, b*ds, p*ls*l, al*p*ds*l, al*p*ls*d, p*ds*d,
                     h*ll, be*h*dl, be*h*ld, h*dd, e*ld, ga*e*dd, ga*e*ll, e*dl])

def rhs(x, P): return N @ rates(x, P)

def jac(x, P):
    l, d, ls, ds, ll, dd, dl, ld = x
    a, b, p, al, h, be, e, ga = (P[k] for k in "a b p al h be e ga".split())
    R = np.zeros((16, 8)); I = IDX
    R[0, I['l']] = a; R[1, I['d']] = a; R[2, I['ls']] = b; R[3, I['ds']] = b
    R[4, I['ls']] = p*l; R[4, I['l']] = p*ls
    R[5, I['ds']] = al*p*l; R[5, I['l']] = al*p*ds
    R[6, I['ls']] = al*p*d; R[6, I['d']] = al*p*ls
    R[7, I['ds']] = p*d; R[7, I['d']] = p*ds
    R[8, I['ll']] = h; R[9, I['dl']] = be*h; R[10, I['ld']] = be*h; R[11, I['dd']] = h
    R[12, I['ld']] = e; R[13, I['dd']] = ga*e; R[14, I['ll']] = ga*e; R[15, I['dl']] = e
    return N @ R

def M1(x, P):  # paper's Methods matrix, transcribed from image M1 (valid for b=beta=gamma=0)
    l, d, ls, ds, ll, dd, dl, ld = x
    a, p, al, h, e = P["a"], P["p"], P["al"], P["h"], P["e"]
    xl, xd, yl, yd = p*ls, p*ds, p*l, p*d
    return np.array([
        [-a-xl-al*xd, 0, -yl, -al*yl, 2*h, 0, 0, 0],
        [0, -a-xd-al*xl, -al*yd, -yd, 0, 2*h, 0, 0],
        [a-xl, -al*xl, -yl-al*yd, 0, 0, 0, 0, 0],
        [-al*xd, a-xd, 0, -yd-al*yl, 0, 0, 0, 0],
        [xl, 0, yl, 0, -h, 0, e, 0],
        [0, xd, 0, yd, 0, -h, 0, e],
        [al*xd, 0, 0, al*yl, 0, 0, -e, 0],
        [0, al*xl, al*yd, 0, 0, 0, 0, -e]])

def P_(**k):
    P = dict(a=1., b=0., p=1., al=.1, h=1., be=0., e=1., ga=0.); P.update(k); return P

S0 = lambda **k: P_(**k)
S1 = lambda **k: P_(**{**dict(al=.1, be=.01, ga=.01), **k})
S2 = lambda **k: P_(**{**dict(a=1e-8, b=5e-4, p=2e-2, al=.35, h=1e-7, be=.2, e=1e-7, ga=.3), **k})

# ---- self-checks: M1 equality, finite-difference Jacobian, mirror symmetry
rng = np.random.default_rng(0)
worst = 0
for _ in range(200):
    x = rng.random(8); P = P_(a=rng.random()*3, p=rng.random()*3, al=rng.random()*5, h=rng.random()*3, e=rng.random()*3)
    worst = max(worst, np.abs(jac(x, P) - M1(x, P)).max())
print("P9 max|J - M1| over 200 random states (b=beta=gamma=0):", worst); assert worst < 1e-12
P = S2(); x = rng.random(8); eps = 1e-6
fd = np.column_stack([(rhs(x+eps*np.eye(8)[i], P)-rhs(x-eps*np.eye(8)[i], P))/(2*eps) for i in range(8)])
print("max|J - finite diff| (S2, random state):", np.abs(fd-jac(x, P)).max()/np.abs(fd).max(), "(relative)")
MIR = np.zeros((8, 8))
for u, v in [('l', 'd'), ('ls', 'ds'), ('ll', 'dd'), ('dl', 'ld')]:
    MIR[IDX[u], IDX[v]] = MIR[IDX[v], IDX[u]] = 1
x = rng.random(8); P = S2()
print("mirror symmetry |rhs(Px)-P rhs(x)|:", np.abs(rhs(MIR@x, P)-MIR@rhs(x, P)).max())
# antisymmetric / symmetric bases
BA = np.zeros((8, 4)); BS = np.zeros((8, 4))
for k, (u, v) in enumerate([('l', 'd'), ('ls', 'ds'), ('ll', 'dd'), ('dl', 'ld')]):
    BA[IDX[u], k], BA[IDX[v], k] = 1, -1; BS[IDX[u], k], BS[IDX[v], k] = 1, 1

# ---- racemic fixed point: unique root in m = [L] = [D] (total residue is monotone in m)
def racemic(c, P):
    a, b, p, al, h, be, e, ga = (P[k] for k in "a b p al h be e ga".split())
    Dn = h*(e + be*h + ga*be*e)
    def parts(m):
        s = a*m/(b + (1+al)*p*m); q = p*s*m
        return s, q*(e + be*h + al*e)/Dn, q*(al*h + al*ga*e + ga*e)/Dn   # s, u=[LL]=[DD], v=[DL]=[LD]
    f = lambda lm: (lambda m: 2*m + 2*parts(m)[0] + 4*sum(parts(m)[1:]) - c)(np.exp(lm))
    lo, hi = -700., np.log(c/2)
    if f(lo) >= 0 or f(hi) < 0: return None
    m = np.exp(brentq(f, lo, hi, xtol=1e-14, rtol=1e-14, maxiter=500)); s, u, v = parts(m)
    x = np.array([m, m, s, s, u, u, v, v])
    return dict(x=x, res=np.abs(rhs(x, P)).max()/max(np.abs(N*rates(x, P)).max(), 1e-300), cerr=abs(W@x-c)/c)

def classify(c, P):
    r = racemic(c, P)
    if r is None: return "none", None, None
    J = jac(r['x'], P)
    la = np.linalg.eigvals(BA.T @ J @ BA / 2); ls_ = np.linalg.eigvals(BS.T @ J @ BS / 2)
    tol = 1e-10*np.abs(J).max()
    k = la.real.argmax(); lam = la[k]
    if lam.real > tol: cls = "U" if abs(lam.imag) < 1e-9*abs(lam.real) + 1e-300 else "O"   # real (symmetry-breaking) / oscillatory
    else: cls = "S"
    return cls, lam, ls_.real.max()

def bis(fcls, lo, hi, n=60):  # log-bisect between grid points with different class
    cl, ch = fcls(lo), fcls(hi)
    for _ in range(n):
        mid = np.sqrt(lo*hi); cm = fcls(mid)
        if cm == cl: lo = mid
        else: hi = mid
    return np.sqrt(lo*hi), cl, ch

def scan(var, grid, base, c, name):
    """Print class along a log grid of `var` (a parameter name or 'c'), and refined edges."""
    def cls(v):
        P = dict(base); cc = c
        if var == 'c': cc = v
        else: P[var] = v
        return classify(cc, P)[0]
    cl = [cls(v) for v in grid]
    runs = []
    for i, k in enumerate(cl):
        if not runs or runs[-1][0] != k: runs.append([k, grid[i], grid[i]])
        else: runs[-1][2] = grid[i]
    print(f"  {name}: classes along grid [{grid[0]:.3g}..{grid[-1]:.3g}] ({len(grid)} pts): " +
          "; ".join(f"{k} on [{a_:.4g},{b_:.4g}]" for k, a_, b_ in runs))
    for i in range(len(cl)-1):
        if cl[i] != cl[i+1]:
            e_, a_, b_ = bis(cls, grid[i], grid[i+1])
            print(f"     edge {a_}->{b_} at {var} = {e_:.6g}")
    return cl

def integ(c, P, ee0=0.01, tend=None, n=3000, x0=None):
    x0 = np.array([c/2*(1+ee0), c/2*(1-ee0), 0, 0, 0, 0, 0, 0]) if x0 is None else x0
    y0 = x0/c
    f = lambda t, y: rhs(y*c, P)/c
    ts = np.linspace(0, tend, n)
    sol = solve_ivp(f, (0, tend), y0, method='Radau', t_eval=ts, rtol=1e-10, atol=1e-16, jac=lambda t, y: jac(y*c, P))
    Y = sol.y*c; ee = (Y[0]+Y[2]+2*Y[4]-Y[1]-Y[3]-2*Y[5])/c
    return sol.t, ee, Y

def horizon(c, P, mult=25, cap=1e14):
    r = racemic(c, P)
    if r is None: return 200/max(P['a'], P['h'], P['e'], 1e-30)
    la = np.linalg.eigvals(BA.T @ jac(r['x'], P) @ BA / 2)
    return min(mult/max(abs(la.real.max()), 1e-30), cap)

def run_ee(label, c, P, ee0=0.01, mult=25):
    tend = horizon(c, P, mult); t, ee, Y = integ(c, P, ee0, tend)
    cross = int(np.sum(np.diff(np.sign(ee[np.abs(ee) > 1e-12])) != 0))
    print(f"    [int] {label}: c={c:g} tend={tend:.3g} ee: start {ee[0]:+.4f} max|ee| {np.abs(ee).max():.4f} end {ee[-1]:+.5f} sign-changes {cross}; residue drift {abs(W@Y[:, -1]-c)/c:.1e}")
    return ee

if __name__ == "__main__":
    print("\n=== A. VALIDATION against the paper's own figures ===")
    print("Fig.3F existence boundary c=2a/(p(1+alpha)) [S0]; computed edge of racemic existence vs formula")
    for al in (0.01, 0.1, 0.5, 1.0):
        cl = scan('c', np.geomspace(0.3, 5, 120), S0(al=al), None, f"S0 alpha={al}")
        print(f"     formula 2/(1+alpha) = {2/(1+al):.6f}")
    print("Fig.3 marked points (S0, ee0=0.1 for A; 0.01 otherwise):")
    for tag, al, c in [("A dead", .1, 1), ("B asym", .1, 3), ("C sym", 1.1, 3.2), ("D sym+osc", 50, 2), ("E unstable-osc", 100, 2)]:
        cls, lam, sl = classify(c, S0(al=al)); print(f"  {tag}: alpha={al} c={c}: class={cls} lam_a={lam} lam_sym_max={sl}")
        run_ee(tag, c, S0(al=al), ee0=0.1 if tag[0] == 'A' else 0.01)
    print("alpha_c(c) (Hopf edge, S0) vs Fig.3F read by eye (c=1:~12, 2:~17, 3:~21, 5:~27):")
    for c in (1.0, 2.0, 3.0, 5.0):
        scan('al', np.geomspace(2, 200, 100), S0(), c, f"c={c}")
    print("Fig.5 check (S2, ee0=0.01) final ee vs figure (1M:0.695, 400mM:~0.655, 200mM:~0.585 @2e10 s, 100mM:~0.09 @2e10 s, 20mM flat):")
    for c in (1.0, 0.4, 0.2, 0.1, 0.02):
        P = S2(); t, ee, Y = integ(c, P, 0.01, 2e10, n=400)
        i = np.searchsorted(t, 2e10) - 1
        t2, ee2, _ = integ(c, P, 0.01, 1e12, n=400)
        print(f"   c={c*1000:6.0f} mM: ee(2e10 s)={ee[-1]:+.4f}  ee(1e12 s)={ee2[-1]:+.4f}")

    print("\n=== B. c_min: smallest c at which the racemic state is unstable ===")
    for name, P in [("S0 a=1 alpha=.1 b=0 beta=gamma=0", S0()), ("S1 alpha=.1 beta=gamma=.01 b=0", S1()), ("S2 (Fig.5 set)", S2())]:
        print(f" {name}")
        grid = np.geomspace(1e-3, 10, 160) if P['b'] == 0 else np.geomspace(1e-4, 10, 220)
        scan('c', grid, P, None, "c scan")
        if P['b'] == 0:
            c0 = 2*P['a']/(P['p']*(1+P['al']))
            cls, lam, sl = classify(c0*(1+1e-9), P)
            print(f"     at c0(1+1e-9), c0={c0:.6f}: class={cls} lam_a={lam} lam_sym_max={sl}")
    print(" S1 with c = 1 M (Fig.4 caption) vs 2 M (Methods):", classify(1.0, S1())[0], "vs", classify(2.0, S1())[0])

    print("\n=== C. drive window (activation a) ===")
    print(" S0 c=2 alpha=.1 b=0; formula upper edge p c (1+alpha)/2 = 1.1")
    scan('a', np.geomspace(1e-6, 10, 200), S0(), 2.0, "a scan")
    print(" S1 c=2 b=0")
    scan('a', np.geomspace(1e-6, 10, 200), S1(), 2.0, "a scan")
    for c in (1.0, 0.2, 0.1):
        print(f" S2 c={c} M, b=5e-4 fixed, a scanned")
        scan('a', np.geomspace(1e-16, 1e2, 300), S2(), c, "a scan")
    print(" S2 c=1 M: drive as a/b (b scanned at a=1e-8 fixed) -- same window in the other variable")
    scan('b', np.geomspace(1e-12, 1e2, 300), S2(), 1.0, "b scan")
    print(" integration checks across the S2 window (c=1 M): a = 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2")
    for a in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2):
        cls, lam, _ = classify(1.0, S2(a=a)); print(f"   a={a:g} class={cls} lam_a_max={lam.real if lam is not None else None}")
        run_ee(f"S2 a={a:g}", 1.0, S2(a=a), mult=12)

    print("\n=== D. sign tests ===")
    print(" S0 c=2 beta=gamma=0, alpha scan:")
    for al in (0.01, 0.1, 0.5, 0.9, 1.0, 1.1, 2, 5, 10, 17, 25, 50, 100):
        cls, lam, sl = classify(2.0, S0(al=al)); print(f"   alpha={al:<5g} class={cls} lam_a={lam:.5g}  sym-block max Re={sl:.3g}")
    for al in (0.1, 10, 100):
        run_ee(f"S0 alpha={al}", 2.0, S0(al=al), mult=40)
    print(" full model c=2, b=0, a=h=e=p=1: class / lam_a over (alpha,beta,gamma) in {0.1,1,10}^3, then integration")
    for al in (0.1, 1, 10):
        for be in (0.1, 1, 10):
            for ga in (0.1, 1, 10):
                P = S0(al=al, be=be, ga=ga); cls, lam, sl = classify(2.0, P)
                print(f"   (alpha,beta,gamma)=({al:g},{be:g},{ga:g}) class={cls} lam_a_max={lam if lam is None else format(lam, '.4g')} sym-block max Re={sl if sl is None else format(sl, '.2g')}")
                if cls != "none": run_ee("   ", 2.0, P, mult=30)
    print(" S2 sign flips (c=1 M): each of alpha,beta,gamma moved to 1/x or 1x or 10x of its paper value")
    for k in ("al", "be", "ga"):
        for f in (0.1, 1, 10, 100):
            P = S2(); P[k] = P[k]*f; cls, lam, _ = classify(1.0, P)
            print(f"   {k} = {P[k]:.4g}: class={cls} lam_a_max={lam if lam is None else format(lam, '.4g')}")
    print(" S2 with beta=gamma=1 and alpha=0.35 (polymerisation selectivity only): ", classify(1.0, S2(be=1., ga=1.))[:2])
    print(" S2 with alpha=1, beta=.2, gamma=.3 (hydrolysis/epimerisation selectivity only): ", classify(1.0, S2(al=1.))[:2])
    print(" S2 alpha=1, beta=gamma=1: ", classify(1.0, S2(al=1., be=1., ga=1.))[:2])
    print("\n=== E. follow-ups added after run 1 (tolerance artefact, finer sign boundary, Hopf region) ===")
    # E1: classify() uses tol = 1e-10*max|J|; for S2 (rates ~1e-8) that hides tiny positive eigenvalues.
    print(" E1 S2 c=1 M b=5e-4: lam_a_max(a) and lam_a/a, sign only (no tolerance); max|J| shown")
    for a_ in (1e-16, 1e-14, 1e-12, 1e-11, 1e-10, 1e-8, 1e-7, 5e-7, 9e-7, 9.9e-7, 1.0e-6, 1.1e-6):
        P = S2(a=a_); r = racemic(1.0, P); J = jac(r['x'], P)
        la = np.linalg.eigvals(BA.T @ J @ BA / 2); k = la.real.argmax()
        print(f"   a={a_:.3g}: lam_a_max={la[k]:.4g}  lam/a={la[k].real/a_:.4g}  max|J|={np.abs(J).max():.3g}  tol_used={1e-10*np.abs(J).max():.3g}")
    def sign_upper(a_):
        P = S2(a=a_); r = racemic(1.0, P); la = np.linalg.eigvals(BA.T @ jac(r['x'], P) @ BA / 2); return "U" if la.real.max() > 0 else "S"
    e_, c1_, c2_ = bis(sign_upper, 1e-8, 1e-5)
    print(f"   sign-only upper edge (U->S) at a = {e_:.6g} ({c1_}->{c2_}); a/b = {e_/5e-4:.4g}")
    # E2: finer alpha=0.1 boundary in (beta, gamma), full model c=2, b=0, a=h=e=p=1
    print(" E2 full model c=2, alpha=0.1: class on a (beta, gamma) grid (U = symmetry-breaking real instability)")
    gr = [0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2]
    print("   beta\\gamma " + " ".join(f"{g:>5g}" for g in gr))
    for be in gr:
        print(f"   {be:>9g} " + " ".join(f"{classify(2.0, S0(al=0.1, be=be, ga=g))[0]:>5}" for g in gr))
    print(" E2b beta=gamma=x, alpha=0.1: class for x in", [0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
    print("   ", [classify(2.0, S0(al=0.1, be=x, ga=x))[0] for x in (0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)])
    print(" E2c alpha scan at beta=gamma=0.1, c=2:", [(al, classify(2.0, S0(al=al, be=0.1, ga=0.1))[0]) for al in (0.01, 0.1, 0.5, 0.9, 1, 1.5, 3, 10)])
    print(" E2d alpha=0.1, c=2: gamma at which U->S flips, for fixed beta (sign-based, bisection) and beta*gamma there")
    for be in (0.1, 0.2, 0.3, 0.5, 1.0, 2.0):
        f_ = lambda g, be=be: classify(2.0, S0(al=0.1, be=be, ga=g))[0]
        e_, c1_, c2_ = bis(f_, 0.03, 2.0)
        print(f"   beta={be:g}: {c1_}->{c2_} at gamma={e_:.4g}; beta*gamma={be*e_:.4g}")
    # E3: oscillatory (Hopf) region in S0 beta=gamma=0: first alpha with Re(lam)>0, as a function of c
    print(" E3 S0 (beta=gamma=0): smallest alpha on a log grid [1,1e5] with Re(lam_antisym)>1e-9, vs c (None = not found)")
    for c_ in (0.1, 0.5, 1, 1.2, 1.5, 1.8, 2, 3, 5):
        first = None
        for al in np.logspace(0, 5, 400):
            P = S0(al=al); r = racemic(c_, P)
            if r is None: continue
            la = np.linalg.eigvals(BA.T @ jac(r['x'], P) @ BA / 2)
            if la.real.max() > 1e-9: first = al; break
        print(f"   c={c_:g}: alpha_c(Hopf) = {first if first is None else round(first, 2)}")
    # F: ADDED AFTER UNBLINDING (after reading aped.py/aped.log); like-for-like probes against what aped.py reports
    print("\n=== F. added AFTER unblinding: like-for-like probes ===")
    for c_ in (2.0, 5.0):
        print(f" F1 selectivities (alpha,beta,gamma)=(0.35,0.2,0.3), S0 rates a=h=e=p=1, b=0, c={c_:g}:", classify(c_, S0(al=.35, be=.2, ga=.3))[:2])
    print(" F1 alpha=beta=gamma=0.1, S0 rates, c=2:", classify(2.0, S0(al=.1, be=.1, ga=.1))[:2])
    def lam_kR(c_, P, kR):  # monomer racemisation L<->D at kR (absent in the 2004 model): symmetric state unchanged, adds -2kR on [l-d]
        r = racemic(c_, P)
        if r is None: return None
        M_ = (BA.T @ jac(r['x'], P) @ BA / 2).copy(); M_[0, 0] -= 2*kR
        return np.linalg.eigvals(M_).real.max()
    print(" F2 PNAS model + monomer racemisation kR, S0 rates, alpha=beta=gamma=0.1, b=0, c=2: a-range with lambda>0 (400-pt log grid 1e-8..1.1)")
    for kR in (0.0, 1e-6, 1e-4):
        grid = np.logspace(-8, np.log10(1.1), 400); s_ = [lam_kR(2.0, S0(a=a_, al=.1, be=.1, ga=.1), kR) for a_ in grid]
        pos = [(g_, v_) for g_, v_ in zip(grid, s_) if v_ is not None and v_ > 0]
        print(f"   kR={kR:g}: lambda>0 for a in [{pos[0][0]:.3g}, {pos[-1][0]:.3g}]; max lambda {max(v_ for _, v_ in pos):.3g}; lambda(a=1e-8)={s_[0]:.3g}" if pos else f"   kR={kR:g}: never positive")
    print("\nDONE")
