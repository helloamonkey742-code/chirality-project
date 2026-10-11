"""Check that frank.ant (Antimony, converts to SBML) is the same chemistry as frank.py. (2026-10-11)

frank.ant writes frank.py's Frank-type scheme in a standard model language so other people can run it in their own
tools. This script converts it to SBML (writes frank.xml), simulates it deterministically with libRoadRunner, and
compares the crossing point (where k*A = kd) with frank.py's own code run with the noise switched off
(Omega = 1e30, one run).

Not in run_all.sh: needs antimony and libroadrunner, which the main environment lacks.
    python -m venv v && v/bin/pip install antimony libroadrunner numpy scipy && v/bin/python frank_sbml.py
Not checked: the stochastic part. frank.py adds counting noise with the chemical Langevin equation; an exact
Gillespie run of this SBML at frank.py's Omega = 1e4 would take ~1e7 events per run, too slow for 2000 runs.

PREDICTION (written 2026-10-11, before the first run):
P1. Antimony converts frank.ant to SBML without errors, and libRoadRunner loads the SBML.
P2. With noise off, the crossing values (A, x = (L+D)/2 and dA/dt at the crossing) from the SBML agree with
    frank.run within 1% for g = 0.003 and g = 0 (tau = 800); frank.py steps with Euler at dt = 0.01, which is
    the main source of any difference.
P3. With g = 0.003, L ends above D; with g = 0, L equals D to 1e-9 (relative).
"""
import antimony
import numpy as np
import roadrunner

import frank


def sbml():
    assert antimony.loadAntimonyFile("frank.ant") >= 0, antimony.getLastError()
    s = antimony.getSBMLString("frank")
    open("frank.xml", "w").write(s)
    return s


def crossing(rr, g, tau):
    """A, x, dA/dt where k*A first reaches kd, plus final L, D (same reading as frank.run)."""
    rr.resetAll()
    rr.g, rr.tau = g, tau
    rr.integrator.relative_tolerance, rr.integrator.absolute_tolerance = 1e-10, 1e-14
    t, A, L, D = rr.simulate(0, tau, int(tau / 0.01) + 1, ["time", "A", "L", "D"]).T
    i = int(np.argmax(rr.k * A >= rr.kd))
    w = max(len(t) // 50, 2)
    dadt = np.polyfit(t[i - w:i + w], A[i - w:i + w], 1)[0]
    return (A[i], (L[i] + D[i]) / 2, dadt), L[-1], D[-1]


def main():
    rr = roadrunner.RoadRunner(sbml())
    print(f"P1. SBML written to frank.xml; libRoadRunner loaded {len(rr.model.getReactionIds())} reactions")
    worst = 0
    for g in (0.003, 0.0):
        mine, L, D = crossing(rr, g, 800)
        _, theirs = frank.run(g, 1e30, 800, runs=1)
        rel = [abs(a - b) / abs(b) for a, b in zip(mine, theirs)]
        worst = max(worst, *rel)
        print(f"P2. g = {g}: SBML  A, x, dA/dt = " + ", ".join(f"{v:.5g}" for v in mine))
        print(f"    g = {g}: frank A, x, dA/dt = " + ", ".join(f"{v:.5g}" for v in theirs) +
              f"   (max rel diff {max(rel):.2e})")
        print(f"P3. g = {g}: final L = {L:.6g}, D = {D:.6g}, (L-D)/L = {(L - D) / L:.2e}")
        assert (L > D) if g else abs(L - D) / L < 1e-9
    assert worst < 0.01, "frank.ant does not reproduce frank.py"


if __name__ == "__main__":
    main()
