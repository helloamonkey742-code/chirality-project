"""POST HOC arithmetic (no simulation, no pre-registration): the share of runs that rule B (journals/droplet_prereg.md Section 2 item 1) calls
NOT SCORABLE only because more than 10% of the 128 control drops are read wrong, i.e. P(k > 12.8) = P(k >= 13) for k ~ Binomial(128, e),
when every control is misread with the same true probability e.  Written after droplet_range_b_ns.log showed NS 0.04-0.05 at e = 0.06.
Compare with misread_sim.log, whose policy B has no such cap (its STEP shares at e = 0.08 and 0.10 are conditional on no cap).

Run: python3 droplet_range_cap.py > droplet_range_cap.log
"""
from math import comb


def tail(n, e, kmin):
    return sum(comb(n, k) * e ** k * (1 - e) ** (n - k) for k in range(kmin, n + 1))


assert abs(tail(2, 0.5, 1) - 0.75) < 1e-12 and abs(tail(3, 0.5, 3) - 0.125) < 1e-12      # hand cases
print("true misread e -> P(NOT SCORABLE from the control cap), 128 controls, k >= 13")
for e in (0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10):
    print(f"  e = {e:.2f}: {tail(128, e, 13):.3f}")
