"""Crystal check: do published NaClO3 crystallization runs behave like the model's "decision count"?

The model's key quantity is how many independent units decide a hand in one volume.
For crystals it can be read straight off the data: if a beaker holds n_eff independent
decisions, each run's order O = (N_L - N_R)/(N_L + N_R) has mean square <O^2> = 1/n_eff
(plus bias^2). n_eff = 1 means the whole beaker picked one hand; n_eff = crystal count
means every crystal decided on its own.

This is a POST-HOC reanalysis of numbers quoted in the text of two papers. It is not a
pre-registered test (the student's own experiment would be). Expectations written
before running, after reading the papers:
  E1 Unstirred: n_eff well above 1 but possibly below the crystal count (clustering).
  E2 Stirred with the bar crushing crystals: n_eff ~ 1.
  E3 Floating bar (little crushing): n_eff of a few. Electro-convection (flow, no
     crushing): no symmetry breaking, i.e. transport without an amplifier does nothing.
     This is the model's point that BOTH amplification and transport are needed.
  E4 A small per-crystal bias is easy to see over thousands of independent crystals but
     invisible in ~20 one-decision runs (the model's sqrt(N) counting).

Sources (numbers quoted in text; figure data would need digitizing):
  K90: Kondepudi, Kaufman & Singh 1990, Science 250:975, doi:10.1126/science.250.4983.975
  M96: Martin, Tharrington & Wu 1996, PRL 77:2826, doi:10.1103/PhysRevLett.77.2826
       (NASA NTRS 19970000425 conference version, which carries the numbers below)
"""
from math import comb, sqrt


def n_eff(o_rms, bias=0.0):
    return 1 / (o_rms**2 - bias**2)


def binom_tail_le(k, n, p):
    return sum(comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k + 1))


rows = []
# K90 unstirred: 17 runs, 1000 crystals, 525 L. Per-run counts "statistically equal".
o = (525 - 475) / 1000
rows.append(("K90 unstirred", 1000 / 17, None, f"pooled O={o:+.3f}, {o*sqrt(1000):.1f} sigma"))
# K90 stirred 100 rpm: 32 runs, 11,829 crystals, every run >99% one hand.
rows.append(("K90 stirred 100 rpm", 11829 / 32, n_eff(0.99), "32/32 runs >99% one hand"))
# M96 quiescent: 94 runs, 5395 crystals, 2872 L; ensemble rms order ~0.25.
bias = (2872 - 2523) / 5395
rows.append(("M96 quiescent", 5395 / 94, n_eff(0.25, bias), f"pooled O={bias:+.3f}, {bias*sqrt(5395):.1f} sigma"))
# M96 5 cc, 1 Hz, bar crushes on the floor: 20 runs, 8 at O~+0.99, 12 at O~-0.95.
o_rms = sqrt((8 * 0.99**2 + 12 * 0.95**2) / 20)
rows.append(("M96 5 cc 1 Hz (crushing)", None, n_eff(o_rms), "8 L / 12 R"))
# M96 7.5 cc, 1 Hz, bar afloat: 92 runs, PDF ~flat -> <O^2> = 1/3, O_avg ~0.6.
rows.append(("M96 7.5 cc 1 Hz (afloat)", None, n_eff(sqrt(1 / 3)), "flat PDF"))

print(f"{'condition':28s} {'crystals/run':>12s} {'n_eff':>7s}  note")
for name, c, n, note in rows:
    print(f"{name:28s} {('%.0f' % c) if c else '-':>12s} {('%.1f' % n) if n else '-':>7s}  {note}")

# E4: the 6.5% left bias M96 saw over 5395 independent crystals would set each
# one-decision run's winner with p ~ 0.53. Over 20 runs that predicts ~10.6 L; 8 seen.
p = 0.5 + bias / 2
print(f"\nE4: p(L)={p:.3f}; P(<=8 L of 20)={binom_tail_le(8, 20, p):.2f} (p=0.5: {binom_tail_le(8, 20, 0.5):.2f})")
print(f"    runs needed to see this bias at 3 sigma: {9 / (bias**2):.0f}")

n_q, n_s, n_f = rows[2][2], rows[3][2], rows[4][2]
assert 5 < n_q < 5395 / 94, "E1: unstirred is many decisions but fewer than crystals"
assert n_s < 1.2, "E2: crushing + stirring gives one decision per beaker"
assert 2 < n_f < 5, "E3: floating bar gives a few decisions"
print("\nE1-E4 consistent.")
