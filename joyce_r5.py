"""Joyce et al. (1984) Nature 310:602, scorecard R5 (fast enough), scored against joyce_r5_predictions.txt.

Predictions hashed before the paper was opened:
  sha256 2610b2a109952e05acf57e4ea0107f66179cdec72b9e14b2e3ae48bbf348ac1f, 2026-10-20 23:57 PDT (2026-10-21 06:57 UTC)
Numbers entered by a separate reader (blind to the predictions), then checked here against the PDF text layer.
Table 1 (p. 602): % of 14C-labelled G_D residues in oligomers of length >= 4, net of template-free controls.
All runs on poly(C_D), 0 C, 7 days (Fig. 2 and Fig. 3 legends).  System I: 2-MeImpG, Mg2+.  System II: ImpG, Mg2+, Zn2+.
n >= 4 is the only printed fraction, so it is f (the predictions' conservative rule): k2 is a lower bound.
DL concentration is taken as total (not stated in the paper); per-enantiomer would halve k2.

Run: python3 joyce_r5.py > joyce_r5.log
"""
import hashlib
import math

PRED_SHA = "2610b2a109952e05acf57e4ea0107f66179cdec72b9e14b2e3ae48bbf348ac1f"
FLOOR = 3e-7                  # scorecard.py spec, pond and lake, /M/s
T = 7 * 86400                 # s
TABLE1 = [                    # system, monomer, c (M), f (n >= 4)
    ("I", "D", 0.025, 0.74), ("I", "D", 0.05, 0.75), ("I", "DL", 0.05, 0.32), ("I", "DL", 0.10, 0.35),
    ("II", "D", 0.01, 0.30), ("II", "D", 0.02, 0.19), ("II", "DL", 0.01, 0.02), ("II", "DL", 0.02, 0.02),
]


def k2(c, f):
    return -math.log(1 - f) / (c * T)


def main():
    with open("joyce_r5_predictions.txt", "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == PRED_SHA, "predictions file changed"
    assert abs(k2(0.05, 0.75) - 4.58e-5) < 1e-7  # hand check: ln 4 / (0.05 * 604800)
    k = {"D": [], "DL": []}
    for sysn, mono, c, f in TABLE1:
        k[mono].append(k2(c, f))
        print(f"System {sysn:2} {mono:2} {c:5.3f} M  f {f:.2f}  k2 {k2(c, f):.2e} /M/s  floor x{k2(c, f) / FLOOR:5.0f}")
    kD, kR = min(k["D"]), min(k["DL"])
    letter = "Y" if kD >= FLOOR and kR >= FLOOR else "P" if kD >= FLOOR else "N"
    print(f"\nmin k2(D) {kD:.2e}, min k2(DL) {kR:.2e} /M/s; floor {FLOOR:.0e}; per-enantiomer DL basis: {kR / 2:.2e}")
    cs = [c for _, _, c, _ in TABLE1]
    res = [
        ("P1", True, "no rate constant, half-life or time course printed (one 7-day point)"),
        ("P2", all(0.02 <= c <= 0.2 for c in cs), f"0 C, 7 d, Mg2+ held; c spans {min(cs)}-{max(cs)} M, System II 0.01 M is below 0.02"),
        ("P3", all(1e-6 <= x <= 1e-3 for x in k["D"]), f"k2(D) {min(k['D']):.1e}-{max(k['D']):.1e}"),
        ("P4", all(max(f for s, m, _, f in TABLE1 if s == sy and m == "DL") < min(f for s, m, _, f in TABLE1 if s == sy and m == "D")
                   for sy in ("I", "II")) and kR >= FLOOR,
         f"DL yields 32-35% vs D 74-75% (I), 2% vs 19-30% (II); min k2(DL) {kR:.1e}"),
        ("P5", letter == "Y", f"R5 = {letter}"),
        ("P6", True, "only pure D, pure L and 1:1 DL monomer used; no starting excess, no amplification reported"),
    ]
    for name, ok, why in res:
        print(f"{name} {'HELD  ' if ok else 'FAILED'} {why}")
    print(f"{sum(ok for _, ok, _ in res)} of {len(res)} held.  R5 = {letter} (scorecard.py not changed).")


if __name__ == "__main__":
    main()
