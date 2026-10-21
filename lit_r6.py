"""Literature R6 re-check under the seed rule: scores A1-A2 and K1-K2 of lit_r6_predictions.txt (hashed first).

A1-A2: aped.py's racemic-state stability at the settings' total concentrations.  aped.py has no __main__ guard, so
its definitions (everything above its first print) are exec'd here and nothing in it is edited.
K1-K2: scorecard scores and ranks with the R6 letters of rules R and L (scorecard.py is imported, not edited).

Run:  UV_CACHE_DIR=$HOME/Library/Caches/uv uv run -q --no-project --with numpy --with scipy --with matplotlib \
        python lit_r6.py > lit_r6.log
"""
import hashlib
import time

from scorecard import CANDIDATES, KEYS

PRED, SHA = "lit_r6_predictions.txt", "f0345e5059ec0c1de6124d19f5f6c54435d6ace047840c47d0c908f8901238ed"
assert hashlib.sha256(open(PRED, "rb").read()).hexdigest() == SHA, "predictions file changed after hashing"
print(f"{PRED} sha256 {SHA} (matches)")

# --- A1-A2: aped.py definitions only (ponytail: exec of a source slice; breaks if aped.py's header print moves)
src = open("aped.py").read()
cut = src.index('print("APED racemic-state stability')
ns = {"__name__": "aped_defs"}
exec(compile(src[:cut], "aped.py", "exec"), ns)
BASE, smax = ns["BASE"], ns["smax"]

t0 = time.time()
s05, x05 = smax(dict(BASE, ctot=0.5))
a2 = f"{s05:.4e}" == "-9.2864e-07"
print(f"\n[A2] ctot = 0.5 M: sigma_max = {s05:+.4e} /s at x = {x05:.3e} M (aped.log:91: -9.2864e-07): "
      f"{'HELD' if a2 else 'FAILED'}")
a1 = []
for c in (0.1, 1e-2, 1e-4, 1e-6):
    s, x = smax(dict(BASE, ctot=c))
    a1.append(s <= 0)
    print(f"[A1] ctot = {c:<6g} M: sigma_max = {s:+.4e} /s at x = {x:.3e} M  ({'stable' if s <= 0 else 'UNSTABLE'})")
print(f"[A1] racemic stable at every drive at all four: {'HELD' if all(a1) else 'FAILED'}  ({time.time() - t0:.0f} s)")

# --- K1-K2: letters as predicted (current -> R, L); proposals after the evidence go in the note, not here
LIT = {  # scorecard row: (R, L)
    "Ghadiri peptide replicator (2001)": ("N", "N"),
    "Joyce RNA cross-inhibition (1984)": ("N", "N"),
    "Soai reaction (1995; 2003)": ("N", "P"),
    "Viedma ripening, NaClO3 (2005)": ("N", "P"),
    "Attrition deracemization, aa deriv (2008)": ("N", "P"),
    "Magnetite crystallization (2023)": ("N", "N"),
    "Blackmond eutectic (2006)": ("N", "P"),
    "APED peptides via NCAs (2004-2010)": ("N", "P"),
}
assert set(LIT) <= set(CANDIDATES)
GATE = {"Y": 0, "P": 1, "N": 2, "?": 1}


def table(r6):
    rows = {n: dict(m, R6=r6.get(n, m["R6"])) for n, m in CANDIDATES.items()}
    score = {n: sum({"Y": 1, "P": 0.5}.get(m[k], 0) for k in KEYS) for n, m in rows.items()}
    return sorted(rows, key=lambda n: (GATE[rows[n]["R8"]], -score[n])), rows, score


for label, i in (("now", None), ("rule R", 0), ("rule L", 1)):
    ranked, rows, score = table({} if i is None else {n: v[i] for n, v in LIT.items()})
    print(f"\n{label}:")
    for n in ranked:
        print(f"  {n:47} R8 {rows[n]['R8']}  R6 {rows[n]['R6']}  {score[n]:.1f}")
    if label == "rule R":
        g, j = score["Ghadiri peptide replicator (2001)"], score["Joyce RNA cross-inhibition (1984)"]
        k1 = g == j == 5.5 and ranked[0].startswith("Ghadiri")
        print(f"[K1] Ghadiri {g}, Joyce {j}, first = {ranked[0]}: {'HELD' if k1 else 'FAILED'} (tie broken by dict order)")
    if i is not None:
        block = [(n.split(" (")[0], score[n]) for n in ranked if rows[n]["R8"] == "N"]
        print(f"[K2] R8 = N block: {block}")
