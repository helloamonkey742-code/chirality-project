"""Figure 5 (seed-fixed): make_fig_sensitivity2.py's tornado plot with the Enceladus rows corrected.

rerun_base.log is not rewritten. seedfix.log Section F lists every changed rerun_base.log line as an
"old:"/"new:" pair; this script checks each "old:" against rerun_base.log, swaps in the "new:" text, and
plots the result with make_fig_sensitivity2's own parser and self-checks. No new runs.
"""
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_fig_sensitivity2 as m2  # noqa: E402

ROOT = os.path.join(HERE, "..")


def patched_lines():
    log = open(os.path.join(ROOT, "rerun_base.log")).read().splitlines()
    sf = open(os.path.join(ROOT, "seedfix.log")).read().splitlines()
    start = next(i for i, l in enumerate(sf) if l.startswith("=== F. rerun_base"))
    n = 0
    for i in range(start, len(sf)):
        m = re.match(r"^  \[.*\] rerun_base\.log:(\d+)$", sf[i])
        if not m:
            continue
        k = int(m.group(1)) - 1
        old, new = sf[i + 1], sf[i + 2]
        assert old.startswith("    old: ") and new.startswith("    new: "), sf[i]
        assert log[k] == old[9:], (k + 1, log[k], old)
        log[k] = new[9:]
        n += 1
    return log, n


def main():
    log, n = patched_lines()
    assert n == 95, n  # 90 Enceladus + 4 aggregate + 1 flagged line (seedfix.log Section F header)
    with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as f:
        f.write("\n".join(log) + "\n")
    m2.LOG, m2.OUT = f.name, os.path.join(HERE, "fig_sensitivity3.png")
    w = m2.parse(m2.LOG)
    assert w["Enceladus", "main"]["base"] == 26.4 and w["Europa", "main"]["base"] == 39.6
    m2.main()
    os.unlink(f.name)
    print("patched", n, "rerun_base.log lines from seedfix.log Section F")


if __name__ == "__main__":
    main()
