# Fresh-clone reproducibility check (2026-10-13)

Question: can a reviewer reproduce the paper from a fresh clone of
https://github.com/helloamonkey742-code/chirality-project.git?

**Short answer:** `./run_all.sh` finishes from a fresh clone, and every number it prints matches the committed
logs. The same is true for `--full`. But run_all.sh does not run the scripts behind the paper's main results
(`sweep2`, `sweep3`, `multiwall`, `seedfix`, `rerun_base`, `tidal_mix`, ...). Those need `numba`, which the
README and the paper's code-availability statement do not list. So "run_all.sh reproduces the paper" is not
true yet. The fixes are listed at the end.

## Setup

- Clone: `git clone https://github.com/helloamonkey742-code/chirality-project.git` into a new scratch directory.
  HEAD `ff06793` is the same as the local `~/Documents/chirality-project` HEAD and `origin/main`.
- Interpreter: `/opt/anaconda3/bin/python` (Python 3.13.5, numpy 2.1.3, scipy 1.15.3, matplotlib 3.10.0,
  numba 0.61.0; `gsw`, `antimony` and `roadrunner` are not installed). Nothing was pip-installed.
- What a reviewer is told to install: README ("Needs numpy, scipy, matplotlib.") and the paper's Code availability
  ("It needs Python with numpy, scipy and matplotlib."). The repo has no requirements.txt or environment file,
  and no versions are pinned.

## Commands and wall time

Fast set (the acceptance run):

```
PATH="/opt/anaconda3/bin:$PATH" OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
  NUMBA_NUM_THREADS=2 MPLBACKEND=Agg nice -n 10 ./run_all.sh > run_all_fast.out 2>&1
```

Exit 0, `== all checks passed (35 scripts)`, **wall time 132 s**. The README estimate is 4–7 min.

Slow (`--full`) set: run one script at a time with the same env, because the README estimates add up to more than
60 min. Command: `nice -n 10 /opt/anaconda3/bin/python <script>`. Total **3024 s (50 min)**. Nothing was skipped.

Paper figure scripts: `journals/make_fig_schematic2.py`, `journals/make_fig_sensitivity3.py`,
`make_fig_shares3.py`, `make_fig_threshold2.py`. All exit 0, take about 1 s each, and pass their self-checks.

## Per-script results

How outputs were compared: run_all's stdout was split by `== script.py` header (Warning lines removed, as run_all
does) and diffed against the committed `<script>.log`. Files the scripts write themselves (PNGs, `soai_check.log`,
`viedma_check.log`) were compared with `git status` / `git diff` in the clone.

**Fast set (35 scripts): all passed.**

| Result | Scripts |
|---|---|
| Identical to committed `.log` | assumed_rules, breakeven, crystal_check, earth_breakeven, impact_melt, jupiter_moons, nudge_realism, ocean_worlds, robust_ranking, soai_check, strat_sweep, titan_95k, titan_rates, uncertainty, viedma_check |
| Identical to their section in the committed `run_all.log` | sim, ocean, domains, patches, scenarios, inherit, dual, scorecard |
| Numeric difference (sign of zero only) | quench: committed `final +0.000` vs new `final -0.000`, on 2 lines (20% seed, slow=1 and slow=0.01). The fade ratios (1.6x, 1.2x) are the same. |
| Text changed since `run_all.log` was written; the shared numbers are the same | design: the headers changed (`smallest detectable k2 ...` became `=== Baseline Detection Limits ...`) and the Arm F sections were added. The detection limits and floors 3e-07 / 3e-09 are the same. |
| No committed log to compare against (passed their own asserts) | dna_network, stress_test, chiral_search, sketch, prebiotic_search, test_candidate_chemistry, extensive_testing, lethal_stress_tests, planetary_chirality_audit, resolve_challenges |

All 13 PNGs that the fast set rewrites came out byte-identical (`git status` is clean after the run).

**Slow set (12 scripts): all passed (exit 0).**

| Script | Time | Diff vs committed log |
|---|---|---|
| review.py | 124 s | identical except the trailing `exit 0` line, which the author's wrapper added |
| k3.py | 79 s | same (only `exit 0`) |
| frank.py | 150 s | same (only `exit 0`) |
| frank2.py | 126 s | same (only `exit 0`) |
| objections.py | 66 s | identical |
| moon_nudges.py | 24 s | identical |
| coarsen_fig.py | 160 s | same (only `exit 0`); `coarsening.png` byte-identical |
| aniso3d.py | 186 s | same (only `exit 0` and a hand-added dated note) |
| domains3d.py | 238 s | numbers identical; hand-annotated lines differ (`spread 31% (recomputed ...)` vs `spread 31%`, and an `exit 0 (verify.log ...)` line) |
| coarsen3d.py | 435 s | same (only `exit 0`) |
| shell3d.py | 1033 s | same (only `exit 0`) |
| sphere_layers.py | 403 s | sections 1, 2 and 3: every number identical; only the `[NNs]` timings and the `real/user/sys` lines differ. **Section 2b, which is in the committed log, is not produced**: with no arguments the script runs `["1","2","3"]`, and run_all.sh passes none. |

Paper figures: `fig_schematic2.png`, `fig_sensitivity3.png`, `shares3.png` and `threshold2.png` were regenerated
byte-identical. Note that these scripts re-plot numbers parsed from the committed logs. They do not recompute them.

**Errors: none. Missing outputs:** sphere_layers section 2b, plus everything in the "not run" list below.

## Other checks

1. **The paper's main results are not in run_all.sh, even with `--full`.** 30 scripts named in the Code
   availability table of `journals/ija_draft.md` are never run by run_all.sh: sweep2, sweep3, multiwall, seedfix,
   rerun_base, tidal_mix, tidal_mix_g005, front2d, pin_onset, gradient_ramp, kz_patches, kz_large, handoff,
   deephand_scale, rho2side, two_zone, polar_gap, currents, realistic, earth, field_walls, field_flips, wobble,
   wobble_sweep, verify_seeds, verify_sstar, plus the 4 figure scripts. `multiwall.log` (32 citations) and
   `seedfix.log` (28 citations) are the logs the paper cites most.
2. **Undeclared dependencies.**
   - `multiwall.py:110` runs `from numba import njit` at top level. deephand_scale, handoff, rerun_base,
     rho2side, seedfix, tidal_mix and tidal_mix_g005 all `import multiwall`, so they need numba too.
   - `verify_sstar.py` needs `gsw` (TEOS-10). It is not installed here, and its docstring points at
     `~/research-env/bin/python`.
   - `frank_sbml.py` needs `antimony` + `libroadrunner`. The README already says so.
3. **Hard-coded paths.**
   - `run_all.sh:7-8` prefers `/opt/anaconda3/bin/python` over `python3` whenever it exists. On this machine that
     overrides the caller's PATH. This is harmless for a reviewer, but `PYTHON=...` is the only override and the
     README does not mention it.
   - `/opt/anaconda3/bin/python` appears in the docstrings/usage lines of kz_large, seedfix, tidal_mix,
     tidal_mix_g005, make_fig_shares2, make_fig_shares3, make_fig_threshold2, check_refs, multiwall,
     journals/make_fig_schematic*.py. Everything still runs, since these are comments only.
   - No code path reads `/Users/...`.
4. **Files referenced but not in the repo.** `journals/ija_draft.md:169` cites `refcheck.py` (marked "scratch"),
   which is not in the repo (the repo has `check_refs.py`). `OUTREACH_EMAILS.md` is gitignored but only mentioned
   in AUTONOMOUS_LOG.md. Every other `*.py/.log/.png/.md` named in README.md, ija_draft.md and run_all.sh exists.
5. **README vs run_all.sh.**
   - The committed `run_all.log` (2026-09-23) covers only 9 scripts. run_all.sh now runs 35.
   - README:263 says `--full` takes "about 22 minutes". The comment at run_all.sh:2 says "~30 min", and run_all.sh:17
     says sphere_layers alone takes ~22 min. Measured here: 50 min for the slow set, with 2 threads.
   - README:669 says "fast (~4–7 min)". Measured here: 2.2 min.
6. **Fragility.** In `run_all.sh:28`, `"$PY" "$s" 2>&1 | grep -v Warning` under `set -o pipefail` aborts the whole
   run if a script ever prints only Warning lines, because grep then exits 1. It also hides any real output line
   that contains "Warning".

## Fixes for the author

1. **Add `requirements.txt`** with `numpy`, `scipy`, `matplotlib`, `numba` (pin the versions you used:
   numpy==2.1.3, scipy==1.15.3, matplotlib==3.10.0, numba==0.61.0; Python 3.13). Optionally add `gsw` for
   verify_sstar, and `antimony` and `libroadrunner` for frank_sbml.
2. **README.md:671**: change "Needs numpy, scipy, matplotlib." to "Needs numpy, scipy, matplotlib and numba
   (multiwall.py and the scripts that import it); gsw for verify_sstar.py; see requirements.txt."
3. **journals/ija_draft.md:782** (Code availability): change "It needs Python with numpy, scipy and matplotlib" to
   "numpy, scipy, matplotlib and numba (gsw for the TEOS-10 check)". Also remove the claim that
   `./run_all.sh --full` covers the slow runs, or make it true (fix 4).
4. **run_all.sh**: add a third tier, e.g. `paper=(sweep2.py sweep3.py multiwall.py seedfix.py rerun_base.py
   tidal_mix.py tidal_mix_g005.py front2d.py pin_onset.py gradient_ramp.py kz_patches.py kz_large.py handoff.py
   deephand_scale.py rho2side.py two_zone.py polar_gap.py verify_seeds.py ...)` behind `--paper`, using the section
   arguments each script's docstring documents (e.g. `multiwall.py <n>`). Add the 4 figure scripts to the fast
   set, since they take ~1 s each. If some of these are too slow, list them in run_all.sh with their run time and
   the exact command, the way `shell_res.py` is listed.
5. **run_all.sh:17 / sphere_layers.py**: call `sphere_layers.py 1 2 2b 3` so the committed log's section 2b is
   regenerated.
6. **run_all.log**: regenerate it with the current 35-script fast set. The committed copy is from 2026-09-23 and
   has 9 scripts.
7. **README.md:263**: change "about 22 minutes" to "about 50 minutes (sphere_layers alone ~6–22 min)". Make
   run_all.sh:2 say the same.
8. **quench.py**: print `final {abs(x) if abs(x) < 5e-4 else x:+.3f}` (or add `+ 0.0`) so that `-0.000` and
   `+0.000` do not differ between machines.
9. **run_all.sh:28**: replace `| grep -v Warning` with `PYTHONWARNINGS=ignore "$PY" "$s"` (or `|| true` on the
   grep only), so the pipeline cannot fail on an all-warnings output.
10. **journals/ija_draft.md:169**: either commit `refcheck.py` or point to `check_refs.py`.
11. Optional: document `PYTHON=/path/to/python ./run_all.sh` in the README. Also drop the `/opt/anaconda3` and
    `~/research-env` paths from the docstrings, or note that they are the author's local paths.
