# Small open items H12, H26, H28 (2026-10-12)

Three items from the SUBMISSION NOTES of `journals/ija_draft.md` (H12, H26, H28). This pass did not edit `ija_draft.md`. The caption and text changes below are proposals and have not been applied. Every run used `/opt/anaconda3/bin/python` under `nice -n 10` with `OMP_NUM_THREADS=2` and `NUMBA_NUM_THREADS=2`.

## H12. Kramers flip-size message in review.py R4

**What changed.** Only the closing message of `r4_longrun()` in `review.py` was edited. The old line printed that a 1 µm³ droplet at 1 µM (about 600 molecules) has "flips possible". The new lines print the flip criterion, the bound at three values of kτ, and the expected number of flips for N = 600, with an assert that this number is below 10⁻¹⁰⁰. The criterion is (√2/2π)·e^(−N/2)·kτ ≥ 1, so N ≤ 2 ln(0.225 kτ).

The rerun output is `review.log`; the old log is kept as `review_pre_r4fix.log`. `diff review_pre_r4fix.log review.log` shows:

- **R4 (intended):** line 42 is replaced by six lines:
  - N ≤ 1.6 at kτ = 10
  - N ≤ 20.0 at kτ = 10⁵
  - N ≤ 57.4 at kτ = 1.3×10¹³
  - "2e-118 expected flips at k tau = 1.3e13 (never flips)"
  - "only pools of tens of molecules can flip; a droplet, pond or ocean never does."

  The Kramers rate check is unchanged: measured 3.4e-03 vs 3.5e-03, and 1.4e-03 vs 1.5e-03.
- **R6 (not intended, cosmetic):** line 60 ends in `final -0.000` where the old log had `final +0.000`. This is the sign of a closed-network ee that rounds to zero at 3 decimals. The R6 code was not touched, so the change comes from floating-point differences between this run and the 2026-10-01 run. The R6 assertions still pass.
- **Format only:** the old log ended with an `exit 0` line written by the run wrapper. The same line was appended to the new log after this run exited 0.

**Proposed text.** Submission-note item 16 already contains the replacement sentence for Section 3.8. Once it is applied, change its source from `review.log:42` to `review.log:42–47`, and mark H12 resolved.

## H28. Larger 1D Kibble–Zurek system

**What changed.** `kz_large.py` is a new, pre-registered script. It imports `kz_patches` unchanged (`check`, `batch`, `report`, `fit_probit`) and runs one 1D periodic system:
- size: L = 1024ξ, 6476 cells, about 80 patches, which is twice the largest earlier run;
- model: the same as `kz_patches.sec_1d` (D = 4, γ = 0.01, ε = 10⁻⁴ per cell, dt = 0.05);
- Δ design: Δ_total = 1.0, 1.5 and 2.0, with 1200 runs each, chosen near the probit information peak;
- seeds 400–405, run in 2 threads.

The nine predictions were committed before the run (`a8def72`, "kz_large: predictions before run"). The results and the log were committed afterwards (`9de3aef`). The RESULTS block is appended to the script's docstring.

**Result.** The end-of-sweep value is c = 0.918 ± 0.021 (the early value is 0.930 ± 0.022). All 9 predictions held:

| # | Prediction | Outcome |
|---|---|---|
| P1 | 3-cell uniform limit within 3 SE | +1.13 SE (0.8478 ± 0.0057 vs 0.8413) |
| P2 | within 2 combined SE of the 0.879 ± 0.021 plateau | +1.31 |
| P3 | this system alone excludes 0.798 by > 3 SE | 5.65 SE |
| P4 | within 2 SE of the Gaussian-field 0.914 | +0.19 SE |
| P5 | early = end majority ≥ 98% | 99.22% |
| P6 | no single-hand runs | 0.00 at every Δ |
| P7 | each P(Δ) within 2 SE of Φ(cΔ) | −0.73, +0.32, +0.30 SE |
| P8 | inverse-variance pool of 128/512/1024ξ excludes 0.798 by > 4 SE | c = 0.898 ± 0.015, 6.77 SE |
| P9 | under 20 min | 619 s (10:19 wall, 196% CPU) |

**What it means.** The H28 concern was that 0.80 was excluded at 3.95 SE only through a pooling step that was not pre-registered. It is now excluded at 5.65 SE by one pre-registered system with no pooling, and at 6.77 SE by a pool that was pre-registered.

The three sizes give 0.876, 0.882 and 0.918. That sequence passes both P2 (plateau) and P4 (theory), so this run cannot tell a flat plateau from a slow drift toward 0.914. The 1024ξ value by itself sits on the Gaussian-field theory.

The pool in P8 combines two Δ designs, those of `kz_patches.py` and of this script; P3 does not depend on the pool. For H20 (propagating c² into the sweep), a c nearer 1 makes the c² correction smaller (c² = 0.84 against 0.77 before).

**Proposed text (not applied).**
- *Table 2, patch-sum row (draft line 519):* change "c = 0.88 ± 0.02 (1D)" to "c = 0.92 ± 0.02 (1D, 1024ξ; 0.88 ± 0.02 at 128–512ξ)".
- *Section 2.3 / 4, wherever c in 1D is quoted:* "In 1D the patch-sum constant is c = 0.918 ± 0.021 for a single 1024ξ system (about 80 patches; `kz_large.py`, pre-registered). This excludes the box patch-sum value 0.80 by 5.7 SE and agrees with the Gaussian-field value 0.914. Smaller systems gave 0.876 ± 0.029 (128ξ) and 0.882 ± 0.029 (512ξ). The pre-registered pool of all three sizes is 0.898 ± 0.015."
- *Submission-note item 49 (line 135):* replace "4.0 SE above 0.80 when pooled (pooling not pre-registered), 2.9 SE at 512 ξ alone" with "5.65 SE above 0.80 from a single pre-registered 1024ξ system (`kz_large.log`, P3); 6.77 SE for the pre-registered 128/512/1024ξ pool (P8)".
- *Submission-note item 51 (pre-registration tallies, Section 2.7):* add "kz_large 9/9 (`kz_large.py` RESULTS)".
- *Submission notes:* mark H28 resolved. In the Kibble–Zurek line of the status list (line 216), change "Open: H20, H28." to "Open: H20."

## H26. Figures 1, 3 and 6 redrawn

New scripts write new PNGs; `threshold.png`, `rnaas_fig.png` and `journals/fig_schematic.png` are untouched. All three new PNGs are 300 dpi. Each script asserts that its plotted numbers match the log lines. The two shares and threshold checks were confirmed to fail when one number is altered.

| Figure | Script | Output | Self-check |
|---|---|---|---|
| 3 | `make_fig_threshold2.py` (imports `ocean.py` curves) | `threshold2.png` | grey band = `sweep2.log` `g_min 4.40e-18` / `g_max 2.64e-16`; each curve at both band edges equals `ocean.c_size`/`c_time` and gives Δ ≥ 2 |
| 6 | `make_fig_shares2.py` | `shares2.png` | all 14 Table 3 worlds re-read from `multiwall.log` Section 4 ("P_pved %, no lid / final": multiwall and cascade ×10 columns), exact to 0.1 |
| 1 | `journals/make_fig_schematic2.py` (copy of `make_fig_schematic.py`) | `journals/fig_schematic2.png` | g range = `sweep2.log` g_min/g_max, and `n_needed(g, 10)` at both ends = the log's 1.84e+34 and 5.12e+30 |

**Text inside the images.**
- Figure 3: the title "Where physics, not chance, picks life's hand" becomes "Minimum concentration for PVED selection (k2 = 1e-3 /M/s)". The band label becomes "published amino-acid PVED range (0.10–6.2 aeV)".
- Figure 6: the x-axis "…hand wins (%)" becomes "Draws in which the PVED-favoured hand is selected (%)". Each bar runs from the wall-time ×10 share to the as-computed share. Lid worlds (Europa, Enceladus, Dione, Mimas, Miranda) get a hatched no-lid bar. Each world is labelled with Table 3's rounded range, rounded half-up, so Triton reads 27–37%.
- Figure 1: panel (c) now spans 4.4×10⁻¹⁸ to 2.64×10⁻¹⁶. The panel already had no anthropomorphic words ("outcome set here", "favoured hand").

None of the three scripts contains "decides", "wins", "chooses", "takes over", "helps" or "picks".

**Proposed caption updates (not applied).**
- *Figure 1, panel (c):* replace "over the g range of the uniform-D_z base model (4×10⁻¹⁹ to 4×10⁻¹⁶; not yet redrawn for the published amino-acid range of Table 1, 4.4×10⁻¹⁸ to 2.64×10⁻¹⁶)" with "over the published amino-acid range of Table 1 (4.4×10⁻¹⁸ to 2.64×10⁻¹⁶)". Change the image path to `fig_schematic2.png` and the last sentence to: "Drawn by `journals/make_fig_schematic2.py`, which uses `review.py: n_needed()` and asserts the 3.6×10³³ and 6×10²⁶ values and the g range and N values of `sweep2.log`."
- *Figure 3:* replace "The grey band is the g range of the uniform-D_z base model (4×10⁻¹⁹ to 4×10⁻¹⁶), not yet redrawn for Table 1's range." with "The grey band is the published amino-acid range of Table 1 (4.4×10⁻¹⁸ to 2.64×10⁻¹⁶, 0.10–6.2 aeV)." Change the image path to `../threshold2.png`.
- *Figure 6:* replace the caption with: "**Figure 6.** Percentage of 200,000 plausible input draws in which the PVED-favoured hand is selected, for each world in Table 3 (`multiwall.log`, Section 4). Each bar runs from the share with wall rise and sink times multiplied by 10 (left end) to the share with those times as computed (right end). Solid bars include the lid test. Hatched bars show the share if no lid formed, for the five worlds where a lid can form. Colour shows the ocean status. Europa is shown with 20 km of ice. Drawn by `make_fig_shares2.py`." Change the image path to `../shares2.png`.
- *Code availability table:* Figure 1 → `journals/make_fig_schematic2.py`; Figure 3 → `make_fig_threshold2.py` (new row); Figure 6 → `make_fig_shares2.py`.
- *Submission notes:* mark H26 resolved, and in item 11 strike the threshold.png and rnaas_fig.png entries (coarsening.png remains). The Figures row of the journal-requirements table still says Figure 1 is drawn at 600 dpi practice; the new PNGs are 300 dpi, which meets the halftone and colour requirement but not the 600 dpi mixed-art one. Raise `dpi` if the journal treats these as mixed art.
