# Verification record — 2026-10-06

What was checked for the numbers in RNAAS_DRAFT.md / PAPER_DRAFT.md §3.14–3.15.

| Check | How | Result |
|---|---|---|
| Reproducible | Re-ran `strat_sweep.py`, `jupiter_moons.py`, `ocean_worlds.py` from scratch; diffed against committed logs | Identical |
| Not a lucky random seed | `verify_seeds.py`: seeds 1–5 | Every world within 0.3 points of the seed-0 value; all ordering checks pass |
| Salt thresholds (S\*) | `verify_sstar.py`: recomputed from TEOS-10 (gsw 3.6.23) without project code | Match to rounding of ice pressure (Europa 3.2 vs 3.0 g/kg at 24.2 MPa, the pressure under 20 km of ice; S\* falls steeply here, to 2.55 at the ~25 MPa that strat_sweep.py's comment rounds to; Dione 5.75 vs 5.6); no S\* above 28–29 MPa, as claimed. Needs a Python with gsw (not /opt/anaconda3) |
| Scripts outside `run_all.sh` (2026-10-09) | Re-ran earth, field_flips, field_walls, currents, polar_gap, two_zone, wobble, wobble_sweep, realistic, make_rnaas_fig, verify_seeds, verify_sstar; diffed against committed logs | All identical apart from run times; `rnaas_fig.png` byte-identical; verify_seeds: all 5 seeds ok, every world within 0.3 points |
| Core model suite | `./run_all.sh` (fast set, now 22 scripts; needs `PYTHON=/opt/anaconda3/bin/python`) | All self-checks pass |
| References | Crossref lookup of each DOI in the RNAAS draft | All resolve to the cited paper |
| References, all documents (2026-10-09) | Crossref lookup of all 149 DOIs in PAPER_DRAFT, RNAAS_DRAFT, README, EXPERIMENT, vertical_mixing_lit and script docstrings; arXiv API for the three 2026 preprints cited without a DOI | 148 resolve to the cited paper. One was wrong: Charette & Smith 2010 was cited as 10.5670/oceanog.2010.09, which is Yoder et al. 2010; fixed to 10.5670/oceanog.2010.51. Five older references now carry their checked DOIs. `pdf/chirality_paper_draft.pdf` keeps the wrong DOI until it is rebuilt |
| Ranking vs newest measurements | `robust_ranking.py`: Europa ice from Juno (29 ± 10 km, Levin et al. 2025); Enceladus salinity 5–20, 7–30, 17–40 g/kg (Kang et al. 2022) | Predictions held. Europa 26% vs Ganymede 28% is a tie, not a clear lead; Enceladus is last under every estimate (0–18%) |
| What the unknowns would have to be (2026-10-09) | `breakeven.py`, `assumed_rules.py`, `earth_breakeven.py`: predictions written in each docstring before the first run; each first reproduces `uncertainty.py` / `earth.py` exactly (asserted) | All 13 predictions held. Moons pass 50% only through vertical mixing (D_z ≥ 2.7×10⁻⁵ / 6.1×10⁻⁶ m²/s); the two assumed rules move a moon ≤ 2.2 points; a 10–100× bias adds 2–6 points; Earth's field must make < ~2×10⁻¹⁵ of the deciding molecules for the weak force to keep deciding |

Now in `run_all.sh`: the three ocean-world scripts (previously not run by it), `robust_ranking.py`, the three lab-data checks (`crystal_check.py`, `soai_check.py`, `viedma_check.py`) and the three studies above.

## Independent rebuild from the paper text alone (2026-10-06)

A separate agent rebuilt the model from PAPER_DRAFT.md and README.md without opening any code or log.
- **Matched:** the molecule-count rule (§3.0: 3.6×10³³, 18/1.8/0.4 km cubes, ponds ~7 orders short, temperature scaling); the layer-merge shares (33%/34%); the layered shares (96%/92%); Europa's lid share (24%); Enceladus 0% with the lid; the never-merge bounds (3%/7%).
- **Did not match:** Enceladus 22% → 29.5%, Europa 29% → 33.8%, Ganymede 28% → 32.8%. Early Earth ranged 32–96% depending on guesses.
- **Cause:** the text never gave the concentration range, the bias range, Earth's mixing ranges, Europa's age range or the exact win rule. With c = 10⁻¹² – 10⁻³ M it got 22.8%/29.4%, so the model is consistent and only the description was incomplete.
- **Fix:** §2.2 now lists every input and the exact win rule.
- **Second blind rebuild, from the paper with the table only:** Enceladus 21.8% (paper 22%), Europa 28.6% (29%), Early Earth 74.6% (75%), Europa with lid 21.9% (22%), Enceladus with lid 0.0% (0%). All within 0.4 points; changing the order of the random draws moves them by at most 0.2. The only missing details were units, now stated.
- **Third pass (2026-10-07), after the §3.15 input table and §2.2 diagnostic definitions were added:** every §3.9 figure and every §3.15 world reproduced from the text within 1 point. Examples: merge shares 33.4/33.6 (paper 33/34), D_z halves 0/43, 0/57, 52/97 (exact), concentration +21/+10, Ganymede 28.5 (22.0–32.7), Mimas 4.7 with 68% pinned. S\* values came from the paper's table and were not recomputed in this pass (see the TEOS-10 row above). The last two ambiguities it found (the volume factor for these worlds, and which depth gives which end of the range) are now stated.

## Not verified (cannot be by re-running code)

- That the model's assumptions are right: mixing rates in these oceans are unmeasured; the chemistry ranges are Europa's for every world.
- That a weak-force-amplifying prebiotic reaction exists (the pre-registered wet-lab test, §5).
- Independent expert review. RNAAS is not peer reviewed.
