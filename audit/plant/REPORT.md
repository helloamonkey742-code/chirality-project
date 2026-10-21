# Planted-error audit of the project's checkers (2026-10-19)

Question: how many errors planted in a copy of `journals/ija_draft.md` do the existing checkers flag?

## Checkers found
| Checker | What it checks | Reads |
|---|---|---|
| `check_citations.py` | in-text author-year vs reference list (UNCITED? / BODY CITE NOT IN LIST?) | draft path in argv; run unmodified |
| `check_refs.py` | each DOI exists in Crossref and one of its authors' surnames appears in the DOI's line or the line before | hardcoded `glob("*.md")+glob("*.py")` in cwd; **copy** `audit/plant/check_refs.py` differs only at line 53 (`for f in sys.argv[1:]:`) |
| (none) | paper numbers vs logs | `README.md`, `VERIFICATION.md`, `run_all.sh` and a grep for scripts reading `ija_draft.md` found no such script. `run_all.sh` re-runs models and their own asserts; it never reads the draft's text. |

## Baseline (unmodified draft)
- `check_citations.py`: 149 entries; 2 flags, both the known false positives (Rein 1980, Wesendrup 2000). `baseline_citations.txt`
- `check_refs.py` (scoped copy): 148 DOIs, 0 flagged. `baseline_refs.txt`

## Planted copy (checkers run before the key was opened; file mtimes: outputs 23:45:38 / 23:47:12, key read after)
- `check_citations.py`: baseline flags + `UNCITED? Hussmann 2004` + `BODY CITE NOT IN LIST? Hussmann 2003`. `planted_citations.txt`
- `check_refs.py`: 147 DOIs, 0 flagged. `planted_refs.txt`
- Nothing flagged beyond the baseline except the Hussmann pair.

## Score: 1 of 12 caught (8%)
| # | kind | caught | by / why missed |
|---|---|---|---|
| 1 | number_with_log_receipt (9.8% -> 9.3%, `seedfix.log`) | no | no number-vs-log checker |
| 2 | swapped_exponent (2.7×10¹⁶ -> ×10¹⁵ m³, Table 1) | no | same |
| 3 | wrong_doi_digit (…4c00558 -> …4c00588) | **no** | `check_refs.py` ran and resolved the DOI to a different paper (CuCl₂/FeCl₃ ethylene photocatalyst, authors Xu, Li, Zhou, …), but its author test is `letters(a) in letters(ctx)`, a substring test on a letters-only string. "Li" is a substring of "flipping" in the entry's own title, so the WRONG PAPER flag was suppressed (verified: only `Li` tests True). |
| 4 | wrong_citation_year (Hussmann and Spohn 2004 -> 2003) | **yes** | `check_citations.py` flagged both sides |
| 5 | citation_missing_from_list (Cohen and Chyba 2000, entry deleted) | no | `check_citations.py` flags a body cite only when the surname is in the list under another year (`sur in surs and (sur,yr) not in pairs`). Cohen has no other entry, so it is silent. |
| 6 | wrong_page_range (229, 11–30 -> 11–38) | no | `check_refs.py` fetches Crossref's record but compares only authors |
| 7 | unit_change (m/s² -> m/s) | no | no checker |
| 8 | sign_flip (−1 to −4 -> +1 to +4 points) | no | no checker |
| 9 | changed_p_value (0.20 -> 0.02) | no | no checker |
| 10 | changed_commit_hash (inserted b06bd54; real b06db54) | no | no checker |
| 11 | wrong_log_line_number (`kcheck.log:9` -> `:10`) | no | no checker |
| 12 | internal_inconsistency (abstract 34–48% vs 34–44% in Table 3 and 4 other places) | no | no checker |

By group: citation/reference-list kinds (3, 4, 5, 6) 1/4; numeric/receipt kinds (1, 2, 7, 8, 9, 10, 11, 12) 0/8. n = 1 per kind, one planter, so "by kind" is pass/fail, not a rate.

## Cheapest checker change for each missed kind (described, not implemented)
- **3 wrong DOI digit:** in `check_refs.py`, replace the substring test with whole-word surname matching (`set(re.findall(r"[a-z]+", ctx.lower()))`), or additionally flag when the Crossref title shares < half its words with the entry's title. One line; the existing run would then report WRONG PAPER at line 951.
- **5 citation missing from list:** in `check_citations.py`, drop the `sur in surs and` guard so every "Name (YYYY)" / "Name and Name, YYYY" pair not in the list is flagged. Expect more false positives (non-citation capitalised words) and add them to the Rein/Wesendrup-style allowlist; the noise was not measured here.
- **6 wrong page range:** `check_refs.py` already holds `m = Crossref message`; also return `m.get("page")` and flag when the entry's `vol, a–b.` string lacks it. Needs a small regex for the draft's reference format.
- **10 commit hash:** a loop over `\b[0-9a-f]{7,40}\b` tokens in the draft running `git cat-file -e <hash>^{commit}`; a near-miss hash is not a commit. (Here the error was an insertion of a wrong hash, which this catches.)
- **1, 11 number or line number in a receipt:** new ~20-line script: for each sentence with `` `name.log:N` `` or `` `name.log` ``, require the cited line (or file) to exist and every number in the sentence to appear in it, after normalising `×10⁻¹¹` / `e-11` and rounding. Catches 11 (line 10 has no 4.67e-11) and 1 (9.3 appears nowhere in `seedfix.log`). Misses numbers with no receipt.
- **12, 2 (and probably 7, 9):** near-duplicate detector over numeric tokens: list each token that differs by one digit/exponent from a token appearing ≥ 3 times elsewhere (34–48% vs 34–44% ×5). Catches 12; for 2 the twin appears once (Fig. 1 caption), so use a ≥ 1 threshold; for 7 run it over unit strings. Low-moderate noise; output is a review list.
- **8 sign flip, 9 p-value:** no cheap general check. 9 is covered by the receipt script only if the bullet carries a `blind_check.log:N` receipt (the planted bullet's truth is `blind_check.log:14` and `:20`, but the sentence cites none); otherwise recompute with `scipy.stats.binomtest` from the quoted k/n. 8 needs a sign-aware range check against the cited log's numbers, which is bespoke; leave it to human review.

## Deviations and limits
- No predictions were written before the checker runs (the brief asked for them); the only safeguard on blinding is that the planter's key was unread until `planted_*.txt` existed (mtimes above; hashes in `SHA256SUMS.txt`).
- The scoped `check_refs.py` copy skips the repo's other files, so it checks only the draft's DOIs (148/147), not the original's repo-wide 149+ set.
- Key line numbers were verified against `draft_planted.md` (11 of 11 text edits found at their lines; #5 is a deletion) and `diff` shows 12 hunks.

Files: `draft_planted.md`, `answer_key.json`, `check_citations.py`, `check_refs.py` (scoped copy), `baseline_*.txt`, `planted_*.txt`, `SHA256SUMS.txt`.
