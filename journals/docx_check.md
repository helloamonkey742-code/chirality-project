# ija_manuscript.docx: conversion check (2026-10-13)

**Output:** `journals/ija_manuscript.docx` (1.3 MB, untracked).
**Source:** `journals/ija_draft.md` (commit ff06793, sha256 11f2e7030b13…, not edited), from the line `# How large a mixed volume…` (line 252) to the end. The submission notes (lines 1–250) are not included.

## How it was built

Script: `/private/tmp/claude-501/…/scratchpad/docx/build_docx.py`, run with `/opt/anaconda3/bin/python build_docx.py`. It does the following:

1. Copies the manuscript part to `scratchpad/docx/manuscript.md` and changes only that copy:
   - `\tag{n}` becomes `\qquad\text{(n)}`, because texmath cannot handle `\tag`. The number stays inside the Word equation as upright "(n)".
   - It adds a hard line break after the author name.
2. Builds a reference document from pandoc's default. Fonts are Times New Roman 12 pt, body text is double-spaced, and headings are black. Table cells stay single-spaced.
3. Runs pandoc 2.12:
   ```
   cd journals && pandoc manuscript.md -f markdown-implicit_figures -t docx \
     --resource-path journals --reference-doc reference.docx -o raw.docx
   ```
   pandoc gave no warnings. `-implicit_figures` is used because the bold **Figure N.** paragraph under each image is already the caption. Without it, every figure would get a second caption taken from the alt text. The alt text is kept as the image description, which helps meet IJA's accessibility rules.
4. Edits the final section settings in `document.xml`:
   - US Letter page with 1-inch margins
   - continuous line numbers (`w:lnNumType`)
   - a centred page-number footer

The checks are in `scratchpad/docx/check_docx.py`. Every hard check passes.

## Counts (source → docx)

| Item | Source | docx |
|---|---|---|
| Image links / `word/media` files / `a:blip` refs | 6 | 6 / 6 / 6 |
| Display equations (`$$`) / `m:oMathPara` / `m:oMath` | 23 | 23 / 23 / 23 |
| Equation numbers | 23 | 23, each once, in this order: 1 2a 2b 2c 3 4 … 21 |
| Tables (pipe) / `w:tbl` | 5 | 5 |
| Raw TeX commands in body text or inside OMML | – | none |
| `[ADD …]` placeholders | 9 | 9 (expected) |
| `[CHECK FULL TEXT…]` markers | 9 | 9 (these are in the manuscript body, not leaked from the notes) |
| `[CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]` | 1 | 1 (Acknowledgements) |
| Note text (`END OF SUBMISSION NOTES`, `To do before submission`, `UNCONFIRMED`, E-numbers) | 0 | 0 |

- **Images:** all six are 5.83 in wide, which fits the 6.5 in text width.
- **XML:** every part of the package parses.

**Equations that failed to convert:** none. I converted the docx back to markdown with pandoc, and all 23 equations came back as TeX that matches the source. The only cosmetic change is that `H'` becomes `H\prime`.

Eqs. 2a–2c are three separate Word equations, numbered (2a), (2b) and (2c). This answers submission note H9.

**Visual check:** limited. LibreOffice (`soffice`) is not installed, and neither are Word or Pages, so I could not render a PDF. A Quick Look thumbnail of page 1 shows the title block and the abstract correctly. Quick Look does not show line numbers, double spacing or equations, so I did not look at the equations on screen. Open the file in Word once.

## IJA requirements (source: [Preparing your materials](https://www.cambridge.org/core/journals/international-journal-of-astrobiology/information/author-instructions/preparing-your-materials))

- **Format:** Word or LaTeX. There is no Cambridge Word template; the template linked from that page is LaTeX only. The Word rule is literal and old: PC files must be named `manuscript.doc` and Mac files saved as `manuscript.rtf`. I kept .docx because ScholarOne normally accepts it. Rename the file, or save it as .doc/.rtf, if the upload rejects it. Cambridge sets equations made in Word again during production.
- **Figures:** embed them for review (done). On acceptance, supply separate EPS or TIFF files: 1200 dpi for line art, 300 dpi for colour, 600 dpi for mixed. Name them by figure number.
- **Equations:** numbered (1), (2), … and set in Equation Editor, not as images. Done, using Word's native OMML.
- **References:** Harvard (author, date), alphabetical, with article titles.
- **Abstract:** at most 300 words, one paragraph. The notes give 298.
- **Title page:** title, authors, full institutional address, and the corresponding author's address, phone and e-mail.
- **Competing interests:** a statement is required. The manuscript has one.
- **AI use:** declare it with tool name, version and dates. The placeholders are still open.
- **Line numbers and double spacing:** not required by IJA. I added both.
- **Keywords:** the page contradicts itself. It asks for 5–10 keywords on the title page, but also says to pick 2 (+1) during online submission instead. The manuscript has none.
- **Submission:** through ScholarOne, at `mc.manuscriptcentral.com/submit_ija`.

## What the author must still do

1. Fill the 9 `[ADD …]` fields:
   - title page: address ×2, telephone, e-mail, ORCID
   - Sections 2.8 and Acknowledgements: AI model version(s) ×2 and date range ×2
2. Resolve the 9 `[CHECK FULL TEXT]` markers (8 in Section 1.1 and 1 in Section 1.2; note H7). Delete the explanatory sentence "Items marked [CHECK FULL TEXT]…" in Section 1.1.
3. Get Dr. McKay's agreement and delete `[CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]`. Ask Prof. Nimmo about the "pers. comm." (H31a).
4. Decide on keywords: either add 5–10 under the abstract, or pick them only in ScholarOne.
5. Number the reaction-network table in Section 2.2.1 or leave it inline (H10). The code-availability table is also unnumbered.
6. Re-render the `coarsening.png` legend for Figure 4 (H11). The embedded image still shows "l = 5.5·√(D t)".
7. Open the file in Word once and check the equations, figure placement and table borders. Pandoc's default table style has light borders.
8. Rename the file or save it as .doc/.rtf if ScholarOne rejects .docx.
9. Prepare EPS/TIFF figure files for the final version.

To rebuild after editing `ija_draft.md`, run `/opt/anaconda3/bin/python3 -I journals/build_manuscript.py`, then `journals/check_manuscript.py`, from the repo root. They replaced the scratchpad `build_docx.py` and `check_docx.py` (2026-10-20), and also build `ija_manuscript.pdf` and `ija_draft.pdf` (pandoc to HTML with `manuscript.css`, then headless Chrome).
