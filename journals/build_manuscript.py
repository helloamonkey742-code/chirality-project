"""Rebuild journals/ija_manuscript.docx, ija_manuscript.pdf and ija_draft.pdf from journals/ija_draft.md.

Recipe in docx_check.md.  ija_draft.md is read only; edited copies and intermediates go to a temp dir.
The manuscript files start at the title line; ija_draft.pdf keeps the submission notes.
Run from the repo root: /opt/anaconda3/bin/python3 -I journals/build_manuscript.py, then journals/check_manuscript.py
Needs pandoc 2.12 (/opt/anaconda3/bin/pandoc) and Google Chrome.
"""
import io
import re
import subprocess
import tempfile
import zipfile
from pathlib import Path

J = Path(__file__).resolve().parent
PANDOC = "/opt/anaconda3/bin/pandoc"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
TITLE = "# How large a mixed volume"


def prep(s):
    """\\tag{n} -> \\qquad\\text{(n)} (texmath cannot do \\tag); hard break after the author name."""
    s, n = re.subn(r"\\tag\{([^}]*)\}", r"\\qquad\\text{(\1)}", s)
    assert n == 23, n
    a = "**Saptarshi Ghosh**\nCupertino"
    assert s.count(a) == 1
    return s.replace(a, "**Saptarshi Ghosh**\\\nCupertino")


def run(*a):
    r = subprocess.run(a, capture_output=True)
    err = "\n".join(l for l in (r.stdout + r.stderr).decode(errors="replace").splitlines()
                    if "deprecated" not in l and "bytes written" not in l)
    assert r.returncode == 0, err
    return r.stdout, err


def reference_docx(path):
    """pandoc's default reference.docx: Times New Roman 12 pt, double-spaced body, black headings,
    single-spaced table cells."""
    zi = zipfile.ZipFile(io.BytesIO(run(PANDOC, "--print-default-data-file", "reference.docx")[0]))
    st = zi.read("word/styles.xml").decode()
    tnr = '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" w:cs="Times New Roman" />'
    st = re.sub(r"<w:rFonts [^>]*Theme[^>]*/>", tnr, st)
    st = st.replace('<w:spacing w:after="200" />', '<w:spacing w:after="200" w:line="480" w:lineRule="auto" />', 1)
    st = st.replace('<w:spacing w:before="36" w:after="36" />',
                    '<w:spacing w:before="36" w:after="36" w:line="240" w:lineRule="auto" />')
    st = re.sub(r'(<w:style [^>]*w:styleId="(?:Heading\d|Title|Subtitle)".*?</w:style>)',
                lambda m: re.sub(r'<w:color w:val="[0-9A-Fa-f]{6}" />', '<w:color w:val="000000" />', m.group(1)),
                st, flags=re.S)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zo:
        for it in zi.infolist():
            zo.writestr(it, st if it.filename == "word/styles.xml" else zi.read(it))


FOOTER = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
          '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
          '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:fldChar w:fldCharType="begin"/></w:r>'
          '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r>'
          '<w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>')
# US Letter, 1-inch margins, continuous line numbers, centred page-number footer
SECT = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdPgNum"/><w:pgSz w:w="12240" w:h="15840"/>'
        '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>'
        '<w:lnNumType w:countBy="1" w:restart="continuous"/></w:sectPr>')


def finish_docx(raw, out):
    with zipfile.ZipFile(raw) as zi, zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zo:
        for it in zi.infolist():
            d = zi.read(it.filename)
            if it.filename == "word/document.xml":
                t = d.decode()
                k = t.rindex("<w:sectPr")
                assert 'xmlns:r=' in t
                d = (t[:k] + SECT + t[t.index("</w:body>", k):]).encode()
            elif it.filename == "word/_rels/document.xml.rels":
                d = d.decode().replace("</Relationships>",
                    '<Relationship Id="rIdPgNum" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
                    'Target="footer1.xml"/></Relationships>').encode()
            elif it.filename == "[Content_Types].xml":
                d = d.decode().replace("</Types>",
                    '<Override PartName="/word/footer1.xml" '
                    'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>').encode()
            zo.writestr(it, d)
        zo.writestr("word/footer1.xml", FOOTER)


def main():
    src = (J / "ija_draft.md").read_text()
    with tempfile.TemporaryDirectory() as tmp:
        T = Path(tmp)
        (T / "full.md").write_text(prep(src))
        (T / "manu.md").write_text(prep(src[src.index("\n" + TITLE) + 1:]))

        # -implicit_figures: the bold "Figure N." paragraph under each image is the caption; the alt text stays as
        # the image description instead of becoming a second caption.
        reference_docx(T / "reference.docx")
        print("pandoc docx warnings:", run(PANDOC, str(T / "manu.md"), "-f", "markdown-implicit_figures", "-t", "docx",
                                           "--resource-path", str(J), "--reference-doc", str(T / "reference.docx"),
                                           "-o", str(T / "raw.docx"))[1] or "none")
        finish_docx(T / "raw.docx", J / "ija_manuscript.docx")

        for md, name in (("full.md", "ija_draft"), ("manu.md", "ija_manuscript")):
            html = T / f"{name}.html"
            print(f"pandoc {name}.html warnings:", run(
                PANDOC, str(T / md), "-f", "markdown-implicit_figures", "-t", "html5", "-s", "--mathml", "--self-contained",
                "--metadata", f"pagetitle={name}", "--resource-path", str(J), "-c", str(J / "manuscript.css"),
                "-o", str(html))[1] or "none")
            run(CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer", f"--user-data-dir={T / 'chrome'}",
                f"--print-to-pdf={J / (name + '.pdf')}", html.as_uri())
    print("built: ija_manuscript.docx ija_manuscript.pdf ija_draft.pdf")


if __name__ == "__main__":
    main()
