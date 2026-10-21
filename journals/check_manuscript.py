"""Check ija_manuscript.docx / .pdf and ija_draft.pdf (from build_manuscript.py) against journals/ija_draft.md.

Run from the repo root: /opt/anaconda3/bin/python3 -I journals/check_manuscript.py   (exit 1 on any FAIL)
"""
import re
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.dom import minidom

J = Path(__file__).resolve().parent
src = (J / "ija_draft.md").read_text()
manu = src[src.index("\n# How large a mixed volume") + 1:]
fails = []


def ok(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'} {name} {detail}")
    if not cond:
        fails.append(name)


def norm(s):
    """Whitespace, soft hyphens and pandoc's smart quotes."""
    s = s.replace("\u00ad", "").translate(str.maketrans("\u201c\u201d\u2018\u2019", "\"\"''"))
    return re.sub(r"\s+", " ", s).strip()


def headings(md):
    return [norm(re.sub(r"[*`]", "", h)) for h in re.findall(r"^#{1,4} (.+)$", md, re.M)]


tags = re.findall(r"\\tag\{([^}]*)\}", manu)
heads = headings(manu)
marks = {m: manu.count(m) for m in ("[ADD ", "[CHECK FULL TEXT", "[CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]")}
leaks = ("END OF SUBMISSION NOTES", "To do before submission", "UNCONFIRMED", "SUBMISSION NOTES")
must = ("Jang et al., 2024", "blind_jang.py", "5055–5064", "Five of seven predictions held", "p = 0.048")
n_img = manu.count("![")
n_tbl = len(re.findall(r"^\|[ :|-]+\|\s*$", manu, re.M))
print(f"source manuscript: {len(tags)} tags, {n_img} images, {n_tbl} tables, {len(heads)} headings, marks {marks}")

print("\n== ija_manuscript.docx")
z = zipfile.ZipFile(J / "ija_manuscript.docx")
bad = []
for n in z.namelist():
    if n.endswith((".xml", ".rels")):
        try:
            minidom.parseString(z.read(n))
        except Exception as e:
            bad.append(f"{n}: {e}")
ok("every XML part parses", not bad, "; ".join(bad))
doc = z.read("word/document.xml").decode()
ok("images", len([n for n in z.namelist() if n.startswith("word/media/")]) == n_img == doc.count("<a:blip "), f"{doc.count('<a:blip ')}")
ok("display equations", doc.count("<m:oMathPara>") == len(tags), f"{doc.count('<m:oMathPara>')} oMathPara, {doc.count('<m:oMath>')} oMath")
paras = re.findall(r"<m:oMathPara>.*?</m:oMathPara>", doc, re.S)
nums = [re.findall(r"\(([0-9]+[a-z]?)\)", "".join(re.findall(r"<m:t[^>]*>([^<]*)</m:t>", p)))[-1:] for p in paras]
ok("equation numbers in order", [x[0] if x else None for x in nums] == tags, " ".join(x[0] if x else "?" for x in nums))
ok("tables", doc.count("<w:tbl>") == n_tbl, f"{doc.count('<w:tbl>')}")
text = "".join(re.findall(r"<(?:w|m):t(?: [^>]*)?>([^<]*)</(?:w|m):t>", doc)).replace("&amp;", "&")
raw = re.findall(r"\\[A-Za-z]{2,}", text)
ok("no raw TeX", not raw, str(sorted(set(raw))[:10]))
for m, c in marks.items():
    ok(f"marker {m!r}", text.count(m) == c, f"{text.count(m)} of {c}")
for l in leaks:
    ok(f"no notes text {l!r}", l not in text)
for m in must:
    ok(f"has {m!r}", m in text)
missing = [h for h in heads if h not in norm(text)]
ok("every heading present", not missing, str(missing[:5]))
sect = doc[doc.rindex("<w:sectPr"):]
ok("letter, 1in margins, line numbers, page footer",
   all(s in sect for s in ('w:w="12240"', 'w:top="1440"', "<w:lnNumType", "footerReference")) and "word/footer1.xml" in z.namelist())
st = z.read("word/styles.xml").decode()
ok("Times New Roman, double spacing", "Times New Roman" in st and 'w:line="480"' in st and "Theme=" not in re.sub(r"<w:latentStyles.*?</w:latentStyles>", "", st, flags=re.S))

for name, part in (("ija_manuscript.pdf", manu), ("ija_draft.pdf", src)):
    print(f"\n== {name}")
    p = J / name
    info = subprocess.run(["pdfinfo", str(p)], capture_output=True, text=True).stdout
    print("  pages " + re.search(r"Pages:\s+(\d+)", info).group(1) + ", " + re.search(r"Page size:\s+(.+)", info).group(1))
    fonts = subprocess.run(["pdffonts", str(p)], capture_output=True, text=True).stdout.splitlines()[2:]
    ok("all fonts embedded", all(" yes " in f[50:] for f in fonts), f"{len(fonts)} fonts")
    t = subprocess.run(["pdftotext", "-enc", "UTF-8", str(p), "-"], capture_output=True, text=True).stdout
    nt = norm(t).replace(" ", "")
    raw = sorted(set(re.findall(r"\\(?:qquad|text|tag|frac|mathrm|left|right|begin)\b|\$\$", t)))
    ok("no raw TeX", not raw, str(raw))
    imgs = len(subprocess.run(["pdfimages", "-list", str(p)], capture_output=True, text=True).stdout.splitlines()[2:])
    ok("images", imgs >= part.count("!["), f"{imgs} image objects for {part.count('![')} figures")
    ok("equation numbers", all(f"({x})" in t for x in tags), f"{sum(f'({x})' in t for x in tags)} of {len(tags)}")
    ph = headings(part)
    missing = [h for h in ph if h.replace(" ", "") not in nt]
    ok("every heading present", not missing, f"{len(ph) - len(missing)} of {len(ph)} {missing[:3]}")
    for m in must:
        ok(f"has {m!r}", m.replace(" ", "") in nt)
    ok("notes " + ("included" if part is src else "excluded"), ("SUBMISSION NOTES" in t) == (part is src))
    print(f"  doi.org links in text: {nt.count('doi.org/')} (source {part.count('doi.org/')})")

print("\nFAILED:", fails if fails else "none")
sys.exit(1 if fails else 0)
