"""Check every DOI cited in the project against Crossref. (2026-10-10)

Asking an AI "is this reference real?" checks nothing; this asks the registry. For each DOI in the .md and .py files it
reports NOT FOUND (Crossref has no such DOI) or WRONG PAPER (none of the paper's authors is named next to the DOI, so
the DOI most likely points at a different paper). Needs internet, so it is not in run_all.sh: run it before sending
the paper or the review to anyone.
    /opt/anaconda3/bin/python check_refs.py            # exit 1 if anything is flagged
A flag is a lead, not a verdict: open https://doi.org/<doi> and look.
"""
import glob
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

DOI = re.compile(r"10\.\d{4,9}/[^\s\"'`,\]]+")


def clean(d):
    """Strip markdown escapes and trailing punctuation; keep parentheses that belong to the DOI."""
    d = d.replace("\\<", "<").replace("\\>", ">").rstrip(".:;*>")
    while d.endswith(")") and d.count(")") > d.count("("):
        d = d[:-1].rstrip(".:;")
    return d


def letters(s):
    return re.sub(r"[^a-z]", "", unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower())


def crossref(d):
    """(title, [author surnames]) or None if Crossref has no such DOI."""
    url = "https://api.crossref.org/works/" + urllib.parse.quote(d, safe="/()<>:;")
    for wait in (1, 4, 16):
        try:
            m = json.load(urllib.request.urlopen(url, timeout=30))["message"]
            return (m.get("title") or [""])[0], [a.get("family", "") for a in m.get("author", [])]
        except urllib.error.HTTPError as e:
            if e.code == 404:  # not Crossref's; arXiv/DataCite DOIs are only in the DOI handle registry
                h = json.load(urllib.request.urlopen("https://doi.org/api/handles/" + urllib.parse.quote(d), timeout=30))
                return ("(registered outside Crossref)", []) if h.get("responseCode") == 1 else None
            time.sleep(wait)  # 429 / 5xx: back off and retry
    raise RuntimeError(f"Crossref kept failing for {d}")


def cites():
    """{doi: (file, line number, the citing line and the one before)} for the first place each DOI appears."""
    out = {}
    for f in sys.argv[1:]:
        lines = open(f, encoding="utf-8").read().splitlines()
        for i, line in enumerate(lines):
            for d in map(clean, DOI.findall(line)):
                out.setdefault(d, (f, i + 1, (lines[i - 1] if i else "") + " " + line))
    return out


def main():
    first = lambda text: clean(DOI.findall(text)[0])
    assert first("doi:10.1016/0012-8252(76)90011-8).") == "10.1016/0012-8252(76)90011-8"
    assert first("(10.1021/ja7106349):") == "10.1021/ja7106349"
    assert first("10.1175/1520-0485(1980)010\\<0083:EOTLRO\\>2.0.CO;2;") == "10.1175/1520-0485(1980)010<0083:EOTLRO>2.0.CO;2"
    assert letters("Castillo‐Rogez") == letters("Castillo-Rogez") == "castillorogez"

    refs, bad = cites(), 0
    for d, (f, n, ctx) in refs.items():
        time.sleep(0.3)
        hit = crossref(d)
        if hit is None:
            print(f"NOT FOUND   {f}:{n}  {d}")
        elif hit[1] and not any(letters(a) and letters(a) in letters(ctx) for a in hit[1]):
            print(f"WRONG PAPER {f}:{n}  {d}  -> really: {hit[0][:80]} ({', '.join(hit[1][:2])})")
        else:
            continue
        bad += 1
    print(f"{len(refs)} DOIs checked, {bad} flagged")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
