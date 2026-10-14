"""Cross-check in-text citations against the reference list of a paper draft (2026-10-13).
Prints UNCITED? for a list entry never cited in the text, and BODY CITE NOT IN LIST? for a text citation with no entry.
A flag is a lead: three-author citations (e.g. Hegstrom, Rein and Sandars 1980) flag as a false positive.
    /opt/anaconda3/bin/python check_citations.py journals/ija_draft.md
"""
import re,sys
s=open(sys.argv[1]).read()
body=s[s.index("\n# How large"):s.index("\n## References")]
refs=s[s.index("\n## References"):]
entries=[l for l in refs.split("\n") if re.match(r"^[A-ZÀ-ſ][^ ]*, ",l)]
def key(e):
    sur=e.split(",")[0].strip()
    m=re.search(r" (\d{4}[ab]?),",e)
    return sur,m.group(1) if m else None
keys=[key(e) for e in entries]
print("entries",len(entries))
# entries cited?
for sur,yr in keys:
    y=yr
    pat=re.escape(sur)+r"[^()]{0,80}?\(?"+r"[^)]{0,40}?"+(y[:4])
    found=False
    for m in re.finditer(re.escape(sur),body):
        seg=body[m.start():m.start()+160]
        if y[:4] in seg and (len(y)==4 or y in seg or (y[:4]+", "+y[:4]) in seg or re.search(y[:4]+r"[ab]?, "+y[:4]+y[4:],seg) or re.search(y[:4]+"a, "+y[4:],seg)):
            found=True;break
    if not found: print("UNCITED?",sur,yr)
# body citations in list?
surs={k[0] for k in keys}
pairs={(k[0],k[1][:4]) for k in keys}
for m in re.finditer(r"([A-ZÀ-Ž][a-zà-ž'’\-]+(?: [A-Z][a-z]+)?)(?: et al\.| and [A-ZÀ-Ž][\w\-]+| [A-Z]?[a-z]*)?,? \(?((?:1[89]|20)\d{2})[ab]?\b",body):
    sur,yr=m.group(1),m.group(2)
    if sur in surs and (sur,yr) not in pairs: print("BODY CITE NOT IN LIST?",sur,yr,body[m.start():m.start()+50].replace("\n"," "))
