"""Independent TEOS-10 check of the S* values in strat_sweep.py / ocean_worlds.py (2026-10-06). Needs gsw: ~/research-env/bin/python."""
import gsw
from scipy.optimize import brentq
def alpha(S, p):  # thermal expansion at the freezing point; p in dbar
    return gsw.alpha(S, gsw.CT_freezing(S, p, 0), p)
def sstar(p_mpa):
    p = p_mpa * 100
    if alpha(0.01, p) > 0: return None  # fresh water already sinks when cooled
    return brentq(alpha, 0.01, 60, args=(p,))
claims = {"Mimas":(1.5,22.6),"Enceladus top":(2,22.5),"Miranda":(2.2,21.9),"Dione":(21,5.6),
          "Umbriel":(21,5.6),"Europa (20 km ice)":(25,3.0),"Ariel":(25,2.7)}
print(gsw.__version__)
for k,(p,c) in claims.items():
    s=sstar(p); print(f"{k:20} {p:5} MPa  paper {c:5}  independent {s:6.2f}")
for p in (26,27,28,29,30,50,200):
    s=sstar(p); print(f"  {p} MPa -> {'none' if s is None else round(s,2)}")
