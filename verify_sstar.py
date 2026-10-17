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
          "Umbriel":(21,5.6),"Europa (20 km ice)":(24.2,3.0),"Europa at ~25 MPa":(25,3.0),"Ariel":(25,2.7)}
# Europa: 917 kg/m3 ice * 1.315 m/s2 * 20 km = 24.2 MPa; strat_sweep.py rounds this to ~25. S* falls steeply here.
# CORRECTION (2026-10-16; the line above is left as written): 24.2 MPa is 920 kg/m3 * 1.315 * 20 km, the density the main
# scripts use; 917 would give 24.1. The 24.2 used below is unchanged.
print(gsw.__version__)
for k,(p,c) in claims.items():
    s=sstar(p); print(f"{k:20} {p:5} MPa  paper {c:5}  independent {s:6.2f}")
for p in (26,27,28,29,30,50,200):
    s=sstar(p); print(f"  {p} MPa -> {'none' if s is None else round(s,2)}")
