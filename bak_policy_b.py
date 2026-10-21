"""POST HOC.  Bak et al. (2023) Fig. 2a under the D2 rule now in journals/droplet_prereg.md (policy B of
misread_sim.py: both fits use the misread rate, c* = -2.21).  Bak had no control drops, so the rate is unknown and
is tried at 0-10%.  Bak's data were read long before this was written; this is a re-fit, not a test.
n at the P = 1 markers is not printed; tried at 3, 6, 10 as in exclusion_bak.py.
sharp = 2 (ll best - ll best STEP with s <= 0.15); > 3.84 rejects a sharp step.

Run: UV_CACHE_DIR=~/Library/Caches/uv uv run --no-project --with numpy --with scipy --with numba python bak_policy_b.py > bak_policy_b.log
"""
import numpy as np

from blind_units import BAK_MAIN, BAK_TINC90
from droplet_power import scorable
from droplet_sim2 import EPS_FIT, classify, curves, stats
from exclusion_bak import bak

C_B = -2.21                                   # misread_sim.log, c*_B(128), frozen
EPS = (0.0, 0.02, 0.04, 0.06, 0.08, 0.10)


def with_eps(cur, e):
    """curves() with the misread blend e in place of EPS_FIT."""
    out = {}
    for fam, (P, _, _, par) in cur.items():
        obs = np.clip(e + (1 - 2 * e) * P, 1e-12, 1 - 1e-12)
        out[fam] = (P, np.log(obs).T, np.log1p(-obs).T, par)
    return out


def best(K, N, cur, fam, keep):
    _, lp, l1, par = cur[fam]
    return (K @ lp + (N - K) @ l1)[:, keep(par)].max(1)[0]


def self_check():
    x, _, _ = bak(BAK_TINC90, 3)
    c, d = curves(x), with_eps(curves(x), EPS_FIT)
    assert all(np.allclose(c[f][i], d[f][i]) for f in c for i in (1, 2))


def main():
    np.seterr(over="ignore", divide="ignore", invalid="ignore")
    self_check()
    print(f"POST HOC.  Re-fit of published points, not new data.  Frozen c* = {C_B}; misread rate unknown, tried 0-10%.")
    for name, series in (("Tinc 90 C", BAK_TINC90), ("main", BAK_MAIN)):
        for n1 in (3, 6, 10):
            x, k, n = bak(series, n1)
            K, N = k[None], n[None]
            row = []
            for e in EPS:
                cur = with_eps(curves(x), e)
                r = stats(K, N, x, cur)
                r["ns"] = np.array([not scorable(k, n)])
                lab = classify(r, C_B)[0]
                top = max(best(K, N, cur, "count", lambda p: p > 0), best(K, N, cur, "step", lambda p: p > 0))
                sharp = 2 * (top - best(K, N, cur, "step", lambda p: p <= 0.15))
                row.append(f"e {e:.2f}: L {r['L'][0]:+5.2f} sharp {sharp:5.1f} {lab}")
            print(f"{name:9} n={n1:<2} | " + " | ".join(row))


if __name__ == "__main__":
    main()
