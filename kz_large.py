"""H28 (journals/ija_draft.md): a larger 1D Kibble-Zurek system, to test c without the unregistered pooling.

kz_patches.py (A1) found c = 0.876 +- 0.029 at L = 128 xi and 0.882 +- 0.029 at 512 xi; 0.80 (box patch sum) was
excluded by 2.90 SE at 512 xi alone and by 3.95 SE only after pooling the two sizes, which was not pre-registered.
This script runs ONE 1D system at L = 1024 xi (2x the largest so far, 6476 cells, M = L / 12.8 xi ~ 80 patches),
with the same model, D = 4, gamma = 0.01, eps = 1e-4 per cell and dt = 0.05 as kz_patches.sec_1d (it imports
kz_patches.batch / report / check unchanged). Only the Delta design differs: Delta_total = 1.0, 1.5, 2.0 with 1200
runs each (the probit information on c per run peaks near c*Delta ~ 1.6, so this gives SE ~ 0.02 at 3600 runs).
Runs are done in 6 chunks of 200 runs per Delta (seeds 400..405), two chunks at a time in 2 threads.
Fit: P_sim = Phi(c * Delta_total) by maximum likelihood (kz_patches.fit_probit), at lam = 0.4 ("early") and at the
end of the sweep ("end"). c below means c_end at 1024 xi.

PREDICTIONS (written 2026-10-12, before any result-producing run):
 P1. kz_patches.check() (3-cell uniform limit reproduces sim.py Phi(Delta)) passes within 3 SE.
 P2. c(1024 xi) agrees with the earlier plateau 0.879 +- 0.021 (128 and 512 xi) within 2 combined SE: the plateau
     extends from ~10 to ~80 patches.
 P3. This single system alone excludes the box patch-sum value 0.798 by more than 3 SE: (c - 0.798) / SE > 3.
 P4. c(1024 xi) is consistent with the Gaussian-field sign-majority theory 0.914 within 2 SE.
 P5. Walls frozen: the majority at lam = 0.4 equals the majority at the end of the sweep in >= 98% of runs.
 P6. No run reaches a single hand by the end of the sweep (single-hand fraction 0.00 at every Delta).
 P7. Probit shape: at each Delta, P_sim lies within 2 SE of Phi(c * Delta) using the fitted c.
 P8. Pre-registered pool (inverse-variance, c_end at 128, 512 and 1024 xi): excludes 0.798 by more than 4 SE.
 P9. Wall time under 20 minutes (nice -n 10, 2 threads).

Usage: OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 nice -n 10 /opt/anaconda3/bin/python kz_large.py > kz_large.log
"""
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.stats import norm

import kz_patches as K

D, DT, L_XI, DELTAS, CHUNK, NCHUNK = 4.0, 0.05, 1024, (1.0, 1.5, 2.0), 200, 6
C_PREV = {128: (0.876, 0.029), 512: (0.882, 0.029)}      # kz_patches.log, end of sweep
C_PLATEAU, C_BOX, C_THEORY = (0.879, 0.021), K.C_PATCH, K.C_GAUSS[1]


def main():
    t0 = time.time()
    K.check()
    n = int(round(L_XI * K.xi(D)))
    print(f"\n1D periodic, L={L_XI}xi = {n} cells (M={L_XI / K.C_KZ:.1f}), D={D}, gamma={K.GAMMA}, eps={K.EPS}/cell, "
          f"dt={DT}; {NCHUNK} chunks x {CHUNK} runs per Delta, 2 threads", flush=True)
    with ThreadPoolExecutor(2) as ex:
        parts = list(ex.map(lambda i: K.batch(n, 1, D, DELTAS, CHUNK, DT, seed=400 + i), range(NCHUNK)))
    rows = {key: [(d, np.concatenate([p[key][j][1] for p in parts]), np.concatenate([p[key][j][2] for p in parts]))
                  for j, d in enumerate(DELTAS)] for key in ("early", "end")}
    fits = K.report(f"L={L_XI}xi ({n} cells)", rows, ("early", "end"))
    c, se = fits["end"]
    agree = np.mean(np.concatenate([e[1] == f[1] for e, f in zip(rows["early"], rows["end"])]))
    single = [float(o.mean()) for _, _, o in rows["end"]]
    print(f"  early==end majority in {agree:.4f} of runs; single-hand fraction at end {single}")
    print(f"P2 c(1024xi)={c:.3f}+-{se:.3f} vs plateau {C_PLATEAU[0]}+-{C_PLATEAU[1]}: "
          f"{(c - C_PLATEAU[0]) / np.hypot(se, C_PLATEAU[1]):+.2f} combined SE")
    print(f"P3 (c - {C_BOX:.3f}) / SE = {(c - C_BOX) / se:.2f}")
    print(f"P4 (c - {C_THEORY}) / SE = {(c - C_THEORY) / se:+.2f}")
    print(f"P5 early==end {agree:.4f} (need >= 0.98)")
    print(f"P6 single-hand at end max {max(single):.4f}")
    for d, s, _ in rows["end"]:
        p, sp = K.pse(s)
        print(f"P7 Delta={d}: P={p:.4f}+-{sp:.4f} vs Phi(c Delta)={norm.cdf(c * d):.4f}: {(p - norm.cdf(c * d)) / sp:+.2f} SE")
    cs = [C_PREV[128], C_PREV[512], (c, se)]
    w = np.array([1 / s ** 2 for _, s in cs])
    cp, sp = float(np.dot(w, [v for v, _ in cs]) / w.sum()), float(w.sum() ** -0.5)
    print(f"P8 pool 128/512/1024xi: c={cp:.3f}+-{sp:.3f}; (c - {C_BOX:.3f}) / SE = {(cp - C_BOX) / sp:.2f}")
    print(f"P9 total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
