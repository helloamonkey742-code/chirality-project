# Silva-Dias & Weber (2026) against Eq. 10–13

Paper: L. Silva-Dias, C.A. Weber, "How fluctuations control chiral symmetry breaking in reacting mixtures of chiral molecules", arXiv:2610.13523 (v1, 2026-10-09), doi:10.48550/arXiv.2610.13523. Read in full (16 pp.) on 2026-10-17. The PDF was kept in a scratch directory outside the repo and is not committed.

**What they do.** They study the Kondepudi–Nelson network with production δ, decay γ = k2/k1, autocatalysis σ, a quadratic loss ω and antagonism ϵ. Their route is: chemical master equation → Kramers–Moyal Fokker–Planck (Eqs. 13–14) → adiabatic elimination of the total density ϕ (Eq. 15) → a 1D Fokker–Planck in the imbalance ψ → splitting probability (Eq. 28), mean first-passage time and Kramers rate. The parameters are held fixed inside the bistable regime (σ > γ). The only bias is the initial excess ψ0. Gillespie runs (50,000 realizations) check the analysis.

## Q1. Does it reduce to Φ(a·ε·√N), or to Eq. 11?

It reduces to Φ(a·ε·√N), not to Eq. 11.

- Eq. 28 is the exact splitting probability for the reduced 1D process.
- Eq. 30 is its near-racemic form, Π ≈ 1/2 + √κ ψ0 / (√π erf(ψ+√κ)), with κ = a′(0)/b(0) = cN.
- Eq. 32 is the large-N limit, Π ≈ 1/2 + √(cN/π) ψ0. They read the control parameter as ψ0/δψ with δψ ∝ N^(−1/2).
- They never write Φ. For a linear unstable drift with constant diffusion, the full result is ½ + ½ erf(ψ0√κ) = Φ(√(2c) · ψ0 · √N), and Eq. 32 is its first-order term. The bias is ε = ψ0, the initial excess. The prefactor is a = √(2c) = √(2a′(0)/(N b(0))), a ratio of network constants with no universal number in it.
- Nothing matches Eq. 11's π^(1/4), its rate bias g or its (kτ)^(1/4).
- Their N is a system size: x = X/N is a density. Our N is a molecule count.

The draft's Eq. 8 gives their limit directly. Take a sudden quench (λ constant, initial excess α0, g = 0, noise ε): then Δ = α0 √(2λ/ε). Setting λ = a′(0) and ε = b(0) recovers their Eq. 32. So their result is the quench, initial-bias limit of the same Gaussian bookkeeping. Eq. 10–11 is the slow-sweep, persistent-bias case. At their Fig. 7 parameters, c = 2.3×10⁻⁴ per unit N.

## Q2. Sweep rate (kτ), or only an initial bias?

Only an initial bias. There is no ramp: σ, γ and δ are constant in time and past the bifurcation from t = 0, and the rates carry no systematic bias. Their γ is a decay ratio. It is not our sweep rate γ = 1/(kτ), so keep the two apart if both papers are discussed together. Time enters only through the mean first-passage time and the Kramers rate, both of which grow with N. Their conclusion (ii) says selection has to happen within the lifetime of the nonequilibrium conditions, a statement about duration, not about sweep rate. The (kτ)^(1/4) in Eq. 11 and the √(kτ) in Eq. 13 have no counterpart in their paper, so that part of the draft gets no support from them and no contradiction either.

## Q3. Anything against frank.py's F ≈ 0.7 (README Part 7)?

No. Their ψ diffusion (Eq. 14) is

β_ψψ = N⁻¹ (2δ + (σ+γ)ϕ + (ω/2)(ϕ² + ψ²)),

which counts production, decay and autocatalysis together, as draft Eq. 5 and frank.py do. Antagonism (X + Y → ·) leaves ψ unchanged and drops out. The amplifying step contributes σϕ. The other steps add noise and no amplification, which is the reason F < 1.

They do not define F. The analog in their model is F² = σϕ* / (2δ + (σ+γ)ϕ* + ωϕ*²/2) at ψ = 0. At the Fig. 7 parameters (ω = 0.01, γ = 0.9, σ = 1, ϵ = 0.1, δ = 2), ϕ* = 9.49 from Eq. A41, and F = 0.65. As δ → 0 and γ → σ (copying ≈ outflow) this tends to 1/√2 = 0.71, which is also the limit of frank.py's √(x/(k0+2x)). The two have the same structure, and a value of 0.65 against 0.70 reflects their larger production term. This is my calculation, not a number from the paper.

## Draft change

- One sentence after Eq. 12 in Section 2.2.4, citing them as the quench, initial-excess counterpart.
- A reference entry.
- Log entry E103.

No number in the draft changes.
