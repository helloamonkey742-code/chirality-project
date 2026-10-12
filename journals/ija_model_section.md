# 2 The model

<!-- Replacement for the opening of "2 Methods" in ija_draft.md (Sections 2.1-2.4 there). Written 2026-10-11.
     Every number below is taken from ija_draft.md, PAPER_DRAFT.md, VERIFICATION.md, a script or its log.
     Each equation names the script and function that implements it. Conflicts between sources are listed
     at the end under "Source conflicts" and are not resolved here. -->

This section describes the physical effect I model, step by step: where the bias comes from and how large it is (2.1); the chemistry that amplifies it and the equation that chemistry reduces to (2.2); how many molecules take part and why mixing sets that number (2.3); what happens after the transition, when regions of opposite hand meet (2.4); the test for a stratified layer under the ice (2.5); how the uncertain inputs become the reported percentages (2.6); and which checks support each step (2.7). Figure 1 summarises the picture.

![Figure 1. Schematic of the model.](fig_schematic.png)

**Figure 1.** Schematic of the model. Nothing in this figure is simulated or measured data. (a) Fixed points of the noise-free normal form, λα − α^3 + g = 0 (Eq. 6), against the control parameter λ. Solid lines are stable states and dashed lines unstable ones. The dotted grey lines are the symmetric pitchfork (g = 0). The bias is drawn as g = 0.05 so that it is visible; the real bias is about 10^-17. Before the transition (λ < 0) the only stable state is offset from racemic by about g/|λ|. After it (λ > 0) the branch connected to that state leads to the PVED-favoured hand. The grey band marks, schematically, the window |λ| ≲ γ^(1/2) in which the outcome is set (Eq. 9). (b) Cartoon of an ice-covered ocean, not to scale: ice shell, a stratified lid that forms when the salinity S is below the threshold S* (Section 2.5), and a convecting interior of depth H with horizontal and vertical eddy diffusivities D_h and D_z. A wall separates a region dominated by L (α ≈ +1) from one dominated by D (α ≈ −1); its width is of order √(D/k). (c) The scale rule, Eq. 13: the number of molecules N needed for Δ = 2 (a 97.7% probability that the PVED-favoured hand is selected) over the range of g sampled in Table 1, for kτ = 10 (solid) and kτ = 10^5 (dashed). The dot is the reference case, g = 10^-17 and kτ = 10, which needs 3.6×10^33 molecules. Horizontal lines are N = c·V·1000·N_A for a 1 m^3 pond at 1 M and for the Enceladus (2.7×10^16 m^3) and Europa (3.3×10^18 m^3) ocean volumes of Table 1 at 1 µM, each treated as fully mixed. Drawn by `journals/make_fig_schematic.py`, which uses `review.py: n_needed()` and asserts the 3.6×10^33 and 6×10^26 values.

## 2.1 The bias: the parity-violating energy difference

The weak nuclear force does not conserve parity. As a result, the two mirror-image forms (enantiomers) of a chiral molecule are not exactly equal in energy. The difference is the parity-violating energy difference, ΔE_pv. It is set by the weak interaction between electrons and nuclei, and it grows steeply with nuclear charge. For molecules made of light atoms (H, C, N, O) it is of order attoelectronvolts (aeV) to femtoelectronvolts (feV); with sulfur, chlorine or bromine it can reach about 200 aeV (Quack et al., 2022). It has never been measured directly, and its sign for the amino acids in water is not settled (Quack et al., 2022).

I express the bias as the dimensionless ratio of this energy to the thermal energy:

$$ g = \frac{\Delta E_{pv}}{k_B T}. \tag{1}$$

Here k_B is Boltzmann's constant and T is the temperature. In the chemistry of Section 2.2, g enters as a fractional difference between the copying rate constants of the two hands: the L hand copies at k(1 + g) and the D hand at k(1 − g). This follows Kondepudi and Nelson (1985), who proposed that the PVED acts as a small bias on the rates of a self-amplifying reaction.

The reference value is g = 10^-17, taken at 0 °C (T = 273 K) in `review.py: r1_scale()`, which sets ΔE_pv = 10^-17 × R × 273 J/mol. In the uncertainty sampling (Section 2.6), g is drawn log-uniformly from 4×10^-19 to 4×10^-16 (`uncertainty.py: sweep()`), which is about 0.01 to 10 aeV per molecule. This range lies at or below the low end of the values given by Quack et al. (2022) for light-atom molecules, so it does not overstate the bias. Because ΔE_pv is a fixed energy, g falls as 1/T; Section 2.2.4 gives the effect of temperature on the result.

## 2.2 The chemistry and the equation it reduces to

### 2.2.1 A Frank-type network

A bias of 10^-17 has no effect on a single reaction. It can matter only in a reaction that copies its own product, so that a small excess grows, and in which the two hands remove each other, so that only one hand can become dominant (Frank, 1953). I use the following network in a flow reactor (`frank.py: run()`):

| Reaction | Rate (M s^-1) | Physical role |
|---|---|---|
| feed → A | f (a₀(t) − a) | supply of achiral feedstock A; a₀ is raised slowly |
| A → L, A → D | k₀ a each | slow, unbiased, uncatalysed production |
| A + L → 2L | k (1 + g) a L | self-copying of L, with the weak-force bias |
| A + D → 2D | k (1 − g) a D | self-copying of D |
| L + D → inactive product | k_i L D | mutual removal (cross-inhibition) |
| L → out, D → out | k_d L, k_d D | outflow |

Symbols and units: a, L and D are the concentrations of A and of the two hands (M = mol/L); a₀ is the feed concentration (M); f is the flow rate (s^-1); k₀ is the uncatalysed rate constant (s^-1); k is the self-copying rate constant (M^-1 s^-1; this is k₂ in the ocean model); k_i is the cross-inhibition rate constant (M^-1 s^-1); k_d is the outflow rate (s^-1); g is dimensionless (Eq. 1). The rate equations are

$$\frac{da}{dt} = f\,(a_0 - a) - 2k_0 a - k a\left[(1+g)L + (1-g)D\right], \tag{2a}$$
$$\frac{dL}{dt} = k_0 a + k(1+g)\,aL - k_i LD - k_d L, \tag{2b}$$
$$\frac{dD}{dt} = k_0 a + k(1-g)\,aD - k_i LD - k_d D. \tag{2c}$$

In the simulations each reaction also carries its own counting noise (the chemical Langevin equation): a reaction with rate r adds a random change of standard deviation √(r dt/Ω) to the concentrations it changes in a time step dt, where Ω is the number of molecules per unit concentration in the reactor (molecules per M; Ω = 1000·N_A·V for a volume V in m^3). This is the term `n = sqrt(r*dt/omega)*...` in `frank.py: run()`. The simulations use dimensionless constants k = k_d = k_i = 1, k₀ = 10^-3 and f = 1, with a₀ ramped from 0.8 to 1.6 over a time τ.

A second network (`frank2.py: run()`) adds the reverse reactions L → A, D → A (rate constant 0.2) and L + L → A + L, D + D → A + D (rate constant 0.5). I use it only as a check (Section 2.7).

### 2.2.2 Reduction near the symmetry-breaking point

Write the sum and the difference of the two hands as s = L + D and u = L − D. Subtracting Eq. 2c from Eq. 2b gives, exactly,

$$\frac{du}{dt} = (ka - k_d)\,u + g\,k a\, s. \tag{3}$$

Cross-inhibition removes one L and one D together, so it does not appear in Eq. 3. The first term says that a difference between the hands grows when copying outpaces outflow (ka > k_d) and decays otherwise. The second term is the bias.

On the racemic branch (u = 0, g = 0), adding Eqs. 2b and 2c and setting the result to zero gives 2k₀a + (ka − k_d)s − (k_i/2)s^2 = 0. At the crossing point, where ka = k_d, each hand has concentration x = s/2 = √(k₀a/k_i). With the constants above this is x = √(10^-3) ≈ 0.032, as used in `frank.py: main()`.

Because a₀ is raised slowly, a rises and the growth rate of u, λ_u = ka − k_d, passes through zero at a time t_c. Near t_c it is linear in time:

$$\lambda_u(t) \approx \gamma\,(t - t_c), \qquad \gamma = k\,\frac{da}{dt} \quad (\text{s}^{-2}). \tag{4}$$

The noise on u comes from every reaction that changes L − D by one molecule: the two uncatalysed productions, the two copying steps and the two outflows. Cross-inhibition changes L − D by zero and contributes nothing. At L = D = x the variance of u grows at the rate

$$\varepsilon_u = \frac{2k_0 a + 2x\,(ka + k_d)}{\Omega} \quad (\text{M}^2\,\text{s}^{-1}). \tag{5}$$

Near t_c, therefore, du = [γ(t − t_c)u + g_eff] dt + √ε_u dW, with g_eff = 2g·x·ka. These are the three quantities computed in `frank.py: predicted()`. For frank2.py the same bookkeeping gives λ = ka − k₀r − k_d − 2k_r x and ε = (2k₀a + 2k₀r x + 2kax + 2k_r x^2 + 2k_d x)/Ω (`frank2.py: mapped()`).

Once |u| is no longer small, feedstock depletion and cross-inhibition limit its growth. Because the equations are unchanged when L and D are exchanged together with g → −g, the leading saturating term in the excess is cubic. Writing α = (L − D)/(L + D), measuring time in units of 1/k and rescaling α so that the cubic coefficient is 1 gives the normal form

$$d\alpha = \left(\lambda(t)\,\alpha - \alpha^3 + g\right)dt + \sqrt{\varepsilon}\,dW, \qquad \lambda = \gamma t. \tag{6}$$

This is implemented in `sim.py: simulate()`. Its symbols are dimensionless: λ measures how far the system is past the transition; g is the bias; ε is the noise strength; γ = 1/(kτ) is the sweep rate, with k = k₂c the copying rate (k₂ in M^-1 s^-1 times the reactant concentration c in M) and τ the time over which conditions pass through the transition; W is a Wiener process. For λ < 0 the only stable state is near α = 0; for λ > 0 there are two stable states near α = +1 and α = −1 (Figure 1a). I did not derive the cubic coefficient for the network; the result below depends only on the linear stage, and Section 2.7 tests it against the full network directly.

### 2.2.3 Which state is selected

The outcome is set while α is still small, so the cubic term can be dropped. Equation 6 becomes linear,

$$d\alpha = (\gamma t\,\alpha + g)\,dt + \sqrt{\varepsilon}\,dW, \tag{7}$$

and its solution is

$$\alpha(t) = e^{\gamma t^2/2}\left[\alpha(t_0)\,e^{-\gamma t_0^2/2} + \int_{t_0}^{t} e^{-\gamma s^2/2}\left(g\,ds + \sqrt{\varepsilon}\,dW_s\right)\right]. \tag{8}$$

The first term in the bracket vanishes when the sweep starts well before the transition (t₀ ≪ −γ^(-1/2)). For t ≫ γ^(-1/2), the prefactor grows without bound and the bracket tends to a fixed random number Z. The sign of Z sets the final state, because the cubic term only stops α at ±1 and cannot change its sign. Z is Gaussian, with

$$\langle Z\rangle = g\int_{-\infty}^{\infty} e^{-\gamma s^2/2}\,ds = g\sqrt{2\pi/\gamma}, \qquad \operatorname{Var} Z = \varepsilon\int_{-\infty}^{\infty} e^{-\gamma s^2}\,ds = \varepsilon\sqrt{\pi/\gamma}. \tag{9}$$

Both integrals are dominated by |t| ≲ γ^(-1/2), that is |λ| ≲ γ^(1/2): this is the window in which the outcome is set (grey band in Figure 1a). The probability that the PVED-favoured state is selected is P = Pr(Z > 0):

$$P = \Phi(\Delta), \qquad \Delta = \frac{\langle Z\rangle}{\sqrt{\operatorname{Var} Z}} = \sqrt{2}\,\pi^{1/4}\,g\,\gamma^{-1/4}\,\varepsilon^{-1/2}, \tag{10}$$

where Φ is the standard normal cumulative distribution (`sim.py: delta()`). Δ = 2 gives P = 97.7%. Physically, Δ compares the excess that the bias builds up during the window with the excess that counting noise builds up in the same window. A slower sweep (smaller γ) widens the window, and both grow, but the bias grows faster, which is why Δ rises as γ^(-1/4).

### 2.2.4 The scale rule

For a pool of N molecules that react together, the counting noise in α is ε ≈ 1/N, and with time in units of 1/k the sweep rate is γ = 1/(kτ). Equation 10 becomes

$$\Delta = \sqrt{2}\,\pi^{1/4}\,g\,(k\tau)^{1/4}\sqrt{N}, \qquad k = k_2 c. \tag{11}$$

For a body of water of volume V (m^3), of which a fraction f is mixed together during the window,

$$N = c\,V f \cdot 1000 \cdot N_A, \tag{12}$$

where N_A is Avogadro's number and the factor 1000 converts m^3 to litres. Equations 11 and 12 together are `ocean.py: delta_at()`.

The condition for the bias to control the outcome is Δ ≥ 2: the bias accumulated over the window must exceed twice the standard deviation of the noise accumulated over the same window. This is the balance between bias and noise that underlies the proposal of Kondepudi and Nelson (1985). Solving Eq. 11 for N step by step:

1. Square Eq. 11: Δ^2 = 2π^(1/2) g^2 (kτ)^(1/2) N.
2. Set Δ = 2: 4 = 2π^(1/2) g^2 (kτ)^(1/2) N.
3. Solve for N:

$$N \ge \left(\frac{2}{\sqrt{2}\,\pi^{1/4}\,g\,(k\tau)^{1/4}}\right)^2 = \frac{2}{\sqrt{\pi}\;g^2\,(k\tau)^{1/2}}. \tag{13}$$

This is `review.py: n_needed()`. The reaction must also run to completion within the available time, which I require as kτ = k₂cτ ≥ 10 (`ocean.py: c_time()`). With g = 10^-17 and kτ = 10, Eq. 13 gives N ≥ 3.6×10^33 molecules, or 6×10^9 mol (`review.py: r1_scale()`). At 1 µM that is a well-mixed cube of water 18 km on a side; at 1 mM, 1.8 km; at 0.1 M, 0.4 km. A 1 m^3 pond at 1 M holds 6×10^26 molecules, about 7 orders of magnitude too few. The rule depends weakly on everything except g: N scales as g^-2, so as T^2 for a fixed ΔE_pv (3.6×10^33 at 0 °C, 4.3×10^33 at 27 °C, 6.7×10^33 at 100 °C), and as (kτ)^(-1/2), so a sweep 10,000 times slower lowers N by 100 times.

The full network adds noise that the simple count ε = 1/N leaves out: uncatalysed production and outflow change L − D as much as copying does (Eq. 5). Dividing Eq. 10 for the network by Eq. 11 for the same molecules gives Δ_network = F·Δ, with F = √(x/(k₀ + 2x)) = 0.70 at the constants of `frank.py` (`frank.py: main()`); the reversible network gives F = 0.71 (`frank2.py: mapped()`). Since Δ ∝ c^(3/4) at fixed V (Eqs. 11, 12), the required concentration rises by F^(-4/3) = 1.60. The uncertainty sampling draws F from U(0.3, 0.7) (Table 1), whose upper end matches these values and whose lower end allows for chemistry with more wasted turnover.

## 2.3 Why mixing sets the number of molecules

N in Eq. 11 is not the number of molecules in a body of water. It is the number whose counting noise is averaged together while the outcome is set. Two regions that do not exchange molecules during the window each have their own noise, their own Z in Eq. 9, and so they can end in opposite states. The relevant N is therefore the number of molecules in a volume that is connected and mixed on the time scale of the reaction.

To describe this I let α vary in space. Mixing by ocean eddies acts on large scales like diffusion, with a horizontal eddy diffusivity D_h and a vertical one D_z (both m^2 s^-1). Equation 6 becomes the field equation

$$\partial_t \alpha = \lambda\alpha - \alpha^3 + g + D_h\nabla_h^2\alpha + D_z\,\partial_z^2\alpha + \text{noise}, \tag{14}$$

implemented in 1D in `domains.py: sweep()`, in 3D with one D in `domains3d.py: sweep3d()` and `coarsen3d.py: step()`, with separate D_h and D_z in `aniso3d.py: run()`, and on a spherical shell in `sphere_layers.py`.

Ocean mixing is strongly anisotropic: D_z is many orders of magnitude smaller than D_h (Table 1). Changing the vertical coordinate to z' = z·√(D_h/D_z) turns D_z ∂_z^2 into D_h ∂_{z'}^2, so Eq. 14 becomes the equal-mixing equation in a box whose depth is stretched from H to

$$H' = H\sqrt{D_h/D_z}. \tag{15}$$

Weak vertical mixing therefore acts like a deep, narrow ocean. I take the ocean depth as H = V/(4πR^2) with R = L/π, and the horizontal size as half the circumference, L = πR (`uncertainty.py: depth()`). Comparing H' with L sorts every input set into one of two cases:

- H' ≤ L (the stretched ocean is wider than it is deep). After the transition the ocean forms a patchwork of regions of each hand that heals into one hand (Section 2.4). The whole volume counts: f = 1 in Eq. 12.
- H' > L (the stretched ocean is a tall column). The patchwork freezes into stacked horizontal layers of opposite hand (Section 2.4). On a sphere the top layer later takes over the whole depth, so the outcome is set by the molecules of the top layer alone. I assume that this layer has a thickness of about L in stretched units, so it holds a share f = L/H' of the molecules. This share is an assumption, based on the 3-5 layers seen in `aniso3d.py`; two `sphere_layers.py` runs gave top-layer shares of 0.52 and 0.07 against L/H' = 0.25, which neither supports nor rules it out (`uncertainty.py` docstring).

The column case begins when D_z < D_h·(H/L)^2, which is D_z < D_h/333 for Enceladus and D_z < D_h/1386 for Europa (`aniso3d.log`).

## 2.4 Walls, coarsening and why walls stop

**What a wall is.** After the transition (λ > 0), each region of the ocean sits near α = +1 or α = −1. Where two such regions meet, α passes from one value to the other across a thin layer: a wall, or domain boundary. Its width is set by the balance between the reaction, which pushes α towards ±1, and mixing, which smooths it; it is of order √(D/k). Because Eq. 14 (with g = 0) is a gradient flow, a wall carries an energy per unit area σ, which scales as σ ∝ √D (`two_zone.py` docstring).

**Why walls move.** In this kind of gradient flow, a wall moves along its normal at a speed equal to the diffusivity times its mean curvature (Allen and Cahn, 1979):

$$v_n = D\,\kappa, \tag{16}$$

where κ is the sum of the two principal curvatures. A closed patch therefore shrinks and disappears, and the typical patch size grows as

$$l = A\sqrt{D\,t}. \tag{17}$$

In 3D runs from a random start I measured exponents of 0.45 and 0.49 (0.5 expected) and A = 5.48 and 5.60; the sampling uses A = 5.5 (`coarsen3d.py: coarsen()`, `coarsen3d.log`). Because the dynamics do not conserve the amount of either hand, the hand that fills more of the volume takes over: in a 3D run starting with 55% of the volume in one hand, that hand filled 93.0% by t = 400 (`coarsen3d.log`). A patchwork of horizontal size L therefore heals into one hand in a time (L/A)^2/D_h. In the 1D runs, the patch size just after the transition follows l ∝ √D·γ^(-1/4) (`domains.py: main()`, check A), but in 3D this law fails (31% spread, `domains3d.log`) because curved walls keep moving after the transition; I therefore use Eq. 17 instead.

The bias also moves a flat wall: in 1D the region of the PVED-favoured hand advances at

$$v = \frac{3}{\sqrt{2}}\,g\sqrt{D} \tag{18}$$

(`domains.py: wall_speed()`, matched within 10% in two runs). For g ≈ 10^-17 this speed is negligible against curvature-driven motion, so I neglect it in the ocean model.

**Why walls stop.** Equation 16 gives zero speed for a flat wall. In a flat box, the horizontal walls between stacked layers therefore do not move, and the layers persist (`aniso3d.py`: 7 of 8 column runs froze into 3-5 layers, while a cube healed in 4 of 4 runs from 55%). On a sphere, a wall at constant depth is itself a sphere of radius r, with κ = 2/r, and its normal is vertical, so only vertical mixing moves it:

$$v = \frac{2D_z}{r}. \tag{19}$$

I measured 0.991, 0.991 and 0.984 times 2D_z/r, with 0.0% change when D_h was raised four times (`sphere_layers.py: wall_speed()`, `sphere_layers.log`). Each layer wall sinks, the top layer grows, and it takes over the full depth in about H·R/(2D_z), with R = L/π. At D_z = 10^-10 m^2 s^-1 this is 1.5×10^12 yr for Enceladus and 2.8×10^13 yr for Europa, far longer than the age of either ocean; at larger D_z it is short.

A wall also stops where the mixing changes with depth. Its energy E = 4πr^2 σ(r) must fall as it moves. Sinking shrinks its area, which lowers E, but sinking into faster-mixing water raises σ, which raises E. The wall can pass only if

$$\left|\frac{d\ln\sigma}{dr}\right| < \frac{2}{r}. \tag{20}$$

In radial runs (`two_zone.py: run()`) with a slow-mixing layer on top or in the middle, the wall stopped about 0.08 depth units above the lower edge of the slow layer and did not move again, for mixing ratios of 10 and 100, run to 3 times the expected time. The pinning begins at a mixing ratio between 1.02 (the wall passed, 1.378 times slower than without the layer) and 1.05 (it stopped). A bias can push the wall through, but it needs h ≈ 0.20 in model units (`polar_gap.py`, Section 3), while g ≤ 4×10^-16, about 15 orders of magnitude too weak. With the slow layer at the bottom the wall passed, 0.67-0.71 times the predicted time (`two_zone.log`).

## 2.5 The stratified-lid test

Seawater near its freezing point normally becomes denser as it cools, so cooling at the top drives convection. Below a threshold salinity S*, the thermal expansion coefficient at the freezing point changes sign: cooled water becomes lighter, and the ocean stacks a stratified upper layer on a convecting interior (Zeng and Jansen, 2021, 2024). Such a layer has a much smaller D_z than the water below it, so by Eq. 20 it traps the wall between the top layer and the interior, and the ocean keeps two hands.

I define S* as the salinity at which the thermal expansion coefficient at the freezing temperature is zero,

$$\alpha_T\!\left(S^*,\; T_f(S^*, p),\; p\right) = 0, \qquad p = \rho_{ice}\, g_s\, h_{ice}, \tag{21}$$

computed from the TEOS-10 equation of state (McDougall and Barker, 2011) at the pressure p at the top of the ocean, where h_ice is the ice thickness, g_s the surface gravity and ρ_ice = 920 kg/m^3 (`verify_sstar.py: sstar()` recomputes it independently). S* is 22.5 g/kg at the top of the Enceladus ocean (18.0 at its floor), and 3.0 g/kg under 20 km of ice on Europa (19.5 under 4 km). The sampling uses S* = 20 g/kg for Enceladus and 3.0 g/kg for Europa (`strat_sweep.py: SALT`).

There is a pressure cutoff. Under more than about 28 MPa of ice, water at its freezing point sinks as it cools at any salinity, so no lid can form (`verify_sstar.py` finds no S* above 28-29 MPa). On Europa this corresponds to an ice shell of about 23 km.

Each input set is given a salinity drawn log-uniformly (Enceladus 5-20 g/kg, the measured range of Postberg et al. (2009); Europa 1-100 g/kg, an assumed prior). A set with S < S* is counted as pinned, and pinned sets are not counted as successes (`strat_sweep.py: main()`). Salinity is drawn independently of all other inputs. I report the lid results beside the lid-free results, not in place of them.

## 2.6 From the equations to the reported percentages

Every input to Eqs. 11-21 is uncertain, most by orders of magnitude (Table 1 of the paper). For each world I draw 200,000 input sets, log-uniformly over each range unless noted, with a fixed random seed of 0 (`uncertainty.py: logu()`, `sweep()`). An input set is counted as a success, meaning that the PVED-favoured hand is selected across the whole ocean with probability of at least 97.7%, only if all the conditions below hold:

1. the reaction finishes: k₂cτ ≥ 10;
2. either (i) H' ≤ L, the patchwork heals horizontally in time, (L/5.5)^2/D_h ≤ τ, and F·Δ(f = 1) ≥ 2;
   or (ii) H' > L, (L/5.5)^2/D_h ≤ τ, the layers merge in time, H·(L/π)/(2D_z) ≤ τ, and F·Δ(f = min(1, L/H')) ≥ 2;
3. for the lid test only: S ≥ S*.

In case (i), vertical healing, (H/5.5)^2/D_z ≤ τ, follows from H' ≤ L and does not need a separate test. The reported percentage is the share of input sets that meet all the conditions. It is not the probability that the bias set life's hand. It is the share of plausible conditions under which the bias, rather than chance, controls the outcome.

Without the lid test, the shares are 22% for Enceladus, 29% for Europa and 75% for early Earth's open ocean (`uncertainty.log`). If layers never merged (the flat-box bound), they would be 3%, 7% and 40%. With the lid test, Enceladus falls from 21.8% to 0.0% (every set pinned) and Europa from 28.6% to 21.7% (23.9% pinned) (`strat_sweep.log`). Splitting each input range in half shows which input matters most: the share in the lower and upper halves of the D_z range is 0% and 43% for Enceladus, 0% and 57% for Europa, and 52% and 97% for early Earth (`uncertainty.py: sweep()`).

## 2.7 Why the results should be believed, and what is assumed

Each step above has been tested separately. Table M1 lists the checks; VERIFICATION.md gives the commands.

**Table M1.** Checks that support each step.

| Step | Check | Result |
|---|---|---|
| Eq. 10 (normal form, linear analysis) | 4,000 brute-force stochastic runs of Eq. 6 per point (`sim.py: simulate()`) | agree within ±0.005 (draft); the script asserts agreement within 0.03 |
| Eqs. 3-5 (reduction from chemistry) | full Frank network with chemical Langevin noise (`frank.py`), 2,000 runs per point | predicted/simulated P: 0.500/0.485, 0.737/0.733, 0.858/0.854, 0.999/0.974, 0.909/0.898 for slow sweeps (τ ≥ 800); a fast sweep (τ = 200) gave 0.628/0.680, so the formula is conservative there |
| Same, second chemistry | reversible network (`frank2.py`) | 0.500/0.493, 0.802/0.785, 0.967/0.952, 0.951/0.948; F = 0.71, within 1.7 points |
| Selected state stays selected | reaction held at saturation for 100 times longer than the selection window | share on the favoured hand stayed at 0.730; the flip rate matches Kramers' escape rate (√2/2π)·e^(−1/(2ε)) within 5% (Kramers, 1940) |
| Eq. 17 (coarsening) | 3D runs (`coarsen3d.py`) | A = 5.48 and 5.60; halving the grid spacing changed patch sizes by < 5% |
| Eq. 19 (layer walls on a sphere) | `sphere_layers.py` | 0.991, 0.991, 0.984 × 2D_z/r; independent of D_h |
| Eq. 20 (pinning) | `two_zone.py`, predictions written before the runs | no-layer control at 1.043 times the predicted time (1.023 at half the time step); pinning between ratios 1.02 and 1.05 |
| Eq. 21 (S*) | independent TEOS-10 recomputation (`verify_sstar.py`, gsw 3.6.23) | matches to the rounding of the ice pressure (Europa 3.2 vs 3.0 g/kg at 24.2 MPa; Dione 5.75 vs 5.6) |
| Sampling (Section 2.6) | seeds 1-5 (`verify_seeds.py`); two blind rebuilds from the paper text alone | every world within 0.3 points of seed 0; rebuilds within 0.4 points |
| Counting rule at lab scale | reanalysis of published amplifiers (`soai_check.py`, `crystal_check.py`, `viedma_check.py`) | Soai 46 of 91 and crystal 37 of 73 split as chance predicts; 12 of 25 Soai pairs opposite (P = 0.018); stirred NaClO₃ behaves as n_eff ≈ 1.0-1.1 pools, unstirred about 17 |

The numerical scheme for all spatial runs is explicit Euler-Maruyama time stepping with the 7-point Laplacian, and the time step meets the 3D stability limit D·dt/dx^2 ≤ 1/6 (dt = min(0.05, 0.1/D) in `domains3d.py: sweep3d()` and `coarsen3d.py: coarsen()`). Thirteen predictions written before the runs in `breakeven.py`, `assumed_rules.py` and `earth_breakeven.py` all held (VERIFICATION.md).

The lab checks support the counting rule (a bias controls the outcome only when it exceeds the chance excess of about 1/√N) at scales of 10^6-10^12 units. They do not test the weak force, the 10^33-molecule threshold or any ocean.

**Main assumptions.**

1. *An amplifying chemistry exists.* No measured prebiotic reaction is known to amplify handedness in this way. The model asks what such a reaction would need, not whether it exists.
2. *The bias is the PVED, as a rate bias.* Eq. 1 identifies the fractional rate bias with ΔE_pv/kT. How the energy difference of the molecules enters the rate constants is not derived here; a factor of order 1 in this identification is small against the three orders of magnitude sampled. The sign for amino acids is unsettled.
3. *Mixing is diffusive.* Eddies are represented by D_h and D_z. Currents that carry a wall are not in the main model; separate runs with an overturning current are reported in Section 3.
4. *The top layer has thickness about L in stretched units* (Section 2.3, case ii). The two rules that are assumed rather than measured move a moon's share by at most 2.2 points (`assumed_rules.py`).
5. *Uniform priors in log space.* The percentages depend on the input ranges. The input ranges and the reason for each are given in Table 1; the vertical range is the one modelled for icy-moon oceans by Zeng and Jansen (2021).
6. *Salinity and lid are independent of the other inputs,* and a pinned set never counts as a success.

---

### Source conflicts (for the author; remove before submission)

These are differences between sources that I found while writing this section. I did not resolve them; the text above uses the value in ija_draft.md where one exists.

1. Coarsening constant: the code uses A = 5.5; the measured values are 5.48 and 5.60 (already noted in the draft).
2. `uncertainty.py` docstring lists a vertical healing rule, (H/A)^2/D_z ≤ τ, that `sweep()` does not apply. It follows from H' ≤ L in case (i), so it has no numerical effect.
3. Column onset for Enceladus: `aniso3d.log` prints D_h/333; the code asserts 334 within 5%.
4. Europa S*: 3.0 g/kg (strat_sweep.py, the draft) vs 3.2 at 24.2 MPa and 2.55 at about 25 MPa (VERIFICATION.md). The strat_sweep.py docstring says "~25 MPa, 20 km ice", while verify_sstar.py gives 24.2 MPa for 20 km.
5. Ice density: the draft uses 920 kg/m^3; `verify_sstar.py` uses 917 kg/m^3 with g_s = 1.315 m/s^2 for Europa. The Europa gravity value is not stated in the draft.
6. Lid cutoff: "about 28 MPa" (draft) vs "no S* above 28-29 MPa" (VERIFICATION.md) vs "about 23 km of ice" for Europa (draft Section 3.6). These agree if g_s ≈ 1.3 m/s^2, but the draft does not state it.
7. Agreement of Eq. 10 with simulation: "±0.005" (draft) vs the asserted tolerance of 0.03 in `sim.py: main()`.
8. Range of g: Table 1 uses 4×10^-19 to 4×10^-16; `review.py: r1_scale()` prints "~1e-19..1e-16 in g" as the modern range.
9. PAPER_DRAFT.md Section 3.3 reports 93% at t = 400 (128^3 grid) for the 55% start; its Figure 3 caption reports 99% (96^3 grid). This section quotes 93.0% from coarsen3d.log.
10. The phrase "10^33-10^34 molecules" (abstract) is not tied in the draft to a stated range of g and kτ; the 3.6×10^33 to 6.7×10^33 values come from T = 0-100 °C at g = 10^-17 (273 K) and kτ = 10.
