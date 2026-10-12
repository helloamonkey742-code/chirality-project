<!--
=====================================================================
SUBMISSION NOTES (remove before submitting)
=====================================================================

Target: The Planetary Science Journal (PSJ), AAS / IOP. Checked 2026-10-11.

A. ARTICLE TYPE AND FORMAT
- Article type chosen: regular PSJ research Article.
- Word limit: NONE FOUND for PSJ Articles (UNCONFIRMED). The 3,500-word cap belongs to ApJ Letters
  and does not apply. Length is paid for by the word (section C), so shorter is cheaper.
  Source: https://journals.aas.org/manuscript-preparation/ (no PSJ article-type list or limit on that page).
- Abstract: one paragraph, no more than 250 words. This draft's abstract: see count at the end of these notes.
  Source: https://journals.aas.org/manuscript-preparation/
- Running head: 44 characters or fewer. This draft: "Weak-force handedness in ocean worlds" (37).
  Source: https://journals.aas.org/manuscript-preparation/
- Keywords: Unified Astronomy Thesaurus (UAT) terms are collected at acceptance and need not be in the
  manuscript. Source: https://journals.aas.org/manuscript-preparation/
- Citations are author-year. Reference title strings are required.
  Sources: https://journals.aas.org/manuscript-preparation/ ; https://journals.aas.org/psj-policies/
- Figures: EPS preferred; PNG/JPG/TIFF acceptable; at least 1000 px wide; no multi-page files; each figure
  at most 3 printed pages. Source: https://journals.aas.org/manuscript-preparation/
  All four PNGs used here pass the width rule: sketch.png 1500 px, threshold.png 1050 px,
  coarsening.png 1800 px, rnaas_fig.png 1950 px.
- PSJ is prepared like other AAS articles (AASTeX). This markdown must be converted to AASTeX before
  upload [TODO]. Source: https://journals.aas.org/psj-policies/
- PSJ is fully gold open access (AAS + DPS), CC-BY 4.0, generally two reviewers.
  Sources: https://journals.aas.org/planetary-science-journal/ ; https://journals.aas.org/psj-policies/

B. DUAL-ANONYMOUS REVIEW (identity risk)
- Dual-anonymous review is the PSJ standard. Authors may opt out in the "Manuscript Comment" box.
  The editorial office removes names, affiliations and acknowledgments.
  Source: https://journals.aas.org/psj-policies/
- RISK: the GitHub URL (user name, licence file), the planned Zenodo record and any self-citation of the RNAAS
  note reveal who I am. Choose one before submitting: (a) opt out of anonymity in the Manuscript Comment box,
  or (b) give reviewers an anonymized code link. [DECIDE]

C. CHARGES
- 2026 PSJ rates: $69.51 per 350-word quantum, plus $61.40 per figure, table or online component.
  "Publication fees are required." Contact: journals.manager@aas.org.
  Source: https://journals.aas.org/article-charges-and-copyright/
- Estimate for this draft: see "FEE ESTIMATE" at the end of these notes.
- Waiver: the AAS Publication Support Fund lists only AJ, ApJ, ApJS and ApJL; PSJ is NOT listed. Requests are
  made at submission, and priority goes to authors without funding.
  Sources: https://journals.aas.org/waivers/ ; https://journals.aas.org/support
  An AAS post (4 April 2022) mentions waivers of "5-100%" depending on the journal:
  https://journals.aas.org/?p=2914
  DOAJ says PSJ has no waiver policy and an APC up to $5,315 (record reviewed Jan 2025, may be stale):
  https://doaj.org/toc/2632-3338
  PSJ WAIVER STATUS: UNCONFIRMED. Before submitting, email journals.manager@aas.org (or OAquestions@aas.org)
  and ask whether an unfunded high-school author can get a PSJ fee waiver.

D. AI POLICY
- AAS author guidelines on LLMs (posted 9 September 2026): https://journals.aas.org/author-llm-guidelines
  Announcement: https://journals.aas.org/?p=4629 ("required AI-usage disclosure statement for all
  submissions", in force "for new submissions later this fall").
- Key rules: disclose at submission; the examples put the statement at the end of the article; AI cannot be
  an author; LLMs may not generate any portion of the text unsupervised; the authors are fully responsible.
  Example wording used as the model: "AI tools were used to generate part of the manuscript text. The GAI tool
  used was ..." and "Responsibility for the final manuscript text lies entirely with the authors."
- Illustrations made with AI tools must say so in the caption (tool, and prompt if possible). My figures were
  drawn by matplotlib scripts that an AI assistant helped write; they are not AI-generated images. The policy
  gives no rule on AI-written code. UNCONFIRMED whether captions need a note; I put the disclosure in the
  Acknowledgments only. Ask the editors if unsure.
- [TODO] Fill in the exact Claude model name and version in the Acknowledgments.

E. SOFTWARE AND DATA
- Software policy: cite the software and a DOI for the version used; an open licence is "highly
  recommended"; archiving on Zenodo/FigShare for a DOI is recommended, not required; use the \software
  macro below the acknowledgments. Source: https://journals.aas.org/policy-statement-on-software/
- Data guide: "Software and notebooks should be published, placed in a persistent archive, and cited".
  Use the AAS Journals Zenodo community; a licence is required and the record description must not be blank.
  Source: https://journals.aas.org/data-guide/
  No separate data-availability-statement rule found (UNCONFIRMED); contact data-editors@aas.org.
- [TODO] Archive the repository on Zenodo (AAS Journals community) and put the DOI in the paper.
- Software list = packages actually imported by the scripts used in this paper: numpy, scipy, matplotlib,
  gsw (TEOS-10 check, verify_sstar.py only), libRoadRunner + Antimony (SBML cross-check, frank_sbml.py only).
  Their own software papers are not cited yet [TODO: add the numpy/scipy/matplotlib/gsw/libRoadRunner papers
  with verified DOIs, as the AAS software policy asks]. Not added now because the task allows references
  only from PAPER_DRAFT.md, RNAAS_DRAFT.md and limits_lit.md, and none of them lists these papers.

F. PRIOR RNAAS NOTE
- No explicit AAS policy found on turning a Research Note into a full article (UNCONFIRMED). Pages checked:
  https://journals.aas.org/research-notes/ (Research Notes: up to 1,500 words, not peer reviewed)
  https://journals.aas.org/scope-of-research-notes-of-the-american-astronomical-society/
  https://journals.aas.org/research-note-preparation-guidelines/
  https://journals.aas.org/policies/
- The ethics page asks authors to certify the manuscript "has not been previously published elsewhere ... in
  any form", and it counts "repeating text from previously published papers by the author" as
  self-plagiarism: https://journals.aas.org/professional-and-ethical-standards-for-the-aas-journals/
- Advice: if the RNAAS note is published first, cite it, do not reuse its sentences word for word, say so in
  the cover letter, and ask the PSJ editors whether this is OK before submitting. RNAAS_DRAFT.md is still a
  draft; whether it has been submitted is unknown to me. Several abstract and ranking sentences here follow
  the RNAAS draft closely and must be reworded if the note is published. [DECIDE]

G. NUMBER CONFLICTS BETWEEN SOURCES (not resolved silently; RNAAS_DRAFT.md governs where it speaks)
1. Europa S*: 3.0 g/kg (PAPER, strat_sweep.py) vs 3.2 (verify_sstar.py, 24.2 MPa). Dione S*: 5.6 (PAPER)
   vs 5.75 (verify_sstar.py). Draft uses 3.0 and 5.6 (PAPER/RNAAS side).
2. Magnetite worst case for icy moons: RNAAS "2-4 x 10^-15" vs log values 1.6e-15 and 3.5e-15. Draft uses RNAAS.
3. Titan inherited excess that overrides the weak force: limits_lit "~1e-16" vs impact_melt.log 1e-15 /
   2e-16 / 6e-17 (1 km crater / Selk / Menrva). Draft gives the three log values.
4. Effect of sideways mixing D_h: uncertainty.log -3 / -2 / -9 (Enceladus/Europa/Earth) vs PAPER "-2 to -9".
   Draft uses the log values.
5. Enceladus range: "0-8%" (lid choices, PAPER section 3.14) vs "0-18%" (all published salinities, RNAAS).
   Consistent but different scopes; both are stated with their scope.
6. Rounding of world scores: PAPER gives Callisto 28.5%, Triton 28%, Mimas 5%; ocean_worlds.log gives 28.5,
   27.7, 4.7. The table uses the log decimals; the text uses RNAAS rounding (Ganymede 28%).
7. Ozturk et al. 2023: two papers. PAPER's ~60% excess cites the Sci. Adv. paper (2023a); RNAAS's magnetized
   surfaces sentence cites the Nat. Commun. paper (2023b). Draft keeps that split.
8. PAPER abstract contains a garbled sentence ("Across every known or suspected ocean world, Among confirmed
   oceans ..."). Not copied.
9. PAPER's "Multi-Plume Handedness" test names Europa Clipper SUDA, but PAPER section 3.10 says SUDA cannot
   measure chirality (mirror molecules have equal mass). Draft drops the SUDA claim and keeps the 3.10 statement.
10. Coarsening figure caption in PAPER says l = 5.5 sqrt(Dt); the corrected text gives A = 5.48 / 5.60.
    Draft caption says l = A sqrt(Dt) with exponent 0.45 (coarsen_fig.log) and gives A in the text.
11. NOT a conflict (confirmed by the coordinator): the 99% at t = 400 (coarsen_fig.log 0.989) is the 96^3 figure
    run; PAPER's 93% is a separate 128^3 run. The draft uses 99% with Figure 3 and does not quote the 93% run.
12. Hochberg et al. 2023 wording comes from RNAAS; limits_lit has no Hochberg text.

H. REFERENCES
- DOIs checked against Crossref with curl -s https://api.crossref.org/works/<doi> (no mailto), title matched:
  every DOI in the reference list below resolves with a matching title. One item is not in Crossref:
  Kang & Zhang 2026 (arXiv DOI 10.48550/arXiv.2603.22602 is DataCite); verified on arXiv (title "Subsurface
  ocean salinity and dissipation rate inferred from Enceladus ice shell morphology") and cited by arXiv ID.
  McDougall & Barker 2011 has no DOI (SCOR/IAPSO WG127 manual).
- Dropped or not used, with reason:
  * Wakita et al. 2023 (Selk melt volumes): no DOI or full reference in the allowed sources. Selk is named
    without a citation [TODO: find and verify the reference].
  * Peter et al. (HCN in the Enceladus plume): Crossref year 2023, PAPER says 2024. Year conflict; dropped.
  * Zhang, Kang & Marshall 2026 and all personal communications: not needed / not citable.
  * Kempf et al. 2025 (SUDA): only supported the dropped SUDA claim.
  * Everything in FRANK_1953_REVIEW.md and test_against_actual_findings.py: excluded (fabricated references).
  * PAPER's other references (Truong, Choblet, Daval, Zega, Connolly, etc.): not needed for this shorter paper.
- PAPER_DRAFT.md credits A. Brandenburg, K. Soderlund and W. Kang via personal communications. They are NOT
  acknowledged here because the task allowed only Dr. McKay "for discussion". [DECIDE whether to ask them
  and acknowledge them.]
- [CONFIRM WITH DR. MCKAY BEFORE SUBMISSION] the acknowledgment wording.

I. OTHER TODOs
- [TODO] Corresponding-author email and ORCID.
- [TODO] Convert to AASTeX, add \software, check figure placement, generate BibTeX with title fields.

COUNTS AND FEE ESTIMATE (counted with a whitespace split, 2026-10-11)
- Abstract: 230 words (limit 250).
- Abstract through Figures (no references): 6,945 words. With references: 8,242 words. Whole file with these
  notes: 9,775 words.
- Items: 5 tables (Tables 1-4 + the section-to-script table) + 4 figures = 9.
- Fee if references are not counted: ceil(6945/350) = 20 quanta x $69.51 = $1,390.20, plus 9 x $61.40 = $552.60,
  total about $1,942.80.
- Fee if references are counted: 24 x $69.51 = $1,668.24 + $552.60 = about $2,220.84.
  Which words AAS counts is UNCONFIRMED; the true fee is set at acceptance.
- Cheapest cuts: drop the section-to-script table into the Zenodo README (-$61.40); trim Sections 4.2 and 6.
- Ask journals.manager@aas.org about a waiver BEFORE submitting (section C).
=====================================================================
-->

# Which Ocean World Could Let the Weak Force Choose Life's Handedness? Stratification, Vertical Mixing and a Ranking of Solar System Oceans

**Saptarshi Ghosh**
Cupertino High School, Cupertino, CA, USA
[TODO: corresponding-author email; ORCID]

Running head: Weak-force handedness in ocean worlds

UAT keywords: collected at acceptance [TODO]

## Abstract

The weak nuclear force makes one mirror-image form of a chiral molecule slightly lower in energy, by about 10⁻¹⁷ kT. A slow, self-amplifying reaction can turn this parity-violating energy difference (PVED) into a one-handed result, but only if about 10³³–10³⁴ molecules decide together. That needs kilometres of well-mixed water, so the question becomes one of ocean physics. I couple a stochastic model of chiral amplification to an ocean in which regions of opposite hand coarsen by eddy diffusion. Three physical results set the answer. First, weak vertical mixing freezes the ocean into stacked layers of opposite hand. On a round moon these merge in time only if the vertical diffusivity exceeds H·R/(2t): about 1.5×10⁻⁷ to 1.5×10⁻⁴ m² s⁻¹ for Enceladus and 7×10⁻⁷ to 3×10⁻⁵ m² s⁻¹ for Europa. Published estimates straddle this bar. Second, a stably stratified layer under the ice pins the boundary between hands. From the TEOS-10 equation of state, such a lid forms below a threshold salinity, and it cannot form under more than about 28 MPa of ice. Third, without a lid, thinner oceans merge their layers faster. I sweep every uncertain input over 200,000 draws per world for fourteen ocean worlds. Of the confirmed oceans, Ganymede (28% of plausible cases) and Europa (22–26%) score highest, and Enceladus scores lowest (0–18%). Early Earth's open ocean scores 75%. The result is a ranking, not an absolute probability.

## 1. Introduction

Life uses left-handed amino acids and right-handed sugars. Why one hand and not the other is still open. One candidate is the weak nuclear force. It makes the two mirror forms of a chiral molecule differ in energy by a tiny amount, the parity-violating energy difference (PVED), about 10⁻¹⁷ kT. The PVED has never been measured directly, and whether it or chance chose life's hand "remains completely open" (Quack et al. 2022).

Kondepudi & Nelson (1985) proposed that the PVED can select life's hand in a large, slow, self-amplifying reaction. The classic amplifier is the Frank (1953) scheme: each hand copies itself, and the two hands destroy each other. Stochastic models show that a small bias selects the hand when it exceeds the noise (Hochberg et al. 2023). Lente (2007) found that the PVED cannot be amplified when self-copying starts from a few chiral molecules. My model instead assumes an open, fed reaction in which a whole population chooses at once. Lunine et al. (1999) discussed abiotic homochirality on Titan. Their chapter asks whether the PVED decides the hand during amplification and proposes Titan as a natural laboratory. It also notes that a visible PVED effect needs a system much larger than a laboratory. I consider liquid-water oceans and ask how large and how well mixed they must be.

The answer turns on ocean physics. A bias of 10⁻¹⁷ beats molecular counting noise only when about 10³³–10³⁴ molecules decide together. That is kilometres of water. Ponds and rock pores are far too small. Icy-moon oceans are large enough, but they are not well mixed. They are covered by ice, may be stratified, and mix slowly in the vertical. So the chemistry is the input, and the result depends on how these oceans mix.

In this paper I do three things. I derive the number of molecules that must decide together (Section 3.1). I show how weak vertical mixing and a stratified lid control whether an ocean ends with one hand (Sections 3.2–3.3), and I compare the needed vertical mixing with published estimates. I then sweep every uncertain input for fourteen Solar System ocean worlds and rank them (Sections 3.4–3.6). Section 4 reconciles the model with racemic samples from asteroids Bennu and Ryugu. Section 5 lists the limitations, and Section 6 describes a planned laboratory test.

## 2. Model

### 2.1. The normal form of chiral amplification

Let α = (L − D)/(L + D) be the excess of one hand. Near the point where a self-amplifying reaction switches on, any such reaction reduces to one equation, the normal form (Figure 1):

$$d\alpha = \left(\lambda(t)\,\alpha - \alpha^3 + g\right)dt + \sqrt{\varepsilon}\,dW, \qquad \lambda = \gamma t . \tag{1}$$

Here λ measures how far past the switch-on point the system is. For λ < 0 the only stable state is 50/50; for λ > 0 there are two stable states, all-L and all-D. The bias is g = ΔE_pv/kT, where ΔE_pv is the PVED per molecule. The noise from counting N molecules is ε ≈ 1/N. The sweep rate is γ = 1/(kτ), where k = k₂c is the self-copying rate (rate constant k₂ times concentration c) and τ is the time the sweep takes. Time is in units of 1/k.

Linearising around λ = 0 gives the probability that the hand favoured by the weak force wins:

$$P = \Phi(\Delta), \qquad \Delta = \sqrt{2}\,\pi^{1/4}\, g\, \gamma^{-1/4}\, \varepsilon^{-1/2} = \sqrt{2}\,\pi^{1/4}\, g\,(k\tau)^{1/4}\sqrt{N}, \tag{2}$$

where Φ is the normal cumulative distribution. Δ = 2 means a 97.7% win. I checked Equation (2) against 4,000 brute-force runs per point; the agreement is ±0.005.

### 2.2. Molecule count and the win rule

For a body of water, N = c·V·f·1000·N_A, where V is the volume in m³, f is the fraction that is well mixed on the decision time, and N_A is Avogadro's number. The weak force "wins" a draw only when two conditions both hold: the size condition Δ ≥ 2, and the time condition k₂cτ ≥ 10 (the reaction actually finishes).

### 2.3. Space and mixing

With mixing, α becomes a field α(x, t):

$$\partial_t \alpha = \lambda\alpha - \alpha^3 + g + D_h\nabla_h^2\alpha + D_z\,\partial_z^2\alpha + \text{noise}, \tag{3}$$

where D_h and D_z are the horizontal and vertical eddy diffusivities. Rescaling depth by √(D_h/D_z) turns Equation (3) exactly into the equal-mixing equation in a box that is taller by that factor. So weak vertical mixing acts like a deep, narrow ocean, with stretched depth H' = H·√(D_h/D_z). After the decision, walls between regions of opposite hand move at D times their curvature (Allen & Cahn 1979). Closed patches shrink, patch size grows as l = A√(Dt), and the majority hand takes over. Flat walls between horizontal layers do not move.

All spatial runs use explicit Euler time stepping (Euler–Maruyama when noise is included) and a 7-point Laplacian. Units are 1/k for time and √(D/k) for length, with grid spacing 1 and time step 0.05. This meets the 3D stability limit D·dt/dx² ≤ 1/6. Sidewalls are periodic; top and bottom are closed.

### 2.4. The ocean sweep

For each world I draw 200,000 input sets. Every range is log-uniform unless noted (Table 1). Ocean depth is H = V/(4π(L/π)²). A draw wins when k₂cτ ≥ 10 and either

- (i) H' ≤ L, (L/5.5)²/D_h ≤ τ and F·Δ(1) ≥ 2; or
- (ii) H' > L, (L/5.5)²/D_h ≤ τ, H·(L/π)/(2D_z) ≤ τ and F·Δ(min(1, L/H')) ≥ 2,

where Δ(f) = √2 π^¼ g (k₂cτ)^¼ (c·V·f·1000·N_A)^½. Case (i) is an ocean that heals sideways to one hand. Case (ii) is a layered ocean. Its layers must merge in time, and the top layer decides with its own molecules, a share L/H' of the ocean (an assumption tested in Section 3.4). The stratified lid of Section 2.6 then removes draws with salinity below a threshold S\*, drawn independently. All quantities are in SI units, with 1 yr = 3.156×10⁷ s.

**Table 1.** Sweep inputs. Ranges are log-uniform unless noted.

| Input | Enceladus | Europa | Early Earth ocean |
|---|---|---|---|
| bias g = PVED/kT | 4×10⁻¹⁹ – 4×10⁻¹⁶ | same | same |
| rate k₂ (M⁻¹ s⁻¹) | 10⁻⁶ – 1 | same | same |
| concentration c (M) | 10⁻¹² – 10⁻⁴ | 10⁻¹² – 10⁻⁴ | 10⁻⁷ – 3×10⁻⁴ |
| horizontal D_h (m² s⁻¹) | 10⁻² – 10² | 10⁻² – 10² | 10² – 10⁴ |
| vertical D_z (m² s⁻¹) | 10⁻¹⁰ – 10⁻³ | 10⁻¹⁰ – 10⁻³ | 10⁻⁶ – 10⁻⁴ |
| time τ (yr) | 10⁶ – 10⁹ | 10⁸ – 4×10⁹ | 10⁷ – 5×10⁸ |
| volume V (m³), each ×U(0.7, 1.3) | 2.7×10¹⁶ | 3.3×10¹⁸ | 1.332×10¹⁸ |
| length L = πR (R in km) | R = 232 | R = 1460 | R = 6371 |
| chemistry factor F | U(0.3, 0.7) | same | same |
| salinity (g kg⁻¹), S\* | 5–20, S\* = 20 | 1–100, S\* = 3.0 | no lid |

The Enceladus volume is derived from ice-shell and core geometry (Čadek et al. 2016).

### 2.5. Real chemistry

I checked the normal form against a full Frank-type network in a flow reactor. An achiral feedstock A is fed in and ramped slowly upward. A turns into L or D at a slow, unbiased rate k₀a. Each hand copies itself from A at k(1 ± g)·a·(L or D). L and D destroy each other at k_i·LD, and both wash out at k_d. Every reaction carries its own counting noise (the chemical Langevin equation). The network follows Equation (2) to within about 1 percentage point in the slow, natural limit. Its selection strength is lower by a factor F ≈ 0.7, because every reaction that adds or removes a chiral molecule adds noise. A second scheme with reverse reactions gives F = 0.71. The sweep therefore draws F from U(0.3, 0.7).

### 2.6. The stratified lid

Below a salinity S\*, water at its freezing point gets lighter as it cools, and the ocean stacks a stratified upper layer on a convecting interior (Zeng & Jansen 2021, 2024). I compute S\* from the TEOS-10 equation of state (McDougall & Barker 2011) at the pressure of the ocean top, 920 kg m⁻³ × g × ice thickness. Under more than about 28 MPa of ice, freezing-point water sinks as it cools at any salinity, so no thermal lid can form. For Europa this happens under about 23 km of ice. For Ganymede and Callisto, TEOS-10 is extrapolated past its 100 MPa range.

## 3. Results

### 3.1. How many molecules must decide together?

Setting Δ = 2 in Equation (2) gives the number of molecules that must react as one well-mixed pool:

$$N \ge \left(\frac{2}{\sqrt{2}\,\pi^{1/4} g\,(k\tau)^{1/4}}\right)^2 . \tag{4}$$

For g = 10⁻¹⁷ (the middle of the modern PVED range) and a reaction that only just finishes (kτ = 10), N ≥ 3.6×10³³ molecules, or 6×10⁹ mol. At 1 µM that is a well-mixed cube of water 18 km on a side; at 1 mM, 1.8 km; at 0.1 M, 0.4 km. A 1 m³ pond at 1 M holds only 6×10²⁶ molecules, about 7 orders of magnitude too few.

Cold helps only a little. The PVED is a fixed energy, so g = ΔE/kT is larger in the cold, and the needed N scales as T²: 3.6×10³³ at 0 °C, 4.3×10³³ at 27 °C and 6.7×10³³ at 100 °C. Time also helps only weakly, as (kτ)^¼: a sweep 10,000× slower lowers N by 100×. Concentration alone cannot rescue a small setting. Only ocean-scale water that is mixed within the decision time pools enough molecules. Figure 2 shows the minimum concentration for a small lake and for three oceans.

### 3.2. Layers, and the vertical-mixing bar

**Patches heal sideways.** In 3D, curved walls between patches of opposite hand straighten out, so small patches shrink even with no bias. Patch size grows as l = A√(Dt), with exponents 0.45 and 0.49 (theory 0.5) and A = 5.48 and 5.60 at two mixing strengths. Starting from a 55% majority, a 96³ cube reaches 99% of the favoured hand by t = 400 (Figure 3, top). Sideways healing within a moon's age needs only D_h ≥ ~6×10⁻⁴ m² s⁻¹ (Enceladus) or ~2×10⁻⁴ m² s⁻¹ (Europa). That is about 180–450× below the ~0.1 m² s⁻¹ sideways eddy mixing estimated for Enceladus by Zhang et al. (2024). Sideways mixing therefore never limits the result.

**Weak vertical mixing makes layers.** The stretched ocean is deeper than it is wide (H' > L) in 92–96% of the plausible range. In such tall boxes the patchwork froze into 3–5 stacked layers in 7 of 8 runs and did not change from t = 500 to t = 1000. A cube healed to one hand in 4 of 4 runs from a 55% start. The frozen layers in Figure 3 (bottom) do not change from t = 400 to t = 5000.

**On a sphere, layers merge only through vertical mixing.** On a round moon a layer boundary is a sphere, and "down" follows the radius. So a boundary at constant depth feels only vertical mixing and sinks at 2·D_z/R. I measured sinking speeds of 0.991, 0.991 and 0.984 times 2·D_z/R, with no change when D_h was multiplied by 4. Layers merge in time when H·R/(2D_z) ≤ τ. This holds in 33% (Enceladus) and 34% (Europa) of draws. At D_z = 10⁻¹⁰ m² s⁻¹ the merge would take up to 1.5×10¹² yr (Enceladus) or 2.8×10¹³ yr (Europa).

Inverting the merge condition gives a bar on the vertical diffusivity, D_z ≥ H·R/(2t):

- Enceladus (H = 39.9 km, R = 232 km): D_z ≳ 1.5×10⁻⁷ m² s⁻¹ for an ocean 1 Gyr old, up to 1.5×10⁻⁴ for one 1 Myr old.
- Europa (H = 123 km, R = 1460 km): D_z ≳ 7×10⁻⁷ m² s⁻¹ for 4 Gyr, up to 3×10⁻⁵ for 100 Myr.

**What the literature says.** No one has measured vertical mixing in an icy-moon ocean. Table 2 lists the published values. The one primary estimate (Zeng & Jansen 2021) works back from the tidal and libration energy that might stir the ocean, 3×10⁻¹² to 3×10⁻⁵ W m⁻². That energy is uncertain by about 7 powers of ten, so κ_z spans about 3×10⁻¹⁰ to 3×10⁻³ m² s⁻¹. The sweep range of 10⁻¹⁰ to 10⁻³ matches it. Dissipation converts to diffusivity through the Osborn relation κ_v = Γε/N², with Γ ≈ 0.2 (Osborn 1980). Jansen et al. (2023) call the energy supply "highly uncertain".

**Table 2.** Published vertical diffusivities for icy-moon oceans, compared with the merge bar.

| Source | κ_z (m² s⁻¹) | Kind |
|---|---|---|
| Zeng & Jansen (2021), Enceladus | ≈ 3×10⁻¹⁰ to 3×10⁻³ | scaling from tidal and libration energy |
| same paper | molecular floor: 10⁻⁷ (thermal), 10⁻⁹ (tracer) | physical floor |
| same paper | 5×10⁻⁵ | model input for their low-salinity run |
| Ames et al. (2025) | "~10⁻⁷ to 10⁻³" | not independent (cites Zeng & Jansen 2021); runs use 10⁻⁵, 10⁻⁴, 10⁻³ |
| Jansen et al. (2023) | calls 5×10⁻⁵ "relatively large" | discussion of other papers' inputs |
| Kang et al. (2022a), as quoted by Jansen et al. (2023) | 5×10⁻³ | model input (secondary source) |
| Zhang et al. (2025) | none; seafloor-to-ice transit "tens of years or more" | convection from uneven bottom heating |
| Kang & Zhang (2026, preprint) | < 10⁻³ unless salinity is around 10 psu | upper bound from ice-shell shape |
| Merge bar, this work | Enceladus 1.5×10⁻⁷ – 1.5×10⁻⁴; Europa 7×10⁻⁷ – 3×10⁻⁵ | D_z ≥ H·R/(2t) |

No Europa-specific value exists, so I transfer the Enceladus range to Europa. This is an assumption. Values below 10⁻⁹ m² s⁻¹ (the bottom ~14% of the sampled range) really mean "only molecular spreading". This changes no result, because all those draws already fail the merge bar. The published range straddles the bar. The values that modellers use in simulations (5×10⁻⁵ to 5×10⁻³) clear it for any ocean age, but they are partly chosen to make the simulations run. The literature therefore does not settle the answer. It confirms that vertical mixing is the unknown that decides it.

Other ocean models support this picture. Salinity sets how stratified these oceans are (Kang et al. 2022a; Zeng & Jansen 2024). In a stratified Enceladus ocean, material takes at least hundreds of years to travel from the floor to the ice (Kang et al. 2022b; Ames et al. 2025). An unstratified ocean heated from below is the fast-mixing end (Bire et al. 2022). None of these gives a vertical diffusivity directly, so I keep the full range.

### 3.3. A stratified lid pins the wall

A real ocean is not one mixing zone. It can have a convecting interior with fast vertical mixing under a stratified layer with slow mixing. I ran a set of experiments on this, with a written prediction before each run.

**A slow layer on top pins the wall.** The model is a gradient flow in which a wall's energy per unit area grows as √D. A wall sinking toward the centre shrinks in area, which helps it. Entering faster-mixing water raises its energy, which stops it. It can pass only if |d ln σ/dr| < 2/r. With the slow layer at the bottom, or with no layer, the wall sinks to the floor and one hand takes over. The control matched the sinking rate 2D/r within 4%. With the slow layer on top or in the middle, the wall stopped about 0.08 depth units above the layer's lower edge and never moved again (mixing ratios 10 and 100, run to 3× the expected time). Pinning starts at a mixing ratio between 1.02 and 1.05, so any real stratification is far past it. A bias that favours one hand can push the wall through, but it needs h ≈ 0.20 in model units. The weak-force bias is g ≤ 4×10⁻¹⁶, about 15 orders of magnitude too weak.

**When is the top stratified?** From TEOS-10, S\* is 22.5 g kg⁻¹ at the top of Enceladus' ocean (18.0 at its floor). It is 3.0 g kg⁻¹ under 20 km of ice on Europa (19.5 under 4 km). Enceladus' measured salinity, 0.5–2% (5–20 g kg⁻¹; Postberg et al. 2009), lies below S\*. Even its poles are stratified under a steady ice shell (Ames et al. 2025). Europa's ocean may be only weakly stratified (Ashkenazy & Tziperman 2021; Kang 2022), so I give it a log-uniform salinity of 1–100 g kg⁻¹ (an assumption). Counting stratified draws as losses gives Table 3.

**Table 3.** Effect of the stratified lid.

| | No lid | With salinity | Range across choices tested |
|---|---|---|---|
| Enceladus | 22% | 0% (stratified in 100% of draws) | 0–8% (S\* = 18, or salinity allowed up to 40–100 g kg⁻¹) |
| Europa | 29% | 22% (stratified in 24%) | 10–29% (S\* from 1.5 to 19.5; salinity floor 1 or 3 g kg⁻¹) |

For Enceladus, S\* = 18 gives 1.6%, a salinity range of 5–40 g kg⁻¹ gives 7.4%, and 1–100 g kg⁻¹ gives 7.7%. For Europa, S\* = 1.5 gives 26.0%, S\* = 19.5 gives 10.1%, and a salinity range of 3–100 g kg⁻¹ gives 28.5%.

**Do holes in the lid help?** Earth's thermocline opens at the poles, where deep water forms. In an axisymmetric hemisphere, a polar gap frees the wall only if its half-width exceeds θ_c, with

$$\tan\theta_c = \frac{h}{R}\sqrt{\frac{D_{\rm side}}{D_{\rm lid}}}, \tag{5}$$

where h is the depth below the lid and R the radius. I measured θ_c between 25° and 30° for a lid with 10× slower mixing (predicted 27°), and between 55° and 70° for a 100× slower lid (predicted 58°). For real oceans this gives 60–87° for Earth, ≥70° for Europa and ≥80° for Enceladus. These are far wider than any real polar convection zone. Mixing alone therefore pins Earth too.

**Do currents carry the wall?** I added one overturning cell that sinks at the pole and rises through the lid at low latitudes. The lid Péclet number Pe = W·(lid thickness)/D_lid compares the current across the lid with mixing across it. For lids 10× and 30× slower, the wall stays pinned at Pe 0.04 and 0.4 and is freed at Pe ≳ 1 (Pe 1.2, 4 and 12). Halving the time step changed the takeover time by 4%. A lid 100× slower stayed pinned up to Pe 80, the strongest current the grid resolves. Water flows through the pinned wall and changes hand as it crosses, so the pin acts like static friction. Once freed, the ocean went to the hand holding more water (the deep hand, ~62%).

For early Earth, the thermocline mixes at about 10⁻⁵ m² s⁻¹ (Ledwell et al. 1993) against about 10⁻⁴ in the abyss (Munk 1966; Waterhouse et al. 2014), a ratio of ~10. Munk's upwelling of ~10⁻⁷ m s⁻¹ across a ~1 km thermocline gives Pe ~5–10, past the freeing threshold. Enceladus' stratified lid sits over a convecting interior, a contrast of ≥1000×. That is far beyond the 100× lid that held up to Pe 80, so the Enceladus result (~0%) stands.

**Closer to a real ocean.** Three changes, each predicted before it ran:

1. Thinner walls grip harder. Halving the wall width relative to the lid raised the current needed to free it from Pe 0.4–1.2 to Pe 1.2–4.
2. A fast-mixed surface layer above the thermocline raised it to Pe 12–40. With a closed lid, Earth's Pe of 5–10 would not be enough.
3. Earth's thermocline is not closed: it reaches the surface near the poles, where deep water sinks. With a 20° polar opening, a current sinking at the pole freed the wall at Pe 1.2, also with walls two and four times thinner. The takeover times were 580, 650 and 740 diffusion-time units, and the top hand took the whole ocean.

So in this model Earth's result holds because its overturning sinks where the thermocline opens. These runs are idealized: one flow shape, axisymmetric, a radius only 4× the depth, and walls about 0.35× the lid thickness. Even the thinnest walls run are still about 1–12× too thick relative to a real thermocline, so real pinning may be stronger.

### 3.4. The sweep for Enceladus, Europa and early Earth

Before the lid, the weak force wins in 22% of plausible cases for Enceladus, 29% for Europa and 75% for early Earth's open ocean. If layers never merged, the moons would fall to 3% and 7%. Vertical mixing is the most decisive input. In the low half of its range the win rates are 0%, 0% and 52% (Enceladus, Europa, Earth); in the high half, 43%, 57% and 97%. The next inputs, measured as the win rate in the upper half of each range minus the lower half, are concentration (+21 and +10 points for the moons), time (+15 and +11; +32 for Earth), bias g (+8 and +4) and rate k₂ (+5). Faster sideways mixing lowers the win rate slightly (−3, −2 and −9). F and volume change the moons by 2 points or less.

**What would the unknowns have to be?** With every other input still drawn over its range, a moon passes 50% only through vertical mixing: D_z ≥ 2.7×10⁻⁵ m² s⁻¹ for Enceladus and ≥ 6.1×10⁻⁶ m² s⁻¹ for Europa. No slice of any other input's range takes either moon above 38%. Even in the fastest tenth of vertical mixing, the moons level off at 66% and 86%. Early Earth stays above 50% unless D_z < 2.9×10⁻⁶ m² s⁻¹ or the time is under 1.4×10⁷ yr. It wins 87% near the measured 1.1×10⁻⁵ m² s⁻¹ of Ledwell et al. (1993).

**Assumed rules.** Two rules were assumed rather than measured: the deciding top layer's share, and the Δ ≥ 2 bar. Scaling the top layer's share by 0.1, 0.28 or 2.1, letting the whole ocean decide, or moving the bar to Δ ≥ 1 or Δ ≥ 3 changes a moon by at most 2.2 points and never changes the order.

**Bias size.** The bias range used lies at or below the low end of current calculations (Quack et al. 2022). A bias 10× larger raises Enceladus to 25.7% and Europa to 30.7%. A bias 100× larger gives 27.8% and 31.2%. Even an unlimited bias stops at 28.6% and 31.2%, because the reaction must still finish and the layers must still merge. Earth is 74.7% in every case.

**Swinging conditions.** Tides, orbital cycles and climate make the amplifier's drive swing while it ramps. A swing only weakens the effective bias. In the full sweep, a mild swing gives 21%/28%/75% (Enceladus/Europa/Earth), a plausible one 20%/26%/75%, and a harsh one 17%/24%/75%.

### 3.5. Ranking the Solar System's ocean worlds

I ran the same sweep and lid test for every moon or dwarf planet thought to have, or to have had, a liquid-water ocean. Only five properties change between worlds: body size, ocean depth, ice thickness, age and (where a lid is possible) salinity. Every world uses Europa's chemistry and mixing ranges. This is an assumption, because no measurement sets these worlds apart on those inputs. Each world except Enceladus and Europa uses V = 4πR²H (times the same U(0.7, 1.3) factor) and L = πR, with R the body radius. Where a lid is possible and no salinity is measured, salinity is log-uniform 1–100 g kg⁻¹ (an assumption). Table 4 gives the inputs and results, and Figure 4 shows them.

**Table 4.** Ocean worlds: inputs and the share of 200,000 draws in which the weak force wins. "Depth sweep" re-runs the world at each end of the depth range; the shallower ocean always scores higher. Age is 10⁸–4×10⁹ yr unless noted.

| World | Ocean status | R (km) | H (km) | S\* (g kg⁻¹) | No lid | Pinned | Final | Depth sweep |
|---|---|---|---|---|---|---|---|---|
| Umbriel | possible (Castillo-Rogez et al. 2023) | 585 | 20 | none | 38.7% | 0% | 38.7% | 41% (10 km) / 37% (30 km) |
| Ariel | possible (Castillo-Rogez et al. 2023) | 579 | 30 | none | 37.0% | 0% | 37.0% | 44% (5 km) / 37% (30 km) |
| Oberon | possible (Castillo-Rogez et al. 2023) | 761 | 30 | none | 36.5% | 0% | 36.5% | 41% (10 km) / 34% (50 km) |
| Titania | possible (Castillo-Rogez et al. 2023) | 789 | 30 | none | 36.4% | 0% | 36.4% | 41% (10 km) / 34% (50 km) |
| Pluto | disputed (Johnson et al. 2016) | 1188 | 100 | none | 30.2% | 0% | 30.2% | 33% (50 km) / 27% (200 km) |
| Callisto | likely (Vance et al. 2018; Cascioli et al. 2026) | 2410 | 100 | none | 28.5% | 0% | 28.5% | 38% (20 km) / 22% (300 km) |
| Ganymede | confirmed (Saur et al. 2015) | 2631 | 100 | none | 28.3% | 0% | 28.3% | 33% (50 km) / 22% (300 km) |
| Triton | likely (Hansen et al. 2021) | 1353 | 150 | none | 27.7% | 0% | 27.7% | 36% (30 km) / 26% (200 km) |
| Titan | disputed (Petricca et al. 2025); depth from Sohl et al. (2014) | 2575 | 200 | none | 24.3% | 0% | 24.3% | 33% (50 km) / 20% (400 km) |
| Europa | confirmed | 1460 | Table 1 | 3.0 | 28.6% | 24% | 21.7% | — |
| Dione | possible (Beuthe et al. 2016) | 561 | 65 | 5.6 | 33.6% | 38% | 21.0% | 23% (30 km) / 20% (100 km) |
| Miranda | past only, 100–500 Myr ago (Strom et al. 2024); age 10⁸–5×10⁸ yr | 236 | 100 | 21.9 | 26.9% | 67% | 8.9% | 9% (100 km) / 8% (150 km) |
| Mimas | likely, < 25 Myr old (Lainey et al. 2024; Denton & Rhoden 2025); age 2×10⁶–2.5×10⁷ yr | 198 | 45 | 22.6 | 14.7% | 68% | 4.7% | 5% (35 km) / 4% (55 km) |
| Enceladus | confirmed | 232 | Table 1 | 20 | 21.8% | 100% | 0.0% | — |

Two things set the number. Where the ice is thin and the water only mildly salty, a stratified lid can form, and it pins the wall (Europa, Enceladus, Dione, Mimas, Miranda). Without a lid, the ocean first freezes into stacked layers of opposite hand that must merge, and thinner oceans merge faster. Every lid-free world falls between 24% and 39%. Ganymede's ice exerts about 200 MPa, which rules out a top lid; Pluto has about 90 km of ice, and the Uranian moons 184–222 km. The higher numbers for Uranus' moons rest on oceans that are model predictions only, and those models disagree on how thick the oceans are today (Castillo-Rogez et al. 2023; Tyler 2025). No percentage here includes the chance that the ocean exists.

**Against the newest measurements.** The Europa lid result assumed 20 km of ice. Juno's microwave radiometer measured a 29 ± 10 km shell (Levin et al. 2025). Drawing the thickness from that measurement, the ice is thinner than 23 km in 27% of draws, Europa is pinned in 9.1% of draws, and it scores 26.1%. That lies inside Ganymede's own depth spread (22–33%), so Ganymede and Europa tie. For Enceladus, the 5–20 g kg⁻¹ salinity gives 0.0% at both S\* = 20 and S\* = 22.5. Higher salinity estimates collected by Kang et al. (2022a) give 6.0% and 4.2% for 7–30 g kg⁻¹, and 17.7% and 14.7% for 17–40 g kg⁻¹. Enceladus is last among confirmed oceans under every estimate. A low-salinity stratified top would mix in hundreds of years or more (Zeng & Jansen 2021), but silica nanoparticles suggest that the Enceladus ocean mixes within years (Hsu et al. 2015). If so, its salinity is high and its score sits near the upper end.

Of the three confirmed oceans, then, Ganymede (28%) and Europa (22–26%) are the most favourable, and Enceladus (0–18%) is the least favourable.

### 3.6. Early Earth and its magnetic field

Early Earth scores 75%, well above the icy moons. Its ocean is shallow relative to its width (volume 1.33×10¹⁸ m³, mean depth 3.7 km), and its horizontal eddies (~10³ m² s⁻¹; Abernathey & Marshall 2013) and even its vertical mixing clear the one-ocean bar.

Early Earth had one chiral push that the icy moons lack. Magnetite magnetized by the planet's field gives an RNA precursor a ~60% excess (Ozturk et al. 2023a). The field points opposite ways in the two hemispheres, however, so the two halves get opposite hands. With the weak force alone, the favoured hand wins 75%, the other 0%, and 25% of draws end with no single hand. Adding the magnetic field gives 35%/35%/31%. The field pins two hands in 7% of draws. Where it acts it outweighs the weak force, and because its polarity is random it picks L or D equally often. Adding a meteorite L-excess (with the weak force favouring D) gives 66% L, 4% D and 31% no single hand. Meteorites fell on every Solar System body, however, so that push is not unique to Earth.

The weak force keeps deciding Earth's hand only if the share of deciding molecules made on magnetized magnetite, f_mag, is below 1.9×10⁻¹⁵. With field reversals every 10⁵–10⁸ yr this is 1.7–2.7×10⁻¹⁵; with reversals every 10⁴ yr it is 5.8×10⁻¹⁵. For the icy moons, if a moon's own field aligned magnetite (Ozturk et al. 2023b), the worst case allows less than about 2–4 × 10⁻¹⁵ of the deciding chemistry there (Europa, Enceladus). These are properties of the model, not measurements of real magnetite.

## 4. Discussion

### 4.1. Racemic Bennu and Ryugu

Samples from asteroid Bennu, plausibly a fragment of a wet parent body, carry chiral non-protein amino acids that are racemic or nearly so (Glavin et al. 2025). Ryugu's amino acids are also racemic (Parker et al. 2023; Naraoka et al. 2023). This tests the model, because the model assumes that some amplifying chemistry exists. With assumed ranges that bracket published values, c·τ on Bennu's parent body was 3×10⁸ to 9×10¹¹ M·s. Any amplifier with k₂ ≥ 10⁻¹¹ to 3×10⁻⁸ M⁻¹ s⁻¹ would have finished there. Taken at face value, a racemic Bennu cuts the weak force's win rate to 2% (Enceladus), 4% (Europa) and 15% (Earth).

But that assumes an amplified excess lasts once it is made. That holds only while something keeps driving the system. At equilibrium both hands are equal, so a closed, reversible network must end up racemic after a temporary "chiral excursion" (Blanco et al. 2011; Stich et al. 2016; Higgs & Blackmond 2025). I reproduced this with my own network, made reversible with rates that obey detailed balance and run closed. From a 0.1% seed, the excess peaks at 70% and fades after about 7×10³ reaction times when reverse steps are half the forward ones. With reverse steps 100× slower, it peaks at 100% and fades after about 2×10⁶. The same network kept open holds 99.8% for as long as it is driven.

Bennu's parent body was closed, and its liquid water lasted only while short-lived radioactive heat did. If its excess rose and then faded, the conflict mostly disappears. Given a racemic Bennu, the win rates become 22%/29%/75% and 23%/30%/75% for the original k₂ range, and 15%/21%/57% and 17%/22%/55% for a range extended down to 10⁻¹² M⁻¹ s⁻¹. If instead the water left before the network relaxed, the excess froze in. Then 97%, 73%, 64% or 43% of Bennu-like draws are racemic, depending on the reverse rates and the k₂ range. With slow reverse steps, a racemic Bennu is about 1.5–1.7× less likely. The reverse rates, which are unmeasured, decide how much Bennu counts against the model. One point stays open: isovaline cannot racemize, so the fade rule does not explain Bennu's racemic isovaline. I report both the 2%/4%/15% and the 22%/29%/75% results as possible.

The resolution makes a prediction. Samples from a driven ocean should carry a lasting excess. Enceladus' plume H₂ points to ongoing water–rock reactions (Waite et al. 2017). Without an amplifier, models expect the Enceladus ocean to be racemic, because racemization outpaces production (Steel et al. 2017). A lasting excess in Enceladus' plume would therefore support the idea, and a racemic plume would count against it. Capillary-electrophoresis instruments being built for ocean worlds can separate mirror-image molecules in salty samples (Mora et al. 2022). A mass spectrometer alone cannot, because mirror-image molecules have the same mass.

### 4.2. Smaller and shorter-lived waters

**Connected compartments.** Sixteen compartments that exchange material faster than the decision time act as one pool of their combined volume (P = 0.724 against 0.675 when isolated; P = 0.686 and 0.714 at intermediate and fast exchange). In a vent mound the decision takes about 56 yr. A diffusion-connected mound behaves like a ~2 m³ pond (P = 0.500). A flow-connected 100 m mound gives P = 0.73, and a flow-connected cubic kilometre gives P = 1.000.

**Long runs.** Once decided, a pool keeps its hand: the favoured share stayed at 0.730 for 10⁴ time units after the decision. Flips follow Kramers' rate within 5% and fall off as e^(−N/2), so a pond or ocean never flips.

**Titan impact-melt pools.** A 1 km impactor could melt crater ice on Titan for thousands of years (Lunine et al. 1999). Running the same model on such pools, the weak force wins in 33% of plausible cases for a 1 km-impactor crater, 51% for the Selk melt pool and 74% for the Menrva melt zone. Concentration decides most. Any inherited excess above 1×10⁻¹⁵, 2×10⁻¹⁶ or 6×10⁻¹⁷ (respectively) would override the weak force.

## 5. Limitations

- All results are order-of-magnitude estimates. The real-chemistry factor F ≈ 0.7 is checked for two reaction schemes only.
- No measured prebiotic reaction yet performs the needed self-amplification. A recent prebiotic network could (Higgs & Blackmond 2025). The rate constant k₂ is therefore swept, not measured.
- Vertical mixing is unmeasured in every icy-moon ocean. It is the single largest reason the icy-moon answer is not settled (Section 3.2).
- Ocean turbulence is treated as eddy diffusion. The lid and current experiments use one idealized, axisymmetric overturning cell with walls thicker than in a real ocean. This is not a full circulation model.
- The time τ is set to the ocean's age, the most favourable case.
- The model leaves out salt-driven layering deep in large oceans (Vance et al. 2014). If such layers pinned the boundary half of the time, Ganymede would drop to 14%.
- Ammonia, which is likely in the Titan, Triton, Pluto and Uranian oceans, is missing from TEOS-10.
- Every world beyond Enceladus and Europa uses Europa's chemistry and mixing ranges.
- The PVED has never been measured, and its sign for amino acids in water is unsettled (Quack et al. 2022).
- Meteorite seeding predicts the same left-handed outcome throughout the Solar System as the weak force. A left-handed icy-moon ocean could not by itself tell the two apart.
- The Bennu reconciliation depends on unmeasured reverse rates and on Bennu's parent body behaving as a closed system (Section 4.1).

## 6. A planned laboratory test

No known prebiotic reaction meets all the requirements for self-amplifying handedness, so I have pre-registered a laboratory search for one. The smallest rate that a run can detect is

$$k_{2,\min} = \frac{\ln(1 + 3\sqrt{2}\,\sigma/e_0)}{cT}, \tag{6}$$

where σ is the per-sample noise in the enantiomeric excess, e₀ the seed excess, c the concentration and T the run length. The realistic arm uses 1 mM, a 20% seed, 3 years and 9 independent vials per condition. It reaches 1.8×10⁻⁶ M⁻¹ s⁻¹ (2.9×10⁻⁶ with 3 vials), below the pond floor of 3×10⁻⁶. The ocean floor, 3×10⁻⁹, is out of reach of any laboratory run at natural concentration, so an accelerated arm measures how the rate scales with concentration. Mirror-seeded (D) controls guard against contamination by biological L-amino acids. This experiment has not been run.

## 7. Conclusions

1. For the weak force to beat chance, about 10³³–10³⁴ molecules must decide together. Ponds and rock pores fall short by many orders of magnitude. Ocean-scale water clears the bar.
2. In an icy-moon ocean the answer is set by ocean physics. Weak vertical mixing freezes the ocean into layers of opposite hand, which merge only if D_z ≥ H·R/(2t). Published estimates straddle that bar.
3. A stratified lid under the ice pins the boundary between hands. From TEOS-10 a lid forms below S\* and cannot form under more than about 28 MPa of ice.
4. Of the confirmed oceans, Ganymede (28%) and Europa (22–26%) are the most favourable, and Enceladus (0–18%) the least. Early Earth's open ocean scores 75%.
5. The testable result is a ranking, not an absolute probability. A deep, unstratified ocean (Ganymede) or a thick-shelled Europa is a better place than Enceladus to look for a hand chosen by the weak force. A lasting excess in a driven ocean, with racemic samples from extinct parent bodies, would support the idea.

## Acknowledgments

I thank Dr. Chris McKay (NASA Ames Research Center) for discussion. [CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]

AI-usage disclosure: An AI assistant (Claude, Anthropic; [TODO: model name and version]) was used to help write the analysis code and to draft part of the manuscript text. I directed the work, and I checked every result, number and reference against the code output and the cited sources. Responsibility for the final manuscript text lies entirely with me.

**Software:** NumPy, SciPy, Matplotlib, GSW-Python (TEOS-10; S\* check only), libRoadRunner and Antimony (SBML cross-check of the reaction network only). [TODO: add software citations and version DOIs]

## Software and data availability

All results come from Python scripts in a public repository: https://github.com/helloamonkey742-code/chirality-project (MIT licence). [TODO: archive repo on Zenodo] Each script prints and asserts its own self-checks. `./run_all.sh` reruns the fast subset in 4–7 minutes, and `./run_all.sh --full` adds the slow 3D and reaction-network runs.

| Section | Scripts |
|---|---|
| 2.1, 3.1 (normal form, molecule threshold) | sim.py, review.py |
| 2.4, 3.4 (ocean sweep, robustness) | ocean.py, uncertainty.py, breakeven.py, assumed_rules.py, wobble.py, wobble_sweep.py |
| 2.3, 3.2 (patches, layers, sphere) | domains3d.py, coarsen3d.py, aniso3d.py, sphere_layers.py |
| 2.5 (reaction network) | frank.py, frank2.py, frank_sbml.py |
| 2.6, 3.3 (stratified lid) | two_zone.py, strat_sweep.py, polar_gap.py, currents.py, realistic.py, verify_sstar.py |
| 3.5 (ocean worlds) | jupiter_moons.py, ocean_worlds.py, robust_ranking.py |
| 3.6 (early Earth) | earth.py, earth_breakeven.py, field_walls.py, field_flips.py |
| 4.1 (Bennu, Ryugu) | review.py, quench.py |
| 4.2 (compartments, long runs, Titan) | review.py, impact_melt.py |
| 6 (experiment sizing) | design.py, review.py |
| Figures 1–4 | sketch.py, ocean.py, coarsen_fig.py, make_rnaas_fig.py |

## Figures

**Figure 1** (sketch.png). The landscape that the excess α moves in. The weak force tilts it very slightly, and molecular noise shakes it. Which valley α falls into is decided in a short window around the switch-on point.

**Figure 2** (threshold.png). Minimum reactant concentration needed for the weak force to win, against the bias g, at k₂ = 10⁻³ M⁻¹ s⁻¹. Lines show a small lake (1 km² × 10 m) and the Enceladus, early-Earth and Europa oceans. The shaded band is the modern PVED range; horizontal lines mark 1 µM and 1 nM.

**Figure 3** (coarsening.png). Coarsening. Top: slices of a 96³ cube starting from a 55% majority; patches round off and shrink, and the majority reaches 99% by t = 400. Bottom: side view of a tall box, equivalent to weak vertical mixing; the patchwork freezes into flat layers of opposite hand, unchanged to t = 5000. Right: patch size grows as l = A√(Dt), with exponent 0.45.

**Figure 4** (rnaas_fig.png). Percentage of 200,000 plausible input draws in which the hand favoured by the weak force wins, for each world. Solid bars include the stratified-lid test. Hatched bars show the value if no lid formed. Colour shows the ocean status.

## References

Abernathey, R. P., & Marshall, J. 2013, Global Surface Eddy Diffusivities Derived from Satellite Altimetry, JGRC, 118, 901, doi:10.1002/jgrc.20066

Allen, S. M., & Cahn, J. W. 1979, A Microscopic Theory for Antiphase Boundary Motion and Its Application to Antiphase Domain Coarsening, AcMet, 27, 1085, doi:10.1016/0001-6160(79)90196-2

Ames, F., Ferreira, D., Czaja, A., & Masters, A. 2025, Ocean Stratification Impedes Particulate Transport to the Plumes of Enceladus, CmEE, 6, 63, doi:10.1038/s43247-025-02036-3

Ashkenazy, Y., & Tziperman, E. 2021, Dynamic Europa Ocean Shows Transient Taylor Columns and Convection Driven by Ice Melting and Salinity, NatCo, 12, 6376, doi:10.1038/s41467-021-26710-0

Beuthe, M., Rivoldini, A., & Trinh, A. 2016, Enceladus's and Dione's Floating Ice Shells Supported by Minimum Stress Isostasy, GeoRL, 43, 10088, doi:10.1002/2016GL070650

Bire, S., Kang, W., Ramadhan, A., Campin, J.-M., & Marshall, J. 2022, Exploring Ocean Circulation on Icy Moons Heated From Below, JGRE, 127, e2021JE007025, doi:10.1029/2021JE007025

Blanco, C., Stich, M., & Hochberg, D. 2011, Temporary Mirror Symmetry Breaking and Chiral Excursions in Open and Closed Systems, CPL, 505, 140, doi:10.1016/j.cplett.2011.02.032

Čadek, O., et al. 2016, Enceladus's Internal Ocean and Ice Shell Constrained from Cassini Gravity, Shape, and Libration Data, GeoRL, 43, 5653, doi:10.1002/2016GL068634

Cascioli, G., et al. 2026, Updated Interior Structure of Callisto From a Reanalysis of Galileo Data, GeoRL, 53, e2025GL120812, doi:10.1029/2025GL120812

Castillo-Rogez, J., et al. 2023, Compositions and Interior Structures of the Large Moons of Uranus and Implications for Future Spacecraft Observations, JGRE, 128, e2022JE007432, doi:10.1029/2022JE007432

Denton, C. A., & Rhoden, A. R. 2025, Herschel Impact Basin Probes the Onset of a Young Ocean on Mimas, E&PSL, 671, 119614, doi:10.1016/j.epsl.2025.119614

Frank, F. C. 1953, On Spontaneous Asymmetric Synthesis, BBA, 11, 459, doi:10.1016/0006-3002(53)90082-1

Glavin, D. P., et al. 2025, Abundant Ammonia and Nitrogen-Rich Soluble Organic Matter in Samples from Asteroid (101955) Bennu, NatAs, 9, 199, doi:10.1038/s41550-024-02472-9

Hansen, C. J., et al. 2021, Triton: Fascinating Moon, Likely Ocean World, Compelling Destination!, PSJ, 2, 137, doi:10.3847/PSJ/abffd2

Higgs, P. G., & Blackmond, D. G. 2025, Autocatalytic Symmetry Breaking and Chiral Amplification in a Feedback Network Combining Amino Acid Synthesis and Ligation, PNAS, 122, e2423683122, doi:10.1073/pnas.2423683122

Hochberg, D., Buhse, T., Micheau, J.-C., & Ribó, J. M. 2023, Chiral Selectivity vs. Noise in Spontaneous Mirror Symmetry Breaking, PCCP, 25, 31583, doi:10.1039/d3cp03311b

Hsu, H.-W., et al. 2015, Ongoing Hydrothermal Activities within Enceladus, Natur, 519, 207, doi:10.1038/nature14262

Jansen, M. F., Kang, W., Kite, E. S., & Zeng, Y. 2023, Energetic Constraints on Ocean Circulations of Icy Ocean Worlds, PSJ, 4, 117, doi:10.3847/PSJ/acda95

Johnson, B. C., Bowling, T. J., Trowbridge, A. J., & Freed, A. M. 2016, Formation of the Sputnik Planum Basin and the Thickness of Pluto's Subsurface Ocean, GeoRL, 43, 10068, doi:10.1002/2016GL070694

Kang, W. 2022, Different Ice-shell Geometries on Europa and Enceladus due to Their Different Sizes: Impacts of Ocean Heat Transport, ApJ, 934, 116, doi:10.3847/1538-4357/ac779c

Kang, W., Mittal, T., Bire, S., Campin, J.-M., & Marshall, J. 2022a, How Does Salinity Shape Ocean Circulation and Ice Geometry on Enceladus and Other Icy Satellites?, SciA, 8, eabm4665, doi:10.1126/sciadv.abm4665

Kang, W., Marshall, J., Mittal, T., & Bire, S. 2022b, Ocean Dynamics and Tracer Transport over the South Pole Geysers of Enceladus, MNRAS, 517, 3485, doi:10.1093/mnras/stac2882

Kang, W., & Zhang, Y. 2026, Subsurface Ocean Salinity and Dissipation Rate Inferred from Enceladus Ice Shell Morphology, arXiv:2603.22602 (preprint)

Kondepudi, D. K., & Nelson, G. W. 1985, Weak Neutral Currents and the Origin of Biomolecular Chirality, Natur, 314, 438, doi:10.1038/314438a0

Lainey, V., et al. 2024, A Recently Formed Ocean inside Saturn's Moon Mimas, Natur, 626, 280, doi:10.1038/s41586-023-06975-9

Ledwell, J. R., Watson, A. J., & Law, C. S. 1993, Evidence for Slow Mixing across the Pycnocline from an Open-Ocean Tracer-Release Experiment, Natur, 364, 701, doi:10.1038/364701a0

Lente, G. 2007, The Effect of Parity Violation on Kinetic Models of Enantioselective Autocatalysis, PCCP, 9, 6134, doi:10.1039/b711546f

Levin, S. M., et al. 2025, Europa's Ice Thickness and Subsurface Structure Characterized by the Juno Microwave Radiometer, NatAs, 10, 84, doi:10.1038/s41550-025-02718-0

Lunine, J. I., Beauchamp, J., Smith, M. A., & Nikolaev, E. 1999, The Abiotic Generation of Homochirality on Saturn's Moon Titan, in Advances in BioChirality, ed. G. Pályi, C. Zucchi, & L. Caglioti (Amsterdam: Elsevier), 257, doi:10.1016/B978-008043404-9/50017-6

McDougall, T. J., & Barker, P. M. 2011, Getting Started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, SCOR/IAPSO WG127

Mora, M. F., Kok, M. G. M., Noell, A., & Willis, P. A. 2022, Detection of Biosignatures by Capillary Electrophoresis Mass Spectrometry in the Presence of Salts Relevant to Ocean Worlds Missions, AsBio, 22, 914, doi:10.1089/ast.2021.0091

Munk, W. H. 1966, Abyssal Recipes, DSRA, 13, 707, doi:10.1016/0011-7471(66)90602-4

Naraoka, H., et al. 2023, Soluble Organic Molecules in Samples of the Carbonaceous Asteroid (162173) Ryugu, Sci, 379, eabn9033, doi:10.1126/science.abn9033

Osborn, T. R. 1980, Estimates of the Local Rate of Vertical Diffusion from Dissipation Measurements, JPO, 10, 83, doi:10.1175/1520-0485(1980)010<0083:EOTLRO>2.0.CO;2

Ozturk, S. F., Liu, Z., Sutherland, J. D., & Sasselov, D. D. 2023a, Origin of Biological Homochirality by Crystallization of an RNA Precursor on a Magnetic Surface, SciA, 9, eadg8274, doi:10.1126/sciadv.adg8274

Ozturk, S. F., et al. 2023b, Chirality-Induced Avalanche Magnetization of Magnetite by an RNA Precursor, NatCo, 14, 6351, doi:10.1038/s41467-023-42130-8

Parker, E. T., et al. 2023, Extraterrestrial Amino Acids and Amines Identified in Asteroid Ryugu Samples Returned by the Hayabusa2 Mission, GeCoA, 347, 42, doi:10.1016/j.gca.2023.02.017

Petricca, F., et al. 2025, Titan's Strong Tidal Dissipation Precludes a Subsurface Ocean, Natur, 648, 556, doi:10.1038/s41586-025-09818-x

Postberg, F., et al. 2009, Sodium Salts in E-ring Ice Grains from an Ocean below the Surface of Enceladus, Natur, 459, 1098, doi:10.1038/nature08046

Quack, M., Seyfang, G., & Wichmann, G. 2022, Perspectives on Parity Violation in Chiral Molecules: Theory, Spectroscopic Experiment and Biomolecular Homochirality, ChSci, 13, 10598, doi:10.1039/d2sc01323a

Saur, J., et al. 2015, The Search for a Subsurface Ocean in Ganymede with Hubble Space Telescope Observations of Its Auroral Ovals, JGRA, 120, 1715, doi:10.1002/2014JA020778

Sohl, F., et al. 2014, Structural and Tidal Models of Titan and Inferences on Cryovolcanism, JGRE, 119, 1013, doi:10.1002/2013JE004512

Steel, E. L., Davila, A., & McKay, C. P. 2017, Abiotic and Biotic Formation of Amino Acids in the Enceladus Ocean, AsBio, 17, 862, doi:10.1089/ast.2017.1673

Stich, M., Ribó, J. M., Blackmond, D. G., & Hochberg, D. 2016, Necessary Conditions for the Emergence of Homochirality via Autocatalytic Self-Replication, JChPh, 145, 074111, doi:10.1063/1.4961021

Strom, C., Nordheim, T. A., Patthoff, D. A., & Fieber-Beyer, S. K. 2024, Constraining Ocean and Ice Shell Thickness on Miranda from Surface Geological Structures and Stress Modeling, PSJ, 5, 226, doi:10.3847/PSJ/ad77d7

Tyler, R. H. 2025, Ocean Worlds Maintained by Ocean Tidal Heat Expected on Several Uranian Moons, GeoRL, 52, e2025GL114911, doi:10.1029/2025GL114911

Vance, S., Bouffard, M., Choukroun, M., & Sotin, C. 2014, Ganymede's Internal Structure Including Thermodynamics of Magnesium Sulfate Oceans in Contact with Ice, P&SS, 96, 62, doi:10.1016/j.pss.2014.03.011

Vance, S. D., et al. 2018, Geophysical Investigations of Habitability in Ice-Covered Ocean Worlds, JGRE, 123, 180, doi:10.1002/2017JE005341

Waite, J. H., et al. 2017, Cassini Finds Molecular Hydrogen in the Enceladus Plume: Evidence for Hydrothermal Processes, Sci, 356, 155, doi:10.1126/science.aai8703

Waterhouse, A. F., et al. 2014, Global Patterns of Diapycnal Mixing from Measurements of the Turbulent Dissipation Rate, JPO, 44, 1854, doi:10.1175/JPO-D-13-0104.1

Zeng, Y., & Jansen, M. F. 2021, Ocean Circulation on Enceladus with a High- versus Low-salinity Ocean, PSJ, 2, 151, doi:10.3847/PSJ/ac1114

Zeng, Y., & Jansen, M. F. 2024, The Effect of Salinity on Ocean Circulation and Ice-Ocean Interaction on Enceladus, PSJ, 5, 13, doi:10.3847/PSJ/ad0cba

Zhang, Y., Kang, W., & Marshall, J. 2024, Ocean Weather Systems on Icy Moons, with Application to Enceladus, SciA, 10, eadn6857, doi:10.1126/sciadv.adn6857

Zhang, Y., et al. 2025, Long Transit Time from the Seafloor to the Ice Shell on Enceladus, MNRAS, 541, 859, doi:10.1093/mnras/staf1008
