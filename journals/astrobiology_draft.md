<!--
==========================================================================================
SUBMISSION NOTES (remove before submitting)
==========================================================================================

Checked 2026-10-11. The Sage pages block scripted fetches, so I read them in a normal browser
session. Anything not seen on an official page is marked UNCONFIRMED.

A. ASTROBIOLOGY AUTHOR INSTRUCTIONS
   Source: https://journals.sagepub.com/author-instructions/ASB
   - Header of the page: "Access: Subscription; Accepts preprints? Yes; Identity transparency: Single anonymized".
   - Article type chosen: RESEARCH ARTICLE. Limit: 7,500 words, "No longer than 15 journal pages".
     Abstract: unstructured, at most 200 words.
     The word limits do NOT count the abstract, disclosure statements, author contribution
     statements, funding information, acknowledgments, tables, figure legends or references.
   - Other types, for comparison:
       Rapid Communication: about 2,000 words, including tables and figures.
       Review: 5,000 words.
       Hypothesis Article: 3,000 words. This is a possible alternative, but it would need cutting
         to about 40% of this draft.
       Perspective: 2,500 words, 150-word abstract, at most 2 figures/tables.
       Forum: 2,000 words. Education: 3,000 words. Letter: 500 words.
   - Draft length: the main text (Sections 1-6, without the abstract, tables, figure legends,
     declarations or references) is counted in the reply that delivered this file. It is under 7,500.
   - Figures:
       No figure-count limit is stated for Research Articles beyond the 15-page cap. UNCONFIRMED.
       Formats are TIFF, PNG, EPS or JPG. Line art needs at least 600 dpi.
       ACTION: rnaas_fig.png is 300 dpi; sketch.png and coarsening.png are 150 dpi. Re-export
       the same plots at 600 dpi or more from make_rnaas_fig.py, sketch.py and coarsen_fig.py
       before submitting. Do not change their content.
   - Format:
       A format-neutral first submission is allowed. Word is preferred and LaTeX is accepted.
       Equations must be editable, so convert the LaTeX blocks below to Word equations, or
       submit in LaTeX.
   - Keywords: at least 4, listed after the abstract. This draft has 7.
   - Reference style is SAGE HARVARD (https://www.sagepub.com/reference-style-sage-harvard-msg,
     which serves a PDF guide).
       In text: (Author, year), (Author and Author, year), (Author et al., year). Separate
       multiple citations with semicolons, in alphabetical order.
       Personal communications are cited in the text only and do not go in the reference list.
       List: alphabetical, with up to three authors and then "et al.". Initials have no full
       stops. Journal names are written in full.
       DOIs are given as "DOI: 10....". Citing preprints is discouraged, so this draft cites none.
   - Declarations order:
       Acknowledgements, then Author Contributions, then "Statements and Declarations".
       The latter has these sub-headings: Ethical considerations, Consent to participate,
       Consent for publication, Declaration of conflicting interest, Funding statement and
       Data availability. Write "Not applicable" where relevant.
       A funding statement is required even with no funding.
   - Code: Sage "strongly recommend" GitHub and Zenodo links.
       ACTION (optional): archive the GitHub release on Zenodo and add the Zenodo DOI to
       the Code availability statement.
   - ORCID is required for the submitting author. ACTION: create one and add it to the title block.
   - Peer review: "To ensure the integrity of the peer review process we assign reviewers and
     cannot accept author recommendations." So do not suggest reviewers.
   - Editors: Chris McKay is listed on the editorial board as Forum Article/Deputy Editor
     (seen on the same page). See RISK 1 below.

B. FEES
   Source: https://journals.sagepub.com/author-instructions/ASB
   - Page charges are $90 per printed page, billed on acceptance.
     At the 15-page cap that is about $1,350. This is my own calculation, not a quote.
   - Colour is free online.
     Colour in print costs $800 for the first figure and $200 for each one after that.
     Choose online-only colour to avoid this.
   - Sage Choice open access is optional. The page shows no price (UNCONFIRMED).
     Quote: "Open access fees do not cover page or color charges."
   - FREE ROUTE: none found. The page charges appear to be mandatory, and the page names no
     waiver. UNCONFIRMED.
     ACTION: before submitting, ask the Managing Editor (Karen Brattain, karen.edits@gmail.com,
     the address listed on the author-instructions page) whether a waiver exists for an
     unfunded high-school author.

C. PRIOR PUBLICATION (the RNAAS note)
   Source: https://www.sagepub.com/journals/publication-ethics-policies/prior-publication
   - Sage asks authors to disclose any prior distribution to the Editor at submission and
     within the manuscript, and to cite it.
     Quote: "If a substantial portion of your manuscript has been previously published, your
     manuscript will generally not be acceptable".
     Preprints and conference abstracts are generally acceptable.
   - Research Notes of the AAS (RNAAS) is not named on the page. Whether it counts as prior
     publication is UNCONFIRMED.
     This paper reuses the RNAAS ranking and several of its sentences. The RNAAS note has
     ~900 words, against the full methods, robustness tests, mission discussion, Bennu
     analysis and experiment here.
   - ACTION:
       (1) If the RNAAS note is published, add it to the reference list and to the
           placeholder in Section 1.
       (2) Disclose it in the cover letter.
       (3) Email the Managing Editor before submitting to ask whether a published RNAAS
           note blocks a full Research Article.
       (4) If RNAAS has not yet published the note, decide which venue goes first.

D. AI POLICY
   Source: https://www.sagepub.com/journals/publication-ethics-policies/artificial-intelligence-policy
   - Assistive tools that only edit language need no disclosure.
     Quote: "AI generated (disclosure required): The primary or partial use of AI tools and/or
     LLMs that produce content such as references, text, images, or any other content that
     directly impacts the research methodology, analysis, results and/or conclusions must be
     disclosed upon submission".
     The disclosure goes "within the methods or acknowledgements".
   - Uses Sage lists as needing disclosure include "Using software or code generated by AI to
     undertake additional research" and "Code enhanced or checked for errors by AI
     assistance".
   - AI cannot be an author.
   - RISK 2: Sage lists as inappropriate "Generation of data or submission using a series of
     prompts" and "Fabricated references or falsified claims".
       All code and text here were written with heavy AI help. Disclose this fully: see the
       statement in Acknowledgements and the copy in Methods (Section 2.6).
       Every DOI below was checked against Crossref (section F).
       The submission form will also ask about AI use. Answer it the same way.

E. CONFLICTS BETWEEN SOURCES (rule followed: RNAAS_DRAFT.md wins; otherwise PAPER_DRAFT.md)
   1. Ganymede:
        RNAAS and PAPER give 28%. The logs give 28.3/28.5%.
        Used: 28%. Callisto stays at 28.5%, from PAPER (RNAAS has no value; Fig. 3 shows 28%).
   2. Magnetite share the icy moons can tolerate:
        RNAAS gives "less than about 2-4 x 10^-15". limits_lit gives ~10^-15. A log gives 1.6-3.5e-15.
        Used: the RNAAS value.
        Earth's 1.9 x 10^-15 (1.7-2.7 x 10^-15) is from PAPER section 3.13 and is a separate number.
   3. Europa S* under 20 km of ice:
        PAPER/strat_sweep gives 3.0 g/kg. VERIFICATION.md recomputes 3.2 g/kg at 24.2 MPa
        (a rounding of the ice pressure).
        Used: 3.0. The score is unchanged either way; VERIFICATION calls it a match to rounding.
   4. Coarsening figure:
        The PAPER caption says the majority "reaches 99% by t = 400". The PAPER text gives 93%
        at t = 400 from a 55% start. These may be different runs (a 96^3 box in the caption,
        a 128^3 box in the text), but this is not confirmed.
        RESOLVED (main session): not a conflict. 99% (0.989 at t = 400) is the 96³ figure run in coarsen_fig.log; 93% is the separate 128³ run in the text. Caption keeps 99%.
   5. Ozturk et al. 2023 is two papers:
        RNAAS cites Nature Communications 14:6351. PAPER cites Science Advances 9:eadg8274
        for the ~60% excess and for experiment arm D.
        Used: both, as 2023a (Sci. Adv.) and 2023b (Nat. Commun.).
   6. Bennu quench step:
        The PAPER section 3.10 text gives 16/21/57% and 17/23/57%; the step-rule table gives
        15/21/57% and 17/22/55%. PAPER says they agree within 2 points.
        Used: the table values only.
   7. Experiment arms:
        PAPER section 5 lists arms A-E. EXPERIMENT.md also has an Arm F.
        Used: A-E. The accelerated arm is 0.01 M and 0.1 M; EXPERIMENT.md's "10 mM" is the
        same value.
   8. RNAAS cites "Ghosh, in prep." That is this paper, so the citation is not reused.
   9. "We" in PAPER section 3.10 is rewritten as "I" (style only).

F. REFERENCES AND DOI CHECK
   Every DOI was fetched with `curl -s https://api.crossref.org/works/<doi>` (no mailto).
   Result: 64 checked, 64 resolved, and every title matches the source file. 0 dropped for a bad DOI.
   - McDougall & Barker 2011 (TEOS-10 manual) has no DOI and was not checked.
     The source files give no initials, so "TJ" and "PM" are marked [CONFIRM].
   - Lunine et al. 1999 is a book chapter. Crossref gives the publisher (Elsevier) and pages
     257-270 but no editors. Editors and city are marked [CONFIRM].
   - Volume/issue/page details not in the source files were taken from the same Crossref records.
   - Left out on purpose (not dropped for errors):
       Wakita et al. 2023: limits_lit names it with no full citation. The Titan impact-melt
         results (33/51/74%) that depend on it are left out too.
       Preprints (Kang & Zhang 2026; Higgins et al. 2026; Oberg et al. 2026; Jannat et al.
         2026): Sage discourages citing preprints.
       The Glavin et al. 2012 per-specimen values: PAPER marks them "[confirm ... Table 1]".
       Anything from FRANK_1953_REVIEW.md or test_against_actual_findings.py.
       The "candidate prebiotic amplifier", stress_test and lethal_stress_tests material: these
         come from untested or questionable scripts.
   - Personal communications are cited in the text only, per Sage Harvard:
       K. M. Soderlund (2026-10-06), on the two-zone ocean.
       W. Kang (2026-10-06), on the lack of a lower bound on vertical mixing.
     ACTION: confirm that both are happy to be cited this way.

G. RISKS
   1. Dr. McKay is a Deputy Editor of Astrobiology.
        Acknowledging him "for discussion" is fine only if he agrees.
        He should not handle the manuscript, so say so in the cover letter.
        He is also a co-author of a cited paper (Steel, Davila and McKay, 2017).
        [CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]
   2. AI use is heavy. See section D.
   3. The RNAAS overlap is unconfirmed. See section C.
   4. The page charges (about $1,350) have no waiver found. See section B.
   5. Figure resolution is below the 600 dpi rule. See section A.
   6. Sole high-school author with no institutional affiliation for research. This is
        allowed, but the cover letter should say so plainly.
==========================================================================================
END OF SUBMISSION NOTES
==========================================================================================
-->

# Which Ocean World Could Let the Weak Force Choose Life's Handedness?

**Article type:** Research Article

**Saptarshi Ghosh**
Cupertino High School, Cupertino, California, USA

Corresponding author: Saptarshi Ghosh, Cupertino High School, Cupertino, California, USA. Email: [ADD]. ORCID: [ADD; required for the submitting author]

**Running title:** Weak-force handedness in ocean worlds

## Abstract

Life uses one hand of its chiral molecules. The weak nuclear force makes one mirror-image form slightly lower in energy. A large, slow, self-amplifying reaction could turn this tiny bias into a one-handed result, but only if about 10³³–10³⁴ molecules decide together. Ponds and rock pores are far too small for this; ocean worlds are not. I apply a stochastic model of the process to fourteen Solar System worlds that have, or may have, liquid-water oceans, using 200,000 plausible input sets per world. Two factors set the result. First, under a thin ice shell, a mildly salty ocean can form a stable layer below the ice that traps the boundary between regions of opposite hand. Second, without this layer, a thinner ocean merges its layers of opposite hand faster. Of the confirmed oceans, Ganymede (28% of plausible cases) and Europa (22–26%) score highest, and Enceladus scores lowest (0–18%). The Uranian moons score 36–39%, but their oceans are not confirmed. Vertical mixing is the main unknown. The model predicts that a driven ocean keeps an amplified excess while a closed parent body loses it, which fits racemic Bennu samples. Testing this needs a lander or returned samples with chiral separation.

**Keywords:** Homochirality; Parity violation; Ocean worlds; Europa; Enceladus; Ganymede; Biosignatures

## 1. Introduction

Ocean worlds are leading targets in the search for habitable environments beyond Earth. Europa, Ganymede and Enceladus have confirmed liquid-water oceans below their ice. Several other moons and dwarf planets may have one (Table 2). Enceladus is the best studied. Its plume carries molecular hydrogen from ongoing water–rock reactions (Waite et al., 2017), HCN, which is a feedstock for amino acids (Peter et al., 2024), and phosphate from the ocean (Postberg et al., 2023). Missions will soon sample these worlds. Europa Clipper's dust analyser, SUDA, will look for organic molecules (Kempf et al., 2025). Dragonfly will land near Selk crater on Titan (Lorenz et al., 2021). The Enceladus Orbilander is a lander mission concept (MacKenzie et al., 2021).

One property these missions could look for is handedness. Life uses left-handed amino acids and right-handed sugars. Without an amplifier, models expect the Enceladus ocean to be racemic (equal amounts of each hand), because racemization outpaces production (Steel et al., 2017). A lasting excess of one hand in an ocean world's organics would therefore be a clean sign that something amplified it. To read such a measurement, one needs to know which oceans could amplify a bias in the first place, and what would set the hand.

One candidate for setting the hand is the weak nuclear force. It makes one mirror-image form of a chiral molecule slightly lower in energy than the other. This parity-violating energy difference (PVED) is about 10⁻¹⁷ kT, where kT is the thermal energy. Kondepudi and Nelson (1985) proposed that the PVED can select life's hand in a large, slow, self-amplifying reaction. Whether the PVED or chance chose life's hand "remains completely open" (Quack et al., 2022). The PVED has never been measured directly.

Several studies bear on this idea. Stochastic models show that a small bias selects the hand when it exceeds the noise (Hochberg et al., 2023). Lente (2007) found that the PVED cannot be amplified when self-copying starts from a few chiral molecules. My model instead assumes an open, fed reaction in which the whole population chooses at once. Lunine et al. (1999) discussed abiotic homochirality on Titan. They asked whether the PVED decides the hand during amplification, and noted that seeing such an effect may need a system far larger than a laboratory. I consider liquid-water oceans only.

Here I ask a question that mission planners can use: if the weak force did choose life's hand, which ocean world gives it the best odds? I first derive how many molecules must decide together. I then follow how the decision spreads through a slowly mixed, possibly stratified ocean, and apply the result to fourteen worlds. The headline is a ranking, not an absolute probability. I also show how the model fits the racemic amino acids returned from asteroid Bennu, and I describe a laboratory search for the missing amplifier. A short summary of the ranking was prepared as a research note [CITE THE RNAAS NOTE HERE IF PUBLISHED; OTHERWISE DELETE THIS SENTENCE].

## 2. Methods

### 2.1 The core model

Let α = (L − D)/(L + D) be the excess of one hand. Near the point where a self-amplifying reaction switches on, any such reaction reduces to the same equation (its normal form; Fig. 1):

$$d\alpha = \left(\lambda(t)\,\alpha - \alpha^3 + g\right)dt + \sqrt{\varepsilon}\,dW, \qquad \lambda = \gamma t . \tag{1}$$

Here λ is how far past the switch-on point the system is. Before it (λ < 0), the only stable state is 50/50. After it (λ > 0), there are two stable states, all-L and all-D. The bias is g = ΔE_pv/kT, where ΔE_pv is the PVED per molecule. The noise is ε ≈ 1/N, where N is the number of molecules that react together. The sweep rate is γ = 1/(kτ), where k = k₂c is the self-copying rate (rate constant k₂ times concentration c) and τ is the time the sweep takes. Time is measured in units of 1/k.

Linearising around λ = 0 gives the chance that the hand favoured by the weak force wins:

$$P = \Phi(\Delta), \qquad \Delta = \sqrt{2}\,\pi^{1/4}\, g\, \gamma^{-1/4}\, \varepsilon^{-1/2} = \sqrt{2}\,\pi^{1/4}\, g\,(k\tau)^{1/4}\sqrt{N}, \tag{2}$$

where Φ is the normal cumulative distribution. Δ = 2 means a 97.7% win. I checked the formula against 4,000 brute-force runs per point, and it agreed within ±0.005.

### 2.2 From a reaction to an ocean

For a body of water, N = c·V·f·1000·N_A. Here V is the volume in m³, f is the fraction that is well mixed on the decision time, and N_A is Avogadro's number. The weak force "wins" only if two conditions both hold. The size condition is Δ ≥ 2. The time condition is k₂cτ ≥ 10, which means the reaction actually finishes.

Real chemistry adds noise from reactions that do not copy. I checked the normal form against a full Frank-type network in a flow reactor. In this network, a feedstock makes both hands slowly, each hand copies itself, and the two hands destroy each other. Every reaction carries its own counting noise (the chemical Langevin equation). At the switch-on point, the network maps onto Eq. 1. Its extra noise lowers Δ by a factor F ≈ 0.7. A second scheme with back-reactions matches the formula within 1.7 points and gives F = 0.71 (against 0.70).

### 2.3 Space and mixing

Real oceans mix slowly, so the decision happens in many patches at once. With mixing, α becomes a field α(x, t):

$$\partial_t \alpha = \lambda\alpha - \alpha^3 + g + D_h\nabla_h^2\alpha + D_z\,\partial_z^2\alpha + \text{noise}. \tag{3}$$

D_h and D_z are the horizontal and vertical eddy diffusivities. Rescaling depth by √(D_h/D_z) turns Eq. 3 exactly into the equal-mixing equation in a box that is taller by that factor. Weak vertical mixing therefore acts like a deep, narrow ocean.

After the decision, the walls between regions of opposite hand move at D × curvature (Allen and Cahn, 1979). Closed patches shrink, so the patch size grows as l = A√(Dt). I measured growth exponents of 0.45 and 0.49 (theory 0.5) and A = 5.48 and 5.60 at two mixing strengths, so I use A ≈ 5.5. Starting from a 55% majority, the majority hand grows to 66% at t = 50, 81% at t = 200 and 93% at t = 400 (Fig. 2).

Flat walls between horizontal layers do not move. The stretched ocean is deeper than wide for 92–96% of the plausible input range. In such tall boxes, the patchwork froze into 3–5 stacked layers in 7 of 8 runs. On a round moon, a layer boundary of constant depth feels only vertical mixing. It sinks at 2·D_z/R, where R is the radius; measured sink speeds were 0.991, 0.991 and 0.984 × 2·D_z/r. So layers merge only as fast as vertical mixing allows. At D_z = 10⁻¹⁰ m²/s, a takeover would need up to 1.5×10¹² yr on Enceladus or 2.8×10¹³ yr on Europa.

All spatial runs use explicit Euler–Maruyama time steps and second-order centred differences (the 7-point Laplacian). The time step is 0.05 and the grid spacing 1, in units of 1/k and √(D/k). This satisfies the 3D stability limit D·dt/dx² ≤ 1/6.

### 2.4 The win rule and the input sweep

For each world, I draw 200,000 input sets. Each input is drawn log-uniform over its plausible range (Table 1). The inputs are the bias, rate constant, concentration, horizontal and vertical mixing, time and volume, plus the chemistry factor F. The ocean depth is H = V/(4π(L/π)²), where L = πR is the horizontal length scale. The stretched depth is H' = H·√(D_h/D_z). With Δ(f) = √2 π^¼ g (k₂cτ)^¼ (c·V·f·1000·N_A)^½, a draw wins when k₂cτ ≥ 10 and either:

- (i) H' ≤ L, (L/5.5)²/D_h ≤ τ and F·Δ(1) ≥ 2; or
- (ii) H' > L, (L/5.5)²/D_h ≤ τ, H·(L/π)/(2D_z) ≤ τ and F·Δ(min(1, L/H')) ≥ 2.

In case (ii), the top layer decides with its own molecules (a share L/H' of the ocean). This is an assumed rule, and Section 3.4 tests it. All quantities are in SI units.

### 2.5 The stratified lid

A real ocean need not be one mixing zone. It can have a convecting part with fast vertical mixing and a stratified part with slow mixing (K. M. Soderlund, personal communication, 6 October 2026). A stably stratified layer at the top of the ocean (a "lid") stops the wall between hands completely.

A one-dimensional radial model shows why. The model is a gradient flow in which a wall's energy per area grows as √D. A wall sinking into faster-mixing water raises its energy, so a slow layer on top pins it. In this model, pinning starts at a mixing ratio between 1.02 and 1.05, so any real stratification is far past it. A bias can push the wall through, but it needs h ≈ 0.20 in model units. The weak-force bias is g ≤ 4×10⁻¹⁶, about 15 orders of magnitude too weak.

A lid forms when the salinity is below a threshold S\*. Below S\*, water at its freezing point becomes lighter as it cools (Zeng and Jansen, 2021, 2024). I test for a lid with the TEOS-10 equation of state (McDougall and Barker, 2011), at the pressure of the ocean top (920 kg/m³ × surface gravity × ice thickness). S\* is 22.5 g/kg at the top of Enceladus' ocean (18.0 at its floor). It is 3.0 g/kg under 20 km of ice on Europa (19.5 under 4 km). Under more than about 28 MPa of ice, water at its freezing point sinks as it cools at any salinity, so no lid can form. On Europa, that means no lid under more than about 23 km of ice.

Holes in the lid and ocean currents could free a pinned wall. A polar gap frees the wall only if it is wider than tan θ_c = (h/R)·√(D_side/D_lid). Here h is the depth below the lid, D_side is the vertical mixing below the lid and D_lid is the mixing within it. This gives critical gaps of 60–87° for Earth, ≥70° for Europa and ≥80° for Enceladus, far wider than any real polar convection zone. An overturning current frees the wall when the lid Péclet number exceeds about 1 (Pe = current speed × lid thickness / D_lid). A lid 100× slower than the water below stayed pinned up to Pe 80, the strongest current the grid resolves. Enceladus' stratified lid sits over a convecting interior, with a contrast of ≥1000×, so its pinning stands. For early Earth, the overturning sinks where the thermocline opens near the poles. With a 20° polar opening, a current sinking at the pole freed the wall at Pe 1.2, and the top hand took the whole ocean. Earth's result therefore holds.

### 2.6 The fourteen worlds

Only five properties change between worlds: body radius, ocean depth, ice thickness, ocean age and (where a lid is possible) salinity (Table 2). Every world except Enceladus and Europa uses V = 4πR²H, times a U(0.7, 1.3) factor, with R the body radius. All worlds use Europa's chemistry and mixing ranges. This is an assumption: no measurement tells the worlds apart on those inputs. Ocean depths come from published interior models. For each world, I re-ran the sweep at both ends of a plausible depth range. Worlds with no measured salinity, where a lid can form, get a log-uniform prior of 1–100 g/kg, and a draw is pinned when its salinity is below S\*. For Enceladus, the measured salinity is 5–20 g/kg (Postberg et al., 2009), and I use S\* = 20 g/kg.

An AI assistant (Claude, Anthropic) helped write the simulation and analysis code. I directed the work and checked all results (see Acknowledgements).

## 3. Results

### 3.1 How many molecules must decide together?

Setting Δ = 2 in Eq. 2 gives the number of molecules that must react as one well-mixed pool:

$$N \ge \left(\frac{2}{\sqrt{2}\,\pi^{1/4}\, g\,(k\tau)^{1/4}}\right)^2 . \tag{4}$$

Take g = 10⁻¹⁷, the middle of the modern PVED range, and a reaction that only just finishes (kτ = 10). Then N ≥ 3.6×10³³ molecules, or 6×10⁹ mol. At 1 µM, that is a well-mixed cube of water 18 km on a side. At 1 mM, the cube is 1.8 km, and at 0.1 M it is 0.4 km. A 1 m³ pond at 1 M holds only 6×10²⁶ molecules, about 7 orders of magnitude too few.

Cold and time help only a little. The PVED is a fixed energy, so g is larger in the cold, and the needed N scales as T². It is 3.6×10³³ at 0 °C, 4.3×10³³ at 27 °C and 6.7×10³³ at 100 °C. Time enters as (kτ)^¼, so a sweep 10,000× slower lowers N by only 100×. Concentration alone cannot rescue a small setting. Only ocean-scale water, mixed within the decision time, can pool enough molecules.

Icy-moon oceans are big enough. At g = 10⁻¹⁷ and k₂ = 10⁻³ M⁻¹ s⁻¹, the minimum reactant concentration for the weak force to win is about 2 mM for a small lake. It is only about 50 pM for Enceladus, 4 pM for Earth and 1 pM for Europa. Size is therefore not the limit for icy-moon oceans. What limits them is time, when chemistry is slow, and mixing.

### 3.2 The ranking of fourteen worlds

Table 3 and Fig. 3 give the share of plausible draws in which the hand favoured by the weak force wins. Worlds without a lid score between 24% and 39%. Ocean depth sets most of this order, because a thinner ocean has fewer layers to merge. Depth uncertainty moves each value by up to about ten points (Ganymede: 22–33% for depths of 50–300 km).

Five worlds form a lid in some draws. Their scores drop as follows, with the share of draws that have a lid in brackets:

- Europa: 29% to 22% (24%)
- Dione: 34% to 21% (38%)
- Miranda: 27% to 9% (67%)
- Mimas: 15% to 5% (68%)
- Enceladus: 22% to 0% (100%)

The measured salinity of Enceladus, 5–20 g/kg (Postberg et al., 2009), is below S\* (18–22.5 g/kg under its thin shell). Across the choices I tested, its score stays at 0–8%. Even its poles are stratified under a steady ice shell (Ames et al., 2025). Europa's ocean may be only weakly stratified (Ashkenazy and Tziperman, 2021). Ganymede's ice shell exerts about 200 MPa, which rules out a top lid.

The Uranian moons score highest of all (36–39%), but their oceans are model predictions. Those models disagree on how thick the oceans are today (Castillo-Rogez et al., 2023; Tyler, 2025). Titan's ocean is disputed (Petricca et al., 2025). No percentage here includes the chance that the ocean exists.

### 3.3 Updating with the newest measurements

I then replaced two inputs with newer measurements.

- **Europa.** The lid result assumed 20 km of ice. Juno measured 29 ± 10 km (Levin et al., 2025). No lid can form under more than about 23 km. With the Juno value, Europa has a lid in only 9% of draws and scores 26%. This is within Ganymede's own spread (22–33%), so the two tie.
- **Enceladus.** Higher salinity estimates collected by Kang et al. (2022) give 4–6% for 7–30 g/kg and 15–18% for 17–40 g/kg.

Of the three confirmed oceans, Ganymede (28%) and Europa (22–26%) are the most favourable. Enceladus (0–18%) is the least favourable under every salinity estimate. Silica nanoparticles suggest that the Enceladus ocean mixes within years (Hsu et al., 2015), which may rule out a stratified top. If so, its salinity is high and it scores near the upper end of this range.

### 3.4 What controls the odds

For Enceladus and Europa without the lid test, and for early Earth's open ocean, the sweep gives 22%, 29% and 75%. If layers never merged, the two moons would fall to 3% and 7%. The merge-in-time condition holds in 33% and 34% of draws.

I wrote my predictions into the scripts before the first run, and all held. To find which input matters, I split the draws at each input's median and compared the two halves. Vertical mixing is the most decisive input. The win rate is 0%, 0% and 52% (Enceladus, Europa, Earth) in the low half of its range, against 43%, 57% and 97% in the high half. Concentration (+21 and +10 points for the two moons) and time (+15 and +11; +32 for Earth) come next, then the bias g (+8 and +4) and the rate k₂ (+5).

With every other input still drawn over its range, a moon passes 50% only through vertical mixing. That needs D_z ≥ 2.7×10⁻⁵ m²/s for Enceladus and D_z ≥ 6.1×10⁻⁶ m²/s for Europa. No slice of any other input's range takes either moon above 38%.

The rules I assumed rather than measured matter little. I scaled the deciding top layer's share, let the whole ocean decide, and moved the bar to Δ ≥ 1 or Δ ≥ 3. None of these changed a moon by more than 2.2 points, and none changed the order.

The bias range I used (≈0.01–10 aeV) lies at or below the low end of current calculations, which give "aeV to feV" for molecules of light atoms (Quack et al., 2022). A bias 10× or 100× larger raises Enceladus to 25.7% or 27.8% and Europa to 30.7% or 31.2%. Even an unlimited bias stops at 28.6% and 31.2%, because the reaction must still finish and the layers must still merge in time. Earth is 74.7% in every case.

Real oceans are not steady. Tides, orbital cycles and climate make the amplifier's drive swing while it slowly ramps. I added a periodic swing to the model and found that it only weakens the effective bias. Across the swings I tested, from mild to harsh, the icy moons lose at most 5 points and Earth is unchanged.

### 3.5 Early Earth for comparison

Early Earth's open ocean scores 75%. It has one chiral push the icy moons lack. Magnetite magnetized by the planet's field gives an RNA precursor a ~60% excess (Ozturk et al., 2023a). The field points opposite ways in the two hemispheres, however, so the two halves get opposite hands. When I added this push, the field outweighed the weak force wherever it acted. Its polarity at the time is random, so it picked L or D equally often (35% each, with 31% of draws ending without one hand).

The weak force keeps deciding Earth's hand only if fewer than about 1.9×10⁻¹⁵ of the deciding molecules form on magnetized magnetite. For the icy moons, in the worst case, the weak force keeps its lead only if less than about 2–4 × 10⁻¹⁵ of the deciding chemistry occurs on such surfaces (Europa, Enceladus; Ozturk et al., 2023a, 2023b).

## 4. Discussion

### 4.1 What missions could test

The model's testable output is a ranking. A deep ocean with no lid (Ganymede) or a thick-shelled Europa is a better place than Enceladus to look for a hand chosen by the weak force. The ranking rests on a few inputs that missions and models can measure, each listed below with what it does in the model.

- **Ice-shell thickness.** It decides whether a lid can form. On Europa, more than about 23 km of ice rules out a lid. Juno's 29 ± 10 km measurement already moved Europa from 22% to 26%.
- **Ocean salinity.** Below S\*, the lid forms and pins the wall. Enceladus' score ranges from 0% to 18% across published salinity estimates.
- **Vertical mixing.** This is the main unknown. Its range of 10⁻¹⁰–10⁻³ m²/s spans a 0% to 57% swing in the moons' win rates (Section 3.4).

**Measuring handedness.** Europa Clipper's SUDA is a time-of-flight mass spectrometer that will look for organic molecules (Kempf et al., 2025). It cannot tell the hands apart, because mirror-image molecules have the same mass. The measurement needs a lander or returned samples with a chiral separation step. Capillary-electrophoresis instruments being built for ocean worlds can separate enantiomers in salty samples (Mora et al., 2022). Radiation does not rule out such a measurement. Amino acids survive radiolysis anywhere on Enceladus' surface, while on Europa samples would need to come from about 20 cm down (Pavlov et al., 2024). A lander such as the Enceladus Orbilander concept could test the plume predictions below. Its published abstract names no chiral measurement, so a chiral separation step would have to be added (MacKenzie et al., 2021).

**Titan and Dragonfly.** My model covers liquid-water oceans only, and Titan's ocean is disputed. Dragonfly cannot reach that ocean. It does target Selk crater, where melt water mixed with organics (Barnes et al., 2021). It first lands at a dune site near Selk (Lorenz et al., 2021). Its mass spectrometer is being adapted to separate the two hands (Boulesteix et al., 2024). Dragonfly could therefore measure whether organics that formed in liquid water on Titan carry an excess.

**Predictions.**

1. A driven ocean should carry a lasting excess (Section 4.2). Samples from Enceladus' plume should show one if an amplifier works there. A racemic Enceladus would count against the idea.
2. Where a lid keeps an ocean split into regions of opposite hand, as the model predicts for Enceladus, plume sources that tap different parts of the ocean could carry opposite hands. Delivery from meteorites predicts the same excess at every source. This prediction is a hypothesis drawn from the model and has not been simulated directly.
3. A left-handed result on an icy moon would not by itself separate the weak force from seeding by meteorites, which also deliver an L-excess throughout the Solar System.

### 4.2 Reconciling the model with racemic Bennu and Ryugu

Bennu is plausibly a fragment of a wet parent body. Its chiral non-protein amino acids, including isovaline, which cannot racemize once formed, are racemic or nearly so (Glavin et al., 2025; Baczynski et al., 2026). Ryugu's amino acids are also racemic or nearly so (Parker et al., 2023). This is a direct test of the model, because the model assumes that some amplifying chemistry exists.

I assumed one chiral amino acid at 10–300 nmol/g of rock, a water-to-rock mass ratio of 0.1–1 and aqueous alteration lasting 1–10 Myr. These ranges bracket published values (Glavin and Dworkin, 2009; Glavin et al., 2025; McCoy et al., 2025). Then any amplifier with k₂ ≥ 10⁻¹¹ to 3×10⁻⁸ M⁻¹ s⁻¹ would have finished on Bennu's parent body. If an amplified excess lasted forever, Bennu would rule out the whole rate range of 10⁻⁶–1 M⁻¹ s⁻¹. If the range is extended down to 10⁻¹² M⁻¹ s⁻¹, 23% of draws keep Bennu racemic. Among those, the weak force wins in only 2% (Enceladus), 4% (Europa) and 15% (early Earth) of cases.

That conflict assumes an amplified excess lasts. It lasts only while something keeps driving the system. At equilibrium, a closed system must end up racemic. In a closed, reversible Frank network, the excess can rise for a while and then return to racemic (Blanco et al., 2011; Higgs and Blackmond, 2025). I reproduced this with my own network, made reversible with rates that obey detailed balance and run closed from a 0.1% seed:

- with each reverse rate half its forward one, the excess peaks at 70% and falls back below a tenth of its peak after about 7×10³ reaction times;
- with reverse steps 100× slower, it peaks at 100% and fades after about 2×10⁶ reaction times;
- the same network kept open (fed and flushed) holds 99.8% for as long as it is driven.

Bennu's parent body was closed, and its liquid water lasted only while short-lived radioactive heat did. Its excess therefore stays at zero if the excess never rose, or if it rose and then faded. The second case removes most of the conflict (Table 4).

Enceladus, by contrast, is still driven. Its plume H₂ points to ongoing water–rock reactions (Waite et al., 2017). Tidal heating of a porous core can keep hydrothermal activity going for tens of Myr to Gyr (Choblet et al., 2017). Serpentinization alone can supply the observed H₂ for at most about 500 Myr (Daval et al., 2022), which is still far longer than the decision times in this model.

If the water left before the network relaxed, the excess would freeze in. I read the excess off the full closed-network curve at the moment the water left, and counted |ee| < 5.2% as racemic. In those runs, the weak force's odds given a racemic Bennu matched the step rule within 2 points. With slow reverse steps, however, a racemic Bennu is about 1.5–1.7× less likely. Measuring the reverse rates therefore decides how much Bennu counts against the model.

The meteorites fit if the rule applies to excesses that stay dissolved and can racemize. In Murchison and Murray, the excesses survive only in α-methyl amino acids such as isovaline (Pizzarello and Cronin, 2000). The large L-excesses of aspartic and glutamic acid in Tagish Lake (43–59%) are in amino acids that crystallize as separate L and D crystals, which fits a crystallization model (Glavin et al., 2012; Klussmann et al., 2006). One point remains open. Isovaline cannot racemize, so the fade rule does not explain Bennu's racemic isovaline. If the closed-system argument is set aside, the weak force's odds fall back to 2%, 4% and 15%. I report both results.

### 4.3 Limitations

- All results are order-of-magnitude estimates. The factor F ≈ 0.7 was checked for only two reaction schemes. The second scheme gives its two making steps different equilibrium constants, which violates detailed balance in the open reactor. The closed-system test used one constant for both routes.
- Vertical mixing in icy-moon oceans is unmeasured across 10⁻¹⁰–10⁻³ m²/s (Zeng and Jansen, 2021). It is the single biggest reason the answer is not settled. Sideways mixing (~0.1 m²/s; Zhang et al., 2024) never limits healing. Ocean models agree that salinity sets how stratified these oceans are (Kang et al., 2022; Zeng and Jansen, 2024), but none gives a vertical diffusivity directly. There is no lower bound at high salinity, where vertical mixing can be very weak (W. Kang, personal communication, 6 October 2026). I keep the full range.
- The time τ is set to the ocean's age, the most favourable case possible.
- No measured prebiotic reaction performs the needed self-amplification yet. The Soai reaction does it in the lab but is not prebiotic, and lab tests of it need biases far larger than the PVED (Hawbaker and Blackmond, 2019). The closest prebiotic candidates are peptide-ligation reactions that amplify an excess (Deng et al., 2024) and a model network that breaks symmetry on its own (Higgs and Blackmond, 2025). The rate k₂ is therefore swept, not measured. A peptide-catalysed transamination step implies rough rates of 6×10⁻⁶–1×10⁻⁴ M⁻¹ s⁻¹ (my estimate from Yu et al., 2024). This is inside the swept range, but it is a conversion rate, not the amplification rate.
- The PVED has never been measured, and its sign for amino acids in water is unsettled (MacDermott et al., 2009; Quack et al., 2022).
- The model does not include:
  - salt-driven layering deep in large oceans (Vance et al., 2014). If it pinned the wall half of the time, Ganymede would drop to 14%;
  - ammonia, which is likely in outer-Solar-System oceans but missing from TEOS-10;
  - full ocean circulation. The lid and current tests use one idealized, axisymmetric flow, with walls far thicker relative to the lid than in a real ocean.
- Random chemical flipping between hands only adds noise. In simulation, it moves the win chance toward 50% (72% to 67% as the flipping rate rises to 0.45 of the amplifier rate).
- The Enceladus ocean volume (2.7×10¹⁶ m³) comes indirectly from ice-shell and core geometry (Čadek et al., 2016). Inputs for Callisto, Mimas and the Uranian moons come from their cited interior models (Cascioli et al., 2026; Denton and Rhoden, 2025; Lainey et al., 2024; Castillo-Rogez et al., 2023).
- All Bennu inputs are assumed, including that the same k₂ applies in both settings.

## 5. Planned laboratory test

No known prebiotic reaction meets the model's requirements. I have therefore pre-registered a wet-lab search for one. The hypothesis is that at least one prebiotic, water-based, peptide-forming system will show seed-following growth of enantiomeric excess (ee).

The smallest detectable rate is k₂,min = ln(1 + 3√2 σ/e₀)/(cT). Here σ is the ee noise per sample, e₀ the seed excess, c the concentration and T the run length. At realistic concentrations, the rate floor for a pond is 3×10⁻⁶ M⁻¹ s⁻¹ (1 mM, 100 yr). For an ocean it is 3×10⁻⁹ M⁻¹ s⁻¹ (1 µM, 10⁸ yr). No lab run can reach the ocean floor at natural concentration, so the design has two tiers:

1. **A realistic arm** at 1 mM, with a 20% seed (similar to Murchison isovaline), runs for 3 years in 9 independent vials per condition at 2.6% ee noise per vial. It reaches 1.8×10⁻⁶ M⁻¹ s⁻¹, below the pond floor.
2. **An accelerated arm** at 0.01 M and 0.1 M is plainly not natural. It exists to find any amplifier at all and to measure how the rate of ee growth scales with concentration.

There are five candidate chemistries:

- (A) amino acids with carbonyl sulfide (Leman et al., 2004);
- (B) amino acids with hydroxy acids under wet–dry cycling (Forsythe et al., 2015);
- (C) cysteine peptides that catalyse their own joining (Foden et al., 2020);
- (D) an RNA precursor on magnetite (Ozturk et al., 2023a);
- (E) peptide-catalysed transamination coupled to peptide ligation, run fed rather than closed (Deng et al., 2024; Higgs and Blackmond, 2025; Yu et al., 2024).

Arms E and C are the priority. A Viedma (2005) grinding positive control checks that the pipeline can detect real amplification. Each arm gets +L, +D (mirror) and racemic seedings, plus three negative controls: no activator, sterile-filtered, and seed alone. Analysts are blinded. An arm counts as an amplifier only if four conditions hold:

- both seeded sets rise by more than the detection threshold;
- the L- and D-seeded sets change by equal and opposite amounts;
- the racemic vials stay at zero;
- the controls stay flat.

The D-seed arm guards against contamination by biological L-amino acids. As of this draft, this experiment has not been run.

## 6. Conclusions

- For the weak force to beat chance, about 10³³–10³⁴ molecules must decide together. Ponds and rock pores fall short by many orders of magnitude. Icy-moon oceans and early Earth's ocean clear the bar.
- Among confirmed oceans, Ganymede (28%) and Europa (22–26%) give the weak force its best odds, and Enceladus (0–18%) its worst. A stratified lid under thin ice explains most of the gap. The Uranian moons score highest (36–39%), but their oceans are unconfirmed.
- Vertical mixing is the main unknown. No other input moves a moon above 38%.
- An amplified excess fades in a closed rock and holds in a driven ocean, which fits racemic Bennu. The reverse reaction rates decide how strongly Bennu counts against the model.
- The decisive measurement is chiral: a lasting excess in plume or near-surface amino acids, measured by a lander or sample return with chiral separation. A racemic Enceladus would count against the idea. A bench-top search for the missing amplifier is designed and ready to run.

## Figure legends

**Figure 1.** (file: `sketch.png`) The landscape the excess α moves in. The weak force tilts it very slightly; molecular noise shakes it. Which valley α falls into is decided in a short window around the switch-on point.

**Figure 2.** (file: `coarsening.png`) Coarsening. Top: slices of a 96³ cube starting from a 55% majority; patches round off and shrink, and the majority reaches 99% by t = 400. Bottom: side view of a tall box (equivalent to weak vertical mixing); the patchwork freezes into flat layers of opposite hands, unchanged to t = 5000. Right: patch size grows as l = 5.5√(Dt) (exponent 0.45).

**Figure 3.** (file: `rnaas_fig.png`) Percentage of 200,000 plausible input draws in which the hand favoured by the weak force wins, for each world. Solid bars include the lid test. Hatched bars show the value if no lid formed. Colour shows the ocean status.

## Tables

**Table 1.** Sweep inputs. Every range is log-uniform unless noted. The other worlds use Europa's ranges for c, D_h, D_z, g, k₂ and F (Table 2 gives their size, depth and age).

| Input | Enceladus | Europa | Early Earth ocean |
|---|---|---|---|
| bias g = PVED/kT | 4×10⁻¹⁹ – 4×10⁻¹⁶ | same | same |
| rate k₂ (M⁻¹ s⁻¹) | 10⁻⁶ – 1 | same | same |
| concentration c (M) | 10⁻¹² – 10⁻⁴ | 10⁻¹² – 10⁻⁴ | 10⁻⁷ – 3×10⁻⁴ |
| horizontal D_h (m²/s) | 10⁻² – 10² | 10⁻² – 10² | 10² – 10⁴ |
| vertical D_z (m²/s) | 10⁻¹⁰ – 10⁻³ | 10⁻¹⁰ – 10⁻³ | 10⁻⁶ – 10⁻⁴ |
| time τ (yr) | 10⁶ – 10⁹ | 10⁸ – 4×10⁹ | 10⁷ – 5×10⁸ |
| volume V (m³) | 2.7×10¹⁶, ×U(0.7, 1.3) | 3.3×10¹⁸, ×U(0.7, 1.3) | 1.332×10¹⁸, ×U(0.7, 1.3) |
| length L = π·R (R in km) | R = 232 | R = 1460 | R = 6371 |
| chemistry factor F | U(0.3, 0.7) | same | same |
| salinity (g/kg), S\* | 5–20, S\* = 20 | 1–100, S\* = 3.0 | no lid |

**Table 2.** Inputs for the other twelve worlds. The depth sweep re-runs each world with H at each end of the range. S\* is the lid threshold at the ocean-top pressure ("none": no lid possible).

| World | R (km) | H (km) | depth sweep (km) | age (yr) | S\* (g/kg) |
|---|---|---|---|---|---|
| Ganymede | 2631 | 100 | 50, 300 | 10⁸ – 4×10⁹ | none |
| Titan | 2575 | 200 | 50, 400 | 10⁸ – 4×10⁹ | none |
| Callisto | 2410 | 100 | 20, 300 | 10⁸ – 4×10⁹ | none |
| Triton | 1353 | 150 | 30, 200 | 10⁸ – 4×10⁹ | none |
| Pluto | 1188 | 100 | 50, 200 | 10⁸ – 4×10⁹ | none |
| Titania | 789 | 30 | 10, 50 | 10⁸ – 4×10⁹ | none |
| Oberon | 761 | 30 | 10, 50 | 10⁸ – 4×10⁹ | none |
| Umbriel | 585 | 20 | 10, 30 | 10⁸ – 4×10⁹ | none |
| Ariel | 579 | 30 | 5, 30 | 10⁸ – 4×10⁹ | none |
| Dione | 561 | 65 | 30, 100 | 10⁸ – 4×10⁹ | 5.6 |
| Miranda | 236 | 100 | 100, 150 | 10⁸ – 5×10⁸ | 21.9 |
| Mimas | 198 | 45 | 35, 55 | 2×10⁶ – 2.5×10⁷ | 22.6 |

**Table 3.** Share of plausible draws in which the hand favoured by the weak force wins. No value includes the chance that the ocean exists.

| World | ocean | ocean depth (km) | stratified top | weak force wins | across ocean depths |
|---|---|---|---|---|---|
| Ganymede | confirmed | 100 (Saur et al., 2015) | no (~200 MPa) | **28%** | 22–33% (50–300 km) |
| Europa | confirmed | Table 1 | 24% of draws (9% with Juno ice thickness) | **22%** (26% with Juno) | 10–29% across lid choices |
| Enceladus | confirmed | Table 1 | 100% of draws | **0%** | 0–18% across salinity estimates |
| Callisto | likely | 100 (Vance et al., 2018; Cascioli et al., 2026) | no | 28.5% | 22–38% (20–300 km) |
| Triton | likely | 150 (Hansen et al., 2021) | no | 28% | 26–36% (30–200 km) |
| Mimas | likely, < 25 Myr old | 45 (Lainey et al., 2024; Denton and Rhoden, 2025) | 68% of draws | 5% | 4–5% (35–55 km) |
| Titan | disputed (Petricca et al., 2025) | 200 (Sohl et al., 2014) | no | 24% | 20–33% (50–400 km) |
| Pluto | disputed | 100 (Johnson et al., 2016) | no (~90 km ice) | 30% | 27–33% (50–200 km) |
| Umbriel, Ariel, Oberon, Titania | possible (Castillo-Rogez et al., 2023) | 20–30 | no (184–222 km ice) | 36–39% | 34–44% (5–50 km) |
| Dione | possible (Beuthe et al., 2016) | 65 | 38% of draws | 21% | 20–23% (30–100 km) |
| Miranda | past only, 100–500 Myr ago (Strom et al., 2024) | 100 | 67% of draws | 9% | 8–9% (100–150 km) |

**Table 4.** The weak force's win rate given a racemic Bennu, when an amplified excess lasts only for a limited number of reaction times in a closed parent body.

| k₂ range (M⁻¹ s⁻¹) | excess lasts (reaction times) | Enceladus / Europa / Early Earth |
|---|---|---|
| 10⁻⁶ – 1 | 8 to 7×10³ | 22% / 29% / 75% |
| 10⁻⁶ – 1 | 0.6 to 2×10⁶ | 23% / 30% / 75% |
| 10⁻¹² – 1 | 8 to 7×10³ | 15% / 21% / 57% |
| 10⁻¹² – 1 | 0.6 to 2×10⁶ | 17% / 22% / 55% |

## Acknowledgements

I thank Dr. Chris McKay (NASA Ames Research Center) for discussion. [CONFIRM WITH DR. MCKAY BEFORE SUBMISSION]

Use of AI: An AI assistant (Claude, Anthropic) helped write the code and text of this study. That includes the simulation and analysis scripts, the drafting of the manuscript, and the search for and summary of the literature. I directed the work and checked all results and references. I checked every DOI in the reference list against its Crossref record. The AI assistant is not an author.

## Author Contributions

S.G. is the sole author. He designed the study, wrote and ran the code with the AI assistance described above, analysed the results, and wrote and revised the manuscript.

## Statements and Declarations

**Ethical considerations:** Not applicable. This study used no human participants, animals or human data.

**Consent to participate:** Not applicable.

**Consent for publication:** Not applicable.

**Declaration of conflicting interest:** The author declared no potential conflicts of interest with respect to the research, authorship, and/or publication of this article.

**Funding statement:** The author received no financial support for the research, authorship, and/or publication of this article.

**Data availability:** All results are produced by the code below. No external data sets were used beyond the published values cited in the text.

**Code availability:** All code is public at https://github.com/helloamonkey742-code/chirality-project (MIT licence). Each script prints and asserts its own self-checks. `./run_all.sh` runs a fast subset in 4–7 minutes, and `./run_all.sh --full` adds the slow 3D and real-chemistry runs. [Optional: add a Zenodo archive DOI.] The table maps the sections to the scripts.

| Section | Scripts |
|---|---|
| 2.1 Core model, Eq. 2 check | `sim.py` |
| 2.2 Chemistry factor F | `frank.py`, `frank2.py` |
| 2.3 Coarsening, layers, sphere sinking, numerics | `coarsen3d.py`, `k3.py`, `aniso3d.py`, `sphere_layers.py`, `shell_res.py` |
| 2.4, 3.4 Input sweep, win rule, robustness | `uncertainty.py`, `breakeven.py`, `assumed_rules.py`, `wobble.py`, `wobble_sweep.py` |
| 2.5 Stratified lid, polar gaps, currents | `two_zone.py`, `strat_sweep.py`, `polar_gap.py`, `currents.py`, `realistic.py` |
| 2.6, 3.2, 3.3 Fourteen worlds, newest measurements | `jupiter_moons.py`, `ocean_worlds.py`, `robust_ranking.py` |
| 3.1 Molecule threshold, minimum concentrations | `review.py` (R1), `ocean.py` |
| 3.5 Earth's magnetic field | `earth.py`, `field_walls.py`, `field_flips.py`, `earth_breakeven.py` |
| 4.2 Bennu, closed-system fade, freezing | `review.py` (R2, R6), `quench.py` |
| 4.3 Random flipping between hands | `review.py` (R7) |
| 5 Experiment sizing | `review.py` (R5), `design.py` |
| Figures 1–3 | `sketch.py`, `coarsen_fig.py`, `make_rnaas_fig.py` |

## References

Allen SM and Cahn JW (1979) A microscopic theory for antiphase boundary motion and its application to antiphase domain coarsening. *Acta Metallurgica* 27(6): 1085–1095. DOI: 10.1016/0001-6160(79)90196-2.

Ames F, Ferreira D, Czaja A, et al. (2025) Ocean stratification impedes particulate transport to the plumes of Enceladus. *Communications Earth & Environment* 6: 63. DOI: 10.1038/s43247-025-02036-3.

Ashkenazy Y and Tziperman E (2021) Dynamic Europa ocean shows transient Taylor columns and convection driven by ice melting and salinity. *Nature Communications* 12: 6376. DOI: 10.1038/s41467-021-26710-0.

Baczynski AA, Mcintosh OM, Simkus DN, et al. (2026) Multiple formation pathways for amino acids in the early Solar System based on carbon and nitrogen isotopes in asteroid Bennu samples. *Proceedings of the National Academy of Sciences* 123(8): e2517723123. DOI: 10.1073/pnas.2517723123.

Barnes JW, Turtle EP, Trainer MG, et al. (2021) Science goals and objectives for the Dragonfly Titan rotorcraft relocatable lander. *The Planetary Science Journal* 2(4): 130. DOI: 10.3847/psj/abfdcf.

Beuthe M, Rivoldini A and Trinh A (2016) Enceladus's and Dione's floating ice shells supported by minimum stress isostasy. *Geophysical Research Letters* 43(19): 10088. DOI: 10.1002/2016GL070650.

Blanco C, Stich M and Hochberg D (2011) Temporary mirror symmetry breaking and chiral excursions in open and closed systems. *Chemical Physics Letters* 505(4–6): 140–147. DOI: 10.1016/j.cplett.2011.02.032.

Boulesteix D, Buch A, Samson J, et al. (2024) Dimethylformamide dimethyl acetal reagent for in situ chiral analyses of organic molecules on Titan with the Dragonfly mass spectrometer space instrument (Dragonfly mission). *Journal of Chromatography A* 1722: 464860. DOI: 10.1016/j.chroma.2024.464860.

Čadek O, Tobie G, Van Hoolst T, et al. (2016) Enceladus's internal ocean and ice shell constrained from Cassini gravity, shape, and libration data. *Geophysical Research Letters* 43(11): 5653–5660. DOI: 10.1002/2016GL068634.

Cascioli G, Petricca F, Mazarico E, et al. (2026) Updated interior structure of Callisto from a reanalysis of Galileo data. *Geophysical Research Letters* 53(15): e2025GL120812. DOI: 10.1029/2025GL120812.

Castillo-Rogez J, Weiss B, Beddingfield C, et al. (2023) Compositions and interior structures of the large moons of Uranus and implications for future spacecraft observations. *Journal of Geophysical Research: Planets* 128(1): e2022JE007432. DOI: 10.1029/2022JE007432.

Choblet G, Tobie G, Sotin C, et al. (2017) Powering prolonged hydrothermal activity inside Enceladus. *Nature Astronomy* 1(12): 841–847. DOI: 10.1038/s41550-017-0289-8.

Daval D, Choblet G, Sotin C, et al. (2022) Theoretical considerations on the characteristic timescales of hydrogen generation by serpentinization reactions on Enceladus. *Journal of Geophysical Research: Planets* 127(2): e2021JE006995. DOI: 10.1029/2021JE006995.

Deng M, Yu J and Blackmond DG (2024) Symmetry breaking and chiral amplification in prebiotic ligation reactions. *Nature* 626(8001): 1019–1024. DOI: 10.1038/s41586-024-07059-y.

Denton C and Rhoden A (2025) Herschel impact basin probes the onset of a young ocean on Mimas. *Earth and Planetary Science Letters* 671: 119614. DOI: 10.1016/j.epsl.2025.119614.

Foden CS, Islam S, Fernández-García C, et al. (2020) Prebiotic synthesis of cysteine peptides that catalyze peptide ligation in neutral water. *Science* 370(6518): 865–869. DOI: 10.1126/science.abd5680.

Forsythe JG, Yu S, Mamajanov I, et al. (2015) Ester-mediated amide bond formation driven by wet–dry cycles: a possible path to polypeptides on the prebiotic Earth. *Angewandte Chemie International Edition* 54(34): 9871–9875. DOI: 10.1002/anie.201503792.

Glavin DP and Dworkin JP (2009) Enrichment of the amino acid L-isovaline by aqueous alteration on CI and CM meteorite parent bodies. *Proceedings of the National Academy of Sciences* 106(14): 5487–5492. DOI: 10.1073/pnas.0811618106.

Glavin DP, Dworkin JP, Alexander CMO, et al. (2025) Abundant ammonia and nitrogen-rich soluble organic matter in samples from asteroid (101955) Bennu. *Nature Astronomy* 9(2): 199–210. DOI: 10.1038/s41550-024-02472-9.

Glavin DP, Elsila JE, Burton AS, et al. (2012) Unusual nonterrestrial L-proteinogenic amino acid excesses in the Tagish Lake meteorite. *Meteoritics & Planetary Science* 47(8): 1347–1364. DOI: 10.1111/j.1945-5100.2012.01400.x.

Hansen CJ, Castillo-Rogez J, Grundy W, et al. (2021) Triton: fascinating moon, likely ocean world, compelling destination! *The Planetary Science Journal* 2(4): 137. DOI: 10.3847/PSJ/abffd2.

Hawbaker NA and Blackmond DG (2019) Energy threshold for chiral symmetry breaking in molecular self-replication. *Nature Chemistry* 11(10): 957–962. DOI: 10.1038/s41557-019-0321-y.

Higgs PG and Blackmond DG (2025) Autocatalytic symmetry breaking and chiral amplification in a feedback network combining amino acid synthesis and ligation. *Proceedings of the National Academy of Sciences* 122(20): e2423683122. DOI: 10.1073/pnas.2423683122.

Hochberg D, Buhse T, Micheau JC, et al. (2023) Chiral selectivity vs. noise in spontaneous mirror symmetry breaking. *Physical Chemistry Chemical Physics* 25(46): 31583–31595. DOI: 10.1039/d3cp03311b.

Hsu HW, Postberg F, Sekine Y, et al. (2015) Ongoing hydrothermal activities within Enceladus. *Nature* 519(7542): 207–210. DOI: 10.1038/nature14262.

Johnson BC, Bowling TJ, Trowbridge AJ, et al. (2016) Formation of the Sputnik Planum basin and the thickness of Pluto's subsurface ocean. *Geophysical Research Letters* 43(19): 10068. DOI: 10.1002/2016GL070694.

Kang W, Mittal T, Bire S, et al. (2022) How does salinity shape ocean circulation and ice geometry on Enceladus and other icy satellites? *Science Advances* 8(29): eabm4665. DOI: 10.1126/sciadv.abm4665.

Kempf S, Tucker S, Altobelli N, et al. (2025) SUDA: a SUrface Dust Analyser for compositional mapping of the Galilean moon Europa. *Space Science Reviews* 221(1): 10. DOI: 10.1007/s11214-025-01134-0.

Klussmann M, Iwamura H, Mathew SP, et al. (2006) Thermodynamic control of asymmetric amplification in amino acid catalysis. *Nature* 441(7093): 621–623. DOI: 10.1038/nature04780.

Kondepudi DK and Nelson GW (1985) Weak neutral currents and the origin of biomolecular chirality. *Nature* 314(6010): 438–441. DOI: 10.1038/314438a0.

Lainey V, Rambaux N, Tobie G, et al. (2024) A recently formed ocean inside Saturn's moon Mimas. *Nature* 626(7998): 280–282. DOI: 10.1038/s41586-023-06975-9.

Leman L, Orgel L and Ghadiri MR (2004) Carbonyl sulfide-mediated prebiotic formation of peptides. *Science* 306(5694): 283–286. DOI: 10.1126/science.1102722.

Lente G (2007) The effect of parity violation on kinetic models of enantioselective autocatalysis. *Physical Chemistry Chemical Physics* 9(46): 6134. DOI: 10.1039/b711546f.

Levin SM, Zhang Z, Bolton SJ, et al. (2025) Europa's ice thickness and subsurface structure characterized by the Juno microwave radiometer. *Nature Astronomy* 10(1): 84–91. DOI: 10.1038/s41550-025-02718-0.

Lorenz RD, MacKenzie SM, Neish CD, et al. (2021) Selection and characteristics of the Dragonfly landing site near Selk crater, Titan. *The Planetary Science Journal* 2(1): 24. DOI: 10.3847/psj/abd08f.

Lunine J, Beauchamp J, Smith M, et al. (1999) The abiotic generation of homochirality on Saturn's moon Titan. In: [EDITORS: CONFIRM] (eds) *Advances in BioChirality*. [CITY: CONFIRM]: Elsevier, pp.257–270. DOI: 10.1016/B978-008043404-9/50017-6.

MacDermott AJ, Fu T, Hyde GO, et al. (2009) Electroweak parity-violating energy shifts of amino acids: the "conformation problem". *Origins of Life and Evolution of Biospheres* 39(5): 407–437. DOI: 10.1007/s11084-009-9161-x.

MacKenzie SM, Neveu M, Davila AF, et al. (2021) The Enceladus Orbilander mission concept: balancing return and resources in the search for life. *The Planetary Science Journal* 2(2): 77. DOI: 10.3847/psj/abe4da.

McCoy TJ, Russell SS, Zega TJ, et al. (2025) An evaporite sequence from ancient brine recorded in Bennu samples. *Nature* 637(8048): 1072–1077. DOI: 10.1038/s41586-024-08495-6.

McDougall TJ and Barker PM [INITIALS: CONFIRM] (2011) *Getting Started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox*. SCOR/IAPSO WG127.

Mora MF, Kok MG, Noell A, et al. (2022) Detection of biosignatures by capillary electrophoresis mass spectrometry in the presence of salts relevant to ocean worlds missions. *Astrobiology* 22(8): 914–925. DOI: 10.1089/ast.2021.0091.

Ozturk SF, Liu Z, Sutherland JD, et al. (2023a) Origin of biological homochirality by crystallization of an RNA precursor on a magnetic surface. *Science Advances* 9(23): eadg8274. DOI: 10.1126/sciadv.adg8274.

Ozturk SF, Bhowmick DK, Kapon Y, et al. (2023b) Chirality-induced avalanche magnetization of magnetite by an RNA precursor. *Nature Communications* 14: 6351. DOI: 10.1038/s41467-023-42130-8.

Parker ET, McLain HL, Glavin DP, et al. (2023) Extraterrestrial amino acids and amines identified in asteroid Ryugu samples returned by the Hayabusa2 mission. *Geochimica et Cosmochimica Acta* 347: 42–57. DOI: 10.1016/j.gca.2023.02.017.

Pavlov AA, McLain H, Glavin DP, et al. (2024) Radiolytic effects on biological and abiotic amino acids in shallow subsurface ices on Europa and Enceladus. *Astrobiology* 24(7): 698–709. DOI: 10.1089/ast.2023.0120.

Peter JS, Nordheim TA and Hand KP (2024) Detection of HCN and diverse redox chemistry in the plume of Enceladus. *Nature Astronomy* 8(2): 164–173. DOI: 10.1038/s41550-023-02160-0.

Petricca F, Vance SD, Parisi M, et al. (2025) Titan's strong tidal dissipation precludes a subsurface ocean. *Nature* 648(8094): 556–561. DOI: 10.1038/s41586-025-09818-x.

Pizzarello S and Cronin J (2000) Non-racemic amino acids in the Murray and Murchison meteorites. *Geochimica et Cosmochimica Acta* 64(2): 329–338. DOI: 10.1016/S0016-7037(99)00280-X.

Postberg F, Kempf S, Schmidt J, et al. (2009) Sodium salts in E-ring ice grains from an ocean below the surface of Enceladus. *Nature* 459(7250): 1098–1101. DOI: 10.1038/nature08046.

Postberg F, Sekine Y, Klenner F, et al. (2023) Detection of phosphates originating from Enceladus's ocean. *Nature* 618(7965): 489–493. DOI: 10.1038/s41586-023-05987-9.

Quack M, Seyfang G and Wichmann G (2022) Perspectives on parity violation in chiral molecules: theory, spectroscopic experiment and biomolecular homochirality. *Chemical Science* 13(36): 10598–10643. DOI: 10.1039/d2sc01323a.

Saur J, Duling S, Roth L, et al. (2015) The search for a subsurface ocean in Ganymede with Hubble Space Telescope observations of its auroral ovals. *Journal of Geophysical Research: Space Physics* 120(3): 1715–1737. DOI: 10.1002/2014JA020778.

Sohl F, Solomonidou A, Wagner FW, et al. (2014) Structural and tidal models of Titan and inferences on cryovolcanism. *Journal of Geophysical Research: Planets* 119(5): 1013–1036. DOI: 10.1002/2013JE004512.

Steel EL, Davila A and McKay CP (2017) Abiotic and biotic formation of amino acids in the Enceladus ocean. *Astrobiology* 17(9): 862–875. DOI: 10.1089/ast.2017.1673.

Strom C, Nordheim TA, Patthoff DA, et al. (2024) Constraining ocean and ice shell thickness on Miranda from surface geological structures and stress modeling. *The Planetary Science Journal* 5(10): 226. DOI: 10.3847/PSJ/ad77d7.

Tyler RH (2025) Ocean worlds maintained by ocean tidal heat expected on several Uranian moons. *Geophysical Research Letters* 52(20): e2025GL114911. DOI: 10.1029/2025GL114911.

Vance S, Bouffard M, Choukroun M, et al. (2014) Ganymede's internal structure including thermodynamics of magnesium sulfate oceans in contact with ice. *Planetary and Space Science* 96: 62–70. DOI: 10.1016/j.pss.2014.03.011.

Vance SD, Panning MP, Stähler S, et al. (2018) Geophysical investigations of habitability in ice-covered ocean worlds. *Journal of Geophysical Research: Planets* 123(1): 180–205. DOI: 10.1002/2017JE005341.

Viedma C (2005) Chiral symmetry breaking during crystallization: complete chiral purity induced by nonlinear autocatalysis and recycling. *Physical Review Letters* 94(6): 065504. DOI: 10.1103/PhysRevLett.94.065504.

Waite JH, Glein CR, Perryman RS, et al. (2017) Cassini finds molecular hydrogen in the Enceladus plume: evidence for hydrothermal processes. *Science* 356(6334): 155–159. DOI: 10.1126/science.aai8703.

Yu J, Darú A, Deng M, et al. (2024) Prebiotic access to enantioenriched amino acids via peptide-mediated transamination reactions. *Proceedings of the National Academy of Sciences* 121(7): e2315447121. DOI: 10.1073/pnas.2315447121.

Zeng Y and Jansen MF (2021) Ocean circulation on Enceladus with a high- versus low-salinity ocean. *The Planetary Science Journal* 2(4): 151. DOI: 10.3847/PSJ/ac1114.

Zeng Y and Jansen MF (2024) The effect of salinity on ocean circulation and ice–ocean interaction on Enceladus. *The Planetary Science Journal* 5(1): 13. DOI: 10.3847/psj/ad0cba.

Zhang Y, Kang W and Marshall J (2024) Ocean weather systems on icy moons, with application to Enceladus. *Science Advances* 10(45): eadn6857. DOI: 10.1126/sciadv.adn6857.
