# lit_r6: the eight literature R6 letters under the seed rule

Written 2026-10-21. Applied the same day with the author's OK: the rule R letters (ija_draft change log E114) and
side findings 1-4 (E115).

## What ran

| file | role |
|---|---|
| `lit_r6_predictions.txt` | rules R and L, predicted letters, evidence predictions F1-F4, A1-A2, K1-K2; written and hashed before any evidence was gathered (sha256 `f0345e5059ec0c1de6124d19f5f6c54435d6ace047840c47d0c908f8901238ed`, asserted in the script) |
| `lit_r6.py` | A1-A2 (aped.py definitions, unedited) and K1-K2 (scorecard.py imported, unedited); 1 s |
| `lit_r6.log` | output; `lit_r6.err` is empty |
| papers | F1-F4 read by three agents, open text only; no paywall, Cloudflare challenge or TLS check was bypassed. Downloads are in the session scratchpad, not the repo |

Rule R is the approved seed rule (arm_e_r6.py:262), applied with three kinds of evidence: (a) a model at measured rates;
(b) a deciding count, where a hand decided among N_d units cannot follow a seed below ~2/sqrt(N_d), so following at e_need
needs N_d >= 1.2e23 (pond), 3.8e27 (lake), 5.1e35 (ocean); (c) a seed that decays, is only passed through, or does not set
the hand. Rule L is a lenient alternative: P if a seed of <= 1% grows and sets the hand in the lab or in a model of the
mechanism. Arms E and G are N under both rules.

## Predictions, scored as written: 7 of 8 held

| | result | evidence |
|---|---|---|
| F1 Saghatelian 2001 | **not held** (2 of 3 clauses) | Held: no run gives a product ee above the seed ee; the seeds are enantiopure T_LL with single-hand fragments, and there is no enantiomeric cross-inhibition (p. 799). Held with caveat: homochiral products from racemic fragments are plotted together, so "equal amounts" is implied, not measured (the metric is diastereomeric excess). Failed: the paper states no growth order, so "sub-exponential" has no support in it |
| F2 Ozturk 2023 | held | The hand is set by the surface magnetization direction. The start is 65 mM racemic; no seed of <= 1% ee is used. The Nat. Commun. avalanche paper starts enantiopure |
| F3 Noorduin 2008 | held | The smallest initial ee in JACS 130:1158 is 2-3%, and the same hand is reached (thesis Fig. 2.1, SI p. 56). Crystals decide; no crystal count is given, only mass (4 g). The Angew 47:6445 seed ee is not stated in the thesis chapter |
| F4 Klussmann 2006 | held | Proline in DMSO: solution ee 47-54% across 20-80% overall ee (Supp. Table 1). The solid + solution mass balance reproduces each employed ee, so overall ee is only partitioned, not grown. The sign follows the overall excess, but only L-excess and 0% rows exist (0% gives -0.9, 0.2, -0.3) |
| A1 | held | sigma_max < 0 at ctot 0.1, 1e-2, 1e-4, 1e-6 M (lit_r6.log:4-8). The grid is 8 points per decade, so a very narrow window could be missed |
| A2 | held | -9.2864e-07 at 0.5 M reproduces aped.log:91 (lit_r6.log:3) |
| K1 | held | Ghadiri and Joyce tie at 5.5 under R; Ghadiri stays first only by dict order (lit_r6.log:33) |
| K2 | held | R8 = N block: R 5.0, 4.0, 3.5, 2.5; L 5.5, 4.5, 3.5, 3.0 (lit_r6.log:34, 47) |

Sources: Noorduin numbers come from the author's PhD thesis (reprints of both papers), not the journal PDFs. Klussmann
numbers come from the public Nature SI; the article body is paywalled. The Saghatelian text was read from a third-party
copy, so cite doi:10.1038/35057238 only.

## Proposed letters (not applied)

| row | now | R | L | deciding fact |
|---|---|---|---|---|
| Ghadiri peptide replicator (2001) | Y | N | N | (c) No run follows a seed; nothing amplifies ee (F1) |
| Joyce RNA cross-inhibition (1984) | N | N | N | (c) No self-copying, so nothing to seed |
| Soai reaction (1995; 2003) | Y | N | P | R: (b) N_d <~ 3.6e9 (soai_check.py S4). L: a 5e-7 seed sets the hand (S5) |
| Viedma ripening, NaClO3 (2005) | Y | N | P | R: (b) crystals decide. L: viedma_big.py follows ~0.4% (Delta = 3 at 256,000 crystals) |
| Attrition deracemization, aa deriv (2008) | Y | N | P | R: (b) crystals decide (F3). L: the grinding model, plus 0.3-0.35% ee in a later imine from the same group (thesis ch. 5) |
| Magnetite crystallization (2023) | Y | N | N | (c) The field-set magnetization sets the hand (F2) |
| Blackmond eutectic (2006) | P | N | P* | R: (b) crystals carry the excess. L: see below |
| APED peptides via NCAs (2004-2010) | P | N | P* | R: (c) racemic is stable at every setting's ctot, so a seed decays (A1). L: see below |

\* Two L letters are judgment calls:
- **Eutectic.** The smallest overall ee tested is 5%, not <= 1%. The solution ee rises only because the excess
  partitions between phases; the overall ee never grows (F4). L gives P only if partitioning counts as growth.
  Otherwise it gives N.
- **APED.** aped.py grows a 1e-6 seed only at ctot 2 M (aped.log:47-49). At every setting's ctot the seed decays
  (A1). That decay is the clause that keeps arm E at N under L, so reading APED by the arms' standard gives N.

Recommendation: rule R. It is the approved rule, it treats arms and literature alike, and all eight rows land on
evidence of type (b) or (c) with no judgment call.

## Score and rank effects (lit_r6.log)

- **Under R, every row loses its R6 credit:** Ghadiri 6.5 -> 5.5, Soai 4.5 -> 3.5, Attrition 6.0 -> 5.0,
  Viedma 5.0 -> 4.0, Magnetite 4.5 -> 3.5, APED 5.0 -> 4.5, Eutectic 3.0 -> 2.5. Joyce stays at 5.5.
- **Order is unchanged except at the top.** Ghadiri now ties Joyce at 5.5. The scorecard.py:87 assert
  (`ranked[0] == "Ghadiri ..."`) still passes, but only because Ghadiri comes first in the dict. The same tie holds
  under L.
- **Fix if R is applied:** either assert that Ghadiri is among the top-scoring routes, or state the tie in the text.
  README:464 (the table) would also need the new letters. ija_draft does not name the top route.
- **The R8 = N block keeps its order under both rules,** so the README:535 claim stands: crystal and surface routes
  fall below every solution route.

## Side findings (all four fixed 2026-10-21, E115)

1. README:526 and scorecard.py:60 cite "growth is sub-exponential (rho ~0.5-0.65)" to Saghatelian 2001. Neither
   appears in that paper's main text. The SI was not read, so the source is unknown: either cite it or remove the
   clause. The R9 = N inference itself still holds, since it rests only on own-hand catalysis without
   cross-inhibition, which the paper states. **Fixed:** clause removed; README and scorecard.py now say the main
   text gives no growth order.
2. README:798 and scorecard.py:63 give Ozturk 2023 as "racemic -> ~60% ee". That is the best single run. Repeat
   runs average about 35-38% (down 59, 36, 18; up -10, -22, -41, -54), so say "up to ~60%". **Fixed** in README
   (Part 6 and the reference), scorecard.py and ija_draft Section 1.1 ("about 35% on average", the authors' figure).
3. In Noorduin's JACS work, every racemic start went to the R solid. The authors suspect a trace chiral impurity,
   and 0.1 mol% of (S)- or (R)-phenylglycine sets the opposite hand. This is direct support for (c): the crystals
   follow whatever biases them, not the pool. **Added** to the README R6 paragraph and the scorecard.py attrition
   comment, with the JACS 130:1158 reference added to README.
4. scorecard.py:62 says the eutectic row used "CHCl3". The Klussmann SI shows DMSO, water and DMF. The article body
   was not read, so check before changing it. **Fixed:** the abstract (open) names no solvent and the SI lists only
   DMSO, water and DMF, so the comment now says "proline in DMSO or water (Klussmann 2006 SI)".
