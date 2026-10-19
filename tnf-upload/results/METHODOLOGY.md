# Track 3 — pH-selective TNF-α binders: methodology

## Hypothesis

pH selectivity cannot be filtered with a structure predictor. AF2, Boltz and
Chai all ignore protonation state, so a folding model scores a design
identically at pH 7.4 and pH 6.0. We therefore treated the two requirements as
separate problems solved by separate tools: **binding** was designed and
filtered on GPU, **pH selectivity** was designed rationally and scored by
thermodynamic linkage on CPU.

The switch mechanism is histidine. With pKa ~6.0–6.5 it is neutral at pH 7.4
and cationic at pH 6.0. Placed at the interface beside a **target arginine or
lysine**, protonation creates charge–charge repulsion and the complex weakens.
Placed beside a target **aspartate or glutamate** it does the opposite —
protonation forms a favourable interaction and *strengthens* binding at low pH
— so carboxylate proximity was treated as disqualifying, not merely neutral.

For an ionisable group with pKa_free in the unbound binder and pKa_bound in the
complex:

    ΔΔG = ΔG_bind(6.0) − ΔG_bind(7.4)
        = RT·ln[(1+10^(pKa_b−7.4))/(1+10^(pKa_f−7.4))]
          − RT·ln[(1+10^(pKa_b−6.0))/(1+10^(pKa_f−6.0))]

ΔΔG > 0 means weaker binding at pH 6.0, which is the objective. It occurs when
pKa_bound < pKa_free: the complex suppresses protonation, so lowering the pH
costs binding energy. At ~1.4 kcal/mol per 10× Kd shift, "no detectable
binding" needs roughly 2.8–4.2 kcal/mol, i.e. **two to three cooperating
histidines**, not one.

## Epitope

The receptor-binding groove at the protomer–protomer interface of the trimer.
Track 3 imposes no mouse cross-reactivity constraint, which puts the apex loop
back in play — and that is where the cationic residues are. Hotspots:
**B31/B32 (an exposed arginine pair, SASA 200.7 and 93.0)**, B33, A87, A145,
A146. That arginine pair is the entire reason this epitope was chosen: it is
the counter-charge the switch histidines push against.

Target: **1TNF chains A+B**, two-protomer crop. Position 143 was patched
Leu→Asp to match UniProt P01375 and the assay construct TNA-H4211; the crop was
verified character-identical to the 157-residue construct over residues 6–157.
Ligand-bound TNF entries (9OJO, 9OJS, 6X81, 2AZ5) were deliberately avoided —
they carry smTNF modulators that distort the trimer at precisely the receptor
site being designed against.

## Pipeline

1. **Backbone generation — RFdiffusion.** 383 backbones against the target
   cropped to a 15 Å shell around the hotspot patch (173 of 304 target
   residues). The crop is a throughput measure: the full two-protomer target
   ran at 1.84 min/backbone versus 1.09 min/backbone cropped. Binder length
   sampled 70–110 aa, noise scales zeroed.

2. **Sequence design — ProteinMPNN (soluble weights).** Binder chain designed
   with the target as fixed context, 4 sequences per backbone, sampling
   temperature 0.1, cysteine forbidden at design time. Two pools were run in
   parallel: an unbiased pool and one with a mild histidine bias
   (`--bias_AA_jsonl {"H": 1.5}`) to seed switch positions during design rather
   than grafting them only post hoc. 12 of the final 20 come from the biased
   pool, 8 from the unbiased one.

3. **Binding filter — Boltz-2.** 1906 designs folded as binder + both target
   protomers. Two points materially affected the scores:

   * **The target was pinned to the crystal geometry via template.** In
     single-sequence mode without it, Boltz folds each protomer well
     (chain pTM ~0.85) but does **not dock them** — the target-internal pair
     ipTM was 0.291. Since the receptor groove sits *at* that interface, the
     epitope was effectively absent from the prediction. With the template it
     rises to ~0.71 and binder scores roughly double. Predictions where the
     protomers still failed to assemble (~8%) were discarded.
   * **Only binder↔target chain pairs were scored**, taking the minimum of the
     two. Global ipTM on a homotrimer crop is meaningless because the
     target-internal pair contributes regardless of the binder.

   The cutoff (pair ipTM ≥ 0.285) was calibrated in-house against
   composition- and length-matched shuffled sequences, whose scores topped out
   at 0.235. No published threshold transfers — AF2's `pae_interaction < 10` is
   a different model, metric and scale.

4. **pH scoring — PROPKA3 thermodynamic linkage.** Every interface position
   near a target cation and away from target carboxylates was substituted to
   histidine, the side chain built as a real imidazole (grafted on the N/CA/CB
   frame, χ1 chosen by clash), and pKa computed in the complex and in the free
   binder. Combinations were built additively, with two constraints:

   * **Histidines in a combination must be >8 Å apart (Cβ–Cβ).** Closer than
     that they couple electrostatically and the additive sum overestimates the
     pair — the dangerous direction, since it promotes designs that will not
     actually switch. Coupled combinations were rejected rather than scored.
   * **Pre-existing histidines are included in the linkage sum**, both the
     target's own (TNF carries three per protomer; B15 sits 6.6 Å from the
     hotspot patch) and any the biased MPNN pool installed. Each is referenced
     against the correct unbound state — free binder for binder His, free
     target for target His. Across 177 combinations this contribution was
     positive in 162 and negative in 9 (median +0.77 kcal/mol); in the final
     slate it contributes a median of +0.98 kcal/mol.

5. **pH 7.4 survival — Boltz-2 re-fold.** Every histidine combination was
   rebuilt as a sequence and re-folded. A design that switches beautifully but
   no longer binds at pH 7.4 is worthless, so each variant had to clear the
   same absolute pair-ipTM cutoff. A variant comfortably strong in absolute
   terms is exempt from the "lost too much versus parent" test, because
   otherwise designs descending from unusually strong parents are penalised
   for where they started rather than where they ended. 66 of 170 variants
   survived.

6. **Developability and diversity.** Rejected: free cysteine, N-X-S/T sequons,
   homopolymer runs >4, and a **total histidine ceiling of 8 (or 9% of
   length)**. The ceiling counts pre-existing histidines, not just designed
   substitutions, and removed 20 otherwise-surviving variants carrying 10–13
   histidines — beyond that the binder is broadly pH-labile rather than
   switched at a designed site. Survivors were clustered by TM-score at 0.6,
   capped at 2 variants per parent binder (so slots go to distinct binders
   rather than near-identical variants of one) and 8 per structural cluster.

## The submitted slate

20 designs, 77–108 aa, from **15 distinct parent binders across 6 structural
clusters**.

| metric | range | median |
|---|---|---|
| ΔΔG(6.0 vs 7.4) | 2.84 – 4.05 kcal/mol | 3.26 |
| implied Kd shift | 123× – 950× | ~250× |
| binder↔target pair ipTM (after His installation) | 0.327 – 0.650 | 0.403 |
| designed His substitutions | 2 – 4 | 3 |

Ranking is by ΔΔG first, with post-substitution interface confidence as the
tiebreak.

## Honest limitations

* **No positive control in the ipTM calibration.** The cutoff separates designs
  from shuffled sequences, so it means "above noise" — not "comparable to a
  known binder". Absolute ipTM values here are modest.
* **Structural diversity is the weakest axis.** Six clusters across 20 designs.
  Prior Adaptyv rounds show in-silico scores predict real binding at roughly
  chance (ROC-AUC 0.5–0.64), which makes scaffold concentration the main
  avoidable failure mode, and we are more concentrated than we would like.
* **The additive linkage model ignores coupling** between sites closer than
  8 Å. We rejected those combinations rather than approximating them, but the
  model is still additive among the ones we kept.
* **PROPKA pKa values are empirical**, computed on predicted rather than
  experimental complex structures, and the histidine side chains are modelled
  rotamers rather than repacked.
* **Track 3 asks for *no detectable* binding at 6.0, not a large fold-change.**
  A design shifting 50 nM → 25 µM moves 500× and may still register on BLI. We
  favoured designs with both a strong 7.4 interface and a large shift, but the
  absolute Kd at pH 6.0 is not something these methods establish.

## Design methods used

| method | role | link |
|---|---|---|
| RFdiffusion | backbone generation | https://github.com/RosettaCommons/RFdiffusion |
| ProteinMPNN (soluble weights) | sequence design | https://github.com/dauparas/ProteinMPNN |
| Boltz-2 | complex prediction, binding filter, re-fold | https://github.com/jwohlwend/boltz |
| PROPKA3 | pKa prediction for thermodynamic linkage | https://github.com/jensengroup/propka |
| TM-align (tmtools) | structural diversity clustering | https://github.com/jvkersch/tmtools |
| Biopython | structure handling, SASA, mutagenesis | https://biopython.org |

No PyRosetta and no `dl_binder_design`; the AF2 initial-guess /
`pae_interaction` filter was replaced by the calibrated Boltz-2 pair metric
described above.
