# TNF-α Track 3 — pH-selective binder

Target: bind human TNF-α at pH 7.4, no detectable binding at pH 6.0. 20 designs.
Deadline: Mon 19 Oct 22:59 AEDT (23:59 AoE Sun 18).
Compute: 1× A100-80GB, ~13 GPU-h budget.

## Strategy

Two phases. GPU makes binders; CPU makes them pH-sensitive.

**Phase A (GPU).** Standard de novo campaign against the receptor-binding groove —
RFdiffusion → SolubleMPNN → AF2 initial-guess. Goal is *volume*: 300–500 designs
passing pAE_interaction < 10, not 20 perfect ones. The His scan downstream kills
most of what goes into it, so the funnel has to be wide at the top.

**Phase B (CPU, free).** Install histidine switches into the surviving interfaces and
score them by thermodynamic linkage. This is the part that actually wins Track 3 and
it costs no GPU time, so it can run on the instance's cores in parallel with Phase A's
later stages.

## Why these methods

- **RFdiffusion over Chai-2/BoltzGen.** Chai-2 beat TNFα in the BoltzGen benchmark and
  BoltzGen missed it, but plain binding is not the bottleneck here — surviving His
  substitution is. RFdiffusion yields the most backbones per GPU-hour and runs locally.
- **No BindCraft.** ~34 GPU-h per accepted design. Does not fit.
- **SolubleMPNN over vanilla ProteinMPNN.** Better expression in Adaptyv's cell-free
  pipeline, which is where these get tested.

## Epitope

Receptor-binding groove at the protomer–protomer interface. Unlike Track 2 there is no
mouse-conservation constraint, so the apex loop is back in play and that is where the
cationic residues are.

Hotspots (mature numbering, two-protomer crop A+B):
- Chain A bottom/E-F: D143, F144, A145, E146, S147 + S86, Y87
- Chain B apex loop: **R31, R32**, A33, N34, A35

**Switch-placement rule.** Designed His must sit near target **Arg/Lys** (R31, R32, K65,
K90, K112) — protonation at pH 6.0 then creates charge–charge repulsion. Keep His away
from target **Asp/Glu** (D143, E146, E23, E127): a protonated His beside a carboxylate
*strengthens* binding at low pH, which is the opposite of the objective.

## pH scoring

For an ionizable His with pKa_free (unbound binder) and pKa_bound (complex):

    ΔΔG = ΔG_bind(6.0) − ΔG_bind(7.4)
        = RT·ln[(1+10^(pKa_b−7.4))/(1+10^(pKa_f−7.4))]
        − RT·ln[(1+10^(pKa_b−6.0))/(1+10^(pKa_f−6.0))]

Want ΔΔG > 0 — weaker binding at pH 6.0. This happens when pKa_bound < pKa_free, i.e.
the complex suppresses protonation, so dropping the pH costs binding energy.

Targets: ~1.4 kcal/mol per 10× Kd shift. "No detectable binding" realistically needs
100–1000×, so **2.8–4.2 kcal/mol**, which means **2–3 cooperating His**, not one.
Single-His designs will not clear the bar.

Constraint pulling the other way: each His costs affinity at 7.4. Cap at 4 substitutions
and re-fold every variant to confirm the 7.4 interface survives.

## Schedule (hours from 12:30 Sun)

| h | Phase | Step |
|---|---|---|
| 0–1.5 | setup | Target prep, numbering map, hotspot verification |
| 1.5–6 | GPU | RFdiffusion, ~1000 backbones |
| 6–7 | GPU | SolubleMPNN, 8 seq/backbone (+ His-biased second pass) |
| 7–10 | GPU | AF2 initial-guess filter → expect 300–500 survivors |
| 10–14 | CPU | His scan + PROPKA3 linkage scoring (parallel with above tail) |
| 14–18 | GPU | Re-fold top His variants, confirm pH 7.4 interface intact |
| 18–22 | CPU | Developability filter, diversity clustering, pick 20 |
| 22–26 | — | Writeup, submission |
| 26–34 | — | Buffer |

## Ranking for final 20

1. Predicted ΔΔG(6.0 vs 7.4) ≥ 2.8 kcal/mol — the objective
2. AF2 pAE_interaction < 10 **after** His installation — binding at 7.4 must survive
3. His environment sanity — each switch His within 6 Å of a target Arg/Lys, none within
   6 Å of a target Asp/Glu
4. Composition sanity — reject poly-Ala/poly-Gly, free Cys, N-X-S/T sequons
5. Structural diversity — cluster at TM-score 0.6, cap 3 per cluster

Submit a diverse slate. Prior Adaptyv rounds show in-silico scores predict binding at
roughly chance (ROC-AUC ≈ 0.5–0.64), so concentrating all 20 on one scaffold family is
the main avoidable failure mode.

## Known traps

- **Numbering.** Three conventions coexist: mature (1 = Val), 1TNF crystal (6–157 modeled,
  1–5 disordered), and whatever the generator renumbers to. Build the map in step 01 and
  never hand-write a hotspot index.
- **ipTM on a homotrimer is meaningless.** Three of six chain pairs are target-internal
  and score 0.85+ regardless of binder. Compute per-pair binder↔target only.
- **1TNF's 29–36 loop is poorly ordered** and that is exactly the apex loop this track
  needs. Use 9OJO (1.36 Å) as the primary target, or graft the ordered loop from 4TSV.
