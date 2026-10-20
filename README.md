# tnf-a-comp2 — pH-switchable TNF-α binders

De novo single-chain binders designed to bind human TNF-α at **pH 7.4** and show
**no detectable binding at pH 6.0**.

Submitted to the Anthropic × Adaptyv Protein Design Challenge 02, Track 3
(October 2026). 20 designs, 77–108 aa, selected from 383 RFdiffusion backbones
and 2464 folded complexes. Everything here — scripts, target preparation,
intermediate scores, predicted structures and the submitted CSV — is the actual
run, not a cleaned-up reconstruction.

| | |
|---|---|
| ΔΔG(6.0 vs 7.4) | **2.84 – 4.56 kcal/mol** (median 3.27) → ~120–2200× Kd shift |
| binder↔target pair ipTM, after His installation | 0.294 – 0.650 (median 0.411) |
| diversity | 15 distinct parent binders, 6 structural clusters (TM ≥ 0.6) |
| hardware | **1 × NVIDIA A100-SXM4-80GB**, 12 vCPU, ~22 GPU-hours |

---

## The core idea

**pH selectivity cannot be filtered with a structure predictor.** AF2, Boltz and
Chai all ignore protonation state, so a folding model scores a design
identically at pH 7.4 and pH 6.0. Treating the two requirements as one problem
is the main trap in this track.

So the problem is split. **The GPU makes binders; the CPU makes them
pH-sensitive.**

The switch is histidine — pKa ~6.0–6.5, so neutral at 7.4 and cationic at 6.0.
Placed at the interface next to a **target arginine or lysine**, protonation
creates charge–charge repulsion and the complex weakens. Placed next to a target
**aspartate or glutamate** it does the opposite: protonation forms a favourable
interaction and *strengthens* binding at low pH, which is backwards. Carboxylate
proximity is therefore disqualifying, not merely neutral.

Selectivity is scored by thermodynamic linkage, on CPU, for free. For an
ionisable group with pKa_free in the unbound binder and pKa_bound in the complex:

```
ΔΔG = ΔG_bind(6.0) − ΔG_bind(7.4)
    = RT·ln[(1+10^(pKa_b−7.4))/(1+10^(pKa_f−7.4))]
      − RT·ln[(1+10^(pKa_b−6.0))/(1+10^(pKa_f−6.0))]
```

ΔΔG > 0 means weaker binding at pH 6.0 — the objective. It happens when
pKa_bound < pKa_free: the complex suppresses protonation, so dropping the pH
costs binding energy. At ~1.4 kcal/mol per 10× Kd shift, "no detectable binding"
needs roughly 2.8–4.2 kcal/mol, i.e. **two to three cooperating histidines**,
not one.

## Epitope

The receptor-binding groove at the protomer–protomer interface of the trimer.
Track 3 imposes no mouse cross-reactivity constraint, which puts the apex loop
back in play — and that is where the cationic residues are.

Hotspots: **B31/B32** (an exposed arginine pair, SASA 200.7 and 93.0), B33, A87,
A145, A146. That arginine pair is the counter-charge the switch histidines push
against, and the entire reason this epitope was chosen.

Target is `target/tnf_AB.pdb` — 1TNF chains A+B, two-protomer crop, position 143
patched Leu→Asp to match UniProt P01375 and Adaptyv's assay construct TNA-H4211,
verified character-identical to that construct over residues 6–157. Ligand-bound
TNF entries (9OJO, 9OJS, 6X81, 2AZ5) are deliberately avoided: they carry smTNF
modulators that distort the trimer at precisely the receptor site being designed
against.

## Pipeline

| stage | script | what it does |
|---|---|---|
| 01 | `scripts/01_prep_target.py` | fetch 1TNF, crop A+B, patch 143, build numbering map |
| 02 | `scripts/02_crop_target.py` | 15 Å shell around the hotspots for RFdiffusion |
| 02 | `rfdiffusion.sh` | backbone generation via the Docker image |
| 03 | `scripts/03_mpnn.py` | ProteinMPNN, soluble weights, two pools (plain + His-biased) |
| 04 | `scripts/04_boltz.py` | Boltz-2 complex prediction and binding filter |
| 04b/c | `04b_make_controls.py`, `04c_calibrate.py` | shuffled-sequence negatives, cutoff calibration |
| 05 | `scripts/05_ph_scan.py` | PROPKA3 thermodynamic linkage scan (CPU, parallel) |
| 06 | `scripts/06_his_refold.py` | re-fold every His variant, confirm pH 7.4 binding survives |
| 07 | `scripts/07_submit.py` | developability, TM-score clustering, ranked CSV |

Driven by `run_all.sh` per stage, or `round2.sh` / `drive_*.sh` for unattended
multi-stage runs.

```bash
source env.sh
./run_all.sh 02      # backbones (tmux, long)
./run_all.sh 03      # sequences
./run_all.sh 04      # Boltz screen
./run_all.sh 05      # pH linkage scan (CPU)
CUTOFF=0.285 ./run_all.sh 06
./run_all.sh 07      # submission CSV
```

`setup/install.sh` builds the environments. RFdiffusion runs from the
`rosettacommons/rfdiffusion` Docker image; ProteinMPNN and Boltz-2 run in conda
environments; the pH scan runs in a venv with biopython + propka.

## Measured throughput

All on the single A100-80GB.

| operation | rate |
|---|---|
| RFdiffusion, full two-protomer target (~382 res) | 1.84 min/backbone |
| RFdiffusion, 15 Å cropped target (~263 res) | **1.09 min/backbone** |
| Boltz-2, templated complex, batched | **~20 s/design** |
| PROPKA linkage scan, 10 processes on 12 cores | ~6 s/design |

The crop exists because 1000 backbones against the full target is 30.7 GPU-h
against a ~30 h deadline. Cropping to a 15 Å shell keeps 173 of 304 target
residues and cuts that by ~40%.

## Four things that are easy to get wrong

Each of these was a real bug in this run, and each produced plausible-looking
output while being wrong.

**1. Boltz does not assemble the target on its own.** In single-sequence mode
without a template, it folds each TNF protomer well (chain pTM ~0.85) but does
**not dock them to each other** — the target-internal pair ipTM was 0.291. Since
the receptor groove *is* the protomer–protomer interface, the epitope was
effectively absent from every prediction. Pinning both target chains to the
crystal via template (`scripts/make_template_cif.py` → `target/tnf_AB.cif`)
raises it to ~0.71 and roughly doubles binder scores.

**2. A histidine without its imidazole is invisible to PROPKA.** An early version
of the scan truncated substituted residues to CB. PROPKA identifies titratable
histidines by their ring atoms, so it did not register them as His at all and
silently omitted them from its output — the scan returned **zero candidates for
every design**, which looked like the designs being bad. The scan now grafts a
real imidazole on the N/CA/CB frame with χ1 chosen by clash.

**3. Pre-existing histidines have to be counted.** Both the target's own (TNF has
three per protomer; B15 sits 6.6 Å from the hotspot patch) and any the
His-biased MPNN pool installed. Each needs the correct unbound reference — free
binder for binder His, free target for target His. Across 177 combinations this
contribution was positive in 162 and negative in 9 (median +0.77 kcal/mol), so
omitting it is not a small correction.

**4. Score binder↔target chain pairs only.** Global ipTM on a two-protomer crop
is dominated by the target-internal pair and means nothing. We take the minimum
of the two binder–target pairs: a binder that grips one protomer and ignores the
other is not binding the groove.

Also enforced: histidines within a combination must be **>8 Å apart (Cβ–Cβ)**.
Closer than that they couple electrostatically and the additive ΔΔG
*overestimates* the pair — the dangerous direction, since it promotes designs
that will not actually switch. Coupled combinations are rejected rather than
approximated.

## Calibration

No published threshold transfers here — AF2's `pae_interaction < 10` is a
different model, metric and scale. The cutoff (pair ipTM ≥ 0.285) was calibrated
against composition- and length-matched **shuffled** sequences, which topped out
at 0.235.

This calibrates the null only. The cutoff means **"above noise"**, not
"comparable to a known binder". `scripts/04b_make_controls.py --positives`
accepts a FASTA of known binders to anchor the top end; we did not run it.

## Repo layout

```
scripts/            pipeline stages 01–07 plus calibration helpers
run_all.sh          per-stage driver
round2.sh, drive_*.sh   unattended multi-stage drivers
setup/install.sh    environment setup
target/             design target, crop, template CIF, numbering, hotspots
results/
  submission.csv            the 20 submitted designs
  submission_report.tsv     per-design ΔΔG, ipTM, His count, mutations, cluster
  METHODOLOGY.md            full write-up
  ph_scan.tsv               every candidate His site scored individually
  ph_combos.tsv             combinations and their linkage ΔΔG
  06_refold.tsv             every variant re-folded, with pH 7.4 verdict
  structures/
    complex/      20 predicted complexes (chain X binder, A+B target)
    binder_only/  chain X extracted
```

Earlier slates are kept (`submission_20_uncorrected.csv`,
`submission_pre_round2.csv`, …) so the effect of each correction is traceable
rather than asserted.

**Note on the structures:** the target chains were template-pinned, so their
coordinates are close to crystal geometry *by construction* and should not be
read as a prediction. Only the binder pose is predicted. See
`results/structures/README.txt`.

## Limitations

- **No positive control in the calibration.** Absolute ipTM values are modest
  and the cutoff only separates designs from shuffled sequences.
- **Structural diversity is the weakest axis** — 6 clusters across 20 designs.
  Prior Adaptyv rounds show in-silico scores predict real binding at roughly
  chance (ROC-AUC 0.5–0.64), which makes scaffold concentration the main
  avoidable failure mode, and this slate is more concentrated than intended.
- **The linkage model is additive** among the sites kept, and PROPKA pKa values
  are empirical, computed on predicted structures with modelled rotamers rather
  than repacked side chains.
- **Track 3 asks for *no detectable* binding at pH 6.0, not a large fold-change.**
  A design shifting 50 nM → 25 µM moves 500× and may still register on BLI.
  These methods do not establish the absolute Kd at pH 6.0.
- **Two designs carry unusual composition** (37% Ala, 41% Glu) — a known
  ProteinMPNN bias on idealised de novo helices, with some molten-globule risk.
  They were kept because filtering them cost a whole structural cluster.

## Methods used

[RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) ·
[ProteinMPNN](https://github.com/dauparas/ProteinMPNN) (soluble weights) ·
[Boltz-2](https://github.com/jwohlwend/boltz) ·
[PROPKA3](https://github.com/jensengroup/propka) ·
[TM-align / tmtools](https://github.com/jvkersch/tmtools) ·
[Biopython](https://biopython.org)

No PyRosetta and no `dl_binder_design`; the AF2 initial-guess /
`pae_interaction` filter is replaced by the calibrated Boltz-2 pair metric.

## License

MIT — see [LICENSE](LICENSE). Model weights and third-party tools carry their
own licences. `target/raw/1tnf.pdb` is from the RCSB PDB.
