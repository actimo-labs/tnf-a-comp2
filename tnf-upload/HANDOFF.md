# TNF-α Track 3 — handoff

Read this plus `PLAN.md` before touching anything. Written 2026-10-18 ~16:30 AEDT.

## Deadline

**Mon 19 Oct 22:59 AEDT** (23:59 AoE Sun 18). ~30 hours from writing.

## The task

Anthropic × Adaptyv Challenge 02, Track 3: a de novo single-chain binder,
10–250 aa, that binds human TNF-α at pH 7.4 and shows **no detectable binding
at pH 6.0**. Submit 20. De novo only. Adaptyv validates by cell-free
expression + BLI/SPR.

Not Track 2 — mouse cross-reactivity is NOT required. (Earlier work in the
chat assumed Track 2; that constraint is dropped and the epitope choice
inverts because of it, see below.)

## Environment

GCP A100-80GB, `~/tnf-track3`, user HarjasG.

| Thing | Where | State |
|---|---|---|
| pH scan venv (CPU) | `venv-phscan/` | works — biopython + propka |
| RFdiffusion | Docker `rosettacommons/rfdiffusion`, wrapper `rfdiffusion.sh` | image pulled, weights in `RFdiffusion/models/` |
| ProteinMPNN | `ProteinMPNN/` | cloned, torch 2.5.1+cu121 cuda True, soluble weights present |
| Boltz-2 | conda env `boltz` | installed, **not yet wired into any script** |
| dl_binder_design | not installed | needs PyRosetta licence — deliberately skipped |

`source env.sh` first. It sources conda by absolute path and makes the output dirs.

## What is done

- `target/tnf_AB.pdb` — 1TNF chains A+B, two-protomer crop. **This is the target.**
- `target/numbering.json` — numbering map, epitope sets, charge annotations.
- `target/hotspots.txt` — `B31,B32,B33,A87,A145,A146`
- 1TNF position 143 was Leu in the crystal but Asp in P01375 and in Adaptyv's
  assay construct TNA-H4211. **Already patched to ASP in both chains.** Verify
  with `grep -c "ASP A 143\|ASP B 143" target/tnf_AB.pdb` → expect 2+.
  Numbering verified clean 136–152 otherwise, no insertion codes.

## The strategy, and why

**pH switching is the whole game.** Structure predictors ignore protonation —
AF2, Boltz, Chai all do. So pH selectivity cannot be filtered with a folding
model. It is designed rationally and scored by thermodynamic linkage on CPU.
GPU makes binders; CPU makes them pH-sensitive.

**Mechanism.** Histidine, pKa ~6.0–6.5: neutral at 7.4, cationic at 6.0. Place
His at the interface next to a **target Arg/Lys** so protonation creates
charge–charge repulsion and the complex falls apart. B31/B32 are an exposed
arginine pair (SASA 200.7 and 93.0) sitting in the receptor groove — that pair
is the entire reason the chosen epitope works. Keep His **away** from target
Asp/Glu (A146, A143, E23): protonated His beside a carboxylate *strengthens*
binding at low pH, which is backwards.

**Scoring.** For His with pKa_free (unbound binder) and pKa_bound (complex):

    ddG = dG_bind(6.0) − dG_bind(7.4)

Positive ddG = weaker at pH 6.0 = the objective. Happens when
pKa_bound < pKa_free. ~1.4 kcal/mol per 10× Kd shift; "no detectable binding"
needs 2.8–4.2 kcal/mol, so **2–3 cooperating His**, not one. Implemented in
`scripts/05_ph_scan.py` via PROPKA3.

**Tension to manage:** every His costs affinity at 7.4. Cap at 4. Re-fold every
variant to confirm the 7.4 interface survives. A design that switches
beautifully but no longer binds at 7.4 is worthless.

**Why RFdiffusion and not BindCraft/Chai-2/BoltzGen.** His installation
destroys most designs, so the funnel must be wide — hundreds of backbones, not
fifteen excellent ones. RFdiffusion gives the most backbones per GPU-hour and
runs locally. BindCraft is ~34 GPU-h per accepted design; doesn't fit.

## Current step — RFdiffusion

Run from `~/tnf-track3` inside tmux. `$PWD` mounts at `/work`, so use
container-absolute paths:

```bash
./rfdiffusion.sh \
  inference.output_prefix=/work/rfdiff/tnf \
  inference.input_pdb=/work/target/tnf_AB.pdb \
  'contigmap.contigs=[A6-157/0 B6-157/0 70-110]' \
  "ppi.hotspot_res=[$(cat target/hotspots.txt)]" \
  inference.num_designs=1000 \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 \
  hydra.run.dir=/work/.hydra \
  2>&1 | tee logs/02_rfdiffusion.log
```

`hydra.run.dir` is required — the container runs as your UID and cannot write
Hydra's default output dir inside root-owned `/app`.

~4 GPU-h for 1000 backbones. First run spends several minutes on
`Calculating IGSO3` before output; not a hang.

## What still needs building

`scripts/run_pipeline.sh` is **stale** — written before the installer, expects
`dl_binder_design` for MPNN and AF2. Do not run it. Stages needed:

**03 — ProteinMPNN.** `ProteinMPNN/protein_mpnn_run.py`, design the binder
chain only with the target as fixed context: `--pdb_path_chains` for the binder
chain, `--use_soluble_model`, `--num_seq_per_target 8`, `--sampling_temp 0.1`.
Also worth a second pass with a mild His bias via `--bias_AA_jsonl` ({"H": 1.5})
to seed switch positions during design rather than relying only on the
post-hoc scan.

**04 — Boltz-2 validation.** BUILT: `scripts/04_boltz.py`. Score binder↔target
chain pairs ONLY — but see the correction below, the reason is not the one
originally given here. Cutoff calibrated by `04c_calibrate.py`: **0.285**.

**05 — pH scan.** FIXED. The >8 Å rule was documented in the docstring but
never implemented; it is now enforced. It also could never have produced
output at all — see the correction below.

**06 — BUILT:** `scripts/06_his_refold.py`. Gates on both the absolute cutoff
and fractional drop vs parent. **07 — BUILT:** `scripts/07_submit.py`,
developability + TM-score clustering + the ranked submission CSV.


## Corrections from the 18 Oct session — read these before trusting the above

Three claims in this document are wrong, and two of them were hiding silent
failures. Measured, not estimated.

**1. The target-internal A–B pair does not score ~0.9.** It scores **0.291**.
Boltz folds each protomer well (ptm 0.85) but does not dock them to each
other. Since the receptor-binding groove sits *at* the protomer–protomer
interface, the epitope was absent from every prediction — the binder was being
scored against a surface that did not exist. Fix: pin both target chains to the
crystal with a template (`scripts/make_template_cif.py` →
`target/tnf_AB.cif`). A–B then rises to **0.708** and binder pair scores roughly
double. Scoring per-pair is still right; it was just not sufficient.
A `target_pair_iptm` sanity gate now rejects predictions that fail to assemble
(~8% of them).

**2. Stage 05 returned zero candidates for every design, always.**
`mutate_to_his` truncated the side chain to CB. PROPKA identifies titratable
histidines by their ring atoms, so with no imidazole it did not register the
residue as His at all and omitted it from SUMMARY; the caller read back `None`
and dropped the site. Silent — it looked like the designs were bad. Fixed by
grafting a real imidazole (superimposed on the N/CA/CB frame, chi1 chosen by
clash). A test design now yields sites at 1.14 and 0.51 kcal/mol.

**3. The Asp143 side chain is truncated to CB in `target/tnf_AB.pdb`**, but
this does NOT matter for scoring, contrary to what it looks like. Stages 05
and 06 run on Boltz output, where the target is rebuilt from the full 157-aa
sequence and comes out complete. It affects only the template input. Do not
"fix" the protected target file for this.

### Target verified against the real assay construct

Pulled from the Proteinbase challenge page. TNA-H4211 is **157 aa and
genuinely tag-free** (ends `...DFAESGQVYFGIIAL`, no polyhistidine — a His tag
on the target would have titrated at pH 6.0 and confounded the switch).
`target/tnf_AB.pdb` residues 6–157 match it with **zero mismatches**, and
position 143 is D, so the Leu→Asp patch was correct. Mouse TNF is not needed;
Track 3 drops cross-reactivity.

Submission format: ranked CSV, columns `name,sequence,molecule_class`
(="protein"). Deadline confirmed: Sun 18 Oct 23:59 AoE = **Mon 19 Oct 22:59
AEDT**.

### Native target histidines matter

TNF has three His per protomer (15, 73, 78) and **B15 sits 6.6 Å from the
hotspot patch**. It titrates between 7.4 and 6.0 regardless of what the binder
does. Stage 05 now includes native target His within 10 Å in the linkage sum
(`native_his_ddg` column) — the contribution can be negative, eating into the
designed switch. A binder carboxylate placed near B15 would reverse the switch.

### Environment fixes

- RFdiffusion Docker died on `PermissionError` creating
  `/app/RFdiffusion/schedules` — it runs as the host UID and `/app` is
  root-owned. `rfdiffusion.sh` now bind-mounts a writable `schedules/`.
  This was never an entrypoint problem.
- The `boltz` conda env had `torch 2.14.1+cu130` against driver 550
  (CUDA 12.4): GPU completely unusable. Reinstalled `torch==2.14.1+cu126`
  (same torch version, no dependency churn). Also needed
  `cuequivariance-torch` + `cuequivariance-ops-torch-cu12` for the triangle
  kernels, else `--no_kernels` is required.
- `tmtools` installed in `venv-phscan` for TM-score clustering; there is no
  TMalign/USalign binary on this box.

### Measured throughput

| stage | rate |
|---|---|
| RFdiffusion, full 2-protomer target | 1.84 min/backbone |
| RFdiffusion, 15 Å cropped target | **1.09 min/backbone** |
| Boltz-2, batched, templated | **20.7 s/design** |

The crop (`scripts/02_crop_target.py`, 173/304 target residues) exists because
1000 backbones against the full target is 30.7 GPU-h against a 29.6 h
deadline. `PLAN.md`'s "~4 GPU-h for 1000 backbones" is wrong by ~7×.

### Still true, still worth heeding

The funnel guidance holds: under ~150 survivors into the His scan is too thin.
At the calibrated cutoff of 0.285, 46% of designs pass, so screen ~500
sequences to clear it.

The calibration has **no positive control** — the cutoff is "clearly above
noise", not "as good as a known binder". Adding known TNF binders via
`04b_make_controls.py --positives` would anchor the top end.

### Known-stale in PLAN.md

- "Use 9OJO as the primary target, or graft the ordered loop from 4TSV"
  contradicts the CLAUDE.md non-negotiable. Ignored; `target/tnf_AB.pdb`
  is authoritative. `target_4tsv/raw/` is empty.
- AF2 `pae_interaction < 10` appears in Phase A and ranking criterion 2.
  `dl_binder_design` is not installed. Replaced by the Boltz pair cutoff.
- The schedule counts hours from 12:30 Sun and is ~5 h behind reality.

## Ranking the final 20

1. ddG(6.0 vs 7.4) ≥ 2.8 kcal/mol
2. Interface confidence **after** His installation — 7.4 binding must survive
3. Each switch His within 6 Å of target Arg/Lys, none within 6 Å of Asp/Glu
4. Composition sanity — no poly-Ala/poly-Gly, no free Cys, no N-X-S/T sequons
5. Diversity — cluster at TM-score 0.6, cap 3 per cluster

Prior Adaptyv rounds: in-silico scores predict real binding at roughly chance
(ROC-AUC ≈ 0.5–0.64). **Submit a diverse slate.** Concentrating all 20 on one
scaffold family is the main avoidable failure mode.

Track 3 asks for *no detectable* binding at 6.0, not a large fold-change. A
design going 50 nM → 25 µM shifts 500× but may still register on BLI. Favour
strong 7.4 affinity **and** a large shift, not shift alone.

## Traps already hit — do not rediscover

- Biopython `PDBList` uses a retired RCSB endpoint, 403s on everything. Fetch
  from `https://files.rcsb.org/download/<ID>.pdb` directly. Already fixed.
- A real PDB file has hundreds of header lines before the first ATOM record —
  don't validate a download by scanning only its first few KB.
- **Do not use ligand-bound TNF entries** (9OJO, 9OJS, 6X81, 2AZ5) as the design
  target. They carry smTNF modulators that distort the trimer at the receptor
  site — the exact surface being designed against. Adaptyv assays apo trimer.
  9OJO is 1.36 Å and tempting; it is still wrong.
- The RFdiffusion image has an ENTRYPOINT. Pass Hydra overrides only — no
  `python`, no script path, no `-w` override.
- If N surviving designs after stage 04 is **under ~150**, the funnel is too
  thin: the His scan keeps only 10–20%. Raise `inference.num_designs`, do not
  loosen the confidence cutoff — a design failing confidence won't survive His
  substitution either.

## Open decision

PyRosetta licence. `dl_binder_design` needs one (free academic, paid
commercial; Actimo has paying clients). Taking it restores the validated AF2
`pae_interaction < 10` filter and removes the Boltz calibration step. Current
path assumes no licence.
