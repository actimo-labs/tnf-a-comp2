# TNF-α Track 3 binder design

**Read `HANDOFF.md` first — full state, strategy and rationale. `PLAN.md` has the
schedule and ranking criteria.**

Deadline: **Mon 19 Oct 22:59 AEDT**. Hard. Competition submission.

## What this is

Anthropic × Adaptyv Challenge 02, Track 3. Design 20 de novo single-chain
binders (10–250 aa) that bind human TNF-α at pH 7.4 and show **no detectable
binding at pH 6.0**. Not Track 2 — mouse cross-reactivity is not required.

## Environment

- Work in `~/tnf-track3`. `source env.sh` before anything.
- GPU: one A100-80GB, shared across all stages. Budget-constrained — do not
  launch speculative GPU jobs.
- RFdiffusion runs via Docker wrapper `./rfdiffusion.sh`. Pass **Hydra
  overrides only** (the image has an ENTRYPOINT) and **container-absolute
  paths** under `/work`. Include `hydra.run.dir=/work/.hydra`.
- Python: `venv-phscan/bin/python` for anything using biopython/propka.
  conda envs `mpnn` and `boltz` for those stages. There is no system biopython.
- `scripts/run_pipeline.sh` is **stale** — it expects `dl_binder_design`, which
  is not installed. Do not run it. Rewrite stages 03/04 against ProteinMPNN and
  Boltz-2, or delete it.

## Non-negotiables

- **`target/tnf_AB.pdb` is the design target.** Do not regenerate it from a
  different PDB entry. Ligand-bound TNF structures (9OJO, 9OJS, 6X81, 2AZ5)
  carry smTNF modulators that distort the trimer at the receptor site. Position
  143 is already patched Leu→Asp to match the assay construct.
- **pH selectivity cannot be scored by a folding model.** AF2, Boltz and Chai
  ignore protonation. It is scored by thermodynamic linkage on CPU in
  `scripts/05_ph_scan.py`. Do not substitute a confidence score for it.
- **Switch histidines go near target Arg/Lys (B31, B32), never near Asp/Glu
  (A143, A146).** A protonated His beside a carboxylate strengthens binding at
  low pH — backwards.
- **Score binder↔target chain pairs only.** The target is a homotrimer crop;
  global ipTM includes the target-internal A–B pair at ~0.9 and is meaningless.
- **Every His variant must be re-folded to confirm the pH 7.4 interface
  survives.** A design that switches but no longer binds at 7.4 is worthless.

## Working style

- Long GPU jobs run in tmux. Never start one in a bare foreground shell.
- Check in before anything that costs more than ~30 GPU-minutes.
- Say when a number comes from a run versus an estimate.
- Flag it when the plan in `PLAN.md` stops matching reality rather than
  quietly improvising around it.
