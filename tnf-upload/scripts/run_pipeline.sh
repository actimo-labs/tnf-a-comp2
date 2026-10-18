#!/usr/bin/env bash
# TNF-alpha Track 3 pipeline. Single A100-80GB, ~13 GPU-h budget.
# Each stage is idempotent and checkpointed -- a crash resumes, it does not restart.
set -euo pipefail

RFD=${RFD:-$HOME/RFdiffusion}
DLBD=${DLBD:-$HOME/dl_binder_design}
NBACKBONES=${NBACKBONES:-1000}
SEQS_PER=${SEQS_PER:-8}
BINDER_LEN=${BINDER_LEN:-"70-110"}
PAE_CUT=${PAE_CUT:-10}

mkdir -p target rfdiff mpnn af2 results logs

stage_done() { [[ -f "logs/$1.done" ]]; }
mark_done()  { touch "logs/$1.done"; }

# ---------------------------------------------------------------- 01 target
if ! stage_done 01; then
  echo "=== 01 target prep ==="
  python scripts/01_prep_target.py --pdb 9OJO --chains A B --outdir target
  mark_done 01
fi
HOTSPOTS=$(cat target/hotspots.txt)
echo "hotspots: $HOTSPOTS"

# ------------------------------------------------------------ 02 rfdiffusion
# Two-protomer crop. noise_scale=0 tightens PPI interfaces.
# ~4 GPU-h for 1000 backbones at this length.
if ! stage_done 02; then
  echo "=== 02 RFdiffusion ($NBACKBONES backbones) ==="
  "$RFD/scripts/run_inference.py" \
    inference.output_prefix=rfdiff/tnf \
    inference.input_pdb=target/tnf_AB.pdb \
    "contigmap.contigs=[A6-157/0 B6-157/0 ${BINDER_LEN}]" \
    "ppi.hotspot_res=[${HOTSPOTS}]" \
    inference.num_designs="$NBACKBONES" \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0 \
    2>&1 | tee logs/02_rfdiffusion.log
  mark_done 02
fi

# -------------------------------------------------------------- 03 mpnn
# SolubleMPNN -- better cell-free expression, which is how Adaptyv tests these.
# Two passes: unbiased for affinity, His-biased to seed switch positions directly
# rather than relying only on the post-hoc scan in step 05.
if ! stage_done 03; then
  echo "=== 03 SolubleMPNN ==="
  python "$DLBD/mpnn_fr/dl_interface_design.py" \
    -pdbdir rfdiff/ \
    -outpdbdir mpnn/unbiased/ \
    -seqs_per_struct "$SEQS_PER" \
    -checkpoint_path "$DLBD/mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt" \
    2>&1 | tee logs/03a_mpnn.log

  # His bias: nudges histidine into interface positions during design.
  # Keep it mild -- too much and affinity at 7.4 collapses before step 05 sees it.
  python "$DLBD/mpnn_fr/dl_interface_design.py" \
    -pdbdir rfdiff/ \
    -outpdbdir mpnn/hisbias/ \
    -seqs_per_struct 4 \
    -checkpoint_path "$DLBD/mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt" \
    -bias_AA "H:1.5" \
    2>&1 | tee logs/03b_mpnn_his.log
  mark_done 03
fi

# ---------------------------------------------------------- 04 af2 filter
# AF2 initial-guess. pae_interaction < 10 is the headline filter.
# IMPORTANT: this target is a homotrimer crop. Global ipTM is meaningless here --
# the A-B pair is target-internal and scores ~0.9 regardless of the binder.
# Score binder-vs-target pairs only.
if ! stage_done 04; then
  echo "=== 04 AF2 initial-guess ==="
  for src in unbiased hisbias; do
    python "$DLBD/af2_initial_guess/predict.py" \
      -pdbdir "mpnn/$src/" \
      -outpdbdir "af2/$src/" \
      -scorefilename "af2/${src}_scores.sc" \
      2>&1 | tee "logs/04_af2_${src}.log"
  done
  mark_done 04
fi

echo "=== survivors at pae_interaction < $PAE_CUT ==="
python - <<PY
import glob, csv
total = 0
for f in glob.glob("af2/*_scores.sc"):
    rows = [l.split() for l in open(f) if l.strip() and not l.startswith("SCORE:  ")]
    hdr = rows[0]
    try:
        i = hdr.index("pae_interaction")
    except ValueError:
        continue
    n = sum(1 for r in rows[1:] if len(r) > i and float(r[i]) < $PAE_CUT)
    print(f"  {f}: {n}")
    total += n
print(f"  TOTAL: {total}")
if total < 150:
    print("  WARNING: thin funnel. The His scan typically keeps 10-20% of inputs.")
PY

# ------------------------------------------------------- 05 pH scan (CPU)
# No GPU. Run this on all cores in parallel with anything else queued.
if ! stage_done 05; then
  echo "=== 05 pH scan ==="
  python scripts/05_ph_scan.py \
    --complex af2/unbiased/*.pdb af2/hisbias/*.pdb \
    --binder-chain C \
    --target-json target/numbering.json \
    --out results/ph_scan.tsv \
    --combos-out results/ph_combos.tsv \
    2>&1 | tee logs/05_ph_scan.log
  mark_done 05
fi

echo
echo "=== top pH-switch candidates ==="
head -25 results/ph_combos.tsv
echo
echo "Next: build the His combination mutants, repack, re-fold with AF2, and"
echo "confirm pae_interaction < $PAE_CUT still holds at pH 7.4. Anything that"
echo "loses the 7.4 interface is dead regardless of its switch score."
