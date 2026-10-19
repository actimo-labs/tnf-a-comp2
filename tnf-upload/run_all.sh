#!/usr/bin/env bash
# TNF-alpha Track 3 pipeline.
#
# Replaces scripts/run_pipeline.sh, which was written against dl_binder_design
# (ProteinMPNN + AF2 initial-guess, needing a PyRosetta licence). That toolchain
# is not installed and is not being installed; this drives ProteinMPNN and
# Boltz-2 directly instead.
#
# Stages are independent and resumable -- run one at a time rather than firing
# the lot. Stage 02 and the Boltz stages are the GPU consumers; 05 and 07 are
# CPU-only and can run while the GPU is busy.
#
#   ./run_all.sh 02          # RFdiffusion backbones (tmux, long)
#   ./run_all.sh 03          # ProteinMPNN sequences
#   ./run_all.sh 04          # Boltz-2 screen
#   ./run_all.sh 05          # pH linkage scan (CPU)
#   ./run_all.sh 06          # re-fold His variants
#   ./run_all.sh 07          # submission CSV
set -euo pipefail

cd "$(dirname "$0")"
# shellcheck disable=SC1091
source env.sh

STAGE="${1:-}"
NDESIGNS="${NDESIGNS:-500}"
CUTOFF="${CUTOFF:-}"
TAG="${TAG:-screen}"

need_cutoff() {
  if [[ -z "$CUTOFF" ]]; then
    echo "CUTOFF is unset. Calibrate it from boltz/calib/scores.json first:" >&2
    echo "  the negatives (NEG_shuf_*) give the floor; set CUTOFF above it." >&2
    exit 1
  fi
}

case "$STAGE" in
  02)
    # Cropped target: ~1.1 min/backbone vs ~1.8 for the full two-protomer
    # crop. See scripts/02_crop_target.py for why that matters.
    [[ -f target/tnf_AB_crop.pdb ]] || \
      venv-phscan/bin/python scripts/02_crop_target.py
    tmux new-session -d -s rfdiff \
      "cd '$PWD' && ./rfdiffusion.sh \
        inference.output_prefix=/work/rfdiff/tnf \
        inference.input_pdb=/work/target/tnf_AB_crop.pdb \
        \"contigmap.contigs=\$(cat target/crop_contigs.txt)\" \
        'ppi.hotspot_res=[$(cat target/hotspots.txt)]' \
        inference.num_designs=$NDESIGNS \
        denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
        hydra.run.dir=/work/.hydra > logs/02_rfdiffusion.log 2>&1"
    echo "started in tmux 'rfdiff'; watch: ls rfdiff/tnf_*.pdb | wc -l"
    ;;

  03)
    conda run -n mpnn python scripts/03_mpnn.py "${@:2}"
    ;;

  04)
    [[ -f target/tnf_AB.cif ]] || \
      conda run -n boltz python scripts/make_template_cif.py
    conda run -n boltz python scripts/04_boltz.py \
      --fasta-dir mpnn/plain/seqs mpnn/hisb/seqs --tag "$TAG" "${@:2}"
    ;;

  05)
    # propka3 lives in the venv; the script shells out to it by name.
    PATH="$PWD/venv-phscan/bin:$PATH" \
      venv-phscan/bin/python scripts/05_ph_scan.py \
        --binder-chain X \
        --complex boltz/"$TAG"/pred/*/predictions/*/*_model_0.pdb \
        --out results/ph_scan.tsv --combos-out results/ph_combos.tsv "${@:2}"
    ;;

  06)
    need_cutoff
    conda run -n boltz python scripts/06_his_refold.py \
      --parent-tag "$TAG" --cutoff "$CUTOFF" "${@:2}"
    ;;

  07)
    venv-phscan/bin/python scripts/07_submit.py "${@:2}"
    ;;

  *)
    sed -n '2,20p' "$0"
    exit 1
    ;;
esac
