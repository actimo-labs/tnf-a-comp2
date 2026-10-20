#!/usr/bin/env bash
# TNF-alpha Track 3 -- environment setup for a single A100-80GB GCP instance.
#
# Idempotent and staged. Re-running skips completed stages. Every stage writes
# a marker to .setup/<stage>.done; delete a marker to force a redo.
#
#   bash setup/install.sh                 # Boltz-2 validation path (no PyRosetta)
#   VALIDATOR=af2ig bash setup/install.sh # dl_binder_design path (needs PyRosetta licence)
#   RFD_MODE=conda bash setup/install.sh  # conda instead of docker for RFdiffusion
#
set -euo pipefail

WORKDIR=${WORKDIR:-$HOME/tnf-track3}
VALIDATOR=${VALIDATOR:-boltz}     # boltz | af2ig
RFD_MODE=${RFD_MODE:-docker}      # docker | conda
MARKERS="$WORKDIR/.setup"

mkdir -p "$WORKDIR" "$MARKERS"
cd "$WORKDIR"

log()  { echo -e "\n=== $* ==="; }
done_() { [[ -f "$MARKERS/$1.done" ]]; }
mark() { touch "$MARKERS/$1.done"; }

# ---------------------------------------------------------------- preflight
log "preflight"

if ! command -v nvidia-smi >/dev/null; then
  echo "FATAL: no nvidia-smi. This needs to run on the GPU instance." >&2
  exit 1
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
if (( GPU_MEM < 70000 )); then
  echo "WARNING: ${GPU_MEM}MB VRAM. Plan assumes 80GB." >&2
fi

AVAIL_GB=$(df -BG --output=avail "$WORKDIR" | tail -1 | tr -dc '0-9')
echo "disk available: ${AVAIL_GB}GB"
if (( AVAIL_GB < 60 )); then
  echo "FATAL: need ~60GB (RFdiffusion weights ~3GB, AF2 params 5.3GB, outputs)." >&2
  exit 1
fi

# ------------------------------------------------------------- 00 miniforge
if ! done_ 00_conda; then
  log "00 miniforge"
  if ! command -v conda >/dev/null; then
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    rm Miniforge3-Linux-x86_64.sh
  fi
  mark 00_conda
fi
# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || \
  source "$(conda info --base)/etc/profile.d/conda.sh"

# ------------------------------------------------ 01 pH scan env (CPU only)
# Deliberately a plain venv, isolated from every GPU env. This is the only
# stage that must work -- it is where Track 3 is actually won, and it has no
# CUDA, no PyRosetta, and no licence surface.
if ! done_ 01_phscan; then
  log "01 pH scan venv (CPU)"
  python3 -m venv "$WORKDIR/venv-phscan"
  "$WORKDIR/venv-phscan/bin/pip" install -q --upgrade pip
  "$WORKDIR/venv-phscan/bin/pip" install -q \
    "biopython>=1.83" "propka>=3.5" numpy scipy pandas
  "$WORKDIR/venv-phscan/bin/python" -c "
import Bio, propka
from Bio.PDB.SASA import ShrakeRupley
print('  biopython', Bio.__version__, '/ propka ok')
"
  "$WORKDIR/venv-phscan/bin/propka3" --version 2>/dev/null \
    || echo "  note: propka3 CLI present via venv/bin/propka3"
  mark 01_phscan
fi

# ------------------------------------------------------------ 02 RFdiffusion
if ! done_ 02_rfdiffusion; then
  log "02 RFdiffusion ($RFD_MODE)"
  [[ -d RFdiffusion ]] || git clone -q https://github.com/RosettaCommons/RFdiffusion.git

  # Weights (~3GB). Only the two needed for PPI binder design plus the beta
  # model, which gives non-helical topologies -- worth having, since the
  # default model is overwhelmingly helical and topology diversity matters
  # when in-silico scores predict binding at close to chance.
  mkdir -p RFdiffusion/models && pushd RFdiffusion/models >/dev/null
  base=http://files.ipd.uw.edu/pub/RFdiffusion
  declare -A W=(
    [Base_ckpt.pt]=6f5902ac237024bdd0c176cb93063dc4
    [Complex_base_ckpt.pt]=e29311f6f1bf1af907f9ef9f44b8328b
    [Complex_beta_ckpt.pt]=f572d396fae9206628714fb2ce00f72e
  )
  for f in "${!W[@]}"; do
    if [[ ! -s "$f" ]]; then
      echo "  fetching $f"
      wget -q --show-progress "$base/${W[$f]}/$f" -O "$f" || {
        echo "  FAILED $f -- check network egress from this instance" >&2
        rm -f "$f"; exit 1; }
    fi
  done
  popd >/dev/null

  if [[ "$RFD_MODE" == "docker" ]]; then
    # SE3nv.yml pins CUDA 11.1; a GCP A100 image ships CUDA 12.x drivers and
    # the combination breaks more often than it works. The maintained image
    # sidesteps it.
    command -v docker >/dev/null || { echo "FATAL: docker not installed" >&2; exit 1; }
    docker pull rosettacommons/rfdiffusion:latest
    cat > "$WORKDIR/rfdiffusion.sh" <<'EOF'
#!/usr/bin/env bash
# Wrapper so run_pipeline.sh can call RFdiffusion identically in either mode.
exec docker run --rm --gpus all \
  -v "$PWD:/work" -v "$HOME/tnf-track3/RFdiffusion/models:/app/RFdiffusion/models" \
  -w /work rosettacommons/rfdiffusion:latest \
  python /app/RFdiffusion/scripts/run_inference.py "$@"
EOF
    chmod +x "$WORKDIR/rfdiffusion.sh"
    echo "  wrapper: $WORKDIR/rfdiffusion.sh"
  else
    pushd RFdiffusion >/dev/null
    conda env create -f env/SE3nv.yml -y 2>/dev/null || echo "  (env exists)"
    conda activate SE3nv
    pushd env/SE3Transformer >/dev/null
    pip install -q --no-cache-dir -r requirements.txt
    python setup.py install
    popd >/dev/null
    pip install -q -e .
    conda deactivate
    popd >/dev/null
    cat > "$WORKDIR/rfdiffusion.sh" <<EOF
#!/usr/bin/env bash
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate SE3nv
exec python "$WORKDIR/RFdiffusion/scripts/run_inference.py" "\$@"
EOF
    chmod +x "$WORKDIR/rfdiffusion.sh"
  fi
  mark 02_rfdiffusion
fi

# ------------------------------------------------------------ 03 ProteinMPNN
# Standalone, MIT, no PyRosetta. Includes the soluble weights, which give
# better cell-free expression -- which is how Adaptyv tests these.
if ! done_ 03_mpnn; then
  log "03 ProteinMPNN"
  [[ -d ProteinMPNN ]] || git clone -q https://github.com/dauparas/ProteinMPNN.git
  conda create -n mpnn python=3.10 -y 2>/dev/null || true
  conda activate mpnn
  pip install -q torch --index-url https://download.pytorch.org/whl/cu121
  pip install -q numpy biopython
  python -c "import torch; print('  torch', torch.__version__, 'cuda', torch.cuda.is_available())"
  conda deactivate
  ls ProteinMPNN/soluble_model_weights/ >/dev/null && echo "  soluble weights present"
  mark 03_mpnn
fi

# ------------------------------------------------------------ 04 validator
if ! done_ 04_validator; then
  log "04 validator ($VALIDATOR)"

  if [[ "$VALIDATOR" == "boltz" ]]; then
    # No PyRosetta and no licence question. Cost: the published
    # pae_interaction<10 threshold does not
    # transfer -- calibrate on ipTM against a handful of known binders first.
    conda create -n boltz python=3.11 -y 2>/dev/null || true
    conda activate boltz
    pip install -q boltz
    python -c "import boltz; print('  boltz ok')"
    conda deactivate
    echo "  NOTE: point BOLTZ_MSA_SERVER at your private MSA server before running."
    echo "  NOTE: score binder-vs-target chain pairs ONLY. Global ipTM on this"
    echo "        homotrimer crop is meaningless -- the A-B pair is target-internal"
    echo "        and scores ~0.9 regardless of the binder."

  else
    # Requires PyRosetta. Free for academics; COMMERCIAL USE REQUIRES A LICENCE.
    echo "  PyRosetta licence required for commercial use -- https://els2.comotion.uw.edu/product/pyrosetta"
    read -rp "  Confirm you have a valid PyRosetta licence [y/N] " ok
    [[ "$ok" == "y" ]] || { echo "  aborting; re-run with VALIDATOR=boltz" >&2; exit 1; }

    [[ -d dl_binder_design ]] || git clone -q https://github.com/nrbennet/dl_binder_design.git
    grep -q pyrosetta ~/.condarc 2>/dev/null || cat >> ~/.condarc <<'EOF'
channels:
  - https://conda.graylab.jhu.edu
EOF
    pushd dl_binder_design >/dev/null
    # af2_binder_design.yml is the CUDA 12 env; dl_binder_design.yml is CUDA 11.1
    conda env create -f include/af2_binder_design.yml -y 2>/dev/null || echo "  (env exists)"
    mkdir -p model_weights/params && pushd model_weights/params >/dev/null
    if [[ ! -s params_model_1_ptm.npz ]]; then
      echo "  fetching AF2 params (5.3GB)"
      wget -q --show-progress https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
      tar -xf alphafold_params_2022-12-06.tar params_model_1_ptm.npz
      rm alphafold_params_2022-12-06.tar
    fi
    popd >/dev/null
    [[ -d mpnn_fr/ProteinMPNN ]] || git submodule update --init --recursive
    popd >/dev/null
  fi
  mark 04_validator
fi

# ----------------------------------------------------------------- 04b viz
# Optional: PyMOL for rendering designed complexes (scripts/render_complex.py).
# Not on the design path, so a failure here must not fail the install -- but
# without it the structures in results/ can only be inspected by hand.
if ! done_ 04b_viz; then
  log "04b visualisation (pymol, optional)"
  if conda create -n viz -c conda-forge pymol-open-source pillow -y >/dev/null 2>&1 \
     && conda run -n viz python -c "from pymol import cmd" >/dev/null 2>&1; then
    echo "  pymol OK -- render with:"
    echo "    conda run -n viz pymol -cq scripts/render_complex.py -- \\"
    echo "        results/structures/complex/tnf3_01_complex.pdb figures/tnf3_01 \"Q78H,K16H,E31H\""
    mark 04b_viz
  else
    echo "  pymol install failed -- skipping. Rendering will be unavailable;"
    echo "  nothing else is affected."
  fi
fi

# ----------------------------------------------------------------- 05 env
if ! done_ 05_env; then
  log "05 pipeline env file"
  cat > "$WORKDIR/env.sh" <<EOF
# source this before run_pipeline.sh
export WORKDIR="$WORKDIR"
export RFD_WRAPPER="$WORKDIR/rfdiffusion.sh"
export MPNN="$WORKDIR/ProteinMPNN"
export DLBD="$WORKDIR/dl_binder_design"
export VALIDATOR="$VALIDATOR"
export PHSCAN_PY="$WORKDIR/venv-phscan/bin/python"
export NBACKBONES=\${NBACKBONES:-1000}
export PAE_CUT=\${PAE_CUT:-10}
source "\$(conda info --base)/etc/profile.d/conda.sh"
EOF
  mark 05_env
fi

# ---------------------------------------------------------------- smoke test
log "smoke test"
source "$WORKDIR/env.sh"

echo "-- target prep (fast, no GPU)"
if [[ -f "$WORKDIR/scripts/01_prep_target.py" ]]; then
  "$PHSCAN_PY" "$WORKDIR/scripts/01_prep_target.py" --pdb 9OJO --chains A B \
    --outdir "$WORKDIR/target" || echo "  target prep FAILED -- fix before proceeding"
else
  echo "  scripts/01_prep_target.py not found. Copy the repo files into $WORKDIR first."
fi

echo "-- RFdiffusion (3 backbones, ~2 min)"
mkdir -p "$WORKDIR/smoke"
if [[ -f "$WORKDIR/target/hotspots.txt" ]]; then
  HS=$(cat "$WORKDIR/target/hotspots.txt")
  "$RFD_WRAPPER" \
    inference.output_prefix="$WORKDIR/smoke/test" \
    inference.input_pdb="$WORKDIR/target/tnf_AB.pdb" \
    'contigmap.contigs=[A6-157/0 B6-157/0 70-110]' \
    "ppi.hotspot_res=[${HS}]" \
    inference.num_designs=3 \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    && echo "  RFdiffusion OK ($(ls "$WORKDIR"/smoke/*.pdb 2>/dev/null | wc -l) backbones)" \
    || echo "  RFdiffusion FAILED -- try RFD_MODE=conda, or check docker --gpus all"
else
  echo "  skipped: no target/hotspots.txt"
fi

log "setup complete"
cat <<EOF
  workdir:   $WORKDIR
  validator: $VALIDATOR
  rfdiff:    $RFD_MODE

  source $WORKDIR/env.sh
  bash scripts/run_pipeline.sh

  First RFdiffusion run spends several minutes on "Calculating IGSO3" before
  producing anything. That is normal, not a hang.
EOF
