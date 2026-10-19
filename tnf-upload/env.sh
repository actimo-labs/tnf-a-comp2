# source this before running anything
# Replaces the version install.sh generated, which assumed conda was on PATH.

export WORKDIR="$HOME/tnf-track3"
export RFD_WRAPPER="$WORKDIR/rfdiffusion.sh"
export MPNN="$WORKDIR/ProteinMPNN"
export VALIDATOR="${VALIDATOR:-boltz}"
export PHSCAN_PY="$WORKDIR/venv-phscan/bin/python"
export NBACKBONES="${NBACKBONES:-1000}"
export BINDER_LEN="${BINDER_LEN:-70-110}"

# Absolute path -- `conda info --base` needs conda already on PATH, which it is
# not in a fresh shell unless miniforge wrote to .bashrc.
for _c in "$HOME/miniforge3/etc/profile.d/conda.sh" \
          "$HOME/miniconda3/etc/profile.d/conda.sh" \
          "/opt/conda/etc/profile.d/conda.sh"; do
  if [[ -f "$_c" ]]; then
    # shellcheck disable=SC1090
    source "$_c"
    break
  fi
done
unset _c

if ! command -v conda >/dev/null; then
  echo "WARNING: conda still not found. MPNN and Boltz stages will fail." >&2
fi

mkdir -p "$WORKDIR"/{rfdiff,mpnn,boltz,results,logs}
