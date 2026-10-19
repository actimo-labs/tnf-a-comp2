#!/usr/bin/env bash
# Unattended driver: stages 03 -> 07, no operator required.
#
# Why this exists: the stage hand-offs were being done by hand. Over a long
# overnight run against a hard deadline that is a single point of failure --
# if the controlling session stops, the GPU job finishes and then nothing
# happens. This script makes the remaining pipeline self-driving, so the only
# thing that has to stay up is the VM.
#
# Run it detached:
#   tmux new-session -d -s drive './drive_overnight.sh'
#   tail -f logs/drive.log
#
# Safe to run while RFdiffusion is still generating: it waits for a backbone
# target, and stage 04 is resumable, so re-running re-folds nothing.
#
# The existing results/submission.csv is treated as a floor. It is replaced
# only if a later pass yields strictly more designs.
set -uo pipefail

cd "$(dirname "$0")"
WORK="$PWD"
LOG="$WORK/logs/drive.log"
mkdir -p "$WORK/logs" "$WORK/results"

TARGET_BB="${TARGET_BB:-230}"       # stop waiting for backbones at this count
BB_DEADLINE_MIN="${BB_DEADLINE_MIN:-240}"  # ...or after this long, whichever first
CUTOFF="${CUTOFF:-0.285}"
NSEQ="${NSEQ:-4}"
JOBS="${JOBS:-10}"
STOP_RFDIFF="${STOP_RFDIFF:-1}"     # free the GPU for Boltz once generating is done

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

count_bb() { ls "$WORK"/rfdiff/tnf_*.pdb 2>/dev/null | wc -l; }

log "=== driver start: target ${TARGET_BB} backbones, cutoff ${CUTOFF} ==="

# ---------------------------------------------------------------- wait for backbones
waited=0
while :; do
  n=$(count_bb)
  if (( n >= TARGET_BB )); then
    log "backbone target reached: $n"; break
  fi
  if (( waited >= BB_DEADLINE_MIN )); then
    log "backbone wait capped at ${BB_DEADLINE_MIN} min with $n backbones; proceeding"
    break
  fi
  if ! tmux has-session -t rfdiff 2>/dev/null; then
    log "rfdiff session gone with $n backbones; proceeding"
    break
  fi
  sleep 300; waited=$((waited+5))
  (( waited % 30 == 0 )) && log "waiting for backbones: $n/$TARGET_BB (${waited}m)"
done

if [[ "$STOP_RFDIFF" == "1" ]] && tmux has-session -t rfdiff 2>/dev/null; then
  log "stopping RFdiffusion to give the GPU to Boltz"
  tmux kill-session -t rfdiff 2>/dev/null
  docker ps -q | xargs -r docker stop >/dev/null 2>&1
  # Wait for the GPU to actually be released before loading Boltz onto it.
  for _ in $(seq 1 60); do
    used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    (( used < 2000 )) && break
    sleep 5
  done
  log "GPU free"
fi

NBB=$(count_bb)
log "proceeding with $NBB backbones"

# ---------------------------------------------------------------- 03 ProteinMPNN
log "stage 03: ProteinMPNN ($NSEQ seq/backbone, plain + hisb)"
if ! conda run -n mpnn python scripts/03_mpnn.py --num-seq "$NSEQ" \
      --pools plain,hisb >>"$LOG" 2>&1; then
  log "STAGE 03 FAILED -- stopping"; exit 1
fi
log "stage 03 done: $(ls mpnn/plain/seqs/*.fa 2>/dev/null | wc -l) backbones x $NSEQ x 2 pools"

# ---------------------------------------------------------------- 04 Boltz screen
log "stage 04: Boltz screen (resumable; already-folded designs are skipped)"
if ! conda run -n boltz python scripts/04_boltz.py \
      --fasta-dir mpnn/plain/seqs mpnn/hisb/seqs \
      --per-backbone "$NSEQ" --tag screen >>"$LOG" 2>&1; then
  log "STAGE 04 FAILED -- stopping"; exit 1
fi
log "stage 04 done: $(ls boltz/screen/pred/*/predictions 2>/dev/null | wc -l) predictions"

# ---------------------------------------------------------------- select passers
venv-phscan/bin/python - <<EOF >>"$LOG" 2>&1
import json, pathlib
rows = json.load(open('boltz/screen/scores.json'))
keep = [r['name'] for r in rows
        if (r.get('target_pair_iptm') or 0) >= 0.5 and r['iptm_min'] >= $CUTOFF]
root = pathlib.Path('boltz/screen/pred')
paths = []
for n in keep:
    hits = sorted(root.rglob(f"{n}_model_0.pdb"))
    if hits:
        paths.append(str(hits[0]))
pathlib.Path('results/screen_passers.txt').write_text("\n".join(paths) + "\n")
print(f"{len(paths)} passers of {len(rows)} scored")
EOF
NPASS=$(wc -l < results/screen_passers.txt)
log "passers at cutoff $CUTOFF: $NPASS"
if (( NPASS == 0 )); then log "no passers -- stopping"; exit 1; fi

# ---------------------------------------------------------------- 05 pH scan (CPU)
log "stage 05: pH linkage scan on $NPASS passers (${JOBS} jobs)"
# shellcheck disable=SC2046
if ! PATH="$WORK/venv-phscan/bin:$PATH" venv-phscan/bin/python scripts/05_ph_scan.py \
      --complex $(tr '\n' ' ' < results/screen_passers.txt) \
      --binder-chain X --jobs "$JOBS" \
      --out results/ph_scan.tsv --combos-out results/ph_combos.tsv >>"$LOG" 2>&1; then
  log "STAGE 05 FAILED -- stopping"; exit 1
fi
NCOMBO=$(( $(wc -l < results/ph_combos.tsv) - 1 ))
log "stage 05 done: $NCOMBO combinations >= 2.8 kcal/mol"
if (( NCOMBO <= 0 )); then log "no combinations -- stopping"; exit 1; fi

# ---------------------------------------------------------------- 06 re-fold
log "stage 06: re-folding $NCOMBO His variants"
if ! conda run -n boltz python scripts/06_his_refold.py \
      --combos results/ph_combos.tsv --parent-tag screen \
      --cutoff "$CUTOFF" --tag refold >>"$LOG" 2>&1; then
  log "STAGE 06 FAILED -- stopping"; exit 1
fi
NSURV=$(cut -f13 results/06_refold.tsv | grep -c '^1$')
log "stage 06 done: $NSURV variants keep their pH 7.4 interface"

# ---------------------------------------------------------------- 07 submission
# Keep the previous CSV as a floor; only replace it if this pass did better.
PREV=0
[[ -f results/submission.csv ]] && PREV=$(( $(wc -l < results/submission.csv) - 1 ))
cp -f results/submission.csv "results/submission_prev.csv" 2>/dev/null || true

log "stage 07: developability, clustering, submission CSV"
if ! venv-phscan/bin/python scripts/07_submit.py \
      --refold results/06_refold.tsv --pred-root boltz/refold >>"$LOG" 2>&1; then
  log "STAGE 07 FAILED -- previous CSV left in place"; exit 1
fi

NEW=$(( $(wc -l < results/submission.csv) - 1 ))
if (( NEW < PREV )); then
  log "new slate ($NEW) smaller than previous ($PREV) -- restoring previous"
  cp -f results/submission_prev.csv results/submission.csv
  NEW=$PREV
fi

log "=== DONE: results/submission.csv holds $NEW designs (was $PREV) ==="
log "review: results/submission_report.tsv"
