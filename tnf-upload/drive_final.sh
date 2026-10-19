#!/usr/bin/env bash
# Final driver: bank a corrected slate, then fold in the round-2 backbones.
#
# Written after drive_round2.sh died at stage 03 with "conda: command not
# found" and left the GPU idle for 13.8 h. Two changes follow from that:
#
#   1. `source env.sh` up front. A tmux session is a non-login shell with no
#      conda on PATH.
#   2. Stage completion is signalled by sentinel files, not by polling whether
#      a tmux session is alive. "Session gone" is indistinguishable from
#      "session crashed", which is how the earlier run nearly consumed
#      uncorrected data as though it were finished.
#
# Scope is deliberately trimmed for the remaining budget: ProteinMPNN still
# emits 4 sequences per backbone (unchanged, so the round-1 sequences it
# regenerates are byte-identical and Boltz correctly skips them), but only the
# first 2 per backbone are screened. Changing --num-seq instead would risk
# MPNN emitting different sequences under names Boltz has already folded,
# which would silently score the wrong sequence.
set -uo pipefail
cd "$(dirname "$0")"
# shellcheck disable=SC1091
source env.sh
WORK="$PWD"
LOG="$WORK/logs/drive_final.log"
SENT="$WORK/.sentinels"
mkdir -p "$SENT" "$WORK/logs" "$WORK/results"

CUTOFF="${CUTOFF:-0.285}"
PER_BACKBONE="${PER_BACKBONE:-2}"
NSEQ="${NSEQ:-4}"
JOBS="${JOBS:-10}"
PER_PARENT="${PER_PARENT:-2}"
PER_CLUSTER="${PER_CLUSTER:-6}"
# Stop starting new work after this; a partial slate already on disk beats
# being halfway through a stage when the deadline lands.
BAIL_UTC="${BAIL_UTC:-2026-10-19 09:30:00}"

log() { echo "[$(date -u '+%H:%M:%S')] $*" | tee -a "$LOG"; }
past_bail() { (( $(date -u +%s) >= $(date -u -d "$BAIL_UTC" +%s) )); }

snapshot() {  # keep every slate we produce, newest wins only if larger
  local n; n=$(( $(wc -l < results/submission.csv) - 1 ))
  cp -f results/submission.csv "results/submission_${n}designs_$(date -u +%H%M).csv"
  log "snapshot: ${n}-design slate saved"
}

select_and_bank() {
  local prev; prev=$(( $(wc -l < results/submission.csv) - 1 ))
  cp -f results/submission.csv results/submission_prev.csv
  if venv-phscan/bin/python scripts/07_submit.py \
        --refold results/06_refold.tsv --pred-root boltz/refold \
        --per-parent "$PER_PARENT" --per-cluster "$PER_CLUSTER" >>"$LOG" 2>&1; then
    local new; new=$(( $(wc -l < results/submission.csv) - 1 ))
    if (( new < prev )); then
      log "slate $new < previous $prev -- restoring previous"
      cp -f results/submission_prev.csv results/submission.csv
    else
      log "slate improved: $prev -> $new designs"
      snapshot
    fi
  else
    log "stage 07 FAILED -- previous CSV intact"
    cp -f results/submission_prev.csv results/submission.csv
  fi
}

log "=== final driver start (bail at $BAIL_UTC) ==="

# ---------------------------------------------------------------- 1. wait for the running re-fold
while pgrep -f "06_his_refold" >/dev/null 2>&1; do
  log "waiting for stage 06 on the corrected combos..."
  sleep 120
done
if [[ ! -s results/06_refold.tsv ]]; then
  log "no refold output -- aborting, submission.csv untouched"; exit 1
fi
NSURV=$(cut -f13 results/06_refold.tsv | grep -c '^1$')
log "stage 06 done: $NSURV survivors from $(( $(wc -l < results/ph_combos.tsv) - 1 )) combos"

# ---------------------------------------------------------------- 2. bank the corrected slate NOW
log "stage 07: banking corrected round-1 slate"
select_and_bank
touch "$SENT/round1_banked"

# ---------------------------------------------------------------- 3. round-2 backbones
if past_bail; then log "past bail time -- stopping with current slate"; exit 0; fi

log "stage 03: ProteinMPNN over $(ls rfdiff/*.pdb 2>/dev/null | wc -l) backbones"
if ! conda run -n mpnn python scripts/03_mpnn.py --num-seq "$NSEQ" \
      --pools plain,hisb >>"$LOG" 2>&1; then
  log "stage 03 FAILED -- keeping banked slate"; exit 1
fi
touch "$SENT/stage03"
log "stage 03 done"

log "stage 04: screening $PER_BACKBONE seq/backbone/pool (resumable)"
if ! conda run -n boltz python scripts/04_boltz.py \
      --fasta-dir mpnn/plain/seqs mpnn/hisb/seqs \
      --per-backbone "$PER_BACKBONE" --tag screen >>"$LOG" 2>&1; then
  log "stage 04 FAILED -- keeping banked slate"; exit 1
fi
touch "$SENT/stage04"
log "stage 04 done: $(ls boltz/screen/pred/*/predictions 2>/dev/null | wc -l) predictions"

if past_bail; then log "past bail -- stopping with banked slate"; exit 0; fi

# ---------------------------------------------------------------- 4. new passers only
cp -f results/ph_combos.tsv results/ph_combos_round1.tsv
venv-phscan/bin/python - <<EOF >>"$LOG" 2>&1
import json, pathlib
rows = json.load(open('boltz/screen/scores.json'))
keep = [r['name'] for r in rows
        if (r.get('target_pair_iptm') or 0) >= 0.5 and r['iptm_min'] >= $CUTOFF]
root = pathlib.Path('boltz/screen/pred')
done = {pathlib.Path(x).stem.replace('_model_0','')
        for x in pathlib.Path('results/screen_passers.txt').read_text().split()}
allp, newp = [], []
for n in keep:
    hits = sorted(root.rglob(f"{n}_model_0.pdb"))
    if not hits: continue
    allp.append(str(hits[0]))
    if n not in done: newp.append(str(hits[0]))
pathlib.Path('results/screen_passers.txt').write_text("\n".join(allp) + "\n")
pathlib.Path('results/passers_new.txt').write_text("\n".join(newp) + "\n")
print(f"passers {len(allp)} total, {len(newp)} new")
EOF
NNEW=$(grep -c . results/passers_new.txt 2>/dev/null || echo 0)
log "new passers: $NNEW"

if (( NNEW > 0 )) && ! past_bail; then
  log "stage 05: scanning $NNEW new passers"
  # shellcheck disable=SC2046
  if PATH="$WORK/venv-phscan/bin:$PATH" venv-phscan/bin/python scripts/05_ph_scan.py \
        --complex $(tr '\n' ' ' < results/passers_new.txt) \
        --binder-chain X --jobs "$JOBS" \
        --out results/ph_scan_r2.tsv --combos-out results/ph_combos_r2.tsv \
        >>"$LOG" 2>&1; then
    cp results/ph_combos_round1.tsv results/ph_combos.tsv
    tail -n +2 results/ph_combos_r2.tsv >> results/ph_combos.tsv
    log "combined combos: $(( $(wc -l < results/ph_combos.tsv) - 1 ))"
    touch "$SENT/stage05"
  else
    log "stage 05 FAILED -- keeping banked slate"; exit 1
  fi

  if ! past_bail; then
    log "stage 06: re-folding combined set"
    if conda run -n boltz python scripts/06_his_refold.py \
          --combos results/ph_combos.tsv --parent-tag screen \
          --cutoff "$CUTOFF" --tag refold >>"$LOG" 2>&1; then
      log "stage 06 done: $(cut -f13 results/06_refold.tsv | grep -c '^1$') survivors"
      select_and_bank
    else
      log "stage 06 FAILED -- keeping banked slate"
    fi
  fi
fi

log "=== DONE: submission.csv holds $(( $(wc -l < results/submission.csv) - 1 )) designs ==="
