#!/usr/bin/env bash
# Round-2 driver: cap generation, screen the new backbones, merge, finish.
#
# Differs from drive_overnight.sh in one important way: stage 05 is run only on
# passers that are NEW this round, and its combos are concatenated onto the
# round-1 results. Re-scanning all ~1000 passers would cost ~4 h against a
# budget that does not have it, and the scan is per-design independent so the
# concatenation is exact rather than an approximation.
#
# Waits for the corrected round-1 pH scan to land before touching its output.
#
#   tmux new-session -d -s drive2 './drive_round2.sh'
set -uo pipefail
cd "$(dirname "$0")"
# A tmux session starts from a non-login shell with no conda on PATH. Without
# this every `conda run` in the script dies with "conda: command not found" --
# which is exactly how the first round-2 attempt failed, silently, leaving the
# GPU idle for hours. env.sh sources conda by absolute path.
# shellcheck disable=SC1091
source env.sh
WORK="$PWD"
LOG="$WORK/logs/drive2.log"
mkdir -p "$WORK/logs" "$WORK/results"

TARGET_BB="${TARGET_BB:-150}"        # round-2 backbones (r2_*.pdb)
HARD_STOP_UTC="${HARD_STOP_UTC:-}"   # optional ISO time to stop generating
CUTOFF="${CUTOFF:-0.285}"
NSEQ="${NSEQ:-4}"
JOBS="${JOBS:-10}"
PER_PARENT="${PER_PARENT:-2}"
PER_CLUSTER="${PER_CLUSTER:-6}"

log() { echo "[$(date -u '+%H:%M:%S')] $*" | tee -a "$LOG"; }
count_r2() { ls "$WORK"/rfdiff/r2_*.pdb 2>/dev/null | wc -l; }

log "=== round 2: target ${TARGET_BB} new backbones ==="

# ---- 1. wait for the corrected round-1 scan (it owns results/ph_combos.tsv)
while tmux has-session -t phscan 2>/dev/null; do
  log "waiting for corrected round-1 pH scan..."
  sleep 120
done
if [[ ! -s results/ph_combos.tsv ]]; then
  log "round-1 combos missing -- aborting rather than overwriting good state"
  exit 1
fi
cp results/ph_combos.tsv results/ph_combos_round1.tsv
cp results/screen_passers.txt results/passers_round1.txt
R1=$(( $(wc -l < results/ph_combos_round1.tsv) - 1 ))
log "round-1 corrected combos banked: $R1"

# ---- 2. cap generation
while :; do
  n=$(count_r2)
  (( n >= TARGET_BB )) && { log "round-2 target reached: $n"; break; }
  tmux has-session -t rfdiff2 2>/dev/null || { log "rfdiff2 gone at $n"; break; }
  if [[ -n "$HARD_STOP_UTC" ]] && \
     (( $(date -u +%s) >= $(date -u -d "$HARD_STOP_UTC" +%s) )); then
    log "hard stop reached with $n backbones"; break
  fi
  sleep 180
done

tmux kill-session -t rfdiff2 2>/dev/null
docker ps -q | xargs -r docker stop >/dev/null 2>&1
for _ in $(seq 1 60); do
  used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
  (( used < 2000 )) && break
  sleep 5
done
log "generation stopped; GPU free; r2 backbones: $(count_r2)"

# ---- 3. sequences over every backbone (round 1 regenerates identically)
log "stage 03: ProteinMPNN over $(ls rfdiff/*.pdb | wc -l) backbones"
conda run -n mpnn python scripts/03_mpnn.py --num-seq "$NSEQ" --pools plain,hisb \
  >>"$LOG" 2>&1 || { log "STAGE 03 FAILED"; exit 1; }

# ---- 4. screen (resumable: round-1 designs are already folded and skipped)
log "stage 04: Boltz screen, resumable"
conda run -n boltz python scripts/04_boltz.py \
  --fasta-dir mpnn/plain/seqs mpnn/hisb/seqs \
  --per-backbone "$NSEQ" --tag screen >>"$LOG" 2>&1 \
  || { log "STAGE 04 FAILED"; exit 1; }

# ---- 5. passers, split into new vs already scanned
venv-phscan/bin/python - <<EOF >>"$LOG" 2>&1
import json, pathlib
rows = json.load(open('boltz/screen/scores.json'))
keep = [r['name'] for r in rows
        if (r.get('target_pair_iptm') or 0) >= 0.5 and r['iptm_min'] >= $CUTOFF]
root = pathlib.Path('boltz/screen/pred')
old = set()
for line in pathlib.Path('results/passers_round1.txt').read_text().split():
    old.add(pathlib.Path(line).stem.replace('_model_0', ''))
allp, newp = [], []
for n in keep:
    hits = sorted(root.rglob(f"{n}_model_0.pdb"))
    if not hits: continue
    allp.append(str(hits[0]))
    if n not in old: newp.append(str(hits[0]))
pathlib.Path('results/screen_passers.txt').write_text("\n".join(allp) + "\n")
pathlib.Path('results/passers_new.txt').write_text("\n".join(newp) + "\n")
print(f"passers total {len(allp)}, new this round {len(newp)}")
EOF
NNEW=$(wc -l < results/passers_new.txt)
log "passers: $(wc -l < results/screen_passers.txt) total, $NNEW new"

# ---- 6. scan only the new ones, then merge onto round 1
if (( NNEW > 0 )); then
  log "stage 05: scanning $NNEW new passers"
  PATH="$WORK/venv-phscan/bin:$PATH" venv-phscan/bin/python scripts/05_ph_scan.py \
    --complex $(tr '\n' ' ' < results/passers_new.txt) \
    --binder-chain X --jobs "$JOBS" \
    --out results/ph_scan_r2.tsv --combos-out results/ph_combos_r2.tsv \
    >>"$LOG" 2>&1 || { log "STAGE 05 FAILED"; exit 1; }
  cp results/ph_combos_round1.tsv results/ph_combos.tsv
  tail -n +2 results/ph_combos_r2.tsv >> results/ph_combos.tsv
  tail -n +2 results/ph_scan_r2.tsv  >> results/ph_scan.tsv
else
  log "no new passers; keeping round-1 combos"
  cp results/ph_combos_round1.tsv results/ph_combos.tsv
fi
NCOMBO=$(( $(wc -l < results/ph_combos.tsv) - 1 ))
log "combined combinations: $NCOMBO (round 1 had $R1)"

# ---- 7. re-fold and select
log "stage 06: re-folding $NCOMBO variants (cached ones are free)"
conda run -n boltz python scripts/06_his_refold.py \
  --combos results/ph_combos.tsv --parent-tag screen \
  --cutoff "$CUTOFF" --tag refold >>"$LOG" 2>&1 \
  || { log "STAGE 06 FAILED"; exit 1; }
NSURV=$(cut -f13 results/06_refold.tsv | grep -c '^1$')
log "stage 06: $NSURV survivors"

PREV=$(( $(wc -l < results/submission.csv) - 1 ))
cp -f results/submission.csv results/submission_prev.csv

log "stage 07: selection"
venv-phscan/bin/python scripts/07_submit.py \
  --refold results/06_refold.tsv --pred-root boltz/refold \
  --per-parent "$PER_PARENT" --per-cluster "$PER_CLUSTER" >>"$LOG" 2>&1 \
  || { log "STAGE 07 FAILED -- previous CSV intact"; exit 1; }

NEW=$(( $(wc -l < results/submission.csv) - 1 ))
if (( NEW < PREV )); then
  log "new slate $NEW < previous $PREV -- restoring previous"
  cp -f results/submission_prev.csv results/submission.csv
  NEW=$PREV
fi
log "=== DONE: submission.csv holds $NEW designs (was $PREV) ==="
