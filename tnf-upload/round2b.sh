#!/usr/bin/env bash
# Round 2, linear. No waiting on other jobs -- every previous driver failed in
# its wait loop (tmux-session polling, then pgrep matching the launcher's own
# command line). This script only runs its own stages, in order, and each one
# either succeeds or stops the script.
#
# submission.csv currently holds a validated 20-design corrected slate. It is
# replaced only if this run produces a strictly larger one, and a copy is kept
# either way.
set -uo pipefail
cd "$(dirname "$0")"
source env.sh                    # conda is not on PATH in a tmux shell
LOG="logs/round2b.log"

CUTOFF=0.285
PER_BACKBONE=1                   # screen 2 seq/backbone/pool; MPNN still makes 4
NSEQ=4                           # unchanged, so round-1 sequences regenerate identically
JOBS=10
BAIL_UTC="2026-10-19 10:30:00"   # 21:30 AEDT

log() { echo "[$(date -u '+%H:%M:%S')] $*" | tee -a "$LOG"; }
past_bail() { (( $(date -u +%s) >= $(date -u -d "$BAIL_UTC" +%s) )); }
keep() { cp -f results/submission.csv "results/submission_$1.csv"; }

keep "pre_round2"
BASE=$(( $(wc -l < results/submission.csv) - 1 ))
log "=== round 2 start; banked slate = $BASE designs; bail $BAIL_UTC ==="

log "stage 03: SKIPPED -- 383 fastas already present from the earlier run"

log "stage 04: screening (resumable; round-1 designs skipped)"
conda run -n boltz python scripts/04_boltz.py \
  --fasta-dir mpnn/plain/seqs mpnn/hisb/seqs \
  --per-backbone "$PER_BACKBONE" --tag screen >>"$LOG" 2>&1 \
  || { log "stage 04 FAILED -- slate untouched"; exit 1; }
log "stage 04 done: $(ls boltz/screen/pred/*/predictions | wc -l) predictions"

past_bail && { log "past bail -- stopping with banked slate"; exit 0; }

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
(( NNEW == 0 )) && { log "no new passers -- keeping banked slate"; exit 0; }

past_bail && { log "past bail -- stopping with banked slate"; exit 0; }

log "stage 05: scanning $NNEW new passers"
PATH="$PWD/venv-phscan/bin:$PATH" venv-phscan/bin/python scripts/05_ph_scan.py \
  --complex $(tr '\n' ' ' < results/passers_new.txt) \
  --binder-chain X --jobs "$JOBS" \
  --out results/ph_scan_r2.tsv --combos-out results/ph_combos_r2.tsv \
  >>"$LOG" 2>&1 || { log "stage 05 FAILED -- slate untouched"; exit 1; }
cp results/ph_combos_round1.tsv results/ph_combos.tsv
tail -n +2 results/ph_combos_r2.tsv >> results/ph_combos.tsv
log "combined combos: $(( $(wc -l < results/ph_combos.tsv) - 1 ))"

past_bail && { log "past bail -- stopping with banked slate"; exit 0; }

log "stage 06: re-folding combined set (cached variants are free)"
conda run -n boltz python scripts/06_his_refold.py \
  --combos results/ph_combos.tsv --parent-tag screen \
  --cutoff "$CUTOFF" --tag refold >>"$LOG" 2>&1 \
  || { log "stage 06 FAILED -- slate untouched"; exit 1; }
log "stage 06 done: $(cut -f13 results/06_refold.tsv | grep -c '^1$') survivors"

log "stage 07: selection"
cp -f results/submission.csv results/submission_prev.csv
if venv-phscan/bin/python scripts/07_submit.py \
      --refold results/06_refold.tsv --pred-root boltz/refold \
      --per-parent 2 --per-cluster 8 >>"$LOG" 2>&1; then
  NEW=$(( $(wc -l < results/submission.csv) - 1 ))
  if (( NEW < BASE )); then
    log "new slate $NEW < banked $BASE -- restoring banked"
    cp -f results/submission_prev.csv results/submission.csv
  else
    log "slate: $BASE -> $NEW designs"
    keep "round2_final"
  fi
else
  log "stage 07 FAILED -- restoring banked slate"
  cp -f results/submission_prev.csv results/submission.csv
fi

log "=== DONE: submission.csv holds $(( $(wc -l < results/submission.csv) - 1 )) designs ==="
