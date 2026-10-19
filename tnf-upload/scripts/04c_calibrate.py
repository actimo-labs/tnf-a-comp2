#!/usr/bin/env python
"""Pick the stage-04 pair-ipTM cutoff from the calibration run.

Handoff.md: no published threshold transfers here. AF2's `pae_interaction < 10`
is a different model, a different metric and a different scale, so the number
has to come from this setup's own behaviour.

What this reads: scores.json from a stage-04 run that included the controls
built by 04b_make_controls.py. Negatives are composition-matched shuffles, so
their scores are what this pipeline hands a sequence that cannot bind. The
cutoff goes above that floor.

Two numbers are reported:
  neg_max     -- the highest score any shuffle achieved
  neg_p95     -- the 95th percentile of the shuffle distribution

Default cutoff is max(neg_max, neg_p95) plus a margin. Using neg_max rather
than the mean is deliberate: a cutoff at the mean of the null admits half the
noise, and with ~50% of designs expected to fail anyway the cost of a loose
cutoff is a stage-06 queue full of designs that were never binding.
"""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path,
                    default=ROOT / "boltz/calib/scores.json")
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--min-target-pair", type=float, default=0.5,
                    help="reject predictions where the target protomers "
                         "failed to assemble; the epitope is then absent")
    args = ap.parse_args()

    if not args.scores.exists():
        raise SystemExit(f"no scores at {args.scores}")
    rows = json.loads(args.scores.read_text())

    bad_target = [r for r in rows
                  if (r.get("target_pair_iptm") or 0) < args.min_target_pair]
    rows = [r for r in rows
            if (r.get("target_pair_iptm") or 0) >= args.min_target_pair]

    negs = [r for r in rows if r["name"].startswith("NEG_")]
    poss = [r for r in rows if r["name"].startswith("POS_")]
    designs = [r for r in rows if not r["name"].startswith(("NEG_", "POS_"))]

    if bad_target:
        print(f"WARNING: {len(bad_target)} prediction(s) dropped -- target "
              f"protomers did not assemble (pair ipTM < {args.min_target_pair})")
    if not negs:
        raise SystemExit("no NEG_ controls in scores -- rerun stage 04 with "
                         "--extra-seqs boltz/controls.fasta")

    nv = sorted(r["iptm_min"] for r in negs)
    dv = sorted((r["iptm_min"] for r in designs), reverse=True)

    neg_max = nv[-1]
    neg_p95 = nv[int(0.95 * (len(nv) - 1))]
    cutoff = round(max(neg_max, neg_p95) + args.margin, 3)

    print(f"negatives  n={len(nv):3d}  "
          f"min={nv[0]:.3f} median={statistics.median(nv):.3f} "
          f"p95={neg_p95:.3f} max={neg_max:.3f}")
    if poss:
        pv = sorted(r["iptm_min"] for r in poss)
        print(f"positives  n={len(pv):3d}  "
              f"min={pv[0]:.3f} median={statistics.median(pv):.3f} "
              f"max={pv[-1]:.3f}")
    else:
        print("positives  none -- cutoff is 'above noise', NOT 'as good as a "
              "known binder'. Treat it as a floor.")
    if dv:
        print(f"designs    n={len(dv):3d}  "
              f"min={dv[-1]:.3f} median={statistics.median(dv):.3f} "
              f"max={dv[0]:.3f}")

    print(f"\nsuggested CUTOFF = {cutoff}")
    for c in (cutoff, cutoff + 0.05, cutoff + 0.10):
        n_pass = sum(1 for v in dv if v >= c)
        rate = 100.0 * n_pass / len(dv) if dv else 0.0
        print(f"  cutoff {c:.3f}: {n_pass}/{len(dv)} designs pass ({rate:.0f}%)")

    print("\nHandoff rule of thumb: under ~150 survivors going into the His "
          "scan means the funnel is too thin. Raise backbone count rather "
          "than lowering this cutoff.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
