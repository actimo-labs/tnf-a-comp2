#!/usr/bin/env python
"""Stage 03 -- ProteinMPNN sequence design on RFdiffusion backbones.

Designs the binder chain only, with the target held as fixed context.

RFdiffusion writes the binder as chain A (renumbered from 1) and collapses the
whole cropped target into chain B, so "design chain A, fix chain B" is exactly
right here -- but it is derived from the output layout, not assumed, and the
script refuses to run on a backbone that doesn't match.

Two pools are produced and kept separate:
  plain  -- soluble weights, no bias. The control.
  hisb   -- soluble weights with a mild His bias, to seed switch positions
            during design rather than only grafting them post hoc at stage 05.

Keeping them apart means stage 05 can show whether the bias actually bought
anything instead of us assuming it did.

Runs under the `mpnn` conda env:
  conda run -n mpnn python scripts/03_mpnn.py --limit 50
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MPNN = ROOT / "ProteinMPNN"
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"


def chains_of(pdb: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[line[21]] = out.get(line[21], 0) + 1
    return out


def stage_inputs(src: Path, dest: Path, limit: int | None) -> list[Path]:
    """Copy backbones into a clean dir; MPNN's parser globs a whole directory.

    Skips anything whose chain layout isn't binder=A/target=B, since a silent
    mismatch here would design the target instead of the binder.
    """
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)

    pdbs = sorted(src.glob("*.pdb"))
    if limit:
        pdbs = pdbs[:limit]

    staged, skipped = [], []
    for p in pdbs:
        ch = chains_of(p)
        if set(ch) != {BINDER_CHAIN, TARGET_CHAIN} or ch[BINDER_CHAIN] < 10:
            skipped.append((p.name, ch))
            continue
        shutil.copy(p, dest / p.name)
        staged.append(p)

    if skipped:
        print(f"skipped {len(skipped)} backbone(s) with unexpected chains:", file=sys.stderr)
        for name, ch in skipped[:5]:
            print(f"  {name}: {ch}", file=sys.stderr)
    return staged


def run(cmd: list[str]) -> None:
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rfdiff", type=Path, default=ROOT / "rfdiff")
    ap.add_argument("--out", type=Path, default=ROOT / "mpnn")
    ap.add_argument("--num-seq", type=int, default=8)
    ap.add_argument("--temp", default="0.1")
    ap.add_argument("--his-bias", type=float, default=1.5)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pools", default="plain,hisb")
    ap.add_argument("--seed", type=int, default=37)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    staged_dir = args.out / "_backbones"
    staged = stage_inputs(args.rfdiff, staged_dir, args.limit)
    if not staged:
        print(f"no usable backbones in {args.rfdiff}", file=sys.stderr)
        return 1
    print(f"staged {len(staged)} backbones")

    parsed = args.out / "parsed.jsonl"
    assigned = args.out / "assigned.jsonl"
    run([sys.executable, MPNN / "helper_scripts/parse_multiple_chains.py",
         f"--input_path={staged_dir}", f"--output_path={parsed}"])
    run([sys.executable, MPNN / "helper_scripts/assign_fixed_chains.py",
         f"--input_path={parsed}", f"--output_path={assigned}",
         f"--chain_list={BINDER_CHAIN}"])

    bias_path = args.out / "bias_his.jsonl"
    bias_path.write_text(json.dumps({"H": args.his_bias}) + "\n")

    for pool in [p.strip() for p in args.pools.split(",") if p.strip()]:
        folder = args.out / pool
        # Clear prior output: MPNN appends into seqs/ and stale fastas from an
        # earlier, smaller backbone set would silently ride along downstream.
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, MPNN / "protein_mpnn_run.py",
               "--jsonl_path", parsed,
               "--chain_id_jsonl", assigned,
               "--out_folder", folder,
               "--num_seq_per_target", args.num_seq,
               "--sampling_temp", args.temp,
               "--seed", args.seed,
               "--batch_size", 1,
               # Free Cys is an automatic reject at stage 07 (disulfide
               # scrambling in cell-free expression), and MPNN emits it
               # readily on these backbones. Cheaper to forbid here than to
               # throw away finished designs later.
               "--omit_AAs", "CX",
               "--use_soluble_model"]
        if pool == "hisb":
            cmd += ["--bias_AA_jsonl", bias_path]
        run(cmd)
        n = len(list((folder / "seqs").glob("*.fa")))
        print(f"pool {pool}: {n} fasta files -> {folder/'seqs'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
