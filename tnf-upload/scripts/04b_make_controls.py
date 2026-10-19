#!/usr/bin/env python
"""Build calibration controls for the stage-04 Boltz cutoff.

Handoff.md is right that AF2's `pae_interaction < 10` does not transfer to a
Boltz pair-ipTM, so the cutoff has to be set from data rather than imported.

Negatives are composition- and length-matched shuffles of real designs. A
shuffle keeps amino-acid content and chain length but destroys the designed
surface, so whatever ipTM it still earns is the floor this scoring setup gives
to a sequence that cannot bind. The cutoff goes above that floor.

Caveat worth stating plainly: this calibrates the *null*, not the positive end.
Without a known-binder positive the cutoff is "clearly better than noise",
which is weaker than "comparable to a real binder". Add --positives with a
FASTA of known TNF binders to anchor the top end too.
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def read_mpnn_fastas(fasta_dir: Path, per_backbone: int) -> list[tuple[str, str]]:
    out = []
    for fa in sorted(fasta_dir.glob("*.fa")):
        recs, name, buf = [], None, []
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        for i, (_, seq) in enumerate(recs[1 : 1 + per_backbone]):
            out.append((f"{fa.stem}_d{i}", seq.split("/")[0]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", type=Path, default=ROOT / "mpnn/plain/seqs")
    ap.add_argument("--out", type=Path, default=ROOT / "boltz/controls.fasta")
    ap.add_argument("--n-negatives", type=int, default=12)
    ap.add_argument("--per-backbone", type=int, default=1)
    ap.add_argument("--positives", type=Path, default=None)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    designs = read_mpnn_fastas(args.fasta_dir, args.per_backbone)
    if not designs:
        raise SystemExit(f"no sequences in {args.fasta_dir}")
    rng.shuffle(designs)

    lines = []
    for i, (name, seq) in enumerate(designs[: args.n_negatives]):
        chars = list(seq)
        rng.shuffle(chars)
        lines.append(f">NEG_shuf_{i:02d}_{name}")
        lines.append("".join(chars))

    n_pos = 0
    if args.positives and args.positives.exists():
        name, buf = None, []
        for line in args.positives.read_text().splitlines():
            if line.startswith(">"):
                if name:
                    lines += [f">POS_{name}", "".join(buf)]
                    n_pos += 1
                name, buf = line[1:].strip().replace(" ", "_"), []
            elif line.strip():
                buf.append(line.strip())
        if name:
            lines += [f">POS_{name}", "".join(buf)]
            n_pos += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.n_negatives} negatives + {n_pos} positives -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
