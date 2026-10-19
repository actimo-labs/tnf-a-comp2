#!/usr/bin/env python
"""Crop the TNF target to a shell around the hotspot patch, for RFdiffusion only.

Why: RFdiffusion cost is ~quadratic in total denoised length. The full
two-protomer crop is 304 target residues; with a 70-110 aa binder that is ~382
residues and 1.84 min/design measured on the A100, i.e. 30.7 GPU-h for 1000
backbones against a 29.6 h deadline. A 15 A shell around the hotspots keeps
~162 target residues and cuts that to roughly 12 GPU-h.

This does NOT replace target/tnf_AB.pdb. That file stays the design target and
is what stages 04 and 06 validate against -- the crop exists only to make
backbone generation affordable, and the coordinates are copied verbatim from
it, so no new PDB entry enters the pipeline.

Writes:
  target/tnf_AB_crop.pdb   cropped coordinates, original chain IDs + numbering
  target/crop_contigs.txt  the contigmap.contigs string for the cropped input
  target/crop_report.json  segments kept, counts, hotspot check
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select

ROOT = Path(__file__).resolve().parent.parent


def parse_hotspots(path: Path) -> list[tuple[str, int]]:
    """target/hotspots.txt is RFdiffusion's own format: B31,B32,A87,..."""
    raw = path.read_text().strip()
    out = []
    for tok in raw.replace("[", "").replace("]", "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        out.append((tok[0], int(tok[1:])))
    return out


def segments(nums: list[int], max_gap: int, min_len: int) -> list[tuple[int, int]]:
    """Contiguous runs, merging gaps <= max_gap, dropping runs < min_len.

    Gaps are bridged rather than cut because a backbone with a one- or
    two-residue hole denoises worse than one carrying the extra residues, and
    the extra residues are nearly free.
    """
    if not nums:
        return []
    nums = sorted(nums)
    segs = [[nums[0], nums[0]]]
    for n in nums[1:]:
        if n - segs[-1][1] <= max_gap + 1:
            segs[-1][1] = n
        else:
            segs.append([n, n])
    return [(a, b) for a, b in segs if b - a + 1 >= min_len]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=ROOT / "target/tnf_AB.pdb")
    ap.add_argument("--hotspots", type=Path, default=ROOT / "target/hotspots.txt")
    ap.add_argument("--out", type=Path, default=ROOT / "target/tnf_AB_crop.pdb")
    ap.add_argument("--radius", type=float, default=15.0)
    ap.add_argument("--max-gap", type=int, default=4)
    ap.add_argument("--min-len", type=int, default=5)
    ap.add_argument("--binder-len", default="70-110")
    args = ap.parse_args()

    model = PDBParser(QUIET=True).get_structure("t", str(args.pdb))[0]
    hotspots = parse_hotspots(args.hotspots)

    missing = [h for h in hotspots if (" ", h[1], " ") not in model[h[0]]]
    if missing:
        raise SystemExit(f"hotspot residues absent from {args.pdb}: {missing}")

    atoms = [a for a in model.get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms)
    hot_atoms = [a for ch, num in hotspots for a in model[ch][(" ", num, " ")]]

    near: dict[str, set[int]] = {ch.id: set() for ch in model}
    for a in hot_atoms:
        for nb in ns.search(a.coord, args.radius):
            res = nb.get_parent()
            near[res.get_parent().id].add(res.id[1])

    # Hotspots themselves are in by construction, but assert it rather than assume.
    for ch, num in hotspots:
        near[ch].add(num)

    kept: dict[str, list[tuple[int, int]]] = {}
    for ch_id, nums in near.items():
        segs = segments(sorted(nums), args.max_gap, args.min_len)
        if segs:
            kept[ch_id] = segs

    # A hotspot dropped by min_len would silently remove the thing we are
    # designing against, so re-check after segment filtering.
    lost = [
        (ch, num)
        for ch, num in hotspots
        if not any(a <= num <= b for a, b in kept.get(ch, []))
    ]
    if lost:
        raise SystemExit(f"hotspots lost to segment filtering: {lost}")

    keep_set = {
        (ch_id, n) for ch_id, segs in kept.items() for a, b in segs for n in range(a, b + 1)
    }

    class Keep(Select):
        def accept_residue(self, residue):
            ch = residue.get_parent().id
            return (ch, residue.id[1]) in keep_set

        def accept_atom(self, atom):
            return atom.element != "H" and atom.get_altloc() in (" ", "A")

    io = PDBIO()
    io.set_structure(model)
    io.save(str(args.out), select=Keep())

    # Contig: segments of one chain stay in that chain separated by "/", and
    # "/0" ends a chain. Binder is the trailing length range, as in stage 02.
    parts = []
    for ch_id in sorted(kept):
        parts.append("/".join(f"{ch_id}{a}-{b}" for a, b in kept[ch_id]) + "/0")
    contig = "[" + " ".join(parts) + f" {args.binder_len}]"
    (ROOT / "target/crop_contigs.txt").write_text(contig + "\n")

    n_kept = len(keep_set)
    n_total = sum(1 for _ in model.get_residues())
    report = {
        "source_pdb": str(args.pdb),
        "radius": args.radius,
        "max_gap": args.max_gap,
        "min_len": args.min_len,
        "segments": {k: [list(s) for s in v] for k, v in kept.items()},
        "residues_kept": n_kept,
        "residues_total": n_total,
        "hotspots": [f"{c}{n}" for c, n in hotspots],
        "contig": contig,
    }
    (ROOT / "target/crop_report.json").write_text(json.dumps(report, indent=2))

    print(f"kept {n_kept}/{n_total} target residues at {args.radius} A")
    for ch_id in sorted(kept):
        spans = ", ".join(f"{a}-{b}" for a, b in kept[ch_id])
        total = sum(b - a + 1 for a, b in kept[ch_id])
        print(f"  chain {ch_id}: {total:3d} res  [{spans}]")
    print(f"contig: {contig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
