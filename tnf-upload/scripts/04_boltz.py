#!/usr/bin/env python
"""Stage 04 -- Boltz-2 validation of designed binders against TNF-alpha.

Scoring rule (CLAUDE.md non-negotiable): the target is a homotrimer crop, so
global ipTM is meaningless -- the target-internal A-B pair scores ~0.9 no
matter what the binder does, and including it inflates every design. This
script reads `pair_chains_iptm` and keeps ONLY the binder-vs-target entries,
dropping A-B on the floor. The reported score is the *minimum* of the two
binder-target pair values, not the mean: a binder that grips protomer A and
ignores B is not binding the receptor groove we designed against.

The target is supplied as two separate chains carrying the full 157-residue
assay construct (TNA-H4211, verified tag-free and identical to
target/tnf_AB.pdb over 6-157), so the chain pairs are real chain pairs.

Single-sequence mode by default: the target MSA would be identical across
every design, and the binders are de novo so they have no meaningful MSA.

  conda run -n boltz python scripts/04_boltz.py --fasta-dir mpnn/plain/seqs
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Mature human TNF-alpha, Adaptyv assay construct TNA-H4211, 157 aa, tag-free.
# Verified character-identical to target/tnf_AB.pdb over residues 6-157.
TNF = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYS"
       "QVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYL"
       "GGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")

BINDER_ID = "X"
TARGET_IDS = ("A", "B")


def read_binders(fasta_dir: Path, per_backbone: int | None) -> dict[str, str]:
    """MPNN fastas: record 0 is the input backbone, records 1..n are designs."""
    out: dict[str, str] = {}
    for fa in sorted(fasta_dir.glob("*.fa")):
        stem = fa.stem
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

        designs = recs[1:]
        if per_backbone:
            designs = designs[:per_backbone]
        for hdr, seq in designs:
            m = re.search(r"sample=(\d+)", hdr)
            idx = m.group(1) if m else str(len(out))
            # Chain A only is emitted by stage 03, but guard anyway: a "/" would
            # mean the target leaked into the design record.
            seq = seq.split("/")[0]
            out[f"{stem}_s{idx}"] = seq
    return out


def write_yaml(path: Path, binder: str, target: str, template: Path | None) -> None:
    """One complex: binder + two target protomers, target pinned to the crystal.

    The template matters more than it looks. Without it, single-sequence Boltz
    folds each protomer fine (ptm ~0.85) but does NOT dock them to each other
    (A-B ipTM ~0.29) -- and since the receptor-binding groove sits at the
    protomer-protomer interface, the epitope we designed against simply is not
    present in the prediction. Pinning both target chains to target/tnf_AB.pdb
    restores the real geometry and leaves Boltz predicting only the binder pose,
    which is the thing we actually want measured.
    """
    chains = [(BINDER_ID, binder)] + [(cid, target) for cid in TARGET_IDS]
    lines = ["version: 1", "sequences:"]
    for cid, seq in chains:
        lines += ["  - protein:",
                  f"      id: {cid}",
                  f"      sequence: {seq}",
                  "      msa: empty"]
    if template:
        lines += ["templates:",
                  f"  - cif: {template}",
                  f"    chain_id: [{', '.join(TARGET_IDS)}]",
                  f"    template_id: [{', '.join(TARGET_IDS)}]"]
    path.write_text("\n".join(lines) + "\n")


def parse_confidence(pred_root: Path, name: str) -> dict | None:
    """Pull binder-target pair ipTM out of a Boltz confidence JSON."""
    hits = sorted(pred_root.rglob(f"confidence_{name}_model_0.json"))
    if not hits:
        return None
    data = json.loads(hits[0].read_text())

    # Boltz keys pair_chains_iptm by chain *index* as a string, in the order
    # the chains appear in the YAML: binder is 0, targets are 1 and 2.
    pair = data.get("pair_chains_iptm") or {}
    ids = [BINDER_ID, *TARGET_IDS]

    def pair_val(i: int, j: int) -> float | None:
        """Symmetrised: the matrix is not exactly symmetric (0.187 vs 0.215)."""
        got = []
        for a, b in ((i, j), (j, i)):
            try:
                got.append(float(pair[str(a)][str(b)]))
            except (KeyError, TypeError, ValueError):
                continue
        return sum(got) / len(got) if got else None

    vals = {tgt: pair_val(0, ids.index(tgt)) for tgt in TARGET_IDS}
    good = [v for v in vals.values() if v is not None]
    if not good:
        return None

    chains_ptm = data.get("chains_ptm") or {}
    return {
        "name": name,
        "pair_iptm": vals,
        "iptm_min": min(good),        # the score that matters
        "iptm_mean": sum(good) / len(good),
        # Sanity gate, not a rank key: if the target protomers did not assemble
        # then the epitope is absent and the binder score is meaningless.
        "target_pair_iptm": pair_val(1, 2),
        "binder_ptm": chains_ptm.get("0"),
        "global_iptm": data.get("iptm"),          # recorded, never used to rank
        "complex_plddt": data.get("complex_plddt"),
        "complex_iplddt": data.get("complex_iplddt"),
        "ptm": data.get("ptm"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", type=Path, nargs="*", default=None,
                    help="one or more MPNN seqs dirs; omit to score only "
                         "--extra-seqs. Names are prefixed with the pool dir "
                         "because the plain and hisb pools reuse backbone "
                         "names and would otherwise overwrite each other.")
    ap.add_argument("--out", type=Path, default=ROOT / "boltz")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--per-backbone", type=int, default=None,
                    help="cap designs taken per backbone")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--recycling-steps", type=int, default=3)
    ap.add_argument("--diffusion-samples", type=int, default=1)
    ap.add_argument("--sampling-steps", type=int, default=200)
    ap.add_argument("--extra-seqs", type=Path, default=None,
                    help="FASTA of extra sequences (calibration controls)")
    ap.add_argument("--devices", type=int, default=1)
    ap.add_argument("--override", action="store_true",
                    help="recompute predictions that already exist. Off by "
                         "default so the same --tag can be re-run as new "
                         "backbones arrive and only new designs are folded.")
    ap.add_argument("--template", type=Path, default=ROOT / "target/tnf_AB.cif",
                    help="pin target chains to this geometry; '' to disable")
    args = ap.parse_args()

    binders: dict[str, str] = {}
    for d in (args.fasta_dir or []):
        pool = d.parent.name if d.name == "seqs" else d.name
        for name, seq in read_binders(d, args.per_backbone).items():
            binders[f"{pool}_{name}"] = seq

    if args.limit:
        binders = dict(list(binders.items())[: args.limit])

    if args.extra_seqs:
        name = None
        buf: list[str] = []
        for line in args.extra_seqs.read_text().splitlines():
            if line.startswith(">"):
                if name:
                    binders[name] = "".join(buf)
                name, buf = re.sub(r"\W+", "_", line[1:].strip()), []
            elif line.strip():
                buf.append(line.strip())
        if name:
            binders[name] = "".join(buf)

    if not binders:
        print("no sequences to score", file=sys.stderr)
        return 1

    tmpl = args.template if args.template and str(args.template) else None
    if tmpl and not tmpl.exists():
        print(f"template not found: {tmpl}", file=sys.stderr)
        return 1

    work = args.out / args.tag
    yaml_dir = work / "yaml"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    for name, seq in binders.items():
        write_yaml(yaml_dir / f"{name}.yaml", seq, TNF, tmpl)
    print(f"wrote {len(binders)} yaml -> {yaml_dir}")

    pred_dir = work / "pred"
    cmd = ["boltz", "predict", str(yaml_dir),
           "--out_dir", str(pred_dir),
           "--model", "boltz2",
           "--accelerator", "gpu",
           "--devices", str(args.devices),
           "--recycling_steps", str(args.recycling_steps),
           "--diffusion_samples", str(args.diffusion_samples),
           "--sampling_steps", str(args.sampling_steps),
           "--output_format", "pdb"]
    if args.override:
        cmd.append("--override")
    print("+", " ".join(cmd), flush=True)
    t0 = time.time()
    subprocess.run(cmd, check=True)
    dt = time.time() - t0

    rows = [r for r in (parse_confidence(pred_dir, n) for n in binders) if r]
    rows.sort(key=lambda r: r["iptm_min"], reverse=True)
    out_json = work / "scores.json"
    out_json.write_text(json.dumps(rows, indent=2))

    print(f"\n{len(rows)}/{len(binders)} scored in {dt/60:.1f} min "
          f"({dt/max(len(binders),1):.1f} s/design)")
    for r in rows[:15]:
        pa, pb = r["pair_iptm"]["A"], r["pair_iptm"]["B"]
        tp = r["target_pair_iptm"]
        print(f"  {r['name']:<28} min={r['iptm_min']:.3f} "
              f"(A={pa:.3f} B={pb:.3f})  targetAB={tp:.3f}  "
              f"binder_ptm={r['binder_ptm']:.3f}")
    print(f"scores -> {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
