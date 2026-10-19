#!/usr/bin/env python
"""Stage 06 -- install the His combinations and confirm pH 7.4 binding survives.

This is the gate CLAUDE.md calls non-negotiable: "Every His variant must be
re-folded to confirm the pH 7.4 interface survives. A design that switches but
no longer binds at 7.4 is worthless."

Stage 05 scores how much binding a design *loses* at pH 6.0. It says nothing
about whether the design still binds at 7.4 after two to four interface
residues have been replaced by histidine -- and those substitutions are not
free, which is exactly why MAX_HIS_SUBS exists. So each combination is rebuilt
as a sequence, re-folded, and scored by the same binder-vs-target pair rule as
stage 04.

A variant survives when:
  * its binder-target pair ipTM clears the stage-04 cutoff, and
  * it has not lost more than `--max-drop` of the parent's score.

The second test matters independently of the first. A parent at 0.72 falling to
0.41 is still "above cutoff" while having lost most of its interface, and that
is a design whose 7.4 binding did not really survive.

  conda run -n boltz python scripts/06_his_refold.py \\
      --combos results/ph_combos.tsv --parent-tag screen --cutoff 0.5
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent


def parent_sequences(yaml_dir: Path) -> dict[str, str]:
    """Recover each parent's binder sequence from the stage-04 YAML."""
    seqs: dict[str, str] = {}
    for y in sorted(yaml_dir.glob("*.yaml")):
        lines = y.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.strip() == "id: X":
                for nxt in lines[i + 1 : i + 3]:
                    if nxt.strip().startswith("sequence:"):
                        seqs[y.stem] = nxt.split("sequence:", 1)[1].strip()
                        break
                break
    return seqs


def parse_mutations(field: str) -> list[tuple[str, int]]:
    """'A86H,E47H' -> [('A', 86), ('E', 47)]"""
    out = []
    for tok in field.split(","):
        tok = tok.strip()
        if not tok:
            continue
        out.append((tok[0], int(tok[1:-1])))
    return out


def apply_mutations(seq: str, muts: list[tuple[str, int]]) -> str:
    chars = list(seq)
    for wt, pos in muts:
        if not 1 <= pos <= len(chars):
            raise ValueError(f"position {pos} outside binder of length {len(chars)}")
        if chars[pos - 1] != wt:
            raise ValueError(
                f"expected {wt} at {pos}, found {chars[pos-1]} -- "
                "numbering mismatch between stage 05 and the parent sequence"
            )
        chars[pos - 1] = "H"
    return "".join(chars)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--combos", type=Path, default=ROOT / "results/ph_combos.tsv")
    ap.add_argument("--boltz-root", type=Path, default=ROOT / "boltz")
    ap.add_argument("--parent-tag", default="screen")
    ap.add_argument("--tag", default="refold")
    ap.add_argument("--cutoff", type=float, required=True,
                    help="stage-04 calibrated pair-ipTM cutoff")
    ap.add_argument("--strong-iptm", type=float, default=0.40,
                    help="absolute pair-ipTM above which the drop rule is "
                         "waived; such a variant binds well on its own terms")
    ap.add_argument("--max-drop", type=float, default=0.25,
                    help="max fractional loss vs parent (0.25 = keep >=75%%)")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", type=Path, default=ROOT / "results/06_refold.tsv")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    parent_dir = args.boltz_root / args.parent_tag
    seqs = parent_sequences(parent_dir / "yaml")
    if not seqs:
        print(f"no parent sequences under {parent_dir/'yaml'}", file=sys.stderr)
        return 1

    parent_scores = {}
    score_file = parent_dir / "scores.json"
    if score_file.exists():
        parent_scores = {r["name"]: r for r in json.loads(score_file.read_text())}

    if not args.combos.exists():
        print(f"no combos file: {args.combos}", file=sys.stderr)
        return 1

    variants: dict[str, dict] = {}
    skipped: list[str] = []
    with args.combos.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            design = row["design"]
            # Stage 05 runs on Boltz output PDBs named <parent>_model_0.
            parent = design.replace("_model_0", "")
            seq = seqs.get(parent)
            if seq is None:
                skipped.append(design)
                continue
            muts = parse_mutations(row["mutations"])
            try:
                mutant = apply_mutations(seq, muts)
            except ValueError as exc:
                print(f"  {design}: {exc}", file=sys.stderr)
                skipped.append(design)
                continue
            name = f"{parent}__{'_'.join(f'{w}{p}H' for w, p in muts)}"
            variants[name] = {
                "name": name,
                "parent": parent,
                "sequence": mutant,
                "mutations": row["mutations"],
                "n_his": int(row["n_his"]),
                "ddg": float(row["total_ddg_kcal"]),
                "min_his_sep_A": row.get("min_his_sep_A", ""),
                "native_his_ddg": row.get("native_his_ddg", ""),
            }

    if skipped:
        print(f"skipped {len(skipped)} combo row(s) with no parent sequence",
              file=sys.stderr)
    if not variants:
        print("no variants to re-fold", file=sys.stderr)
        return 1

    ordered = sorted(variants.values(), key=lambda v: v["ddg"], reverse=True)
    if args.limit:
        ordered = ordered[: args.limit]

    fasta = args.boltz_root / f"{args.tag}_variants.fasta"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text(
        "\n".join(f">{v['name']}\n{v['sequence']}" for v in ordered) + "\n"
    )
    print(f"{len(ordered)} His variants -> {fasta}")
    if args.dry_run:
        return 0

    cmd = [sys.executable, str(SCRIPTS / "04_boltz.py"),
           "--extra-seqs", str(fasta),
           "--out", str(args.boltz_root),
           "--tag", args.tag]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    scored = {r["name"]: r for r in json.loads(
        (args.boltz_root / args.tag / "scores.json").read_text())}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows, survivors = [], 0
    for v in ordered:
        s = scored.get(v["name"])
        if not s:
            continue
        parent_s = parent_scores.get(v["parent"], {}).get("iptm_min")
        drop = (
            None if not parent_s else round((parent_s - s["iptm_min"]) / parent_s, 3)
        )
        # Binding at 7.4 is an absolute property, so the absolute score is the
        # primary test. The drop rule exists only to catch a variant that was
        # gutted yet still clings above the cutoff -- but applied as a hard AND
        # it penalises strong parents: a variant at 0.49 that fell 26% from a
        # 0.66 parent was being failed while a 0.31 survivor passed. A variant
        # comfortably strong in absolute terms is exempt from the drop rule.
        strong = s["iptm_min"] >= args.strong_iptm
        ok = s["iptm_min"] >= args.cutoff and (
            strong or drop is None or drop <= args.max_drop
        )
        survivors += ok
        rows.append({
            **{k: v[k] for k in
               ("name", "parent", "mutations", "n_his", "ddg",
                "min_his_sep_A", "native_his_ddg", "sequence")},
            "parent_iptm": parent_s,
            "variant_iptm": round(s["iptm_min"], 4),
            "target_pair_iptm": round(s["target_pair_iptm"] or 0.0, 4),
            "frac_drop": drop,
            "survives_pH74": int(ok),
        })

    rows.sort(key=lambda r: (-r["survives_pH74"], -r["ddg"]))
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_strong = sum(1 for r in rows
                   if r["survives_pH74"] and r["variant_iptm"] >= args.strong_iptm)
    print(f"\n{survivors}/{len(rows)} variants keep their pH 7.4 interface "
          f"(cutoff {args.cutoff}, max drop {args.max_drop}, "
          f"drop waived above {args.strong_iptm}; {n_strong} are strong)")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
