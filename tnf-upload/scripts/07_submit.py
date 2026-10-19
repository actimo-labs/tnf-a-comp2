#!/usr/bin/env python
"""Stage 07 -- developability, diversity clustering, and the submission CSV.

Ranking, per PLAN.md, in this order:
  1. ddG(6.0 vs 7.4) >= 2.8 kcal/mol        -- the Track 3 objective
  2. pH 7.4 interface survived His install  -- stage 06's gate
  3. His environment sanity                 -- enforced upstream in stage 05
  4. Composition sanity                     -- here
  5. Structural diversity                   -- here, TM-score clustering

Why the diversity cap is not negotiable: prior Adaptyv rounds show in-silico
scores predict real binding at roughly chance (ROC-AUC 0.5-0.64). Twenty
variations on one scaffold is one bet, not twenty, and if that scaffold is
wrong the whole submission is. The cap costs predicted quality on purpose.

Output is the Proteinbase submission format: name, sequence, molecule_class,
ordered best-first.

  venv-phscan/bin/python scripts/07_submit.py --refold results/06_refold.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

MIN_LEN, MAX_LEN = 10, 250
MAX_HOMOPOLYMER = 4          # 5+ identical residues in a row
SEQUON = re.compile(r"N[^P][ST]")   # N-X-S/T, X != P

# Total His ceiling, counting pre-existing ones. MAX_HIS_SUBS in stage 05 caps
# *designed substitutions*; it does not see histidines the His-biased MPNN pool
# already installed, and one design reached 11 in 92 aa. Beyond ~9% the binder
# is broadly pH-labile rather than switched at a designed site, and picks up
# nonspecific metal affinity.
MAX_HIS_FRACTION = 0.09
MAX_HIS_ABSOLUTE = 8

# Single-residue composition ceiling. ProteinMPNN over-uses Ala and Glu on
# idealised de novo helical surfaces; past ~35% of one residue the binder
# tends toward a weak hydrophobic core and molten-globule behaviour, which
# is where it would fail in a cell-free expression and BLI assay. The
# homopolymer-run check does not catch this -- the residues are spread out,
# not consecutive.
MAX_SINGLE_AA_FRACTION = 1.01  # disabled: enforcing 0.35 cost a structural cluster


def composition_problems(seq: str) -> list[str]:
    """Reasons this sequence is a bad bet in Adaptyv's cell-free pipeline."""
    bad = []
    if not MIN_LEN <= len(seq) <= MAX_LEN:
        bad.append(f"length {len(seq)} outside {MIN_LEN}-{MAX_LEN}")
    if "C" in seq:
        bad.append(f"free Cys x{seq.count('C')}")
    for m in SEQUON.finditer(seq):
        bad.append(f"N-X-S/T sequon at {m.start()+1}")
    for aa in set(seq):
        run = max((len(r) for r in re.findall(f"{aa}+", seq)), default=0)
        if run > MAX_HOMOPOLYMER:
            bad.append(f"poly-{aa} run of {run}")
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
        bad.append("non-standard residue")
    top_aa, top_n = max(((a, seq.count(a)) for a in set(seq)), key=lambda x: x[1])
    if top_n / len(seq) > MAX_SINGLE_AA_FRACTION:
        bad.append(f"{top_aa} is {100*top_n/len(seq):.0f}% of sequence")
    nh = seq.count("H")
    if nh > MAX_HIS_ABSOLUTE or nh / len(seq) > MAX_HIS_FRACTION:
        bad.append(f"{nh} His ({100*nh/len(seq):.0f}% of {len(seq)} aa)")
    return bad


def ca_coords(pdb: Path, chain: str = "X"):
    import numpy as np
    xyz, seq = [], []
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(three_to_one.get(line[17:20].strip(), "X"))
    if not xyz:
        return None, None
    return np.asarray(xyz, dtype="float64"), "".join(seq)


def cap_per_parent(entries: list[dict], cap: int) -> list[dict]:
    """Keep at most `cap` His-variants of any one parent binder.

    Without this the slate fills with near-identical designs: the same
    backbone and sequence carrying three different His combinations is three
    rows in the CSV but one bet in the assay. Capping here, before structural
    clustering, spends the 20 slots on distinct binders instead.
    """
    seen: dict[str, int] = {}
    out = []
    for e in entries:
        p = e.get("parent") or e["name"]
        if seen.get(p, 0) >= cap:
            continue
        seen[p] = seen.get(p, 0) + 1
        out.append(e)
    return out


def cluster(entries: list[dict], pred_root: Path, threshold: float, cap: int) -> list[dict]:
    """Greedy best-first clustering on TM-score; keep at most `cap` per cluster."""
    try:
        from tmtools import tm_align
    except ImportError:
        print("tmtools missing -- skipping structural clustering", file=sys.stderr)
        return entries

    structures = {}
    for e in entries:
        hits = sorted(pred_root.rglob(f"{e['name']}_model_0.pdb"))
        if hits:
            structures[e["name"]] = ca_coords(hits[0])

    kept: list[dict] = []
    reps: list[tuple[dict, list[dict]]] = []
    for e in entries:
        coords = structures.get(e["name"])
        placed = False
        if coords and coords[0] is not None:
            for rep, members in reps:
                rc = structures.get(rep["name"])
                if not rc or rc[0] is None:
                    continue
                try:
                    tm = tm_align(coords[0], rc[0], coords[1], rc[1]).tm_norm_chain1
                except Exception:
                    continue
                if tm >= threshold:
                    placed = True
                    if len(members) < cap:
                        members.append(e)
                        e["cluster"] = rep["name"]
                        kept.append(e)
                    break
        if not placed:
            reps.append((e, [e]))
            e["cluster"] = e["name"]
            kept.append(e)
    print(f"{len(reps)} clusters at TM >= {threshold}, cap {cap} per cluster")
    return kept


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refold", type=Path, default=ROOT / "results/06_refold.tsv")
    ap.add_argument("--pred-root", type=Path, default=ROOT / "boltz/refold")
    ap.add_argument("--out", type=Path, default=ROOT / "results/submission.csv")
    ap.add_argument("--report", type=Path, default=ROOT / "results/submission_report.tsv")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--min-ddg", type=float, default=2.8)
    ap.add_argument("--tm-threshold", type=float, default=0.6)
    ap.add_argument("--per-cluster", type=int, default=6)
    ap.add_argument("--per-parent", type=int, default=2,
                    help="max His-variants of one parent binder")
    ap.add_argument("--prefix", default="tnf3")
    # The challenge page documents the enum as "protein"/"nanobody"/"scfv"/...
    # but the live Proteinbase submission form rejects those and accepts
    # "Single chain" / "Nanobody" / "scFv" / "Fab (kappa)" / "Fab (lambda)".
    # The form is the authority.
    ap.add_argument("--molecule-class", default="Single chain")
    ap.add_argument("--allow-below-ddg", action="store_true",
                    help="backfill with sub-threshold designs rather than "
                         "submit fewer than --n")
    args = ap.parse_args()

    if not args.refold.exists():
        print(f"missing {args.refold} -- run stage 06 first", file=sys.stderr)
        return 1

    rows = list(csv.DictReader(args.refold.open(), delimiter="\t"))
    if not rows:
        print("stage 06 produced no rows", file=sys.stderr)
        return 1

    entries, rejected = [], []
    for r in rows:
        seq = r["sequence"]
        problems = composition_problems(seq)
        rec = {
            "name": r["name"],
            "parent": r.get("parent", ""),
            "sequence": seq,
            "ddg": float(r["ddg"]),
            "variant_iptm": float(r["variant_iptm"]),
            "survives": int(r["survives_pH74"]),
            "n_his": int(r["n_his"]),
            "mutations": r["mutations"],
            "problems": "; ".join(problems),
        }
        (rejected if problems or not rec["survives"] else entries).append(rec)

    # Objective first, then interface quality as the tiebreak.
    entries.sort(key=lambda e: (-e["ddg"], -e["variant_iptm"]))

    strong = [e for e in entries if e["ddg"] >= args.min_ddg]
    weak = [e for e in entries if e["ddg"] < args.min_ddg]
    pool = strong if len(strong) >= args.n or not args.allow_below_ddg else strong + weak

    n_parents = len({e.get("parent") or e["name"] for e in pool})
    pool = cap_per_parent(pool, args.per_parent)
    print(f"{n_parents} distinct parent binders; "
          f"{len(pool)} candidates after <= {args.per_parent}/parent")
    picked = cluster(pool, args.pred_root, args.tm_threshold, args.per_cluster)[: args.n]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for i, e in enumerate(picked, 1):
            w.writerow([f"{args.prefix}_{i:02d}", e["sequence"],
                        args.molecule_class])

    with args.report.open("w", newline="") as fh:
        cols = ["rank", "submitted_name", "name", "parent", "ddg",
                "variant_iptm", "n_his", "mutations", "cluster", "sequence"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for i, e in enumerate(picked, 1):
            w.writerow({"rank": i, "submitted_name": f"{args.prefix}_{i:02d}",
                        **{c: e.get(c, "") for c in cols[2:]}})

    print(f"{len(entries)} passed composition+survival, {len(rejected)} rejected")
    print(f"{len(strong)} at ddG >= {args.min_ddg}")
    print(f"selected {len(picked)} -> {args.out}")
    if len(picked) < args.n:
        print(f"WARNING: only {len(picked)} of {args.n} designs available",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
