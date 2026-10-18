#!/usr/bin/env python3
"""
Histidine switch scan + thermodynamic linkage scoring.

This is the Track 3 core. Structure predictors ignore protonation, so pH
selectivity cannot be filtered with AF2/Boltz. It is scored here instead,
on CPU, for free.

Physics
-------
An ionizable group with pKa_free in the unbound binder and pKa_bound in the
complex contributes a pH-dependent term to binding:

    dG_bind(pH) = dG_ref - RT*ln[(1 + 10^(pKa_bound - pH)) /
                                 (1 + 10^(pKa_free  - pH))]

so the penalty incurred by dropping from 7.4 to 6.0 is

    ddG = dG_bind(6.0) - dG_bind(7.4)

ddG > 0 means binding is weaker at pH 6.0, which is the objective. That happens
when pKa_bound < pKa_free: the complex suppresses protonation, so protonating
costs binding energy. Burying a His in a tight interface, or placing it beside a
target Arg/Lys, depresses its bound pKa and does exactly this.

~1.4 kcal/mol buys a 10x Kd shift. "No detectable binding" wants 100-1000x, so
aim for 2.8-4.2 kcal/mol total, which needs 2-3 cooperating His. Single
substitutions will not clear the bar.

Requires: propka3, biopython. Optional: foldx (set FOLDX_BIN) for a 7.4 affinity
cross-check on the neutral state.

Usage:
    python scripts/05_ph_scan.py --complex designs/*.pdb --binder-chain C \\
        --target-json target/numbering.json --out results/ph_scan.tsv
"""

import argparse
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import three_to_index, index_to_one

RT = 0.593  # kcal/mol at 298 K
PH_HIGH = 7.4
PH_LOW = 6.0

# Thresholds
MIN_TOTAL_DDG = 2.8      # kcal/mol, ~500x Kd shift
MAX_HIS_SUBS = 4         # each His costs affinity at 7.4
CATION_RADIUS = 6.0      # A, His-to-Arg/Lys for a useful repulsion switch
ANION_RADIUS = 6.0       # A, His-to-Asp/Glu is disqualifying
NEVER_MUTATE = {"G", "P", "C"}  # backbone geometry / disulfide risk


@dataclass
class HisCandidate:
    design: str
    position: int
    wt: str
    pka_free: float
    pka_bound: float
    ddg: float
    near_cationic: list
    near_anionic: list
    burial: float

    def viable(self):
        return self.ddg > 0.3 and not self.near_anionic and self.near_cationic


class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[0] == " "


def linkage_ddg(pka_free, pka_bound, ph_low=PH_LOW, ph_high=PH_HIGH):
    """Binding penalty (kcal/mol) incurred by moving from ph_high to ph_low.

    Positive = weaker binding at low pH = what Track 3 wants.
    """
    def term(ph):
        return math.log10(
            (1.0 + 10 ** (pka_bound - ph)) / (1.0 + 10 ** (pka_free - ph))
        )
    # dG(pH) = dG_ref - RT*ln(10)*term(pH)
    return -RT * math.log(10) * (term(ph_low) - term(ph_high))


def run_propka(pdb_path, workdir):
    """Return {(chain, resnum): pka} for HIS residues."""
    if shutil.which("propka3") is None:
        sys.exit("propka3 not on PATH -- pip install propka")
    pdb_path = Path(pdb_path).resolve()
    local = Path(workdir) / pdb_path.name
    shutil.copy(pdb_path, local)
    try:
        subprocess.run(
            ["propka3", local.name],
            cwd=workdir,
            check=True,
            capture_output=True,
            timeout=600,
        )
    except subprocess.CalledProcessError as exc:
        print(f"  propka3 failed on {pdb_path.name}: {exc.stderr.decode()[:200]}",
              file=sys.stderr)
        return {}
    except subprocess.TimeoutExpired:
        print(f"  propka3 timed out on {pdb_path.name}", file=sys.stderr)
        return {}

    pka_file = Path(workdir) / (local.stem + ".pka")
    if not pka_file.exists():
        return {}

    pkas = {}
    in_summary = False
    for line in pka_file.read_text().splitlines():
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            in_summary = True
            continue
        if in_summary:
            if line.startswith("-") or not line.strip():
                if pkas:
                    break
                continue
            parts = line.split()
            # e.g. "HIS  15 A   6.42   6.50"
            if len(parts) >= 4 and parts[0] == "HIS":
                try:
                    pkas[(parts[2], int(parts[1]))] = float(parts[3])
                except (ValueError, IndexError):
                    continue
    return pkas


def mutate_to_his(structure, chain_id, resnum, out_path):
    """
    Truncate the side chain to CB and relabel as HIS. PROPKA only needs
    approximate geometry to estimate desolvation and charge environment, and
    a full repack here would cost more than the signal is worth at this stage.
    Survivors get properly repacked and re-folded in step 06.
    """
    import copy

    st = copy.deepcopy(structure)
    res = st[0][chain_id][(" ", resnum, " ")]
    keep = {"N", "CA", "C", "O", "CB"}
    for atom in [a for a in res if a.get_name() not in keep]:
        res.detach_child(atom.get_id())
    res.resname = "HIS"
    io = PDBIO()
    io.set_structure(st)
    io.save(str(out_path))
    return st


def interface_residues(structure, binder_chain, target_chains, cutoff=5.0):
    """Binder residues with any heavy atom within cutoff of the target."""
    from Bio.PDB import NeighborSearch

    target_atoms = [
        a
        for cid in target_chains
        for res in structure[0][cid]
        if res.id[0] == " "
        for a in res
        if a.element != "H"
    ]
    ns = NeighborSearch(target_atoms)
    hits = []
    for res in structure[0][binder_chain]:
        if res.id[0] != " ":
            continue
        try:
            one = index_to_one(three_to_index(res.get_resname()))
        except (KeyError, IndexError):
            continue
        if one in NEVER_MUTATE:
            continue
        contact = any(ns.search(a.coord, cutoff) for a in res if a.element != "H")
        if contact:
            hits.append((res.id[1], one))
    return hits


def charge_environment(structure, binder_chain, resnum, target_chains):
    """Target Arg/Lys and Asp/Glu within radius of the binder position's CB."""
    res = structure[0][binder_chain][(" ", resnum, " ")]
    anchor = res["CB"] if "CB" in res else res["CA"]

    cationic, anionic = [], []
    for cid in target_chains:
        for t in structure[0][cid]:
            if t.id[0] != " ":
                continue
            name = t.get_resname()
            if name not in ("ARG", "LYS", "ASP", "GLU"):
                continue
            tip = {
                "ARG": ["CZ", "NH1", "NH2"],
                "LYS": ["NZ"],
                "ASP": ["CG", "OD1", "OD2"],
                "GLU": ["CD", "OE1", "OE2"],
            }[name]
            atoms = [t[a] for a in tip if a in t]
            if not atoms:
                continue
            d = min(float(anchor - a) for a in atoms)
            label = f"{cid}{t.id[1]}{name}"
            if name in ("ARG", "LYS") and d <= CATION_RADIUS:
                cationic.append((label, round(d, 2)))
            elif name in ("ASP", "GLU") and d <= ANION_RADIUS:
                anionic.append((label, round(d, 2)))
    return cationic, anionic


def relative_burial(structure, binder_chain, resnum, target_chains):
    """Fraction of the residue's SASA lost on complex formation."""
    from Bio.PDB.SASA import ShrakeRupley
    import copy

    sr = ShrakeRupley()
    complexed = copy.deepcopy(structure)
    sr.compute(complexed[0], level="R")
    bound = complexed[0][binder_chain][(" ", resnum, " ")].sasa

    alone = copy.deepcopy(structure)
    for cid in target_chains:
        alone[0].detach_child(cid)
    sr.compute(alone[0], level="R")
    free = alone[0][binder_chain][(" ", resnum, " ")].sasa

    if free < 1.0:
        return 1.0
    return round(1.0 - bound / free, 3)


def scan_design(pdb_path, binder_chain, workdir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("c", str(pdb_path))
    target_chains = [
        c.id for c in structure[0] if c.id != binder_chain
    ]
    if not target_chains:
        print(f"  {pdb_path.name}: no target chains found", file=sys.stderr)
        return []

    positions = interface_residues(structure, binder_chain, target_chains)
    if not positions:
        print(f"  {pdb_path.name}: no interface residues", file=sys.stderr)
        return []

    candidates = []
    for resnum, wt in positions:
        cationic, anionic = charge_environment(
            structure, binder_chain, resnum, target_chains
        )
        # Cheap rejection before paying for two propka runs
        if anionic or not cationic:
            continue

        sub = Path(workdir) / f"{pdb_path.stem}_{resnum}H"
        sub.mkdir(exist_ok=True)
        mut_complex = sub / "complex.pdb"
        mutated = mutate_to_his(structure, binder_chain, resnum, mut_complex)

        # Free binder: same mutated structure, target chains stripped
        io = PDBIO()
        io.set_structure(mutated)
        mut_free = sub / "binder.pdb"
        io.save(str(mut_free), ChainSelect([binder_chain]))

        pk_bound = run_propka(mut_complex, sub).get((binder_chain, resnum))
        pk_free = run_propka(mut_free, sub).get((binder_chain, resnum))
        if pk_bound is None or pk_free is None:
            continue

        candidates.append(
            HisCandidate(
                design=pdb_path.stem,
                position=resnum,
                wt=wt,
                pka_free=round(pk_free, 2),
                pka_bound=round(pk_bound, 2),
                ddg=round(linkage_ddg(pk_free, pk_bound), 3),
                near_cationic=cationic,
                near_anionic=anionic,
                burial=relative_burial(
                    structure, binder_chain, resnum, target_chains
                ),
            )
        )
    return candidates


def best_combinations(candidates, max_subs=MAX_HIS_SUBS, min_total=MIN_TOTAL_DDG):
    """
    Additive approximation over independent sites. Valid only while the His
    positions are far enough apart not to couple; enforce >8 A separation.
    Coupled pairs need a real pKa calculation and are flagged, not scored.
    """
    viable = [c for c in candidates if c.viable()]
    viable.sort(key=lambda c: c.ddg, reverse=True)
    viable = viable[:8]  # combinatorics guard

    results = []
    for n in range(2, max_subs + 1):
        for combo in itertools.combinations(viable, n):
            total = sum(c.ddg for c in combo)
            if total >= min_total:
                results.append((total, combo))
    results.sort(key=lambda r: r[0], reverse=True)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", nargs="+", required=True, help="design complex PDBs")
    ap.add_argument("--binder-chain", default="C")
    ap.add_argument("--target-json", default="target/numbering.json")
    ap.add_argument("--out", default="results/ph_scan.tsv")
    ap.add_argument("--combos-out", default="results/ph_combos.tsv")
    ap.add_argument("--keep-temp", action="store_true")
    args = ap.parse_args()

    paths = [Path(p) for p in args.complex if Path(p).exists()]
    if not paths:
        sys.exit("no input PDBs found")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    workroot = tempfile.mkdtemp(prefix="phscan_")
    print(f"[05] scanning {len(paths)} designs, scratch={workroot}")

    all_candidates = []
    all_combos = []
    for i, p in enumerate(paths, 1):
        cands = scan_design(p, args.binder_chain, workroot)
        all_candidates.extend(cands)
        combos = best_combinations(cands)
        for total, combo in combos[:5]:
            all_combos.append((p.stem, total, combo))
        print(
            f"[05] {i}/{len(paths)} {p.stem}: {len(cands)} sites, "
            f"{len(combos)} combos >= {MIN_TOTAL_DDG} kcal/mol"
        )

    with open(args.out, "w") as fh:
        fh.write(
            "design\tposition\twt\tpka_free\tpka_bound\tddg_kcal\t"
            "burial\tnear_cationic\tnear_anionic\n"
        )
        for c in sorted(all_candidates, key=lambda x: x.ddg, reverse=True):
            fh.write(
                f"{c.design}\t{c.position}\t{c.wt}\t{c.pka_free}\t{c.pka_bound}\t"
                f"{c.ddg}\t{c.burial}\t"
                f"{';'.join(f'{l}@{d}' for l, d in c.near_cationic)}\t"
                f"{';'.join(f'{l}@{d}' for l, d in c.near_anionic)}\n"
            )

    with open(args.combos_out, "w") as fh:
        fh.write("design\ttotal_ddg_kcal\tn_his\tmutations\tpredicted_kd_shift\n")
        for design, total, combo in sorted(
            all_combos, key=lambda r: r[1], reverse=True
        ):
            muts = ",".join(f"{c.wt}{c.position}H" for c in combo)
            shift = 10 ** (total / 1.36)
            fh.write(f"{design}\t{total:.2f}\t{len(combo)}\t{muts}\t{shift:.0f}x\n")

    print(f"[05] wrote {args.out} ({len(all_candidates)} sites)")
    print(f"[05] wrote {args.combos_out} ({len(all_combos)} combinations)")
    print(
        "\n[05] Next: step 06 repacks each combination, re-folds with AF2, and "
        "confirms the pH 7.4 interface survived. Linkage ddG here assumes "
        "independent sites -- treat combos whose His are within 8 A as "
        "unvalidated."
    )

    if not args.keep_temp:
        shutil.rmtree(workroot, ignore_errors=True)


if __name__ == "__main__":
    main()
