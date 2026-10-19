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

# Two His closer than this couple electrostatically: protonating one raises the
# pKa penalty on the other, so summing their individually-computed ddG
# overestimates the pair. The additive model is only valid beyond it.
MIN_HIS_SEP = 8.0        # A, CB-CB

# Native target His titrate too. TNF has three (15, 73, 78) and B15 sits 6.6 A
# from the hotspot patch, close enough to contribute to the linkage -- and
# close enough that a binder carboxylate placed beside it would REVERSE the
# switch, strengthening binding at pH 6.0.
TARGET_HIS_RADIUS = 10.0  # A, how close a target His must be to count


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
    cb: tuple = (0.0, 0.0, 0.0)   # CB coords, for the MIN_HIS_SEP check

    def viable(self):
        return self.ddg > 0.3 and not self.near_anionic and self.near_cationic

    def distance_to(self, other) -> float:
        return math.dist(self.cb, other.cb)


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
    # Callers sometimes already wrote the file into workdir.
    if local.resolve() != pdb_path:
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


HIS_SIDECHAIN = ("CG", "ND1", "CD2", "CE1", "NE2")


def _his_template(structure):
    """A complete native HIS to graft from. TNF has three per protomer."""
    for chain in structure[0]:
        for res in chain:
            if res.get_resname() == "HIS" and all(
                a in res for a in ("N", "CA", "CB", *HIS_SIDECHAIN)
            ):
                return res
    return None


def _rotate_about_axis(coords, origin, axis, angle):
    """Rodrigues rotation of coords about the line origin->origin+axis."""
    import numpy as np

    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    c, s = math.cos(angle), math.sin(angle)
    v = np.asarray(coords, dtype=float) - np.asarray(origin, dtype=float)
    rot = (v * c
           + np.cross(axis, v) * s
           + np.outer(np.dot(v, axis), axis) * (1.0 - c))
    return rot + np.asarray(origin, dtype=float)


def mutate_to_his(structure, chain_id, resnum, out_path, template=None):
    """
    Replace a residue with a real histidine, imidazole ring included.

    This cannot be a CB stub. PROPKA identifies titratable histidines by their
    ring atoms (ND1/NE2); given only N/CA/C/O/CB it does not register the
    residue as a His at all and silently omits it from the SUMMARY block. The
    caller then reads back `None` for the pKa and drops the site, which made
    the whole scan return zero candidates for every design.

    The ring is grafted from a complete native His by superimposing on the
    N/CA/CB frame, then spun about the CA-CB axis to the least-clashing chi1.
    Approximate, but it is a genuine imidazole in a plausible rotamer, which
    is what PROPKA's desolvation and charge terms need.
    """
    import copy

    import numpy as np
    from Bio.PDB import Superimposer

    st = copy.deepcopy(structure)
    res = st[0][chain_id][(" ", resnum, " ")]

    if template is None:
        template = _his_template(structure)
    if template is None:
        raise ValueError("no complete HIS in the structure to graft from")
    if "CB" not in res:
        raise ValueError(f"residue {chain_id}{resnum} has no CB to build from")

    # Superimpose the template's backbone+CB frame onto this residue's.
    frame = ("N", "CA", "CB")
    sup = Superimposer()
    sup.set_atoms([res[a] for a in frame], [template[a] for a in frame])

    moved = copy.deepcopy(template)
    sup.apply(list(moved.get_atoms()))

    # Drop everything past CB, then graft the ring.
    for atom in [a for a in res if a.get_name() not in ("N", "CA", "C", "O", "CB")]:
        res.detach_child(atom.get_id())

    ring = [copy.deepcopy(moved[a]) for a in HIS_SIDECHAIN]

    # Pick the chi1 that clashes least with everything else nearby.
    others = [a for a in st[0].get_atoms()
              if a.get_parent() is not res and a.element != "H"]
    kd_origin = np.asarray(res["CB"].coord, dtype=float)
    axis = kd_origin - np.asarray(res["CA"].coord, dtype=float)
    near = [a for a in others
            if np.linalg.norm(np.asarray(a.coord) - kd_origin) < 12.0]

    best, best_score = None, None
    for step in range(12):
        angle = step * math.pi / 6.0
        coords = _rotate_about_axis([a.coord for a in ring], kd_origin, axis, angle)
        score = 0.0
        for c in coords:
            for a in near:
                d = float(np.linalg.norm(np.asarray(a.coord) - c))
                if d < 3.2:
                    score += (3.2 - d) ** 2
        if best_score is None or score < best_score:
            best, best_score = coords, score

    for atom, coord in zip(ring, best):
        atom.set_coord(coord.astype("f"))
        res.add(atom)
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
        # A position that is already His is not a substitution. Scoring it as
        # one let "H52H" appear in a mutation list, consuming a slot of the
        # MAX_HIS_SUBS budget while changing nothing, and double-counted it
        # against native_his_ddg which now covers pre-existing binder His.
        if wt == "H":
            continue
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
                cb=tuple(float(x) for x in structure[0][binder_chain][
                    (" ", resnum, " ")]["CB"].coord),
            )
        )
    return candidates, native_his_ddg(structure, binder_chain, target_chains,
                                      workdir, pdb_path)


def native_his_ddg(structure, binder_chain, target_chains, workdir, pdb_path):
    """Linkage contribution of histidines that are already present.

    Two sources, both of which titrate between 7.4 and 6.0 regardless of what
    we substitute:

    * The target's own. TNF carries three His per protomer and B15 sits 6.6 A
      from the hotspot patch.
    * The binder's own. The His-biased ProteinMPNN pool installs ~4 His per
      binder before this scan runs, and one design reached 8. Leaving those
      out made the reported ddG for every `hisb` design incomplete.

    Only His near the interface are counted; a His 25 A away does not change
    its protonation because a binder arrived. The contribution can be
    negative, eating into the designed switch rather than adding to it.
    """
    import numpy as np

    binder_atoms = [a for a in structure[0][binder_chain].get_atoms()
                    if a.element != "H"]
    if not binder_atoms:
        return 0.0, []
    binder_xyz = np.array([a.coord for a in binder_atoms], dtype=float)
    target_xyz = np.array(
        [a.coord for cid in target_chains
         for a in structure[0][cid].get_atoms() if a.element != "H"],
        dtype=float,
    )

    near = []
    for cid in target_chains:
        for res in structure[0][cid]:
            if res.get_resname() != "HIS" or "CB" not in res:
                continue
            d = float(np.min(np.linalg.norm(
                binder_xyz - np.asarray(res["CB"].coord, dtype=float), axis=1)))
            if d <= TARGET_HIS_RADIUS:
                near.append((cid, res.id[1], round(d, 2)))

    # Pre-existing binder His, measured against the target rather than itself.
    if len(target_xyz):
        for res in structure[0][binder_chain]:
            if res.get_resname() != "HIS" or "CB" not in res:
                continue
            d = float(np.min(np.linalg.norm(
                target_xyz - np.asarray(res["CB"].coord, dtype=float), axis=1)))
            if d <= TARGET_HIS_RADIUS:
                near.append((binder_chain, res.id[1], round(d, 2)))

    if not near:
        return 0.0, []

    sub = Path(workdir) / f"{Path(pdb_path).stem}_native"
    sub.mkdir(exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    complex_pdb = sub / "complex.pdb"
    io.save(str(complex_pdb))
    target_pdb = sub / "target.pdb"
    io.save(str(target_pdb), ChainSelect(target_chains))
    binder_pdb = sub / "binder.pdb"
    io.save(str(binder_pdb), ChainSelect([binder_chain]))

    bound = run_propka(complex_pdb, sub)
    # The unbound reference differs by which side the His sits on: a target
    # His is unperturbed in the free target, a binder His in the free binder.
    free_target = run_propka(target_pdb, sub)
    free_binder = run_propka(binder_pdb, sub)

    total, detail = 0.0, []
    for cid, resnum, dist in near:
        free = free_binder if cid == binder_chain else free_target
        pb, pf = bound.get((cid, resnum)), free.get((cid, resnum))
        if pb is None or pf is None:
            continue
        ddg = linkage_ddg(pf, pb)
        total += ddg
        detail.append((f"{cid}{resnum}HIS", round(dist, 2), round(ddg, 3)))
    return round(total, 3), detail


def best_combinations(candidates, max_subs=MAX_HIS_SUBS, min_total=MIN_TOTAL_DDG,
                      min_sep=MIN_HIS_SEP, base_ddg=0.0):
    """
    Additive approximation over independent sites.

    The additivity only holds while the sites are electrostatically
    independent. Two His within ~8 A couple: protonating the first makes
    protonating the second more costly, so the individually-computed ddG
    values do not sum -- the additive total overestimates the pair, and
    overestimating is the dangerous direction here because it promotes
    designs that will not actually switch.

    Combinations containing any pair closer than `min_sep` are therefore
    rejected outright rather than scored. They are not necessarily bad
    designs, but scoring them would need a real coupled pKa calculation
    (constant-pH MD or a Poisson-Boltzmann titration) and that does not fit
    in the time available.

    `base_ddg` carries the contribution of native target His, which titrate
    regardless of what the binder does and so apply to every combination.
    """
    viable = [c for c in candidates if c.viable()]
    viable.sort(key=lambda c: c.ddg, reverse=True)
    viable = viable[:8]  # combinatorics guard

    results, rejected = [], 0
    for n in range(2, max_subs + 1):
        for combo in itertools.combinations(viable, n):
            closest = min(
                (a.distance_to(b) for a, b in itertools.combinations(combo, 2)),
                default=float("inf"),
            )
            if closest < min_sep:
                rejected += 1
                continue
            total = base_ddg + sum(c.ddg for c in combo)
            if total >= min_total:
                results.append((total, combo, round(closest, 2)))

    if rejected:
        print(f"  rejected {rejected} combination(s) with His closer than "
              f"{min_sep} A (additive ddG invalid)", file=sys.stderr)
    results.sort(key=lambda r: r[0], reverse=True)
    return results


def _scan_worker(job):
    """Module-level so ProcessPoolExecutor can pickle it.

    Each design gets its own scratch subdir, so designs are independent. One
    failing design must not take the whole batch down with it.
    """
    path, binder_chain, workroot = job
    try:
        return path, scan_design(path, binder_chain, workroot), None
    except Exception as exc:
        return path, ([], (0.0, [])), exc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", nargs="+", required=True, help="design complex PDBs")
    ap.add_argument("--binder-chain", default="C")
    ap.add_argument("--target-json", default="target/numbering.json")
    ap.add_argument("--out", default="results/ph_scan.tsv")
    ap.add_argument("--combos-out", default="results/ph_combos.tsv")
    ap.add_argument("--keep-temp", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                    help="designs scanned concurrently")
    args = ap.parse_args()

    paths = [Path(p) for p in args.complex if Path(p).exists()]
    if not paths:
        sys.exit("no input PDBs found")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    workroot = tempfile.mkdtemp(prefix="phscan_")
    print(f"[05] scanning {len(paths)} designs, scratch={workroot}")

    all_candidates = []
    all_combos = []

    if args.jobs > 1 and len(paths) > 1:
        # Processes, not threads. scan_design spends most of its time inside
        # Biopython (deepcopy, ShrakeRupley SASA) holding the GIL, so a thread
        # pool pinned the whole scan to ~1 core -- 107% CPU with --jobs 10.
        # Fork keeps module globals, so the worker needs no setup.
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            scanned = list(pool.map(
                _scan_worker,
                [(p, args.binder_chain, workroot) for p in paths],
                chunksize=1,
            ))
    else:
        scanned = [_scan_worker((p, args.binder_chain, workroot))
                   for p in paths]

    for i, (p, (cands, (base_ddg, base_detail)), err) in enumerate(scanned, 1):
        if err is not None:
            print(f"[05] {i}/{len(paths)} {p.stem}: FAILED ({err})",
                  file=sys.stderr)
            continue
        all_candidates.extend(cands)
        combos = best_combinations(cands, base_ddg=base_ddg)
        for total, combo, sep in combos[:5]:
            all_combos.append((p.stem, total, combo, sep, base_ddg))
        native = ",".join(f"{n}:{d:+.2f}" for n, _, d in base_detail) or "none"
        print(
            f"[05] {i}/{len(paths)} {p.stem}: {len(cands)} sites, "
            f"{len(combos)} combos >= {MIN_TOTAL_DDG} kcal/mol "
            f"| native His {base_ddg:+.2f} ({native})"
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
        fh.write("design\ttotal_ddg_kcal\tn_his\tmutations\t"
                 "min_his_sep_A\tnative_his_ddg\tpredicted_kd_shift\n")
        for design, total, combo, sep, base in sorted(
            all_combos, key=lambda r: r[1], reverse=True
        ):
            muts = ",".join(f"{c.wt}{c.position}H" for c in combo)
            shift = 10 ** (total / 1.36)
            fh.write(f"{design}\t{total:.2f}\t{len(combo)}\t{muts}\t"
                     f"{sep}\t{base:+.2f}\t{shift:.0f}x\n")

    print(f"[05] wrote {args.out} ({len(all_candidates)} sites)")
    print(f"[05] wrote {args.combos_out} ({len(all_combos)} combinations)")
    print(
        "\n[05] Next: step 06 repacks each combination, re-folds with Boltz-2, and "
        "confirms the pH 7.4 interface survived. Linkage ddG here assumes "
        "independent sites -- combos with His closer than 8 A were rejected outright."
    )

    if not args.keep_temp:
        shutil.rmtree(workroot, ignore_errors=True)


if __name__ == "__main__":
    main()
