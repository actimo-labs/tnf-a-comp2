#!/usr/bin/env python3
"""
Prepare the human TNF-alpha two-protomer target crop for binder design.

Produces:
  target/tnf_AB.pdb        two adjacent protomers, renumbered contiguously
  target/numbering.json    mature <-> crystal <-> crop numbering map
  target/hotspots.txt      hotspot residues in CROP numbering, for RFdiffusion

Why a two-protomer crop: the receptor-binding site sits in the cleft between
adjacent subunits, so both contribute. Keeping all three chains lets the binder
straddle a symmetry mate and wastes GPU time.

Usage:
    python scripts/01_prep_target.py --pdb 9OJO --chains A B
    python scripts/01_prep_target.py --pdb 1TNF --chains A B   # fallback
"""

import argparse
import json
import sys
from pathlib import Path

from Bio.PDB import PDBList, PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

# Mature human TNF-alpha soluble domain, UniProt P01375 residues 77-233.
# Mature numbering: position 1 = V.
MATURE_SEQ = (
    "VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLL"
    "THTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYF"
    "GIIAL"
)

# Receptor-binding groove, mature numbering. Chain roles:
#   "bottom" = E-F loop side, contributed by one protomer
#   "apex"   = 29-35 loop side, contributed by the adjacent protomer
EPITOPE = {
    "bottom": [86, 87, 143, 144, 145, 146, 147, 148],
    "apex": [29, 30, 31, 32, 33, 34, 35],
}

# Track 3 switch chemistry. Designed His near these is good (repulsion when
# protonated at pH 6.0); near the anionic set is actively harmful.
CATIONIC = [31, 32, 65, 90, 98, 112, 131, 138]
ANIONIC = [23, 110, 116, 127, 130, 135, 140, 143, 146]


class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        # standard amino acids only; drop waters, ligands, ions
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def fetch(pdb_id, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(outdir), file_format="pdb")
    if not Path(path).exists():
        sys.exit(f"could not fetch {pdb_id}")
    return Path(path)


def build_numbering_map(structure, chains):
    """
    Crystal numbering in TNF structures is already mature numbering, but
    residues 1-5 are usually disordered and chains may have gaps. Map every
    observed residue to its mature index by sequence, not by assuming the
    author numbering is continuous.
    """
    from Bio.PDB.Polypeptide import three_to_index, index_to_one

    mapping = {}
    for chain_id in chains:
        chain = structure[0][chain_id]
        entries = []
        for res in chain:
            if res.id[0] != " ":
                continue
            try:
                one = index_to_one(three_to_index(res.get_resname()))
            except (KeyError, IndexError):
                continue
            author_num = res.id[1]
            # crystal numbering in 1TNF/9OJO == mature numbering
            mature = author_num
            ok = (
                1 <= mature <= len(MATURE_SEQ)
                and MATURE_SEQ[mature - 1] == one
            )
            entries.append(
                {
                    "author": author_num,
                    "mature": mature,
                    "resname": one,
                    "seq_match": ok,
                }
            )
        mismatches = [e for e in entries if not e["seq_match"]]
        if mismatches:
            print(
                f"  WARNING chain {chain_id}: {len(mismatches)} residues do not "
                f"match the mature sequence at their author number. "
                f"First: {mismatches[:3]}",
                file=sys.stderr,
            )
        mapping[chain_id] = entries
    return mapping


def compute_sasa(structure, chains):
    sr = ShrakeRupley()
    sr.compute(structure[0], level="R")
    sasa = {}
    for chain_id in chains:
        for res in structure[0][chain_id]:
            if res.id[0] == " ":
                sasa[(chain_id, res.id[1])] = round(res.sasa, 1)
    return sasa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="9OJO", help="9OJO (1.36 A) preferred; 1TNF fallback")
    ap.add_argument("--chains", nargs=2, default=["A", "B"])
    ap.add_argument("--outdir", default="target")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    raw = fetch(args.pdb, outdir / "raw")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(args.pdb, str(raw))

    present = {c.id for c in structure[0]}
    missing = [c for c in args.chains if c not in present]
    if missing:
        sys.exit(f"chains {missing} not in {args.pdb} (present: {sorted(present)})")

    print(f"[01] {args.pdb}, chains {args.chains}")
    numbering = build_numbering_map(structure, args.chains)
    sasa = compute_sasa(structure, args.chains)

    # Write the cropped target
    io = PDBIO()
    io.set_structure(structure)
    out_pdb = outdir / "tnf_AB.pdb"
    io.save(str(out_pdb), ChainSelect(args.chains))
    print(f"[01] wrote {out_pdb}")

    # Resolve hotspots and sanity-check each one is actually observed and exposed.
    bottom_chain, apex_chain = args.chains
    hotspots = []
    report = []
    for role, chain_id in (("bottom", bottom_chain), ("apex", apex_chain)):
        observed = {e["mature"]: e for e in numbering[chain_id]}
        for mature in EPITOPE[role]:
            entry = observed.get(mature)
            if entry is None:
                report.append(f"  MISSING {chain_id}{mature} ({role}) - disordered in {args.pdb}")
                continue
            area = sasa.get((chain_id, entry["author"]), 0.0)
            flag = ""
            if mature in CATIONIC:
                flag = "  <- cationic, good His neighbour"
            elif mature in ANIONIC:
                flag = "  <- anionic, keep His away"
            if area < 10:
                report.append(
                    f"  BURIED  {chain_id}{mature}{entry['resname']} SASA={area} - skipping{flag}"
                )
                continue
            hotspots.append(f"{chain_id}{entry['author']}")
            report.append(
                f"  ok      {chain_id}{mature}{entry['resname']} SASA={area}{flag}"
            )

    print("[01] hotspot check:")
    print("\n".join(report))

    (outdir / "hotspots.txt").write_text(",".join(hotspots) + "\n")
    print(f"[01] hotspots ({len(hotspots)}): {','.join(hotspots)}")

    with open(outdir / "numbering.json", "w") as fh:
        json.dump(
            {
                "pdb": args.pdb,
                "chains": args.chains,
                "bottom_chain": bottom_chain,
                "apex_chain": apex_chain,
                "numbering": numbering,
                "epitope_mature": EPITOPE,
                "cationic_mature": CATIONIC,
                "anionic_mature": ANIONIC,
                "hotspots_author": hotspots,
            },
            fh,
            indent=2,
        )
    print(f"[01] wrote {outdir / 'numbering.json'}")

    if any("MISSING" in r for r in report):
        print(
            "\n[01] Some apex-loop residues are disordered. 1TNF's 29-36 loop is "
            "poorly ordered and that loop matters for this track. Re-run with "
            "--pdb 9OJO, or graft the ordered loop from 4TSV.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
