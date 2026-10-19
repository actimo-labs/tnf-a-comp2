#!/usr/bin/env python
"""Convert the target PDB into an mmCIF Boltz-2 will accept as a template.

Two things have to be true or Boltz's template parser fails:

  1. mmCIF, not PDB. parse_pdb routes through the mmCIF path anyway and
     trips over a Biopython-written PDB.
  2. Each polymer entity must carry a non-empty `full_sequence`. gemmi's
     setup_entities() leaves it empty for a file with no SEQRES records,
     and boltz's parse_polymer then indexes `sequence[j]` off the end and
     dies with a bare `IndexError: list index out of range`.

Run under the boltz env:
  conda run -n boltz python scripts/make_template_cif.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=ROOT / "target/tnf_AB.pdb")
    ap.add_argument("--out", type=Path, default=ROOT / "target/tnf_AB.cif")
    args = ap.parse_args()

    st = gemmi.read_structure(str(args.pdb))
    st.setup_entities()
    st.assign_label_seq_id()

    # Boltz matches `template_id` against the mmCIF label_asym_id, and gemmi
    # invents ids like "Axp"/"Bxp" for those. Force them to the author chain
    # names so the YAML can say template_id: [A, B] and mean it.
    for model in st:
        for chain in model:
            for res in chain:
                res.subchain = chain.name
    for ent in st.entities:
        ent.subchains = [s for s in {c.name for m in st for c in m}
                         if s in set(ent.subchains) or s == ent.name]

    # Fill full_sequence from the observed residues of each entity's subchain.
    for ent in st.entities:
        if ent.entity_type != gemmi.EntityType.Polymer or ent.full_sequence:
            continue
        names: list[str] = []
        for model in st:
            for chain in model:
                for sub in chain.subchains():
                    if sub.subchain_id() in ent.subchains:
                        names.extend(r.name for r in sub)
                if names:
                    break
            if names:
                break
        ent.full_sequence = names

    bad = [e.name for e in st.entities
           if e.entity_type == gemmi.EntityType.Polymer and not e.full_sequence]
    if bad:
        raise SystemExit(f"entities still without a sequence: {bad}")

    st.make_mmcif_document().write_file(str(args.out))

    for ent in st.entities:
        print(f"entity {ent.name}: {len(ent.full_sequence)} residues")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
