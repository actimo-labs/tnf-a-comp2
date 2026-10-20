#!/usr/bin/env python
"""Render a designed binder:TNF-alpha complex with PyMOL.

Shows what the design is actually doing: the binder sitting in the
receptor-binding groove at the protomer-protomer interface, with the switch
histidines picked out against the target arginine pair they push against when
protonated at pH 6.0.

Run headless:
  conda run -n viz pymol -cq scripts/render_complex.py -- \
      results/structures/complex/tnf3_01_complex.pdb figures/tnf3_01 "Q78H,K16H,E31H"
"""
import sys
import os

from pymol import cmd


def autocrop(path, pad=40):
    """Trim the white border PyMOL leaves around the subject.

    PyMOL frames on the bounding box of the whole loaded object, so a wide,
    flat complex ends up with a third of the canvas empty no matter how the
    zoom is set. Cropping afterwards is simpler than fighting the camera.
    """
    try:
        from PIL import Image, ImageChops
    except ImportError:
        return
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if not bbox:
        return
    l, t, r, b = bbox
    l, t = max(0, l - pad), max(0, t - pad)
    r, b = min(im.width, r + pad), min(im.height, b + pad)
    im.crop((l, t, r, b)).save(path)


def render(pdb_path, out_stem, mutations=""):
    os.makedirs(os.path.dirname(out_stem) or ".", exist_ok=True)

    cmd.reinitialize()
    cmd.load(pdb_path, "cx")
    cmd.remove("solvent or hydrogens")

    # Chain X is the binder; A and B are the two TNF protomers.
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.set("cartoon_transparency", 0.0)

    # Target muted and cool, binder warm and saturated, so the eye goes to the
    # design first and the two protomers stay distinguishable from each other.
    cmd.color("grey70", "chain A")
    cmd.color("slate", "chain B")
    cmd.color("orange", "chain X")
    cmd.set("cartoon_transparency", 0.25, "chain A or chain B")

    # Switch histidines on the binder.
    his_resi = [m.strip()[1:-1] for m in mutations.split(",") if m.strip()]
    if his_resi:
        sel = "chain X and resi " + "+".join(his_resi)
        cmd.show("sticks", f"({sel}) and not (name C+N+O)")
        cmd.color("yellow", f"({sel}) and elem C")
        cmd.set("stick_radius", 0.22, sel)

    # The target cations each switch His actually sits beside. These are found
    # from the structure rather than assumed: the target is a symmetric
    # two-protomer crop, so both chains carry residues 31/32, and a design can
    # engage either copy. Hardcoding the chain B pair pointed the figure at
    # arginines 28-42 A away from the binder.
    if his_resi:
        partners = cmd.get_model(
            f"(chain A or chain B) and resn ARG+LYS and sidechain "
            f"and (byres ((chain X and resi {'+'.join(his_resi)}) around 7))"
        )
        near = sorted({(a.chain, a.resi) for a in partners.atom})
        if near:
            sel = " or ".join(f"(chain {c} and resi {r})" for c, r in near)
            cmd.show("sticks", f"({sel}) and not (name C+N+O)")
            cmd.color("marine", f"({sel}) and elem C")
            cmd.set("stick_radius", 0.22, f"({sel})")
            print("cation partners: " + ", ".join(f"{c}{r}" for c, r in near))

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.18)
    cmd.set("direct", 0.6)

    # View 1: the interface. Orient on the binder plus whatever target residues
    # it actually touches, so the framing follows the real epitope rather than
    # a guessed one, and centre on that selection to avoid dead space.
    cmd.select("epitope", "byres ((chain A or chain B) within 5 of chain X)")
    cmd.orient("chain X or epitope")
    cmd.center("chain X or epitope")
    cmd.zoom("chain X or epitope", 3)
    # Wide frame: a binder lying along a groove is a broad, flat subject, and
    # a 4:3 canvas leaves a third of the image empty.
    cmd.ray(1800, 1000)
    cmd.png(f"{out_stem}_interface.png", dpi=300)

    # View 2: whole complex, so the binder is seen in context of the trimer crop.
    cmd.orient("cx")
    cmd.zoom("cx", 2)
    cmd.ray(1600, 1200)
    cmd.png(f"{out_stem}_overview.png", dpi=300)

    # View 3: same as 2, rotated 90 degrees about y.
    cmd.turn("y", 90)
    cmd.ray(1600, 1200)
    cmd.png(f"{out_stem}_overview_y90.png", dpi=300)

    for suffix in ("_interface", "_overview", "_overview_y90"):
        autocrop(f"{out_stem}{suffix}.png")

    print(f"wrote {out_stem}_interface.png, _overview.png, _overview_y90.png")


# PyMOL execs scripts with __name__ == "pymol", not "__main__", so guarding on
# "__main__" silently does nothing and the run exits 0 having rendered nothing.
if __name__ in ("__main__", "pymol"):
    args = sys.argv[1:]
    pdb = args[0] if args else "results/structures/complex/tnf3_01_complex.pdb"
    stem = args[1] if len(args) > 1 else "figures/tnf3_01"
    muts = args[2] if len(args) > 2 else ""
    render(pdb, stem, muts)
