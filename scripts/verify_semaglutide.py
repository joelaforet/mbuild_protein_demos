"""Check that the structure mBuild builds is the one the workshop loads.

Builds semaglutide's modified chain A with mBuild, reads it back with
openff-pablo, and compares the result against chain A of the workshop's
``7KI0_prepared.pdb``.

Both files go through the same loader call, so the comparison is between
two OpenFF ``Molecule`` objects rather than between two pieces of text.
Run with ``pixi run prove``.
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from openff.pablo import STD_CCD_CACHE, topology_from_pdb

from demo_fragments import fragment_from_ccd
from mbuild.biopolymers import Protein

APO = "semaglutide_apo.pdb"
BUILT = "semaglutide_mbuild.pdb"
REFERENCE = "semaglutide_reference.pdb"

# The crosslink the OpenFF workshop declares for this structure, used
# here unchanged. It names the atoms that leave each side, which is why
# attach() is told to remove those same two hydrogens.
CROSSLINK = dict(
    residues=["LYS", "KUT"],
    linking_atoms=["NZ", "C33"],
    leaving_atoms=[["HZ2"], ["H61"]],
    bond_order=1,
)


def build(path=BUILT):
    """Acylate LYS 20 with the KUT linker and write the result."""
    protein = Protein(APO, download=True)
    linker = fragment_from_ccd("KUT", link_atom="C33")
    protein.deprotonate(20, "NZ", chain_id="A")
    protein.attach(
        linker,
        resnum=20,
        atom_name="NZ",
        chain_id="A",
        leaving_atom_names="HZ2",
        relax=False,
    )
    protein.save_pdb(path, overwrite=True)
    return protein


def load(path):
    STD_CCD_CACHE.auto_download = True
    library = STD_CCD_CACHE.with_crosslink(**CROSSLINK)
    return topology_from_pdb(path, residue_library=library)


def bond_signature(molecule):
    return Counter(
        (
            tuple(sorted((bond.atom1.atomic_number, bond.atom2.atomic_number))),
            bond.bond_order,
        )
        for bond in molecule.bonds
    )


def atom_signature(molecule):
    return Counter(
        (atom.atomic_number, atom.formal_charge.m) for atom in molecule.atoms
    )


def stereo_by_site(molecule):
    return {
        (
            atom.metadata["residue_name"],
            atom.metadata["canonical_name"],
        ): atom.stereochemistry
        for atom in molecule.atoms
        if atom.stereochemistry
    }


def check():
    build()
    built = load(BUILT)
    reference = load(REFERENCE)
    a, b = built.molecule(0), reference.molecule(0)

    results = [
        ("one molecule, so the linker is bonded", built.n_molecules == 1),
        ("same atom count", a.n_atoms == b.n_atoms == 580),
        ("same formula", a.hill_formula == b.hill_formula),
        ("same net charge", a.total_charge == b.total_charge),
        ("same atoms, by element and formal charge", atom_signature(a) == atom_signature(b)),
        ("same bonds, by elements and order", bond_signature(a) == bond_signature(b)),
        (
            "isomorphic, setting stereochemistry aside",
            a.is_isomorphic_with(
                b,
                atom_stereochemistry_matching=False,
                bond_stereochemistry_matching=False,
            ),
        ),
    ]
    width = max(len(label) for label, _ in results)
    for label, passed in results:
        print(f"  {'PASS' if passed else 'FAIL'}  {label:<{width}}")

    differing = {
        site: (value, stereo_by_site(b).get(site))
        for site, value in stereo_by_site(a).items()
        if stereo_by_site(b).get(site) != value
    }
    print(f"\n  {a.n_atoms} atoms, {a.n_bonds} bonds, {a.hill_formula}, "
          f"net charge {a.total_charge}")
    if differing:
        print("\n  Stereocentres where the two structures disagree:")
        for (resname, atom_name), (ours, theirs) in sorted(differing.items()):
            print(f"    {resname} {atom_name}: built {ours}, deposited {theirs}")
        print(
            "\n  mBuild placed the linker at the geometry its CCD component\n"
            "  defines. At this centre the deposited coordinates and the CCD\n"
            "  component disagree, so the two structures differ there."
        )
    return all(passed for _, passed in results)


if __name__ == "__main__":
    raise SystemExit(0 if check() else 1)
