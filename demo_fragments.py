"""Build an mBuild fragment from a Chemical Component Dictionary code.

``mbuild.biopolymers.CCDLibrary`` reads a CCD component into a residue
template: atom names, elements, formal charges, bonds and bond orders.
It is a matcher, so it keeps no coordinates. The component file it
reads does carry them, in the ``pdbx_model_Cartn_*_ideal`` columns, and
this module picks them up to build the fragment.

The point of going through the CCD rather than through SMILES is the
atom names. A fragment built from SMILES gets invented names, so a
downstream residue library needs a hand-written definition to recognise
it. A fragment built from its CCD component already has the names that
library uses, and needs nothing.

This is a thin layer over ``CCDLibrary``, and the obvious next thing to
add to ``mbuild.biopolymers`` itself.
"""

from mbuild import Compound
from mbuild.biopolymers import CCDLibrary
from mbuild.biopolymers.ccd import _parse_cif_blocks, _cif_category_rows
from mbuild.biopolymers.residue import Residue

# CIF coordinates are in Angstrom; mBuild works in nanometres.
ANGSTROM_TO_NM = 0.1


def ideal_coordinates(code, library=None):
    """Return ``{atom name: (x, y, z)}`` in nm for one CCD component."""
    library = library or CCDLibrary(download=True)
    text = _component_cif(code, library)
    keys, loops = _parse_cif_blocks(text)
    rows = _cif_category_rows(keys, loops, "_chem_comp_atom")
    return {
        row["atom_id"]: tuple(
            float(row[f"pdbx_model_Cartn_{axis}_ideal"]) * ANGSTROM_TO_NM
            for axis in "xyz"
        )
        for row in rows
    }


def _component_cif(code, library):
    """Return the raw CIF text for one component, downloading if needed."""
    library[code]  # populates the cache, and raises if the code is unknown
    for directory in library._paths:
        path = directory / f"{code.upper()}.cif"
        if path.exists():
            return path.read_text()
    raise KeyError(f"No cached CCD file for {code!r}.")


def fragment_from_ccd(code, link_atom, library=None):
    """Return a fragment Residue for one CCD component.

    Parameters
    ----------
    code : str
        The component's CCD code, for example ``"KUT"``.
    link_atom : str
        Name of the atom that will bond to the protein. It is recorded
        on the residue, so ``Protein.attach`` finds it without being
        told again.
    library : CCDLibrary, optional
        The template library to read from.

    Returns
    -------
    Residue
        Atoms named as the CCD names them, carrying the component's
        formal charges, bonds and bond orders, placed at its ideal
        coordinates.
    """
    library = library or CCDLibrary(download=True)
    template = library[code][0]
    positions = ideal_coordinates(code, library)

    residue = Residue(resname=code.upper(), resnum=1, hetatm=True)
    particles = {}
    for atom in template.atoms:
        particle = Compound(
            name=atom.name, element=atom.element, pos=positions[atom.name]
        )
        particles[atom.name] = particle
        residue.add(particle)
    for bond in template.bonds:
        residue.add_bond(
            (particles[bond.atom1], particles[bond.atom2]),
            bond_order=float(bond.order),
        )

    residue.atom_formal_charges = {
        atom.name: atom.formal_charge
        for atom in template.atoms
        if atom.formal_charge
    }
    residue.formal_charge = sum(residue.atom_formal_charges.values())
    if link_atom not in particles:
        raise KeyError(
            f"{code} has no atom named {link_atom!r}. Its atoms are "
            f"{sorted(particles)}."
        )
    residue.link_atoms = {"1": link_atom}
    return residue
