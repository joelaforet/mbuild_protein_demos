"""Split the partial charges of a modified protein between two charge models.

The protein atoms keep the Amber ff14SB library charges. The attached
fragment, and the residue it is bonded to, take graph charges computed on a
small capped model of the modification site. The two sets are stitched into
one per-atom array on the conjugate molecule, and the array is handed to
Interchange as preset charges.

The split is needed because the ff14SB library charges cover standard
residues only. The conjugate contains one non-standard residue, so the
all-or-nothing coverage test of the LibraryCharges handler fails and Sage
2.3.0 falls back to its own NAGLCharges handler for the whole protein.

The charge model for the fragment is NAGL am1bcc graph charges
(``openff-gnn-am1bcc-0.1.0-rc.3.pt``). It is not AshGC.

This module supports the demonstration notebooks. It is not part of mBuild.
"""

import numpy as np
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.units import unit

#: Handlers removed from ff14SB before the reference charges are computed.
#: Only LibraryCharges and Electrostatics are needed, and the removal cuts
#: the run to under one second.
_UNUSED_HANDLERS = (
    "Bonds",
    "Angles",
    "ProperTorsions",
    "ImproperTorsions",
    "Constraints",
    "vdW",
)

#: Bond length of the hydrogen atoms that cap the local model, in angstrom.
_CAP_BOND_LENGTH = 1.09


def atom_key(atom):
    """Return chain id, residue number, insertion code and name of an atom.

    The key pairs the atoms of two topologies. A positional pairing is
    wrong, because the fragment atoms shift every later index. Pablo writes
    no ``atom_name`` metadata, so the name comes from the atom.

    Parameters
    ----------
    atom : openff.toolkit.topology.Atom
        Atom with the PDB metadata that openff-pablo writes.

    Returns
    -------
    tuple
        The four values, in that order.
    """
    metadata = atom.metadata
    return (
        metadata["chain_id"],
        metadata["residue_number"],
        metadata["insertion_code"],
        atom.name,
    )


def library_charge_map(topology):
    """Return the ff14SB library charge of every atom in a protein topology.

    Parameters
    ----------
    topology : openff.toolkit.Topology
        Topology of standard residues only.

    Returns
    -------
    dict
        Map from :func:`atom_key` to charge in elementary charge. Two atoms
        with one key raise ``ValueError``.
    """
    force_field = ForceField("ff14sb_off_impropers_0.0.4.offxml")
    for handler in _UNUSED_HANDLERS:
        if handler in force_field.registered_parameter_handlers:
            force_field.deregister_parameter_handler(handler)

    charges = np.zeros(topology.n_atoms)
    for key, charge in force_field.create_interchange(topology)[
        "Electrostatics"
    ].charges.items():
        charges[key.atom_indices[0]] = charge.m_as(unit.elementary_charge)

    charge_map = {}
    for index, atom in enumerate(topology.atoms):
        key = atom_key(atom)
        if key in charge_map:
            raise ValueError(f"two atoms share the key {key}")
        charge_map[key] = charges[index]
    return charge_map


def _selected_indices(conjugate, site, fragment_resname):
    """Return the atom indices of the modified residue and the fragment."""
    residue = set()
    fragment = set()
    for index, atom in enumerate(conjugate.atoms):
        metadata = atom.metadata
        if (metadata["chain_id"], metadata["residue_number"]) == (
            site["chain_id"],
            site["residue_number"],
        ):
            residue.add(index)
        if metadata["residue_name"] == fragment_resname:
            fragment.add(index)
    return residue, fragment


def build_local_model(conjugate, site, fragment_resname):
    """Cut a capped model of the modification site out of the conjugate.

    The selection holds the modified residue and the fragment residue. It
    grows by two bonds through the molecule, so that the atoms of the
    modification keep their real first and second neighbours. Every bond
    that leaves the selection is replaced by a hydrogen atom on the bond
    vector.

    Parameters
    ----------
    conjugate : openff.toolkit.Molecule
        The modified protein, with one conformer.
    site : dict
        ``chain_id`` and ``residue_number`` of the modified residue.
    fragment_resname : str
        Residue name of the fragment.

    Returns
    -------
    local : openff.toolkit.Molecule
        The capped model, with one conformer.
    kept_indices : list of int
        Index in ``conjugate`` of every atom of ``local`` that was cut out,
        in the order of the atoms of ``local``. The cap hydrogens follow
        these atoms and have no counterpart in ``conjugate``.

    Raises
    ------
    ValueError
        If the selection is empty, or if a bond of order greater than one
        leaves the selection.
    """
    residue, fragment = _selected_indices(conjugate, site, fragment_resname)
    if not residue:
        raise ValueError(f"no atom of the conjugate matches the site {site}")
    if not fragment:
        raise ValueError(f"the conjugate has no residue named {fragment_resname}")

    neighbours = [set() for _ in range(conjugate.n_atoms)]
    for bond in conjugate.bonds:
        first, second = bond.atom1_index, bond.atom2_index
        neighbours[first].add(second)
        neighbours[second].add(first)

    kept = residue | fragment
    for _ in range(2):
        kept = kept | {other for index in kept for other in neighbours[index]}
    kept_indices = sorted(kept)
    local_index = {index: order for order, index in enumerate(kept_indices)}

    conformer = conjugate.conformers[0].m_as(unit.angstrom)
    local = Molecule()
    positions = []
    for index in kept_indices:
        atom = conjugate.atom(index)
        local.add_atom(
            atom.atomic_number,
            atom.formal_charge.m_as(unit.elementary_charge),
            atom.is_aromatic,
            name=atom.name,
        )
        positions.append(conformer[index])

    for bond in conjugate.bonds:
        ends = (bond.atom1_index, bond.atom2_index)
        inside = [index for index in ends if index in kept]
        if len(inside) == 2:
            local.add_bond(
                local_index[bond.atom1_index],
                local_index[bond.atom2_index],
                bond.bond_order,
                bond.is_aromatic,
            )
        elif len(inside) == 1:
            if bond.bond_order != 1:
                raise ValueError(
                    f"the selection cuts a bond of order {bond.bond_order} "
                    f"between atoms {bond.atom1_index} and {bond.atom2_index}"
                )
            (anchor,) = inside
            outside = (
                bond.atom2_index if anchor == bond.atom1_index else bond.atom1_index
            )
            cap = local.add_atom(1, 0, False, name="HC")
            local.add_bond(local_index[anchor], cap, 1, False)
            vector = conformer[outside] - conformer[anchor]
            direction = vector / np.linalg.norm(vector)
            positions.append(conformer[anchor] + _CAP_BOND_LENGTH * direction)

    local.add_conformer(np.array(positions) * unit.angstrom)
    return local, kept_indices


def assign_split_charges(
    conjugate,
    unmodified_topology,
    fragment_resname,
    site,
    nagl_model="openff-gnn-am1bcc-0.1.0-rc.3.pt",
    nagl_scope="fragment_and_residue",
    tolerance=1e-6,
):
    """Set the partial charges of the conjugate from the two charge models.

    Parameters
    ----------
    conjugate : openff.toolkit.Molecule
        The modified protein. Its ``partial_charges`` are overwritten.
    unmodified_topology : openff.toolkit.Topology
        The protein before the modification. It supplies the ff14SB charges.
    fragment_resname : str
        Residue name of the fragment.
    site : dict
        ``chain_id`` and ``residue_number`` of the modified residue.
    nagl_model : str, optional
        NAGL model file. The default is the am1bcc graph model.
    nagl_scope : {'fragment_and_residue', 'fragment'}, optional
        Atoms that take graph charges. The default also moves the modified
        residue, because the acylation changes its electron distribution.
    tolerance : float, optional
        Largest accepted difference between the net charge and the formal
        charge, in elementary charge.

    Returns
    -------
    openff.toolkit.Molecule
        The conjugate, with ``partial_charges`` set.

    Raises
    ------
    ValueError
        If ``nagl_scope`` is unknown, if an out-of-scope atom has no ff14SB
        charge, if the residual per atom is larger than 0.005 e, or if the
        net charge misses the formal charge by more than ``tolerance``.
    """
    if nagl_scope not in ("fragment_and_residue", "fragment"):
        raise ValueError(f"nagl_scope must name a known scope, not {nagl_scope!r}")

    local, kept_indices = build_local_model(conjugate, site, fragment_resname)
    NAGLToolkitWrapper().assign_partial_charges(
        local, partial_charge_method=nagl_model
    )
    graph_charge = {
        index: local.partial_charges[order].m_as(unit.elementary_charge)
        for order, index in enumerate(kept_indices)
    }

    residue, fragment = _selected_indices(conjugate, site, fragment_resname)
    scope = fragment | residue if nagl_scope == "fragment_and_residue" else fragment

    reference = library_charge_map(unmodified_topology)
    keys = {}
    charges = np.zeros(conjugate.n_atoms)
    for index, atom in enumerate(conjugate.atoms):
        key = atom_key(atom)
        if key in keys:
            raise ValueError(f"two atoms of the conjugate share the key {key}")
        keys[key] = index
        if index in scope:
            charges[index] = graph_charge[index]
            continue
        if key not in reference:
            raise ValueError(
                f"atom {atom.name} of residue "
                f"{atom.metadata['residue_name']}{atom.metadata['residue_number']} "
                f"has no ff14SB charge; the unmodified topology does not hold it"
            )
        charges[index] = reference[key]

    formal = conjugate.total_charge.m_as(unit.elementary_charge)
    residual = formal - charges.sum()
    smear = residual / len(scope)
    if abs(smear) > 0.005:
        raise ValueError(
            f"the residual of {residual:+.6f} e over {len(scope)} atoms leaves "
            f"{smear:+.6f} e per atom, which is more than the limit of 0.005 e"
        )
    for index in scope:
        charges[index] += smear

    net = charges.sum()
    if abs(net - formal) > tolerance:
        raise ValueError(
            f"the net charge of {net:.9f} e misses the formal charge of "
            f"{formal:.1f} e by more than {tolerance:g} e"
        )
    for key, index in keys.items():
        if index in scope:
            continue
        if charges[index] != reference[key]:
            raise ValueError(
                f"atom {key} is out of scope but its charge moved by "
                f"{charges[index] - reference[key]:+.9f} e"
            )

    label = f"{site['chain_id']}{site['residue_number']}"
    print(f"charges of {label} and {fragment_resname}")
    print(f"{'atom':>6} {'ff14SB':>9} {'NAGL':>9} {'final':>9}")
    for index in sorted(residue):
        atom = conjugate.atom(index)
        key = atom_key(atom)
        library = reference.get(key)
        print(
            f"{atom.name:>6} "
            f"{'--' if library is None else format(library, '+.4f'):>9} "
            f"{graph_charge[index]:>+9.4f} {charges[index]:>+9.4f}"
        )
    print(f"residual {residual:+.6f} e over {len(scope)} atoms")
    print(f"per-atom smear {smear:+.6f} e | net charge {net:.9f} e")

    conjugate.partial_charges = charges * unit.elementary_charge
    return conjugate


def build_interchange(
    conjugate_topology,
    charged_conjugate,
    forcefield_names=("ff14sb_off_impropers_0.0.4.offxml", "openff-2.3.0.offxml"),
):
    """Parameterize a topology and keep the preset charges of the conjugate.

    Parameters
    ----------
    conjugate_topology : openff.toolkit.Topology
        Topology that holds the conjugate, and any other molecule.
    charged_conjugate : openff.toolkit.Molecule
        The conjugate with ``partial_charges`` set.
    forcefield_names : tuple of str, optional
        Force field files, in the order Interchange reads them.

    Returns
    -------
    openff.interchange.Interchange
        The parameterized system.

    Raises
    ------
    ValueError
        If the conjugate has no partial charges, if the topology does not
        hold it, or if the charges of the system differ from the preset
        array. The last check proves that the NAGLCharges handler of Sage
        2.3.0 did not write over the split charges.
    """
    if charged_conjugate.partial_charges is None:
        raise ValueError("the conjugate carries no partial charges")

    interchange = ForceField(*forcefield_names).create_interchange(
        conjugate_topology, charge_from_molecules=[charged_conjugate]
    )

    start = None
    for molecule in conjugate_topology.molecules:
        same_size = molecule.n_atoms == charged_conjugate.n_atoms
        if molecule is charged_conjugate or same_size:
            start = conjugate_topology.atom_index(molecule.atom(0))
            break
    if start is None:
        raise ValueError("the topology holds no molecule of the size of the conjugate")

    preset = charged_conjugate.partial_charges.m_as(unit.elementary_charge)
    charges = interchange["Electrostatics"].charges
    written = np.array(
        [
            charges[key].m_as(unit.elementary_charge)
            for key in sorted(charges, key=lambda key: key.atom_indices[0])
        ]
    )[start : start + charged_conjugate.n_atoms]
    if not np.array_equal(written, preset):
        raise ValueError(
            "the force field overwrote the preset charges; the largest "
            f"difference is {np.abs(written - preset).max():.6f} e"
        )
    return interchange
