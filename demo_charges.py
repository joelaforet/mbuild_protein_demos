"""Split the partial charges of a modified protein between two charge models.

The protein atoms keep the Amber ff14SB library charges. The attached
fragment, and the residue it is bonded to, take graph charges computed on a
small capped model of the modification site. The two sets are stitched into
one per-atom array on the conjugate molecule, and the array is handed to
Interchange as preset charges.

The split is needed because the ff14SB library charges cover standard
residues only. The conjugate holds one non-standard residue. The coverage
test, which needs a library charge for every atom, fails. Sage 2.3.0 then
falls back to its own NAGLCharges handler for the whole protein.

The charge model for the fragment is NAGL am1bcc graph charges
(``openff-gnn-am1bcc-0.1.0-rc.3.pt``).

For the remaining parameters, ``parameterize_with_preset_charges`` copies
ff14SB terms wholly outside the site from the unmodified reference. All
terms touching the modified residue or fragment use Sage, including terms
across the peptide boundaries.

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


def atom_key(atom):
    """Return the chain id, residue number, insertion code and name of an atom.

    The key pairs the atoms of two topologies. A positional pairing is wrong,
    because the fragment atoms shift every later index.
    """
    metadata = atom.metadata
    return (
        metadata["chain_id"],
        metadata["residue_number"],
        metadata["insertion_code"],
        atom.name,
    )


def library_charge_map(topology):
    """Return the ff14SB library charge of every atom of a protein topology.

    The topology holds standard residues only. The map goes from the key of
    :func:`atom_key` to the charge in elementary charge. Two atoms with one
    key raise ``ValueError``.
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


def _label(atom):
    """Return the chain, residue and atom name of an atom, for a message."""
    data = atom.metadata
    return (
        f"chain {data['chain_id']} "
        f"{data['residue_name']}{data['residue_number']} {atom.name}"
    )


def build_local_model(conjugate, fragment_resname, resnum, chain_id):
    """Cut a capped model of the modification site out of the conjugate.

    The selection holds the modified residue and the fragment. It grows by
    two bonds through the molecule, so that the atoms of the modification
    keep their real first and second neighbours. A hydrogen atom replaces
    every bond that leaves the selection. The model carries no conformer,
    because NAGL reads the molecular graph only.

    Parameters
    ----------
    conjugate : openff.toolkit.Molecule
        The modified protein.
    fragment_resname : str
        Residue name of the fragment.
    resnum, chain_id : int, str
        Residue number and chain of the modified residue. Its insertion
        code must be blank.

    Returns
    -------
    local : openff.toolkit.Molecule
        The capped model. The cap hydrogens follow the cut-out atoms.
    kept_indices : list of int
        Index in ``conjugate`` of every cut-out atom, in the order of the
        atoms of ``local``.
    residue, fragment : set of int
        Indices in ``conjugate`` of the modified residue and of the fragment.

    Raises
    ------
    ValueError
        On an empty selection, a fragment bond outside the modified residue,
        or a cut bond of order above one.
    """
    residue = set()
    fragment = set()
    for index, atom in enumerate(conjugate.atoms):
        data = atom.metadata
        if (
            data["chain_id"] == chain_id
            and data["residue_number"] == resnum
            and str(data["insertion_code"]).strip() == ""
        ):
            residue.add(index)
        if data["residue_name"] == fragment_resname:
            fragment.add(index)
    if not residue:
        raise ValueError(
            f"no atom of the conjugate lies in chain {chain_id} residue {resnum}"
        )
    if not fragment:
        raise ValueError(f"the conjugate has no residue named {fragment_resname}")

    neighbours = [set() for _ in range(conjugate.n_atoms)]
    for bond in conjugate.bonds:
        neighbours[bond.atom1_index].add(bond.atom2_index)
        neighbours[bond.atom2_index].add(bond.atom1_index)

    attached = {other for index in fragment for other in neighbours[index]} - fragment
    if not attached <= residue:
        outside = ", ".join(
            _label(conjugate.atom(index)) for index in sorted(attached - residue)
        )
        raise ValueError(
            f"the fragment {fragment_resname} is bonded to {outside}, outside "
            f"chain {chain_id} residue {resnum}; the site residue number or the "
            "fragment residue name is wrong"
        )

    kept = residue | fragment
    for _ in range(2):
        kept = kept | {other for index in kept for other in neighbours[index]}
    kept_indices = sorted(kept)
    local_index = {index: order for order, index in enumerate(kept_indices)}

    local = Molecule()
    for index in kept_indices:
        atom = conjugate.atom(index)
        local.add_atom(
            atom.atomic_number,
            atom.formal_charge.m_as(unit.elementary_charge),
            atom.is_aromatic,
            name=atom.name,
        )

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
            cap = local.add_atom(1, 0, False, name="HC")
            local.add_bond(local_index[inside[0]], cap, 1, False)

    return local, kept_indices, residue, fragment


def _print_table(conjugate, residue, fragment, reference, graph_charge, charges):
    """Print one line per atom of the site, and name the largest change."""
    print(
        f"{'residue':>8} {'atom':>6} {'ff14SB':>9} "
        f"{'NAGL':>9} {'final':>9} {'delta':>9}"
    )
    largest = None
    for index in sorted(residue) + sorted(fragment):
        atom = conjugate.atom(index)
        data = atom.metadata
        label = f"{data['residue_name']}{data['residue_number']}"
        library = reference.get(atom_key(atom))
        if library is None:
            library_text = delta_text = "-"
        else:
            delta = graph_charge[index] - library
            library_text = format(library, "+.4f")
            delta_text = format(delta, "+.4f")
            if largest is None or abs(delta) > abs(largest[1]):
                largest = (f"{label} {atom.name}", delta)
        print(
            f"{label:>8} {atom.name:>6} {library_text:>9} "
            f"{graph_charge[index]:>+9.4f} {charges[index]:>+9.4f} {delta_text:>9}"
        )
    if largest is not None:
        print(f"largest change {largest[1]:+.4f} e on {largest[0]}")


def assign_split_charges(
    conjugate,
    unmodified_topology,
    fragment_resname,
    resnum,
    chain_id,
    nagl_model="openff-gnn-am1bcc-0.1.0-rc.3.pt",
    tolerance=1e-6,
):
    """Set the partial charges of the conjugate from the two charge models.

    The site holds the fragment and the modified residue. Its atoms take the
    graph charges of the local model, and every other atom takes its ff14SB
    library charge. The seam between the two models lies on the peptide bonds
    of the modified residue. The modified residue takes graph charges too,
    because the attachment changed its chemistry and the ff14SB charges of
    the standard residue no longer describe it.

    Parameters
    ----------
    conjugate : openff.toolkit.Molecule
        The modified protein. Its ``partial_charges`` are overwritten.
    unmodified_topology : openff.toolkit.Topology
        The protein before the modification. It supplies the ff14SB charges.
    fragment_resname : str
        Residue name of the fragment.
    resnum, chain_id : int, str
        Residue number and chain of the modified residue.
    nagl_model : str, optional
        NAGL model file. The default is the am1bcc graph model.
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
        On a shared atom key, a missing ff14SB charge, a non-integer charge
        sum outside the site, a residual per atom of the site above 0.005 e,
        or a net charge off the formal charge by more than ``tolerance``.
    """
    local, kept_indices, residue, fragment = build_local_model(
        conjugate, fragment_resname, resnum, chain_id
    )
    NAGLToolkitWrapper().assign_partial_charges(local, partial_charge_method=nagl_model)
    graph_charge = {
        index: local.partial_charges[order].m_as(unit.elementary_charge)
        for order, index in enumerate(kept_indices)
    }

    site = residue | fragment
    reference = library_charge_map(unmodified_topology)
    seen = set()
    charges = np.zeros(conjugate.n_atoms)
    for index, atom in enumerate(conjugate.atoms):
        key = atom_key(atom)
        if key in seen:
            raise ValueError(f"two atoms of the conjugate share the key {key}")
        seen.add(key)
        if index in site:
            charges[index] = graph_charge[index]
        elif key in reference:
            charges[index] = reference[key]
        else:
            raise ValueError(
                f"{_label(atom)} has no ff14SB charge; the unmodified topology "
                "does not hold it"
            )

    outside = np.ones(conjugate.n_atoms, dtype=bool)
    outside[sorted(site)] = False
    outside_charge = charges[outside].sum()
    if abs(outside_charge - round(outside_charge)) > 1e-4:
        raise ValueError(
            f"the charges outside the site sum to {outside_charge:+.6f} e, which "
            "is not an integer within 1e-4 e; an atom of a standard residue is "
            "missing from the charge map, or the site is wrong"
        )

    formal = conjugate.total_charge.m_as(unit.elementary_charge)
    residual = formal - charges.sum()
    smear = residual / len(site)
    if abs(smear) > 0.005:
        raise ValueError(
            f"the residual of {residual:+.6f} e over {len(site)} atoms leaves "
            f"{smear:+.6f} e per atom, which is more than the limit of 0.005 e; "
            "check the leaving atoms of the attachment and the atoms of the site"
        )
    for index in site:
        charges[index] += smear

    net = charges.sum()
    if abs(net - formal) > tolerance:
        raise ValueError(
            f"the net charge of {net:.9f} e misses the formal charge of "
            f"{formal:.1f} e by more than {tolerance:g} e"
        )

    site_name = conjugate.atom(min(residue)).metadata["residue_name"]
    print(f"charges of {site_name}{resnum} and {fragment_resname}")
    _print_table(conjugate, residue, fragment, reference, graph_charge, charges)
    print(f"residual {residual:+.6f} e over {len(site)} atoms")
    print(f"per-atom smear {smear:+.6f} e | net charge {net:.9f} e")

    conjugate.partial_charges = charges * unit.elementary_charge
    return conjugate


def parameterize_with_preset_charges(
    conjugate_topology,
    charged_conjugate,
    unmodified_topology,
    *,
    fragment_resname,
    resnum,
    chain_id,
    protein_forcefield="ff14sb_off_impropers_0.0.4.offxml",
    site_forcefield="openff-2.3.0.offxml",
):
    """Use ff14SB outside the site and Sage for every term touching the site.

    The site is the complete modified residue plus its attached fragment.
    Bonds, angles, and torsions crossing its boundary use Sage. Unmodified
    protein terms are copied from an independently parameterized reference,
    including every Fourier term and improper permutation. Solvent uses Sage's
    solvent parameters. Constraints on bonds to hydrogen follow Sage's constraint
    policy with the equilibrium lengths of the selected bond force field.

    Charges must already be assigned by ``assign_split_charges``. The
    reference topology and explicit site selection are required so force-field
    precedence cannot silently change the intended split. Returns Interchange.
    """
    if charged_conjugate.partial_charges is None:
        raise ValueError("the conjugate carries no partial charges")
    _, _, residue, fragment = build_local_model(
        charged_conjugate, fragment_resname, resnum, chain_id
    )
    site = residue | fragment

    # Locate the actual conjugate, not merely a molecule with the same size.
    wanted = [atom_key(atom) for atom in charged_conjugate.atoms]
    candidates = []
    for molecule in conjugate_topology.molecules:
        if molecule.n_atoms != charged_conjugate.n_atoms:
            continue
        if [atom_key(atom) for atom in molecule.atoms] != wanted:
            continue
        if molecule.is_isomorphic_with(charged_conjugate):
            candidates.append(molecule)
    if len(candidates) != 1:
        raise ValueError(
            "expected one conjugate with matching chemistry and atom order"
        )
    start = conjugate_topology.atom_index(candidates[0].atom(0))
    outside = {start + i for i in range(charged_conjugate.n_atoms) if i not in site}
    target_indices = {key: start + i for i, key in enumerate(wanted)}
    if len(target_indices) != len(wanted):
        raise ValueError("duplicate atom keys in the conjugate")

    reference_indices = {}
    seen = set()
    for i, atom in enumerate(unmodified_topology.atoms):
        key = atom_key(atom)
        if key in seen:
            raise ValueError(f"two reference atoms share the key {key}")
        seen.add(key)
        if key in target_indices:
            reference_indices[i] = target_indices[key]
    if not outside <= set(reference_indices.values()):
        raise ValueError("an unmodified protein atom is missing from the reference")

    sage = ForceField(site_forcefield)
    amber = ForceField(protein_forcefield)
    # The Amber port omits a Constraints handler. Apply the same bonds-to-hydrogen
    # constraint policy, resolving unspecified distances from Amber bonds.
    if "Constraints" not in amber.registered_parameter_handlers:
        amber.register_parameter_handler(sage["Constraints"])
    interchange = sage.create_interchange(
        conjugate_topology, charge_from_molecules=[charged_conjugate]
    )
    reference = amber.create_interchange(unmodified_topology)

    # Both force fields use Lorentz-Berthelot and the same exclusion/1-4
    # conventions (the Amber port rounds 5/6 to six decimal places).
    for name in ("vdW", "Electrostatics"):
        for attr in ("scale_12", "scale_13", "scale_14", "scale_15"):
            if not np.isclose(
                getattr(interchange[name], attr),
                getattr(reference[name], attr),
                atol=1e-6,
                rtol=0,
            ):
                raise ValueError(f"incompatible {name} {attr} between force fields")
    if interchange["vdW"].mixing_rule != reference["vdW"].mixing_rule:
        raise ValueError("incompatible Lennard-Jones mixing rules")

    for name in (
        "Bonds",
        "Angles",
        "ProperTorsions",
        "ImproperTorsions",
        "vdW",
        "Constraints",
    ):
        target = interchange[name]
        source = reference[name]
        # Remove all Sage terms wholly within the unmodified protein. Do not
        # replace torsions one key at a time: multiplicities can differ.
        removed = {key for key in target.key_map if set(key.atom_indices) <= outside}
        for key in removed:
            del target.key_map[key]
        copied = set()
        for key, potential_key in source.key_map.items():
            mapped = tuple(reference_indices.get(i) for i in key.atom_indices)
            if not set(mapped) <= outside:
                continue
            new_key = key.model_copy(update={"atom_indices": mapped})
            new_potential = potential_key.model_copy(
                update={"id": "ff14SB::" + potential_key.id}
            )
            target.key_map[new_key] = new_potential
            target.potentials[new_potential] = source.potentials[potential_key]
            copied.add(new_key)
        # Bond/angle/proper/vdW coverage must survive the transfer. Improper
        # definitions and optional constraints can differ between models.
        if name not in ("ImproperTorsions", "Constraints"):
            expected = {tuple(sorted(key.atom_indices)) for key in removed}
            actual = {tuple(sorted(key.atom_indices)) for key in copied}
            if actual != expected:
                raise ValueError(f"ff14SB reference does not cover unchanged {name}")
        used = set(target.key_map.values())
        target.potentials = {
            key: value for key, value in target.potentials.items() if key in used
        }
        print(
            f"{name}: {len(copied)} ff14SB terms | {len(target.key_map) - len(copied)} Sage/solvent terms"
        )

    preset = charged_conjugate.partial_charges.m_as(unit.elementary_charge)
    written = np.array(
        [
            charge.m_as(unit.elementary_charge)
            for key, charge in sorted(
                interchange["Electrostatics"].charges.items(),
                key=lambda item: item[0].atom_indices[0],
            )
        ]
    )[start : start + charged_conjugate.n_atoms]
    if not np.array_equal(written, preset):
        raise ValueError("the force field overwrote the preset charges")
    return interchange
