"""Check the charge split of ``demo_charges`` on the octanoyl conjugate.

Run it from the repository root through pixi::

    pixi run python verify_charges.py

The script repeats notebook 02 without a notebook: it deprotonates the lysine
side chain, attaches the octanoyl fragment to the protein in memory, builds
the pablo residue library from the bond record that ``attach`` wrote, loads
the two topologies, runs the functions of ``demo_charges``, and reports the
timing and the result of every check. The site loses its proton first,
because the protonated amine carries no lone pair and is not the reactive
form of the side chain. The neutral amine reacts, so the product is a neutral
secondary amide. The conjugate PDB file goes to a scratch directory, so the
repository files do not change. The script exits with status 1 when a check
fails.
"""

import os
import sys
import time

import numpy as np
from openff.pablo import STD_CCD_CACHE, ResidueDefinition, topology_from_pdb
from openff.toolkit import Molecule
from openff.units import unit
from rdkit import Chem

from demo_charges import (
    assign_split_charges,
    atom_key,
    build_local_model,
    library_charge_map,
    parameterize_with_preset_charges,
)

UNMODIFIED_PDB = "1ubq_protonated.pdb"
FRAGMENT_SMILES = "*C(=O)CCCCCCC"
# The wwPDB Chemical Component Dictionary does not assign the code OC8,
# so no CCD component can take the place of this fragment definition.
FRAGMENT_RESNAME = "OC8"
ATOM_NAME = "NZ"
RESNUM = 63
CHAIN_ID = "A"
SCRATCH_DIR = (
    "/tmp/claude-1000/-home-joelaforet-Shirts-Lab-Linux-mbuild/"
    "e079c359-8ac9-4840-993c-a2a01a7ce492/scratchpad/vcfix"
)
CONJUGATE_PDB = os.path.join(SCRATCH_DIR, "1ubq_octanoyl.pdb")

results = []


def check(name, passed, detail):
    """Record the result of one check and print it."""
    results.append(passed)
    print(f"[{'pass' if passed else 'FAIL'}] {name}: {detail}")


def build_conjugate():
    """Attach the fragment in memory and write the conjugate to scratch.

    ``Protein.deprotonate`` runs before ``attach``, as in notebook 02. The
    lysine side chain carries an ammonium group at pH 7, and that group does
    not attack the acyl carbon. The site therefore loses one proton first.
    The net formal charge of the protein falls from 0 to -1 at that step, and
    the new bond takes the place of an N-H bond, so the nitrogen keeps a
    formal charge of zero and the net charge stays at -1.

    ``relax=False`` matches the notebook: the charges read the molecular
    graph only, so the placement of the fragment does not change them.

    Returns
    -------
    fragment : mbuild.Compound
        The pristine fragment. Its atom names transfer to the pablo residue
        definition by position.
    record : dict
        The first entry of ``Protein.bond_records()``.
    """
    from mbuild.biopolymers import Protein, prepare_fragment

    protein = Protein(UNMODIFIED_PDB)
    protein.deprotonate(RESNUM, ATOM_NAME, chain_id=CHAIN_ID)
    fragment = prepare_fragment(FRAGMENT_SMILES, FRAGMENT_RESNAME)
    protein.attach(
        fragment,
        resnum=RESNUM,
        atom_name=ATOM_NAME,
        chain_id=CHAIN_ID,
        relax=False,
    )
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    protein.save_pdb(CONJUGATE_PDB, overwrite=True)
    return fragment, protein.bond_records()[0]


def residue_library(fragment, record):
    """Return the pablo residue library that reads the octanoyl conjugate.

    The definition comes from the same starred SMILES that mBuild used, so
    that the atom order matches. The names come from the pristine fragment,
    because ``attach`` removes the leaving hydrogen from the copy it bonds to
    the protein. The crosslink comes from the bond record, the way notebook
    02 builds it.
    """
    star = Chem.RWMol(Chem.MolFromSmiles(FRAGMENT_SMILES))
    for atom in star.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
    molecule = star.GetMol()
    Chem.SanitizeMol(molecule)
    offmol = Molecule.from_rdkit(Chem.AddHs(molecule), allow_undefined_stereo=True)
    for atom, particle in zip(offmol.atoms, fragment.particles()):
        atom.name = particle.name

    return STD_CCD_CACHE.with_(
        {
            FRAGMENT_RESNAME: [
                ResidueDefinition.from_molecule(offmol, residue_name=FRAGMENT_RESNAME)
            ]
        }
    ).with_crosslink(
        residues=list(record["residue_names"]),
        linking_atoms=list(record["atom_names"]),
        leaving_atoms=[list(side) for side in record["leaving_atoms"]],
        bond_order=record["bond_order"],
    )


def main():
    start = time.time()
    fragment, record = build_conjugate()
    print(
        f"attachment: {time.time() - start:.1f} s | "
        f"{record['residue_names'][0]}{record['residue_numbers'][0]} "
        f"{record['atom_names'][0]} - {record['residue_names'][1]} "
        f"{record['atom_names'][1]} | leaving {record['leaving_atoms']}"
    )

    start = time.time()
    library = residue_library(fragment, record)
    print(f"residue library: {time.time() - start:.1f} s")

    start = time.time()
    unmodified_topology = topology_from_pdb(UNMODIFIED_PDB)
    conjugate_topology = topology_from_pdb(CONJUGATE_PDB, residue_library=library)
    conjugate = conjugate_topology.molecule(0)
    print(
        f"topologies: {time.time() - start:.1f} s | "
        f"{unmodified_topology.n_atoms} atoms unmodified | "
        f"{conjugate.n_atoms} atoms conjugate, formal charge "
        f"{conjugate.total_charge.m_as(unit.elementary_charge):+.0f}"
    )

    start = time.time()
    local, kept_indices, residue, fragment_atoms = build_local_model(
        conjugate, FRAGMENT_RESNAME, RESNUM, CHAIN_ID
    )
    site = residue | fragment_atoms
    print(
        f"local model: {time.time() - start:.1f} s | {local.n_atoms} atoms "
        f"({len(kept_indices)} cut out, {local.n_atoms - len(kept_indices)} caps) | "
        f"formal charge {local.total_charge.m_as(unit.elementary_charge):+.0f}"
    )
    smiles = local.to_smiles(explicit_hydrogens=False)
    local_charge = local.total_charge.m_as(unit.elementary_charge)
    print(f"local model SMILES: {smiles}")
    check(
        "the site is a neutral secondary amide",
        local_charge == 0.0 and "N+" not in smiles,
        f"{local.n_atoms} atoms, formal charge {local_charge:+.0f} e, "
        "no positively charged nitrogen in the SMILES",
    )

    reference = library_charge_map(unmodified_topology)

    start = time.time()
    charged = assign_split_charges(
        conjugate, unmodified_topology, FRAGMENT_RESNAME, RESNUM, CHAIN_ID
    )
    print(f"split charges: {time.time() - start:.1f} s")

    formal = charged.total_charge.m_as(unit.elementary_charge)
    charges = charged.partial_charges.m_as(unit.elementary_charge)
    check(
        "net charge",
        abs(charges.sum() - formal) <= 1e-6,
        f"{charges.sum():.9f} e against the formal {formal:+.0f} e",
    )

    deltas = [
        abs(charges[index] - reference[atom_key(atom)])
        for index, atom in enumerate(conjugate.atoms)
        if index not in site
    ]
    check(
        "unmodified residues keep the ff14SB charges",
        max(deltas) == 0.0,
        f"largest change {max(deltas):.1e} e over {len(deltas)} atoms",
    )

    start = time.time()
    interchange = parameterize_with_preset_charges(conjugate_topology, charged)
    print(f"interchange: {time.time() - start:.1f} s")
    written = np.array(
        [
            value.m_as(unit.elementary_charge)
            for _, value in sorted(
                interchange["Electrostatics"].charges.items(),
                key=lambda item: item[0].atom_indices[0],
            )
        ]
    )
    check(
        "the force field keeps the preset charges",
        np.array_equal(written[: charged.n_atoms], charges),
        f"{charged.n_atoms} atoms equal to the preset array",
    )

    return 0 if all(results) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:  # a failed check of demo_charges raises
        print(f"[FAIL] {type(error).__name__}: {error}")
        sys.exit(1)
