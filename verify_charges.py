"""Check the charge split of ``demo_charges`` on the octanoyl conjugate.

Run it from the repository root through pixi::

    pixi run python verify_charges.py

The script builds the pablo residue library the way notebook 02 builds it,
loads the unmodified protein and the conjugate, runs the three functions of
``demo_charges``, and reports the timing and the result of every check. It
exits with status 1 when a check fails.
"""

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
    build_interchange,
    build_local_model,
    library_charge_map,
)

UNMODIFIED_PDB = "1ubq_protonated.pdb"
CONJUGATE_PDB = "1ubq_octanoyl.pdb"
FRAGMENT_SMILES = "*C(=O)CCCCCCC"
FRAGMENT_RESNAME = "OCT"
SITE = {"chain_id": "A", "residue_number": 63}

results = []


def check(name, passed, detail):
    """Record the result of one check and print it."""
    results.append(passed)
    print(f"[{'pass' if passed else 'FAIL'}] {name}: {detail}")


def residue_library():
    """Return the pablo residue library that reads the octanoyl conjugate.

    The definition comes from the same starred SMILES that mBuild used, so
    that the atom order matches. The names come from the pristine fragment
    that ``prepare_fragment`` returns, because ``attach`` removes the
    leaving hydrogen from the copy it bonds to the protein.
    """
    from mbuild.biopolymers import prepare_fragment

    star = Chem.RWMol(Chem.MolFromSmiles(FRAGMENT_SMILES))
    for atom in star.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
    molecule = star.GetMol()
    Chem.SanitizeMol(molecule)
    offmol = Molecule.from_rdkit(Chem.AddHs(molecule), allow_undefined_stereo=True)
    for atom, particle in zip(
        offmol.atoms, prepare_fragment(FRAGMENT_SMILES, FRAGMENT_RESNAME).particles()
    ):
        atom.name = particle.name

    return STD_CCD_CACHE.with_(
        {
            FRAGMENT_RESNAME: [
                ResidueDefinition.from_molecule(offmol, residue_name=FRAGMENT_RESNAME)
            ]
        }
    ).with_crosslink(
        residues=("LYS", FRAGMENT_RESNAME),
        linking_atoms=("NZ", "C1"),
        leaving_atoms=[["HZ1"], ["H1"]],
        bond_order=1,
    )


def main():
    start = time.time()
    library = residue_library()
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
    local, kept_indices = build_local_model(conjugate, SITE, FRAGMENT_RESNAME)
    print(
        f"local model: {time.time() - start:.1f} s | {local.n_atoms} atoms "
        f"({len(kept_indices)} cut out, {local.n_atoms - len(kept_indices)} caps) | "
        f"formal charge {local.total_charge.m_as(unit.elementary_charge):+.0f}"
    )
    print(f"local model SMILES: {local.to_smiles(explicit_hydrogens=False)}")

    start = time.time()
    charged = assign_split_charges(
        conjugate, unmodified_topology, FRAGMENT_RESNAME, SITE
    )
    print(f"split charges: {time.time() - start:.1f} s")

    formal = charged.total_charge.m_as(unit.elementary_charge)
    charges = charged.partial_charges.m_as(unit.elementary_charge)
    check(
        "net charge",
        abs(charges.sum() - formal) <= 1e-6,
        f"{charges.sum():.9f} e against the formal {formal:+.0f} e",
    )

    reference = library_charge_map(unmodified_topology)
    moved = {
        index
        for index, atom in enumerate(conjugate.atoms)
        if atom.metadata["residue_name"] == FRAGMENT_RESNAME
        or (atom.metadata["chain_id"], atom.metadata["residue_number"])
        == (SITE["chain_id"], SITE["residue_number"])
    }
    deltas = [
        abs(charges[index] - reference[atom_key(atom)])
        for index, atom in enumerate(conjugate.atoms)
        if index not in moved
    ]
    check(
        "unmodified residues keep the ff14SB charges",
        max(deltas) == 0.0,
        f"largest change {max(deltas):.1e} e over {len(deltas)} atoms",
    )

    start = time.time()
    interchange = build_interchange(conjugate_topology, charged)
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
