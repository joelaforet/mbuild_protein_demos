"""Generate the demo notebooks.

Keeping the notebooks in a script means the prose and the code stay
reviewable in a normal diff. Run with ``pixi run notebooks``, then
``pixi run verify`` to execute them.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")}


def code(text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.strip("\n"),
    }


def notebook(cells):
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


ONE = [
    md(
        """
# Reading a protein into mBuild with its chemistry intact

A PDB file gives you element symbols and coordinates. To assign force
field parameters you need more than that: the bond orders and the formal
charge of every atom. Those are not in the file.

Most readers fill the gap by guessing — bonds from interatomic
distances, charges left at zero or copied from a column that is usually
blank. The guess is invisible until something downstream is wrong.

`mbuild.biopolymers` does not guess. Every residue is matched by atom
name against a template from the wwPDB Chemical Component Dictionary,
and the bonds, bond orders and formal charges come from that template. A
residue no template explains raises an error naming the residue.
"""
    ),
    code(
        """
from mbuild.biopolymers import Protein

protein = Protein("1ubq_protonated.pdb")
print(len(list(protein.residues())), "residues,", protein.n_particles, "atoms")
print("net formal charge:", protein.net_formal_charge)
"""
    ),
    md(
        """
The charge is the part a coordinate reader cannot give you. It comes
from the matched templates, one residue at a time.
"""
    ),
    code(
        """
for resnum in (1, 48, 76):
    residue = protein.get_residue(resnum, chain_id="A")
    print(f"{residue.name} {resnum:>3}  charge {residue.formal_charge:+d}  "
          f"{residue.atom_formal_charges}")
"""
    ),
    md(
        """
Bonds carry orders, so the structure is a chemical graph rather than a
set of connected points. `to_rdkit` is the shortest way to see that: it
refuses to export a bond whose order is unknown, so the fact that it
returns a sanitized molecule at all is the check.
"""
    ),
    code(
        """
from rdkit import Chem

mol = protein.to_rdkit()
print(mol.GetNumAtoms(), "atoms, formal charge", Chem.GetFormalCharge(mol))
"""
    ),
    md(
        """
## Handing it to OpenFF

`save_pdb` writes a file with real residue numbers, chain identifiers,
a TER after each chain, and CONECT records for the bonds that residue
adjacency cannot imply — disulfides here. Peptide bonds are left
implied, because a residue-template reader rejects a CONECT its own
definitions cannot explain.

openff-pablo reads that file with no extra arguments.
"""
    ),
    code(
        """
from openff.pablo import STD_CCD_CACHE, topology_from_pdb

protein.save_pdb("1ubq_prepared.pdb", overwrite=True)
topology = topology_from_pdb("1ubq_prepared.pdb", residue_library=STD_CCD_CACHE)

molecule = topology.molecule(0)
print(molecule.n_atoms, "atoms, net charge", molecule.total_charge)
"""
    ),
    code(
        """
view = topology.visualize()
view.clear_representations()
view.add_representation("cartoon", color="#990000")
view
"""
    ),
    md(
        """
That is the round trip for an unmodified protein. The next notebook
does the part that needed new code: modifying the protein first, and
still handing OpenFF something it can read.
"""
    ),
]


TWO = [
    md(
        """
# Building a modified protein that OpenFF can read

Semaglutide is a 31-residue peptide with a lipid linker acylating one
lysine. The OpenFF post-translational-modification workshop simulates
it, starting from a prepared PDB file that someone had to make first.

This notebook makes that file. It starts from the unmodified peptide,
attaches the linker with mBuild, and reads the result back with
openff-pablo using the workshop's own loader call, unchanged.
"""
    ),
    md(
        """
## The starting structure

An ordinary protonated peptide. Residue 2 is Aib, and residue 20 is the
lysine that will carry the linker — protonated, as it is at neutral pH.

Aib is not one of the twenty standard residues, so its template is
downloaded from the RCSB on first use.
"""
    ),
    code(
        """
from mbuild.biopolymers import Protein

protein = Protein("semaglutide_apo.pdb", download=True)
print(protein.n_particles, "atoms, net formal charge", protein.net_formal_charge)

lysine = protein.get_residue(20, chain_id="A")
print("residue 2:", protein.get_residue(2, chain_id="A").name)
print(f"LYS 20: charge {lysine.formal_charge:+d}")
"""
    ),
    md(
        """
## The linker

The linker is a registered Chemical Component Dictionary component,
`KUT`. Building it from that component rather than from a SMILES string
matters for one reason: its atoms come out with the names the CCD gives
them. A fragment built from SMILES would get invented names, and every
tool downstream would then need to be told what those names mean.
"""
    ),
    code(
        """
from demo_fragments import fragment_from_ccd

linker = fragment_from_ccd("KUT", link_atom="C33")
print(linker.n_particles, "atoms, formal charge", linker.formal_charge)
print("bonds to the protein through:", linker.link_atoms["1"])
"""
    ),
    md(
        """
## Making the bond

Two calls. `deprotonate` takes the proton off the lysine side chain,
because the ammonium ion at neutral pH has no lone pair and does not
acylate; the neutral amine does. `attach` then removes one hydrogen
from each side and forms the bond.

`leaving_atom_names` says which hydrogen goes. The three hydrogens on
that nitrogen are chemically equivalent, so the choice is arbitrary as
chemistry — but the residue library downstream describes the product by
naming the atom that is *absent*, so the file has to agree with it.
"""
    ),
    code(
        """
protein.deprotonate(20, "NZ", chain_id="A")
record = protein.attach(
    linker,
    resnum=20,
    atom_name="NZ",
    chain_id="A",
    leaving_atom_names="HZ2",
    relax=False,
)

print(f"{record.residue1.name} {record.atom1_name} - "
      f"{record.residue2.name} {record.atom2_name}")
print("hydrogens removed:", record.leaving1 + record.leaving2)
print(protein.n_particles, "atoms, net formal charge", protein.net_formal_charge)
"""
    ),
    code(
        """
protein.save_pdb("semaglutide_mbuild.pdb", overwrite=True)
print(protein.bond_records()[-1])
"""
    ),
    md(
        """
## Reading it back

This is the workshop's loader call, copied unchanged. It declares the
crosslink: which two residues, which two atoms bond, and which atom
leaves each side. Those are the same names `attach` was given.

Nothing else is needed — no hand-written residue definition, no sidecar
file. `KUT` is a CCD component, so pablo already knows it.
"""
    ),
    code(
        """
from openff.pablo import STD_CCD_CACHE, topology_from_pdb

STD_CCD_CACHE.auto_download = True
library = STD_CCD_CACHE.with_crosslink(
    residues=["LYS", "KUT"],
    linking_atoms=["NZ", "C33"],
    leaving_atoms=[["HZ2"], ["H61"]],
    bond_order=1,
)

topology = topology_from_pdb("semaglutide_mbuild.pdb", residue_library=library)
molecule = topology.molecule(0)
print(topology.n_molecules, "molecule:", molecule.n_atoms, "atoms,",
      molecule.hill_formula + ",", "net charge", molecule.total_charge)
"""
    ),
    md(
        """
One molecule, not two. That is the whole claim: the linker came through
as part of the peptide, not as a separate thing sitting next to it.
"""
    ),
    code(
        """
view = topology.visualize()
view.clear_representations()
view.add_representation("cartoon", color="#990000")
view.add_representation("licorice", selection="[KUT]", color="#FF7733")
view.center(selection="[KUT]")
view
"""
    ),
    code(
        """
molecule.visualize("rdkit", show_all_hydrogens=False)
"""
    ),
    md(
        """
## Is it the right molecule?

`semaglutide_reference.pdb` is chain A of the workshop's own
`7KI0_prepared.pdb`. Reading both files through the same loader call
turns the question into a comparison of two OpenFF molecules.
"""
    ),
    code(
        """
import sys

sys.path.insert(0, "scripts")
from verify_semaglutide import check

check()
"""
    ),
    md(
        """
Set stereochemistry aside and the two are the same molecule — not
merely the same formula, but the same graph, atom for atom and bond for
bond.

The one stereocentre that does differ is worth a word, because it is the
kind of thing this pipeline exists to surface. mBuild placed the linker
at the geometry its CCD component defines; at that centre the deposited
coordinates disagree with the component. Neither structure is malformed.
They are two different molecules, and a loader that reads chemistry
rather than coordinates is what tells you so.

## Parameters and a short simulation

From here it is the ordinary OpenFF path, the same as the workshop's.
"""
    ),
    code(
        """
from openff.toolkit import ForceField

force_field = ForceField("openff_no_water-3.0.0-alpha0.offxml", "opc3.offxml")
interchange = force_field.create_interchange(topology)
print("parameterized", interchange.topology.n_atoms, "atoms")
"""
    ),
    code(
        """
import openmm
from openmm import unit

simulation = interchange.to_openmm_simulation(
    integrator=openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
    ),
)
simulation.minimizeEnergy(maxIterations=200)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.step(500)
state = simulation.context.getState(getEnergy=True)
print("potential energy:", state.getPotentialEnergy())
"""
    ),
    md(
        """
## What this needed

Reading the starting file demanded bond orders and formal charges that
are not in a PDB, so the loader matches CCD templates instead of
guessing. Writing the modified file demanded residue numbers, chain
identifiers and CONECT records for exactly the bonds a reader cannot
infer. And making the bond demanded control over which hydrogen leaves,
so that the file agrees with the residue library that reads it.

With those three, the handoff is one loader call with no special
pleading — and the structure that comes out is the one the workshop
simulates.
"""
    ),
]


def main():
    for name, cells in (
        ("01_load_and_export.ipynb", ONE),
        ("02_semaglutide.ipynb", TWO),
    ):
        path = HERE / name
        path.write_text(json.dumps(notebook(cells), indent=1) + "\n")
        print(f"wrote {name} ({len(cells)} cells)")


if __name__ == "__main__":
    main()
