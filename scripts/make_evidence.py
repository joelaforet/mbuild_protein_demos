"""Generate the evidence notebooks in ``evidence/``.

Each notebook backs one claim made in the pull request descriptions of
the mBuild biopolymers stack. They are deliberately short: a one-line
claim, then the code that proves it, with every line visible. Run with
``pixi run evidence`` to regenerate and execute them.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
OUT = HERE / "evidence"

CORPUS_URL = (
    "https://raw.githubusercontent.com/openforcefield/openff-pablo/main/"
    "openff/pablo/_tests/data/prepared_pdbs/"
)
CORPUS_FILES = [
    "193l_prepared.pdb",
    "1FLR_prepared.pdb",
    "1a4t_samechain.pdb",
    "1csa_maestro.pdb",
    "1csa_maestro_waterfirst.pdb",
    "1hje_diffchain.pdb",
    "1hje_samechain.pdb",
    "1p3q_noter.pdb",
    "2MUM_blowup.pdb",
    "2MUM_composed_function.pdb",
    "2MUM_discontiguous_resseq.pdb",
    "2MUM_discontiguous_serial.pdb",
    "2MUM_dryrun.pdb",
    "2MUM_icode.pdb",
    "2MUM_letters_in_resseq.pdb",
    "2MUM_letters_in_serial.pdb",
    "2MUM_neutralized.pdb",
    "2MUM_reuse_resseq.pdb",
    "2MUM_reuse_serial.pdb",
    "2hi7_prepared.pdb",
    "2zuq_prepared.pdb",
    "3h34_prepared.pdb",
    "3ip9_dye_solvated.pdb",
    "5eil_fixed.pdb",
    "ions.pdb",
    "polyglycines.pdb",
]


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


QUIET = """
import logging, warnings
from rdkit import RDLogger

warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)
RDLogger.DisableLog("rdApp.*")
"""

DOWNLOAD = f"""
# Pablo's own test corpus of prepared PDB files, fetched into a git-ignored cache.
import urllib.request
from pathlib import Path

CORPUS_URL = "{CORPUS_URL}"
CORPUS_FILES = {json.dumps(CORPUS_FILES, indent=4)}
cache = Path("../assets_cache/pablo_prepared_pdbs")
cache.mkdir(parents=True, exist_ok=True)
for name in CORPUS_FILES:
    if not (cache / name).exists():
        urllib.request.urlretrieve(CORPUS_URL + name, cache / name)
len(list(cache.glob("*.pdb"))), "files"
"""


LOADER = [
    md(
        """
# The existing PDB loader carries no bond orders and no formal charges

`mb.load` reads a PDB file through mdtraj. `Protein` reads the same file
through the residue definitions. Same atoms, same bonds, different chemistry.
"""
    ),
    code(QUIET),
    code(
        """
from collections import Counter

import mbuild as mb
from mbuild.biopolymers import Protein

generic = mb.load("../semaglutide_apo.pdb")
protein = Protein("../semaglutide_apo.pdb")


def bond_orders(compound):
    return dict(Counter(d["bond_order"] for *_, d in compound.bonds(return_bond_order=True)))


print("mb.load :", generic.n_particles, "atoms", generic.n_bonds, "bonds", bond_orders(generic))
print("Protein :", protein.n_particles, "atoms", protein.n_bonds, "bonds", bond_orders(protein))
"""
    ),
    code(
        """
print("mb.load particle charges:", dict(Counter(p.charge for p in generic.particles())))
print("Protein net formal charge:", protein.net_formal_charge)
print({f"{r.name}{r.resnum}": r.formal_charge for r in protein.residues() if r.formal_charge})
"""
    ),
    md("Only one of the two can become an OpenFF `Molecule`."),
    code(
        """
from openff.toolkit import Molecule

molecule = Molecule.from_rdkit(protein.to_rdkit(), allow_undefined_stereo=True)
print("from Protein:", molecule.n_atoms, "atoms, net charge", molecule.total_charge)

try:
    Molecule.from_rdkit(generic.to_rdkit(), allow_undefined_stereo=True)
except Exception as error:
    print("from mb.load:", type(error).__name__, str(error).splitlines()[0])
"""
    ),
]


CORPUS = [
    md(
        """
# Loading openff-pablo's prepared PDB corpus

Every file in Pablo's `prepared_pdbs` test set, loaded with `Protein`.
Templates outside the 34 shipped components download from the RCSB.
"""
    ),
    code(QUIET),
    code(DOWNLOAD),
    code(
        """
from mbuild.biopolymers import Protein

loaded, failed = [], []
for path in sorted(cache.glob("*.pdb")):
    try:
        protein = Protein(path, download=True)
    except Exception as error:
        failed.append((path.name, str(error).splitlines()[0]))
        continue
    loaded.append(path.name)
    print(
        f"{path.name:32s} {len(protein.chains):3d} chains"
        f" {len(list(protein.residues())):5d} residues"
        f" {protein.n_particles:6d} atoms"
        f"  net charge {protein.net_formal_charge:+d}"
        f"  crosslinks {len(protein.bond_records())}"
    )
print()
print(len(loaded), "loaded,", len(failed), "refused")
"""
    ),
    md("Each refusal names the residue and says why."),
    code(
        """
for name, reason in failed:
    print(f"{name:32s} {reason}")
"""
    ),
    md(
        """
Pablo's verdict on the same eight files, for comparison. Three of them Pablo
refuses too without extra input. One is DNA, which `Protein` does not claim.
The other four are two copies each of a cyclic peptide (`1csa`) and of a
peptide with a C-terminal `NH2` cap (`1hje`). Pablo handles both and
`Protein` does not yet.
"""
    ),
    code(
        """
from openff.pablo import STD_CCD_CACHE, topology_from_pdb

STD_CCD_CACHE.auto_download = True
for name, _ in failed:
    try:
        topology = topology_from_pdb(cache / name)
        print(f"{name:32s} pablo loads it, {topology.n_atoms} atoms")
    except Exception as error:
        print(f"{name:32s} pablo refuses it too: {type(error).__name__}")
"""
    ),
]


ROUND_TRIP = [
    md(
        """
# Write a protein back out, and read it with openff-pablo unchanged

Hen lysozyme from Pablo's corpus, which has four disulfides. `save_pdb`
writes a `CONECT` for each disulfide and for nothing else, and Pablo reads
the file with no extra arguments.
"""
    ),
    code(QUIET),
    code(
        """
import urllib.request
from pathlib import Path

cache = Path("../assets_cache")
cache.mkdir(exist_ok=True)
source = cache / "193l_prepared.pdb"
if not source.exists():
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/openforcefield/openff-pablo/main/"
        "openff/pablo/_tests/data/prepared_pdbs/193l_prepared.pdb",
        source,
    )
"""
    ),
    code(
        """
from mbuild.biopolymers import Protein

protein = Protein(source)
print(len(list(protein.residues())), "residues,", protein.n_particles, "atoms,", protein.n_bonds, "bonds, net charge", protein.net_formal_charge)
for record in protein.bond_records():
    print(record["residue_names"], record["residue_numbers"], record["atom_names"], "leaving", record["leaving_atoms"])
"""
    ),
    code(
        """
written = cache / "193l_mbuild.pdb"
protein.save_pdb(written, overwrite=True)

lines = written.read_text().splitlines()
print({kind: sum(line.startswith(kind) for line in lines) for kind in ("ATOM", "HETATM", "TER", "CONECT")})
print(*[line for line in lines if line.startswith("CONECT")], sep="\\n")
"""
    ),
    code(
        """
from openff.pablo import topology_from_pdb

topology = topology_from_pdb(written)
print(topology.n_molecules, "molecule,", topology.n_atoms, "atoms,", topology.n_bonds, "bonds, net charge", topology.molecule(0).total_charge)
assert topology.n_bonds == protein.n_bonds
"""
    ),
]


MULTI = [
    md(
        """
# Attaching a fragment made of many residues

The fragment is the 31-residue semaglutide peptide, taken as a `Chain`
whose children are bonded `Residue` objects. It is attached through the
alpha carbon of its last glycine to lysine 48 of ubiquitin. Every fragment
residue keeps its name, its atoms and its formal charge.
"""
    ),
    code(QUIET),
    code(
        """
from mbuild.biopolymers import Protein

fragment = Protein("../semaglutide_apo.pdb").chains[0]
print(len(fragment.children), "residues:", [(r.name, r.resnum) for r in fragment.children][:5], "...")

protein = Protein("../1ubq_protonated.pdb")
protein.deprotonate(48, "NZ")
bond = protein.attach(fragment, fragment_atom_name="CA", fragment_resnum=31, resnum=48, atom_name="NZ", relax=False)
"""
    ),
    code(
        """
residues = list(protein.residues())
print(len(residues), "residues,", protein.n_particles, "atoms, net charge", protein.net_formal_charge)
print("ubiquitin ends, fragment begins:", [(r.name, r.resnum, r.formal_charge) for r in residues[74:80]])
print("fragment ends:", [(r.name, r.resnum, r.formal_charge) for r in residues[-3:]])
"""
    ),
    md("One bond record describes the modification. It is what a downstream residue library needs."),
    code(
        """
record, = protein.bond_records()
record
"""
    ),
    code(
        """
from pathlib import Path

written = Path("../assets_cache/1ubq_plus_peptide.pdb")
protein.save_pdb(written, overwrite=True)
lines = written.read_text().splitlines()
print(sorted({(line[17:20], int(line[22:26])) for line in lines if line.startswith(("ATOM", "HETATM"))}, key=lambda t: t[1])[74:82])
"""
    ),
    md("Pablo reads the file once the bond record is handed over as a crosslink."),
    code(
        """
from openff.pablo import STD_CCD_CACHE, topology_from_pdb

library = STD_CCD_CACHE.with_crosslink(
    residues=list(record["residue_names"]),
    linking_atoms=list(record["atom_names"]),
    leaving_atoms=[list(side) for side in record["leaving_atoms"]],
    bond_order=record["bond_order"],
)
topology = topology_from_pdb(written, residue_library=library)
print(topology.n_molecules, "molecule,", topology.n_atoms, "atoms,", topology.n_bonds, "bonds, net charge", topology.molecule(0).total_charge)
assert topology.n_bonds == protein.n_bonds
"""
    ),
]


DEPROTONATE = [
    md(
        """
# Why `deprotonate` comes before `attach`

`attach` replaces a hydrogen with a bond, so the anchor atom keeps its formal
charge. A methyl group is attached to lysine 48 of ubiquitin, with and
without deprotonating the amine first.
"""
    ),
    code(QUIET),
    code(
        """
from mbuild.biopolymers import Protein, fragment_from_smiles
from openff.toolkit import Molecule

methyl = fragment_from_smiles("*C", "MEE")


def report(label, protein):
    molecule = Molecule.from_rdkit(protein.to_rdkit(), allow_undefined_stereo=True)
    print(f"{label:28s} LYS48 {protein.get_residue(48).formal_charge:+d}   net {protein.net_formal_charge:+d}   exported {molecule.total_charge}")
"""
    ),
    code(
        """
protein = Protein("../1ubq_protonated.pdb")
report("as loaded", protein)

protein.attach(methyl, resnum=48, atom_name="NZ", relax=False)
report("attach only", protein)
"""
    ),
    code(
        """
protein = Protein("../1ubq_protonated.pdb")
protein.deprotonate(48, "NZ")
report("deprotonate", protein)

protein.attach(methyl, resnum=48, atom_name="NZ", relax=False)
report("deprotonate, then attach", protein)
"""
    ),
]


NOTEBOOKS = {
    "loader_bond_orders.ipynb": LOADER,
    "pablo_corpus.ipynb": CORPUS,
    "pdb_round_trip.ipynb": ROUND_TRIP,
    "multi_residue_fragment.ipynb": MULTI,
    "deprotonate_then_attach.ipynb": DEPROTONATE,
}


def main():
    OUT.mkdir(exist_ok=True)
    for name, cells in NOTEBOOKS.items():
        path = OUT / name
        path.write_text(json.dumps(notebook(cells), indent=1) + "\n")
        print("wrote", path.relative_to(HERE))


if __name__ == "__main__":
    main()
