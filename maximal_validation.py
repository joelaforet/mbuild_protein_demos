"""Validation for the maximal mixed-conjugation demo.

Builds, on ubiquitin (1UBQ):
1. LYS63 NZ: an SBM-NHS-EGP methacrylate trimer, amide-bonded through
   the middle (NHS-derived) monomer. Multi-residue linear polymer.
2. ASN60 ND2: a chitobiose N-glycan loaded from a GLYCAM PDB with
   CONECT records. Multi-residue glycan from PDB.
3. SER20 OG: a three-branch N-glycan-like sugar from SMILES as one
   residue. Branched glycan from SMILES.

Then writes the prepared PDB and loads it through openff-pablo.
"""

import numpy as np
import mbuild as mb
from mbuild.biopolymers import Protein, prepare_fragment
from mbuild_extras import fragment_from_pdb, pablo_crosslink_kwargs

DEMO = "."
GLYCAN_PDB = "glycam_G42666HT_CONECT.pdb"

# Product-form monomer units: backbone written as CH3-C(CH3)(R)-H, so one
# hydrogen on each backbone carbon can leave when the chain links.
SBM_SMILES = "CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
EGP_SMILES = "CC(C)C(=O)OCCOc1ccccc1"
NHS_SMILES = "CC(C)C=O"  # the aldehyde H marks the future amide bond

THREEBRANCH_SMILES = (
    "CC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[CH2]1)CO)"  # anomeric OH -> CH2
    "O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)"
    "O[C@H]3[C@H]([C@H]([C@@H]([C@H](O3)CO[C@@H]4[C@H]([C@H]([C@@H]([C@H](O4)"
    "CO[C@@H]5[C@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O)O)"
    "O[C@@H]6[C@H]([C@H]([C@@H]([C@H](O6)CO)O)O)"
    "O[C@@H]7[C@H]([C@H]([C@@H]([C@H](O7)CO)O)O)O)O)O)"
    "O[C@@H]8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)"
    "O[C@@H]9[C@H]([C@H]([C@@H]([C@H](O9)CO)O)O)"
    "O[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O)O)NC(=O)C)O"
)


def prepare_monomer(smiles, resname):
    """Prepare a monomer unit with fixed backbone atom names CBH/CBT."""
    residue = prepare_fragment(mb.load(smiles, smiles=True), resname)
    particles = list(residue.particles())
    particles[0].name = "CBH"  # head backbone carbon (CH3 in the free unit)
    particles[1].name = "CBT"  # tail backbone carbon (carries R and one H)
    return residue


def carbonyl_carbon(residue):
    """Return the name of the carbon double-bonded to an oxygen."""
    root = residue.root if residue.parent is not None else residue
    for p1, p2, data in root.bonds(return_bond_order=True):
        if data["bond_order"] == 2.0 and {p1.element.symbol, p2.element.symbol} == {"C", "O"}:
            return (p1 if p1.element.symbol == "C" else p2).name
    raise ValueError("no carbonyl found")


def anomeric_ch2(residue):
    """Return the ring CH2 bonded to the ring oxygen (the attachment C)."""
    for particle in residue.particles():
        if particle.element.symbol != "C":
            continue
        neighbors = list(particle.direct_bonds())
        hydrogens = [n for n in neighbors if n.element.symbol == "H"]
        oxygens = [n for n in neighbors if n.element.symbol == "O"]
        carbons = [n for n in neighbors if n.element.symbol == "C"]
        if len(hydrogens) == 2 and len(oxygens) == 1 and len(carbons) == 1:
            ring_oxygen = oxygens[0]
            heavy = [
                n for n in ring_oxygen.direct_bonds() if n.element.symbol != "H"
            ]
            if len(heavy) == 2:
                return particle.name
    raise ValueError("no anomeric CH2 found")


# ----------------------------------------------------------------------
# mBuild side
# ----------------------------------------------------------------------
protein = Protein(f"{DEMO}/1ubq_protonated.pdb")
print("loaded:", protein.n_particles, "atoms, net", protein.net_formal_charge)

# --- 1. trimer at LYS63, reactive monomer in the middle ---------------
sbm_unit = prepare_monomer(SBM_SMILES, "SBM")
egp_unit = prepare_monomer(EGP_SMILES, "EGP")
nhs_unit = prepare_monomer(NHS_SMILES, "NHS")
nhs_site = carbonyl_carbon(nhs_unit)
print("NHS carbonyl atom:", nhs_site)

rec_lys = protein.attach(nhs_unit, nhs_site, resnum=63, atom_name="NZ",
                         chain_id="A", fragment_resname="NHS")
# Amide nitrogen is neutral: remove a second proton from NZ.
lys63 = protein.get_residue(63, chain_id="A")
hz3 = [p for p in lys63.particles() if p.name == "HZ3"][0]
protein.remove(hz3)
lys63.formal_charge -= 1

nhs_res = rec_lys.residue2
rec_sbm = protein.attach(sbm_unit, "CBT", resnum=nhs_res.resnum,
                         atom_name="CBH", chain_id="A")
rec_egp = protein.attach(egp_unit, "CBH", resnum=nhs_res.resnum,
                         atom_name="CBT", chain_id="A")
# Backbone order in the file must be SBM-NHS-EGP for Pablo adjacency.
rec_sbm.residue2.resnum, nhs_res.resnum, rec_egp.residue2.resnum = 77, 78, 79
print("trimer records:", [(r.atom1_name, r.atom2_name, r.leaving1, r.leaving2)
                          for r in (rec_lys, rec_sbm, rec_egp)])

# --- 2. chitobiose from GLYCAM PDB at ASN60 ---------------------------
glycan = fragment_from_pdb(
    GLYCAN_PDB,
    bond_orders={
        ((2, "C2N"), (2, "O2N")): 2,  # N-acetyl C=O of 4YB
        ((3, "C2N"), (3, "O2N")): 2,  # N-acetyl C=O of 0YB
    },
)
resnames = [r.name for r in glycan.children]
print("glycan residues:", resnames)
# The GLYCAM ROH cap is the anomeric OH; it leaves on conjugation.
roh = [r for r in glycan.children if r.name == "ROH"][0]
fyb = [r for r in glycan.children if r.name == "4YB"][0]
c1 = [p for p in fyb.particles() if p.name == "C1"][0]
o1 = [p for p in roh.particles() if p.name == "O1"][0]
direction = (o1.pos - c1.pos) / np.linalg.norm(o1.pos - c1.pos)
glycan.remove(list(roh.particles()))
placeholder = mb.Compound(name="H0", element="H", pos=c1.pos + 0.109 * direction)
fyb.add(placeholder)
glycan.add_bond((c1, placeholder), bond_order=1.0)
# Rename the second sugar's anomeric carbon so the inter-sugar linking
# bond (O4 -> C1A) cannot be confused with 4YB's crosslink atom C1.
oyb = [r for r in glycan.children if r.name == "0YB"][0]
[p for p in oyb.particles() if p.name == "C1"][0].name = "C1A"

rec_glycan = protein.attach(glycan, "C1", resnum=60, atom_name="ND2",
                            chain_id="A", fragment_resnum=2)
print("glycan record:", rec_glycan.atom1_name, rec_glycan.atom2_name,
      rec_glycan.leaving1, rec_glycan.leaving2)

# --- 3. three-branch glycan from SMILES at SER20 ----------------------
ng3 = prepare_fragment(mb.load(THREEBRANCH_SMILES, smiles=True), "NG3")
ng3_site = anomeric_ch2(ng3)
rec_ng3 = protein.attach(ng3, ng3_site, resnum=20, atom_name="OG",
                         chain_id="A")
print("NG3 site:", ng3_site, "record:", rec_ng3.leaving1, rec_ng3.leaving2)

# --- export ------------------------------------------------------------
protein.save_pdb(f"{DEMO}/1ubq_modified.pdb", overwrite=True)
records = protein.bond_records()
crosslink_only = [
    pablo_crosslink_kwargs(record)
    for record in records
    if set(record["atom_names"]) != {"CBH", "CBT"}  # backbone rides linking
]
print("specs for with_crosslink:", crosslink_only)
print("net formal charge:", protein.net_formal_charge)


# ----------------------------------------------------------------------
# Pablo side (demo-notebook glue; no OpenFF code in mBuild)
# ----------------------------------------------------------------------
from rdkit import Chem
from openff.toolkit import Molecule
from openff.pablo import STD_CCD_CACHE, ResidueDefinition, topology_from_pdb
from openff.pablo.residue import BondDefinition

POLYMER_LINK = BondDefinition.with_defaults("CBT", "CBH")
GLYCO_LINK = BondDefinition.with_defaults("O4", "C1A")


def named_offmol_from_smiles(smiles, prepared_residue):
    """OpenFF molecule from SMILES with mBuild's final atom names."""
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    offmol = Molecule.from_rdkit(
        rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True
    )
    particles = list(prepared_residue.particles())
    assert [a.symbol for a in offmol.atoms] == [
        p.element.symbol for p in particles
    ], "atom order mismatch between SMILES and mBuild fragment"
    for atom, particle in zip(offmol.atoms, particles):
        atom.name = particle.name
    return offmol


def monomer_definition(smiles, prepared_residue, resname, leaving_names):
    offmol = named_offmol_from_smiles(smiles, prepared_residue)
    for atom in offmol.atoms:
        if atom.name in leaving_names:
            atom.metadata["leaving_atom"] = True
    return ResidueDefinition.from_molecule(
        offmol, residue_name=resname, linking_bond=POLYMER_LINK
    )


def sugar_definition(residue, resname, extra_hydrogens, linking_bond):
    """Definition from an mBuild fragment residue graph plus leaving Hs.

    Sugars from the GLYCAM PDB are neutral; bond orders come from the
    fragment (single, plus the declared C=O bonds).
    """
    particles = list(residue.particles())
    editable = Chem.RWMol()
    index = {}
    for particle in particles:
        atom = Chem.Atom(particle.element.symbol)
        atom.SetNoImplicit(True)
        index[particle] = editable.AddAtom(atom)
    root = residue.root
    for p1, p2, data in root.bonds(return_bond_order=True):
        if p1 in index and p2 in index:
            bond_type = (
                Chem.BondType.DOUBLE
                if data["bond_order"] == 2.0
                else Chem.BondType.SINGLE
            )
            editable.AddBond(index[p1], index[p2], bond_type)
    names = [particle.name for particle in particles]
    by_name = {particle.name: index[particle] for particle in particles}
    for host_name, hydrogen_name in extra_hydrogens:
        atom = Chem.Atom("H")
        atom.SetNoImplicit(True)
        new_index = editable.AddAtom(atom)
        editable.AddBond(by_name[host_name], new_index, Chem.BondType.SINGLE)
        names.append(hydrogen_name)
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)
    offmol = Molecule.from_rdkit(
        mol, allow_undefined_stereo=True, hydrogens_are_explicit=True
    )
    extra_names = {hydrogen for _, hydrogen in extra_hydrogens}
    for atom, name in zip(offmol.atoms, names):
        atom.name = name
        if name in extra_names:
            atom.metadata["leaving_atom"] = True
    return ResidueDefinition.from_molecule(
        offmol, residue_name=resname, linking_bond=linking_bond
    )


sbm_def = monomer_definition(SBM_SMILES, sbm_unit, "SBM", set(rec_sbm.leaving2))
egp_def = monomer_definition(EGP_SMILES, egp_unit, "EGP", set(rec_egp.leaving2))
nhs_def = monomer_definition(
    NHS_SMILES, nhs_unit, "NHS",
    set(rec_sbm.leaving1) | set(rec_egp.leaving1),
)
fyb_def = sugar_definition(fyb, "4YB", [("O4", "HO4")], GLYCO_LINK)
oyb_def = sugar_definition(oyb, "0YB", [("C1A", "HL")], GLYCO_LINK)
ng3_def_mol = named_offmol_from_smiles(THREEBRANCH_SMILES, ng3)
ng3_def = ResidueDefinition.from_molecule(ng3_def_mol, residue_name="NG3")

library = STD_CCD_CACHE.with_(
    {
        "SBM": [sbm_def],
        "NHS": [nhs_def],
        "EGP": [egp_def],
        "4YB": [fyb_def],
        "0YB": [oyb_def],
        "NG3": [ng3_def],
    }
)
for spec in crosslink_only:
    library = library.with_crosslink(**spec)

top = topology_from_pdb(f"{DEMO}/1ubq_modified.pdb", residue_library=library)
mol = top.molecule(0)
print("PABLO OK — molecules:", top.n_molecules, "atoms:", mol.n_atoms,
      "net charge:", mol.total_charge)
nz = [a for a in mol.atoms
      if a.name == "NZ" and a.metadata.get("residue_number") == 63][0]
print("LYS63 NZ charge:", nz.formal_charge, "bonded:",
      sorted(n.name for n in nz.bonded_atoms))
sbm_charges = {a.name: int(a.formal_charge.m)
               for a in mol.atoms
               if a.metadata.get("residue_name") == "SBM" and a.formal_charge.m}
print("SBM charged atoms:", sbm_charges)
