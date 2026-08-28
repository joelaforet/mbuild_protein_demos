"""Write a multi-state PDB movie of the double-tether placement.

Frames: initial rigid placement (chain dangling), 8 rotation steps,
8 shear steps, then protein-fixed minimization in chunks.
Play in PyMOL: load tether_movie.pdb, press play (set all_states off).
"""

import numpy as np

from mbuild.biopolymers import Protein, prepare_fragment
from mbuild.biopolymers.protein import _atom_in_residue

DEMO = "."
SMILES = "[*:1]CCOCCOCCOCCOCCOCCOCCOCC[*:2]"
SITE1 = dict(resnum=5, atom_name="NZ", chain_id="A")
SITE2 = dict(resnum=12, atom_name="NZ", chain_id="A")

frames = []


def capture(protein):
    """Snapshot the current coordinates as PDB atom lines."""
    path = f"{DEMO}/_frame.pdb"
    protein.save_pdb(path, overwrite=True)
    lines = open(path).read().splitlines()
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER"))]
    footer = [l for l in lines if l.startswith("CONECT")]
    frames.append(atoms)
    return footer


# --- placement: first tether by rigid alignment ------------------------
protein = Protein("6m03_protonated.pdb")
fragment = prepare_fragment(SMILES, "PEG")
record1 = protein.attach(fragment, fragment.link_atoms["1"],
                         relax=False, **SITE1)
peg = record1.residue2

# form the second tether topologically (hydrogens leave, bond forms) so
# the atom count is constant across all movie frames
site2_res = protein.get_residue(SITE2["resnum"], chain_id=SITE2["chain_id"])
site2_atom = protein.get_atom(**SITE2)
frag2_atom = _atom_in_residue(peg, peg.link_atoms["2"])
for hydrogen in (
    *Protein._bonded_hydrogens(site2_atom, site2_res.name, 1),
    *Protein._bonded_hydrogens(frag2_atom, peg.name, 1),
):
    protein.remove(hydrogen)
protein.add_bond((site2_atom, frag2_atom), bond_order=1.0)
footer = capture(protein)  # frame 1: dangling chain
print("frame 1 (placement): tether =",
      round(float(np.linalg.norm(site2_atom.pos - frag2_atom.pos)) * 10, 1), "A")

peg_particles = list(peg.particles())
pivot = _atom_in_residue(peg, peg.link_atoms["1"])

# --- rotation in 8 increments ------------------------------------------
v_from = frag2_atom.pos - pivot.pos
v_to = site2_atom.pos - pivot.pos
u_from = v_from / np.linalg.norm(v_from)
u_to = v_to / np.linalg.norm(v_to)
axis = np.cross(u_from, u_to)
axis /= np.linalg.norm(axis)
angle = float(np.arccos(np.clip(np.dot(u_from, u_to), -1, 1)))
skew = np.array([[0, -axis[2], axis[1]],
                 [axis[2], 0, -axis[0]],
                 [-axis[1], axis[0], 0]])
for step in range(8):
    theta = angle / 8
    rotation = (np.eye(3) + np.sin(theta) * skew
                + (1 - np.cos(theta)) * (skew @ skew))
    center = pivot.pos.copy()
    for particle in peg_particles:
        particle.pos = center + rotation @ (particle.pos - center)
    footer = capture(protein)
print("after rotation: tether =",
      round(float(np.linalg.norm(site2_atom.pos - frag2_atom.pos)) * 10, 1), "A")

# --- shear in 8 increments ----------------------------------------------
axis_vec = frag2_atom.pos - pivot.pos
length = float(np.linalg.norm(axis_vec))
unit = axis_vec / length
gap = site2_atom.pos - frag2_atom.pos
target_shift = gap - 0.15 * gap / np.linalg.norm(gap)
weights = {
    id(p): float(np.clip(np.dot(p.pos - pivot.pos, unit) / length, 0, 1))
    for p in peg_particles
}
for step in range(8):
    for particle in peg_particles:
        particle.pos = particle.pos + weights[id(particle)] * target_shift / 8
    footer = capture(protein)
print("after shear: tether =",
      round(float(np.linalg.norm(site2_atom.pos - frag2_atom.pos)) * 10, 1), "A")

# --- protein-fixed minimization, captured in chunks ---------------------
from mbuild.simulation import OpenMMSimulation

mobile = set(peg_particles)
simulation = OpenMMSimulation(protein, forcefield=None, kick=False)
for index, particle in enumerate(protein.particles()):
    if particle not in mobile:
        simulation.system.setParticleMass(index, 0.0)
for chunk in range(30):
    simulation.minimize(n_steps=25, tolerance=10.0)
    footer = capture(protein)
print("after minimization: tether =",
      round(float(np.linalg.norm(site2_atom.pos - frag2_atom.pos)) * 10, 1), "A")

# --- write the multi-state PDB ------------------------------------------
with open(f"{DEMO}/tether_movie.pdb", "w") as handle:
    for index, atoms in enumerate(frames, start=1):
        handle.write(f"MODEL     {index:4d}\n")
        handle.write("\n".join(atoms) + "\n")
        handle.write("ENDMDL\n")
    handle.write("\n".join(footer) + "\nEND\n")
print(f"wrote tether_movie.pdb with {len(frames)} states")
