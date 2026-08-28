"""Fork-only extras on top of mbuild.biopolymers.

These live outside the upstream mBuild contribution on purpose:
- fragment_from_pdb: ingestion of legacy fragment PDB files (e.g.
  GLYCAM glycans) with CONECT-only bonds and declared bond orders.
  Fragment construction is PolyzyMD territory.
- attach_multi: tether one fragment to several protein sites
  (rotate + shear placement, then protein-fixed relaxation).
- pablo_crosslink_kwargs: format a neutral Protein.bond_records()
  entry into openff-pablo's with_crosslink vocabulary.
"""

import logging

import numpy as np

from mbuild.biopolymers import Protein, Residue, prepare_fragment
from mbuild.biopolymers.protein import (
    InterResidueBond,
    _atom_in_residue,
    _parse_pdb,
)
from mbuild.compound import Compound
from mbuild.exceptions import MBuildError

logger = logging.getLogger(__name__)


def pablo_crosslink_kwargs(record):
    """Format one Protein.bond_records() entry for with_crosslink."""
    names = record["residue_names"]
    atoms = record["atom_names"]
    leaving = record["leaving_atoms"]
    if names[0] == names[1] and atoms[0] == atoms[1] and leaving[0] == leaving[1]:
        return {
            "residues": [names[0]],
            "linking_atoms": [atoms[0]],
            "leaving_atoms": [list(leaving[0])],
            "bond_order": record["bond_order"],
        }
    return {
        "residues": list(names),
        "linking_atoms": list(atoms),
        "leaving_atoms": [list(side) for side in leaving],
        "bond_order": record["bond_order"],
    }


def fragment_from_pdb(filename, bond_orders=None):
    """Load a fragment PDB file into Residue compounds for attach().

    Bonds come only from CONECT records; multiple bonds must be
    declared through ``bond_orders`` (keys: pairs of (resnum, atom
    name) tuples).
    """
    with open(filename) as handle:
        text = handle.read()
    groups, conects, _ = _parse_pdb(text)
    if not conects:
        raise MBuildError(
            f"{filename} has no CONECT records. fragment_from_pdb takes "
            "connectivity only from CONECT records."
        )
    orders = {frozenset(key): float(order) for key, order in (bond_orders or {}).items()}
    fragment = Compound(name="fragment")
    serial_to_particle = {}
    serial_key = {}
    for group in groups:
        residue = Residue(
            resname=group.resname, resnum=group.resnum,
            icode=group.icode, hetatm=True,
        )
        particles = []
        for record in group.records:
            if not record.element:
                raise MBuildError(
                    f"Atom {record.name!r} of {group.label} has no element column."
                )
            particle = Compound(
                name=record.name, element=record.element.capitalize(),
                pos=record.pos,
            )
            particles.append(particle)
            serial_to_particle[record.serial] = particle
            serial_key[record.serial] = (group.resnum, record.name)
        residue.add(particles)
        fragment.add(residue)
    used = set()
    for pair in conects:
        serials = tuple(pair)
        if len(serials) != 2:
            continue
        particles = [serial_to_particle.get(serial) for serial in serials]
        if None in particles:
            raise MBuildError(f"CONECT references unknown serial in {serials}.")
        if not fragment.bond_graph.has_edge(*particles):
            key = frozenset(serial_key[serial] for serial in serials)
            if key in orders:
                used.add(key)
            fragment.add_bond(particles, bond_order=orders.get(key, 1.0))
    unused = set(orders) - used
    if unused:
        raise MBuildError(
            f"bond_orders entries match no CONECT bond: "
            f"{sorted(tuple(sorted(key)) for key in unused)}."
        )
    return fragment


def _rotate_about(particles, pivot, from_vector, to_vector):
    """Rigidly rotate particles about a pivot, aligning two vectors."""
    norm_from = np.linalg.norm(from_vector)
    norm_to = np.linalg.norm(to_vector)
    if norm_from < 1e-8 or norm_to < 1e-8:
        return
    unit_from = from_vector / norm_from
    unit_to = to_vector / norm_to
    axis = np.cross(unit_from, unit_to)
    sine = np.linalg.norm(axis)
    cosine = float(np.dot(unit_from, unit_to))
    if sine < 1e-8:
        return
    axis = axis / sine
    skew = np.array(
        [[0.0, -axis[2], axis[1]],
         [axis[2], 0.0, -axis[0]],
         [-axis[1], axis[0], 0.0]]
    )
    rotation = np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)
    for particle in particles:
        particle.pos = pivot + rotation @ (particle.pos - pivot)


def _stretch_along(particles, pivot, link_atom, target):
    """Shear particles along pivot->link so the link atom hits target."""
    axis = link_atom.pos - pivot
    length = float(np.linalg.norm(axis))
    if length < 1e-8:
        return
    unit = axis / length
    displacement = target - link_atom.pos
    for particle in particles:
        weight = float(np.dot(particle.pos - pivot, unit)) / length
        particle.pos = particle.pos + np.clip(weight, 0.0, 1.0) * displacement


def attach_multi(protein, fragment, sites, fragment_resname=None, separation=0.15, relax=True):
    """Tether one fragment to several protein sites at once.

    sites maps each attachment-point label of the fragment ([*:1],
    [*:2] in its SMILES, or particle tags) to keyword arguments naming
    the protein site. The first site is placed by rigid port alignment;
    each further tether is placed by rotating the fragment about its
    first bond and shearing it onto the site, then a protein-fixed
    relaxation removes the strain.
    """
    if isinstance(fragment, str):
        fragment = prepare_fragment(fragment, fragment_resname or "LIG")
    probe_residues = (
        [fragment] if isinstance(fragment, Residue)
        else [c for c in fragment.successors() if isinstance(c, Residue)]
    )
    available = {
        label: (residue.resnum, residue.link_atoms[label])
        for residue in probe_residues
        for label in residue.link_atoms
    }
    missing = set(sites) - set(available)
    if missing:
        raise MBuildError(
            f"The fragment has no attachment points labeled {sorted(missing)}; "
            f"it has {sorted(available)}."
        )
    labels = list(sites)
    first_resnum, first_atom = available[labels[0]]
    records = [
        protein.attach(
            fragment, first_atom, fragment_resnum=first_resnum,
            fragment_resname=fragment_resname, separation=separation,
            relax=False, **sites[labels[0]],
        )
    ]
    for label in labels[1:]:
        candidates = [
            residue for residue in protein.residues()
            if label in residue.link_atoms and residue.hetatm
        ]
        if len(candidates) != 1:
            raise MBuildError(
                f"Attachment label {label!r} matches {len(candidates)} residues."
            )
        site = dict(sites[label])
        bond_order = int(site.pop("bond_order", 1))
        site_residue = protein.get_residue(
            site["resnum"], chain_id=site.get("chain_id"), icode=site.get("icode", "")
        )
        site_atom = protein.get_atom(
            site["resnum"], site["atom_name"],
            chain_id=site.get("chain_id"), icode=site.get("icode", ""),
        )
        frag_residue = candidates[0]
        frag_atom = _atom_in_residue(frag_residue, frag_residue.link_atoms[label])
        pivot = _atom_in_residue(records[0].residue2, records[0].atom2_name)
        fragment_particles = [
            particle
            for residue in protein.residues()
            if residue.hetatm and set(residue.link_atoms) & set(available)
            for particle in residue.particles()
        ]
        _rotate_about(
            fragment_particles, pivot.pos,
            frag_atom.pos - pivot.pos, site_atom.pos - pivot.pos,
        )
        approach = frag_atom.pos - site_atom.pos
        approach_norm = float(np.linalg.norm(approach))
        if approach_norm > 1e-8:
            target = site_atom.pos + separation * approach / approach_norm
        else:
            target = site_atom.pos + np.array([separation, 0.0, 0.0])
        _stretch_along(fragment_particles, pivot.pos, frag_atom, target)
        site_hydrogens = Protein._bonded_hydrogens(site_atom, site_residue.name, bond_order)
        frag_hydrogens = Protein._bonded_hydrogens(frag_atom, frag_residue.name, bond_order)
        for hydrogen in (*site_hydrogens, *frag_hydrogens):
            protein.remove(hydrogen)
        protein.add_bond((site_atom, frag_atom), bond_order=float(bond_order))
        distance = float(np.linalg.norm(site_atom.pos - frag_atom.pos))
        message = (
            f"Tether {label!r} formed at {distance * 10:.1f} A between "
            f"{site_residue.name} {site_residue.resnum} {site_atom.name} and "
            f"{frag_residue.name} {frag_residue.resnum} {frag_atom.name}."
        )
        if distance > 0.2:
            logger.warning(message + " The fragment cannot reach this site.")
        else:
            logger.info(message)
        record = InterResidueBond(
            residue1=site_residue, residue2=frag_residue,
            atom1_name=site_atom.name, atom2_name=frag_atom.name,
            order=bond_order,
            leaving1=tuple(sorted(h.name for h in site_hydrogens)),
            leaving2=tuple(sorted(h.name for h in frag_hydrogens)),
        )
        protein.cross_bonds.append(record)
        records.append(record)
    if relax:
        protein.relax_fragments(n_steps=2000)
    return records
