"""Build the semaglutide assets the notebooks use.

Reads ``7KI0_prepared.pdb`` from the OpenFF workshop and writes two files:

``semaglutide_reference.pdb``
    Chain A as the workshop has it: the 31-residue peptide plus the KUT
    linker bonded to LYS 20, 580 atoms. This is the target the notebook
    compares against. It is the same structure openff-pablo ships as its
    own test fixture (``openff/pablo/_tests/data/semaglutide.pdb``),
    which labels the KUT residue chain B instead of chain A.

``semaglutide_apo.pdb``
    The same peptide with KUT removed and the lysine side chain put back
    into its neutral-pH form, 470 atoms. This is the notebook's input:
    an ordinary protonated protein, of the kind a preparation tool
    writes.

Run with ``pixi run assets``. Both outputs are committed, so a reader
never needs the 21 MB source file.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
SOURCE = HERE / "assets_cache" / "7KI0_prepared.pdb"
REFERENCE = HERE / "semaglutide_reference.pdb"
APO = HERE / "semaglutide_apo.pdb"

# N-H single bond length in Angstrom, and the tetrahedral angle.
NH_BOND = 1.01
TETRAHEDRAL = np.radians(109.5)


def unit(vector):
    return vector / np.linalg.norm(vector)


def position(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def atom_line(serial, name, resname, chain, resnum, pos, element):
    """Format one ATOM record in the fixed columns of the PDB standard."""
    field = f" {name:<3s}" if len(name) < 4 else name
    return (
        f"ATOM  {serial:5d} {field} {resname:>3s} {chain}{resnum:4d}    "
        f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00"
        f"          {element:>2s}\n"
    )


def split_chain_a(lines):
    """Return chain A's records, its TER lines, and its highest serial."""
    body, highest = [], 0
    for line in lines:
        record = line[:6]
        if record in ("ATOM  ", "HETATM") and line[21] == "A":
            body.append(line)
            highest = max(highest, int(line[6:11]))
        elif record == "TER   " and line[21:22] == "A":
            body.append(line)
    return body, highest


def conects_within(lines, highest):
    """Return the CONECT records whose serials all lie in chain A."""
    kept = []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        serials = [
            int(line[i : i + 5])
            for i in range(6, len(line.rstrip()), 5)
            if line[i : i + 5].strip()
        ]
        if serials and all(serial <= highest for serial in serials):
            kept.append(line)
    return kept


def complete_amine(nz, ce, hz1):
    """Return the two missing hydrogen positions on a lysine NZ.

    NZ keeps one hydrogen in the modified structure. Its other two sit
    where a tetrahedral nitrogen puts them, given the directions to CE
    and to the hydrogen that is present.

    Their exact placement does not matter to the notebook: the workflow
    removes one of them at ``deprotonate`` and the other at ``attach``,
    so neither survives into the file being compared.
    """
    to_ce, to_hz1 = unit(ce - nz), unit(hz1 - nz)
    bisector = -unit(to_ce + to_hz1)
    normal = unit(np.cross(to_ce, to_hz1))
    half = TETRAHEDRAL / 2
    return [
        nz + NH_BOND * unit(bisector * np.cos(half) + sign * normal * np.sin(half))
        for sign in (1, -1)
    ]


def build():
    if not SOURCE.exists():
        sys.exit(
            f"{SOURCE} is missing. Run `pixi run fetch` first; it downloads "
            "the 21 MB workshop file into a directory git ignores."
        )
    lines = SOURCE.read_text().splitlines(keepends=True)
    header = [line for line in lines[:5] if line[:6] in ("CRYST1", "REMARK")]
    body, highest = split_chain_a(lines)

    REFERENCE.write_text(
        "".join(header + body + conects_within(lines, highest)) + "END\n"
    )
    atoms = sum(1 for line in body if line[:6] in ("ATOM  ", "HETATM"))
    print(f"{REFERENCE.name}: {atoms} atoms")

    peptide = [
        line
        for line in body
        if not (line[:6] == "HETATM" and line[17:20] == "KUT")
        and not (line[:6] == "TER   " and line[17:20] == "KUT")
    ]

    def lysine_atom(name):
        return next(
            line
            for line in peptide
            if line[17:20] == "LYS"
            and int(line[22:26]) == 20
            and line[12:16].strip() == name
        )

    nz_line = lysine_atom("NZ")
    hz2, hz3 = complete_amine(
        position(nz_line), position(lysine_atom("CE")), position(lysine_atom("HZ1"))
    )

    # Insert the two hydrogens after HZ1 so the lysine's records stay
    # together, then renumber every serial from one.
    index = peptide.index(lysine_atom("HZ1")) + 1
    peptide[index:index] = [
        atom_line(0, "HZ2", "LYS", "A", 20, hz2, "H"),
        atom_line(0, "HZ3", "LYS", "A", 20, hz3, "H"),
    ]

    out, serial = header[:], 1
    for line in peptide:
        if line[:6] in ("ATOM  ", "HETATM"):
            out.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        elif line[:6] == "TER   ":
            out.append(f"TER   {serial:5d}      {line[17:27].rstrip()}\n")
            serial += 1
    APO.write_text("".join(out) + "END\n")
    print(f"{APO.name}: {serial - 2} atoms")


if __name__ == "__main__":
    build()
