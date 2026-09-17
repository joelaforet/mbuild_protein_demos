# mBuild protein demos

**Build a covalently modified protein with
[mBuild](https://github.com/mosdef-hub/mbuild), and read it straight into
OpenFF with [Pablo](https://github.com/openforcefield/openff-pablo).**

Notebook 02 builds semaglutide's modified chain A — the structure the
OpenFF post-translational-modification workshop simulates — and reads it
back with the workshop's own loader call, unchanged. No hand-written
residue definition, no companion file.

```
pixi install && pixi run setup
pixi run prove
```

```
  PASS  one molecule, so the linker is bonded
  PASS  same atom count
  PASS  same formula
  PASS  same net charge
  PASS  same atoms, by element and formal charge
  PASS  same bonds, by elements and order
  PASS  isomorphic, setting stereochemistry aside

  580 atoms, 585 bonds, C187H289N45O59, net charge -2.0 elementary_charge
```

## The notebooks

| | |
| --- | --- |
| `01_load_and_export.ipynb` | Read a protein with its chemistry intact, and hand it to OpenFF. Why template matching rather than distance-based bond perception. |
| `02_semaglutide.ipynb` | The demo. Attach the lipid linker to a lysine, write the PDB, read it with Pablo, check it against the workshop's structure, parameterize and simulate. |
| `03_extensions.ipynb` | Fragments the CCD does not define, relaxing a fragment that lands in a clash, splitting partial charges across the modification site, and a 106-atom FRET dye. |

## The evidence notebooks

`evidence/` holds five short notebooks, one per claim made in the pull
request descriptions. Each states the claim in a line and proves it in
the code below, with nothing hidden in a helper.

| | |
| --- | --- |
| `loader_bond_orders.ipynb` | `mb.load` gives every bond order 0.0 and no atom a charge. `Protein` gives 53 double bonds and a net charge of -1 on the same file. |
| `pablo_corpus.ipynb` | 18 of the 26 files in openff-pablo's prepared PDB corpus load. Each refusal names the residue, and Pablo's verdict on the same files is shown. |
| `pdb_round_trip.ipynb` | Lysozyme written by `save_pdb` has 8 `CONECT` lines, one per disulfide end, and Pablo reads it with no arguments. |
| `multi_residue_fragment.ipynb` | A 31-residue peptide attached to a lysine keeps every residue's name, atoms and charge, and one bond record is enough for Pablo to read the product. |
| `deprotonate_then_attach.ipynb` | Attaching to a charged lysine keeps the +1. Deprotonating first gives the neutral product. |

## What needed new code

A PDB file carries elements and coordinates. Assigning force field
parameters needs bond orders and formal charges as well, and those are
not in the file. Three things follow.

**Reading.** Every residue is matched by atom name against a wwPDB
Chemical Component Dictionary template, and the chemistry comes from that
template. Nothing is inferred from interatomic distances. A residue no
template explains raises an error naming the residue, instead of
producing a structure that is quietly wrong.

**Writing.** `save_pdb` emits real residue numbers, chain identifiers, a
TER after each chain, and CONECT records for exactly the bonds that
residue adjacency cannot imply. Peptide bonds stay implied, because a
residue-template reader rejects a CONECT its own definitions cannot
explain.

**Bonding.** `attach` takes `leaving_atom_names`, so the caller chooses
which hydrogen leaves. The hydrogens on one nitrogen are chemically
equivalent, so the choice is arbitrary as chemistry — but a residue
library describes the product by naming the atom that is *absent*, and
the file has to agree with it.

## The structures

`semaglutide_apo.pdb` (470 atoms) is the notebook's input: the
31-residue peptide with an ordinary protonated lysine at position 20.
`semaglutide_reference.pdb` (580 atoms) is the answer: chain A of the
workshop's `7KI0_prepared.pdb`, peptide plus linker. Both are committed.

`scripts/build_assets.py` derives them from the workshop file, which
`pixi run fetch` downloads and git ignores. The reference is the same
structure openff-pablo ships as its own semaglutide test fixture.

One stereocentre in the linker differs between the two. mBuild places
the linker at the geometry its CCD component defines, and at that centre
the deposited coordinates disagree with the component. Notebook 02 says
so rather than hiding it: a loader that reads chemistry rather than
coordinates is what makes such a disagreement visible at all.

## Tasks

| | |
| --- | --- |
| `pixi run setup` | Install mBuild and Pablo, and cache the KUT and AIB residue templates. Needs network once. |
| | Installs from `feat/biopolymers-docs`, the tip of the review stack on the fork, which contains every layer. Swap the branch for a commit SHA to pin a talk to an exact build. |
| `pixi run prove` | Build the structure and check it against the reference. No notebook. |
| `pixi run lab` | Open the notebooks. |
| `pixi run verify` | Execute all three notebooks headless. |
| `pixi run dev` | Point the environment at a local mBuild checkout. |
| `pixi run fetch` / `assets` | Re-download the workshop file and rebuild the committed structures. |
| `pixi run notebooks` | Regenerate notebooks 01 and 02 from `scripts/make_notebooks.py`. |
| `pixi run evidence` | Regenerate and execute the evidence notebooks from `scripts/make_evidence.py`. Downloads Pablo's 9 MB corpus once. |

## Force fields

Notebook 02 uses the OpenFF Rosemary alpha
(`openff_no_water-3.0.0-alpha0.offxml`), which covers the protein and the
modification with one model, so no charge surgery is needed.

Notebook 03 keeps the split-charge treatment for the case where that is
not true: Amber ff14SB library charges on the unmodified residues, NAGL
AM1-BCC graph charges on the fragment and the residue it is attached to,
and the small residual spread over the atoms of that site. The two sets
come from different fits, so atoms across the seam are not mutually
polarized. See `docs/charge-splitting.md`.

mBuild owns the structure. The force field assignment is downstream on
purpose: MosDef ships no biopolymer force field, so the OpenFF ecosystem
takes that step.
