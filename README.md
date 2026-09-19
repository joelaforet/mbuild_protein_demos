# mBuild protein demos

**Build a covalently modified protein with
[mBuild](https://github.com/mosdef-hub/mbuild), and read it straight into
OpenFF with [Pablo](https://github.com/openforcefield/openff-pablo).**

Notebook 01 loads a protein with its chemistry intact, notebook 02 puts a
fragment written as SMILES onto it and takes the result through OpenFF to a
short simulation, and notebook 05 rebuilds semaglutide's modified chain A,
the structure the OpenFF post-translational-modification workshop
simulates, and reads it back with the workshop's own loader call,
unchanged. Notebooks 03 and 04 mutate a residue's side chain and put two
dyes on the mutant by reactions given as a string.

```
pixi install && pixi run setup
pixi run lab
```

Open notebook 01 and go in order. Every check a notebook makes is written
out in the notebook itself.

## The notebooks

| | |
| --- | --- |
| `01_load_and_export.ipynb` | Read a protein with its chemistry intact, and hand it to OpenFF. Why template matching rather than distance-based bond perception. |
| `02_modify_a_protein.ipynb` | The `attach` tutorial. A fragment from star-marked SMILES, deprotonate the site, attach, relax the clash as a movie, write the PDB and bond records, build the Pablo definition and crosslink in the open, split the partial charges, simulate. |
| `03_point_mutations.ipynb` | `mutate`: the side chains a residue can take (the 20 canonical residues and any peptide-linking CCD component, drawn with RDKit), the fibronectin S1381AzF/S1500C construct from PDB 1FNF, the same mutation from a SMILES side chain, L and D side by side, export and a Pablo round trip. |
| `04_label_with_reactions.ipynb` | Reaction strings on `attach`: a DBCO donor clicked onto the azide, a maleimide acceptor added to the cysteine, the click product merged into one residue, and Pablo residue definitions for the three new residues built from the mBuild residues with Pablo's public API. |
| `05_semaglutide.ipynb` | The case study. Attach the lipid linker from its CCD component to a lysine, write the PDB, read it with Pablo with no hand-written definition, check it against the workshop's structure, parameterize and simulate. |

## The evidence notebooks

`evidence/` holds six short notebooks, one per claim made in the pull
request descriptions. Each states the claim in a line and proves it in
the code below, with nothing hidden in a helper. Run one to see its
output; the committed copies carry none.

| | |
| --- | --- |
| `loader_bond_orders.ipynb` | `mb.load` gives every bond order 0.0 and no atom a charge. `Protein` gives 53 double bonds and a net charge of -1 on the same file. |
| `pablo_corpus.ipynb` | 18 of the 26 files in openff-pablo's prepared PDB corpus load. Each refusal names the residue, and Pablo's verdict on the same files is shown. |
| `pdb_round_trip.ipynb` | Lysozyme written by `save_pdb` has 8 `CONECT` lines, one per disulfide end, and Pablo reads it with no arguments. |
| `multi_residue_fragment.ipynb` | A 31-residue peptide attached to a lysine keeps every residue's name, atoms and charge, and one bond record is enough for Pablo to read the product. |
| `deprotonate_then_attach.ipynb` | Attaching to a charged lysine keeps the +1. Deprotonating first gives the neutral product. |
| `glycan_from_pdb.ipynb` | Three GLYCAM-Web glycans read with `fragment_from_pdb`, N-linked to ASN60 of ubiquitin with the ROH hydroxyl as the leaving group. The GLYCAM residue names survive into the written file. |

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
| `pixi run lab` | Open the notebooks. |
| `pixi run verify` | Execute every notebook headless into `assets_cache/executed/`. The committed notebooks carry no outputs. |
| `pixi run dev` | Point the environment at a local mBuild checkout. |
| `pixi run fetch` / `assets` | Re-download the workshop file and rebuild the committed structures. |
| `pixi run clean` | Clear notebook outputs in place before committing. |

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

mBuild is responsible for the structure and the topology. Assigning a
force field happens in another package. Here that package is OpenFF.
