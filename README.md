# mBuild Protein Demos

**Prepare and covalently modify proteins with
[mBuild](https://github.com/mosdef-hub/mbuild), then hand off to the
parameterization engine of your choice.**

Two notebooks tell the whole story.

**`01_load_and_export.ipynb` — load and export.** A protonated protein
PDB file becomes an mBuild `Compound` tree with full chemistry: every
atom identified against a Chemical Component Dictionary template, bonds
with orders, residues with formal charges and real PDB numbers. Nothing
is guessed. An atom that the templates cannot explain raises an error
that names the residue and the fix. The same object then exports to
**GMSO**, **ParmEd**, **RDKit**, a prepared **PDB** file, and **OpenFF**.

**`02_modify_and_simulate.ipynb` — modify, relax, simulate.** Seven
steps:

1. Load the protein.
2. Build a fragment from scratch as star-sited SMILES. The `*` marks
   the atom that forms the bond.
3. `protein.deprotonate(...)` takes the acidic proton off the
   attachment site. A protonated lysine side chain has no lone pair, so
   it is not the reactive form; the conjugation runs through the neutral
   free amine. The residue re-matches to the CCD variant of the result,
   and each modified lysine lowers the net formal charge of the protein
   by one.
4. `protein.attach(...)` places the fragment, removes one hydrogen per
   side, and forms the bond. The bond replaces one N-H bond, so the
   product is a neutral secondary amide.
5. Look at the modified site in 3D, then **relax the fragment as a
   movie**: 40 frames of energy minimization with the protein held
   fixed, played in the notebook with a frame slider. Frame 0 is the
   rigid placement, so the movie starts at the steric clash and shows
   the fragment settle. It takes about 15 s to compute.
6. `save_pdb` writes the modified PDB file, and `bond_records()`
   returns one plain dict per new bond.
7. Those two outputs are all that OpenFF Pablo needs. From there the
   notebook assigns charges and parameters and runs a short solvated MD
   simulation in OpenMM.

An appendix repeats the mBuild half with a 106-atom sulfonated cyanine
FRET dye, to show that nothing in the workflow depends on the size of
the fragment.

mBuild owns the coordinates. The force-field assignment is downstream on
purpose: MosDef ships no biopolymer force field yet, so the demo uses
OpenFF for that step.

## Charges and parameters

The notebook splits the partial charges between two models.
`demo_charges.assign_split_charges` gives every atom of a standard
residue its Amber **ff14SB** library charge, read from the unmodified
protein. The fragment and the modified residue take **NAGL am1bcc graph
charges** (`openff-gnn-am1bcc-0.1.0-rc.3`), computed on a capped local
model of the modification site. The seam between the two models lies on
the peptide bonds of the modified residue. The two sets do not sum to
the formal charge exactly, so the small residual spreads over the atoms
of that site.

The parameters come from the Amber **ff14SB** port plus **Sage 2.3.0**
(`openff-2.3.0.offxml`) through Interchange. The split charges go in as
preset charges, and the notebook reads them back out of the Interchange
to prove that the NAGLCharges handler of Sage 2.3.0 did not write over
them.

The split claims no more than that. The two sets come from different
fits, so the atoms across the seam are not mutually polarized, and the
library charges of the neighbouring residues do not respond to the
modification. NAGL am1bcc is the model that runs today. A different
graph model, such as AshGC, would take the place of NAGL in the same
call.

## Files

| File | What it is |
|---|---|
| `1UBQ_testProtein.cleaned.pdb` | The ubiquitin crystal structure. Notebook 01 protonates it once with pdbfixer at pH 7. Substitute your own protein the same way. |
| `1ubq_protonated.pdb` | The protonated input, written by notebook 01. |
| `1ubq_octanoyl.pdb` | The modified protein, written by notebook 02. This is the file that OpenFF Pablo reads back. |
| `demo_charges.py` | The charge split: the capped local model of the modification site, the two charge models, and the checks on the result. |
| `demo_utils.py` | The code around the demo: the NGLView views, the relaxation movie, and the openff-pablo boilerplate. The notebooks call the mBuild API themselves. |

## Getting started

### 1. Install pixi (once per machine)

[Pixi](https://pixi.sh) manages the whole software environment for these
tutorials. No conda setup is needed.

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Then restart your shell (or `source ~/.bashrc`). On Windows, see the
[pixi installation docs](https://pixi.sh/latest/#installation).

### 2. Install the tutorial environment (once per clone)

```bash
git clone https://github.com/joelaforet/mbuild_protein_demos.git
cd mbuild_protein_demos
pixi install      # solves and installs all conda dependencies
pixi run setup    # installs the mBuild biopolymers branch + openff-pablo
```

### 3. Open the tutorials

```bash
pixi run lab
```

This starts JupyterLab **inside the pixi environment**, with both
notebooks open. The default `Python 3 (ipykernel)` kernel of this
JupyterLab *is* the tutorial environment, so you do not need to switch
kernels. (If you use an external Jupyter or VS Code instead, register
the environment as a named kernel once:
`pixi run python -m ipykernel install --user --name mbuild-protein-demos`,
then select `mbuild-protein-demos` in the kernel picker.)

Run notebook 01 first: it writes the protonated input that notebook 02
reads.

### 4. Check that the demo still runs

```bash
pixi run verify
```

This executes both notebooks headless and writes the outputs back into
the files. It takes about 95 s. Exit code 0 with no error output in the
notebooks means the demo runs end to end.

## Notes

- The protein input must be **fully protonated** at your target pH. The
  loader errors, naming the residue and the fix, on anything it cannot
  match. It never guesses chemistry.
- The 3D views are NGLView widgets. They appear in JupyterLab, but not
  in a static rendering of the notebook on GitHub.
- The movie exports to an animated GIF only when the environment
  variable `DEMO_RENDER` is set. NGLView asks the browser for each
  picture, so the export needs a live front end.
- **More demos** live on the
  [`showcase-extras`](../../tree/showcase-extras) branch: one polymer
  tethered at several protein sites, glycans from GLYCAM PDB files,
  branched sugars from SMILES, mixed force fields, and a PyMOL placement
  movie.

> The notebooks track the `feat/scope-trim` branch of
> [joelaforet/mbuild](https://github.com/joelaforet/mbuild), proposed
> upstream to mosdef-hub/mbuild.
