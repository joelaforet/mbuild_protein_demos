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

**`02_modify_and_simulate.ipynb` — modify, relax, simulate.** Six steps:

1. Load the protein.
2. Build a fragment from scratch as star-sited SMILES. The `*` marks
   the atom that forms the bond.
3. `protein.attach(...)` places the fragment, removes one hydrogen per
   side, and forms the bond.
4. Look at the modified site in 3D, then **relax the fragment as a
   movie**: 40 frames of energy minimization with the protein held
   fixed, played in the notebook with a frame slider. Frame 0 is the
   rigid placement, so the movie starts at the steric clash and shows
   the fragment settle. It takes about 15 s to compute.
5. `save_pdb` writes the modified PDB file, and `bond_records()`
   returns one plain dict per new bond.
6. Those two outputs are all that OpenFF Pablo needs. From there the
   notebook assigns charges and parameters and runs a short solvated MD
   simulation in OpenMM.

An appendix repeats the mBuild half with a 106-atom sulfonated cyanine
FRET dye, to show that nothing in the workflow depends on the size of
the fragment.

mBuild owns the coordinates. The force-field assignment is downstream on
purpose: MosDef ships no biopolymer force field yet, so the demo uses
OpenFF for that step.

## Charges and parameters

The notebook assigns **NAGL am1bcc graph charges**
(`openff-gnn-am1bcc-0.1.0-rc.3`) to the whole conjugate in one call. The
parameters come from the Amber **ff14SB** port plus **Sage 2.3.0**
(`openff-2.3.0.offxml`) through Interchange.

One call over the whole conjugate is the short path, not the physically
best one: it replaces the ff14SB library charges on every protein atom
that ff14SB already describes. Splitting the two — ff14SB charges on the
protein, graph charges on the modified site — is the next step for this
notebook. A different graph model, such as AshGC, would take the place
of NAGL there.

## Files

| File | What it is |
|---|---|
| `1UBQ_testProtein.cleaned.pdb` | The ubiquitin crystal structure. Notebook 01 protonates it once with pdbfixer at pH 7. Substitute your own protein the same way. |
| `1ubq_protonated.pdb` | The protonated input, written by notebook 01. |
| `1ubq_octanoyl.pdb` | The modified protein, written by notebook 02. This is the file that OpenFF Pablo reads back. |
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
