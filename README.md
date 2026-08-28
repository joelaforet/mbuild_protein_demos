# mBuild Protein Demos

**Prepare and covalently modify proteins with
[mBuild](https://github.com/mosdef-hub/mbuild), then hand off to the
parameterization engine of your choice.**

Two notebooks tell the whole story:

| Notebook | What it shows |
|---|---|
| `01_load_and_export.ipynb` | Load a protonated protein PDB into an mBuild Compound with full chemistry (formal charges, bond orders, residue identity), and export it to **GMSO**, **ParmEd**, **RDKit**, a prepared **PDB**, and **OpenFF**. |
| `02_modify_and_simulate.ipynb` | Build a fragment from scratch (star-sited SMILES), covalently `attach()` it to a lysine, export the modified PDB + bond record — then watch that output feed the OpenFF ecosystem: Pablo ingestion, NAGL charges, ff14SB + Sage via Interchange, solvation, and a short OpenMM MD run. |

mBuild owns the coordinates; the force-field assignment is deliberately
downstream (MosDef does not ship biopolymer force fields yet, so the
demo uses OpenFF for that step).

> The notebooks track the `feat/scope-trim` branch of
> [joelaforet/mbuild](https://github.com/joelaforet/mbuild), proposed
> upstream to mosdef-hub/mbuild.
> **More demos** (multi-site polymer tethering, glycans from GLYCAM
> PDBs, branched sugars, mixed force fields, a PyMOL placement movie)
> live on the [`showcase-extras`](../../tree/showcase-extras) branch.

## Getting started

### 1. Install pixi (once per machine)

[Pixi](https://pixi.sh) manages the whole software environment for
these tutorials - no conda setup needed.

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
JupyterLab *is* the tutorial environment - you do not need to switch
kernels. (If you instead use an external Jupyter or VS Code, register
the environment as a named kernel once:
`pixi run python -m ipykernel install --user --name mbuild-protein-demos`
and select `mbuild-protein-demos` in the kernel picker.)

Run notebook 01 first (it creates the protonated input), then 02.

## Inputs

`1UBQ_testProtein.cleaned.pdb` — the ubiquitin crystal structure.
Notebook 01 protonates it once with pdbfixer at pH 7; substitute your
own protein the same way.

## Notes

- The protein input must be **fully protonated** at your target pH.
  The loader errors, naming the residue and the fix, on anything it
  cannot match - it never guesses chemistry.
