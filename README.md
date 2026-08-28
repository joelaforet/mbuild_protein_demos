# mBuild Protein Demos

Tutorials for preparing and **covalently modifying proteins with
[mBuild](https://github.com/mosdef-hub/mbuild)**, then handing the
modified structure off to standard parameterization engines:
**GMSO/foyer, ParmEd, RDKit, and OpenFF (Pablo + Interchange)**.

The workflow in one breath: load a protonated protein PDB with full
chemistry (formal charges, bond orders, residue identity), attach any
fragment written as star-sited SMILES / SDF / PDB, get clash-checked and
relaxed coordinates, and export a prepared PDB plus neutral bond
records that downstream tools consume.

> These tutorials track the `feat/scope-trim` branch of
> [joelaforet/mbuild](https://github.com/joelaforet/mbuild), proposed
> upstream to mosdef-hub/mbuild.

## Getting started

### 1. Install pixi (once per machine)

[Pixi](https://pixi.sh) manages the whole software environment for
these tutorials — no conda setup needed.

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

This starts JupyterLab **inside the pixi environment**, with the
quickstart notebook open. The default `Python 3 (ipykernel)` kernel of
this JupyterLab *is* the tutorial environment — you do not need to
switch kernels. (If you instead use an external Jupyter or VS Code,
register the environment as a named kernel once:
`pixi run python -m ipykernel install --user --name mbuild-protein-demos`
and select `mbuild-protein-demos` in the kernel picker.)

Everything also runs without notebooks:

```bash
pixi run showcase   # the full triply-modified-ubiquitin build + OpenFF ingestion
```

## The tutorials

| File | What it shows |
|---|---|
| `quickstart_functionalize.ipynb` | **Start here.** PDB + star-sited SMILES fragment → modified PDB + bond records → OpenFF, RDKit, GMSO, ParmEd handoffs. ~15 lines of user code. |
| `mbuild_pablo_ptm_demo.ipynb` | The maximal showcase: a polymer trimer (amide at LYS63), a GLYCAM glycan from PDB (ASN60), and a branched glycan from SMILES (SER20) on one ubiquitin, loaded through OpenFF Pablo. |
| `maximal_validation.py` | The same showcase as a plain script (`pixi run showcase`). |
| `tether_movie.py` / `tether_movie.pdb` | Generates (and ships) a 47-state PyMOL movie of a polymer being tethered to two lysines: rigid placement → rotation → shear → relaxation. `pymol tether_movie.pdb`, then `mplay`. |
| `mbuild_extras.py` | Fork-only helpers used by the showcase: `fragment_from_pdb` (legacy PDB fragments), `attach_multi` (multi-site tethering), and the OpenFF Pablo formatting glue for bond records. |

## Inputs

`1UBQ_testProtein.cleaned.pdb` (ubiquitin crystal structure),
`1ubq_protonated.pdb` and `6m03_protonated.pdb` (pdbfixer, pH 7), and
`glycam_G42666HT_CONECT.pdb` (a GLYCAM glycan with CONECT records).
The quickstart shows the one-time pdbfixer protonation step, so you
can substitute your own protein.

## Notes

- The protein input must be **fully protonated** at your target pH
  (pdbfixer or reduce). The loader errors, naming the residue and the
  fix, on anything it cannot match — it never guesses chemistry.
- Force-field assignment is deliberately **out of scope for mBuild**:
  the tutorials hand off to OpenFF here, and the same exports feed
  GMSO/foyer or ParmEd workflows.
