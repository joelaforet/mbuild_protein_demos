"""Helpers for the two mBuild protein demo notebooks.

The notebooks stay short because the display code, the relaxation movie,
and the OpenFF Pablo boilerplate live here. Nothing in this file is part
of mBuild: it uses only the public API of ``mbuild.biopolymers`` and
``mbuild.simulation``.

Three groups of functions:

- Display: ``show_protein``, ``show_movie``, ``save_gif``. They build
  NGLView widgets from PDB text. ``mbuild.Compound.visualize`` is not
  used, because its nglview branch raises TypeError on this mBuild
  branch.
- Geometry: ``relax_movie`` minimizes the attached fragments in small
  steps and records one movie frame per step.
- OpenFF handoff: ``functionalize``, ``pablo_residue_library``, and
  ``pablo_crosslink_kwargs`` turn one attachment into the residue
  definition and crosslink declaration that openff-pablo needs.
"""

import logging
import os

import nglview
import numpy as np
from nglview.base_adaptor import Structure, Trajectory

logger = logging.getLogger(__name__)

# Scratch file for one set of coordinates. Protein.save_pdb writes a
# file, so the display helpers write and read this path instead of
# holding a PDB string in memory.
SCRATCH_PDB = "_frame.pdb"


# ----------------------------------------------------------------------
# OpenFF handoff
# ----------------------------------------------------------------------
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


def functionalize(protein, smiles, resname, site, relax=True):
    """Attach a starred-SMILES fragment to one site of a protein.

    Parameters
    ----------
    protein : mbuild.biopolymers.Protein
        The protein to modify in place.
    smiles : str
        SMILES of the fragment with one ``*``. The star marks the atom
        that forms the new bond and becomes the leaving hydrogen.
    resname : str
        Residue name given to the attached fragment. Use it later as
        the Pablo residue-definition key.
    site : dict
        Keyword arguments for ``Protein.attach``, for example
        ``dict(resnum=63, atom_name="NZ", chain_id="A")``.
    relax : bool, optional, default=True
        Passed to ``Protein.attach``. Pass ``False`` when a later
        ``relax_movie`` call is to show the relaxation.

    Returns
    -------
    tuple of (mbuild.Compound, mbuild.biopolymers.protein.InterResidueBond)
        The pristine fragment and the recorded bond. ``attach`` clones
        the fragment and removes one hydrogen from the clone, so the
        returned fragment still holds every atom of the free molecule.
        ``pablo_residue_library`` needs that complete atom list.
    """
    from mbuild.biopolymers import prepare_fragment

    fragment = prepare_fragment(smiles, resname)
    record = protein.attach(fragment, relax=relax, **site)
    return fragment, record


def pablo_residue_library(smiles, fragment, resname, records):
    """Build an openff-pablo residue library for one attached fragment.

    Pablo knows the CCD residues but not the fragment, so it needs one
    named residue definition plus a crosslink declaration per new bond.

    Parameters
    ----------
    smiles : str
        The same starred SMILES that built the fragment.
    fragment : mbuild.Compound
        The pristine fragment, as returned by ``prepare_fragment`` or
        ``functionalize``. Its atom order matches the definition built
        from the SMILES, so the atom names transfer by position. Do not
        pass the fragment residue taken out of the protein: ``attach``
        removed its leaving hydrogen, so the two atom lists differ in
        length and the names shift by one atom.
    resname : str
        Residue name of the fragment in the PDB file.
    records : list of dict
        Output of ``Protein.bond_records()``.

    Returns
    -------
    openff.pablo.ResidueDefinitionLibrary
        Pass it to ``topology_from_pdb(..., residue_library=...)``.
    """
    from openff.pablo import STD_CCD_CACHE, ResidueDefinition
    from openff.toolkit import Molecule
    from rdkit import Chem

    # Replace the star with hydrogen exactly as prepare_fragment does,
    # so the definition and the fragment hold the same atoms in the
    # same order.
    star = Chem.RWMol(Chem.MolFromSmiles(smiles))
    for atom in star.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
    mol = star.GetMol()
    Chem.SanitizeMol(mol)
    offmol = Molecule.from_rdkit(Chem.AddHs(mol), allow_undefined_stereo=True)

    particles = list(fragment.particles())
    if len(particles) != offmol.n_atoms:
        raise ValueError(
            f"The fragment has {len(particles)} atoms but the molecule built "
            f"from {smiles!r} has {offmol.n_atoms}. Pass the pristine "
            "fragment, not the residue that attach() put into the protein."
        )
    for atom, particle in zip(offmol.atoms, particles):
        atom.name = particle.name

    library = STD_CCD_CACHE.with_(
        {resname: [ResidueDefinition.from_molecule(offmol, residue_name=resname)]}
    )
    for record in records:
        library = library.with_crosslink(**pablo_crosslink_kwargs(record))
    return library


# ----------------------------------------------------------------------
# Display
# ----------------------------------------------------------------------
class _TextFrames(Structure, Trajectory):
    """One PDB text plus a coordinate array, as one NGLView component.

    NGLView takes topology from a ``Structure`` and coordinates from a
    ``Trajectory``. One object that is both lets the widget show a
    frame slider over coordinates held in memory. The frames must be in
    Angstrom, the unit of the PDB coordinate columns; mBuild positions
    are in nanometers, so callers multiply by 10.

    Parameters
    ----------
    text : str
        PDB text of one frame, including the CONECT records. It fixes
        the atom order and the bonds.
    frames : array_like
        Coordinates, shape (n_frames, n_atoms, 3), in Angstrom.
    """

    def __init__(self, text, frames):
        Structure.__init__(self)
        Trajectory.__init__(self)
        self.ext = "pdb"
        self.params = {}
        self._text = text
        self._frames = np.asarray(frames, dtype=np.float32)
        if self._frames.ndim != 3 or self._frames.shape[2] != 3:
            raise ValueError(
                "frames must have shape (n_frames, n_atoms, 3), got "
                f"{self._frames.shape}."
            )

    def get_structure_string(self):
        """Return the PDB text of the reference frame."""
        return self._text

    def get_coordinates(self, index):
        """Return the coordinates of frame ``index`` in Angstrom."""
        return self._frames[index]

    @property
    def n_frames(self):
        """Return the number of frames."""
        return len(self._frames)


def _pdb_text(source, scratch=SCRATCH_PDB):
    """Return PDB text from a Protein, a file path, or PDB text."""
    if hasattr(source, "save_pdb"):
        source.save_pdb(scratch, overwrite=True)
        source = scratch
    text = str(source)
    if "\n" not in text:
        with open(text) as handle:
            return handle.read()
    return text


def _highlight_selection(resnums, resnames):
    """Build one NGL selection string from residue numbers and names.

    Residue numbers may be plain integers (``63``) or, when numbers
    repeat across chains, ``"number:chain"`` strings (``"63:A"``).
    """
    parts = [f"[{resname}]" for resname in resnames or ()]
    parts += [str(resnum) for resnum in resnums or ()]
    return " or ".join(parts)


def _style_view(view, highlight, chain_colors):
    """Draw a cartoon backbone plus ball-and-stick on the highlight."""
    view.clear_representations()
    view.add_representation("cartoon", selection="protein", color=chain_colors)
    if highlight:
        view.add_representation("ball_and_stick", selection=highlight)
        view.center(selection=highlight)


def show_protein(
    source,
    highlight_resnums=None,
    highlight_resnames=None,
    chain_colors="chainname",
    width="700px",
    height="500px",
):
    """Show a protein as a cartoon, with chosen residues in atom detail.

    Parameters
    ----------
    source : mbuild.biopolymers.Protein or str
        A protein, the path of a PDB file, or PDB text. A protein is
        written to ``SCRATCH_PDB`` first, so the view shows exactly the
        file that a downstream loader reads.
    highlight_resnums : iterable, optional
        Residue numbers to draw as ball-and-stick. Use ``"63:A"`` to
        name the chain.
    highlight_resnames : iterable of str, optional
        Residue names to draw as ball-and-stick, for example
        ``["OCT"]`` for an attached fragment.
    chain_colors : str, optional, default="chainname"
        NGL color scheme for the cartoon.
    width, height : str, optional
        Widget size as a CSS length.

    Returns
    -------
    nglview.NGLWidget
        Display it as the last expression of a notebook cell.
    """
    text = _pdb_text(source)
    view = nglview.NGLWidget(nglview.TextStructure(text, ext="pdb"))
    highlight = _highlight_selection(highlight_resnums, highlight_resnames)
    _style_view(view, highlight, chain_colors)
    view.layout.width = width
    view.layout.height = height
    return view


def _split_models(text):
    """Split a multi-MODEL PDB text into a frame text and coordinates.

    Returns the PDB text of the first model with the CONECT footer
    appended, plus an (n_frames, n_atoms, 3) array in Angstrom.
    """
    models = []
    footer = []
    current = None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            models.append(current)
            current = None
        elif current is not None:
            current.append(line)
        elif line.startswith("CONECT"):
            footer.append(line)
    if not models:
        raise ValueError(
            "This PDB text holds no MODEL records, so it is not a movie. "
            "Pass the file that relax_movie wrote."
        )
    frames = np.array(
        [
            [
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                for line in model
                if line.startswith(("ATOM", "HETATM"))
            ]
            for model in models
        ],
        dtype=np.float32,
    )
    frame_text = "\n".join(models[0] + footer + ["END"]) + "\n"
    return frame_text, frames


def show_movie(
    source,
    highlight_resnums=None,
    highlight_resnames=None,
    chain_colors="chainname",
    width="700px",
    height="500px",
):
    """Show a relaxation movie with a frame slider.

    Parameters
    ----------
    source : str or tuple
        The path of a multi-MODEL PDB file, or the tuple that
        ``relax_movie`` returns. The tuple form skips re-reading the
        coordinates from the file.
    highlight_resnums, highlight_resnames, chain_colors, width, height
        As in ``show_protein``.

    Returns
    -------
    nglview.NGLWidget
        Press play in the widget to run the movie. Under nbconvert
        there is no front end, so the widget shows as a placeholder and
        the cell still succeeds.
    """
    if isinstance(source, (str, os.PathLike)):
        frame_text, frames = _split_models(_pdb_text(source))
    else:
        path, frames = source[0], np.asarray(source[1])
        frame_text, _ = _split_models(_pdb_text(path))
    view = nglview.NGLWidget(_TextFrames(frame_text, frames))
    highlight = _highlight_selection(highlight_resnums, highlight_resnames)
    _style_view(view, highlight, chain_colors)
    view.layout.width = width
    view.layout.height = height
    return view


def save_gif(view, frames, path, duration=120, loop=0):
    """Write the frames of a movie widget to an animated GIF.

    ``NGLWidget.render_image`` asks the browser for a picture, so it
    returns an empty image without a live front end. A headless run
    (``pixi run verify``) has no front end. This function therefore
    writes a file only when the environment variable ``DEMO_RENDER`` is
    set to a non-empty value. Set it before starting JupyterLab, then
    run the cell with the widget already displayed.

    Parameters
    ----------
    view : nglview.NGLWidget
        A displayed movie widget, as returned by ``show_movie``.
    frames : int or array_like
        The number of frames, or the frame array whose length is that
        number.
    path : str
        Path of the GIF file to write.
    duration : int, optional, default=120
        Time per frame in milliseconds.
    loop : int, optional, default=0
        GIF loop count; 0 repeats forever.

    Returns
    -------
    str or None
        The path written, or None when DEMO_RENDER is unset.
    """
    import time

    if not os.environ.get("DEMO_RENDER"):
        logger.info(
            "DEMO_RENDER is unset, so no GIF was written. Set DEMO_RENDER=1 "
            "in a JupyterLab session to render %s.",
            path,
        )
        return None

    from io import BytesIO

    from PIL import Image

    n_frames = frames if isinstance(frames, int) else len(frames)
    images = []
    for index in range(n_frames):
        view.frame = index
        image = view.render_image()
        deadline = time.monotonic() + 30.0
        while not image.value and time.monotonic() < deadline:
            time.sleep(0.1)
        if not image.value:
            raise RuntimeError(
                f"The front end returned no picture for frame {index} within "
                "30 s. Display the widget in the notebook before calling "
                "save_gif."
            )
        images.append(Image.open(BytesIO(image.value)).convert("P"))
    images[0].save(
        path,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop,
    )
    return path


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------
def _capture(protein, scratch):
    """Return the current atom lines and the CONECT footer."""
    protein.save_pdb(scratch, overwrite=True)
    with open(scratch) as handle:
        lines = handle.read().splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM", "TER"))]
    footer = [line for line in lines if line.startswith("CONECT")]
    return atoms, footer


def _step_schedule(n_frames, steps_per_frame):
    """Return the minimization iterations for each frame after frame 0."""
    n_steps = int(n_frames) - 1
    if n_steps < 1:
        raise ValueError(f"n_frames must be at least 2, got {n_frames}.")
    if isinstance(steps_per_frame, int):
        return [steps_per_frame] * n_steps
    schedule = [int(steps) for steps in steps_per_frame]
    if len(schedule) != n_steps:
        raise ValueError(
            f"steps_per_frame holds {len(schedule)} entries but n_frames="
            f"{n_frames} needs {n_steps}: frame 0 is the placement, and every "
            "later frame follows one minimization. Pass an int to use the "
            "same number of iterations for every frame."
        )
    return schedule


def relax_movie(
    protein,
    out_path,
    residues=None,
    n_frames=30,
    steps_per_frame=10,
    tolerance=1.0,
    platform="CPU",
    scratch=SCRATCH_PDB,
):
    """Minimize the attached fragments and record the motion as frames.

    This is ``Protein.relax_fragments`` split into steps. It builds one
    ``OpenMMSimulation`` with mBuild's generic parameters, gives every
    atom outside the chosen residues zero mass (OpenMM holds a
    zero-mass atom still), then minimizes in short bursts and saves the
    coordinates after each burst. The physics is the same as
    ``relax_fragments``; only the number of minimization calls differs.

    One simulation for the whole movie costs about half the time of one
    ``relax_fragments`` call per frame, because the OpenMM system is
    built once. Measured on ubiquitin plus an octanoyl group: 14.9 s for
    30 frames, against 25.4 s for repeated ``relax_fragments``.

    Most of the motion happens in the first frames. Pass a sequence to
    ``steps_per_frame`` to spend more frames on that part, for example
    ``[2] * 10 + [10] * 29`` for 40 frames.

    Parameters
    ----------
    protein : mbuild.biopolymers.Protein
        The protein to relax. Its coordinates change in place.
    out_path : str
        Path of the multi-MODEL PDB file to write.
    residues : iterable of mbuild.biopolymers.Residue, optional
        The residues allowed to move. Default: every HETATM residue,
        which is every attached fragment. This matches the default of
        ``relax_fragments``.
    n_frames : int, optional, default=30
        Number of frames, counting frame 0. Frame 0 holds the rigid
        placement, before any minimization.
    steps_per_frame : int or sequence of int, optional, default=10
        Maximum minimization iterations between frames. A sequence must
        hold ``n_frames - 1`` entries.
    tolerance : float, optional, default=1.0
        Force tolerance in kJ/mol/nm. Lower than the
        ``relax_fragments`` default of 50, so the short bursts keep
        making progress instead of stopping early.
    platform : str, optional, default="CPU"
        OpenMM platform name.
    scratch : str, optional
        Path of the single-frame file the writer reuses.

    Returns
    -------
    tuple of (str, numpy.ndarray, numpy.ndarray)
        The path written, the frames as an (n_frames, n_atoms, 3) array
        in Angstrom, and the potential energy in kJ/mol after each
        minimization (``n_frames - 1`` values). Pass the whole tuple to
        ``show_movie``.
    """
    from mbuild.simulation import OpenMMSimulation

    targets = (
        list(residues)
        if residues is not None
        else [residue for residue in protein.residues() if residue.hetatm]
    )
    if not targets:
        raise ValueError(
            "This protein has no HETATM residue to relax. Attach a fragment "
            "first, or name the residues to move through the residues "
            "argument."
        )
    schedule = _step_schedule(n_frames, steps_per_frame)

    mobile = set()
    for residue in targets:
        mobile.update(residue.particles())
    simulation = OpenMMSimulation(
        protein, forcefield=None, kick=False, platform=platform
    )
    for index, particle in enumerate(protein.particles()):
        if particle not in mobile:
            simulation.system.setParticleMass(index, 0.0)

    # Frame 0 is the placement, so capture it before the first minimize.
    models = []
    frames = []
    atoms, footer = _capture(protein, scratch)
    models.append(atoms)
    frames.append(protein.xyz * 10.0)
    for steps in schedule:
        simulation.minimize(n_steps=steps, tolerance=tolerance)
        atoms, footer = _capture(protein, scratch)
        models.append(atoms)
        frames.append(protein.xyz * 10.0)

    with open(out_path, "w") as handle:
        for number, atom_lines in enumerate(models, start=1):
            handle.write(f"MODEL     {number:4d}\n")
            handle.write("\n".join(atom_lines) + "\n")
            handle.write("ENDMDL\n")
        # The bonds do not change during a minimization, so the CONECT
        # records apply to every model and appear once, as a footer.
        if footer:
            handle.write("\n".join(footer) + "\n")
        handle.write("END\n")

    energies = np.array(
        [entry["potential_energy"] for entry in simulation.energies], dtype=float
    )
    return out_path, np.array(frames, dtype=np.float32), energies
