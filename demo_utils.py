"""Helpers for the two mBuild protein demo notebooks.

The notebooks call the mBuild API directly, because showing that API is
the point of the demo. What lives here is the code around it: the
NGLView display, the relaxation movie, and the openff-pablo
boilerplate. Nothing in this file is part of mBuild, and nothing here
hides an mBuild call that a reader should see.

Three groups of functions:

- Display: ``show_protein``, ``show_movie``, ``save_gif``. They build
  NGLView widgets from PDB text. ``mbuild.Compound.visualize`` is not
  used, because its nglview branch raises TypeError on this mBuild
  branch.
- Geometry: ``relax_movie`` minimizes the attached fragments in small
  steps and records one movie frame per step.
- OpenFF handoff: ``pablo_residue_library`` and
  ``pablo_crosslink_kwargs`` turn one attachment into the residue
  definition and crosslink declaration that openff-pablo needs. They
  wrap openff-pablo and RDKit calls, not mBuild calls: the notebooks
  call the mBuild API themselves.
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


def pablo_residue_library(smiles, fragment, resname, records):
    """Build an openff-pablo residue library for one attached fragment.

    Pablo knows the CCD residues but not the fragment, so it needs one
    named residue definition plus a crosslink declaration per new bond.

    Parameters
    ----------
    smiles : str
        The same starred SMILES that built the fragment.
    fragment : mbuild.Compound
        The pristine fragment, as returned by ``prepare_fragment``.
        Its atom order matches the definition built
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


# NGL reads a color as a hex value or as a CSS color name. These are hex
# values, so the picture does not change with the color name table of
# the browser.
PROTEIN_COLOR = "#BFBFBF"  # grey
LINK_COLOR = "#00CED1"  # cyan
FRAGMENT_COLOR = "#3CB44B"  # green


def _residue_term(resnum, chain_id):
    """Return the NGL selection term for one residue.

    In the NGL selection language a residue number stands alone
    (``63``) and a chain identifier follows a colon (``63:A``). A
    protein that was read without chain identifiers has none to add.

    Parameters
    ----------
    resnum : int
        The PDB residue sequence number.
    chain_id : str
        The chain identifier, or an empty string.

    Returns
    -------
    str
        One NGL selection term.
    """
    return f"{resnum}:{chain_id}" if chain_id else str(resnum)


def _chain_ids(protein):
    """Return the chain identifiers of a protein, each one once."""
    return list(dict.fromkeys(chain.chain_id for chain in protein.chains))


def _residue_chains(protein):
    """Map each ``(residue name, residue number)`` pair to its chains.

    ``Protein.bond_records()`` names a residue by name and number and
    carries no chain identifier, so a record alone cannot address a
    residue of a protein of several chains. This map supplies the
    chain. The value is a list, because one name and number can occur
    in more than one chain.

    Parameters
    ----------
    protein : mbuild.biopolymers.Protein
        The protein to walk.

    Returns
    -------
    dict
        Map of ``(residue name, residue number)`` to a list of chain
        identifiers.
    """
    mapping = {}
    for chain_id in _chain_ids(protein):
        for residue in protein.residues(chain_id=chain_id):
            key = (residue.name, residue.resnum)
            mapping.setdefault(key, []).append(chain_id)
    return mapping


def _modification_selections(protein):
    """Return the NGL selections of the linkage and the fragment residues.

    ``Residue.hetatm`` marks every attached fragment residue.
    ``Protein.bond_records()`` names the residues on both sides of each
    recorded inter-residue bond; those are the residues that take part
    in a covalent linkage. The fragment residues are taken out of the
    linkage set, so each residue carries one color.

    Parameters
    ----------
    protein : mbuild.biopolymers.Protein
        The protein to read.

    Returns
    -------
    tuple of (str, str)
        The linkage selection and the fragment selection. Either is an
        empty string when the protein holds no such residue.
    """
    chains = _residue_chains(protein)
    fragment_keys = set()
    fragment_terms = []
    for chain_id in _chain_ids(protein):
        for residue in protein.residues(chain_id=chain_id):
            if residue.hetatm:
                fragment_keys.add((residue.name, residue.resnum))
                fragment_terms.append(_residue_term(residue.resnum, chain_id))

    link_terms = []
    seen = set()
    for record in protein.bond_records():
        for key in zip(record["residue_names"], record["residue_numbers"]):
            if key in fragment_keys or key in seen:
                continue
            seen.add(key)
            chain_ids = chains.get(key, [])
            if len(chain_ids) == 1:
                link_terms.append(_residue_term(key[1], chain_ids[0]))
                continue
            # The name and the number match a residue of more than one
            # chain, and the record holds no chain identifier. The bare
            # number selects that residue in every chain.
            logger.warning(
                "Residue %s %s occurs in chains %s, and a bond record "
                "carries no chain identifier. The view draws that residue "
                "in every one of those chains. Pass link_selection to draw "
                "one of them.",
                key[0],
                key[1],
                chain_ids or "(none)",
            )
            link_terms.append(str(key[1]))
    return " or ".join(link_terms), " or ".join(fragment_terms)


def _resolve_selections(protein, link_selection, fragment_selection):
    """Fill the selections the caller left out from the protein itself.

    Parameters
    ----------
    protein : mbuild.biopolymers.Protein or None
        The protein to derive from. None when the source is a PDB file
        or PDB text, which carries no bond record.
    link_selection, fragment_selection : str or None
        The values the caller passed. None asks for the derived value.

    Returns
    -------
    tuple of (str, str)
        The linkage selection and the fragment selection.
    """
    if link_selection is not None and fragment_selection is not None:
        return link_selection, fragment_selection
    derived_link, derived_fragment = "", ""
    if protein is not None:
        derived_link, derived_fragment = _modification_selections(protein)
    return (
        derived_link if link_selection is None else link_selection,
        derived_fragment if fragment_selection is None else fragment_selection,
    )


def _style_view(view, link_selection, fragment_selection):
    """Draw a grey cartoon, cyan linkage sticks and green fragment sticks.

    ``licorice`` is the NGL name of the stick representation. NGLView
    sends a representation name straight to NGL, and NGL draws nothing
    for a name that its registry does not hold. ``ball_and_stick`` is
    such a name: the registered name is ``ball+stick``, and the
    ``add_ball_and_stick`` shortcut of NGLView translates it, but
    ``add_representation`` does not.

    The NGL ``protein`` selection covers the standard polymer residues
    only. An attached fragment carries a residue name outside that set,
    so the cartoon leaves the fragment out. The stick representation on
    the fragment selection is the only drawing of the fragment.

    NGL holds no template for a fragment residue, so it finds the bonds
    inside that residue by interatomic distance.

    Parameters
    ----------
    view : nglview.NGLWidget
        The widget to style.
    link_selection : str
        NGL selection of the linkage residues, or an empty string.
    fragment_selection : str
        NGL selection of the fragment residues, or an empty string.
    """
    view.clear_representations()
    view.add_representation("cartoon", selection="protein", color=PROTEIN_COLOR)
    if link_selection:
        view.add_representation("licorice", selection=link_selection, color=LINK_COLOR)
    if fragment_selection:
        view.add_representation(
            "licorice", selection=fragment_selection, color=FRAGMENT_COLOR
        )
        # The modification is the subject of the picture, so the view
        # centers on it.
        view.center(selection=fragment_selection)


def show_protein(
    source,
    link_selection=None,
    fragment_selection=None,
    width="700px",
    height="500px",
):
    """Show a protein as a grey cartoon with the modification in sticks.

    A ``Protein`` describes its own modifications, so the default view
    needs no further argument. ``Protein.bond_records()`` names the
    residues of every recorded inter-residue bond, and ``Residue.hetatm``
    marks the attached fragment residues. The fragment residues become
    green sticks, the other residues of those bonds become cyan sticks,
    and the view centers on the fragment.

    A PDB path and PDB text hold no bond record, so the derivation does
    not run for those sources. Pass ``link_selection`` and
    ``fragment_selection`` with such a source, or the view shows the
    grey cartoon alone.

    Parameters
    ----------
    source : mbuild.biopolymers.Protein or str
        A protein, the path of a PDB file, or PDB text. A protein is
        written to ``SCRATCH_PDB`` first, so the view shows exactly the
        file that a downstream loader reads.
    link_selection : str, optional
        NGL selection of the residues that take part in a covalent
        linkage, drawn as cyan sticks. In the NGL selection language a
        residue number stands alone (``63``), a chain identifier
        follows a colon (``63:A``), a residue name goes in brackets
        (``[OC8]``), and terms join with ``or`` and ``and``. Default:
        derived from a ``Protein``, empty for any other source.
    fragment_selection : str, optional
        NGL selection of the attached fragment, drawn as green sticks
        and centered in the view. Default: the HETATM residues of a
        ``Protein``, empty for any other source.
    width, height : str, optional
        Widget size as a CSS length.

    Returns
    -------
    nglview.NGLWidget
        Display it as the last expression of a notebook cell.
    """
    protein = source if hasattr(source, "bond_records") else None
    text = _pdb_text(source)
    view = nglview.NGLWidget(nglview.TextStructure(text, ext="pdb"))
    _style_view(view, *_resolve_selections(protein, link_selection, fragment_selection))
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
    protein=None,
    link_selection=None,
    fragment_selection=None,
    width="700px",
    height="500px",
):
    """Show a relaxation movie with a frame slider.

    The picture is the one of ``show_protein``: a grey cartoon, cyan
    sticks on the residues of each covalent linkage, green sticks on the
    attached fragment, and the view centered on the fragment.

    The source of a movie is a multi-MODEL PDB file, which holds no bond
    record. Pass the ``Protein`` that the movie was made from, and this
    function derives the two selections from it, as ``show_protein``
    does. Pass the selections instead when no protein object is at hand.

    Parameters
    ----------
    source : str or tuple
        The path of a multi-MODEL PDB file, or the tuple that
        ``relax_movie`` returns. The tuple form skips re-reading the
        coordinates from the file.
    protein : mbuild.biopolymers.Protein, optional
        The protein whose relaxation the movie shows. It supplies the
        default selections. The atom order of the movie is the atom
        order of that protein, so the selections match the frames.
    link_selection, fragment_selection, width, height
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
    _style_view(view, *_resolve_selections(protein, link_selection, fragment_selection))
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
