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

The OpenFF hand-off, the residue definition and the crosslink declaration
that openff-pablo needs, is written out in the notebooks themselves, so
that the reader sees every step. The notebooks call the mBuild API
themselves too.
"""

import json
import logging
import os

import nglview
import ipywidgets as widgets
import numpy as np
from nglview.base_adaptor import Structure, Trajectory

logger = logging.getLogger(__name__)

# Scratch file for one set of coordinates. Protein.save_pdb writes a
# file, so the display helpers write and read this path instead of
# holding a PDB string in memory.
SCRATCH_PDB = "_frame.pdb"


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
        # The frame goes to PyMOL and NGL for visualization only, so it
        # needs no bond-records file.
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


def _respect_fragment_connections(view, text):
    """Keep explicit fragment crosslinks and remove NGL's proximity guesses.

    NGL supplements CONECT with distance-based bonds between nonstandard
    residues and nearby atoms. Frame zero can contain clashes, so those
    guesses can produce false lysine-fragment bonds. Only fragments covered
    by CONECT are filtered; ordinary protein bond perception is retained.
    """
    hetero = set()
    explicit = set()
    connected = set()
    for line in text.splitlines():
        if line.startswith("HETATM"):
            hetero.add(int(line[6:11]))
        elif line.startswith("CONECT"):
            serials = [int(line[i : i + 5]) for i in range(6, len(line.rstrip()), 5)]
            connected.update(serials)
            explicit.update(
                f"{min(serials[0], other)}:{max(serials[0], other)}"
                for other in serials[1:]
            )
    fragment_atoms = hetero & connected
    if not fragment_atoms:
        return
    # This callback runs after the component loads, before styling. It only
    # changes the view's bond graph; the protein and movie data stay intact.
    view._execute_js_code(
        """
        const structure = this.stage.compList[0].structure;
        const explicit = new Set(EXPLICIT);
        const fragment = new Set(FRAGMENT);
        const store = structure.bondStore;
        const a = structure.getAtomProxy(), b = structure.getAtomProxy();
        const seen = new Set();
        let count = 0;
        for (let i = 0; i < store.count; i++) {
            a.index = store.atomIndex1[i]; b.index = store.atomIndex2[i];
            const serialKey = Math.min(a.serial, b.serial) + ":" + Math.max(a.serial, b.serial);
            const atomKey = Math.min(a.index, b.index) + ":" + Math.max(a.index, b.index);
            const crossesFragment = a.residueIndex !== b.residueIndex &&
                (fragment.has(a.serial) || fragment.has(b.serial));
            if (seen.has(atomKey) || (crossesFragment && !explicit.has(serialKey))) continue;
            seen.add(atomKey);
            store.atomIndex1[count] = a.index;
            store.atomIndex2[count] = b.index;
            store.bondOrder[count] = store.bondOrder[i];
            count++;
        }
        store.count = count;
        structure.finalizeBonds();
        structure.refreshPosition();
    """.replace("EXPLICIT", json.dumps(sorted(explicit))).replace(
            "FRAGMENT", json.dumps(sorted(fragment_atoms))
        )
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

    The PDB CONECT records supply the covalent crosslink. A single stick
    representation includes both ends, with colors assigned by selection.

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
    selections = [term for term in (link_selection, fragment_selection) if term]
    if selections:
        # NGL only draws a bond when both endpoints belong to the same
        # representation. Separate residue representations hide the crosslink.
        colors = []
        if fragment_selection:
            colors.append([FRAGMENT_COLOR, fragment_selection])
        if link_selection:
            colors.append([LINK_COLOR, link_selection])
        colors.append([PROTEIN_COLOR, "*"])
        scheme = nglview.color._ColorScheme(colors, f"modification_{view.model_id}")
        view.add_representation(
            "licorice",
            selection=" or ".join(f"({term})" for term in selections),
            color=scheme,
        )
    if fragment_selection:
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
    _respect_fragment_connections(view, text)
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


class MovieView(widgets.VBox):
    """NGL canvas with permanently visible playback and frame controls."""

    def __init__(self, ngl_widget):
        self.ngl_widget = ngl_widget
        self.play = widgets.Play(min=0, max=ngl_widget.max_frame, interval=100)
        self.slider = widgets.IntSlider(
            min=0,
            max=ngl_widget.max_frame,
            description="Frame",
            continuous_update=True,
            layout=widgets.Layout(width="400px"),
        )
        # Keep links alive and run them in the browser, including playback.
        self._links = [
            widgets.jslink((self.play, "value"), (self.slider, "value")),
            widgets.jslink((self.slider, "value"), (ngl_widget, "frame")),
        ]
        super().__init__([ngl_widget, widgets.HBox([self.play, self.slider])])

    @property
    def frame(self):
        return self.ngl_widget.frame

    @frame.setter
    def frame(self, value):
        self.ngl_widget.frame = value


class RelaxationMovie(tuple):
    """Unpackable movie data that displays a player as a cell's last expression."""

    def __new__(cls, path, frames, energies, protein):
        result = super().__new__(cls, (path, frames, energies))
        result.protein = protein
        result._view = None
        return result

    def _repr_mimebundle_(self, **kwargs):
        if self._view is None:
            self._view = show_movie(self, protein=self.protein)
        return self._view._repr_mimebundle_(**kwargs)


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
        The path of a multi-MODEL PDB file, or a tuple ``(path, frames)``
        as ``relax_movie`` returns, where ``frames`` is an array of shape
        (n_frames, n_atoms, 3) in Angstrom and ``path`` supplies the
        topology; a plain single-frame PDB file serves there too.
    protein : mbuild.biopolymers.Protein, optional
        The protein whose relaxation the movie shows. It supplies the
        default selections. The atom order of the movie is the atom
        order of that protein, so the selections match the frames.
    link_selection, fragment_selection, width, height
        As in ``show_protein``.

    Returns
    -------
    MovieView
        Press play below the canvas or drag the frame slider. Under nbconvert
        there is no front end, so the widget shows as a placeholder and
        the cell still succeeds.
    """
    if protein is None and isinstance(source, RelaxationMovie):
        protein = source.protein
    if isinstance(source, (str, os.PathLike)):
        frame_text, frames = _split_models(_pdb_text(source))
    else:
        # Frames given directly: the path only supplies the topology, so
        # it may be a plain single-frame PDB file as well as a movie.
        path, frames = source[0], np.asarray(source[1])
        text = _pdb_text(path)
        frame_text = _split_models(text)[0] if "\nMODEL" in text or text.startswith("MODEL") else text
    view = nglview.NGLWidget(_TextFrames(frame_text, frames))
    _respect_fragment_connections(view, frame_text)
    _style_view(view, *_resolve_selections(protein, link_selection, fragment_selection))
    view.layout.width = width
    view.layout.height = height
    return MovieView(view)


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
    view : MovieView or nglview.NGLWidget
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

    view = getattr(view, "ngl_widget", view)
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
    # The frame goes to PyMOL and NGL for visualization only, so it needs
    # no bond-records file.
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
    RelaxationMovie
        The path written, the frames as an (n_frames, n_atoms, 3) array
        in Angstrom, and the potential energy in kJ/mol after each
        minimization (``n_frames - 1`` values). Pass the whole tuple to
        ``show_movie``. As the last expression of a live notebook cell, the
        result also displays its own player with play/pause and a slider.
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
    return RelaxationMovie(
        out_path, np.array(frames, dtype=np.float32), energies, protein
    )


# ----------------------------------------------------------------------
# Point-mutation helpers (notebook 3)


def residue_to_rdkit(residue, remove_hs=True):
    """Return an RDKit molecule of one detached Residue, for 2D drawing.

    The residue's particles give the elements, ``atom_formal_charges``
    gives the formal charges, and the bonds carry their orders, so the
    molecule sanitizes without guessing anything. Hydrogens are removed
    by default because a 2D depiction reads better without them.

    Parameters
    ----------
    residue : mbuild.biopolymers.Residue
        A detached residue, such as ``fragment_from_ccd(code, "CB")``.
    remove_hs : bool, optional, default=True
        Drop the hydrogens from the returned molecule.

    Returns
    -------
    rdkit.Chem.Mol
    """
    from rdkit import Chem

    orders = {1.0: Chem.BondType.SINGLE, 2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE}
    mol = Chem.RWMol()
    index = {}
    for particle in residue.particles():
        atom = Chem.Atom(particle.element.symbol)
        atom.SetFormalCharge(residue.atom_formal_charges.get(particle.name, 0))
        atom.SetNoImplicit(True)
        index[particle] = mol.AddAtom(atom)
    for p1, p2, data in residue.bonds(return_bond_order=True):
        mol.AddBond(index[p1], index[p2], orders[float(data["bond_order"])])
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol) if remove_hs else mol


def alpha_carbon_cip(protein, resnum, chain_id):
    """Return RDKit's CIP label (R or S) of one residue's alpha carbon.

    The label is assigned from the 3D coordinates of the exported
    molecule, so it reads the geometry the written file will carry.
    """
    from rdkit import Chem

    mol = protein.to_rdkit()
    Chem.AssignStereochemistryFrom3D(mol)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if (
            info.GetResidueNumber() == resnum
            and info.GetChainId() == chain_id
            and info.GetName().strip() == "CA"
        ):
            return atom.GetPropsAsDict().get("_CIPCode", "none")
    raise KeyError(f"No CA in residue {resnum} of chain {chain_id}.")
