# @Author: Tristan Croll
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
import sys, os, glob
import ctypes

from . import clipper

from chimerax.atomic import molc, structure
# from chimerax.atomic.molc import CFunctions, string, cptr, pyobject, \
#     set_c_pointer, pointer, size_t

CFunctions = molc.CFunctions
string = molc.string
cptr = molc.cptr
pyobject = molc.pyobject
set_c_pointer = molc.set_c_pointer
pointer = molc.pointer
size_t = molc.size_t

dpath = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(dpath, 'lib_symmetry*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])

_symmetry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), libfile))
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

HIDE_ISOLDE = 0x02


BACKBONE_MODE_RIBBON=0
BACKBONE_MODE_CA_TRACE=1

_backbone_mode_descr = {
    BACKBONE_MODE_RIBBON: "ribbon",
    BACKBONE_MODE_CA_TRACE: "CA trace"
}

def _format_sym_tuple(result):
    from chimerax.atomic import ctypes_support as convert
    from chimerax.core.geometry import Places
    primary_atoms = convert.atoms(result[0])
    sym_atoms = convert.atoms(result[1])
    n_sym_atoms = len(sym_atoms)
    sym_coords = result[2].reshape((n_sym_atoms,3))
    atom_syms = result[3]
    sym_bonds = convert.bonds(result[4])
    nbonds = len(sym_bonds)
    bond_positions = Places(opengl_array=result[5].reshape((nbonds*2,4,4)))
    bond_syms = result[6]
    return (primary_atoms, sym_atoms, sym_coords, atom_syms, sym_bonds, bond_positions, bond_syms)


def sym_transforms_in_sphere(atoms, transforms, center, cutoff, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double,
            ctypes.c_bool),
        ret = ctypes.py_object)
    natoms = len(atoms)
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(atoms._c_pointers, natoms, pointer(tf), n_tf, pointer(c), cutoff, visible_only)
    return _format_sym_tuple(result)

def whole_residue_sym_sphere(residues, transforms, center, cutoff, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms_by_residue',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double,
            ctypes.c_bool),
            ret = ctypes.py_object)
    nres = len(residues)
    tf = numpy.empty(transforms.shape,numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(residues._c_pointers, nres, pointer(tf), n_tf, pointer(c), cutoff, visible_only)
    return _format_sym_tuple(result)

def atom_and_bond_sym_transforms_from_sym_atoms(atoms, symops, sym_indices, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms_from_sym_atoms',
        args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.c_bool),
        ret=ctypes.py_object)
    natoms = len(atoms)
    nsym = len(sym_indices)
    tf = numpy.empty(symops.shape, numpy.double)
    tf[:] = symops
    si = numpy.empty(len(sym_indices), numpy.uint8)
    si[:] = sym_indices
    result = f(atoms._c_pointers, pointer(sym_indices), natoms, pointer(tf), nsym, visible_only)
    return _format_sym_tuple(result)


def sym_ribbons_in_sphere(tether_coords, transforms, center, cutoff):
    f = c_function('close_sym_ribbon_transforms',
        args=(ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
              ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
              ctypes.POINTER(ctypes.c_double), ctypes.c_double),
        ret=ctypes.py_object)
    tc = numpy.empty(tether_coords.shape, numpy.double)
    tc[:] = tether_coords
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    result = f(pointer(tc), len(tc), pointer(tf), len(tf), pointer(c), cutoff)
    return result

def symmetry_coords(atoms, sym_matrices, sym_indices):
    unique_indices = numpy.unique(sym_indices)
    coords = atoms.coords
    from chimerax.core.geometry import Place
    for i in unique_indices:
        mask = sym_indices == i
        tf = Place(matrix=sym_matrices[i])
        coords[mask] = tf.moved(coords[mask])
    return coords

from chimerax.core.models import Model, Drawing

def get_symmetry_handler(structure, create=False):
    for m in structure.child_models():
        if isinstance(m, XtalSymmetryHandler):
            return m
    if create:
        return XtalSymmetryHandler(structure)
    return None

class XtalSymmetryHandler(Model):
    '''
    Handles crystallographic symmetry and maps for an atomic model.
    '''
    def __init__(self, model, mtzfile=None, calculate_maps=True, map_oversampling=1.5,
        min_voxel_size = 0.5, spotlight_mode = True, map_scrolling_radius=12,
        atomic_symmetry_radius=15, hydrogens='polar'):
        name = 'Crystal'
        session = self.session = model.session
        super().__init__(name, session)

        self.structure = model
        self._box_center = session.view.center_of_rotation
        self._box_center_grid = None
        self._atomic_sym_radius = atomic_symmetry_radius
        self._stepper = None

        from chimerax.core.triggerset import TriggerSet
        trig = self.triggers = TriggerSet()

        trigger_names = (
            'box center moved',       # Centre of rotation moved to a new grid point

            'map box changed',  # Changed shape of box for map viewing
            'map box moved',    # Changed location of box for map viewing
            'atom box changed', # Changed shape or centre of box for showing symmetry atoms
            'backbone mode changed', # Changed backbone mode from ribbon to CA trace or vice versa
        )
        for t in trigger_names:
            trig.add_trigger(t)

        self.mtzdata=None
        self.xmapset = None

        map_data = None

        if mtzfile is not None:
            from .clipper_mtz import ReflectionDataContainer
            mtzdata = self.mtzdata = ReflectionDataContainer(self.session, mtzfile, shannon_rate = map_oversampling)
            self.add([mtzdata])
            self.cell = mtzdata.cell
            self.spacegroup = mtzdata.spacegroup
            self.grid = mtzdata.grid_sampling
            self.hklinfo = mtzdata.hklinfo

            if len(mtzdata.calculated_data):
                datasets = mtzdata.calculated_data
        else:
            from .crystal import symmetry_from_model_metadata
            self.cell, self.spacegroup, self.grid = symmetry_from_model_metadata(model)

        from .crystal import XmapSet
        xmapset = self.xmapset = XmapSet(session, self,
            datasets = mtzdata.calculated_data,
            live_scrolling = spotlight_mode,
            display_radius = map_scrolling_radius)
        xmapset.pickable = False
        self.add([xmapset])

        cell = self.cell
        spacegroup = self.spacegroup
        grid = self.grid

        self._voxel_size = cell.dim/grid.dim

        ref = model.bounds().center()

        from .main import atom_list_from_sel
        ca = self._clipper_atoms = atom_list_from_sel(model.atoms)

        uc = self._unit_cell = clipper.Unit_Cell(ref, ca, cell, spacegroup, grid)


        self._atomic_symmetry_model = AtomicSymmetryModel(model, self, uc,
            radius = atomic_symmetry_radius, live = spotlight_mode)

        self._update_handler = session.triggers.add_handler('new frame',
            self.update)

        model.add([self])

        display_atoms = model.atoms
        display_atoms.draw_modes = display_atoms.STICK_STYLE
        display_atoms.displays = False
        if hydrogens == 'polar':
            display_atoms = display_atoms[display_atoms.idatm_types != 'HC']
        elif hydrogens == 'none':
            display_atoms = display_atoms[display_atoms.element_names != 'H']
        display_atoms.displays = True


        from .mousemodes import initialize_map_contour_mouse_modes
        initialize_map_contour_mouse_modes(session)
        self.spotlight_mode = spotlight_mode

    @property
    def stepper(self):
        '''
        Provides methods for "stepping" back and forth through the
        model according to secondary structure. For example, each call
        to stepper.step_forward() (with default arguments) will return
        an atom selection corresponding to the next pair of defined
        secondary structure elements plus their flanking loops.
        '''
        if self._stepper is None:
            from .structurestepper import StructureStepper
            self._stepper = StructureStepper(self.session, self.structure)
        return self._stepper


    @property
    def spotlight_mode(self):
        return self._spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, flag):
        from chimerax.std_commands import cofr
        cofr.cofr(self.session, 'centerOfView', show_pivot=True)
        self._atomic_symmetry_model.live_scrolling = flag
        if self.xmapset is not None:
            self.xmapset.live_scrolling = flag
        self._spotlight_mode = flag
        self.update(force=True)

    @property
    def atomic_symmetry_model(self):
        return self._atomic_symmetry_model

    @property
    def atomic_sym_radius(self):
        return self._atomic_symmetry_model.radius

    @atomic_sym_radius.setter
    def atomic_sym_radius(self, radius):
        self._atomic_symmetry_model.radius = radius

    @property
    def map_view_radius(self):
        return self.xmapset.display_radius

    @map_view_radius.setter
    def map_view_radius(self, radius):
        self.xmapset.display_radius = radius

    @property
    def unit_cell(self):
        return self._unit_cell

    def delete(self):
        if self._update_handler is not None:
            self.session.triggers.remove_handler(self._update_handler)
            self._update_handler = None
        super().delete()

    def update(self, *_, force=False):
        v = self.session.view
        cofr = self._box_center = v.center_of_rotation
        cofr_grid = clipper.Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
        if force:
            update_needed=True
        else:
            update_needed = False
            if self._box_center_grid is None:
                update_needed = True
            elif (cofr_grid != self._box_center_grid):
                update_needed = True
        if update_needed:
            self.triggers.activate_trigger('box center moved', (cofr, cofr_grid))
            self._box_center_grid = cofr_grid

    @property
    def atom_sym(self):
        return self._atomic_symmetry_model

    def isolate_and_cover_selection(self, atoms, include_surrounding_residues=5,
        show_context=5, mask_radius=3, extra_padding=0, hide_surrounds=True,
        focus=True, include_symmetry = True):
        '''
        Expand the map(s) (if present) to cover a given atomic selection,  then
        mask them to within a given distance of said atoms to reduce visual
        clutter. Adjust the atomic visualisation to show only the selected
        atoms, plus an optional surrounding buffer zone.
        Args:
          atoms (ChimeraX Atoms object):
            The main selection we're interested in. The existing selection will
            be expanded to include the whole residue for every selected atom.
          include_surrounding_residues (float):
            Any residue with an atom coming within this radius of the primary
            selection will be added to the selection covered by the map. To
            cover only the primary selection, set this value to zero.
          show_context (float):
            Any residue within an atom coming within this radius of the previous
            two selections will be displayed as a thinner stick representation,
            but will not be considered for the map masking calculation.
          mask_radius (float):
            Components of the map more than this distance from any atom will
            be hidden.
          extra_padding (float):
            Optionally, further pad the volume by this distance. The extra
            volume will be hidden, but available for calculations.
          hide_surrounds (bool):
            If true, all residues outside the selection region will be hidden
          focus (bool):
            If true, the camera will be moved to focus on the selection (only
            the atoms in the master model will be considered)
          include_symmetry (bool):
            If true, symmetry atoms will be considered for the purposes of
            show_context.
        '''
        self.spotlight_mode = False
        original_atoms = atoms
        atoms = original_atoms.unique_residues.atoms
        asm = self.atomic_symmetry_model
        maps = self.xmapset
        main_set = asm.sym_select_within(atoms, include_surrounding_residues)
        main_coords = symmetry_coords(*main_set[0:3])
        context_set = asm.sym_select_within(main_set[0], show_context,
            coords=main_coords)
        asm._current_focal_set = context_set
        asm.set_sym_display(context_set[3], *atom_and_bond_sym_transforms_from_sym_atoms(*context_set[0:3]))

        cell = self.cell
        grid = self.grid
        xmaps = self.xmapset
        from .crystal import calculate_grid_padding
        pad = calculate_grid_padding(mask_radius, grid, cell)
        ep = calculate_grid_padding(extra_padding, grid, cell)
        box_bounds_grid = clipper.Util.get_minmax_grid(main_coords, cell, grid) \
            + numpy.array((-pad, pad)) + numpy.array((-ep, ep))
        xmaps.set_box_limits(box_bounds_grid)
        xmaps._surface_zone.update(mask_radius, coords = main_coords)
        xmaps._reapply_zone()
        if focus:
            focus_on_selection(self.session, self.session.view, atoms)

    _cube_pairs = numpy.array([[0,1], [0,2], [0,4], [1,3], [1,5], [2,3], [2,6], [3,7], [4,5], [4,6], [5,7], [6,7]], numpy.int)

    def draw_unit_cell_box(self, offset=None, cylinder_radius = 0.05):
        '''
        Draw a rhombohedral box around one unit cell.
        '''
        m, d = _get_special_positions_model(self)
        model = self.structure
        uc = self.unit_cell

        if offset is None:
            offset = numpy.zeros(3)

        from chimerax.core.geometry import Place, Places
        positions = []
        colors = []
        rgba_edge = numpy.array([0,255,255,128],numpy.uint8)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset\
                        + uc.min.coord_frac(self.grid).uvw

        corners = numpy.array([clipper.Coord_frac(c).coord_orth(self.cell).xyz for c in corners_frac])
        from chimerax.surface.shapes import cylinder_geometry
        d.set_geometry(*cylinder_geometry())
        d.set_color(rgba_edge)
        from chimerax.core.geometry import cylinder_rotations, Places
        xyz0 = numpy.array([corners[p[0]] for p in self._cube_pairs])
        xyz1 = numpy.array([corners[p[1]] for p in self._cube_pairs])
        radii = numpy.ones(12)*cylinder_radius
        p = numpy.empty((12, 4, 4), numpy.float32)
        cylinder_rotations(xyz0, xyz1, radii, p)
        p[:,3,:3] = 0.5*(xyz0+xyz1)
        pl = Places(opengl_array = p)
        d.set_positions(pl)
        m.display = True




    def draw_unit_cell_and_special_positions(self, offset = None):
        '''
        Quick-and-dirty drawing mapping out the special positions
        (positions which map back to themselves by at least one
        non-unity symop) within one unit cell. A sphere will be drawn
        at each grid-point with non-unit multiplicity, and colour-coded
        according to multiplicity:
            2-fold: white
            3-fold: cyan
            4-fold: yellow
            6-fold: magenta

        Ultimately it would be nice to replace this with something more
        elegant, that masks and scrolls continuously along with the model/
        map visualisation.

        Args:
            offset (1x3 numpy array, default = None):
                Optional (u,v,w) offset (in fractions of a unit cell axis)
        '''
        m, d = _get_special_positions_model(self)
        model = self.structure

        ref = model.bounds().center().astype(float)
        frac_coords = clipper.Coord_orth(ref).coord_frac(self.cell).uvw
        if offset is None:
            offset = numpy.array([0,0,0],int)

        xmap = self.xmapset[0].xmap
        uc = self.unit_cell
        spc = numpy.array(xmap.special_positions_unit_cell_xyz(uc, offset))
        from chimerax.surface.shapes import sphere_geometry2
        sphere = numpy.array(sphere_geometry2(80))
        sphere[0]*=0.25
        d.set_geometry(*sphere)
        #d.vertices, d.normals, d.triangles = sphere

        from chimerax.core.geometry import Place, Places
        positions = []
        colors = []
        rgba_corner = numpy.array([255,0,255,128],numpy.uint8)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset\
                        + self.unit_cell.min.coord_frac(self.grid).uvw

        corners = []
        for c in corners_frac:
            co = clipper.Coord_frac(c).coord_orth(self.cell)
            positions.append(Place(axes=numpy.identity(3)*4, origin=co.xyz))
            colors.append(rgba_corner)

        if len(spc):

            coords = spc[:,0:3]
            multiplicity = spc[:,3].astype(int)
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.uint8)
            rgba_3fold = numpy.array([0,255,255,255],numpy.uint8)
            rgba_4fold = numpy.array([255,255,0,255],numpy.uint8)
            rgba_6fold = numpy.array([255,0,0,255],numpy.uint8)

            for coord, mult in zip(coords, multiplicity):
                if mult == 2:
                    positions.append(Place(axes=scale_2fold, origin=coord))
                    colors.append(rgba_2fold)
                elif mult == 3:
                    positions.append(Place(axes=scale_3fold, origin=coord))
                    colors.append(rgba_3fold)
                elif mult == 4:
                    positions.append(Place(axes=scale_4fold, origin=coord))
                    colors.append(rgba_4fold)
                elif mult == 6:
                    positions.append(Place(axes=scale_6fold, origin=coord))
                    colors.append(rgba_6fold)
            for c in corners:
                positions.append(Place(axes=scale_6fold, origin=c))
                colors.append(rgba_corner)

        d.set_positions(Places(positions))
        d.set_colors(numpy.array(colors, numpy.uint8))
        # m.add_drawing(d)
        # model.parent.add([m])
        m.display = True


def _get_special_positions_model(m):
    for cm in m.child_models():
        if cm.name == "Special positions":
            return cm, cm.child_drawings()[0]
    from chimerax.core.models import Model, Drawing
    d = Drawing("Special positions")
    sm = Model("Special positions", m.session)
    sm.add_drawing(d)
    m.add([sm])
    return sm, d


class AtomicSymmetryModel(Model):
    '''
    Finds and draws local symmetry atoms for an atomic structure
    '''
    def __init__(self, atomic_structure, parent, unit_cell, radius = 15,
        dim_colors_to = 0.6, backbone_mode = 'CA trace', live = True):
        self._live_scrolling = False
        self.structure = atomic_structure
        session = self.session = atomic_structure.session
        self.unit_cell = unit_cell
        self.cell = parent.cell
        self.spacegroup =parent.spacegroup
        self.grid = parent.grid
        self.session = atomic_structure.session
        self._color_dim_factor = dim_colors_to

        self._current_focal_set = None

        self.manager = parent

        super().__init__('Atomic symmetry', session)
        parent.add([self])

        self._last_hides = atomic_structure.atoms.hides

        self._box_changed_handler = None
        self._center = session.view.center_of_rotation
        self._radius = radius
        self._box_dim = numpy.array([radius*2, radius*2, radius*2], numpy.double)
        ad = self._atoms_drawing = SymAtomsDrawing('Symmetry atoms')
        self.add_drawing(ad)
        #from chimerax.core.atomic.structure import PickedBonds
        bd = self._bonds_drawing = SymBondsDrawing('Symmetry bonds', PickedSymBond, structure.PickedBonds)
        self.add_drawing(bd)
        rd = self._ribbon_drawing = SymRibbonDrawing('Symmetry ribbons',
            atomic_structure._ribbon_drawing, dim_colors_to)
        self.add_drawing(rd)
        self._current_atoms = atomic_structure.atoms
        self._model_changes_handler = atomic_structure.triggers.add_handler(
                                        'changes', self._model_changed_cb)
        self.live_scrolling = live
        from .crystal import set_to_default_cartoon
        set_to_default_cartoon(session)

    def sym_select_within(self, atoms, cutoff, coords=None, whole_residues = True):
        '''
        Given a set of atoms, return a (atoms, symops, sym_indices) tuple
        giving all atoms and their symmetry operators within the given cutoff
        distance from any atom in the primary set.
        Args:
            atoms:
                The core set of atoms to be surrounded by the new selection.
            cutoff:
                Cutoff distance in Angstroms.
            coords (default None):
                Optionally, you can provide the coordinates for all atoms
                (useful for expanding a selection that already includes
                symmetry atoms).
            whole_residues (default true):
                Whether to expand the selections to whole_residues.
        '''
        if coords is None:
            coords = atoms.coords
        coords = coords.astype(numpy.float32)
        master_atoms = self.structure.atoms
        master_coords = master_atoms.coords.astype(numpy.float32)
        grid_minmax = clipper.Util.get_minmax_grid(coords, self.cell, self.grid)
        from .crystal import calculate_grid_padding
        pad = calculate_grid_padding(cutoff, self.grid, self.cell)
        grid_minmax += numpy.array((-pad, pad))
        min_xyz = clipper.Coord_grid(grid_minmax[0]).coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = grid_minmax[1]-grid_minmax[0]
        symops = self.unit_cell.all_symops_in_box(min_xyz, dim, True)
        symmats = symops.all_matrices_orth(self.cell, format='3x4')
        from chimerax.core.geometry import Place
        target = [(coords, Place().matrix.astype(numpy.float32))]
        search_list = []
        for i, s in enumerate(symmats):
            search_list.append((master_coords, s.astype(numpy.float32)))
        from chimerax.core.geometry import find_close_points_sets
        # We want to respond gracefully if the cutoff is zero, but
        # find_close_points_sets returns nothing if the cutoff is
        # *actually* zero. So just set it to a very small non-zero value.
        if cutoff == 0:
            cutoff = 1e-6
        i1, i2 = find_close_points_sets(search_list, target, cutoff)
        found_atoms = []
        sym_indices = []
        sym_count = 0
        for i, (c, s) in enumerate(zip(i1, symmats)):
            if len(c):
                sel = master_atoms[c]
                if whole_residues:
                    sel = sel.unique_residues.atoms
                found_atoms.append(sel)
                indices = numpy.empty(len(sel), numpy.uint8)
                indices[:] = i
                sym_indices.append(indices)
        if len(found_atoms) > 1:
            from chimerax.atomic import concatenate
            found_atoms = concatenate(found_atoms)
            sym_indices = numpy.concatenate(sym_indices)
        else:
            found_atoms = found_atoms[0]
            sym_indices = sym_indices[0]
        return (found_atoms, symmats, sym_indices, symops)

    def delete(self):
        bh = self._box_changed_handler
        if bh is not None:
            self.manager.triggers.remove_handler(bh)
        mh = self._model_changes_handler
        if mh is not None:
            self.structure.triggers.remove_handler(mh)
        self.unhide_all_atoms()
        super().delete()

    @property
    def backbone_mode(self):
        return _backbone_mode_descr[self._backbone_mode]

    @backbone_mode.setter
    def backbone_mode(self, mode):
        old_mode = self.backbone_mode
        if mode == "ribbon":
            self._backbone_mode = BACKBONE_MODE_RIBBON
        elif mode == "CA trace":
            self._backbone_mode = BACKBONE_MODE_CA_TRACE
        else:
            raise TypeError('Unrecognised mode! Should be one of "ribbon" or "CA trace"')
        if old_mode != mode:
            self.triggers.activate_trigger('backbone mode changed', mode)

    def unhide_all_atoms(self):
        self.structure.atoms.hides &= ~HIDE_ISOLDE

    @property
    def dim_factor(self):
        return self._color_dim_factor

    @dim_factor.setter
    def dim_factor(self, factor):
        self._color_dim_factor = factor
        self._ribbon_drawing.dim_factor = factor

    @property
    def live_scrolling(self):
        return self._live_scrolling

    @live_scrolling.setter
    def live_scrolling(self, flag):
        bh = self._box_changed_handler
        if flag and not self._live_scrolling:
            if bh is None:
                self._box_changed_handler = self.parent.triggers.add_handler('box center moved',
                    self._box_moved_cb)
                self._update_box()
                self.update_graphics()
        elif not flag and self._live_scrolling:
            from chimerax.atomic import concatenate
            res = whole_residue_sym_sphere(self.structure.residues, self._current_tfs, self._center, self._radius, visible_only=False)
            ma = res[0]
            sa = res[1]
            fa = concatenate((ma, sa))
            fs = numpy.concatenate((numpy.zeros(len(ma), numpy.uint8), res[3]))
            self._current_focal_set = (fa, self._current_tfs, fs, self._current_symops)

            if bh is not None:
                self.parent.triggers.remove_handler(bh)
                self._box_changed_handler = None
        self._live_scrolling = flag

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius):
        self._radius = radius
        if self.visible:
            self.update_graphics()

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, flag):
        if not flag:
            self.unhide_all_atoms()
        super().set_display(flag)
        if flag:
            if self.live_scrolling:
                self._center = self.session.view.center_of_rotation
                self._update_box()
            self.update_graphics()

    @property
    def _level_of_detail(self):
        from chimerax.atomic import structure_graphics_updater
        gu = structure_graphics_updater(self.session)
        return gu.level_of_detail


    def _box_moved_cb(self, trigger_name, box_params):
        if not self.visible:
            return
        self._center = box_params[0]
        self._update_box()
        self.update_graphics()

    def _update_box(self):
        from .crystal import find_box_corner
        center = self._center
        radius = self._radius
        box_corner_grid, box_corner_xyz = find_box_corner(center, self.cell, self.grid, radius)
        dim = self._box_dim
        dim[:] = radius*2
        grid_dim = (dim / self.parent._voxel_size).astype(numpy.int32)
        symops = self._current_symops = self.unit_cell.all_symops_in_box(box_corner_xyz, grid_dim, True)
        # Identity symop will always be the first in the list
        tfs = self._current_tfs = symops.all_matrices_orth(self.cell, format='3x4')
        atoms = self.structure.atoms
        atoms.hides |=HIDE_ISOLDE
        self._current_master_atoms, self._current_sym_atoms, self._current_sym_atom_coords, \
        self._current_atom_syms, self._current_bonds, self._current_bond_tfs, \
        self._current_bond_syms = \
                whole_residue_sym_sphere(self.structure.residues, tfs, center, radius)
                #sym_transforms_in_sphere(atoms, tfs[first_symop:], center, radius)
        self._current_master_atoms.hides &= ~HIDE_ISOLDE
        cs = self._current_sym_atoms
        if len(cs):
            crs = self._current_ribbon_syms = numpy.unique(self._current_atom_syms)
            self._current_ribbon_syms = crs[crs!=0]
        else:
            self._current_ribbon_syms = sym_ribbons_in_sphere(self._ribbon_drawing._tether_coords, tfs, center, radius)


    def set_sym_display(self, symops, primary_atoms, sym_atoms, sym_coords, atom_sym_indices,
        sym_bonds, sym_bond_tfs, bond_sym_indices):
        '''
        Manually define the set of symmetry atoms to be displayed.
        Args:
            primary_atoms:
                A ChimeraX Atoms object defining the atoms in the master model
                to be displayed.
            symops:
                A Clipper Symops object defining the symmetry transforms
            sym_atoms:
                A ChimeraX Atoms object listing the symmetry atoms to be displayed
            sym_coords:
                The symmetry coordinates corresponding to sym_atoms
            atom_sym_indices:
                A Numpy integer array of length equal to sym_atoms, giving the
                index of the corresponding symmetry operator in symops
            sym_bonds:
                A ChimeraX Bonds object listing the symmetry bonds to be
                displayed.
            sym_bond_tfs:
                A Places object providing the halfbond cylinder transforms to
                draw the sym_bonds in their correct positions
            bond_sym_indices:
                A Numpy integer array of length equal to sym_bonds, giving the
                index of the corresponding symmetry operator in symops
        '''
        self.live_scrolling = False
        all_atoms = self.structure.atoms
        all_atoms.hides |= HIDE_ISOLDE
        primary_atoms.hides &= ~HIDE_ISOLDE
        self._current_symops = symops
        self._current_tfs = symops.all_matrices_orth(self.cell, format='3x4')
        self._current_master_atoms = primary_atoms
        self._current_sym_atoms = sym_atoms
        self._current_sym_atom_coords = sym_coords
        csym = self._current_atom_syms = atom_sym_indices
        self._current_bonds = sym_bonds
        self._current_bond_tfs = sym_bond_tfs
        self._current_bond_syms = bond_sym_indices
        self._current_ribbon_syms = numpy.unique(csym[csym!=0])
        self.update_graphics()

    def _model_changed_cb(self, trigger_name, changes):
        if not self.visible:
            return
        changes = changes[1]
        update_needed = False
        ribbon_update_needed = False
        if len(changes.created_atoms()):
            update_needed = True
            ribbon_update_needed = True
        if changes.num_deleted_atoms() > 0:
            update_needed = True
            ribbon_update_needed = True
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if 'display changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if  'hide changed' in reasons:
            hides = self.structure.atoms.hides
            # Prevent repeated callback with every display update in spotlight mode
            current_hides = hides&~HIDE_ISOLDE
            if numpy.any(current_hides != self._last_hides):
                update_needed = True
            self._last_hides = current_hides
        if 'color changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if 'ribbon_display changed' in changes.residue_reasons():
            ribbon_update_needed = True
        if (ribbon_update_needed):
            self._ribbon_drawing.delayed_rebuild(self.session)
        if (update_needed):
            if self.live_scrolling:
                self._update_box()
                self.update_graphics()
            else:
                self._update_sym_coords()

    def _update_sym_coords(self):
        focal_set = self._current_focal_set
        # try:
        res = atom_and_bond_sym_transforms_from_sym_atoms(*focal_set[0:3])
        self.set_sym_display(focal_set[3], *res)
        # except:
        #     from chimerax.atomic import Atoms, Bonds
        #     self._current_sym_atoms = Atoms()
        #     self._current_bonds = Bonds()
        #     self.update_graphics()


    def update_graphics(self):
        lod = self._level_of_detail
        self._update_atom_graphics(lod)
        self._update_bond_graphics(lod)
        self._update_ribbon_graphics()

    def _update_atom_graphics(self, lod):
        ad = self._atoms_drawing
        ca = self._current_sym_atoms
        if not len(ca):
            ad.display=False
            return
        syms = self._current_atom_syms
        coords = self._current_sym_atom_coords
        # if not self.live_scrolling:
        #     mask = _atoms_only_hidden_by_clipper(ca)
        #     ca = ca[mask]
        #     syms = syms[mask]
        #     coords = coords[mask]
        ad.visible_atoms = ca
        lod.set_atom_sphere_geometry(ad)

        na = len(ca)
        if na > 0:
            ad.display = True
            xyzr = numpy.empty((na, 4), numpy.float32)
            xyzr[:,:3] = coords
            xyzr[:,3] = self.structure._atom_display_radii(ca)
            from chimerax.core.geometry import Places
            ad.positions = Places(shift_and_scale = xyzr)
            colors = ca.colors.astype(numpy.float32)
            colors[:,:3] *= self._color_dim_factor
            ad.colors = colors.astype(numpy.uint8)
        else:
            ad.display = False


    def _update_bond_graphics(self, lod):
            bd = self._bonds_drawing
            bonds = self._current_bonds
            if not len(bonds):
                bd.display = False
                return
            bsym = self._current_bond_syms
            b_tfs = self._current_bond_tfs
            # if not self.live_scrolling:
            #     mask = _bonds_only_hidden_by_clipper(bonds)
            #     bonds = bonds[mask]
            #     mask = numpy.concatenate((mask, mask))
            #     bsym = bsym[mask]
            #     from chimerax.core.geometry import Places
            #     b_tfs = Places(place_array=b_tfs.array()[mask])
            lod.set_bond_cylinder_geometry(bd)
            bd.visible_bonds = bonds
            nb = len(bonds)
            if nb > 0:
                bd.display = True
                bd.positions = b_tfs
                colors = bonds.half_colors.astype(numpy.float32)
                colors[:,:3] *= self._color_dim_factor
                bd.colors = colors.astype(numpy.uint8)
            else:
                bd.display = False

    def _update_ribbon_graphics(self):
        rd = self._ribbon_drawing
        prd = self.structure._ribbon_drawing
        position_indices = self._current_ribbon_syms
        tfs = self._current_tfs[position_indices]
        from chimerax.core.geometry import Places
        rd.positions = Places(place_array=tfs)
        rd.display = prd.display
        if not rd.display:
            return
        rd.update()

#from chimerax.core.atomic.structure import AtomsDrawing
class SymAtomsDrawing(structure.AtomsDrawing):
    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        if not self.display or self.visible_atoms is None or (exclude and exclude(self)):
            return None
        xyzr = self.positions.shift_and_scale_array()
        coords, radii = xyzr[:,:3], xyzr[:,3]

        from chimerax.core.geometry import closest_sphere_intercept
        fmin, anum = closest_sphere_intercept(coords, radii, mxyz1, mxyz2)
        if fmin is None:
            return None
        atom = self.visible_atoms[anum]
        atom_syms = self.parent._current_atom_syms
        sym = self.parent._current_symops[int(atom_syms[anum])]
        return PickedSymAtom(atom, fmin, sym)

    def planes_pick(self, planes, exclude=None):
        return []

    def bounds(self, positions=True):
        if not positions:
            return self._geometry_bounds()
        cpb = self._cached_position_bounds
        if cpb is not None:
            return cpb
        xyzr = self.positions.shift_and_scale_array()
        if xyzr is None:
            return self._geometry_bounds()
        coords, radii = xyzr[:, :3], xyzr[:,3]
        from chimerax.core.geometry import sphere_bounds
        b = sphere_bounds(coords, radii)
        self._cached_position_bounds = b

    def update_selection(self):
        pass

from chimerax.core.graphics import Pick

class PickedSymAtom(Pick):
    def __init__(self, atom, distance, sym):
        super().__init__(distance)
        self.atom = atom
        self.sym = sym

    def description(self):
        return '({}) {}'.format(self.sym, self.atom)

    def select(self, mode = 'add'):
        pass

#from chimerax.core.atomic.structure import BondsDrawing
class SymBondsDrawing(structure.BondsDrawing):
    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        return None #too-hard basket for now.

        if not self.display or (exclude and exclude(self)):
            return None
        #from chimerax.core.atomic.structure import _bond_intercept
        b, f = structure._bond_intercept(bonds, mxyz1, mxyz2)

    def planes_pick(self, mxyz1, mxyz2, exclude=None):
        return []

    def bounds(self, positions=True):
        #TODO: properly calculate bounds
        return self._geometry_bounds()

    def select(self, mode = 'add'):
        pass

    def update_selection(self):
        pass

class PickedSymBond(Pick):
    def __init__(self, bond, distance, sym):
        super().__init__(distance)
        self.bond = bond
        self.sym = sym
    def description(self):
        return '({}) {}'.format(self.sym, self.bond)

class SymRibbonDrawing(Drawing):
    pickable = False
    def __init__(self, name, master_ribbon, dim_factor):
        super().__init__(name)
        m = self._master = master_ribbon
        self._tether_coords = numpy.array([], numpy.double)
        _copy_ribbon_drawing(m, self, dim_factor)
        self._dim_factor = dim_factor

    @property
    def dim_factor(self):
        return self._dim_factor

    @dim_factor.setter
    def dim_factor(self, factor):
        self._dim_factor = factor
        self.update

    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        return None

    def planes_pick(self, myxz1, mxyz2, exclude=None):
        return []

    def delayed_rebuild(self, session):
        self._rebuild_handler = session.triggers.add_handler('frame drawn', self.rebuild)

    def rebuild(self, *args):
        self.remove_all_drawings()
        _copy_ribbon_drawing(self._master, self, self.dim_factor)
        self.find_tether_coords()
        if args and args[0] == "frame drawn":
            self._rebuild_handler = None
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def find_tether_coords(self):
        tethers = []
        for d in self.child_drawings():
            if 'ribbon_tethers' in d.name:
                tethers.append(d.positions.array()[:,:,3])
        if len(tethers):
            self._tether_coords = numpy.concatenate(tethers)
        else:
            self._tether_coords = numpy.array([], numpy.double)



    def update(self):
        dim_factor = self.dim_factor
        sad, mad = self.all_drawings(), self._master.all_drawings()
        if len(sad) != len(mad):
            self.rebuild()
            return
        for td, md in zip(sad, mad):
            td.set_geometry(md.vertices, md.normals, md.triangles)
            td.display = md.display
            vertex_colors = md.vertex_colors
            if vertex_colors is not None:
                vertex_colors = vertex_colors.astype(numpy.float32)
                vertex_colors[:,:3] *= dim_factor
                td.vertex_colors = vertex_colors.astype(numpy.uint8)
            else:
                colors = md.colors.astype(numpy.float32)
                colors[:,:3] *= dim_factor
                td.colors = colors
        self.find_tether_coords()


def _atoms_only_hidden_by_clipper(atoms):
    hides = atoms.hides
    return numpy.logical_and(atoms.displays,
        numpy.logical_not(numpy.logical_and(hides&HIDE_ISOLDE, hides&~HIDE_ISOLDE)))

def _bonds_only_hidden_by_clipper(bonds):
    atoms = bonds.atoms
    return numpy.logical_and(bonds.displays, numpy.logical_and(
        _atoms_only_hidden_by_clipper(atoms[0]),
        _atoms_only_hidden_by_clipper(atoms[1])
    ))


def focus_on_selection(session, view, atoms, clip = True):
    v = view
    pad = 5.0
    bounds = atoms.scene_bounds
    bounds.xyz_min = bounds.xyz_min - pad
    bounds.xyz_max = bounds.xyz_max + pad
    radius = bounds.radius() + pad
    cofr_method = v.center_of_rotation_method
    v.view_all(bounds)
    v.center_of_rotation = center = bounds.center()
    v.center_of_rotation_method = cofr_method
    cam = v.camera
    vd = cam.view_direction()
    if clip:
        cp = v.clip_planes
        cp.set_clip_position('near', center - radius*vd, cam)
        cp.set_clip_position('far', center + radius*vd, cam)
    session.selection.clear()
    atoms.selected=True


def _copy_ribbon_drawing(master_drawing, target_drawing, dim_factor):
    from chimerax.core.graphics import Drawing
    d = master_drawing
    t = target_drawing
    def recursively_add_drawings(fromd, tod, dim_factor):
        children = fromd.child_drawings()
        for c in children:
            nc = Drawing(c.name)
            v, n, t, p, d = c.vertices, c.normals, c.triangles, c.positions, c.display
            if v is None or len(v) == 0:
                continue
            nc.set_geometry(v, n, t)
            nc.positions, nc.display = p, d
            vc = c.vertex_colors
            if vc is not None:
                vc = vc.astype(numpy.float32)
                vc[:,:3] *= dim_factor
                nc.vertex_colors = vc.astype(numpy.uint8)
            else:
                col = c.colors.astype(numpy.float32)
                col[:, :3] *= dim_factor
                nc.colors = col.astype(numpy.uint8)
            tod.add_drawing(nc)
            recursively_add_drawings(c, nc, dim_factor)
    recursively_add_drawings(d, t, dim_factor)

def _get_ca_pbg(m):
    from chimerax.atomic import PseudobondGroup
    for cm in m.child_models():
        if isinstance(cm, PseudobondGroup) and m.name == "CA trace":
            return cm
    pbg = m.pseudobond_group("CA trace")
    pbg.dashes = 1
    return pbg

def create_ca_trace(m):
    pbg = _get_ca_pbg(m)
    pbg.clear()
    chains = m.polymers(missing_structure_treatment = m.PMS_NEVER_CONNECTS)
    for chain in chains:
        chain = chain[0]
        if not chain[0].polymer_type == Residue.PT_AMINO:
            continue
        for i in range(len(chain)-1):
            r1, r2 = chain[i:i+2]
            try:
                ca1 = r1.atoms[r1.atoms.names=='CA'][0]
                ca2 = r2.atoms[r2.atoms.names=='CA'][0]
                pb = pbg.new_pseudobond(ca1, ca2)
                pb.color = ca1.color
            except:
                continue
