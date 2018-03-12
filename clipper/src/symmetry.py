import numpy
import sys, os, glob
import ctypes

from . import clipper

from chimerax.core.atomic import molc

from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t

from chimerax.core.atomic.molobject import _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map


dpath = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(dpath, '_symmetry.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])

_symmetry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), libfile))
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

HIDE_ISOLDE = 0x02

def sym_transforms_in_sphere(atoms, transforms, center, cutoff):
    f = c_function('atom_and_bond_sym_transforms',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double),
        ret = ctypes.py_object)
    natoms = len(atoms)
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(atoms._c_pointers, natoms, pointer(tf), n_tf, pointer(c), cutoff)
    from chimerax.core.atomic.molarray import _atoms, _bonds
    from chimerax.core.geometry import Places
    natoms = len(result[0])
    atom_coords = result[1].reshape((natoms,3))
    nbonds = len(result[3])
    bond_positions = Places(opengl_array=result[4].reshape((nbonds*2,4,4)))
    return (_atoms(result[0]), atom_coords, result[2], _bonds(result[3]), bond_positions, result[5])

def whole_residue_sym_sphere(residues, transforms, center, cutoff):
    f = c_function('atom_and_bond_sym_transforms_by_residue',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double),
            ret = ctypes.py_object)
    nres = len(residues)
    tf = numpy.empty(transforms.shape,numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(residues._c_pointers, nres, pointer(tf), n_tf, pointer(c), cutoff)
    from chimerax.core.atomic.molarray import _atoms, _bonds
    from chimerax.core.geometry import Places
    natoms = len(result[0])
    atom_coords = result[1].reshape((natoms,3))
    nbonds = len(result[3])
    bond_positions = Places(opengl_array=result[4].reshape((nbonds*2,4,4)))
    return (_atoms(result[0]), atom_coords, result[2], _bonds(result[3]), bond_positions, result[5])

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
    return f(pointer(tc), len(tc), pointer(tf), len(tf), pointer(c), cutoff)

from chimerax.core.models import Model, Drawing

class XtalSymmetryHandler(Model):
    '''
    Handles crystallographic symmetry for an atomic model.
    '''
    def __init__(self, model, mtzfile=None, calculate_maps=True, map_oversampling=1.5,
        live_map_scrolling=True, map_scrolling_radius=12, live_atomic_symmetry=True,
        atomic_symmetry_radius=15):
        name = 'Crystal'
        session = self.session = model.session
        super().__init__(name, session)

        self.structure = model
        self._box_center = session.view.center_of_rotation
        self._box_center_grid = None
        self._atomic_sym_radius = atomic_symmetry_radius

        from chimerax.core.triggerset import TriggerSet
        trig = self.triggers = TriggerSet()

        trigger_names = (
            'box center moved',       # Centre of rotation moved to a new grid point

            'map box changed',  # Changed shape of box for map viewing
            'map box moved',    # Changed location of box for map viewing
            'atom box changed', # Changed shape or centre of box for showing symmetry atoms
        )
        for t in trigger_names:
            trig.add_trigger(t)

        self.mtzdata=None
        self.xmapset = None

        if mtzfile is not None:
            from .clipper_mtz import ReflectionDataContainer
            mtzdata = self.mtzdata = ReflectionDataContainer(self.session, mtzfile, shannon_rate = map_oversampling)
            self.add([mtzdata])
            self.cell = mtzdata.cell
            self.spacegroup = mtzdata.spacegroup
            self.grid = mtzdata.grid_sampling
            self.hklinfo = mtzdata.hklinfo

            if len(mtzdata.calculated_data):
                from .crystal import XmapSet
                xmapset = self.xmapset = XmapSet(session, mtzdata.calculated_data, self)
                xmapset.pickable = False
                self.add([xmapset])
        else:
            from .crystal import symmetry_from_model_metadata
            self.cell, self.spacegroup, self.grid = symmetry_from_model_metadata(model)

        cell = self.cell
        spacegroup = self.spacegroup
        grid = self.grid

        self._voxel_size = cell.dim/grid.dim

        ref = model.bounds().center()

        from .main import atom_list_from_sel
        ca = self._clipper_atoms = atom_list_from_sel(model.atoms)

        uc = self._unit_cell = clipper.Unit_Cell(ref, ca, cell, spacegroup, grid)


        self._atomic_symmetry_model = AtomicSymmetryModel(model, self, uc,
            radius = atomic_symmetry_radius, live = live_atomic_symmetry)

        self._update_handler = session.triggers.add_handler('new frame',
            self.update)

        model.add([self])


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

    def update(self, *_):
        v = self.session.view
        cofr = self._box_center = v.center_of_rotation
        cofr_grid = clipper.Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
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

class AtomicSymmetryModel(Model):
    '''
    Finds and draws local symmetry atoms for an atomic structure
    '''
    def __init__(self, atomic_structure, parent, unit_cell, radius = 15,
        spotlight_mode = True, dim_colors_to = 0.6, live = True):
        self._live = False
        self.structure = atomic_structure
        session = self.session = atomic_structure.session
        self.unit_cell = unit_cell
        self.cell = parent.cell
        self.spacegroup =parent.spacegroup
        self.grid = parent.grid
        self.session = atomic_structure.session
        self._include_identity = spotlight_mode
        self._spotlight_mode = spotlight_mode
        self._color_dim_factor = dim_colors_to

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
        from chimerax.core.atomic.structure import PickedBonds
        bd = self._bonds_drawing = SymBondsDrawing('Symmetry bonds', PickedSymBond, PickedBonds)
        self.add_drawing(bd)
        rd = self._ribbon_drawing = SymRibbonDrawing('Symmetry ribbons',
            atomic_structure._ribbon_drawing, dim_colors_to)
        self.add_drawing(rd)
        self._current_atoms = atomic_structure.atoms
        self._model_changes_handler = atomic_structure.triggers.add_handler(
                                        'changes', self._model_changed_cb)
        self.live = live
        from .crystal import set_to_default_cartoon
        set_to_default_cartoon(session)

    def delete(self):
        bh = self._box_changed_handler
        if bh is not None:
            self.manager.triggers.remove_handler(bh)
        mh = self._model_changes_handler
        if mh is not None:
            self.structure.triggers.remove_handler(mh)
        self.unhide_all_atoms()
        super().delete()

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
    def live(self):
        return self._live

    @live.setter
    def live(self, flag):
        if flag and not self._live:
            self._box_changed_handler = self.parent.triggers.add_handler('box center moved',
                self._box_moved_cb)
            self._update_box()
            self.update_graphics()
        elif not flag and self._live:
            self.parent.triggers.remove_handler(self._box_changed_cb)
        self._live = flag

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius):
        self._radius = radius
        if self.visible:
            self.update_graphics()

    @property
    def spotlight_mode(self):
        return self._spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, flag):
        if flag != self._spotlight_mode:
            self._spotlight_mode = flag
            if not flag:
                self.unhide_all_atoms()

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, flag):
        if not flag:
            self.unhide_all_atoms()
        super().set_display(flag)
        if flag:
            if self.spotlight_mode:
                self._center = self.session.view.center_of_rotation
                self._update_box()
            self.update_graphics()

    @property
    def _level_of_detail(self):
        from chimerax.core.atomic.structure import structure_graphics_updater
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
        if self._include_identity:
            first_symop = 0
        else:
            first_symop = 1
        tfs = self._current_tfs = symops.all_matrices_orth(self.cell, format='3x4')
        atoms = self.structure.atoms
        if self._spotlight_mode:
            atoms.hides |=HIDE_ISOLDE
        self._current_atoms, self._current_atom_coords, self._current_atom_syms, \
            self._current_bonds, self._current_bond_tfs, self._current_bond_syms = \
                whole_residue_sym_sphere(self.structure.residues,tfs[first_symop:], center, radius)
                #sym_transforms_in_sphere(atoms, tfs[first_symop:], center, radius)
        self._current_atom_syms += first_symop
        ca = self._current_atoms
        csym = self._current_atom_syms
        if len(ca):
            crs = self._current_ribbon_syms = numpy.unique(csym)
            if self._spotlight_mode:
                ca[csym==0].hides &= ~HIDE_ISOLDE
                self._current_ribbon_syms = crs[crs!=0]
        else:
            if self._spotlight_mode:
                # Don't want the identity operator in this case
                tfs = self._current_tfs = tfs[1:]
            self._current_ribbon_syms = sym_ribbons_in_sphere(self._ribbon_drawing._tether_coords, tfs, center, radius)

    def _model_changed_cb(self, trigger_name, changes):
        if not self.visible:
            return
        changes = changes[1]
        update_needed = False
        if len(changes.created_atoms()):
            update_needed = True
        if changes.num_deleted_atoms() > 0:
            update_needed = True
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons:
            update_needed = True
        if  'hide changed' in reasons:
            hides = self.structure.atoms.hides
            # Prevent repeated callback with every display update in spotlight mode
            current_hides = hides&~HIDE_ISOLDE
            if numpy.any(current_hides != self._last_hides):
                update_needed = True
            self._last_hides = current_hides
        if 'color changed' in reasons:
            update_needed = True
        if 'ribbon_display changed' in changes.residue_reasons():
            self._ribbon_drawing.delayed_rebuild(self.session)
        if (update_needed):
            self._update_box()
            self.update_graphics()

    def update_graphics(self):
        lod = self._level_of_detail
        self._update_atom_graphics(lod)
        self._update_bond_graphics(lod)
        self._update_ribbon_graphics()

    def _update_atom_graphics(self, lod):
        ad = self._atoms_drawing
        lod.set_atom_sphere_geometry(ad)
        ca = self._current_atoms
        ad.visible_atoms = ca
        syms = self._current_atom_syms
        na = len(ca)
        if na > 0:
            ad.display = True
            xyzr = numpy.empty((na, 4), numpy.float32)
            xyzr[:,:3] = self._current_atom_coords
            xyzr[:,3] = self.structure._atom_display_radii(ca)
            from chimerax.core.geometry import Places
            ad.positions = Places(shift_and_scale = xyzr)
            colors = ca.colors.astype(numpy.float32)
            if self._include_identity:
                colors[syms!=0,:3] *= self._color_dim_factor
            else:
                colors[:,:3] *= self._color_dim_factor
            ad.colors = colors.astype(numpy.uint8)
        else:
            ad.display = False


    def _update_bond_graphics(self, lod):
            bd = self._bonds_drawing
            lod.set_bond_cylinder_geometry(bd)
            bonds = bd.visible_bonds = self._current_bonds
            nb = len(bonds)
            if nb > 0:
                bd.display = True
                bd.positions = self._current_bond_tfs
                colors = bonds.half_colors.astype(numpy.float32)
                if self._include_identity:
                    mask = numpy.concatenate([self._current_bond_syms!=0]*2)
                    colors[mask, :3] *= self._color_dim_factor
                else:
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

from chimerax.core.atomic.structure import AtomsDrawing
class SymAtomsDrawing(AtomsDrawing):
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

from chimerax.core.atomic.structure import BondsDrawing
class SymBondsDrawing(BondsDrawing):
    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        return None #too-hard basket for now.

        if not self.display or (exclude and exclude(self)):
            return None
        from chimerax.core.atomic.structure import _bond_intercept
        b, f = _bond_intercept(bonds, mxyz1, mxyz2)

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
            td.vertices, td.normals, td.triangles, td.display = \
                md.vertices, md.normals, md.triangles, md.display
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
            nc.vertices, nc.normals, nc.triangles, nc.positions, nc.display = v, n, t, p, d
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
