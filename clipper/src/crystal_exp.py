# @Author: Tristan Croll
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
import copy
from collections import defaultdict

from chimerax.core.triggerset import TriggerSet
from chimerax.atomic import AtomicStructure, concatenate
from chimerax.core.geometry import Place, Places
from chimerax.core.geometry import find_close_points, find_close_points_sets
from chimerax.surface import zone
from chimerax.surface.shapes import sphere_geometry

from chimerax.core.models import Model, Drawing
from chimerax.std_commands import camera, cofr, cartoon
from chimerax.core.commands import atomspec
from chimerax.map.data import ArrayGridData
from chimerax.map import Volume, volumecommand
from chimerax.map.volume import VolumeSurface

from .mousemodes import initialize_map_contour_mouse_modes
from .main import atom_list_from_sel
# from . import clipper
from .clipper_mtz import ReflectionDataContainer

DEFAULT_BOND_RADIUS = 0.2

def _available_cores():
    import os
    return max(os.cpu_count()-2, 1)

def set_to_default_cartoon(session, model = None):
    '''
    Adjust the ribbon representation to provide information without
    getting in the way.
    '''
    try:
        if model is None:
            atoms = None
        else:
            arg = atomspec.AtomSpecArg('thearg')
            atoms = arg.parse('#' + model.id_string, session)[0]
        cartoon.cartoon(session, atoms = atoms, suppress_backbone_display=False)
        cartoon.cartoon_style(session, atoms = atoms, width=0.4, thickness=0.1, arrows_helix=True, arrow_scale = 2)
        cartoon.cartoon_tether(session, structures=atoms, opacity=0)
    except:
        return

def _model_volume(model, exclude_hydrogens=True, radius_scale=1):
    from math import pi
    atoms = model.atoms
    if exclude_hydrogens:
        atoms = atoms[atoms.element_names != 'H']
    return sum( 4/3 * pi * (atoms.radii*radius_scale)**3 )


def guess_suitable_contour(volume, model, mask_radius=3, atom_radius_scale = 0.5):
    '''
    Find the contour level that would make the volume inside the contour for the
    region immediately surrounding the model approximately equal to the volume
    of the model's atoms scaled by atom_radius_scale.
    '''
    import numpy
    session = model.session

    from .symmetry import is_crystal_map
    is_xmap = is_crystal_map(volume)
    if is_xmap:
        sh = volume.manager.crystal

        spotlight_mode = sh.spotlight_mode
        # Expand the map to cover the whole model
        sh.isolate_and_cover_selection(model.atoms, focus=False)

    from .util import voxel_volume
    vv = voxel_volume(volume)
    mv = _model_volume(model, radius_scale=atom_radius_scale)

    from chimerax.core.geometry import find_close_points
    grid_coords = volume.grid_points(volume.scene_position)

    data = numpy.array(volume.data.matrix(), order='C')
    close_i = find_close_points(model.atoms.scene_coords, grid_coords, mask_radius)[1]
    close_vals = data.ravel()[close_i]


    target_percentile = (1-mv/(vv*len(close_i)))*100
    level = numpy.percentile(close_vals, target_percentile)

    if is_xmap:
        sh.spotlight_mode = spotlight_mode

    return level



#TODO: update Surface_Zone class to handle symmetry atoms
class Surface_Zone:
    '''
    Add this as a property to a Volume object to provide it with the
    necessary information to update its triangle mask after re-contouring.
    '''
    def __init__(self, distance, atoms = None, coords = None):
        '''
        Args:
          distance (float in Angstroms):
            distance from points to which the map will be masked
          atoms:
            Atoms to mask to (coordinates will be updated upon re-masking)
          coords:
            (x,y,z) coordinates to mask to (will not be updated upon
            re-masking).

          Set both atoms and coords to None to disable automatic re-masking.
        '''
        self.update(distance, atoms, coords)

    def update(self, distance, atoms = None, coords = None):
        self.distance = distance
        self.atoms = atoms
        self.coords = coords

    @property
    def all_coords(self):
        if self.atoms is not None:
            if self.coords is not None:
                return numpy.concatenate(self.atoms.coords, self.coords)
            return self.atoms.coords
        return self.coords

def surface_zones(models, points, distance):
    from chimerax.surface import zone
    for m in models:
        for s in m.surfaces:
            spoints = s.position.inverse(is_orthonormal=True) * points
            zone.surface_zone(s, spoints, distance, auto_update=True)


class XmapSet_Live(Model):
    '''
    Handles creation, deletion, recalculation and visualisation of
    crystallographic maps based on the current model and a set of observed
    reflection data (FOBS, SIGFOBS).
    '''

    STANDARD_LOW_CONTOUR = numpy.array([1.5])
    STANDARD_HIGH_CONTOUR = numpy.array([2.5])
    STANDARD_DIFFERENCE_MAP_CONTOURS = numpy.array([-3.0, 3.0])

    DEFAULT_MESH_MAP_COLOR = [0,1.0,1.0,1.0] # Solid cyan
    DEFAULT_SOLID_MAP_COLOR = [0,1.0,1.0,0.4] # Transparent cyan
    DEFAULT_DIFF_MAP_COLORS = [[1.0,0,0,1.0],[0,1.0,0,1.0]] #Solid red and green


    def __init__(self, session, crystal, model, fsigf_name = 'FOBS, SIGFOBS',
                bsharp_vals=[], exclude_free_reflections=False,
                 fill_with_fcalc=False, display_radius = 12,
                 live_update=True, show_r_factors=True):
        '''
        Prepare the C++ Xtal_mgr object and create a set of crystallographic
        maps. The standard 2mFo-DFc and mFo-DFc maps will always be created,
        while any number of maps with different sharpening/smoothing can be
        created by providing a list of b-factors in bsharp_vals.
        Args:
            session:
                The ChimeraX session
            crystal:
                The parent XtalSymmetryHandler object
            fsigf_name:
                The label of the :class:`F_sigF_float` object containing the
                observed amplitudes
            model:
                A ChimeraX `AtomicStructure` object to be used to
                calculate the maps.
            bsharp_vals:
                For each value in this list, a 2mFo-DFc map will be generated
                with the given B_sharp value. A negative B_sharp yields a
                sharpened map, while a positive B_sharp gives smoothing. As a
                rough rule of thumb, a value of B_sharp=-100 tends to be
                reasonable for a 3.8 Angstrom map, and -50 for 3 Angstroms.
                For maps with resolutions better than ~2 Angstroms, it may be
                useful in some circumstances to apply smoothing. This can
                sometimes bring out weaker, low-resolution features hidden by
                high-resolution noise in standard maps.
            exclude_free_reflections:
                If True, observed amplitudes corresponding to the free set will
                not be used in generating the maps. The values used in their
                place will depend on the value of `fill_with_fcalc`. Note that
                this only affects maps generated for viewing - the MDFF potential
                map is always generated with free reflections excluded.
            fill_with_fcalc:
                If `exclude_free_reflections` is False this argument will be
                ignored. Otherwise, if `fill_with_fcalc` is True then the
                excluded amplitudes will be replaced by sigmaa-weighted F(calc).
                If False, the excluded amplitudes will be set to zero.
            display_radius:
                The radius (in Angstroms) of the display sphere used in
                live scrolling mode.
            live_update:
                If True, maps will be automatically recalculated whenever
                coordinates change
            show_r_factors:
                If True, print new R factors to the status bar on each map
                recalculation
        '''
        Model.__init__(self, 'Live real-space maps', session)
        self.crystal = crystal
        from chimerax.core.triggerset import TriggerSet
        trig = self.triggers = TriggerSet()

        trigger_names = (
            'map box changed',  # Changed shape of box for map viewing
            'map box moved',    # Just changed the centre of the box
            'maps recalculated'
        )
        for t in trigger_names:
            trig.add_trigger(t)


        self._live_update = False
        self._recalc_needed = False
        self._model_changes_handler = None
        self._delayed_recalc_handler = None
        self._show_r_factors = show_r_factors

        self.model = model

        # Since the calculation of crystallographic maps is directly reliant on
        # the atomic model, we need to ensure that deletion of the model deletes
        # this class

        from chimerax.core.models import REMOVE_MODELS
        self._model_removed_handler = session.triggers.add_handler(
            REMOVE_MODELS, self._model_removed_cb)

        atoms = model.atoms
        # Handler for live box update
        self._box_update_handler = None
        # Is the box already initialised?
        self._box_initialized = False
        # Object storing the parameters required for masking (used after
        # adjusting contours)
        self._surface_zone = Surface_Zone(display_radius, None, None)
        # Is the map box moving with the centre of rotation?
        self._live_scrolling = False
        # Radius of the sphere in which the map will be displayed when
        # in live-scrolling mode
        self.display_radius = display_radius
        # Actual box dimensions in (u,v,w) grid coordinates
        self._box_dimensions = None
        # Centre of the box (used when tracking the centre of rotation)
        self._box_center = None
        # Last grid coordinate of the box centre. We only need to update
        # the map if the new centre maps to a different integer grid point
        self._box_center_grid = None
        # Minimum corner of the box in (x,y,z) coordinates. This must
        # correspond to the grid coordinate in self._box_corner_grid
        self._box_corner_xyz = None
        # Minimum corner of the box in grid coordinates
        self._box_corner_grid = None


        self.live_scrolling = True

        self._box_initialized = True

        # The master C++ manager for handling all map operations
        from .clipper_python.ext import Xtal_thread_mgr
        xm = self._xtal_mgr = Xtal_thread_mgr(crystal.hklinfo,
            crystal.mtzdata.free_flags.data, crystal.grid,
            crystal.mtzdata.experimental_data.datasets[fsigf_name].data,
            num_threads=_available_cores())

        from . import atom_list_from_sel
        xm.init(atom_list_from_sel(atoms))
        # Only this map will actually be used as the MDFF potential
        self.add_live_xmap('MDFF potential', 0, is_difference_map=False,
            exclude_missing_reflections=True,
            exclude_free_reflections=True, fill_with_fcalc=True,
            display=False)

        self.add_live_xmap('2mFo-DFc', 0, is_difference_map=False,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc=fill_with_fcalc,
            display=True)

        self.add_live_xmap('mFo-DFc', 0, is_difference_map=True,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc=fill_with_fcalc,
            display=True)

        for b in bsharp_vals:
            if abs(b) < 5:
                continue
            elif b < 0:
                name_str = "2mFo-DFc_smooth_{:.0f}".format(-b)
            else:
                name_str = "2mFo-DFc_sharp_{:.0f}".format(b)
            self.add_live_xmap(name_str, b, is_difference_map=False,
                exclude_free_reflections=exclude_free_reflections,
                fill_with_fcalc=fill_with_fcalc
            )

        self.display=False
        # Apply the surface mask
        self._rezone_pending=False
        self.session.triggers.add_handler('frame drawn', self._first_init_cb)
        self.live_update = live_update

    def _first_init_cb(self, *_):
        self.display = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _model_removed_cb(self, trigger_name, removed_models):
        if self.model in removed_models:
            if not self.deleted:
                if self in self.session.models.list():
                    self.session.models.remove([self])
            self.delete()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def rezone_needed(self):
        if not self._rezone_pending:
            self.session.triggers.add_handler('frame drawn', self._rezone_once_cb)

    def _rezone_once_cb(self, *_):
        self._reapply_zone()
        self._rezone_pending=False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    @property
    def show_r_factors(self):
        return self._show_r_factors

    @show_r_factors.setter
    def show_r_factors(self, flag):
        self._show_r_factors = flag

    @property
    def live_update(self):
        return self._live_update

    @live_update.setter
    def live_update(self, flag):
        if flag == self._live_update:
            return
        if flag:
            if self._model_changes_handler is None:
                self._model_changes_handler = self.model.triggers.add_handler(
                    'changes', self._model_changed_cb
                )
        else:
            if self._model_changes_handler is not None:
                self.model.triggers.remove_handler(
                    self._model_changes_handler
                )
                self._model_changes_handler = None
        self._live_update = flag

    @property
    def rfree(self):
        return self._xtal_mgr.rfree

    @property
    def rwork(self):
        return self._xtal_mgr.rwork

    @property
    def hklinfo(self):
        return self.crystal.hklinfo

    @property
    def spacegroup(self):
        return self.crystal.spacegroup

    @property
    def cell(self):
        return self.crystal.cell

    @property
    def res(self):
        return self.hklinfo.resolution.limit

    @property
    def resolution(self):
        return self.hklinfo.resolution.limit

    @property
    def grid(self):
        return self.crystal.grid

    @property
    def voxel_size(self):
        return self.cell.dim / self.grid.dim

    @property
    def voxel_size_frac(self):
        return 1/ self.grid.dim

    @property
    def unit_cell(self):
        return self.crystal.unit_cell

    @property
    def display_radius(self):
        return self._display_radius

    @display_radius.setter
    def display_radius(self, radius):
        '''Set the radius (in Angstroms) of the live map display sphere.'''
        self._display_radius = radius
        v = self.session.view
        cofr = self._box_center = v.center_of_rotation
        from .clipper_python import Coord_orth
        self._box_center_grid = Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
        dim = self._box_dimensions = \
            2 * calculate_grid_padding(radius, self.grid, self.cell)
        self._box_corner_grid, self._box_corner_xyz = _find_box_corner(
            cofr, self.cell, self.grid, radius)
        self.triggers.activate_trigger('map box changed',
            ((self._box_corner_xyz, self._box_corner_grid, dim), False))
        self._surface_zone.update(radius, coords = numpy.array([cofr]))
        self._reapply_zone()

    @property
    def live_scrolling(self):
        '''Turn live map scrolling on and off.'''
        return self._live_scrolling

    @live_scrolling.setter
    def live_scrolling(self, switch):
        if switch:
            self.position = Place()
            if not self._live_scrolling:
                '''
                Set the box dimensions to match the stored radius.
                '''
                self.display_radius = self._display_radius
            self._start_live_scrolling()
        else:
            self._stop_live_scrolling()

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, switch):
        if switch:
            if self.live_scrolling:
                self._start_live_scrolling()
        Model.display.fset(self, switch)

    def _start_live_scrolling(self):
        if self._box_update_handler is None:
            self._box_update_handler = self.crystal.triggers.add_handler(
                'box center moved', self.update_box)
        if self._box_center is not None:
            self.update_box(None, (self._box_center, self._box_center_grid), force=True)
        self.positions = Places()
        self._live_scrolling = True

    def _stop_live_scrolling(self):
        if self._box_update_handler is not None:
            self.crystal.triggers.remove_handler(self._box_update_handler)
            self._box_update_handler = None
        self._live_scrolling = False

    def items(self):
        return ((m.name, m) for m in self.child_models() if isinstance(m, XmapHandler_Live))

    def __getitem__(self, name_or_index):
        '''Get one of the child maps by name or index.'''
        if type(name_or_index) == str:
            for m in self.child_models():
                if m.name == name_or_index:
                    return m
            raise KeyError('No map with that name!')
        else:
            return self.child_models()[name_or_index]


    def set_box_limits(self, minmax, force_fill = False):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum grid coordinates. Automatically turns off live scrolling.
        '''
        self.live_scrolling = False
        from .clipper_python import Coord_grid
        cmin = Coord_grid(minmax[0])
        cmin_xyz = cmin.coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = (minmax[1]-minmax[0])
        self.triggers.activate_trigger('map box changed',
            ((cmin_xyz, cmin, dim), force_fill))

    def cover_unit_cells(self, nuvw = [1,1,1], offset = [0,0,0]):
        '''
        Expand the map(s) to cover multiple unit cells. In order to
        maintain reasonable performance, this method cheats a little by
        filling just one unit cell and then tiling it using the graphics
        engine. This leaves some minor artefacts at the cell edges, but
        is a worthwhile tradeoff.
        Automatically turns off live scrolling.
        Args:
            nuvw (array of 3 positive integers):
                Number of unit cells to show in each direction.
            offset (array of 3 integers):
                Shifts the starting corner of the displayed volume by
                this number of unit cells in each direction.
        '''
        self.live_scrolling = False
        uc = self.unit_cell
        box_min_grid = uc.min.uvw
        # Add a little padding to the max to create a slight overlap between copies
        from .clipper_python import Coord_grid
        box_max_grid = (uc.max+Coord_grid([2,2,2])).uvw
        minmax = [box_min_grid, box_max_grid]
        self.set_box_limits(minmax)
        self._surface_zone.update(None, None)
        # Tile by the desired number of cells
        places = []
        grid_dim = self.grid.dim
        nu, nv, nw = nuvw
        ou, ov, ow = offset
        for i in range(ou, nu+ou):
            for j in range(ov, nv+ov):
                for k in range(ow, nw+ow):
                    thisgrid = Coord_grid(numpy.array([i,j,k])*grid_dim)
                    thisorigin = thisgrid.coord_frac(self.grid).coord_orth(self.cell).xyz
                    places.append(Place(origin = thisorigin))
        self.positions = Places(places)

    _map_impacting_changes = set ((
        "aniso_u changed",
        "bfactor changed",
        "coord changed",
        "coordset changed",
        "occupancy changed",
    ))
    def _model_changed_cb(self, trigger_name, changes):
        if changes is not None:
            changes = set(changes[1].atom_reasons()).intersection(self._map_impacting_changes)
            if changes:
                self._recalc_needed = True
                if self._delayed_recalc_handler is None:
                    self._delayed_recalc_handler = self.session.triggers.add_handler(
                        'new frame', self._recalculate_maps_if_needed
                    )

    def _recalculate_maps_if_needed(self, *_):
        xm = self._xtal_mgr
        if self._recalc_needed and not xm.thread_running:
            self.recalculate_all_maps(self.model.atoms)
            if self._delayed_recalc_handler is not None:
                self.session.triggers.remove_handler(self._delayed_recalc_handler)
                self._delayed_recalc_handler = None

    def recalculate_all_maps(self, atoms):
        from . import atom_list_from_sel
        from .delayed_reaction import delayed_reaction
        xm = self._xtal_mgr
        delayed_reaction(self.session.triggers, 'new frame',
            xm.recalculate_all_maps, [atom_list_from_sel(atoms)],
            xm.ready,
            self._apply_new_maps, []
            )
        self._recalc_needed = False

    def _apply_new_maps(self):
        self._xtal_mgr.apply_new_maps()
        if self.show_r_factors:
            self.session.logger.status(
                'R-work: {:0.4f}  Rfree: {:0.4f}'.format(
                    self._xtal_mgr.rwork, self._xtal_mgr.rfree
                )
            )
        self.triggers.activate_trigger('maps recalculated', None)

    def add_nxmap(self, volume, is_difference_map=False, color=None,
        style=None, contour=None, delete_original=True, display=True):
        '''
        Add a new real-space map from an existing Volume object. Makes a copy of
        the original Volume's data, and optionally closes the original.
        '''
        from .real_space_map import NXmapHandler
        new_handler = NXmapHandler(self.session, self, volume,
            is_difference_map = is_difference_map)

        if style is None:
            style = 'mesh'
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]

        if contour is None:
            if is_difference_map:
                contour = self.STANDARD_DIFFERENCE_MAP_CONTOURS
            else:
                contour = [guess_suitable_contour(new_handler, self.model)]
        self.add([new_handler])
        new_handler.set_representation(style)
        new_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False})
        if delete_original:
            self.session.models.remove([volume])
            volume.delete()
        self._reapply_zone()

        return new_handler

    def add_live_xmap(self, name, b_sharp,
        is_difference_map=False,
        exclude_missing_reflections=False,
        exclude_free_reflections=True,
        fill_with_fcalc=True,
        color=None, style=None, contour=None, display=True):
        xm = self._xtal_mgr
        xm.add_xmap(name, b_sharp, is_difference_map=is_difference_map,
            exclude_missing_reflections=exclude_missing_reflections,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc = fill_with_fcalc)
        #xmap = xm.get_xmap_ref(name)
        new_handler = self.add_xmap_handler(name, is_difference_map=is_difference_map,
            color=color, style=style, contour=contour)
        if display:
            new_handler.show()
            self._reapply_zone()
        else:
            new_handler.display = False

    def add_xmap_handler(self, name, is_difference_map = None,
                color = None, style = None, contour = None):
        '''
        Add a new XmapHandler_Live based on the given reflections and phases.
        Args:
            name:
                A unique string describing this map
            xmap:
                a Clipper :class:`Xmap_float` object.
            is_difference_map:
                Decides whether this map is to be treated as a difference
                map (with positive and negative contours) or a standard
                map with positive contours only. Leave as None to allow
                this to be determined automatically from the dataset,
                or set it to True or False to force an explicit choice.
            color:
                an array of [r, g, b, alpha] as integers in the range 0-255
                for a positive density map, or
                [[r,g,b,alpha],[r,g,b,alpha]] representing the colours of
                the positive and negative contours respectively if the
                map is a difference map
            style:
                one of 'mesh' or 'surface'
            contour:
                The value(s) (in sigma units) at which to contour the map
                on initial loading. For a standard map this should be a
                single value; for a difference map it should be
                [negative contour, positive contour]
        '''
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        new_handler = XmapHandler_Live(self.session, self, name,
            self._box_corner_xyz, self._box_corner_grid, self._box_dimensions,
            is_difference_map = is_difference_map)
        if style is None:
            style = 'mesh'
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]
        if contour is None:
            if is_difference_map:
                contour = self.STANDARD_DIFFERENCE_MAP_CONTOURS
            else:
                contour = self.STANDARD_LOW_CONTOUR
        elif not hasattr(contour, '__len__'):
                contour = numpy.array([contour])
        else:
            contour = numpy.array(contour)
        contour = contour * new_handler.sigma
        self.add([new_handler])
        new_handler.set_representation(style)
        new_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False})
        # new_handler.update_surface()
        return new_handler

    def update_box(self, trigger_name, new_center, force=True):
        '''Update the map box to surround the current centre of rotation.'''
        center, center_grid = new_center
        self._box_center = center
        self._box_center_grid = center_grid
        if not self.visible and not force:
            # Just store the box parameters for when we're re-displayed
            return
        if self.live_scrolling:
            box_corner_grid, box_corner_xyz = _find_box_corner(center, self.cell, self.grid, self.display_radius)
            self.triggers.activate_trigger('map box moved', (box_corner_xyz, box_corner_grid, self._box_dimensions))
            self._surface_zone.update(self.display_radius, coords = numpy.array([center]))
            self._reapply_zone()

    def _reapply_zone(self):
        '''
        Reapply any surface zone applied to the volume after changing box
        position.
        '''
        coords = self._surface_zone.all_coords
        radius = self._surface_zone.distance
        if coords is not None:
            surface_zones(self, coords, radius)

    def removed_from_session(self, session):
        if self._model_removed_handler is not None:
            session.triggers.remove_handler(self._model_removed_handler)
        super().removed_from_session(session)
        self.crystal.xmapset = None

    def delete(self):
        self.live_scrolling = False
        self.live_update = False
        super().delete()

def map_potential_bsharp(resolution):
    '''
    Return a recommended sharpening/smoothing B-factor for a given map to
    optimise its use as an MDFF potential. For now this is a simple linear
    function of resolution, passing through zero at 2.5 Angstroms (smoothing
    below, sharpening above). Improved methods will be developed over time.
    '''
    # smooth by 30 A**2 at 1.5A res; sharpen by 30A**2 at 3.5A res; 0 at 2.5A res
    bsharp_base = 30
    return bsharp_base*resolution-2.5*bsharp_base

def viewing_bsharp(resolution):
    '''
    For viewing purposes it is also often useful to have a smoothed or
    sharpened visualisation of your map, but the optimal degree of sharpening
    for visualisation is not necessarily the same as that for MDFF.
    '''
    # smooth by 50 A**2 at 1.5A res; sharpen by 50 A**2 at 3.5A res, 0 at 2.5A res
    bsharp_base = 50
    return bsharp_base*resolution-2.5*bsharp_base

def calculate_grid_padding(radius, grid, cell):
    '''
    Calculate the number of grid steps needed on each crystallographic axis
    in order to capture at least radius angstroms in x, y and z.
    '''
    corner_mask = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,0,0],[1,1,1]])
    corners = (corner_mask * radius).astype(float)
    grid_upper = numpy.zeros([8,3], numpy.int)
    grid_lower = numpy.zeros([8,3], numpy.int)
    from .clipper_python import Coord_orth
    for i, c in enumerate(corners):
        co = Coord_orth(c)
        cm = co.coord_frac(cell).coord_map(grid).uvw
        grid_upper[i,:] = numpy.ceil(cm).astype(int)
        grid_lower[i,:] = numpy.floor(cm).astype(int)
    return grid_upper.max(axis=0) - grid_lower.min(axis=0)

def find_box_corner(center, cell, grid, radius=20):
    return _find_box_corner(center, cell, grid, radius)

def _find_box_corner(center, cell, grid, radius = 20):
    '''
    Find the bottom corner (i.e. the origin) of a rhombohedral box
    big enough to hold a sphere of the desired radius.
    '''
    from .clipper_python import Coord_frac, Coord_orth, Coord_grid
    radii_frac = Coord_frac(radius/cell.dim)
    center_frac = Coord_orth(center).coord_frac(cell)
    bottom_corner_grid = center_frac.coord_grid(grid) \
                - Coord_grid(calculate_grid_padding(radius, grid, cell))
    bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
    return bottom_corner_grid, bottom_corner_orth.xyz

def _get_bounding_box(coords, padding, grid, cell):
    '''
    Find the minimum and maximum grid coordinates of a box which will
    encompass the given (x,y,z) coordinates plus padding (in Angstroms).
    '''
    from .clipper_python import Util, Coord_grid
    grid_pad = calculate_grid_padding(padding, grid, cell)
    box_bounds_grid = Util.get_minmax_grid(coords, cell, grid)\
                        + numpy.array((-grid_pad, grid_pad))
    box_origin_grid = box_bounds_grid[0]
    box_origin_xyz = Coord_grid(box_origin_grid).coord_frac(grid).coord_orth(cell)
    dim = box_bounds_grid[1] - box_bounds_grid[0]
    return [box_origin_grid, box_origin_xyz, dim]

from .clipper_python import Xmap_float
class Xmap(Xmap_float):
    def __init__(self, spacegroup, cell, grid_sampling,
                 name = None, hkldata = None, is_difference_map = None):
        super().__init__(spacegroup, cell, grid_sampling)
        self.name = name
        self.is_difference_map = is_difference_map
        if hkldata is not None:
            self.fft_from(hkldata)
        from .clipper_python import Map_stats
        self._stats = Map_stats(self)

    @property
    def stats(self):
        if self._stats is None:
            from .clipper_python import Map_stats
            self._stats = Map_stats(self)
        return self._stats

    @property
    def mean(self):
        return self.stats.mean

    @property
    def std_dev(self):
        return self.stats.std_dev

    @property
    def sigma(self):
        return self.stats.std_dev

    @property
    def min(self):
        return self.stats.min

    @property
    def max(self):
        return self.stats.max

    @property
    def range(self):
        return self.stats.range




class XmapHandler_Live(Volume):
    '''
    An XmapHandler_Live is in effect a resizable window into a periodic
    crystallographic map. The actual map data (a clipper Xmap object) is
    held within, and filled into the XmapHandler_Live.data array as needed.
    Methods are included for both live updating (e.g. tracking and filling
    a box centred on the centre of rotation) and static display of a
    given region.
    '''
    def __init__(self, session, manager, name, origin, grid_origin, dim,
        is_difference_map = False):
        '''
        Args:
            sesssion:
                The ChimeraX session
            manager:
                The XmapSet_Live object this belongs to
            name:
                A descriptive name for this map
            origin:
                The (x,y,z) coordinates of the bottom left corner of the
                volume.
            grid_origin:
                The (u,v,w) integer grid coordinates corresponding to
                origin.
            dim:
                The shape of the box in (u,v,w) grid coordinates.
            is_difference_map:
                Is this a difference map?
        '''
        self.box_params = (origin, grid_origin, dim)
        self.manager = manager
        darray = self._generate_data_array(origin, grid_origin, dim)
        Volume.__init__(self, darray, session)
        self.name = name
        self._fill_volume_data(self._data_fill_target, grid_origin)
        self.is_difference_map = is_difference_map
        self.initialize_thresholds()

        # If the box shape changes while the volume is hidden, the change
        # will not be applied until it's shown again.
        self._needs_update = True
        self.show()
        self._box_shape_changed_cb_handler = self.manager.triggers.add_handler(
            'map box changed', self._box_changed_cb)
        self._box_moved_cb_handler = self.manager.triggers.add_handler(
            'map box moved', self._box_moved_cb)
        self._map_recalc_cb_handler = self.manager.triggers.add_handler(
            'maps recalculated', self._map_recalc_cb
        )

    @property
    def xmap(self):
        return self.manager._xtal_mgr.get_xmap_ref(self.name)

    @property
    def stats(self):
        return self.manager._xtal_mgr.get_map_stats(self.name)

    def show(self, *args, **kwargs):
        if self._needs_update:
            self._swap_volume_data(self.box_params, force_update = True)
            self._needs_update = False
        else:
            # Just set the origin and fill the box with the data for
            # the current location
            origin, grid_origin, ignore = self.box_params
            self._fill_volume_data(self._data_fill_target, grid_origin)
        super().show(*args, **kwargs)

    @property
    def hklinfo(self):
        return self.manager.hklinfo

    @property
    def spacegroup(self):
        return self.manager.spacegroup

    @property
    def cell(self):
        return self.manager.cell

    @property
    def res(self):
        return self.hklinfo.resolution

    @property
    def grid(self):
        return self.manager.grid

    @property
    def voxel_size(self):
        return self.cell.dim / self.grid.dim

    @property
    def voxel_size_frac(self):
        return 1/ self.grid.dim

    @property
    def unit_cell(self):
        return self.manager.unit_cell

    @property
    def _surface_zone(self):
        return self.manager._surface_zone

    def mean_sd_rms(self):
        '''
        Overrides the standard Volume method to give the overall values
        from the Clipper object.
        '''
        s = self.stats
        return (s.mean, s.std_dev, s.std_dev)

    @property
    def sigma(self):
        return self.stats.std_dev

    def add_surface(self, level, rgba=None):
        ses = self.session
        s = FastVolumeSurface(self, level, rgba)
        self._surfaces.append(s)
        if self.id is None:
            self.add([s])
        else:
            ses.models.add([s], parent=self)
        return s

    def force_refill(self):
        self._box_changed_cb(None, (self.box_params, True))

    def _box_changed_cb(self, name, params):
        box_params, force_fill = params
        self.box_params = box_params
        self._needs_update = True
        if not self.display and not force_fill:
            # No sense in wasting cycles on this if the volume is hidden.
            # We'll just store the params and apply them when we show the
            # volume.
            # NOTE: this means we need to over-ride show() to ensure
            # it's updated before re-displaying.
            return
        self._swap_volume_data(box_params)
        self._use_thread = True
        self.data.values_changed()
        # self.call_change_callbacks('data values changed')

    def _box_moved_cb(self, name, params):
        self.box_params = params
        if not self.display:
            return
        self.data.set_origin(params[0])
        self._fill_volume_data(self._data_fill_target, params[1])
        for s in self.surfaces:
            s._use_thread=True
        self.data.values_changed()

    def _map_recalc_cb(self, name, *_):
        for s in self.surfaces:
            s._use_thread=True
        self._fill_volume_data(self._data_fill_target, self.box_params[1])
        self.data.values_changed()


    def delete(self):
        bh = self._box_shape_changed_cb_handler
        if bh is not None:
            self.manager.triggers.remove_handler(bh)
            self._box_shape_changed_cb_handler = None
        bm = self._box_moved_cb_handler
        if bm is not None:
            self.manager.triggers.remove_handler(bm)
            self._box_moved_cb_handler = None
        super().delete()


    def _swap_volume_data(self, params, force_update = False):
        '''
        Replace this Volume's data array with one of a new shape/size
        Args:
            params:
                A tuple of (new_origin, new_grid_origin, new_dim)
        '''
        if not self._needs_update and not force_update:
            # Just store the parameters
            self.box_params = params
            return
        new_origin, new_grid_origin, new_dim = params
        darray = self._generate_and_fill_data_array(new_origin, new_grid_origin, new_dim)
        self._box_dimensions = new_dim
        self.replace_data(darray)
        self.new_region(ijk_min=(0,0,0), ijk_max=darray.size, ijk_step=(1,1,1), adjust_step=False)

    def _generate_and_fill_data_array(self, origin, grid_origin, dim):
        darray = self._generate_data_array(origin, grid_origin, dim)
        self._fill_volume_data(self._data_fill_target, grid_origin)
        return darray

    def _generate_data_array(self, origin, grid_origin, dim):
        data = self._data_fill_target = numpy.empty(dim, numpy.float32)
        order = numpy.array([2,1,0], int)
        darray = ArrayGridData(data.transpose(), origin = origin,
            step = self.voxel_size, cell_angles = self.cell.angles_deg)
        return darray


    def _fill_volume_data(self, target, start_grid_coor):
        #shape = (numpy.array(target.shape)[::-1] - 1)
        #end_grid_coor = start_grid_coor + clipper.Coord_grid(shape)
        #self.data.set_origin(origin_xyz)
        from .clipper_python import Coord_grid
        xmap = self.xmap
        xmap.export_section_numpy(Coord_grid(start_grid_coor), target)

class FastVolumeSurface(VolumeSurface):
    def __init__(self, volume, level, rgba=(1.0, 1.0, 1.0, 1.0)):
        super().__init__(volume, level, rgba)
        self._update_needed = False

    def _postprocess(self, varray, narray, tarray, rendering_options, level):
        ro = rendering_options
        if ro.flip_normals and level < 0:
          from chimerax.surface import invert_vertex_normals
          invert_vertex_normals(narray, tarray)

        # reverse_triangle_vertex_order done in thread

        if ro.subdivide_surface:
            from chimerax.surface import subdivide_triangles
            for i in range(ro.subdivision_levels):
                varray, tarray, narray = subdivide_triangles(varray, tarray, narray)

        if ro.square_mesh:
            from numpy import empty, uint8
            hidden_edges = empty((len(tarray),), uint8)
            from chimerax.map import _map
            _map.principle_plane_edges(varray, tarray, hidden_edges)
        else:
            hidden_edges = None

        if ro.surface_smoothing:
          sf, si = ro.smoothing_factor, ro.smoothing_iterations
          from chimerax.surface import smooth_vertex_positions
          smooth_vertex_positions(varray, tarray, sf, si)
          smooth_vertex_positions(narray, tarray, sf, si)

        # Transforms and normalization done in thread
        return varray, narray, tarray, hidden_edges

    def _use_fast_thread_result(self, show_mesh, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None:
            va, ta, na = sct.get_result()
            va, na, ta, hidden_edges = self._postprocess(va, na, ta, self.volume.rendering_options, self.level)
            self._set_surface(va, na, ta, hidden_edges)
            self._set_appearance(show_mesh, rendering_options)
            self._surf_calc_thread = None
            if self._update_needed:
                # surface properties were changed while the thread was working
                self.update_surface(show_mesh, rendering_options)
            self._update_needed = False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def update_surface(self, show_mesh, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None and not sct.ready():
            self._update_needed = True
            return

        v = self.volume
        level = self.level
        vertex_transform = v.matrix_indices_to_xyz_transform()
        normal_transform = vertex_transform.zero_translation()
        #normal_transform = vertex_transform.inverse().transpose().zero_translation()
        det = vertex_transform.determinant()

        from .delayed_reaction import delayed_reaction
        from .contour_thread import Contour_Thread_Mgr
        sct = self._surf_calc_thread = Contour_Thread_Mgr()
        delayed_reaction(self.volume.session.triggers, 'new frame',
            sct.start_compute, (v.matrix(), level, det, vertex_transform.matrix, normal_transform.matrix, False, True),
            sct.ready,
            self._use_fast_thread_result, (show_mesh, rendering_options))
