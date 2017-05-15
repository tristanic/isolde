import numpy
import copy
from collections import defaultdict

from chimerax.core.triggerset import TriggerSet
from chimerax.core.atomic import AtomicStructure, concatenate
from chimerax.core.geometry import Place, Places
from chimerax.core.geometry import find_close_points, find_close_points_sets
from chimerax.core.surface import zone
from chimerax.core.surface.shapes import sphere_geometry

from chimerax.core.models import Model, Drawing
from chimerax.core.commands import camera, cofr, cartoon, atomspec
from chimerax.core.map.data import Array_Grid_Data
from chimerax.core.map import Volume, volumecommand

from . import initialize_mouse_modes
from .main import atom_list_from_sel
from . import clipper
from .lib import clipper_python_core as clipper_core
from .data_tree import db_levels, DataTree
from .clipper_mtz_new import ReflectionDataContainer

DEFAULT_BOND_RADIUS = 0.2


def move_model(session, model, new_parent):
    '''
    Temporary method until something similar is added to the ChimeraX
    core. Picks up a model from the ChimeraX model tree and transplants
    it (with all its children intact) as the child of a different model.
    '''

    mlist = model.all_models()
    model_id = model.id
    if new_parent in mlist:
        raise RuntimeError('Target model cannot be one of the models being moved!')
    for m in mlist:
        m.removed_from_session(session)
        mid = m.id
        if mid is not None:
            del session.models._models[mid]
            m.id = None
    session.triggers.activate_trigger('remove models', mlist)
    if len(model_id) == 1:
        parent = session.models.drawing
    else:
        parent = session.models._models[model_id[:-1]]
    parent.remove_drawing(model, delete=False)
    parent._next_unused_id = None
    new_parent.add([model])

def symmetry_from_model_metadata(model):
    '''
    Generate Cell, Spacegroup and a default Grid_Sampling from the PDB
    CRYST1 card (or mmCIF equivalent metadata once it's available in
    ChimeraX).
    '''
    if 'CRYST1' in model.metadata.keys():
        cryst1 = model.metadata['CRYST1'][0].split()
        abc = cryst1[1:4]
        angles = cryst1[4:7]

        remarks = model.metadata['REMARK']
        i = 0

        '''
        Get the resolution. We need this to define a Grid_sampling
        for the unit cell (needed even in the absence of a map since
        atomic symmetry lookups are done with integerised symops for
        performance). We want to be as forgiving as possible at this
        stage - we'll use the official resolution if we find it, and
        set a default resolution if we don't. This will be overridden
        by the value from any mtz file that's loaded later.
        '''
        try:
            while 'REMARK   2' not in remarks[i]:
                i += 1
            # The first 'REMARK   2' line is blank by convention, and
            # resolution is on the second line
            i += 1
            line = remarks[i].split()
            res = line[3]
        except:
            res = 3.0

        '''
        The spacegroup identifier tends to be the most unreliable part
        of the CRYST1 card, so it's considered safer to let Clipper
        infer it from the list of symmetry operators at remark 290. This
        typically looks something like the following:

        REMARK 290      SYMOP   SYMMETRY
        REMARK 290     NNNMMM   OPERATOR
        REMARK 290       1555   X,Y,Z
        REMARK 290       2555   -X,-Y,Z+1/2
        REMARK 290       3555   -Y+1/2,X+1/2,Z+1/4
        REMARK 290       4555   Y+1/2,-X+1/2,Z+3/4
        REMARK 290       5555   -X+1/2,Y+1/2,-Z+1/4
        REMARK 290       6555   X+1/2,-Y+1/2,-Z+3/4
        REMARK 290       7555   Y,X,-Z
        REMARK 290       8555   -Y,-X,-Z+1/2

        Clipper is able to initialise a Spacegroup object from a
        string containing a semicolon-delimited list of the symop
        descriptors in the SYMMETRY OPERATOR column, so we need to
        parse those out.
        '''
        # Find the start of the REMARK 290 section
        while remarks[i][0:10] != 'REMARK 290':
            i += 1
        while 'NNNMMM' not in remarks[i]:
            i += 1
        i+=1
        symstr = ''
        thisline = remarks[i]
        while 'X' in thisline and 'Y' in thisline and 'Z' in thisline:
            if len(symstr):
                symstr += ';'
            splitline = thisline.split()
            symstr += splitline[3]
            i+=1
            thisline = remarks[i]




    # TODO: add equivalent lookup for mmCIF metadata once available

    cell_descr = clipper.Cell_descr(*abc, *angles)
    cell = clipper.Cell(cell_descr)
    spgr_descr = clipper.Spgr_descr(symstr)
    spacegroup = clipper.Spacegroup(spgr_descr)
    resolution = clipper.Resolution(res)
    grid_sampling = clipper.Grid_sampling(spacegroup, cell, resolution)
    return cell, spacegroup, grid_sampling





class CrystalStructure(Model):
    '''
    Master container class for a crystal structure, designed to act as the
    head node of a Model tree with the following general format:

    CrystalStructure
      |
      -- master model (AtomicStructure)
      |
      -- symmetry copies (SymModels)
      |    |
      |    -- AtomicStructure
      |    |
      |    -- AtomicStructure
      |    |
      |    -- ...
      |
      -- reciprocal space data (ReflectionDataContainer)
      |    |
      |    -- Free Flags (ReflectionData_FreeFlags)
      |    |
      |    -- Experimental (ReflectionData_Node)
      |    |    |
      |    |    -- F/SigF (ReflectionData_Exp)
      |    |    |
      |    |    -- ...
      |    |
      |    -- Calculated (ReflectionData_Calc)
      |         |
      |         -- 2mFo-DFc (ReflectionData_Calc)
      |         |
      |         -- ...
      |
      |
      -- real-space maps (XMapSet)
           |
           -- 2mFo-Fc (Volume)
           |
           -- ...
    '''

    def __init__(self, session, model, 
                 mtzfile = None, calculate_maps = True,
                 live_map_scrolling = True, map_scrolling_radius = 12,
                 live_atomic_symmetry = True, atomic_symmetry_radius = 15,
                 show_nonpolar_H = False):
        '''
        Create a new crystal structure object from an atomic model and
        (optionally) a set of reciprocal-space data.
        Args:
            session:
                The ChimeraX session.
            model:
                A loaded AtomicStructure model. NOTE: this will be moved from
                its existing place in the session.models tree to become a
                child of this one.
            mtzfile (string):
                The name of an MTZ file containing experimental and/or 
                calculated amplitudes and phases
            calculate_maps (bool, default True):
                If an MTZ file containing map structure factors is provided,
                generate an XmapSet containing the real-space maps.
            live_map_scrolling (bool, default True):
                If maps are generated, initialise live scrolling of a 
                sphere of density around the centre of rotation.
            map_scrolling_radius (float, default 12):
                The radius (in Angstroms) of the live map sphere.
            live_atomic_symmetry (bool, default True):
                If true, switch to a mode in which atoms (including 
                symmetry atoms) within the given radius are displayed as
                sticks, while the remainder of the structure is shown as
                a minimal cartoon. Cartoons of symmetry-related molecules
                will be automatically shown/hidden as they enter/leave 
                the sphere. 
            atomic_symmetry_radius (float, default 15):
                Radius (in Angstroms) of the live atomic symmetry sphere.
            show_nonpolar_H (bool, default False):
                Do you want non-polar hydrogens to be visible in all default
                visualisations?
        '''
        name = 'Crystal (' + model.name +')'
        Model.__init__(self, name, session)

        initialize_mouse_modes(session)
        volumecommand.volume(session, pickable=False)

        self.session.models.add([self])
        move_model(self.session, model, self)
        self.master_model = model
        self.mtzdata = None
        if mtzfile is not None:
            self.mtzdata = ReflectionDataContainer(self.session, mtzfile)
            self.add([self.mtzdata])
            self.cell = self.mtzdata.cell
            self.spacegroup = self.mtzdata.spacegroup
            self.grid = self.mtzdata.grid_sampling
            self.hklinfo = self.mtzdata.hklinfo
        else:
            self.cell, self.spacegroup, self.grid = symmetry_from_model_metadata(model)

        '''
        Named triggers to simplify handling of changes on key events (e.g. live
        updating of maps/symmetry).
        '''
        self.triggers = TriggerSet()
        trigger_names = (
            'map box changed',  # Changed shape of box for map viewing
            'map box moved',    # Changed location of box for map viewing
            'atom box changed', # Changed shape of box for showing symmetry atoms
            'atom box moved',   # Changed location of box for showing symmetry atoms
        )
        for t in trigger_names:
            self.triggers.add_trigger(t)





        self._voxel_size = self.cell.dim / self.grid.dim

        self.show_nonpolar_H = show_nonpolar_H
        
        ref = model.bounds().center()

        # Convert the atoms to a format recognised by Clipper
        self._clipper_atoms = atom_list_from_sel(model.atoms)
        # A C++ object holding the symmetry operations required to pack one
        # unit cell relative to the current atomic model, along with fast
        # functions returning the symops necessary to pack a given box in
        # xyz space.
        self.unit_cell = clipper.Unit_Cell(ref,
                    self._clipper_atoms, self.cell, self.spacegroup, self.grid)
        
        # Container for managing all the symmetry copies
        self._sym_model_container = None
        
        # Container for drawing of special positions
        self._special_positions_model = None
        
        ###
        # LIVE ATOMIC SYMMETRY
        ###
        
        # Do we want to always have the reference model shown?
        self.sym_always_shows_reference_model = True
        # Trigger handler for live display of symmetry
        self._sym_handler = None
        # Centroid of the search space for symmetry equivalents
        self._sym_box_center = None
        # Half-width of the box in which we want to search
        self._sym_box_radius = atomic_symmetry_radius
        # Grid dimensions of the box
        self._sym_box_dimensions = None
        # Is the box already initialised?
        self._sym_box_initialized = False
        # Last grid coordinate for the cofr at which the symmetry was updated
        self._sym_last_cofr_grid = None
        # Factor defining the frequency of steps in the symmetry search
        # (larger number = more steps)
        self._sym_search_frequency = 2
        
        # Do we want to find and show atoms in the search radius at each iteration?
        self._live_atomic_symmetry = None
        self.live_atomic_symmetry = live_atomic_symmetry
        
        ###
        # REAL-SPACE MAPS / LIVE SCROLLING
        ###
        
        self.xmaps = None
        if calculate_maps:        
            if self.mtzdata is not None:
                if self.mtzdata.calculated_data is not None:
                        # We want to add the capability to generate maps from
                        # the atomic model and experimental data here, but for
                        # now we'll only add maps if amplitudes and phases were
                        # provided in the MTZ file.
                        self.xmaps = XmapSet(self.session, 
                            self.mtzdata.calculated_data, self, 
                            live_scrolling = live_map_scrolling,
                            display_radius = map_scrolling_radius)
                        self.add([self.xmaps])


    @property
    def live_atomic_symmetry_radius(self):
        '''
        Set the radius (in Angstroms) of the volume in which to show 
        symmetry atoms.
        '''
        return self._sym_box_radius

    @live_atomic_symmetry_radius.setter
    def live_atomic_symmetry_radius(self, radius):
        self._change_sym_box_radius(radius)

    @property
    def live_atomic_symmetry(self):
        '''
        Turn live display of symmetry atoms on or off.
        '''
        return self._live_atomic_symmetry
    
    @live_atomic_symmetry.setter
    def live_atomic_symmetry(self, switch):
        if switch:
            if not self._sym_box_initialized:
                self.initialize_symmetry_display(self.live_atomic_symmetry_radius)
            set_to_default_cartoon(self.session, self.master_model)
            for key, m in self.items():
                m.bonds.radii = DEFAULT_BOND_RADIUS
            self._start_live_atomic_symmetry()
        else:
            self._stop_live_atomic_symmetry()
        self._live_atomic_symmetry = switch
    

    @property
    def sym_model_container(self):
        if self._sym_model_container is None:
            self._sym_model_container = SymModels(self.session, self)
        return self._sym_model_container
        
    @property
    def display(self):
        return super().display
    
    @display.setter
    def display(self, switch):
        if switch:
            if self.live_atomic_symmetry:
                self._start_live_atomic_symmetry()
            if self.xmaps.live_scrolling:
                self.xmaps._start_live_scrolling()
        else:
            self._stop_live_atomic_symmetry()
            self.xmaps._stop_live_scrolling()
        Model.display.fset(self, switch)
        

    def items(self):
        return ((clipper.RTop_frac.identity(), self.master_model), *self.sym_model_container.items())

    def add_model_to_self(self, model):
        '''
        Transplant a model from
        '''

    def sym_select_within(self, coords, radius):
        '''
        Given an array of (x,y,z) coordinates, return a list of Atoms lists
        (one per symmetry equivalent molecule) covering all atoms within
        radius of any input coordinate.
        Args:
          coords:
            An (n * [x,y,z)) array of coordinates
          radius:
            Search radius in Angstroms
        '''
        c = numpy.empty((len(coords), 3), numpy.float32)
        c[:] = coords
        master_atoms = self.master_model.atoms
        master_coords = master_atoms.coords.astype(numpy.float32)
        grid_minmax = clipper.Util.get_minmax_grid(coords, self.cell, self.grid)
        pad = calculate_grid_padding(radius, self.grid, self.cell)
        grid_minmax += numpy.array((-pad, pad))
        min_xyz = clipper.Coord_grid(grid_minmax[0]).coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = grid_minmax[1] - grid_minmax[0]
        symops = self.unit_cell.all_symops_in_box(min_xyz, dim)
        symmats = symops.all_matrices_orth(self.cell, format = '3x4')
        target = [(c, Place().matrix.astype(numpy.float32))]
        search_list = []
        model_list = []
        for i, s in enumerate(symops):
            search_list.append((master_coords, symmats[i].astype(numpy.float32)))
            model_list.append(self.sym_model_container[s])
        i1, i2 = find_close_points_sets(search_list, target, radius)

        found = []
        for i, c in enumerate(i1):
            if len(c):
                found.append(model_list[i].atoms[c])
        return found

    def initialize_symmetry_display(self, radius = 14):
        '''
        Continually update the display to show all symmetry equivalents of
        the atomic model which enter a box of the given size, centred on the
        centre of rotation.

        Args:
          radius (float):
            The search volume is actually the smallest parallelepiped defined
            by the same angles as the unit cell, which will contain a sphere
            of the given radius.
        '''
        if self._sym_box_initialized:
            raise RuntimeError('''
              The symmetry box is already intialised for this structure. If you
              want to reset it, run sym_box_reset().
              ''')
        camera.camera(self.session, 'ortho')
        cofr.cofr(self.session, 'centerOfView', show_pivot = True)
        self.sym_always_shows_reference_model = True
        self._sym_box_radius = radius
        uc = self.unit_cell
        v = self.session.view
        c = self.cell
        g = self.grid
        self._sym_box_center = v.center_of_rotation
        self._sym_last_cofr_grid = clipper.Coord_orth(self._sym_box_center).coord_frac(c).coord_grid(g)
        box_corner_grid, box_corner_xyz = _find_box_corner(self._sym_box_center, c, g, radius)
        self._sym_box_dimensions = (numpy.ceil(radius / self._voxel_size * 2)).astype(int)
        self._update_sym_box(force_update = True)
        self._sym_box_initialized = True

    def _change_sym_box_radius(self, radius):
        dim = (numpy.ceil(radius / self._voxel_size * 2)).astype(int)
        self._sym_box_dimensions = dim
        self._sym_box_radius = radius

    def _update_sym_box(self, *_, force_update = False):
        v = self.session.view
        uc = self.unit_cell
        cofr = v.center_of_rotation
        cofr_grid = clipper.Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
        if not force_update:
            if self._sym_last_cofr_grid is not None:
                if cofr_grid == self._sym_last_cofr_grid:
                    return
        self._sym_last_cofr_grid = cofr_grid
        self._sym_box_center = cofr
        box_corner_grid, box_corner_xyz = _find_box_corner(
                  self._sym_box_center, self.cell,
                  self.grid, self._sym_box_radius)
        symops = uc.all_symops_in_box(box_corner_xyz,
                  self._sym_box_dimensions.astype(numpy.int32), self.sym_always_shows_reference_model,
                  sample_frequency = self._sym_search_frequency)
        l1 = []
        search_entry = [(numpy.array([cofr], numpy.float32), Place().matrix.astype(numpy.float32))]
        coords = self.master_model.atoms.coords.astype(numpy.float32)
        display_mask = numpy.array([False]*len(coords))
        if not self.show_nonpolar_H:
            h_mask = self.master_model.atoms.idatm_types != 'HC'
        for s in symops:
            this_model = self.sym_model_container[s]
            this_set = (coords, s.rtop_orth(self.cell).mat34.astype(numpy.float32))
            l1.append(this_set)
        i1, i2 = find_close_points_sets(l1, search_entry, self.live_atomic_symmetry_radius)
        for i, indices in enumerate(i1):
            this_model = self.sym_model_container[symops[i]]
            if len(indices):
                display_mask[:] = False
                a = this_model.atoms
                display_mask[a.indices(a[indices].residues.atoms)] = True
                if not self.show_nonpolar_H:
                    display_mask = numpy.logical_and(display_mask, h_mask)
                a.displays = display_mask
                #if this_model is not self.master_model:
                    #a.residues.ribbon_displays = display_mask
                this_model.display = True
            else:
                if this_model is not self.master_model:
                    this_model.display = False
        for s, m in self.sym_model_container.items():
            if s not in symops:
                m.display = False

    def show_large_scale_symmetry(self, box_width = 200):
        '''
        Show the model symmetry over a large volume by tiling the current
        representation of the master model. NOTE: this will automatically
        stop live atomic symmetry updating.
        '''
        self.live_atomic_symmetry = False
        box_center = self.master_model.bounds().center()
        uc = self.unit_cell
        dim = (numpy.ceil(box_width / self._voxel_size)).astype(int)
        box_corner_grid, box_corner_xyz = _find_box_corner(
                  box_center, self.cell,
                  self.grid, box_width/2)
        symops = uc.all_symops_in_box(box_corner_xyz, dim, True)
        num_symops = len(symops)
        sym_matrices = symops.all_matrices_orth(self.cell, format = '3x4')
        self.master_model.positions = Places(place_array=sym_matrices)

    def hide_large_scale_symmetry(self, restart_live_atomic_symmetry = True):
        self.master_model.position = Place()
        self.live_atomic_symmetry = restart_live_atomic_symmetry

    def _start_live_atomic_symmetry(self):
        self.hide_large_scale_symmetry(restart_live_atomic_symmetry = False)
        set_to_default_cartoon(self.session)
        if self._sym_handler is None:
            self._sym_handler = self.session.triggers.add_handler('new frame', self._update_sym_box)

    def _stop_live_atomic_symmetry(self):
        if self._sym_handler is not None:
            self.session.triggers.remove_handler(self._sym_handler)
            self._sym_handler = None
        for key, m in self.sym_model_container.items():
            m.display = False

    def isolate_and_cover_selection(self, atoms, include_surrounding_residues = 5,
                        show_context = 5, mask_radius = 3, hide_surrounds = True, focus = True):
        '''
        Expand the map to cover a given atomic selection, then mask it to
        within a given distance of said atoms to reduce visual clutter. Adjust
        the atomic visualisation to show only the selected atoms, plus an
        optional surrounding buffer zone.
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
          hide_surrounds (bool):
            If true, all residues outside the selection region will be hidden
          focus (bool):
            If true, the camera will be moved to focus on the selection
        '''
        # If we're in live mode, turn it off
        self.live_atomic_symmetry = False
        self.xmaps.live_scrolling = False
        orig_atoms = atoms
        atoms = atoms.residues.atoms
        coords = atoms.coords
        if include_surrounding_residues > 0:
            atoms = concatenate(
              self.sym_select_within(
                  coords, include_surrounding_residues)).residues.atoms
            coords = atoms.coords
        context_atoms = None
        if show_context > 0:
            context_atoms = concatenate(
              self.sym_select_within(
                  coords, show_context)).residues.atoms.subtract(atoms)
        pad = calculate_grid_padding(mask_radius, self.grid, self.cell)
        box_bounds_grid = clipper.Util.get_minmax_grid(coords, self.cell, self.grid) \
                                + numpy.array((-pad, pad))
        self.xmaps.set_box_limits(box_bounds_grid)

        self.xmaps._surface_zone.update(mask_radius, atoms, None)
        self.xmaps._reapply_zone()
        if context_atoms is None:
            found_models = atoms.unique_structures
        else:
            found_models = concatenate((context_atoms, atoms)).unique_structures
        self.master_model.bonds.radii = 0.05
        for key, m in self.items():
            if m not in found_models:
                m.display = False
            else:
                m.display = True
                m.atoms.displays = False
                m.residues.ribbon_displays = False
        self.master_model.atoms[numpy.in1d(
            self.master_model.atoms.names, numpy.array(
                ['N','C','CA']))].displays = True
        if not self.show_nonpolar_H:
            atoms = atoms.filter(atoms.idatm_types != 'HC')
        atoms.displays = True
        atoms.inter_bonds.radii = 0.2
        atoms.residues.ribbon_displays = True
        if context_atoms is not None:
            if not self.show_nonpolar_H:
                context_atoms = context_atoms.filter(context_atoms.idatm_types != 'HC')
            context_atoms.displays = True
            context_atoms.inter_bonds.radii = 0.1
        if focus:
            self.session.view.view_all(atoms.scene_bounds, 0.2)
        # return the original selection in case we want to re-run with modified settings
        return orig_atoms

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
        m = self._special_positions_model
        if not (m is None or m.deleted):
            # Just show the existing model
            m.display = True
            return
        model = self.master_model

        ref = model.bounds().center().astype(float)
        frac_coords = clipper.Coord_orth(ref).coord_frac(self.cell).uvw
        if offset is None:
            offset = numpy.array([0,0,0],int)
            
        positions = []
        colors = []
        rgba_corner = numpy.array([255,0,255,128],numpy.int32)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset\
                        + self.unit_cell.min.coord_frac(self.grid).uvw

        corners = []
        for c in corners_frac:
            co = clipper.Coord_frac(c).coord_orth(self.cell)
            positions.append(Place(axes=numpy.identity(3)*4, origin=co.xyz))
            colors.append(rgba_corner)
    
        
        m = self.special_positions_model = Model('Special Positions',self.session)
        xmap = self.xmaps.child_models()[0].xmap
        uc = self.unit_cell
        spc = numpy.array(xmap.special_positions_unit_cell_xyz(uc, offset))
        d = Drawing('points')
        sphere = numpy.array(sphere_geometry(80))
        sphere[0]*=0.25
        d.vertices, d.normals, d.triangles = sphere

        if len(spc):
            
            coords = spc[:,0:3]
            multiplicity = spc[:,3].astype(int)
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.int32)
            rgba_3fold = numpy.array([0,255,255,255],numpy.int32)
            rgba_4fold = numpy.array([255,255,0,255],numpy.int32)
            rgba_6fold = numpy.array([255,0,0,255],numpy.int32)

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

        d.positions = Places(positions)
        d.colors = numpy.array(colors)
        m.add_drawing(d)
        model.parent.add([m])
        m.display = True




class SymModels(defaultdict):
    '''
    Handles creation, destruction and organisation of symmetry copies
    of an atomic model. Uses the Clipper RTop_frac object for the given
    symop as the dict key. If the key is not found, automatically creates
    a copy of the master model, sets colours, applies the Place transform
    for the symop, and adds the model to the session.
    '''
    def __init__(self, session, parent):
        '''
        Just create an empty dict.
        Args:
          session:
            The ChimeraX session.
          parent:
            The AtomicCrystalStructure describing the master model and symmetry
        '''
        self.session = session
        self.parent = parent
        self.master = parent.master_model

        # Add a sub-model to the master to act as a container for the
        # symmetry copies
        self._sym_container = None

    @property
    def sym_container(self):
        if self._sym_container is None or self._sym_container.deleted:
            self._sym_container = Model('symmetry equivalents', self.session)
            self.parent.add([self._sym_container])
            self.clear()
        return self._sym_container

    def __missing__(self, key):
        if type(key) is not clipper.RTop_frac:
            raise TypeError('Key must be a clipper.RTop_frac!')
        thisplace = Place(matrix=key.rtop_orth(self.parent.cell).mat34)
        if not thisplace.is_identity():
            thismodel = self.master.copy(name=key.format_as_symop)
            atoms = thismodel.atoms
            #thismodel.position = thisplace
            atoms.coords = thisplace.moved(atoms.coords)
            atom_colors = atoms.colors
            atom_colors[:,0:3] = (self.master.atoms.colors[:,0:3].astype(float)*0.6).astype(numpy.uint8)
            atoms.colors = atom_colors
            ribbon_colors = thismodel.residues.ribbon_colors
            ribbon_colors[:,0:3] = (self.master.residues.ribbon_colors[:,0:3].astype(float)*0.6).astype(numpy.uint8)
            thismodel.residues.ribbon_colors = ribbon_colors
            self.sym_container.add([thismodel])
            set_to_default_cartoon(self.session, thismodel)
            self[key] = thismodel
            thismodel.display = False
            return thismodel
        return self.master

    def __getitem__(self, key):
        s = self.sym_container
        m = super(SymModels, self).__getitem__(key)
        return m

def set_to_default_cartoon(session, model = None):
    '''
    Adjust the ribbon representation to provide information without
    getting in the way.
    '''
    if model is None:
        atoms = None
    else:
        arg = atomspec.AtomSpecArg('thearg')
        atoms = arg.parse('#' + model.id_string(), session)[0]
    cartoon.cartoon(session, atoms = atoms, suppress_backbone_display=False)
    cartoon.cartoon_style(session, atoms = atoms, width=0.4, thickness=0.1, arrows_helix=True, arrow_scale = 2)


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
    '''
    Essentially a copy of chimerax.core.surface.zone.surface_zone, but uses
    find_close_points_sets to eke a little extra performance
    '''
    vlist = []
    dlist = []
    ident_matrix = Place().matrix.astype(numpy.float32)
    search_entry = [(numpy.array(points, numpy.float32), Place().matrix.astype(numpy.float32))]
    for m in models:
        for d in m.child_drawings():
            if not d.display:
                continue
            if d.vertices is not None:
                dlist.append(d)
                vlist.append((d.vertices.astype(numpy.float32), ident_matrix))
  
    i1, i2 = find_close_points_sets(vlist, search_entry, distance)
  
    for vp, i, d in zip(vlist, i1, dlist):
        v = vp[0]
        nv = len(v)
        mask = numpy.zeros((nv,), numpy.bool)
        numpy.put(mask, i, 1)
        t = d.triangles
        if t is None:
            return
        tmask = numpy.logical_and(mask[t[:,0]], mask[t[:,1]])
        numpy.logical_and(tmask, mask[t[:,2]], tmask)
        d.triangle_mask = tmask


class XmapSet(Model):
    '''
    Each crystal will typically have multiple maps associated with it -
    the standard 2mFo-DFc and mFo-DFc maps, for a start, but also
    potentially sharpened/smoothed versions of these, anomalous difference
    maps, omit maps, etc. etc. etc.

    XmapSet is designed as a class to contain and organise these, control
    their display and recalculation, etc.
    '''
    
    STANDARD_LOW_CONTOUR = numpy.array([1.5])
    STANDARD_HIGH_CONTOUR = numpy.array([2.0])
    STANDARD_DIFFERENCE_MAP_CONTOURS = numpy.array([-3.0, 3.0])
    
    DEFAULT_MESH_MAP_COLOR = [0,1.0,1.0,1.0] # Solid cyan
    DEFAULT_SOLID_MAP_COLOR = [0,1.0,1.0,0.4] # Transparent cyan
    DEFAULT_DIFF_MAP_COLORS = [[1.0,0,0,1.0],[0,1.0,0,1.0]] #Solid red and green

    
    def __init__(self, session, datasets, crystal,
                 live_scrolling = True, display_radius = 12,
                 atom_selection = None, padding = 3):
        '''
        Args:
            session:
                The ChimeraX session
            datasets:
                An iterable of ReflectionData_Calc objects
            parent:
                The CrystalStructure object this object will be subordinate to.
            live_scrolling:
                If True, the maps will be initialised in live scrolling
                mode, displaying a sphere of density centred on the 
                centre of rotation. This option takes precedence over
                atom_selection. 
            display_radius:
                The radius (in Angstroms) of the display sphere used in
                live scrolling mode.
            atom_selection:
                If live_scrolling is False, this argument must provide a
                selection of atoms around which the maps will be masked.
            padding:
                The radius (in Angstroms) of the mask surrounding each
                atom in atom_selection.
        '''
        if not live_scrolling and not atom_selection:
            raise TypeError('''
            If live_scrolling is False, you must provide a set of atoms
            to mask the maps to!
            ''')
        Model.__init__(self, 'Real-space maps', session)
        self.crystal = crystal
        #############
        # Variables involved in handling live redrawing of maps in a box
        # centred on the cofr
        #############

        # ChimeraX session.triggers handler for live box update
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
        self._last_box_center_grid = None
        # Minimum corner of the box in (x,y,z) coordinates. This must
        # correspond to the grid coordinate in self._box_corner_grid
        self._box_corner_xyz = None
        # Minimum corner of the box in grid coordinates
        self._box_corner_grid = None
        
        self.live_scrolling = live_scrolling

        if not live_scrolling:
            # Get the initial box parameters based on atom_selection and padding
            self._box_corner_grid, self._box_corner_xyz, self._box_dimensions = \
                _get_bounding_box(atom_selection.coords, padding, self.grid, self.cell)
            self._surface_zone.update(padding, atoms = atom_selection)
        
        self._box_initialized = True
        
        for dataset in datasets:
            print('Working on dataset: {}'.format(dataset.name))
            self.add_map_handler(dataset)
        
        # Apply the surface mask
        # self._reapply_zone()


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
        return self.hklinfo.resolution
    
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
        cofr = v.center_of_rotation
        dim = self._box_dimensions = \
            2 * calculate_grid_padding(radius, self.grid, self.cell)
        self._box_corner_grid, self._box_corner_xyz = _find_box_corner(
            cofr, self.cell, self.grid, radius)
        self.crystal.triggers.activate_trigger('map box changed', 
            (self._box_corner_xyz, self._box_corner_grid, dim))
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
        self._live_scrolling = switch

    def _start_live_scrolling(self):
        if self._box_update_handler is None:
            self._box_update_handler = self.session.triggers.add_handler(
                'new frame', self.update_box)
    
    def _stop_live_scrolling(self):
        if self._box_update_handler is not None:
            self.session.triggers.remove_handler(self._box_update_handler)
            self._box_update_handler = None
        
    
    
    
    def __getitem__(self, name):
        '''Get one of the child maps by name.'''
        for m in self.child_models():
            if m.name == name:
                return m
        raise KeyError('No map with that name!')
    

    def set_box_limits(self, minmax):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum grid coordinates. Automatically turns off live scrolling.
        '''
        self.live_scrolling = False
        cmin = clipper.Coord_grid(minmax[0])
        cmin_xyz = cmin.coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = (minmax[1]-minmax[0])
        self.crystal.triggers.activate_trigger('map box changed', 
            (cmin_xyz, cmin, dim))
    
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
        box_max_grid = (uc.max+clipper.Coord_grid([2,2,2])).uvw
        minmax = [box_min_grid, box_max_grid]
        self.set_box_limits(minmax)
        # Tile by the desired number of cells
        places = []
        grid_dim = self.grid.dim
        nu, nv, nw = nuvw
        ou, ov, ow = offset
        for i in range(ou, nu+ou):
            for j in range(ov, nv+ov):
                for k in range(ow, nw+ow):
                    thisgrid = clipper.Coord_grid(numpy.array([i,j,k])*grid_dim)
                    thisorigin = thisgrid.coord_frac(self.grid).coord_orth(self.cell).xyz
                    places.append(Place(origin = thisorigin))
        self.positions = Places(places)
                    
        
        
        
    
    
    def add_map_handler(self, dataset, is_difference_map = None, 
                color = None, style = None, contour = None):
        '''
        Add a new XmapHandler based on the given reflections and phases.
        Args:
            dataset: 
                a ReflectionData_Calc object.
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
        data = dataset.data
        new_xmap = clipper.Xmap(self.spacegroup, self.cell, self.grid, name = dataset.name, hkldata = data)
        if is_difference_map is None:
            is_difference_map = dataset.is_difference_map
        new_xmap.is_difference_map = is_difference_map
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for 
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        new_handler = XmapHandler(self.session, self.crystal, dataset.name, new_xmap, 
            self._box_corner_xyz, self._box_corner_grid, self._box_dimensions)
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
        contour = contour * new_xmap.sigma
        self.add([new_handler])
        new_handler.set_representation(style)
        new_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': True})
        new_handler.update_surface()
        new_handler.show()
    
    def update_box(self, *_, force_update = False):
        '''Update the map box to surround the current centre of rotation.'''
        v = self.session.view
        cofr = v.center_of_rotation
        self._box_center = cofr
        cofr_grid = clipper.Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
        if not force_update:
            if self._last_box_center_grid is not None:
                if cofr_grid == self._last_box_center_grid:
                    return
        self._last_box_center_grid = cofr_grid      
        box_corner_grid, box_corner_xyz = _find_box_corner(cofr, self.cell, self.grid, self.display_radius)
        self.crystal.triggers.activate_trigger('map box moved', (box_corner_xyz, box_corner_grid, self._box_dimensions))
        self._surface_zone.update(self.display_radius, coords = numpy.array([cofr]))
        self._reapply_zone()
     
    def _reapply_zone(self):
        '''
        Reapply any surface zone applied to the volume after changing 
        position or contour level.
        '''
        coords = self._surface_zone.all_coords
        radius = self._surface_zone.distance
        if coords is not None:
            surface_zones(self.child_models(), coords, radius)
            
    def delete(self):
        self.live_scrolling = False
        super(XmapSet, self).delete()

    
def calculate_grid_padding(radius, grid, cell):
    '''
    Calculate the number of grid steps needed on each crystallographic axis
    in order to capture at least radius angstroms in x, y and z.
    '''
    corner_mask = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,0,0],[1,1,1]])
    corners = (corner_mask * radius).astype(float)
    grid_upper = numpy.zeros([8,3], numpy.int)
    grid_lower = numpy.zeros([8,3], numpy.int)
    for i, c in enumerate(corners):
        co = clipper.Coord_orth(c)
        cm = co.coord_frac(cell).coord_map(grid).uvw
        grid_upper[i,:] = numpy.ceil(cm).astype(int)
        grid_lower[i,:] = numpy.floor(cm).astype(int)
    return grid_upper.max(axis=0) - grid_lower.min(axis=0)

def _find_box_corner(center, cell, grid, radius = 20):
    '''
    Find the bottom corner (i.e. the origin) of a rhombohedral box
    big enough to hold a sphere of the desired radius.
    '''
    radii_frac = clipper.Coord_frac(radius/cell.dim)
    center_frac = clipper.Coord_orth(center).coord_frac(cell)
    bottom_corner_grid = center_frac.coord_grid(grid) \
                - clipper.Coord_grid(calculate_grid_padding(radius, grid, cell))
    bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
    return bottom_corner_grid, bottom_corner_orth.xyz

def _get_bounding_box(coords, padding, grid, cell):
    '''
    Find the minimum and maximum grid coordinates of a box which will 
    encompass the given (x,y,z) coordinates plus padding (in Angstroms).
    '''
    grid_pad = calculate_grid_padding(padding, grid, cell)
    box_bounds_grid = clipper.Util.get_minmax_grid(coords, cell, grid)\
                        + numpy.array((-grid_pad, grid_pad))
    box_origin_grid = box_bounds_grid[0]
    box_origin_xyz = clipper.Coord_grid(box_origin_grid).coord_frac(grid).coord_orth(cell)
    dim = box_bounds_grid[1] - box_bounds_grid[0]
    return [box_origin_grid, box_origin_xyz, dim]


class XmapHandler(Volume):
    '''
    An XmapHandler is in effect a resizable window into a periodic 
    crystallographic map. The actual map data (a clipper Xmap object) is
    held within, and filled into the XmapWindow.data array as needed.
    Methods are included for both live updating (e.g. tracking and filling
    a box centred on the centre of rotation) and static display of a 
    given region.
    '''
    def __init__(self, session, crystal, name, xmap, origin, grid_origin, dim):
        '''
        Args:
            sesssion:
                The ChimeraX session
            crystal:
                The CrystalStructure object this belongs to
            name:
                A descriptive name for this map
            xmap:
                A clipper.Xmap
            origin:
                The (x,y,z) coordinates of the bottom left corner of the 
                volume.
            grid_origin:
                The (u,v,w) integer grid coordinates corresponding to
                origin.
            dim:
                The shape of the box in (u,v,w) grid coordinates.
        '''
        self.box_params = (origin, grid_origin, dim)
        self.xmap = xmap
        self.crystal = crystal
        darray = self._generate_and_fill_data_array(origin, grid_origin, dim)
        Volume.__init__(self, darray, session)
        
        self.is_difference_map = xmap.is_difference_map
        self.name = name
        self.initialize_thresholds()
        
        # If the box shape changes while the volume is hidden, the change
        # will not be applied until it's shown again.
        self._needs_update = True
        self.show()
        self._box_shape_changed_cb_handler = self.crystal.triggers.add_handler(
            'map box changed', self._box_changed_cb)
        self._box_moved_cb_handler = self.crystal.triggers.add_handler(
            'map box moved', self._box_moved_cb)

        
    
    def show(self):
        if self._needs_update:
            self._swap_volume_data(self.box_params, force_update = True)
            self._needs_update = False
        else:
            # Just set the origin and fill the box with the data for 
            # the current location
            origin, grid_origin, ignore = self.box_params
            self._fill_volume_data(self.data.array, grid_origin)
        super(XmapHandler, self).show()

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
        return self.hklinfo.resolution
    
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
    def _surface_zone(self):
        return self.parent._surface_zone

    def _box_changed_cb(self, name, params):
        self.box_params = params
        self._needs_update = True
        if not self.display:
            # No sense in wasting cycles on this if the volume is hidden.
            # We'll just store the params and apply them when we show the
            # volume.
            # NOTE: this means we need to over-ride show() to ensure 
            # it's updated before re-displaying.
            return
        self._swap_volume_data(params)
        self.data.values_changed()
        self.show()
        
    def _box_moved_cb(self, name, params):
        self.box_params = params
        if not self.display:
            return
        self.data.set_origin(params[0])
        self._fill_volume_data(self.data.array, params[1])
        self.data.values_changed()
        self.update_surface()
    
    def delete(self):
        if self._box_shape_changed_cb_handler is not None:
            self.triggers.remove_handler(self._box_shape_changed_cb_handler)
            self._box_shape_changed_cb_handler = None
        if self._box_moved_cb_handler is not None:
            self.crystal.triggers.remove_handler(self._box_moved_cb_handler)
            self._box_moved_cb_handler = None
        super(XmapHandler, self).delete()
    
    
        
    
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
        self.new_region((0,0,0), darray.size)
    
    def _generate_and_fill_data_array(self, origin, grid_origin, dim):
        dim = dim[::-1]
        data = numpy.empty(dim, numpy.double)
        self._fill_volume_data(data, grid_origin)
        darray = Array_Grid_Data(data, origin = origin,
            step = self.voxel_size, cell_angles = self.cell.angles_deg)
        return darray
    
        
    def _fill_volume_data(self, target, start_grid_coor):
        #shape = (numpy.array(target.shape)[::-1] - 1)
        #end_grid_coor = start_grid_coor + clipper.Coord_grid(shape)
        #self.data.set_origin(origin_xyz)
        xmap = self.xmap
        xmap.export_section_numpy(start_grid_coor, target = target,  order = 'C', rot = 'zyx')
