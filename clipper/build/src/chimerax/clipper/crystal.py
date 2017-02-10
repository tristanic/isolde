import numpy
from . import clipper
from .clipper_mtz import Clipper_MTZ
from .data_tree import db_levels, DataTree
class Xtal_Project:
  '''
  The master object for handling a crystallographic data set within
  ChimeraX. Contains a Clipper_MTZ object holding all reciprocal-space
  data for one or more crystals, any maps generated for each crystal,
  and all ChimeraX-specific functions for working with them.
  '''
  def __init__(self, session, db, name):
    self.session = session
    self.master_db = db
    from .clipper_mtz import Clipper_MTZ
    self.data = db['Experiment'][name] = Clipper_MTZ(parent=db['Experiment'])
    
  def load_data(self, filename):
    self.data.load_hkl_data(filename)
  
  def add_maps(self, crystal_name):
    maps = self.data[crystal_name]['maps'] = \
          Map_set(self.session, self.data[crystal_name]['dataset'])
    for key in self.data[crystal_name]['dataset']['F_Phi'].keys():
      maps.generate_map_from_f_phi(key)
    return maps

DEFAULT_SOLID_MAP_COLOR = [0,1.0,1.0,0.6] # Transparent cyan
DEFAULT_MESH_MAP_COLOR = [0,1.0,1.0,1.0] # Solid cyan
DEFAULT_DIFF_MAP_COLORS = [[1.0,0,0,1],[1.0,1.0,0,1]] #Solid red and yellow


class Map_set:
  '''
  Each crystal dataset will usually have multiple maps associated with it.
  We need an object to organise these, to control their display, and to
  re-calculate them as necessary.  
  '''
  
  def __init__(self, session, dataset):
    '''
    Parameters:
      session: the current ChimeraX session
      dataset: a node from the Clipper_MTZ tree corresponding to one
               crystallographic dataset
    '''
    self.session = session
    if not isinstance(dataset, DataTree) or dataset.level != db_levels.DATASET:
      raise TypeError('dataset should be the fourth level of a Clipper_MTZ tree!')
    self._data = dataset
    # The Clipper HKL_info object holding the (h,k,l) arrays and cell/symmetry information
    hkl = self.hklinfo = dataset.find_ancestor(db_levels.CRYSTAL_SET)['hkl']
    sg = self.spacegroup = hkl.spacegroup
    cell = self.cell = hkl.cell
    res = self.resolution = hkl.resolution
    grid = self.grid = clipper.Grid_sampling(sg, cell, res)
    
    self.voxel_size = cell.dim / grid.dim
    self._voxel_size_frac = 1 / grid.dim
    
    
    # List of Xmap objects generated and handled by this object
    self.maps = []
    
    # Atomic model associated with this object
    self.atomic_model = None
    
    
    #############
    # Variables involved in handling live redrawing of maps in a box
    # centred on the cofr
    #############
    
    # FIXME: The model box and its parameters should really be handled
    # at a higher level - we really just want one box per session, into
    # which we can pour multiple maps from different crystals as necessary
    
    from chimerax.core.models import Model
    self._master_box_model = Model('Xmap display', self.session)
    
    self._standard_contour = numpy.array([1.5])
    self._standard_difference_map_contours = numpy.array([-3.0, 3.0])
    # Default "radius" (actually half-width of the rhombohedron) of
    # the box in which the maps will be drawn.
    self._box_radius = None
    # Box dimensions in grid points (depends on resolution)
    self._box_dimensions = None
    # Centre point of the box (will typically be set to the centre of
    # rotation).
    self._box_center = None
    # Last grid coordinate of the box centre. We only need to update the
    # map if the new centre changes
    self._last_cofr_grid = None
    # Vector mapping the box centre to its origin
    self._box_origin_offset = None
    # ChimeraX Volume objects to draw the map into
    self._volumes = []
    # Array_Grid_Data objects held by the Volume objects
    self._array_grid_data = []
    # The actual numpy arrays to send the map data to
    self._volume_data = []
    # session.triggers handler for live box update
    self._box_handler = None
    # Is the box already initialised?
    self._box_initialized = False
    # Distance any axis of the cofr has to move before triggering an update.
    # Since we're working on Clipper's integer grid, we only need trigger
    # an update once the cofr moves more than half a voxel in any direction
    # This is a bit of a fudge given that the voxels aren't necessarily on
    # orthogonal axes, but it's a useful compromise in the interest of speed.
    v = self.session.view
    cofr = v.center_of_rotation
    self._cofr_eps = self.session.main_view.pixel_size(cofr)

    
    ################
    # Variables involved in handling other useful annotations
    ################
    
    # Model object to hold Drawings defining the special positions
    self._special_positions_model = None
    
    self._crosshairs = None
    
    self._show_crosshairs = True

  @property
  def show_crosshairs(self):
    return self._show_crosshairs
  
  @show_crosshairs.setter
  def show_crosshairs(self, state):
    if state:
      if self._crosshairs is None or self._crosshairs.deleted:
        if self._box_handler:
          self._crosshairs = draw_crosshairs()
          self._master_box_model.add([self._crosshairs])
        else:
          print('Please initialise map drawing first!')
          return
    self._show_crosshairs = state
    self._crosshairs.visible = state
  
  def generate_map_from_f_phi (self, data_key):
    name = data_key
    data = self._data['F_Phi'][data_key]
    xmap = clipper.Xmap(self.spacegroup, self.cell, self.grid, name = data_key)
    xmap.fft_from(data)
    # FIXME: Ugly fudge for now to find difference maps
    if '2' not in data_key:
      xmap.is_difference_map = True
    self.maps.append(xmap)

  ################
  # Live-updating map box
  ################
  
  def initialize_box_display(self, radius = 15, pad = 2):
    '''
    Generate a Volume big enough to hold a sphere of the given radius,
    plus an optional padding of pad voxels on each side. The session
    view will be changed to use an orthographic camera, with the centre
    of rotation updating to always remain at the centre of the screen.
    The volume will automatically track the centre of rotation and 
    update its position and contents to reflect the local density.
    '''
    if self.maps is None or not len(self.maps):
      print('You must generate at least one map first!')
      return
    if self._box_initialized:
      print(''' The live map box is already initialised for this crystal.
                If you want to reset it, run box_reset().''')
      return
    from chimerax.core.map.data import Array_Grid_Data
    from chimerax.core.map import Volume    
    from chimerax.core.commands import camera, cofr
    camera.camera(self.session, 'ortho')
    cofr.cofr(self.session, 'centerOfView')
    self._box_pad = pad
    self._box_radius = radius
    v = self.session.view
    c = self.cell
    g = self.grid
    self._box_center = v.center_of_rotation
    self._last_cofr_grid = clipper.Coord_orth(self._box_center).coord_frac(c).coord_grid(g)
    if self.show_crosshairs:
      self._crosshairs = draw_crosshairs(self._box_center)
      self._master_box_model.add([self._crosshairs])

    box_corner_grid, box_corner_xyz = self._find_box_corner(self._box_center, radius, pad)
    self.box_origin_offset = box_corner_xyz - self._box_center
    self._box_dimensions = (numpy.ceil(radius / self.voxel_size * 2)+pad*2).astype(int)[::-1]
    for m in self.maps:
      data = numpy.empty(self._box_dimensions, numpy.double)
      self._volume_data.append(data)
      grid_data = Array_Grid_Data(data, origin = box_corner_xyz,
          step = self.voxel_size, cell_angles = c.angles_deg)
      self._array_grid_data.append(grid_data)
      volume = Volume(grid_data, self.session)
      volume.name = m.name
      self._volumes.append(volume)
      self._fill_volume_data(m, data, box_corner_grid)
      volume.initialize_thresholds()
      self._master_box_model.add([volume])
      if not m.is_difference_map:
        contour_val = [self._standard_contour * m.sigma]
        colorset = {'surface_colors': [DEFAULT_SOLID_MAP_COLOR]}
      else:
        contour_val = self._standard_difference_map_contours * m.sigma
        colorset = {'surface_colors': DEFAULT_DIFF_MAP_COLORS}
        volume.set_representation('mesh')
      volume.set_parameters(**{'cap_faces': False})
      self.change_contour(volume, contour_val)
      volume.set_parameters(**colorset)
      volume.update_surface()
      volume.show()
    self.session.models.add([self._master_box_model])
    self.update_box(force_update = True)
    self._box_go_live()

  def change_box_radius(self, radius, pad=1):
    self._box_go_static()
    v = self.session.view
    cofr = v.center_of_rotation
    self._box_radius = radius
    dim = (numpy.ceil(radius / self.voxel_size * 2)+pad*2).astype(int)[::-1]
    box_corner_grid, box_corner_xyz = self._find_box_corner(cofr, radius, pad)
    from chimerax.core.map.data import Array_Grid_Data
    for i, volume in enumerate(self._volumes):
      data = numpy.empty(dim, numpy.double)
      self._volume_data[i] = data    
      darray = Array_Grid_Data(data, origin = box_corner_xyz,
        step = self.voxel_size, cell_angles = self.cell.angles_deg)
      self._array_grid_data[i] = darray
      volume.replace_data(darray)
      volume.new_region((0,0,0), darray.size)
    self._box_pad = pad
    self._box_dimensions = dim
    self.update_box(force_update=True)
    self._box_go_live()
    
    
    
  def update_box(self, *_, force_update = False):
    from chimerax.core.geometry import Place
    from chimerax.core.surface import zone
    v = self.session.view
    cofr = v.center_of_rotation
    self._cofr_eps = self.session.main_view.pixel_size(cofr)
    # We need to redraw the crosshairs if the cofr moves by a pixel...
    if not force_update:
      if numpy.all(abs(self._box_center - cofr) < self._cofr_eps):
        return
    if self.show_crosshairs:
      self._crosshairs.position = Place(origin=cofr)
    # ... and redraw the box if it moves by a grid point
    self._box_center = cofr
    cofr_grid = clipper.Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
    if not force_update:
      if cofr_grid == self._last_cofr_grid:
        return
    self._last_cofr_grid = cofr_grid      
    box_corner_grid, box_corner_xyz = self._find_box_corner(cofr, self._box_radius, self._box_pad)
    for i, volume in enumerate(self._volumes):
      self._array_grid_data[i].set_origin(box_corner_xyz)
      self._fill_volume_data(self.maps[i], self._volume_data[i], box_corner_grid)
      volume.update_surface()
      zone.surface_zone(volume, [cofr], self._box_radius)
    
  def change_contour(self, volume, contour_vals):
    from chimerax.core.map import volumecommand
    volume.set_parameters(**{'surface_levels': contour_vals})
    volume.update_surface()
    

  def _fill_volume_data(self, xmap, target, start_grid_coor):
    shape = (numpy.array(target.shape)[::-1] - 1).tolist()
    end_grid_coor = start_grid_coor + clipper.Coord_grid(shape)
    xmap.export_section_numpy(target, start_grid_coor, end_grid_coor, 'C')
        
  
  def _find_box_corner(self, center, radius = 20, pad = 1):
    '''
    Find the bottom corner (i.e. the origin) of a rhombohedral box
    big enough to hold a sphere of the desired radius, and padded by
    pad voxels in each dimension.
    '''
    cell = self.cell
    grid = self.grid
    voxel_size = self.voxel_size
    radii_in_voxels = radius / voxel_size
    radii_frac = clipper.Coord_frac(radii_in_voxels * self._voxel_size_frac)
    center_ortho = clipper.Coord_orth(center)
    center_frac = center_ortho.coord_frac(cell)
    bottom_corner_frac = center_frac - radii_frac
    bottom_corner_grid = bottom_corner_frac.coord_grid(grid) -\
              clipper.Coord_grid([pad,pad,pad])
    bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
    return bottom_corner_grid, bottom_corner_orth.xyz
    
  def _box_go_live(self):
    if self._box_handler is None:
      self._box_handler = self.session.triggers.add_handler('new frame', self.update_box)
  
  def _box_go_static(self):
    if self._box_handler is not None:
      self.session.triggers.remove_handler(self._box_handler)
      self._box_handler = None


def draw_crosshairs(origin):
  from chimerax.core.surface.shapes import cylinder_geometry
  from chimerax.core.models import Drawing
  from chimerax.core.geometry import Place, Places
  d = Drawing('axis')
  axes = Drawing('axes')
  d.vertices, d.normals, d.triangles = cylinder_geometry(radius = 0.05, height = 0.5)
  p = []
  p.append(Place())
  p.append(Place(axes = [[0,0,1],[0,1,0],[-1,0,0]]))
  p.append(Place(axes = [[1,0,0],[0,0,-1],[0,1,0]]))
  c = [[255,0,0,255],[0,255,0,255],[0,0,255,255]]
  d.positions = Places(p)
  d.colors = c
  axes.add_drawing(d)
  axes.position = Place(origin=origin)
  return axes
  
'''

                  
  def draw_unit_cell_and_special_positions(self, model, offset = None):
    from chimerax.core.models import Model, Drawing
    from chimerax.core.geometry import Place, Places
    from chimerax.core.surface.shapes import sphere_geometry
    import copy
    
    ref = model.bounds().center().astype(float)
    frac_coords = Coord_orth(ref).coord_frac(self.cell()).uvw
    if offset is None:
      offset = numpy.array([0,0,0],int)
    corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset
    corners = []
    for c in corners_frac:
      cf = Coord_frac(c).coord_orth(self.cell())
      corners.append(cf.xyz)
    m = self.special_positions_model
    
    if m is None or m.deleted:
      m = self.special_positions_model = Model('Special Positions',self.session)
      spc = numpy.array(self.special_positions_unit_cell_xyz(offset))
      coords = spc[:,0:3]
      multiplicity = spc[:,3].astype(int)
      sphere = numpy.array(sphere_geometry(80))
      sphere[0]*=0.5
      scale_2fold = numpy.identity(3)
      scale_3fold = numpy.identity(3)* 1.5
      scale_4fold = numpy.identity(3)* 2
      scale_6fold = numpy.identity(3)* 3
      rgba_2fold = numpy.array([255,255,255,255],numpy.int32)
      rgba_3fold = numpy.array([0,255,255,255],numpy.int32)
      rgba_4fold = numpy.array([255,255,0,255],numpy.int32)
      rgba_6fold = numpy.array([255,0,0,255],numpy.int32)
      rgba_corner = numpy.array([255,0,255,255],numpy.int32)
      positions = []
      colors = []
      d = Drawing('points')
      d.vertices, d.normals, d.triangles = sphere
      
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
      self.session.models.add([m])
  
'''    
