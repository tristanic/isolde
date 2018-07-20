# @Author: Tristan Croll
# @Date:   16-Oct-2017
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



from . import clipper_python as clipper
import numpy

# Override logging of Clipper messages to send them to the ChimeraX log
# rather than just printing them to the console.

# _clipper_messages = clipper._clipper_messages
#
# def _log_clipper(func):
#   def func_wrapper(*args, **kwargs):
#     _clipper_messages.clear()
#     func(*args, **kwargs)
#     message_string = _clipper_messages.read_and_clear()
#     if message_string:
#       session.logger.info("CLIPPER WARNING:")
#       session.logger.info(message_string)
#   return func_wrapper

# clipper._clipper.log_clipper = _log_clipper

def add_crystal(session, name):
  from .crystal import Xtal_Project
  return Xtal_Project(session, name)





def apply_b_factors_to_hydrogens(atom_list):
    '''
    When hydrogens are added to a pre-existing structure, their isotropic
    and/or anisotropic B-factors are not necessarily set. This routine simply
    sets them to match their adjacent heavy atoms. The input atom_list can
    be either the hydrogens alone, or the whole structure.
    '''
    hydrogens = atom_list.filter(atom_list.element_names == 'H')
    # Get heavy atoms by taking advantage of the fact that hydrogens are
    # only ever covalently bonded to a single other atom.
    for h in hydrogens:
        b = h.bonds[0]
        for atom in b.atoms:
            if atom.element_name != 'H':
                h.bfactor = atom.bfactor
                h.aniso_u6 = atom.aniso_u6
                break


def atom_list_from_sel(atom_list):
    '''
    Takes a ChimeraX Atoms object, and creates a Clipper Atoms_list object
    from the relevant atom properties.
    '''
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancies
    u_iso = atom_list.bfactors
    # Have to be careful here. Atoms.aniso_u6 returns None if any atom in
    # the array has no aniso_u6 entry.
    u_aniso = atom_list.aniso_u6
    if u_aniso is None:
        u_aniso = numpy.ones([n,6],numpy.float32)*numpy.nan
        u_aniso[atom_list.has_aniso_u] = atom_list.filter(atom_list.has_aniso_u).aniso_u6
        # FIXME Once new ChimeraX build arrives with Atoms.has_aniso_u entry
        #for i, a in enumerate(atom_list):
            #if a.aniso_u6 is not None:
                #u_aniso[i] = a.aniso_u6
    clipper_atom_list = clipper.Atom_list(elements, coords, occupancies, u_iso, u_aniso)
    return clipper_atom_list

def make_unit_cell(model, xmap, draw = True):
    from chimerax.core.geometry import Place, Places
    cell = xmap.cell()
    atoms = model.atoms
    clipper_atoms = atom_list_from_sel(atoms)
    coord = model.bounds().center()
    coord_orth = clipper.Coord_orth(coord)
    coord_frac = coord_orth.coord_frac(cell)
    sg = xmap.spacegroup()
    unit_cell_frac_symops = xmap.unit_cell_symops(coord_frac, clipper_atoms)
    if draw:
        uc_places = []
        for op in unit_cell_frac_symops.symops():
            op_orth = op.rtop_orth(cell)
            uc_places.append(Place(matrix=op_orth.matrix()[0:3,:]))
        ucp = Places(uc_places)
        model.positions = ucp
    return unit_cell_frac_symops

def test_pack_box(model, xmap, size = 50):
    from chimerax.core.geometry import Place, Places
    uc = make_unit_cell(model, xmap, draw = False)
    coord = model.bounds().center()
    box_size = (numpy.ones(3, numpy.int32)*size)
    bo = xmap.all_symops_in_box(coord-box_size/2, box_size, uc)
    p = []
    for b in bo:
        p.append(Place(matrix=b.rtop_orth(xmap.cell()).matrix()[0:3,:]))
    P = Places(p)
    model.positions = P

    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('box corners')
    m = Model('box', session)
    d.set_geometry(*sphere_geometry(80))
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    box_min = clipper.Coord_orth(coord-box_size/2).coord_frac(xmap.cell()).coord_grid(xmap.grid_sampling())
    box_max = box_min + clipper.Coord_grid(*box_size.tolist())
    minmax = [box_min, box_max]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                thiscg = clipper.Coord_grid([minmax[i].u(), minmax[j].v(), minmax[k].w()])
                dp.append(Place(origin=thiscg.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()).xyz))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])

    return bo

def pack_box(model, xmap, box_origin_xyz, size = 100):
    box_size_xyz = numpy.ones(3)*size
    box_origin_xyz = numpy.array(box_origin_xyz)
    model_bounds = model.bounds().box_corners()
    bo = xmap.pack_xyz_box(model_bounds, box_origin_xyz, box_size_xyz);
    from chimerax.core.geometry import Place, Places

    p = []
    for b in bo:
        p.append(Place(matrix=b.rtop_orth(xmap.cell()).matrix()[0:3,:]))
    P = Places(p)
    model.positions = P

    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('box corners')
    m = Model('box', session)
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    d.set_geometry(*sphere_geometry(80))
    minmax = [box_origin_xyz, box_origin_xyz+box_size_xyz]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    return bo

def draw_box(min_corner, max_corner, name='box'):
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('corners')
    m = Model(name, session)
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    d.set_geometry(*sphere_geometry(80))
    minmax = [min_corner, max_corner]
    dp = []
    base_color = numpy.array([255,255,255,255])
    color_increment = numpy.array([0,-32,-32,0])
    colors = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
                colors.append(base_color-color_increment*(i+j+k));
    d.positions = Places(dp)
    d.colors = colors
    m.add_drawing(d)
    session.models.add([m])
    return d




def draw_asu(xmap):
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('asu corners')
    m = Model('asu box', session)
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    d.set_geometry(*sphere_geometry(80))
    asu = xmap.grid_asu()
    grid = xmap.grid_sampling()
    cell = xmap.cell()
    minmax = [asu.min().coord_frac(grid).coord_orth(cell).xyz, asu.max().coord_frac(grid).coord_orth(cell).xyz]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    return d

def box_corners(origin_xyz, size_xyz):
    ret = []
    minmax = [origin_xyz, size_xyz]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ret.append([minmax[i][0],minmax[j][1], minmax[k][2]])
    return ret
