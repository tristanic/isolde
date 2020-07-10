# @Author: Tristan Croll <tic20>
# @Date:   20-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 23-May-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from enum import IntEnum
class map_styles(IntEnum):
    mesh_square             = 0
    mesh_triangle           = 1
    solid_t20               = 2
    solid_t40               = 3
    solid_t60               = 4
    solid_t80               = 5
    solid_opaque            = 6

human_readable_map_display_styles = {
    map_styles.mesh_square: "Mesh (squares)",
    map_styles.mesh_triangle: "Mesh (triangles)",
    map_styles.solid_t20: "Solid (20% opacity)",
    map_styles.solid_t40: "Solid (40% opacity)",
    map_styles.solid_t60: "Solid (60% opacity)",
    map_styles.solid_t80: "Solid (80% opacity)",
    map_styles.solid_opaque: "Solid (opaque)"
    }

# array of settings to apply to the map depending on the chosen
# representation.
map_style_settings = {
    map_styles.mesh_square: {'style': 'mesh', 'square_mesh': True, 'transparency':0, 'flip_normals': True},
    map_styles.mesh_triangle: {'style': 'mesh', 'square_mesh': False, 'transparency':0, 'flip_normals': True},
    map_styles.solid_t20: {'style': 'surface', 'transparency': 0.8},
    map_styles.solid_t40: {'style': 'surface', 'transparency': 0.6},
    map_styles.solid_t60: {'style': 'surface', 'transparency': 0.4},
    map_styles.solid_t80: {'style': 'surface', 'transparency': 0.2},
    map_styles.solid_opaque: {'style': 'surface', 'transparency': 0.0}
    }

def default_atom_visualisation(model):
    session = model.session
    from chimerax.core.objects import Objects
    from chimerax.std_commands import color
    from chimerax.atomic import Atom
    atoms = model.atoms
    ao = Objects(atoms=atoms)
    color.color(session, ao, color='bychain')
    color.color(session, ao, color='byhetero')
    model.bonds.radii = 0.2
    model.bonds.displays=True
    atoms.draw_modes = Atom.STICK_STYLE
    atoms.displays = True
    from chimerax.clipper.util import nonpolar_hydrogens
    atoms[nonpolar_hydrogens(atoms)].displays=False
    model.residues.ribbon_displays = True
