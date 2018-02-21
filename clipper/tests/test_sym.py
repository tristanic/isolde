import numpy
import os
dpath = os.path.dirname(os.path.abspath(__file__))
from chimerax.core.commands import open
m = open.open(session, os.path.join(dpath, '5f4y.pdb')[0]
from chimerax.clipper import CrystalStructure
cs = CrystalStructure(session, m, os.path.join('5f4y_map_coeffs.mtz')
from chimerax.clipper.crystal import _find_box_corner

box_corner_grid, box_corner_xyz = _find_box_corner(cs._sym_box_center, cs.cell, cs.grid, 25)
dim = numpy.array([50,50,50], numpy.int32)
symops = symops = cs.unit_cell.all_symops_in_box(box_corner_xyz, dim, True)
tfs = symops.all_matrices_orth(cs.cell)
from chimerax.clipper import geometry
m.atoms.displays = True
results = geometry.sym_transforms_in_sphere(m.atoms, tfs, session.view.center_of_rotation, 50)
