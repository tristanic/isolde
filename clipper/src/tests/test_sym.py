import numpy
def test_clipper_sym(session, radius=25):
    import os
    from chimerax.core.commands import open
    libdir = os.path.abspath(os.path.dirname(__file__))
    m = open.open(session, os.path.join(libdir,'5f4y.pdb'))[0]
    from chimerax.clipper import CrystalStructure
    cs = CrystalStructure(session, m, os.path.join(libdir, '5f4y_map_coeffs.mtz'))
    return sym_around_current_cofr(session, cs, radius)

def sym_around_current_cofr(session, cs, radius=25):
    m = cs.master_model
    from chimerax.clipper.crystal import _find_box_corner
    center = session.view.center_of_rotation
    box_corner_grid, box_corner_xyz = _find_box_corner(center, cs.cell, cs.grid, radius)
    dim = numpy.array([radius*2,radius*2,radius*2], numpy.int32)
    symops = symops = cs.unit_cell.all_symops_in_box(box_corner_xyz, dim, True)
    tfs = symops.all_matrices_orth(cs.cell)
    from chimerax.clipper import symmetry
    m.atoms.displays = True
    results = symmetry.sym_transforms_in_sphere(m.atoms, tfs, session.view.center_of_rotation, radius)
    return m, results
