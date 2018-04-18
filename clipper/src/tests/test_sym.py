# @Author: Tristan Croll
# @Date:   01-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
def test_clipper_sym(session, radius=15):
    import os
    from chimerax.core.commands import open
    libdir = os.path.abspath(os.path.dirname(__file__))
    m = open.open(session, os.path.join(libdir,'5f4y.pdb'))[0]
    m.atoms.displays = True
    from chimerax.core.commands import cofr
    cofr.cofr(session, method='center of view', show_pivot=True)

    from chimerax.clipper import symmetry
    sym_handler = symmetry.XtalSymmetryHandler(m, os.path.join(libdir,'5f4y_map_coeffs.mtz'),
        atomic_symmetry_radius = radius)
    return sym_handler
