# @Author: Tristan Croll
# @Date:   01-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
def test_clipper_sym(session, radius=12):
    import os
    from chimerax.core.commands import open
    libdir = os.path.abspath(os.path.dirname(__file__))
    m = open.open(session, os.path.join(libdir,'3io0.pdb'))[0]
    m.atoms.displays = True
    from chimerax.std_commands import cofr
    cofr.cofr(session, method='center of view', show_pivot=True)

    from chimerax.clipper import symmetry
    sym_handler = symmetry.get_symmetry_handler(m, create=True)
    sym_handler.spotlight_radius = radius
    sym_handler.map_mgr.add_xmapset_from_mtz(os.path.join(libdir,'3io0_combined.mtz'), oversampling_rate=1.5)
    sym_handler.map_mgr.nxmapset.add_nxmap_handler_from_file(os.path.join(libdir, '3io0_real_space.ccp4'))
    return sym_handler
