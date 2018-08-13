import os

_basedir = os.path.dirname(os.path.abspath(__file__))

def test_xtal_mgr(session):
    from chimerax.core.commands import open
    m = open.open(session, os.path.join(_basedir, '3io0.pdb'))[0]
    from ..clipper_mtz import ReflectionDataContainer
    from .. import atom_list_from_sel
    from ..clipper_python.ext import Xtal_mgr

    rdc = ReflectionDataContainer(session, os.path.join(_basedir, '3io0-sf.mtz'))
    xm = Xtal_mgr(rdc.hklinfo, rdc.free_flags.data, rdc.grid_sampling,
        rdc.experimental_data.datasets['FOBS, SIGFOBS'].data)
    xm.generate_fcalc(atom_list_from_sel(m.atoms))
    print("Rfree: {}".format(xm.rfree))
    print("Rwork: {}".format(xm.rwork))
    return xm, rdc
