# @Author: Tristan Croll <tic20>
# @Date:   10-Jul-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 20-Jul-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



# Managers
from ..molobject import ChiralMgr, ProperDihedralMgr, RamaMgr, RotaMgr

# Manager get functions
from ..session_extensions import get_proper_dihedral_mgr, get_ramachandran_mgr, \
                                 get_rotamer_mgr

# Singular objects
from ..molobject import ChiralCenter, ProperDihedral, Rama, Rotamer

# Collections
from ..molarray import ChiralCenters, ProperDihedrals, Ramas, Rotamers

global _templates_loaded
_templates_loaded = False

def load_default_mmcif_templates():
    global _templates_loaded
    if not _templates_loaded:
        import os
        from chimerax import mmcif
        base_dir = os.path.dirname(os.path.abspath(__file__))
        mmcif.load_mmCIF_templates(os.path.join(base_dir, '..', 'dictionaries', 'core_components.cif'))
        _templates_loaded=True
