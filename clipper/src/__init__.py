# @Author: Tristan Croll
# @Date:   09-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



from .main import *
# General objects
from .clipper_python import (
    Cell,
    Cell_descr,
    CCP4MTZfile,
    CIFfile,
    Coord_frac,
    Coord_grid,
    Coord_orth,
    Grid,
    Grid_range,
    Grid_sampling,
    HKL_info,
    Map_stats,
    Resolution,
    RTop_frac,
    RTop_orth,
    Spacegroup,
    Spgr_descr,
    Symop,
    # Symops,
    Unit_Cell,
    Util,
    Xmap_float as Xmap,
    )


# Singular forms of HKL data
from .clipper_python import Flag, Flag_bool
from .clipper_python.data32 import (
    ABCD_float as ABCD,
    E_sigE_float as E_sigE,
    F_phi_float as F_phi,
    F_sigF_float as F_sigF,
    F_sigF_ano_float as F_sigF_ano,
    I_sigI_float as I_sigI,
    I_sigI_ano_float as I_sigI_ano,
    Phi_fom_float as Phi_fom,
    )

# Array forms of HKL data
from .clipper_python import HKL_data_Flag, HKL_data_Flag_bool
from .clipper_python.data32 import (
    HKL_data_ABCD_float as HKL_data_ABCD,
    HKL_data_E_sigE_float as HKL_data_E_sigE,
    HKL_data_F_phi_float as HKL_data_F_phi,
    HKL_data_F_sigF_float as HKL_data_F_sigF,
    HKL_data_F_sigF_ano_float as HKL_data_F_sigF_ano,
    HKL_data_I_sigI_float as HKL_data_I_sigI,
    HKL_data_I_sigI_ano_float as HKL_data_I_sigI_ano,
    HKL_data_Phi_fom_float as HKL_data_Phi_fom,
    )

from .clipper_mtz import ReflectionDataContainer

from chimerax.core.toolshed import BundleAPI
class _ClipperBundle(BundleAPI):
    from chimerax.core.commands import FloatArg
    from chimerax.atomic import StructureArg
    @staticmethod
    def initialize(session, bundle_info):
        from chimerax.clipper import cmd
        # cmd.register_mtz_file_format(session)

    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        from chimerax.clipper import cmd
        if command_name == 'clipper':
            cmd.register_clipper_cmd(logger)

    @staticmethod
    def open_file(session, path, format_name, structure_model=None,
            over_sampling=1.5):
        if structure_model is None:
            from chimerax.core.errors import UserError
            raise UserError('Must specify a structure model to associate with crystallographic data')
        if format_name == 'mtz':
            from .cmd import open_mtz
            return open_mtz(session, path, structure_model=structure_model,
                over_sampling=over_sampling)

bundle_api = _ClipperBundle()
