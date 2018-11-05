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
    # CIFfile,
    Coord_grid,
    Grid,
    Grid_range,
    Grid_sampling,
    HKL_info,
    Resolution,
    RTop_frac,
    RTop_orth,
    Spacegroup,
    Spgr_descr,
    Symop,
    # Symops,
    Unit_Cell,
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
    Phi_fom_float as Phi_fom,
    )

# Array forms of HKL data
from .clipper_python import HKL_data_Flag, HKL_data_Flag_bool
from .clipper_python.data32 import (
    HKL_data_ABCD_float as HKL_data_F_ABCD,
    HKL_data_E_sigE_float as HKL_data_E_sigE,
    HKL_data_F_phi_float as HKL_data_F_phi,
    HKL_data_F_sigF_float as HKL_data_F_sigF,
    HKL_data_F_sigF_ano_float as HKL_data_F_sigF_ano,
    HKL_data_I_sigI_float as HKL_data_I_sigI,
    HKL_data_Phi_fom_float as HKL_data_Phi_fom,
    )

from .clipper_mtz import ReflectionDataContainer
