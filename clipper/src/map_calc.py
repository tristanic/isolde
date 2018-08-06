
def generate_fcalc(container, atoms, fsigfdata, target=None):
    '''
    Generate a crystallographic map from the given atoms and observed
    structure factors.

    Args:
        * atoms:
            - A :class:`chimerax.Atoms` object containing all atoms to include
              in the calculation.
        * fsigfdata:
            - A :class:`ReflectionData_Exp` object holding the observed
              structure factors
        * target:
            - A :class:`ReflectionData_Calc` object, or None. If not None, any
              data in the existing
    '''
    session = fsigfdata.session
    fsigf = fsigfdata.data
    hkls = container.hklinfo
    from . import atom_list_from_sel
    clipper_atoms = atom_list_from_sel(atoms)

    if target is None:
        from .clipper_mtz import ReflectionData_Calc
        from .clipper_python.data64 import HKL_data_F_phi_double
        target = ReflectionData_Calc('Fcalc', session,
            HKL_data_F_phi_double(hkls), is_difference_map=False)
    fcalc = target.data
    from .clipper_python import SFcalc_obs_bulk_double
    SFcalc_obs_bulk_double(fcalc, fsigf, clipper_atoms)
    container.calculated_data.add([target])

def generate_map_coeffs(container, fsigf, fcalc, free_r_flags):
    from .clipper_python.data64 import HKL_data_F_phi_double, HKL_data_Phi_fom_double
    session = container.session
    hkls = container.hklinfo
    best_coeffs = HKL_data_F_phi_double(hkls)
    diff_coeffs = HKL_data_F_phi_double(hkls)
    phiw = HKL_data_Phi_fom_double(hkls)

    from .clipper_python import SFweight_spline_double
    SFweight_spline_double(best_coeffs, diff_coeffs, phiw, fsigf.data, fcalc.data, free_r_flags.data)

    ret = []
    from .clipper_mtz import ReflectionData_Calc
    ret.append(ReflectionData_Calc('2FOFCWT, PH2FOFCWT', session, best_coeffs, is_difference_map=False))
    ret.append(ReflectionData_Calc('FOFCWT, PHFOFCWT', session, diff_coeffs, is_difference_map=True))
    return ret
    # for b in sharpening_factors:
    #     if b == 0:
    #         ret.append(best_coeffs)
    #     else:












    # if sharpening is None:
    #     if target is not None and hasattr(target, 'sharpening'):
    #         sharpening = target.sharpening
    #     else:
    #         sharpening = 0
    #
    # if sharpening < 0:
    #     name_ext = " sharp {:0.1f}".format(sharpening)
    # elif sharpening > 0:
    #     name_ext = " smooth {:0.1f}".format(sharpening)
    # else:
    #     name_ext = ""
    # name_string = type + name_ext
