
_cif_sf_table_names = (
    'audit',                        # General information about file
    'cell',                         # Crystallographic cell information
    'diffrn',                       # General diffraction information
    'diffrn_radiation_wavelength',  # X-ray wavelength
    'entry',                        # PDB ID
    'exptl_crystal',                # Unique crystal identifier for multi-dataset cifs
    'reflns_scale',                 #
    'symmetry',                     # Space group identification
    'reflns',                       # General information about reflections
    'refln',                        # Actual reflection data
    )

_metadata_tables = ('audit', 'diffrn', 'diffrn_radiation_wavelength', 'entry',
    'exptl_crystal', 'reflns')

_required_columns = {
    'cell':       ('length_a', 'length_b', 'length_c',
                    'angle_alpha', 'angle_beta', 'angle_gamma'),
    'refln':      ('index_h', 'index_k', 'index_l')
}

_space_group_identifiers = (
    'space_group_name_H-M',
    'Int_Tables_number'
)


from .. import (
    HKL_data_ABCD,
    HKL_data_E_sigE,
    HKL_data_F_phi,
    HKL_data_F_sigF,
    HKL_data_F_sigF_ano,
    HKL_data_I_sigI,
    HKL_data_I_sigI_ano,
    HKL_data_Phi_fom,
    HKL_data_Flag,
    HKL_data_Flag_bool
)

_data_columns_to_data_types = {
    ('F_meas_au', 'F_meas_sigma_au'):                   (float, HKL_data_F_sigF),
    ('F_meas', 'F_meas_sigma'):                         (float, HKL_data_F_sigF),
    ('pdbx_F_plus', 'pdbx_F_minus',
     'pdbx_F_plus_sigma', 'pdbx_F_minus_sigma'):        (float, HKL_data_F_sigF_ano),

    ('pdbx_anom_difference',
     'pdbx_anom_difference_sigma'):                     (float, HKL_data_F_sigF),
    ('pdbx_anomalous_diff',
     'pdbx_anomalous_diff_sigma'):                      (float, HKL_data_F_sigF),

    ('intensity_meas', 'intensity_sigma'):              (float, HKL_data_I_sigI),

    ('pdbx_I_plus', 'pdbx_I_minus',
     'pdbx_I_plus_sigma', 'pdbx_I_minus_sigma'):        (float, HKL_data_I_sigI_ano),

    ('pdbx_HL_A_iso', 'pdbx_HL_B_iso',
     'pdbx_HL_C_iso', 'pdbx_HL_D_iso'):                 (float, HKL_data_ABCD),

    ('phase_calc', 'fom'):                              (float, HKL_data_Phi_fom),

    ('pdbx_r_free_flag'):                               (int, HKL_data_Flag),
    ########
    # CALCULATED STRUCTURE FACTORS
    ########

    ('F_calc', 'phase_calc'):                           (float, HKL_data_F_phi), # Fcalc, phiFcalc
    ('F_calc_au', 'phase_calc'):                        (float, HKL_data_F_phi), # Fcalc, phiFcalc
    ('pdbx_DELFWT', 'pdbx_DELPHWT'):                    (float, HKL_data_F_phi), # mFo-DFc
    ('pdbx_FWT', 'pdbx_PHWT'):                          (float, HKL_data_F_phi), # 2mFo-DFc
    ('pdbx_F_calc_part_solvent',
     'pdbx_phase_calc_part_solvent'):                   (float, HKL_data_F_phi), # Solvent contribution to calculated structure factors
    ('pdbx_F_calc_with_solvent',
     'pdbx_phase_calc_with_solvent'):                   (float, HKL_data_F_phi), # Calculated structure factors including bulk solvent

}

_anomalous_data_columns = (
    ('pdbx_F_plus', 'pdbx_F_minus', 'pdbx_F_plus_sigma', 'pdbx_F_minus_sigma'),
    ('pdbx_I_plus', 'pdbx_I_minus', 'pdbx_I_plus_sigma', 'pdbx_I_minus_sigma')
)


def load_cif_sf(filename):
    '''
    Load a set of structure factors from a .cif file.
    '''
    from chimerax.atomic.mmcif import get_cif_tables
    table_list = get_cif_tables(filename, _cif_sf_table_names)

    tables = dict(zip(_cif_sf_table_names, table_list))
    metadata = {l: tables[l] for l in _metadata_tables}
    cell_info = tables['cell']
    try:
        cell_dim = [float(d) for d in cell_info.fields(_required_columns['cell'])[0]]
    except:
        raise TypeError('Could not read cell information from file!')
    symm = tables['symmetry']
    for id in _space_group_identifiers:
        if symm.has_field(id):
            spgr_descriptor = symm.fields((id,))[0][0]
            break
    else:
        raise TypeError('Could not read spacegroup information from file!')

    refln_table = tables['refln']
    hkls = _get_miller_indices(refln_table)

    from ..clipper_python import (
        Cell_descr, Cell, Spgr_descr,
        Spacegroup, Resolution, Grid_sampling
    )
    cell = Cell(Cell_descr(*cell_dim))
    spacegroup = Spacegroup(Spgr_descr(spgr_descriptor))

    refln_info = tables['reflns']
    # Resolution information is explicitly given in the 'reflns' table in a
    # surprisingly small minority of structure factor cif files, but it's nice when
    # we can get it.
    if refln_info.has_field('d_resolution_high'):
        res = float(refln_info.fields(('d_resolution_high',))[0][0])
    else:
        # If this turns out to be too slow, it could fairly easily be pushed
        # back to C++
        from ..clipper_python import HKL
        invresolsq_limit = max(HKL(hkl).invresolsq(cell) for hkl in hkls)
        from math import sqrt
        res = 1/sqrt(invresolsq_limit)
    print('Resolution: {}'.format(res))
    resolution = Resolution(res)

    from .. import HKL_info
    hkl_info = HKL_info(spacegroup, cell, resolution, True)

    # OK, now we have all our vital information. Time to find all the data
    data = _parse_status(refln_table, hkl_info, cell, hkls)
    for column_names, type_spec in _data_columns_to_data_types.items():
        result = _cif_columns_to_clipper(refln_table, hkl_info, cell, hkls, column_names, type_spec)
        if result is not None:
            data[column_names] = result

    return cell, spacegroup, resolution, hkl_info, data


def _get_miller_indices(table):
    import numpy
    headers = ('index_h', 'index_k', 'index_l')
    hkls = numpy.array([[int(i) for i in row] for row in table.fields(headers)],
                        dtype=numpy.int32)
    return hkls

def _cif_columns_to_clipper(table, hkl_info, cell, hkls, field_names, type_spec):
    import numpy
    for field in field_names:
        if not table.has_field(field):
            return None
    scalar_type, array_type = type_spec
    if scalar_type == float:
        data = numpy.array(
            [[float(d) if d!='?' else numpy.nan for d in row] for row in table.fields(field_names)],
            numpy.double
        )
    elif scalar_type == int:
        data = numpy.array(
            [[int(d) for d in row] for row in table.fields(field_names)],
            numpy.double
        )
    if field_names in _anomalous_data_columns:
        # Clipper's anomalous data types have a fifth covariance column, but
        # this is not provided in the .cif file. Just set it to zero.
        padded_data = numpy.empty((len(data), 5), numpy.double)
        padded_data[:,:4] = data
        padded_data[:,4] = 0
        data = padded_data
    clipper_data = array_type(hkl_info, cell)
    try:
        clipper_data.set_data(hkls, data)
    except RuntimeError as e:
        err_string = " Field names: {}; HKL array shape: {}; Data array shape: {}".format(
            field_names, hkls.shape, data.shape
        )
        raise RuntimeError(str(e)+err_string)
    return clipper_data

_refln_status_flags = {
    'systematically_absent':    '-',    # Doesn't exist in the given space group
    'unobserved':               '<',    # Too low to measure, but not systematically absent
    'free':                     'f',    # Include in the R-free set
    'cut_high':                 'h',    # Measured, but above the decided high resolution cutoff
    'cut_low':                  'l',    # Measured, but below the decided low resolution cutoff
    'unreliable':               'x',    # Unreliable measurement. Do not use.
    'default':                  'o'     # Normal
}



def _parse_status(table, hkl_info, cell, hkls):
    import numpy
    status = table.fields(('status',))
    n = len(hkls)
    # Due to a quirk of Clipper's template instantiation, all HKL_data types
    # import/export as double regardless of their actual data type.
    free = numpy.zeros((n,1), numpy.double)
    unreliable = numpy.zeros((n,1), numpy.double)
    systematic_absences = numpy.zeros((n,1), numpy.double)
    for i, s in enumerate(status):
        s = s[0]
        if s == 'o':
            continue
        elif s == '-':
            systematic_absences[i] = 1
        elif s == 'f':
            free[i] = 1
        elif s == 'x':
            unreliable[i] = 1
    c_free = HKL_data_Flag(hkl_info, cell)
    c_unreliable = HKL_data_Flag_bool(hkl_info, cell)
    c_systematic_absences = HKL_data_Flag_bool(hkl_info, cell)
    c_free.set_data(hkls, free)
    c_unreliable.set_data(hkls, unreliable)
    c_systematic_absences.set_data(hkls, systematic_absences)
    return {"Free set":     c_free,
            "Unreliable":   c_unreliable,
            "Systematic absences":  c_systematic_absences
            }
