#
#     Copyright (C) 2015-16 CCP-EM
#
#     This code is distributed under the terms and conditions of the
#     CCP-EM Program Suite Licence Agreement as a CCP-EM Application.
#     A copy of the CCP-EM licence can be obtained by writing to the
#     CCP-EM Secretary, RAL Laboratory, Harwell, OX11 0FA, UK.
#

import clipper
import numpy
import numpy as np
import os
import matplotlib.pyplot as plt


class ClipperMTZ(object):
    def __init__(self,
                 mtz_in_path=None):
        self.mtz = clipper.CCP4MTZfile()
        self.hkl_info = clipper.HKL_info()
        self.mtz_in_path = mtz_in_path
        if self.mtz_in_path is not None:
            assert os.path.exists(self.mtz_in_path)
        self.spacegroup = None
        self.cell = None
        self.column_data = {}

    def import_column_data(self, column_label, get_resolution=True):
        if self.mtz_in_path is not None:
            # Convert to ascii from if unicode
            if isinstance(self.mtz_in_path, unicode):
                self.mtz_in_path = self.mtz_in_path.encode('utf8')
            self.mtz.open_read(self.mtz_in_path)
            self.mtz.import_hkl_info(self.hkl_info)
            self.spacegroup = self.hkl_info.spacegroup()
            self.cell = self.hkl_info.cell()
            f_phi = clipper.HKL_data_F_phi_float(self.hkl_info)
            self.mtz.import_hkl_data(f_phi, '/*/*/' + column_label)
            self.mtz.close_read()
            # Convert to numpy
            f_phi_np = numpy.zeros((f_phi.data_size() * len(f_phi)),
                                      numpy.float)
            f_phi.getDataNumpy(f_phi_np)
            # Reshape and transpose
            f_phi_np = numpy.reshape(f_phi_np, (-1, 2))
            f_phi_np = numpy.transpose(f_phi_np)
            # Convert to rec array to store col names
            names = [n for n in f_phi.data_names().split()]
            f_phi_np = np.core.records.fromarrays(
                f_phi_np,
                names=names,
                formats='float64, float64')
            # Append to dictionary
            self.column_data[column_label] = f_phi_np
            # Get resolution column
            if get_resolution:
                res_np = numpy.zeros(f_phi_np.shape[0])
                for n in xrange(f_phi_np.shape[0]):
                    r = self.hkl_info.invresolsq(n)
                    res_np[n] = r
                self.column_data['resolution_1/Ang^2'] = res_np

