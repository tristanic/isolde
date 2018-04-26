# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



# Adds CMAP corrections for improved AMBER performance in implicit solvent
# Ref: http://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b00662
import os
import numpy
import math
from simtk import openmm

class CMAPLoader:
    _map_index = {
        'GLY': 0,
        'PRO': 1,
        'ALA': 2,
        'CYS': 3,
        'CYX': 3,
        'ASP': 3,
        'ASH': 3,
        'GLU': 3,
        'GLH': 3,
        'PHE': 3,
        'HIS': 3,
        'HIE': 3,
        'HID': 3,
        'HIP': 3,
        'ILE': 3,
        'LYS': 3,
        'LYN': 3,
        'MET': 3,
        'ASN': 3,
        'GLN': 3,
        'SER': 3,
        'THR': 3,
        'VAL': 3,
        'TRP': 3,
        'TYR': 3,
        'LEU': 3,
        'ARG': 3
    }
    _map_names = ('gly', 'pro', 'ala', 'gen')

    def __init__(self, alpha_bias = 1.0, beta_bias = 1.0):
        '''
        Loads and prepares the CMAP information necessary for correction
        of backbone behaviour of AMBER protein atoms in implicit solvent.
        Adapted from MELD (http://github.com/mcconnellab/meld)

        Args:
            alpha_bias: strength of the alpha correction, default = 1.0
            beta_bias: strength of the beta correction, default = 1.0
        '''
        self._alpha_bias = alpha_bias
        self._beta_bias = beta_bias
        self._maps = []
        self._load_maps()

    def prepare_cmap_force(self, system):
        '''
        Prepare and return the CMAP force
        '''
        cmap_force = openmm.CMAPTorsionForce()
        for m in self._maps:
            cmap_force.addMap(m.shape[0], m.flatten())

        return cmap_force

    @property
    def cmaps(self):
        return self._maps

    def __getitem__(self, resname):
        return self._map_index[resname]

    def map_index(self, resname):
        return self._map_index[resname]

    def _load_map(self, stem):
        '''
        Load the maps from file. TODO: reformat these into something
        more human-readable that avoids the need for the below numpy
        rearrangements.
        '''
        basedir = os.path.join(os.path.dirname(__file__), 'amberff', 'amap')
        alpha = (numpy.loadtxt(os.path.join(basedir, '{}_alpha.txt'.format(stem))) *
                self._alpha_bias)
        beta = (numpy.loadtxt(os.path.join(basedir, '{}_beta.txt'.format(stem))) *
                self._beta_bias)
        total = alpha+beta
        assert total.shape[0] == total.shape[1]
        n = int(math.ceil(total.shape[0] / 2.0))
        total = numpy.roll(total, -n, axis=0)
        total = numpy.roll(total, -n, axis=1)
        total = numpy.flipud(total)
        return total

    def _load_maps(self):
        '''Load the maps from disk and apply the alpha and beta biases.'''
        for name in self._map_names:
            self._maps.append(self._load_map(name))
