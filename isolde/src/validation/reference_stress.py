# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Baseline-subtraction layer for the local virial stress signal
(local_virial.LocalVirialCalculator). Even a perfectly-modelled, well-minimised
structure carries a systematic per-(residue type, atom name) "frustration" floor
that is a property of the force field, not the model -- e.g. persistent stress at
arginine CZ, where AMBER's three independently-fitted guanidinium angle equilibria
can't all be satisfied at once, at any geometry. ReferenceStressLibrary subtracts
that floor (in tensor space, before scalar invariants are taken, so an anomalous
stress *direction* still shows up even if its magnitude looks ordinary) so the
residual reflects real modelling errors instead.

The reference tensors are built offline (see tools/build_reference_stress_library.py)
from a curated set of high-resolution, well-refined structures run through ISOLDE's
own load+minimise pipeline (rather than isolated capped monomers -- ISOLDE has no
ligand-free monomer-building utility, but reuses its existing SimConstruct/
SimHandler/minimize machinery unchanged for this), and cached to disk as the
element-wise median tensor per (residue name, atom name).
'''

import numpy as np

from .local_virial import LocalVirialCalculator

_DEFAULT_DATA_FILE = 'data/reference_virial.npz'


def _key(resname, atom_name):
    return f'{resname}:{atom_name}'


class ReferenceStressLibrary:
    '''
    Loads a cached per-(residue name, atom name) baseline virial tensor and uses
    it to compute excess (baseline-subtracted) stress invariants for a live
    per-atom virial.

    Args:
        path: path to a .npz cache written by `save()` (or built via
            `from_samples()`). Defaults to the bundled
            `validation/data/reference_virial.npz`.
    '''

    # Floor on the per-type spread scale, in the virial's kJ/mol units, so a type
    # with a near-zero (or single-sample) spread can't blow the z-score up.
    _MIN_SPREAD = 1.0

    def __init__(self, path=None):
        if path is None:
            import os
            path = os.path.join(os.path.dirname(__file__), _DEFAULT_DATA_FILE)
        self._baseline, self._scales = self._load(path)

    @classmethod
    def _load(cls, path):
        data = np.load(path)
        keys = data['keys']
        tensors = data['tensors']
        baseline = {key: tensors[i] for i, key in enumerate(keys)}
        # 'scales' is optional -- older caches predate variance normalisation.
        if 'scales' in data:
            scales = {key: float(data['scales'][i]) for i, key in enumerate(keys)}
        else:
            scales = {}
        return baseline, scales

    def save(self, path):
        keys = np.array(list(self._baseline.keys()))
        tensors = np.array([self._baseline[k] for k in keys])
        scales = np.array([self._scales.get(k, np.nan) for k in keys])
        np.savez(path, keys=keys, tensors=tensors, scales=scales)

    @classmethod
    def from_samples(cls, samples):
        '''
        Build a library directly from accumulated samples, without going
        through disk. Used by the offline builder script.

        Args:
            samples: dict mapping `"RESNAME:ATOMNAME"` -> list/array of (3,3)
                symmetrised virial tensors collected across the reference
                structure set.

        For each key we store two things:
          * baseline: the element-wise median tensor (the systematic frustration
            floor for that atom type), and
          * scale: the *typical spread* of the deviatoric stress about that
            median -- the median over samples of the von Mises invariant of
            ``sample - baseline``. This is the natural per-type yardstick for the
            excess: an atom type whose stress is intrinsically large *and*
            variable (polar/charged sidechain tips) gets a large scale, so an
            ordinary-for-its-type deviation reports a small z-score, while a
            normally-quiet atom (much of the backbone) that deviates stands out.
        '''
        self = cls.__new__(cls)
        self._baseline = {}
        self._scales = {}
        for key, tensors in samples.items():
            arr = np.asarray(tensors)
            median = np.median(arr, axis=0)
            self._baseline[key] = median
            deviations = arr - median
            vm = LocalVirialCalculator.scalar_measures(deviations)['von_mises']
            self._scales[key] = max(float(np.median(vm)), cls._MIN_SPREAD)
        return self

    def has_baseline(self, resname, atom_name):
        return _key(resname, atom_name) in self._baseline

    def compute_excess(self, virial, atoms):
        '''
        Args:
            virial: (N,3,3) raw per-atom virial, e.g. from
                `LocalVirialCalculator.compute()['virial']`.
            atoms: a ChimeraX Atoms collection in the same order as `virial`'s
                first axis (e.g. `sim_construct.all_atoms`).

        Returns:
            dict of (N,) [or (N,3) for 'principal'] arrays. Keys:
              hydrostatic, principal, von_mises, max_shear -- invariants of the
                residual (live_symmetrised - baseline) tensor (as
                `LocalVirialCalculator.scalar_measures()`);
              scale -- the per-type spread yardstick for each atom (NaN if the
                cache predates variance normalisation);
              von_mises_z -- von_mises / scale, the variance-normalised excess
                (how many typical-for-this-atom-type deviations the excess is).
                This is the discriminating metric: it de-emphasises atoms that are
                intrinsically high-variance (polar/charged tips) and highlights
                atoms far outside their own type's normal spread.
            Atoms with no baseline entry (non-standard residues, novel ligands,
            chain termini) get NaN; von_mises_z is also NaN where no scale exists.
        '''
        n = len(atoms)
        if virial.shape[0] != n:
            raise ValueError(
                f'virial has {virial.shape[0]} atoms but {n} atoms were given -- '
                'these must be the same collection the OpenMM System was built from.')
        live_sym = 0.5 * (virial + np.transpose(virial, (0, 2, 1)))
        residual = np.full((n, 3, 3), np.nan)
        scale = np.full(n, np.nan)
        has_baseline = np.zeros(n, dtype=bool)
        for i, atom in enumerate(atoms):
            key = _key(atom.residue.name, atom.name)
            baseline = self._baseline.get(key)
            if baseline is None:
                continue
            residual[i] = live_sym[i] - baseline
            has_baseline[i] = True
            if key in self._scales:
                scale[i] = self._scales[key]

        measures = {
            'hydrostatic': np.full(n, np.nan), 'von_mises': np.full(n, np.nan),
            'max_shear': np.full(n, np.nan), 'principal': np.full((n, 3), np.nan)}
        if has_baseline.any():
            sub = LocalVirialCalculator.scalar_measures(residual[has_baseline])
            for key, arr in measures.items():
                arr[has_baseline] = sub[key]
        measures['scale'] = scale
        measures['von_mises_z'] = measures['von_mises'] / scale
        return measures
