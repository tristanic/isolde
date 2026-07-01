# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Per-atom virial (local stress) tensor, reconstructed from an OpenMM System's
per-interaction forces before they are summed away into the net per-atom force
(which is ~0 at a minimum and therefore useless for finding strained sites).

Supported terms (mirroring OpenMM's Reference-platform math bit-for-bit):
    * HarmonicBondForce
    * HarmonicAngleForce
    * PeriodicTorsionForce
    * NonbondedForce, or ISOLDE's NonbondedSoftcoreForce/
      NonbondedSoftcoreExceptionForce (see notes on the latter below) -- either
      NoCutoff or CutoffNonPeriodic/CutoffPeriodic with reaction field and LJ
      switching.

Intentionally omitted: GB/OBC implicit solvent (GBSAForce/SoftCoreGBSAGBnForce)
and AMBER CMAP backbone corrections (AmberCMAPForce). Both are smooth,
low-spatial-contrast terms for an error-finding signal, and reimplementing their
force loops in Python is not worthwhile unless real strain analysis needs them.
Restraint and tugging forces (RESTRAINT_FORCE_GROUP) and MDFF map coupling
(MAP_FORCE_GROUP) are excluded on principle: the signal should reflect the
physical force field only, not ISOLDE's interactive/artificial terms.

ISOLDE always runs its live simulation with a soft-core nonbonded potential
(NonbondedSoftcoreForce + NonbondedSoftcoreExceptionForce) rather than a plain
NonbondedForce -- purely as a numerical-stability device so minimisation
doesn't blow up on severe clashes, not as an intended physical change to the
force field. The softcore energy expression provably reduces to plain
LJ + Coulomb as its lambda parameter -> 1, so rather than differentiating the
general softcore expression, this module reads the physical charge/sigma/epsilon
parameters straight out of the softcore force objects and evaluates the plain
classical formula -- i.e. the lambda=1 limit -- using the same per-interaction
math as for a plain NonbondedForce.

Per-atom virial convention
--------------------------
For each interaction with participating atoms {k}, forces {f_k} and positions {r_k}:
    T = sum_k (r_k - r_ref) (x) f_k          (origin-independent since sum_k f_k = 0)
    W[k] += T / N_atoms      for each k in the interaction
For a 2-body term this reduces to the standard W = 1/2 r_ij (x) f_ij on each atom.

`W` has units kJ/mol (= nm * kJ/mol/nm). The mechanical stress is sigma = -W / Omega
for a per-atom volume Omega; we return the raw virial and provide scalar invariants,
leaving volume normalisation to the caller (it only rescales the field).

Sign of the trace: Tr(W) > 0  <=>  local compression (atoms pushed apart),
                   Tr(W) < 0  <=>  local tension     (atoms pulled together).
'''

import warnings

import numpy as np

try:
    import openmm as mm
    from openmm import unit
except ImportError:  # pragma: no cover
    import simtk.openmm as mm
    from simtk import unit

from ..openmm.openmm_interface import CORE_FORCE_GROUPS
# Already a bare float in kJ/mol*nm/e**2 (ISOLDE's `Quantity * Unit` idiom strips
# units once they're dimensionally compatible) -- not an openmm.unit.Quantity.
from ..openmm.custom_forces import ONE_ON_4_PI_EPS0 as ONE_4PI_EPS0


def _strip(value, u):
    '''Return a plain float/ndarray in OpenMM internal units (nm, kJ/mol, rad).'''
    if unit.is_quantity(value):
        return value.value_in_unit(u)
    return value


class LocalVirialCalculator:
    '''
    Parse an OpenMM System once, then map positions -> per-atom virial + forces.

    Args:
        system: an openmm.System
        force_groups: iterable of OpenMM force-group indices to reconstruct.
            Default: CORE_FORCE_GROUPS ({0,1,2,3,4}), ISOLDE's own convention for
            "physical, non-restraint, non-map" forces (see module docstring).
        logger: optional ChimeraX-style logger (anything with a .warning(str)
            method, e.g. session.logger) used to report force types found within
            the included groups that this calculator doesn't know how to
            reconstruct (e.g. GBSAForce, AmberCMAPForce). Falls back to
            warnings.warn() if not given.
    '''

    _RECOGNISED = ('HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce',
                   'NonbondedForce')
    # ISOLDE's NonbondedSoftcoreForce/NonbondedSoftcoreExceptionForce are *pure
    # Python* subclasses of CustomNonbondedForce/CustomBondForce (custom_forces.py)
    # with no distinct C++ type -- OpenMM's bindings hand them back from
    # system.getForces() typed as the generic base class, not the ISOLDE subclass,
    # so class-name dispatch alone can never see "NonbondedSoftcoreForce". They are
    # identified instead by the 'softcore_lambda' global parameter the class always
    # defines -- nothing else in ISOLDE's system shares that class *and* group.
    _SOFTCORE_MARKER = 'softcore_lambda'

    def __init__(self, system, force_groups=None, logger=None):
        self.system = system
        self.n = system.getNumParticles()
        self.force_groups = set(CORE_FORCE_GROUPS if force_groups is None else force_groups)
        self._logger = logger
        self._bonds = None  # (M,2) idx, (M,) r0, (M,) k
        self._angles = None  # (M,3) idx, (M,) theta0, (M,) k
        self._torsions = None  # (M,4) idx, (M,) per, (M,) phase, (M,) k
        self._nb = None  # dict of nonbonded params (general pairwise term)
        self._nb_exceptions = None  # dict of exception/softcore-bond params
        self._nb_reaction_field = False
        self.found_force_groups = {}  # class name -> force-group bitmask (for validation)
        self._parse()

    # ------------------------------------------------------------------ parsing
    def _warn(self, message):
        if self._logger is not None:
            self._logger.warning(message)
        else:
            warnings.warn(message)

    @staticmethod
    def _has_global_parameter(force, param_name):
        for i in range(force.getNumGlobalParameters()):
            if force.getGlobalParameterName(i) == param_name:
                return True
        return False

    def _parse(self):
        nb_force = None
        nb_softcore_force = None
        nb_softcore_exceptions = None
        for force in self.system.getForces():
            if force.getForceGroup() not in self.force_groups:
                continue
            name = force.__class__.__name__
            found_as = None
            if name in self._RECOGNISED:
                found_as = name
            elif name == 'CustomNonbondedForce' and self._has_global_parameter(
                    force, self._SOFTCORE_MARKER):
                found_as = 'NonbondedSoftcoreForce'
                nb_softcore_force = force
            elif name == 'CustomBondForce' and self._has_global_parameter(
                    force, self._SOFTCORE_MARKER):
                found_as = 'NonbondedSoftcoreExceptionForce'
                nb_softcore_exceptions = force
            if found_as is None:
                self._warn(
                    f'LocalVirialCalculator: force "{name}" is in an included force '
                    'group but is not a recognised physical term. Its contribution '
                    'to the per-atom virial will be omitted.')
                continue
            self.found_force_groups.setdefault(found_as, 0)
            self.found_force_groups[found_as] |= (1 << force.getForceGroup())
            if found_as == 'HarmonicBondForce':
                self._parse_bonds(force)
            elif found_as == 'HarmonicAngleForce':
                self._parse_angles(force)
            elif found_as == 'PeriodicTorsionForce':
                self._parse_torsions(force)
            elif found_as == 'NonbondedForce':
                nb_force = force
            # NonbondedSoftcoreForce/NonbondedSoftcoreExceptionForce are captured
            # into the outer variables above; parsed together once the loop ends,
            # since the exception force is optional and may be seen either before
            # or after the general pairwise force.
        # ISOLDE's live sim always uses the softcore pair; fall back to a plain
        # NonbondedForce if that's what's actually present (e.g. an older/
        # non-default configuration).
        if nb_softcore_force is not None:
            self._parse_nonbonded_softcore(nb_softcore_force, nb_softcore_exceptions)
        elif nb_force is not None:
            self._parse_nonbonded(nb_force)

    def _parse_bonds(self, force):
        idx, r0, k = [], [], []
        for i in range(force.getNumBonds()):
            a, b, length, kk = force.getBondParameters(i)
            idx.append((a, b))
            r0.append(_strip(length, unit.nanometer))
            k.append(_strip(kk, unit.kilojoule_per_mole / unit.nanometer**2))
        self._bonds = (np.array(idx, dtype=int).reshape(-1, 2), np.array(r0), np.array(k))

    def _parse_angles(self, force):
        idx, t0, k = [], [], []
        for i in range(force.getNumAngles()):
            a, b, c, angle, kk = force.getAngleParameters(i)
            idx.append((a, b, c))
            t0.append(_strip(angle, unit.radian))
            k.append(_strip(kk, unit.kilojoule_per_mole / unit.radian**2))
        self._angles = (np.array(idx, dtype=int).reshape(-1, 3), np.array(t0), np.array(k))

    def _parse_torsions(self, force):
        idx, per, phase, k = [], [], [], []
        for i in range(force.getNumTorsions()):
            a, b, c, d, n, ph, kk = force.getTorsionParameters(i)
            idx.append((a, b, c, d))
            per.append(float(n))
            phase.append(_strip(ph, unit.radian))
            k.append(_strip(kk, unit.kilojoule_per_mole))
        self._torsions = (np.array(idx, dtype=int).reshape(-1, 4), np.array(per),
                          np.array(phase), np.array(k))

    def _parse_nonbonded(self, force):
        '''Plain openmm.NonbondedForce -- fallback path, not ISOLDE's live default.'''
        n = force.getNumParticles()
        charge = np.empty(n)
        sigma = np.empty(n)
        epsilon = np.empty(n)
        for i in range(n):
            q, s, e = force.getParticleParameters(i)
            charge[i] = _strip(q, unit.elementary_charge)
            sigma[i] = _strip(s, unit.nanometer)
            epsilon[i] = _strip(e, unit.kilojoule_per_mole)
        exc_idx, exc_q, exc_sig, exc_eps = [], [], [], []
        excluded = set()
        for i in range(force.getNumExceptions()):
            p, qd, qprod, s, e = force.getExceptionParameters(i)
            excluded.add((min(p, qd), max(p, qd)))
            exc_idx.append((p, qd))
            exc_q.append(_strip(qprod, unit.elementary_charge**2))
            exc_sig.append(_strip(s, unit.nanometer))
            exc_eps.append(_strip(e, unit.kilojoule_per_mole))
        self._set_nb_common(
            force, charge, sigma, epsilon, excluded, supports_reaction_field=True)
        self._nb_exceptions = dict(
            idx=np.array(exc_idx, dtype=int).reshape(-1, 2), q=np.array(exc_q),
            sig=np.array(exc_sig), eps=np.array(exc_eps))

    def _parse_nonbonded_softcore(self, force, exceptions_force):
        '''
        ISOLDE's NonbondedSoftcoreForce (CustomNonbondedForce) + optional
        NonbondedSoftcoreExceptionForce (CustomBondForce). Every original
        NonbondedForce exception (including nonzero-epsilon 1-4 pairs) is moved
        wholesale onto the exceptions force by
        SimHandler._convert_to_soft_core_potentials(), and fully excluded from
        the general pairwise force -- so, unlike plain NonbondedForce, there is
        no "nonzero exception that must still pass through the general loop"
        case to handle here.
        '''
        n = force.getNumParticles()
        charge = np.empty(n)
        sigma = np.empty(n)
        epsilon = np.empty(n)
        for i in range(n):
            q, s, e = force.getParticleParameters(i)
            charge[i] = _strip(q, unit.elementary_charge)
            sigma[i] = _strip(s, unit.nanometer)
            epsilon[i] = _strip(e, unit.kilojoule_per_mole)
        excluded = set()
        for i in range(force.getNumExclusions()):
            p, qd = force.getExclusionParticles(i)
            excluded.add((min(p, qd), max(p, qd)))
        self._set_nb_common(
            force, charge, sigma, epsilon, excluded, supports_reaction_field=False)
        exc_idx, exc_q, exc_sig, exc_eps = [], [], [], []
        if exceptions_force is not None:
            for i in range(exceptions_force.getNumBonds()):
                p1, p2, params = exceptions_force.getBondParameters(i)
                qprod, s, e = params
                exc_idx.append((p1, p2))
                exc_q.append(_strip(qprod, unit.elementary_charge**2))
                exc_sig.append(_strip(s, unit.nanometer))
                exc_eps.append(_strip(e, unit.kilojoule_per_mole))
        self._nb_exceptions = dict(
            idx=np.array(exc_idx, dtype=int).reshape(-1, 2), q=np.array(exc_q),
            sig=np.array(exc_sig), eps=np.array(exc_eps))

    def _set_nb_common(self, force, charge, sigma, epsilon, excluded, supports_reaction_field):
        '''
        Shared cutoff/switching setup for both NonbondedForce and
        NonbondedSoftcoreForce -- both expose the same generic Force API for
        this (getNonbondedMethod/getCutoffDistance/getUseSwitchingFunction/
        getSwitchingDistance), so this only needs to run once per class.

        `supports_reaction_field`: True for plain openmm.NonbondedForce, where a
        cutoff method (CutoffNonPeriodic/CutoffPeriodic) implicitly adds
        reaction-field electrostatic screening as part of OpenMM's own kernel.
        False for NonbondedSoftcoreForce (CustomNonbondedForce): its cutoff
        methods only spatially truncate/switch the user-specified energy
        expression -- there is no implicit reaction-field term, since the
        softcore Coulomb formula (and its lambda=1 classical limit) has none.
        '''
        method = force.getNonbondedMethod()
        no_cutoff_value = mm.NonbondedForce.NoCutoff
        has_cutoff = (method != no_cutoff_value)
        cutoff = None
        krf = crf = 0.0
        reaction_field = supports_reaction_field and has_cutoff
        if has_cutoff:
            cutoff = _strip(force.getCutoffDistance(), unit.nanometer)
        if reaction_field:
            eps_s = force.getReactionFieldDielectric() if hasattr(
                force, 'getReactionFieldDielectric') else 78.5
            krf = cutoff**-3 * (eps_s - 1.0) / (2.0 * eps_s + 1.0)
            crf = (1.0 / cutoff) * (3.0 * eps_s) / (2.0 * eps_s + 1.0)
        use_switch = force.getUseSwitchingFunction()
        switch_dist = (_strip(force.getSwitchingDistance(), unit.nanometer)
                       if use_switch else None)
        self._nb_reaction_field = reaction_field
        self._nb = dict(
            charge=charge, sigma=sigma, epsilon=epsilon, excluded=excluded,
            method=method, cutoff=cutoff, krf=krf, crf=crf, use_switch=use_switch,
            switch_dist=switch_dist)

    # -------------------------------------------------------------- accumulation
    def _accumulate(self, W, F, idx, forces, pos):
        '''idx:(M,n) atom indices, forces:(M,n,3), pos:(N,3). In-place add.'''
        if idx.size == 0:
            return
        n = idx.shape[1]
        pos_per = pos[idx]  # (M,n,3)
        rel = pos_per - pos_per[:, 0:1, :]  # (M,n,3) origin-independent
        T = np.einsum('mki,mkj->mij', rel, forces)  # (M,3,3)
        contrib = T / n
        for k in range(n):
            np.add.at(W, idx[:, k], contrib)
            np.add.at(F, idx[:, k], forces[:, k, :])

    # ------------------------------------------------------------------- forces
    def _bond_forces(self, pos):
        idx, r0, k = self._bonds
        d = pos[idx[:, 1]] - pos[idx[:, 0]]  # getDeltaR(a,b) = b - a
        r = np.linalg.norm(d, axis=1)
        dEdR = np.where(r > 0, k * (r - r0) / r, 0.0)
        f_a = dEdR[:, None] * d  # force on a (toward b if stretched)
        forces = np.stack([f_a, -f_a], axis=1)  # (M,2,3)
        return idx, forces

    def _angle_forces(self, pos):
        idx, t0, k = self._angles
        a, b, c = idx[:, 0], idx[:, 1], idx[:, 2]
        d0 = pos[b] - pos[a]  # getDeltaR(a,b) = b - a
        d1 = pos[b] - pos[c]  # getDeltaR(c,b) = b - c
        r2_0 = np.einsum('mi,mi->m', d0, d0)
        r2_1 = np.einsum('mi,mi->m', d1, d1)
        p = np.cross(d0, d1)
        rp = np.maximum(np.sqrt(np.einsum('mi,mi->m', p, p)), 1.0e-6)
        dot = np.einsum('mi,mi->m', d0, d1)
        cosine = np.clip(dot / np.sqrt(r2_0 * r2_1), -1.0, 1.0)
        angle = np.arccos(cosine)
        dEdR = k * (angle - t0)  # dE/dtheta
        termA = dEdR / (r2_0 * rp)
        termC = -dEdR / (r2_1 * rp)
        f_a = termA[:, None] * np.cross(d0, p)
        f_c = termC[:, None] * np.cross(d1, p)
        f_b = -(f_a + f_c)
        forces = np.stack([f_a, f_b, f_c], axis=1)  # (M,3,3)
        return idx, forces

    def _torsion_forces(self, pos):
        idx, per, phase, k = self._torsions
        a, b, c, d = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]
        d0 = pos[a] - pos[b]  # getDeltaR(B,A) = A - B
        d1 = pos[c] - pos[b]  # getDeltaR(B,C) = C - B
        d2 = pos[c] - pos[d]  # getDeltaR(D,C) = C - D
        cp0 = np.cross(d0, d1)
        cp1 = np.cross(d1, d2)
        r2_1 = np.einsum('mi,mi->m', d1, d1)
        normBC = np.sqrt(r2_1)
        # robust signed dihedral (matches OpenMM's sign convention; see notes)
        y = np.einsum('mi,mi->m', np.cross(cp0, cp1), d1) / normBC
        x = np.einsum('mi,mi->m', cp0, cp1)
        phi = np.arctan2(y, x)
        deltaAngle = per * phi - phase
        dEdAngle = -k * per * np.sin(deltaAngle)
        normCross1 = np.einsum('mi,mi->m', cp0, cp0)
        normCross2 = np.einsum('mi,mi->m', cp1, cp1)
        ff0 = (-dEdAngle * normBC) / normCross1
        ff3 = (dEdAngle * normBC) / normCross2
        ff1 = np.einsum('mi,mi->m', d0, d1) / r2_1
        ff2 = np.einsum('mi,mi->m', d2, d1) / r2_1
        iF0 = ff0[:, None] * cp0
        iF3 = ff3[:, None] * cp1
        s = ff1[:, None] * iF0 - ff2[:, None] * iF3
        iF1 = iF0 - s
        iF2 = iF3 + s
        # accumulation signs from ReferenceProperDihedralBond
        forces = np.stack([iF0, -iF1, -iF2, iF3], axis=1)  # (M,4,3)
        return idx, forces

    def _nonbonded_pairs(self, pos):
        nb = self._nb
        n = self.n
        if nb['method'] == mm.NonbondedForce.NoCutoff:
            i, j = np.triu_indices(n, k=1)
        else:
            from scipy.spatial import cKDTree
            tree = cKDTree(pos)
            pairs = tree.query_pairs(nb['cutoff'], output_type='ndarray')
            if pairs.size == 0:
                i = j = np.empty(0, dtype=int)
            else:
                i, j = pairs[:, 0], pairs[:, 1]
        if len(nb['excluded']) and len(i):
            lo = np.minimum(i, j)
            hi = np.maximum(i, j)
            keys = lo.astype(np.int64) * n + hi
            exc = np.array([a * n + b for (a, b) in nb['excluded']], dtype=np.int64)
            mask = ~np.isin(keys, exc)
            i, j = i[mask], j[mask]
        return i, j

    def _nonbonded_forces(self, pos):
        nb = self._nb
        i, j = self._nonbonded_pairs(pos)
        idx_list, force_list = [], []
        if len(i):
            f_i = self._pair_force(
                pos, i, j, q=nb['charge'][i] * nb['charge'][j],
                sigma=0.5 * (nb['sigma'][i] + nb['sigma'][j]),
                eps=np.sqrt(nb['epsilon'][i] * nb['epsilon'][j]),
                reaction_field=self._nb_reaction_field, use_switch=nb['use_switch'])
            idx_list.append(np.stack([i, j], axis=1))
            force_list.append(np.stack([f_i, -f_i], axis=1))
        # exceptions (or, for the softcore case, every original exception): plain
        # LJ + Coulomb, no reaction field, no switch
        exc = self._nb_exceptions
        if exc is not None and exc['idx'].size:
            keep = (exc['q'] != 0.0) | (exc['eps'] != 0.0)
            if np.any(keep):
                a, b = exc['idx'][keep, 0], exc['idx'][keep, 1]
                f_a = self._pair_force(
                    pos, a, b, q=exc['q'][keep], sigma=exc['sig'][keep], eps=exc['eps'][keep],
                    reaction_field=False, use_switch=False)
                idx_list.append(np.stack([a, b], axis=1))
                force_list.append(np.stack([f_a, -f_a], axis=1))
        return idx_list, force_list

    def _pair_force(self, pos, i, j, q, sigma, eps, reaction_field, use_switch):
        '''Force on atom i for each pair (i,j). Mirrors ReferenceLJCoulombIxn.'''
        nb = self._nb
        d = pos[i] - pos[j]  # r_i - r_j
        r2 = np.einsum('mi,mi->m', d, d)
        r = np.sqrt(r2)
        invR = 1.0 / r
        # LJ (OpenMM folds 4*epsilon into its combined epsilon; eps here is the
        # physical well depth, so use 4*eps)
        four_eps = 4.0 * eps
        sig_over_r = sigma * invR
        sig6 = sig_over_r**6
        switchValue = np.ones_like(r)
        switchDeriv = np.zeros_like(r)
        if use_switch:
            sd = nb['switch_dist']
            cd = nb['cutoff']
            mask = r > sd
            t = np.where(mask, (r - sd) / (cd - sd), 0.0)
            switchValue = np.where(mask, 1 + t**3 * (-10 + t * (15 - t * 6)), 1.0)
            switchDeriv = np.where(mask, t * t * (-30 + t * (60 - t * 30)) / (cd - sd), 0.0)
        dEdR = switchValue * four_eps * (12.0 * sig6 - 6.0) * sig6
        # Coulomb
        ke_q = ONE_4PI_EPS0 * q
        if reaction_field:
            dEdR = dEdR + ke_q * (invR - 2.0 * nb['krf'] * r2)
        else:
            dEdR = dEdR + ke_q * invR
        dEdR = dEdR * invR * invR
        if use_switch:
            lj_energy = four_eps * (sig6 - 1.0) * sig6
            dEdR = dEdR - lj_energy * switchDeriv * invR
        # force on i is dEdR * (r_i - r_j) (note d already = r_i - r_j)
        return dEdR[:, None] * d

    # ------------------------------------------------------------------- public
    def compute(self, positions):
        '''Return dict with 'virial' (N,3,3) and 'forces' (N,3, kJ/mol/nm).'''
        pos = np.asarray(_strip(positions, unit.nanometer), dtype=float).reshape(-1, 3)
        W = np.zeros((self.n, 3, 3))
        F = np.zeros((self.n, 3))
        if self._bonds is not None:
            self._accumulate(W, F, *self._bond_forces(pos), pos=pos)
        if self._angles is not None:
            self._accumulate(W, F, *self._angle_forces(pos), pos=pos)
        if self._torsions is not None:
            self._accumulate(W, F, *self._torsion_forces(pos), pos=pos)
        if self._nb is not None:
            idx_list, force_list = self._nonbonded_forces(pos)
            for idx, forces in zip(idx_list, force_list):
                self._accumulate(W, F, idx, forces, pos)
        return dict(virial=W, forces=F)

    @staticmethod
    def scalar_measures(W, volume=1.0):
        '''
        Per-atom scalar invariants of the (symmetrised) virial.

        Returns dict of (N,) arrays:
            hydrostatic : Tr(W)/3 / volume   (>0 compressive, <0 tensile)
            principal   : (N,3) sorted eigenvalues / volume
            von_mises   : deviatoric magnitude / volume  (shape/shear strain proxy)
            max_shear   : (lambda_max - lambda_min)/2 / volume
        '''
        Ws = 0.5 * (W + np.transpose(W, (0, 2, 1)))
        ev = np.linalg.eigvalsh(Ws) / volume  # (N,3) ascending
        hydro = ev.mean(axis=1)
        s1, s2, s3 = ev[:, 0], ev[:, 1], ev[:, 2]
        vm = np.sqrt(0.5 * ((s1 - s2)**2 + (s2 - s3)**2 + (s3 - s1)**2))
        return dict(hydrostatic=hydro, principal=ev, von_mises=vm, max_shear=0.5 * (s3 - s1))
