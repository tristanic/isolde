# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Soft-core / fade-out / decoupling forces for the garnet force field.

garnet plugs into ISOLDE's existing AMBER soft-core machinery unchanged: a single
``softcore_lambda`` (minimise/equil ramp) + a per-group ``nb_coupling_table`` fade
the double-exponential vdW, the Coulomb term and (via the reused GBSA soft-core
classes) the implicit-solvent term to zero together, kept balanced — electrostatics
fade as the square of the coupling (``COULOMB_DECOUPLE_POWER``), vdW linearly.

Only the vdW block differs from AMBER, so these classes just override
``_vdw_block`` (via :class:`GarnetVdwMixin`) on top of the existing
``NonbondedSoftcoreForce`` ladder; the Coulomb term, the ``pair_lambda`` group
coupling, the ``SymmetryAwareMixin`` and the GBSA soft-core forces are all reused
verbatim.

**The double-exponential is finite at r→0** (no LJ singularity), so it needs no
``(1-lambda)`` regularisation to stay bounded. But to reproduce the LJ soft-core's
*preferential* softening of the short-range wall (LJ lowers the wall first as
lambda drops, which a bare prefactor cannot), we apply the **same Pham–Shirts
radial transform** LJ uses to the dexp's reduced coordinate: with ``x = r/r0``,
``r0 = 2^(1/6)*(sigma1+sigma2)/2``,

    x_soft = ( softcore_alpha*(1-lambda)^(2b) + x^c )^(1/c)

and evaluate the dexp at ``x_soft`` with a linear well-depth prefactor
``lambda^(1/a)``. At ``lambda=1`` ⇒ ``x_soft = x`` ⇒ the exact garnet dexp
(equilibrium force field unchanged). As lambda drops the wall (small-x, steep
``exp(-alpha*x)``) softens while the tail feels mainly the prefactor; at
``lambda→0`` it fades to zero. The prefactor and the softening floor are *separate*
factors, so a future split of the (currently conflated) lambda into independent
"couple" and "soften" knobs drops in without restructuring.

**vdW/Coulomb floor balance (why 2b, not LJ's b vs 4b).** The purpose of lowering
lambda is to soften the low-r wall so a severe clash (e.g. backbone threaded through
a phenyl ring) can release -- so we *want* the wall to come down; we only need the
vdW term to stay **above** the Coulomb as it does, or an oppositely-charged pair
develops a net attractive contact sink and fuses. LJ softens the vdW floor at
``(1-lambda)^b`` and the Coulomb floor at ``(1-lambda)^(4b)`` -- Coulomb softens 4x
slower, fine when the vdW wall is a towering r^-12 power law that dominates anyway.
For the *finite* dexp wall that ordering is fatal: the wall caves while the Coulomb
stays near-singular, and the Coulomb wins at intermediate lambda. So garnet softens
**both floors at the same rate** -- vdW and Coulomb floor exponents equal at ``2b``
(= 4 at the default b=2), overridden via
:meth:`NonbondedSoftcoreForce._vdw_floor_power` /
:meth:`~NonbondedSoftcoreForce._coulomb_floor_power`. The wall then stays above the
Coulomb as both soften, so it can drop far enough to release a clash without
opposite-charge collapse. (LJ keeps its b / 4b defaults, untouched.)
'''
import math

from openmm import openmm

from ..custom_forces import (NonbondedSoftcoreForce, NBGroupNonbondedSoftcoreForce,
                             SymmetryAwareMixin, ONE_ON_4_PI_EPS0)

# Reduced-coordinate + softened-coordinate expression fragments, shared by the
# CustomNonbondedForce mixin and the CustomBondForce exception below.
_R0 = '((2^(1/6))*((sigma1+sigma2)/2))'


def _xsoft_expr(lam, b, c):
    '''Pham–Shirts softened reduced coordinate x_soft as an OpenMM expression string.'''
    return f'( softcore_alpha*(1-{lam})^{b} + (r/{_R0})^{c} )^(1/{c})'


def _dexp_of(arg):
    '''garnet double-exponential evaluated at reduced coordinate ``arg`` (a variable
    name or expression); uses globals ``alpha``/``beta`` and per-particle/-bond
    ``sigma``/``epsilon``. Matches ``garnet_core.openmm_build._DEXP_ENERGY`` with
    ``r/r0`` replaced by ``arg``.'''
    return ('sqrt(epsilon1*epsilon2)*('
            f'((beta*exp(alpha))/(alpha-beta))*exp(-alpha*{arg})'
            f'-((alpha*exp(beta))/(alpha-beta))*exp(-beta*{arg})'
            ')')


class GarnetVdwMixin:
    '''
    Overrides :meth:`NonbondedSoftcoreForce._vdw_block` with the softened
    double-exponential, and registers the dexp ``alpha``/``beta`` globals. Compose
    **first** in the bases (like :class:`SymmetryAwareMixin`) so its ``__init__``
    runs after the base has built the energy + per-particle params.
    '''
    def __init__(self, *args, dexp_alpha=0.0, dexp_beta=0.0, **kwargs):
        super().__init__(*args, **kwargs)
        # Added after the base built the CustomNonbondedForce; the energy string may
        # reference globals added later (validated only at Context creation). Keeping
        # them after softcore_lambda(0)/softcore_alpha(1) preserves LAMBDA_INDEX=0.
        self.addGlobalParameter('alpha', dexp_alpha)
        self.addGlobalParameter('beta', dexp_beta)

    # Equal vdW/Coulomb floor exponents (2b = 4 at the default b=2). The finite dexp
    # wall would otherwise cave (vdW floor (1-l)^b) while the Coulomb stayed near-
    # singular (floor (1-l)^4b), letting the softened Coulomb overwhelm the wall and
    # fuse oppositely-charged atoms at intermediate lambda. Softening both floors at
    # the SAME rate keeps the vdW term above the Coulomb as lambda drops, so the wall
    # can still soften for clash release without an opposite-charge contact collapse.
    @classmethod
    def _vdw_floor_power(cls, b):
        return b * 2

    @classmethod
    def _coulomb_floor_power(cls, b):
        return b * 2

    @classmethod
    def _vdw_block(cls, a, b, c, lam):
        head = f'{lam}^(1/{a}) * ({_dexp_of("xsoft")})'
        defs = f'xsoft = {_xsoft_expr(lam, cls._vdw_floor_power(b), c)}'
        return head, defs


class GarnetNonbondedSoftcoreForce(GarnetVdwMixin, NonbondedSoftcoreForce):
    '''Plain (no groups) garnet soft-core dexp + softened Coulomb.'''
    pass


class NBGroupGarnetNonbondedSoftcoreForce(GarnetVdwMixin, NBGroupNonbondedSoftcoreForce):
    '''Per-group garnet soft-core (the `isolde decouple` / group-coupling variant).'''
    pass


class SymmetryAwareGarnetNonbondedSoftcoreForce(SymmetryAwareMixin, GarnetVdwMixin,
                                                NBGroupNonbondedSoftcoreForce):
    '''Per-group + crystallographic-symmetry-aware garnet soft-core.'''
    pass


class GarnetNonbondedSoftcoreExceptionForce(openmm.CustomBondForce):
    '''
    garnet 1-4 (and gated coordination 1-3) exceptions as soft-core bonds, faded by
    the global ``softcore_lambda`` only (no per-group coupling, exactly like the AMBER
    :class:`NonbondedSoftcoreExceptionForce`). Per-bond ``vdw_scale`` carries garnet's
    ``vdw14scale`` (1-4) or ``vdw13scale`` (gated 1-3); ``charge_prod`` carries the
    ``coulomb14scale``-scaled product (1-4) or 0.0 (1-3, which adds no Coulomb back).
    At ``softcore_lambda=1`` this reduces to garnet's plain 1-4/1-3 terms exactly.
    '''
    def __init__(self, a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2,
                 dexp_alpha=0.0, dexp_beta=0.0):
        # Equal vdW/Coulomb (1-lambda) floor exponents (2b = 4 at default b=2), matching
        # GarnetVdwMixin -- the finite dexp wall stays above the Coulomb as lambda drops.
        floor = b * 2
        energy = (
            'vdw + coulombic;'
            f'vdw = vdw_scale * softcore_lambda^(1/{a}) * ({_dexp_of("xsoft")});'
            f'xsoft = {_xsoft_expr("softcore_lambda", floor, c)};'
            f'coulombic = {ONE_ON_4_PI_EPS0} * charge_prod * '
                f'( 1 / ( softcore_alpha*(1-softcore_lambda)^({floor}) + r^{c} ) )^(1/{c})'
        )
        super().__init__(energy)
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_alpha', alpha)
        self.addGlobalParameter('alpha', dexp_alpha)
        self.addGlobalParameter('beta', dexp_beta)
        for p in ('charge_prod', 'sigma1', 'sigma2', 'epsilon1', 'epsilon2', 'vdw_scale'):
            self.addPerBondParameter(p)
        self.update_needed = False


def find_garnet_nonbonded_forces(system):
    '''
    Locate garnet's plain nonbonded forces in ``system``:
    ``(coulomb_nb, coulomb_idx, dexp_cnb, dexp_idx, extra_bond_forces)`` where
    ``extra_bond_forces`` is a list of ``(force, index)`` for the scaled 1-4 / gated
    1-3 dexp ``CustomBondForce`` companions (identified by their ``w14vdw``/``w13vdw``
    globals). ``coulomb``/``dexp`` are ``None`` if absent.
    '''
    coulomb = dexp = None
    coulomb_idx = dexp_idx = None
    extras = []
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if type(f) is openmm.NonbondedForce:
            coulomb, coulomb_idx = f, i
        elif isinstance(f, openmm.CustomNonbondedForce):
            names = {f.getGlobalParameterName(k) for k in range(f.getNumGlobalParameters())}
            if {'alpha', 'beta'} <= names:
                dexp, dexp_idx = f, i
        elif isinstance(f, openmm.CustomBondForce):
            names = {f.getGlobalParameterName(k) for k in range(f.getNumGlobalParameters())}
            if names & {'w14vdw', 'w13vdw'}:
                extras.append((f, i))
    return coulomb, coulomb_idx, dexp, dexp_idx, extras


# Representative parameters for the illustrative potential-vs-radius plot (mirrors the
# LJ mirror's hard-coded O–O choice). Not used in any simulation.
_PLOT_SIGMA = 0.30        # nm
_PLOT_EPSILON = 0.50      # kJ/mol


def potential_values(radii, lam, a, b, c, softcore_alpha, dexp_alpha=12.24, dexp_beta=4.37,
                     charge=-0.1):
    '''
    NumPy mirror of the garnet soft-core vdW + Coulomb pairwise energy for the GUI
    potential indicator (`ui/general_tab/nonbonded.py`). ``radii`` in nm; returns
    ``(vdw, coulomb)`` in kJ/mol for a representative like-atom pair. Same math as the
    OpenMM expressions above, so the plotted curve matches what the simulation uses
    (shape + how lambda softens the wall).

    The Coulomb term is **identical** to AMBER's
    :meth:`NonbondedSoftcoreForce.potential_values` (same softened form, same
    representative ``charge``, same ``ONE_ON_4_PI_EPS0``) — garnet and AMBER share the
    soft-core Coulomb verbatim, so the two curves must overlay. Only the vdW form
    differs (dexp vs L-J).
    '''
    import numpy
    r = numpy.asarray(radii, dtype=float)
    al, be = float(dexp_alpha), float(dexp_beta)
    floor = b * 2       # equal vdW/Coulomb (1-lambda) floor exponents (see GarnetVdwMixin)
    r0 = (2.0 ** (1.0 / 6.0)) * _PLOT_SIGMA        # like-atom r0
    x = r / r0
    xsoft = (softcore_alpha * (1.0 - lam) ** floor + x ** c) ** (1.0 / c)
    dexp = _PLOT_EPSILON * (((be * math.exp(al)) / (al - be)) * numpy.exp(-al * xsoft)
                            - ((al * math.exp(be)) / (al - be)) * numpy.exp(-be * xsoft))
    vdw = (lam ** (1.0 / a)) * dexp
    coul = ONE_ON_4_PI_EPS0 * charge ** 2 * (
        1.0 / (softcore_alpha * (1.0 - lam) ** floor + r ** c)) ** (1.0 / c)
    return vdw, coul
