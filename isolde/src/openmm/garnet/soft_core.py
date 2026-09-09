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
                             SymmetryAwareMixin, ONE_ON_4_PI_EPS0, COULOMB_DECOUPLE_POWER)

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


# ======================================================================================
# Per-atom repulsive wall (+ short-range Coulomb guard) -- the r10b-era functional form.
# ======================================================================================
# The wall's repulsive decay is a PER-ATOM ``bee`` (combined to a per-pair ``alph``), and
# only ``beta`` stays a global (``alpha`` is gone). The optional short-range Coulomb guard
# is a positive-definite, attractive-pair-only electrostatics correction that keeps the
# (now soft) wall from being overrun by a Coulomb sink. It is folded into the Coulomb block
# under the SAME soft-core radial floor and coupling as bare Coulomb -- so it softens for
# clash release and fades on decouple/symmetry in lockstep, and is never left un-faded
# (the "guard = electrostatics" design). At lambda=1 both reduce to the exact r10b form.

def _dexp_peratom_of(arg):
    '''Per-atom-wall double-exponential at reduced coordinate ``arg``. The repulsive decay
    ``alph`` is a per-pair intermediate built from per-particle/-bond ``bee`` (see
    :data:`_ALPH_DEF`); ``beta`` is global; ``sigma``/``epsilon`` per-particle/-bond.
    Mirrors ``garnet_core.openmm_build._DEXP_ENERGY_PERATOM`` with ``r/rm`` replaced by
    ``arg`` (the softened reduced coordinate). The expression text is identical for a
    CustomNonbondedForce (``bee1``/``sigma1``...) and a CustomBondForce (per-bond params).'''
    return ('sqrt(epsilon1*epsilon2)*('
            f'((beta*exp(alph))/(alph-beta))*exp(-alph*{arg})'
            f'-((alph*exp(beta))/(alph-beta))*exp(-beta*{arg})'
            ')')


# Per-pair repulsive exponent from the two atoms' per-atom decays (arithmetic mean, biased
# to the softer partner). Matches garnet_core.openmm_build._DEXP_DEFS_PERATOM's ``alph``.
_ALPH_DEF = 'alph = 0.5*((2^(1/6))*(bee1*sigma1+bee2*sigma2))'


class GarnetPeratomVdwMixin:
    '''
    Per-atom repulsive wall: overrides :meth:`_vdw_block` with the per-atom
    double-exponential (per-particle ``bee`` -> per-pair ``alph``; only ``beta`` global),
    and keeps the equal ``2b`` vdW/Coulomb floors. Compose **first** in the bases (like
    :class:`GarnetVdwMixin`) so its ``__init__`` runs after the base builds the energy +
    per-particle params. No guard here -- it inherits the base bare-Coulomb block (with the
    garnet ``2b`` Coulomb floor via the override below).
    '''
    def __init__(self, *args, dexp_beta=0.0, **kwargs):
        super().__init__(*args, **kwargs)
        # Added after softcore_lambda(0)/softcore_alpha(1) so LAMBDA_INDEX=0 is preserved.
        self.addGlobalParameter('beta', dexp_beta)
        self.addPerParticleParameter('bee')

    @classmethod
    def _vdw_floor_power(cls, b):
        return b * 2

    @classmethod
    def _coulomb_floor_power(cls, b):
        return b * 2

    @classmethod
    def _vdw_block(cls, a, b, c, lam):
        head = f'{lam}^(1/{a}) * ({_dexp_peratom_of("xsoft")})'
        defs = f'xsoft = {_xsoft_expr(lam, cls._vdw_floor_power(b), c)};{_ALPH_DEF}'
        return head, defs


class GarnetPeratomGuardVdwMixin(GarnetPeratomVdwMixin):
    '''
    Per-atom wall + short-range Coulomb guard. Adds the guard globals (``cgw``/``cgp``) and
    per-particle ``bg``, and overrides :meth:`_coulomb_block` to fold the guard into the
    electrostatics under the same soft-core radial floor and coupling scale as bare Coulomb.
    '''
    def __init__(self, *args, coulomb_guard_w=0.0, **kwargs):
        super().__init__(*args, **kwargs)
        from garnet_core.energy import COULOMB_GUARD_P
        self.addGlobalParameter('cgw', coulomb_guard_w)
        self.addGlobalParameter('cgp', COULOMB_GUARD_P)
        self.addPerParticleParameter('bg')

    @classmethod
    def _coulomb_block(cls, b, c, lam, cs):
        from garnet_core.energy import COULOMB_GUARD_LAMBDA as CGL
        floor = cls._coulomb_floor_power(b)
        # The singular 1/r of BOTH bare Coulomb and the guard is regularised by the SAME
        # Pham-Shirts floor, so at lambda=1 (floor->0) each is its exact 1/r form and as
        # lambda drops both soften together -- the guard can't become an un-releasable wall.
        soft = f'( 1 / ( softcore_alpha*(1-{lam})^({floor}) + r^{c} ) )^(1/{c})'
        bare = f'{ONE_ON_4_PI_EPS0} * charge1 * charge2 * {soft}'
        guard = f'{ONE_ON_4_PI_EPS0} * cg_pref * cg_gg * {soft}'
        # Guard intermediates (chained reverse-dependency order, per the guard brief):
        # pref = pseudo-Huber of the Coulomb sink; gg = Slater-overlap damping; the pair
        # exponent bgij = M_p(bg1, bg2). qq uses the soft-core force's per-particle charges.
        defs = (
            f'cg_pref = 0.5*{CGL!r}*(sqrt(cg_qq*cg_qq+cgw*cgw)-cg_qq);'
            'cg_gg = exp(-cg_xx)*(1+(11/16)*cg_xx+(3/16)*cg_xx^2+(1/48)*cg_xx^3);'
            'cg_qq = charge1*charge2;'
            'cg_xx = cg_bgij*r;'
            'cg_bgij = (0.5*(bg1^cgp+bg2^cgp))^(1/cgp)'
        )
        # cs (the coupling-scale prefix) multiplies the WHOLE electrostatics term, so the
        # guard fades in lockstep with the Coulomb it corrects on decouple.
        return f'{cs}( {bare} + {guard} )', defs


class GarnetPeratomNonbondedSoftcoreForce(GarnetPeratomVdwMixin, NonbondedSoftcoreForce):
    '''Plain (no groups) per-atom-wall garnet soft-core.'''
    pass


class NBGroupGarnetPeratomNonbondedSoftcoreForce(GarnetPeratomVdwMixin,
                                                 NBGroupNonbondedSoftcoreForce):
    '''Per-group per-atom-wall garnet soft-core.'''
    pass


class SymmetryAwareGarnetPeratomNonbondedSoftcoreForce(SymmetryAwareMixin,
                                                       GarnetPeratomVdwMixin,
                                                       NBGroupNonbondedSoftcoreForce):
    '''Per-group + symmetry-aware per-atom-wall garnet soft-core.'''
    pass


class GarnetPeratomGuardNonbondedSoftcoreForce(GarnetPeratomGuardVdwMixin,
                                               NonbondedSoftcoreForce):
    '''Plain (no groups) per-atom-wall + Coulomb-guard garnet soft-core.'''
    pass


class NBGroupGarnetPeratomGuardNonbondedSoftcoreForce(GarnetPeratomGuardVdwMixin,
                                                      NBGroupNonbondedSoftcoreForce):
    '''Per-group per-atom-wall + Coulomb-guard garnet soft-core.'''
    pass


class SymmetryAwareGarnetPeratomGuardNonbondedSoftcoreForce(SymmetryAwareMixin,
                                                            GarnetPeratomGuardVdwMixin,
                                                            NBGroupNonbondedSoftcoreForce):
    '''Per-group + symmetry-aware per-atom-wall + Coulomb-guard garnet soft-core.'''
    pass


class GarnetPeratomNonbondedSoftcoreExceptionForce(openmm.CustomBondForce):
    '''
    Per-atom-wall (+ optional guard) analogue of
    :class:`GarnetNonbondedSoftcoreExceptionForce` for the garnet 1-4 / gated-1-3 pairs.
    Per-bond ``bee1``/``bee2`` drive the per-pair ``alph``; ``beta`` is global. When the
    guard is on it also carries the scaled 1-4 guard (per-bond raw charge product
    ``q_raw_prod`` + scale ``cg_scale`` = coulomb14scale, and per-bond ``bg1``/``bg2``),
    so at ``softcore_lambda=1`` this reproduces the plain build's 1-4 dexp + 1-4 guard
    exactly. Faded by the global ``softcore_lambda`` only (no per-group coupling), like the
    AMBER exception force.
    '''
    def __init__(self, a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2, dexp_beta=0.0,
                 coulomb_guard_w=None):
        floor = b * 2
        has_guard = coulomb_guard_w is not None
        soft = f'( 1 / ( softcore_alpha*(1-softcore_lambda)^({floor}) + r^{c} ) )^(1/{c})'
        bare = f'{ONE_ON_4_PI_EPS0} * charge_prod * {soft}'
        parts = [
            'vdw + coulombic;',
            f'vdw = vdw_scale * softcore_lambda^(1/{a}) * ({_dexp_peratom_of("xsoft")});',
            f'xsoft = {_xsoft_expr("softcore_lambda", floor, c)};',
            _ALPH_DEF + ';',
        ]
        if has_guard:
            from garnet_core.energy import COULOMB_GUARD_LAMBDA as CGL
            guard = (f'cg_scale * {ONE_ON_4_PI_EPS0} * cg_pref * cg_gg * {soft}')
            parts.append(f'coulombic = {bare} + {guard};')
            parts.append(f'cg_pref = 0.5*{CGL!r}*(sqrt(cg_qq*cg_qq+cgw*cgw)-cg_qq);')
            parts.append('cg_gg = exp(-cg_xx)*(1+(11/16)*cg_xx+(3/16)*cg_xx^2+(1/48)*cg_xx^3);')
            parts.append('cg_qq = q_raw_prod;')
            parts.append('cg_xx = cg_bgij*r;')
            parts.append('cg_bgij = (0.5*(bg1^cgp+bg2^cgp))^(1/cgp)')
        else:
            parts.append(f'coulombic = {bare}')
        super().__init__(''.join(parts))
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_alpha', alpha)
        self.addGlobalParameter('beta', dexp_beta)
        per_bond = ['charge_prod', 'sigma1', 'sigma2', 'epsilon1', 'epsilon2',
                    'bee1', 'bee2', 'vdw_scale']
        if has_guard:
            from garnet_core.energy import COULOMB_GUARD_P
            self.addGlobalParameter('cgw', coulomb_guard_w)
            self.addGlobalParameter('cgp', COULOMB_GUARD_P)
            per_bond += ['q_raw_prod', 'cg_scale', 'bg1', 'bg2']
        for p in per_bond:
            self.addPerBondParameter(p)
        self.update_needed = False


def find_garnet_nonbonded_forces(system):
    '''
    Locate garnet's plain nonbonded forces in ``system``, variant-agnostically:
    ``(coulomb_nb, coulomb_idx, dexp_cnb, dexp_idx, extra_bond_forces, guard_forces)``.

    * ``coulomb`` -- the plain ``NonbondedForce`` (LJ zeroed).
    * ``dexp`` -- the double-exponential ``CustomNonbondedForce``, recognised for BOTH
      forms: the old uniform wall carries an ``alpha`` global, the per-atom wall carries a
      per-particle ``bee`` instead; both carry ``beta`` (and no guard-specific ``cgw``).
    * ``extra_bond_forces`` -- ``(force, index)`` for the scaled 1-4 / gated-1-3 dexp
      ``CustomBondForce`` companions (``w14vdw``/``w13vdw`` globals).
    * ``guard_forces`` -- ``(force, index)`` for the short-range Coulomb guard forces
      (the nonbonded ``cgw`` force + its 1-4 ``w14cg`` companion), or ``[]`` when off.

    ``coulomb``/``dexp`` are ``None`` if absent.
    '''
    coulomb = dexp = None
    coulomb_idx = dexp_idx = None
    extras = []
    guards = []
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if type(f) is openmm.NonbondedForce:
            coulomb, coulomb_idx = f, i
        elif isinstance(f, openmm.CustomNonbondedForce):
            names = {f.getGlobalParameterName(k) for k in range(f.getNumGlobalParameters())}
            pnames = {f.getPerParticleParameterName(k)
                      for k in range(f.getNumPerParticleParameters())}
            if 'cgw' in names:                          # short-range Coulomb guard force
                guards.append((f, i))
            elif 'beta' in names and ('alpha' in names or 'bee' in pnames):
                dexp, dexp_idx = f, i                   # old (alpha+beta) or new (beta+bee)
        elif isinstance(f, openmm.CustomBondForce):
            names = {f.getGlobalParameterName(k) for k in range(f.getNumGlobalParameters())}
            if 'w14cg' in names:                        # 1-4 companion of the guard
                guards.append((f, i))
            elif names & {'w14vdw', 'w13vdw'}:
                extras.append((f, i))
    return coulomb, coulomb_idx, dexp, dexp_idx, extras, guards


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
    coul = lam ** COULOMB_DECOUPLE_POWER * ONE_ON_4_PI_EPS0 * charge ** 2 * (
        1.0 / (softcore_alpha * (1.0 - lam) ** floor + r ** c)) ** (1.0 / c)
    return vdw, coul
