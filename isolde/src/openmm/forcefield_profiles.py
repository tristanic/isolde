# @Author: Tristan Croll
# @Date:   25-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 25-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Per-force-field **profiles**: each selectable force field's default
sim-parameter overrides and which parameters apply to it.

ISOLDE historically had a single flat set of simulation-parameter defaults
(:class:`~chimerax.isolde.openmm.sim_param_mgr.SimParams`). With more than one
force field, that no longer holds: different force fields want different default
values (e.g. GARNET's double-exponential vdW runs at soft-core ``lambda=1`` at
equilibrium where AMBER wants 0.95), and some parameters do not apply at all
under a given force field.

A profile is a small, name-keyed data object (the registry key matches the
force-field name used by
:class:`~chimerax.isolde.openmm.forcefields.ForcefieldMgr`). It deliberately does
**not** live on the force-field object itself, because AMBER has no ISOLDE-side
handle (it is a bare OpenMM ``ForceField``); a name-keyed registry is the only
uniform host. A force field with no registered profile gets an empty one, so its
behaviour is identical to the historical single-defaults path.

This module is intentionally dependency-light (no ChimeraX/OpenMM imports at
module load) so it can be imported from both the parameter manager and the UI
without circular-import risk.
'''


class ForceFieldProfile:
    '''
    Default sim-parameter overrides + applicability for one force field.

    Args:
        name:            the force-field name (matches ``sim_params.forcefield``).
        param_defaults:  ``{param_name: value}`` overriding the static
                         ``SimParams._default_params`` default for this force
                         field. Values are in the parameter's stored convention
                         (they pass through ``set_param``, so units are handled).
        not_applicable:  iterable of parameter names that are meaningless under
                         this force field (the UI hides them; the command rejects
                         them). Empty ⇒ every parameter applies.
    '''
    def __init__(self, name, param_defaults=None, not_applicable=()):
        self.name = name
        self.param_defaults = dict(param_defaults or {})
        self._not_applicable = set(not_applicable)

    def has_default(self, key):
        '''True if this profile overrides the static default for ``key``. (Use
        this rather than testing ``default_for(key) is None`` — ``None`` is a
        legitimate override value, e.g. ``rigid_bonds=None``.)'''
        return key in self.param_defaults

    def default_for(self, key):
        '''This force field's override value for ``key`` (KeyError if none —
        callers should guard with :meth:`has_default`).'''
        return self.param_defaults[key]

    def applies(self, key):
        '''True unless ``key`` is explicitly marked not-applicable here.'''
        return key not in self._not_applicable


# --- registry -------------------------------------------------------------

# The empty profile: no overrides, everything applicable. Returned for any force
# field without a registered profile, so e.g. AMBER is byte-for-byte the old
# single-defaults behaviour.
_EMPTY_PROFILE = ForceFieldProfile('_default')

_PROFILES = {
    # GARNET: the double-exponential vdW has no r->0 singularity, so at
    # equilibrium it runs at the exact dexp (soft-core lambda = 1), unlike AMBER's
    # conflated lambda (0.95). (This supersedes the former module-level
    # _FORCEFIELD_PARAM_DEFAULTS dict in sim_param_mgr.py.)
    'garnet': ForceFieldProfile(
        'garnet',
        param_defaults={
            'nonbonded_softcore_lambda_equil': 1.0,
        },
    ),
}


def get_profile(name):
    '''
    Return the :class:`ForceFieldProfile` for a force-field name, or an empty
    profile (no overrides, all parameters applicable) for any force field without
    a registered one — so an unconfigured force field behaves exactly as under the
    historical single set of defaults.
    '''
    return _PROFILES.get(name, _EMPTY_PROFILE)


def register_profile(profile):
    '''Register (or replace) a force-field profile. Keyed by ``profile.name``.'''
    _PROFILES[profile.name] = profile


def ff_sensitive_params():
    '''
    The set of parameter names that *any* registered profile overrides — i.e. the
    parameters whose default depends on the force field. These are the only
    parameters re-resolved when the force field changes; everything else is
    force-field-agnostic and left untouched by a switch.
    '''
    keys = set()
    for prof in _PROFILES.values():
        keys.update(prof.param_defaults.keys())
    return keys
