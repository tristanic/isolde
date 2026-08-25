# @Author: Tristan Croll
# @Date:   25-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 25-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Type / editing metadata for the *tunable* subset of :class:`SimParams`, plus the
value (de)serialisation used by the ``isolde simparam`` command and the
explicit-save persistence store.

This is the single source of truth for "what kind is this parameter and how do I
turn it to/from text or a JSON-able token" -- consumed by the command's value
parser, by the persistence layer, and (later) by the tuning dashboard's editor
widgets. Only parameters listed here are exposed for tuning; the rest of
``SimParams`` is left alone.

A parameter's *unit* is not stored here -- it is read from
``SimParams._default_params[key][1]`` so there is one source for units. ``kind`` is
one of ``'float' | 'int' | 'bool' | 'enum'`` (a unit-bearing float is still
``'float'``; the unit is applied by ``set_param``/read for serialisation).
'''


class ParamSpec:
    def __init__(self, kind, label, tooltip='', minimum=None, maximum=None,
                 step=None, choices=None):
        self.kind = kind
        self.label = label
        self.tooltip = tooltip
        self.minimum = minimum
        self.maximum = maximum
        self.step = step
        # choices: ordered list of (token, value) for kind == 'enum'
        self.choices = list(choices) if choices else None

    def token_for(self, value):
        '''enum value -> token string.'''
        for tok, val in self.choices:
            if val is value or val == value:
                return tok
        raise ValueError('No token for value {!r}'.format(value))

    def value_for(self, token):
        '''token string -> enum value.'''
        for tok, val in self.choices:
            if tok == token:
                return val
        raise ValueError('Unknown choice {!r} (expected one of {})'.format(
            token, ', '.join(t for t, _ in self.choices)))


def _rigid_bonds_choices():
    from openmm.app import HBonds, AllBonds, HAngles
    return [('none', None), ('h-bonds', HBonds),
            ('all-bonds', AllBonds), ('h-angles', HAngles)]


# The curated tunable set. (Extend as needed; anything not here is simply not
# exposed by the command / dashboard.)
def _build_meta():
    return {
        'temperature': ParamSpec('float', 'Temperature (K)',
            'Simulation temperature.', minimum=0, maximum=500, step=1),
        'friction_coefficient': ParamSpec('float', 'Friction (/ps)',
            'Langevin friction. Higher damps motion (suppresses "ringing" of fast '
            'modes) at the cost of responsiveness.', minimum=0.1, maximum=100, step=1),
        'nonbonded_cutoff_distance': ParamSpec('float', 'Nonbonded cutoff (nm)',
            'Direct-space nonbonded cutoff.', minimum=0.5, maximum=2.0, step=0.1),
        'use_gbsa': ParamSpec('bool', 'Implicit solvent (GBSA)',
            'Whether to include GBSA implicit solvent.'),
        'use_softcore_nonbonded_potential': ParamSpec('bool', 'Soft-core nonbonded',
            'Use the soft-core (fade-out) nonbonded potential.'),
        'nonbonded_softcore_lambda_minimize': ParamSpec('float', 'Soft-core λ (min)',
            'Soft-core coupling during minimisation. Lower softens the wall for '
            'clash escape.', minimum=0.01, maximum=1.0, step=0.01),
        'nonbonded_softcore_lambda_equil': ParamSpec('float', 'Soft-core λ (equil)',
            'Soft-core coupling during dynamics. 1.0 = the exact potential.',
            minimum=0.01, maximum=1.0, step=0.01),
        'nonbonded_softcore_alpha': ParamSpec('float', 'Soft-core α',
            'Shape of the soft-core wall within the clashing region.',
            minimum=0.01, maximum=1.0, step=0.01),
        'symmetry_aware': ParamSpec('bool', 'Symmetry-aware',
            'Crystallographic symmetry copies participate in the simulation.'),
        'rigid_bonds': ParamSpec('enum', 'Rigid bonds',
            'Which bonds are held rigid by constraints. "none" = fully flexible '
            '(needs higher friction / HMR for stability).',
            choices=_rigid_bonds_choices()),
        'rigid_water': ParamSpec('bool', 'Rigid water',
            'Constrain water geometry (SETTLE).'),
        'hmr_hydrogen_mass': ParamSpec('float', 'HMR H mass (amu)',
            'Hydrogen mass repartitioning target (0 = off). Raises H mass toward '
            'this, taking it from the bonded heavy atom, floored per centre.',
            minimum=0.0, maximum=6.0, step=0.5),
    }


PARAM_META = _build_meta()


def spec_for(key):
    return PARAM_META.get(key)


def tunable_params():
    '''Ordered list of the tunable parameter names.'''
    return list(PARAM_META.keys())


# --- value parsing (command) ---------------------------------------------

def parse_value(key, text):
    '''
    Parse a user-supplied text value for parameter ``key`` per its ``kind``.
    Returns a value ready to hand to ``SimParams.set_param`` (a raw number for a
    unit-bearing float -- ``set_param`` attaches the unit). Raises ``ValueError``
    on a bad value.
    '''
    spec = PARAM_META[key]
    t = text.strip()
    if spec.kind == 'bool':
        low = t.lower()
        if low in ('true', 't', '1', 'yes', 'on'): return True
        if low in ('false', 'f', '0', 'no', 'off'): return False
        raise ValueError('Expected a boolean (true/false) for {}'.format(key))
    if spec.kind == 'int':
        return int(t)
    if spec.kind == 'enum':
        return spec.value_for(t.lower())
    v = float(t)
    if spec.minimum is not None and v < spec.minimum:
        raise ValueError('{} must be >= {}'.format(key, spec.minimum))
    if spec.maximum is not None and v > spec.maximum:
        raise ValueError('{} must be <= {}'.format(key, spec.maximum))
    return v


# --- (de)serialisation for the persistent store ---------------------------

def _serialize_value(sim_params, key, value):
    '''SimParams value -> JSON-able token for the settings store.'''
    from openmm.unit import Quantity
    spec = PARAM_META.get(key)
    if spec is not None and spec.kind == 'enum':
        return spec.token_for(value)
    if isinstance(value, Quantity):
        unit = sim_params._default_params[key][1]
        return float(value.value_in_unit(unit))
    return value                                       # bool / float / int / None


def _deserialize_value(sim_params, key, token):
    '''Settings-store token -> a value ready for ``set_param``.'''
    spec = PARAM_META.get(key)
    if spec is not None and spec.kind == 'enum':
        return spec.value_for(token)
    return token                                       # bool / float / int / None
                                                       # (raw float re-acquires its
                                                       #  unit via set_param)


def serialize_overrides(sim_params):
    '''
    ``{ff: {param: token}}`` capturing every force field's current live overrides
    (only the tunable, exposed params), ready to write to the settings store.
    '''
    out = {}
    for ff, params in sim_params._overrides.items():
        ser = {}
        for key, value in params.items():
            if key not in PARAM_META:
                continue                               # only persist exposed knobs
            ser[key] = _serialize_value(sim_params, key, value)
        if ser:
            out[ff] = ser
    return out


def apply_overrides(sim_params, store):
    '''
    Load a ``{ff: {param: token}}`` store into ``sim_params``' per-force-field
    override maps and re-apply the active force field's overrides to the live
    values. Unknown params / bad tokens are skipped (logged by the caller if
    desired), so a stale store can never block startup.
    '''
    if not store:
        return
    for ff, params in store.items():
        target = sim_params._overrides.setdefault(ff, {})
        for key, token in params.items():
            if key not in sim_params._default_params or key not in PARAM_META:
                continue
            try:
                target[key] = _deserialize_value(sim_params, key, token)
            except Exception:
                continue
    active = sim_params._params.get('forcefield')
    sim_params._reapply_overrides(active)


# --- explicit-save persistence (ChimeraX Settings backend) ----------------

def save_to_settings(session, sim_params):
    '''
    Persist every force field's current live overrides to the ISOLDE force-field
    settings store (the explicit "Save current settings" action). Returns True if
    written. The store is *not* registered in ChimeraX's Settings panel -- it is a
    silent persistence backend for the tuning dashboard.
    '''
    from .. import settings
    fs = getattr(settings, 'forcefield_settings', None)
    if fs is None:
        return False
    fs.ff_overrides = serialize_overrides(sim_params)
    fs.save()
    return True


def load_from_settings(session, sim_params):
    '''Load saved per-force-field overrides from the settings store into
    ``sim_params`` (applied to the active force field's live values). A no-op if
    the store is unset/empty. Never raises on a stale store.'''
    from .. import settings
    fs = getattr(settings, 'forcefield_settings', None)
    if fs is None:
        return
    try:
        store = dict(fs.ff_overrides or {})
    except Exception:
        return
    apply_overrides(sim_params, store)
