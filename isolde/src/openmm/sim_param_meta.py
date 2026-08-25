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
                 step=None, choices=None, decimals=None):
        self.kind = kind
        self.label = label
        self.tooltip = tooltip
        self.minimum = minimum
        self.maximum = maximum
        self.step = step
        # decimals: spin-box precision for a float editor (None -> derived from step)
        self.decimals = decimals
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


# The curated set of **force-field-specific** tunable parameters -- and ONLY these.
#
# Membership here is what makes a parameter per-force-field: it is remembered
# per force field, re-resolved when the force field changes, exposed by the
# ``isolde simparam`` command and the tuning dashboard, and included in the
# saved store. **Every other SimParams entry is session-wide** and is never
# touched by a force-field switch -- crucially the hardware/UX ones (platform,
# device index, trajectory smoothing, GUI cadence) and the physical conditions
# and interaction knobs (temperature, cutoff, GBSA on/off, symmetry, rigid water,
# the soft-core on/off toggle, restraint spring constants) which have their own
# controls and do not depend on the force field's functional form.
#
# Only put a parameter here if its *optimal value genuinely differs because of the
# force field* (e.g. GARNET's bounded vdW makes flexible hydrogens viable, wanting
# a different rigid_bonds / HMR / friction than AMBER).
def _build_meta():
    return {
        'nonbonded_softcore_lambda_minimize': ParamSpec('float', 'Soft-core λ (min)',
            'Soft-core coupling during minimisation. Lower softens the wall for '
            'clash escape.', minimum=0.01, maximum=1.0, step=0.01),
        'nonbonded_softcore_lambda_equil': ParamSpec('float', 'Soft-core λ (equil)',
            'Soft-core coupling during dynamics. 1.0 = the exact potential.',
            minimum=0.01, maximum=1.0, step=0.01),
        'nonbonded_softcore_alpha': ParamSpec('float', 'Soft-core α',
            'Shape of the soft-core wall within the clashing region.',
            minimum=0.01, maximum=1.0, step=0.01),
        'rigid_bonds': ParamSpec('enum', 'Rigid bonds',
            'Which bonds are held rigid by constraints. "none" = fully flexible '
            '(needs higher friction / HMR for stability).',
            choices=_rigid_bonds_choices()),
        'hmr_hydrogen_mass': ParamSpec('float', 'HMR H mass (amu)',
            'Hydrogen mass repartitioning target (0 = off). Raises H mass toward '
            'this, taking it from the bonded heavy atom, floored per centre.',
            minimum=0.0, maximum=6.0, step=0.5),
        'friction_coefficient': ParamSpec('float', 'Friction (/ps)',
            'Langevin friction. Higher damps motion (suppresses "ringing" of the '
            'fast modes exposed by flexible hydrogens) at the cost of responsiveness.',
            minimum=0.1, maximum=100, step=1),
        'variable_integrator_tolerance': ParamSpec('float', 'Integrator tolerance',
            'Variable-timestep integrator error tolerance. Tightening (e.g. 1e-6) '
            'improves stability -- helpful for flexible hydrogens -- at a modest '
            'performance cost.', minimum=1e-7, maximum=1e-2, step=1e-6, decimals=7),
    }


PARAM_META = _build_meta()

# The authoritative set of per-force-field parameter names (== the tunable set).
# SimParams consults this to decide which params are remembered per force field and
# re-resolved on a switch; everything else is session-wide.
PER_FF_PARAMS = frozenset(PARAM_META)


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
