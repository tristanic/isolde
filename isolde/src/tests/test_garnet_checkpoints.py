# @Author: Tristan Croll
# @Email: tcroll@altoslabs.com
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless check that ISOLDE can find and use the checkpoints GARNET ships.

WHY. ISOLDE used to pin each ``garnet-{run}`` variant as a path *into garnet's repo*
(``garnetff/trained_models/<name>.pt``) and resolve it by walking up out of the
``garnet_core`` package to a sibling directory. That holds in a repo checkout and by
luck in a plain site-packages install, and fails for an installed wheel that owns its
data, for ``pip --target``, and for anything relocated. Since the wheel is how the force
field now reaches other people, the layout assumption had to go: ISOLDE pins a BARE
FILENAME and ``garnet_core.weights`` resolves it, knowing both layouts.

So what is under test is the seam, not the model: every offered variant resolves to a
file that exists, whichever way garnet was installed. The last check then actually
parameterises a structure, because a resolvable path that cannot be loaded is no use.

Skips cleanly (reporting ALL PASS) when garnet is not installed -- it is opt-in, and
most ISOLDE installs will not have it.
'''
import os

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.isolde.openmm import forcefields as ff
    from chimerax.isolde.openmm.garnet import params_cache as pc

    variants = ff._GARNET_VARIANTS

    # --- the pins are bare filenames, not paths into garnet's source layout ----
    for name, ckpt in variants.items():
        if os.sep in ckpt or '/' in ckpt:
            _fail('variant %s pins a PATH (%r); it must pin a bare filename so garnet '
                  'owns the layout' % (name, ckpt))
    print('PASS: all %d variant pins are bare filenames' % len(variants))

    # --- is garnet actually installed? -----------------------------------------
    try:
        pc.require_garnet_core()
    except ImportError as e:
        # The guard must say what to install; a bare "No module named" is the defect.
        if 'garnet-isolde' not in str(e) or 'pip install' not in str(e):
            _fail('missing-garnet error is not actionable: %s' % e)
        print('PASS: garnet_core absent, and the error names the wheel and how to install it')
        print('SKIP: resolution/parameterisation checks need garnet installed')
        print('ALL PASS')
        return

    import garnet_core
    print('PASS: garnet_core imported from %s' % os.path.dirname(garnet_core.__file__))

    # --- every offered variant resolves to a real file -------------------------
    for name, ckpt in sorted(variants.items()):
        path = pc.resolve_checkpoint_path(ckpt)
        if not os.path.isfile(path):
            _fail('variant %s (%s) resolved to a non-existent file: %s' % (name, ckpt, path))
        print('PASS: %-14s -> %s' % (name, path))

    # the bare `garnet` alias / default
    default = pc.default_checkpoint_path()
    if not os.path.isfile(default):
        _fail('default checkpoint does not exist: %s' % default)
    print('PASS: default checkpoint -> %s' % default)

    # --- a legacy repo-relative pin must still resolve -------------------------
    # Sessions and scripts saved before the seam changed carry the old form.
    for name, ckpt in sorted(variants.items()):
        legacy = os.path.join('garnetff', 'trained_models', ckpt)
        if pc.resolve_checkpoint_path(legacy) != pc.resolve_checkpoint_path(ckpt):
            _fail('legacy pin %r does not resolve to the same file as %r' % (legacy, ckpt))
    print('PASS: legacy repo-relative pins still resolve identically')

    # --- an absolute path is still honoured verbatim ---------------------------
    if pc.resolve_checkpoint_path(default) != default:
        _fail('an absolute checkpoint path must be used verbatim')
    print('PASS: absolute checkpoint paths honoured verbatim')

    # --- end to end: parameterise a real structure -----------------------------
    # A resolvable path that will not load is no use, so actually run garnet. The
    # checkpoint is passed as the BARE PIN, so this exercises the resolution seam
    # through the same entry point a simulation start uses.
    import math
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde.openmm.garnet import get_garnet_parameters

    m = run_cmd(session, 'open %s' % FIXTURE)[0]
    try:
        run_cmd(session, 'addh')        # garnet needs explicit hydrogens
        newest = sorted(variants)[-1]
        params = get_garnet_parameters(m, checkpoint_path=variants[newest],
                                       logger=session.logger)
        if not params.parameterised:
            _fail('%s: parameterise() did not complete' % newest)
        if len(params._per_atom) != m.num_atoms:
            _fail('%s: parameterised %d of %d atoms'
                  % (newest, len(params._per_atom), m.num_atoms))
        bad = [a for a, (q, sig, eps) in params._per_atom.items()
               if not (math.isfinite(q) and math.isfinite(sig) and math.isfinite(eps))]
        if bad:
            _fail('%s: %d atoms with non-finite charge/sigma/epsilon' % (newest, len(bad)))
        if not (params._bonds and params._angles and params._propers):
            _fail('%s: bonded terms missing (bonds=%d angles=%d propers=%d)'
                  % (newest, len(params._bonds), len(params._angles), len(params._propers)))
        net = sum(q for q, _, _ in params._per_atom.values())
        print('PASS: %s parameterised %d atoms of %s (net q %+.3f e; '
              '%d bonds, %d angles, %d propers)'
              % (newest, len(params._per_atom), os.path.basename(FIXTURE), net,
                 len(params._bonds), len(params._angles), len(params._propers)))
    finally:
        session.models.close([m])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
