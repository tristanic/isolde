# @Author: Tristan Croll
# @Date:   20-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
isolde decouple -- interactively soften (decouple) a selection's nonbonded
interactions with the rest of the running simulation, using ISOLDE's per-group
soft-core coupling. A hands-on aid for exploring/verifying the machinery that
rotafit and the future placement engines drive internally.

The decoupling is TRANSIENT and SIM-SCOPED: it exists only for the currently
running simulation and is forgotten the moment that simulation stops (the
per-group state lives on the simulation's force objects, which are rebuilt fully
coupled on the next "isolde sim start"). Nothing is persisted to the session.
'''

# Default soft-core coupling for a decoupled selection (= the minimum of the allowed
# range). Values below 0.1 have no practical use (the selection's coupling to its
# surroundings is already effectively off), and the soft-core potential cannot switch
# off completely, so the range is simply [0.1, 1.0] -- 1.0 being normal full coupling.
DEFAULT_DECOUPLE_LAMBDA = 0.1


def decouple_lambda_arg():
    '''ChimeraX annotation for the coupling value: a float in [0.1, 1.0], bounded at the
    command level. Named ``lambda_decouple`` (CLI ``lambdaDecouple``) rather than
    ``lambda``: it avoids the Python reserved word, still abbreviates to ``lambda`` for
    the user, and keeps ``decouple`` an unambiguous keyword prefix.'''
    from chimerax.core.commands import Bounded, FloatArg
    return Bounded(FloatArg, 0.1, 1.0, inclusive=True, name='a coupling value in [0.1, 1]')


def decouple(session, atoms, state=True, lambda_decouple=DEFAULT_DECOUPLE_LAMBDA):
    '''Soften (``on``, the default) or restore (``off``) the nonbonded coupling of
    ``atoms`` to the rest of the running simulation. ``state`` is a boolean (``on``/
    ``off``/``true``/``false``). See :func:`register_decouple_command`.'''
    from chimerax.core.errors import UserError
    isolde = getattr(session, 'isolde', None)
    if isolde is None or not getattr(isolde, 'simulation_running', False):
        raise UserError('isolde decouple only works while an ISOLDE simulation is '
                        'running. Start one first (e.g. "isolde sim start sel").')
    sh = isolde.sim_handler
    # Per-group soft-core coupling must be available: the soft-core nonbonded potential
    # plus >= 3 provisioned group slots (env=0, selection=1, its symmetry copies=2).
    # Every sim provides these by default (SimParams.nb_groups_max=4).
    if not (getattr(sh, 'nb_groups_enabled', False)
            and getattr(sh, 'nb_groups_count', 1) >= 3):
        raise UserError(
            'isolde decouple needs per-group soft-core coupling, which is not '
            'available for this simulation (it needs the soft-core nonbonded potential '
            'and SimParams.nb_groups_max >= 3, the default). Restart the simulation with '
            'the defaults and retry.')

    if not state:
        # off: recouple EVERYTHING -- return all particles (and their symmetry copies) to
        # the environment group and reset the coupling table to full strength. The
        # selection argument is not used when turning off -- only one selection is ever
        # decoupled at a time (see the "on" path below), so this simply clears it.
        sh.assign_nb_group(sh._atoms, 0)
        n = sh.nb_groups_count
        for a in range(n):
            for b in range(a, n):
                sh.set_nb_coupling(a, b, 1.0)
        session.logger.info('isolde decouple: full coupling restored '
                            '(all selections recoupled).')
        return

    # on (the default). lambda_decouple is bounded to [0.1, 1] at the command level.
    n_in = int((sh._atoms.indices(atoms) != -1).sum()) if len(atoms) else 0
    if n_in == 0:
        raise UserError('isolde decouple: none of the selected atoms are part of the '
                        'running simulation (only atoms in the sim can be decoupled).')
    # Decoupling is EXCLUSIVE: reset any previous decoupling first so exactly this
    # selection is softened against everything else. soften_nb_selection puts the
    # selection in group 1 (its symmetry copies in group 2) and softens every coupling
    # touching it -- including its own crystal self-contact -- to lambda_decouple, while
    # keeping the environment (and the selection's own internal geometry) at full.
    sh.assign_nb_group(sh._atoms, 0)
    sh.soften_nb_selection(atoms, lambda_decouple)
    session.logger.info(
        'isolde decouple: %d atom(s) decoupled from their surroundings at lambda=%.3g '
        '(low=soft, 1=full). Transient -- cleared by "isolde decouple sel off" or when '
        'the simulation stops.' % (n_in, lambda_decouple))


def register_decouple_command(logger):
    from chimerax.core.commands import register, CmdDesc, BoolArg
    from chimerax.atomic import AtomsArg
    desc = CmdDesc(
        required=[('atoms', AtomsArg)],
        optional=[('state', BoolArg)],   # on|off|true|false; defaults to on (decouple)
        keyword=[('lambda_decouple', decouple_lambda_arg())],
        synopsis='Soften/restore a selection\'s nonbonded coupling to the rest of the '
                 'running simulation (per-group soft-core; transient, sim-scoped)'
    )
    register('isolde decouple', desc, decouple, logger=logger)
