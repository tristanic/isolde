
def start_pause_resume(session):
    if not hasattr(session, 'isolde'):
        from chimerax.core.errors import UserError
        raise UserError('Please start the ISOLDE GUI first!')
    if not session.isolde.simulation_running:
        from chimerax.core.commands import run
        run(session, 'isolde sim start sel')
    else:
        session.isolde.pause_sim_toggle()

def toolbar_command(session, name):
    from chimerax.core.commands import run
    if name == 'isolde start':
        run(session, 'isolde start')
    elif name == 'associate map':
        pass
    elif name == 'start-pause':
        start_pause_resume(session)
    elif name == 'checkpoint save':
        session.isolde.checkpoint()
    elif name == 'checkpoint revert':
        session.isolde.revert_to_checkpoint()
    elif name == 'stop-keep':
        run(session, 'isolde sim stop')
    elif name == 'stop-revert':
        run(session, 'isolde sim stop discardTo checkpoint')
    elif name == 'stop-discard':
        run(session, 'isolde sim stop discardTo start')
    elif name == 'flip peptide':
        run(session, 'isolde pepflip sel')
    elif name == 'flip cis-trans':
        run(session, 'isolde cisflip sel')

