# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 13-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



def get_singleton(session, create=True):
    if not session.ui.is_gui:
        return None
    from chimerax.core import tools
    from .tool import ISOLDE_ToolUI
    return tools.get_singleton(session, ISOLDE_ToolUI, 'ISOLDE', create=create)


def isolde_start(session):
    ''' Start the ISOLDE GUI '''
    if not session.ui.is_gui:
        session.logger.warning("Sorry, ISOLDE currently requires ChimeraX to be in GUI mode")
        return

    get_singleton(session)
    return session.isolde

def isolde_sim(session, cmd, atoms=None, discard_to=None):
    '''
    Start, stop or pause an interactive simulation.

    Parameters
    ----------
    cmd: One of 'start', 'pause', 'checkpoint', 'revert' or 'stop'
        'pause' toggles between active and paused
        'checkpoint' saves a snapshot of the current simulation state
        'revert' reverts the simulation to the last saved checkpoint (or the
            simulation starting configuration if no checkpoint was saved)
        'stop' stops the simulation.
    atoms: An optional atomic selection, only used if the cmd is 'start'. All
        atoms must be from the same model. If no selection is given, the model
        currently selected in the ISOLDE panel will be simulated in its
        entirety.
    discard_to: Only applicable when cmd is 'stop'. One of 'start' or
        'checkpoint'. Discards all changes since the chosen state when stopping
        the simulation.
    '''
    valid_commands = ('start', 'pause', 'checkpoint', 'revert', 'stop')
    log = session.logger
    if cmd not in valid_commands:
        raise TypeError('Unrecognised command! Should be one of {}'.format(
            ', '.join(valid_commands)
        ))
    if cmd != 'start' and atoms is not None:
        log.warning('Atoms argument is not required for this command. Ignored.')
    isolde = isolde_start(session)

    if cmd == 'start':
        if atoms is None:
            model = isolde.selected_model
            if model is None:
                raise RuntimeError('You must load a model before starting a simulation!')
            atoms = model.atoms

        # from chimerax.core.atomspec import AtomSpec
        # if isinstance(atoms, AtomSpec):
        #     objects = atoms.evaluate(session)


        us = atoms.unique_structures
        if len(us) != 1:
            raise RuntimeError('All atoms must be from the same model!')
        model = us[0]
        isolde.change_selected_model(model)
        session.selection.clear()
        atoms.selected = True
        isolde.start_sim()
        return

    if not isolde.simulation_running:
        log.warning('No simulation is currently running!')
        return

    if cmd == 'pause':
        isolde.pause_sim_toggle()

    elif cmd == 'checkpoint':
        if not isolde.sim_running:
            log.warning('Checkpointing is only available while a simulation is running!')
        isolde.checkpoint()

    elif cmd == 'revert':
        isolde.revert_to_checkpoint()

    elif cmd == 'stop':
        if discard_to is None:
            isolde.commit_sim()
        else:
            isolde.discard_sim(revert_to=discard_to, warn=False)

def isolde_tutorial(session):
    from chimerax.help_viewer import show_url
    import pathlib
    import os
    root_dir = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(root_dir, 'docs', 'isolde', 'tutorials', 'index.html')
    show_url(session, pathlib.Path(fname).as_uri())


def register_isolde(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf)
    from chimerax.atomic import AtomsArg

    def register_isolde_start():
        desc = CmdDesc(
            synopsis = 'Start the ISOLDE GUI'
        )
        register('isolde start', desc, isolde_start, logger=logger)

    def register_isolde_sim():
        desc = CmdDesc(
            optional=[('atoms', AtomsArg)],
            required=[('cmd', EnumOf(('start', 'pause', 'checkpoint', 'revert', 'stop')))],
            keyword=[('discard_to', EnumOf(('start', 'checkpoint')))],
            synopsis='Start, stop or pause an interactive simulation'
            )
        register('isolde sim', desc, isolde_sim, logger=logger)

    def register_isolde_tutorial():
        desc = CmdDesc(
            synopsis='Load a help page with worked examples using ISOLDE'
        )
        register('isolde tutorial', desc, isolde_tutorial, logger=logger)

    register_isolde_start()
    register_isolde_sim()
    register_isolde_tutorial()
    from chimerax.isolde.remote_control import register_remote_commands
    register_remote_commands(logger)
    from chimerax.isolde.restraints.cmd import register_isolde_restrain
    register_isolde_restrain(logger)
    from chimerax.isolde.manipulations.cmd import register_manip_commands
    register_manip_commands(logger)
