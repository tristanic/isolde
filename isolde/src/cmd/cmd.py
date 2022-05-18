# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 25-Nov-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

from chimerax.core.errors import UserError

def block_if_sim_running(session):
    if hasattr(session, 'isolde') and session.isolde.simulation_running:
        raise UserError('This command is not available when a simulation is running!')

def get_singleton(session, create=True):
    if not session.ui.is_gui:
        return None
    from chimerax.core import tools
    from ..tool import ISOLDE_ToolUI
    return tools.get_singleton(session, ISOLDE_ToolUI, 'ISOLDE', create=create)


def isolde_start(session):
    ''' Start the ISOLDE GUI '''
    if not session.ui.is_gui:
        session.logger.warning("Sorry, ISOLDE currently requires ChimeraX to be in GUI mode")
        return

    get_singleton(session)
    return session.isolde

def isolde_set(session, time_steps_per_gui_update=None, temperature=None,
        gpu_device_index=None):
    isolde=isolde_start(session)
    sp = isolde.sim_params
    if time_steps_per_gui_update is not None:
        sp.sim_steps_per_gui_update=time_steps_per_gui_update
    if temperature is not None:
        sp.temperature=temperature
    if gpu_device_index is not None:
        sp.device_index=gpu_device_index

def isolde_select(session, model):
    isolde = isolde_start(session)
    isolde.selected_model = model

def isolde_sim(session, cmd, atoms=None, discard_to=None):
    '''
    Start, stop or pause an interactive simulation.

    Parameters
    ----------
    cmd: One of 'start', 'pause', 'resume', 'checkpoint', 'revert' or 'stop'
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
    valid_commands = ('start', 'pause', 'resume', 'checkpoint', 'revert', 'stop')
    log = session.logger
    if cmd not in valid_commands:
        raise UserError('Unrecognised command! Should be one of {}'.format(
            ', '.join(valid_commands)
        ))
    if cmd != 'start' and atoms is not None:
        log.warning('Atoms argument is not required for this command. Ignored.')
    isolde = isolde_start(session)

    if cmd == 'start':
        if atoms is None:
            model = isolde.selected_model
            if model is None:
                raise UserError('You must load a model before starting a simulation!')
            atoms = model.atoms

        # from chimerax.core.atomspec import AtomSpec
        # if isinstance(atoms, AtomSpec):
        #     objects = atoms.evaluate(session)


        us = atoms.unique_structures
        if len(us) == 0:
            raise UserError('No atoms selected!')
        if len(us) != 1:
            if isolde.selected_model is not None and isolde.selected_model in us:
                model = isolde.selected_model
            else:
                raise UserError('If selection is not from the current model active in ISOLDE, all atoms must be from the same model!')
        else:
            model = us[0]
        if model != isolde.selected_model:
            isolde.change_selected_model(model)
        session.selection.clear()
        atoms.selected = True
        isolde.start_sim()
        return

    if not isolde.simulation_running:
        log.warning('No simulation is currently running!')
        return

    if cmd == 'pause':
        isolde.pause()
    
    elif cmd == 'resume':
        isolde.resume()

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

def isolde_report(session, report=True, interval=20):
    isolde = isolde_start(session)
    if not isolde.simulation_running:
        raise UserError('This command is only valid when a simulation is running!')
    sm = isolde.sim_manager
    if report:
        sm.start_reporting_performance(interval)
    else:
        sm.stop_reporting_performance()

def isolde_ignore(session, residues=None, ignore=True):
    isolde_start(session)
    if session.isolde.simulation_running:
        session.logger.warning('Changes to the set of ignored residues will not'
            'be applied until you stop the current simulation.')
    if not residues and not ignore:
        # Stop ignoring all residues in the selected model
        m = session.isolde.selected_model
        if m is None:
            return
        session.isolde.ignored_residues = None
    if residues:
        for m in residues.unique_structures:
            mres = m.residues.intersect(residues)
            for r in mres:
                r.isolde_ignore = ignore
            session.logger.info('ISOLDE: currently ignoring {} residues in model {}'.format(
                len([r for r in m.residues if getattr(r, 'isolde_ignore', False)]), m.id_string
            ))

def incr_b_factor(session, b_add, atoms=None):
    B_MAX = 500
    if atoms is None:
        from chimerax.atomic import selected_atoms
        atoms = selected_atoms(session)
    if any(atoms.bfactors+b_add < 0):
        raise UserError('Applying this command would reduce the B-factor of at least one atom to below zero.')
    atoms.bfactors += b_add
    atoms[atoms.bfactors > B_MAX].bfactors = B_MAX 

def isolde_stop_ignoring(session, residues=None):
    isolde_ignore(session, residues, ignore=False)


def isolde_tutorial(session):
    from chimerax.help_viewer import show_url
    import pathlib
    import os
    root_dir = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(root_dir, '..', 'docs', 'user', 'tutorials', 'isolde.html')
    show_url(session, pathlib.Path(fname).as_uri())

_available_demos = {
    'crystal_intro': (
        'load_crystal_demo',
        'crystallographic demo: PDB ID 3io0',
        []),
    'cryo_em_intro': (
        'load_cryo_em_demo',
        'cryo-EM demo: PDB ID 6out, EMDB ID 20205',
        ('model_only',)),
}
def isolde_demo(session, demo_name = None, model_only=False, start_isolde=True):
    if start_isolde:
        isolde_start(session)
    def _load_demo(*_, session=session, demo_name=demo_name, model_only=model_only):        
        demo_info = _available_demos[demo_name]
        load_fn_name, description, kwarg_names = demo_info
        if model_only and 'model_only' not in kwarg_names:
            session.logger.warning("modelOnly argument is only applicable to cryo-EM data. Ignoring.")
        kwargs = {}
        if 'model_only' in kwarg_names:
            kwargs['model_only']=model_only
        from .. import isolde
        load_fn = getattr(isolde, load_fn_name)
        load_fn(session, **kwargs)
        session.logger.info("Loaded " + description)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER
    session.triggers.add_handler('new frame', _load_demo)

def isolde_step(session, residue=None, view_distance=None, interpolate_frames=None,
        polymeric_only=True, select=True):
    from chimerax.atomic import Residues, Residue
    if isinstance(residue, Residues):
        if len(residue) == 0:
            raise UserError('Selection contains no residues!')
        if len(residue) > 1:
            session.logger.warning('Multiple residues selected! Going to the first...')
        residue = residue[0]
        m = residue.structure
    else:
        if hasattr(session, 'isolde'):
            m = session.isolde.selected_model
        else:
            from chimerax.atomic import AtomicStructure
            atomic_structures = session.models.list(type=AtomicStructure)
            if len(atomic_structures):
                m = atomic_structures[0]
            else:
                m = None
        if m is None:
            raise UserError('No atomic structures are open!')
    from ..navigate import get_stepper
    rs = get_stepper(m)
    if view_distance is not None:
        rs.view_distance = view_distance
    if interpolate_frames is not None:
        rs.interpolate_frames = interpolate_frames
    if select:
        session.selection.clear()
    if residue is None:
        rs.incr_residue()
    elif isinstance(residue, Residue):
        rs.step_to(residue)
    else:
        # Assume residue is a string argument
        rarg = residue.lower()
        if rarg=='first':
            rs.first_residue(polymeric_only)
            rs.step_direction='next'
        elif rarg=='last':
            rs.last_residue(polymeric_only)
            rs.step_direction='previous'
        elif rarg in ('next', 'prev'):
            rs.incr_residue(rarg, polymeric_only)
        else:
            raise UserError('Unrecognised residue argument! If specified, must '
                'be either a residue, "first", "last", "next" or "prev"')
    if select:
        atoms = rs.current_residue.atoms
        atoms.selected=True
        atoms.intra_bonds.selected=True


def isolde_jump(session, direction="next"):
    isolde_start(session)
    m = session.isolde.selected_model
    from chimerax.core.errors import UserError
    if m is None:
        raise UserError("Please open a model first!")
    from ..navigate import get_stepper
    rs = get_stepper(m)
    rs.incr_chain(direction)

def register_isolde(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, PositiveIntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf, Or)
    from chimerax.atomic import AtomsArg, ResiduesArg

    def register_isolde_start():
        desc = CmdDesc(
            synopsis = 'Start the ISOLDE GUI'
        )
        register('isolde start', desc, isolde_start, logger=logger)

    def register_isolde_set():
        desc = CmdDesc(
            keyword = [
                ('time_steps_per_gui_update', IntArg),
                ('temperature', FloatArg),
                ('gpu_device_index', IntArg),
            ],
            synopsis = 'Adjust ISOLDE simulation settings',
        )
        register('isolde set', desc, isolde_set, logger=logger)

    def register_isolde_sim():
        desc = CmdDesc(
            optional=[('atoms', AtomsArg)],
            required=[('cmd', EnumOf(('start', 'pause', 'resume', 'checkpoint', 'revert', 'stop')))],
            keyword=[('discard_to', EnumOf(('start', 'checkpoint')))],
            synopsis='Start, stop or pause an interactive simulation'
            )
        register('isolde sim', desc, isolde_sim, logger=logger)

    def register_isolde_select():
        from .argspec import IsoldeStructureArg
        desc = CmdDesc(
            required=[('model', IsoldeStructureArg)],
            synopsis = "Set the specified model as ISOLDE's current selected model"
        )
        register('isolde select', desc, isolde_select, logger=logger)

    def register_isolde_report():
        desc = CmdDesc(
            optional=[('report', BoolArg),],
            keyword=[('interval', IntArg)],
            synopsis='Report the current simulation performance to the status bar'
        )
        register ('isolde report', desc, isolde_report, logger=logger)

    def register_isolde_ignore():
        desc = CmdDesc(
            optional=[('residues', ResiduesArg),],
            synopsis=('Tell ISOLDE to ignore a selection of residues during '
                'simulations. To reinstate them, use "isolde ~ignore"')
        )
        register('isolde ignore', desc, isolde_ignore, logger=logger)

    def register_isolde_stop_ignore():
        desc = CmdDesc(
            optional=[('residues', ResiduesArg),],
            synopsis=('Tell ISOLDE to stop ignoring a set of residues during '
                'simulations.')
        )
        register('isolde ~ignore', desc, isolde_stop_ignoring, logger=logger)


    def register_isolde_tutorial():
        desc = CmdDesc(
            synopsis='Load a help page with worked examples using ISOLDE'
        )
        register('isolde tutorial', desc, isolde_tutorial, logger=logger)

    def register_isolde_demo():
        desc = CmdDesc(
            synopsis='Load a small crystallographic or cryo-EM model for use in interactive tutorials',
            required=[('demo_name', EnumOf(list(_available_demos.keys())))],
            keyword=[('model_only', BoolArg),
                      ('start_isolde', BoolArg)]
        )
        register('isolde demo', desc, isolde_demo, logger=logger)

    def register_isolde_step():
        desc = CmdDesc(
            synopsis=('Step the view forward or backward through the chain, or '
                'jump to a specific residue'
            ),
            optional=[
                ('residue', Or(ResiduesArg, StringArg)),
            ],
            keyword=[
                ('view_distance', FloatArg),
                ('interpolate_frames', PositiveIntArg),
                ('polymeric_only', BoolArg),
                ('select', BoolArg)
            ]
        )
        register('isolde stepto', desc, isolde_step, logger=logger)
    
    def register_isolde_change_b():
        desc = CmdDesc(
            synopsis=('Change the B-factors of a group of atoms by the specified amount'),
            required=[
                ('b_add', FloatArg),
            ],
            optional=[
                ('atoms', AtomsArg),
            ],
        )
        register('isolde adjust bfactors', desc, incr_b_factor, logger=logger)

    def register_isolde_jump():
        desc = CmdDesc(
            synopsis=('Jump the view to the first residue of the next chain or '
                'the last residue of the previous chain.'),
            optional=[
                ('direction', EnumOf(('next', 'prev')))
            ]
        )
        register('isolde jumpto', desc, isolde_jump, logger=logger)

    def register_isolde_shorthand():
        desc = CmdDesc(
            synopsis=('Initialise two-character aliases for commonly-used ISOLDE commands')
        )
        from .shorthand import register_isolde_shorthand_commands
        register('isolde shorthand', desc, register_isolde_shorthand_commands, logger=logger)

    register_isolde_start()
    register_isolde_set()
    register_isolde_select()
    register_isolde_sim()
    register_isolde_report()
    register_isolde_ignore()
    register_isolde_stop_ignore()
    register_isolde_tutorial()
    register_isolde_demo()
    register_isolde_step()
    register_isolde_jump()
    register_isolde_change_b()
    register_isolde_shorthand()
    from chimerax.isolde.remote_control import register_remote_commands
    register_remote_commands(logger)
    from chimerax.isolde.restraints.cmd import register_isolde_restrain
    register_isolde_restrain(logger)
    from chimerax.isolde.manipulations.cmd import register_manip_commands
    register_manip_commands(logger)
    from chimerax.isolde.atomic.building.cmd import register_building_commands
    register_building_commands(logger)
    from chimerax.isolde.openmm.cmd import register_ff_cmd
    register_ff_cmd(logger)
    from chimerax.isolde.output.refinement.phenix.cmd import register_phenix_commands
    register_phenix_commands(logger)
    from chimerax.isolde.output.refinement.ccp4 import register_ccp4_commands
    register_ccp4_commands(logger)
    from chimerax.isolde.openmm.amberff.parameterise import register_isolde_param
    register_isolde_param(logger)
