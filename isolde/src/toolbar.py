
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
    elif name == 'start sim':
        run(session, 'isolde sim start sel')
    elif name == 'pause sim':
        run(session, 'isolde sim pause')
    elif name == 'resume sim':
        run(session, 'isolde sim resume')
    
    
    elif name == 'checkpoint save':
        session.isolde.checkpoint()
    elif name == 'checkpoint revert':
        session.isolde.revert_to_checkpoint()
    elif name == 'stop-keep':
        run(session, 'isolde sim stop')
    elif name == 'stop-discard':
        run(session, 'isolde sim stop discardTo start')
    elif name == 'flip peptide':
        run(session, 'isolde pepflip sel')
    elif name == 'flip cis-trans':
        run(session, 'isolde cisflip sel')
    elif 'rotamer' in name:
        _rota_command(session, name)
    elif name == 'spotlight':
        run(session, 'clipper spotlight')
    elif name == 'mask':
        if hasattr(session, 'isolde'):
            radius = session.isolde.params.map_mask_radius
            focus = session.isolde.params.center_on_sel_when_masking
        else:
            radius = 4.0
            focus = True
        run(session, f'clipper isolate sel mask {radius} focus {focus}')
    elif name == 'step n':
        run(session, 'isolde stepto prev')
    elif name == 'step sel':
        from chimerax.atomic import selected_residues, concise_residue_spec
        run(session, f'isolde stepto {concise_residue_spec(session, selected_residues(session))}')
    elif name == 'step c':
        run(session, 'isolde stepto next')
    elif name == 'tug atom':
        run(session, 'ui mousemode right "isolde tug atom"')
    elif name == 'tug residue':
        run(session, 'ui mousemode right "isolde tug residue"')
    elif name == 'tug selection':
        run(session, 'ui mousemode right "isolde tug selection"')

def _rota_command(session, cmd):
    from chimerax.atomic import selected_residues
    r = selected_residues(session)[0]
    from chimerax.isolde import session_extensions as sx
    rmgr = sx.get_rotamer_mgr(session)
    rrmgr = sx.get_rotamer_restraint_mgr(session.isolde.selected_model)
    rota = rmgr.get_rotamer(r)
    if cmd == 'next rotamer':
        rrmgr.next_preview(rota)
    elif cmd == 'commit rotamer':
        rrmgr.commit_preview(rota)
        if session.isolde.simulation_running:
            session.isolde.sim_handler.push_coords_to_sim()
    elif cmd == 'restrain rotamer':
        rrmgr.set_targets(rota)
        restr = rrmgr.get_restraint(rota)
        restr.enabled = True
        from .isolde import OPENMM_RADIAL_SPRING_UNIT
        restr.set_spring_constant(session.isolde.sim_params.rotamer_spring_constant.value_in_unit(OPENMM_RADIAL_SPRING_UNIT))
    elif cmd == 'release rotamer':
        restr = rrmgr.get_restraint(rota)
        restr.enabled=False
    session._isolde_tb._update_rotamer_buttons()



class ToolbarButtonMgr:
    # (tab, section, name, display_name)
    all_buttons = {
        'Start ISOLDE': ('ISOLDE', 'Main', 'isolde start', 'Start ISOLDE'),
        'Start simulation': ('ISOLDE', 'Control', 'start sim', 'Start simulation'),
        'Pause simulation': ('ISOLDE', 'Control', 'pause sim', 'Pause simulation'),
        'Resume simulation': ('ISOLDE', 'Control', 'resume sim', 'Resume simulation'),
        'Store checkpoint': ('ISOLDE', 'Control', 'checkpoint save', 'Store checkpoint'),
        'Revert to checkpoint': ('ISOLDE', 'Control', 'checkpoint revert', 'Revert to checkpoint'),
        'Stop (keep)': ('ISOLDE', 'Control', 'stop-keep', 'Stop (keep)'),
        'Stop (discard)': ('ISOLDE', 'Control', 'stop-discard', 'Stop (discard)'),

        'flip peptide': ('ISOLDE', 'Peptide bond', 'flip peptide', 'Flip peptide'),
        'flip cis-trans': ('ISOLDE', 'Peptide bond', 'flip cis-trans', 'Flip cis<->trans'),

        'rotamer preview': ('ISOLDE', 'Rotamer', 'next rotamer', 'Preview next'),
        'rotamer commit':   ('ISOLDE', 'Rotamer', 'commit rotamer', 'Set coords'),
        'rotamer restrain': ('ISOLDE', 'Rotamer', 'restrain rotamer', 'Restrain'),
        'rotamer release':  ('ISOLDE', 'Rotamer', 'release rotamer', 'Release'),

        'spotlight': ('ISOLDE', 'Map', 'spotlight', 'Spotlight mode'),
        'mask': ('ISOLDE', 'Map', 'mask', 'Mask to selection'),

        'step n': ('ISOLDE', 'Navigate', 'step n', 'Step back'),
        'step sel': ('ISOLDE', 'Navigate', 'step sel', 'Step from here'),
        'step c': ('ISOLDE', 'Navigate', 'step c', 'Step forward'),

        'tug atom': ('ISOLDE', 'Tugging mode', 'tug atom', 'Tug atom'),
        'tug residue': ('ISOLDE', 'Tugging mode', 'tug residue', 'Tug residue'),
        'tug selection': ('ISOLDE', 'Tugging mode', 'tug selecction', 'Tug selection'),
    }

    enable_if_single_peptide_selected=('flip peptide', 'flip cis-trans')

    def set_enabled(self, key, enabled):
        tab, section, name, display_name = self.all_buttons[key]
        self.session.toolbar.set_enabled(enabled, tab, section, display_name)

    def __init__(self, session):
        self.session = session
        tb = session.toolbar
        from collections import defaultdict
        self._handlers = defaultdict(list)
        self._last_selected_rotamer = None
        session.triggers.add_handler('new frame', self._set_button_starting_states)
        self._initialize_callbacks()

    def _set_button_starting_states(self, *_):
        try:
            for key, (tab, section, name, display_name) in self.all_buttons.items():
                if key == 'Start ISOLDE':
                    self.set_enabled(key, True)
                else:
                    self.set_enabled(key, False)
        except ValueError:
            # Temporary hack because in ChimeraX 1.4 the Toolbar manager doesn't update on bundle install,
            # leading to a traceback when ISOLDE is first installed via the GUI.
            from .dialog import choice_warning
            result = choice_warning('You will need to restart ChimeraX before using ISOLDE '
                'for the first time. Would you like to close ChimeraX now?')
            if result:
                self.session.ui.exit()

        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def isolde_started(self):
        isolde = self.session.isolde
        self.set_enabled('Start ISOLDE', False)
        isolde.triggers.add_handler(isolde.ISOLDE_CLOSED, self._isolde_close_cb)
        isolde.triggers.add_handler(isolde.SIMULATION_STARTED, self._sim_start_cb)
        isolde.triggers.add_handler(isolde.SIMULATION_PAUSED, self._sim_pause_cb)
        isolde.triggers.add_handler(isolde.SIMULATION_RESUMED, self._sim_resume_cb)
        isolde.triggers.add_handler(isolde.SIMULATION_TERMINATED, self._sim_end_cb)
        isolde.triggers.add_handler(isolde.SELECTED_MODEL_CHANGED, self._update_map_buttons)
        isolde.triggers.add_handler(isolde.SELECTED_MODEL_CHANGED, self._update_stepper_buttons)

    def _sim_start_cb(self, *_):
        self.set_enabled('Start simulation', False)
        self.set_enabled('Pause simulation', True)
        self.set_enabled('Resume simulation', False)
        self.set_enabled('Store checkpoint', True)
        self.set_enabled('Revert to checkpoint', True)
        self.set_enabled('Stop (keep)', True)
        self.set_enabled('Stop (discard)', True)
        self.set_enabled('tug atom', True)
        self.set_enabled('tug residue', True)
        self.set_enabled('tug selection', True)
        self._update_map_buttons()
        bd = self.all_buttons['Pause simulation']
        def _cb(*_):
            self.session.toolbar.show_group_button(bd[0],bd[1],bd[3])
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('new frame', _cb)


    def _sim_pause_cb(self, *_):
        self.set_enabled('Pause simulation', False)
        self.set_enabled('Resume simulation', True)
        bd = self.all_buttons['Resume simulation']
        def _cb(*_):
            self.session.toolbar.show_group_button(bd[0],bd[1],bd[3])
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('new frame', _cb)
    
    def _sim_resume_cb(self, *_):
        self.set_enabled('Resume simulation', False)
        self.set_enabled('Pause simulation', True)
        bd = self.all_buttons['Pause simulation']
        def _cb(*_):
            self.session.toolbar.show_group_button(bd[0],bd[1],bd[3])
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('new frame', _cb)

    def _sim_end_cb(self, *_):
        self.set_enabled('Pause simulation', False)
        self.set_enabled('Resume simulation', False)
        self._enable_sim_start_button_if_necessary()
        self.set_enabled('Store checkpoint', False)
        self.set_enabled('Revert to checkpoint', False)
        self.set_enabled('Stop (keep)', False)
        self.set_enabled('Stop (discard)', False)
        self.set_enabled('tug atom', False)
        self.set_enabled('tug residue', False)
        self.set_enabled('tug selection', False)
        self._update_map_buttons()
        bd = self.all_buttons['Start simulation']
        def _cb(*_):
            self.session.toolbar.show_group_button(bd[0],bd[1],bd[3])
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('new frame', _cb)


    def _update_map_buttons(self, *_):
        enable_spotlight = False
        enable_mask = False
        if _selected_model_has_maps(self.session):
            isolde = self.session.isolde
            if not isolde.simulation_running:
                enable_spotlight = True
                if isolde.selected_model.residues.selected.sum() > 0:
                    enable_mask = True
        self.set_enabled('spotlight', enable_spotlight)
        self.set_enabled('mask', enable_mask)

    def _update_stepper_buttons(self, *_):
        enable=False
        if hasattr(self.session, 'isolde') and self.session.isolde.selected_model is not None:
            enable=True
        for button in ('step n', 'step c'):
            self.set_enabled(button, enable)
        if enable:
            sel = self.session.isolde.selected_atoms
            if sel is not None and len(sel.unique_residues)==1:
                self.set_enabled('step sel', True)
            else:
                self.set_enabled('step sel', False)
        else:
            self.set_enabled('step sel', False)


    def _update_rotamer_buttons(self, *_):
        enableds = {
            'rotamer preview': False,
            'rotamer commit':  False,
            'rotamer restrain': False,
            'rotamer release': False
        }
        session = self.session
        if hasattr(session, 'isolde') and session.isolde.selected_model is not None:
            m = session.isolde.selected_model
            if m is not None and not m.deleted:
                from chimerax.isolde import session_extensions as sx
                rmgr = sx.get_rotamer_mgr(session)
                rrmgr = sx.get_rotamer_restraint_mgr(m)
                sel_res = m.residues[m.residues.selected]
                if len(sel_res)==1:
                    r = sel_res[0]
                    rota = rmgr.get_rotamer(r)
                    if rota != self._last_selected_rotamer:
                        rrmgr.remove_preview()
                        self._last_selected_rotamer = rota
                    if rota is not None:
                        enableds['rotamer preview'] = True
                        rr = rrmgr.get_restraint(rota)
                        if rr is not None:
                            if rr.enabled:
                                enableds['rotamer release']=True
                        if rrmgr.showing_preview():
                            enableds['rotamer commit'] = True
                            enableds['rotamer restrain'] = True
                else:
                    rrmgr.remove_preview()
                    self._last_selected_rotamer = None
        for key, flag in enableds.items():
            self.set_enabled(key, flag)
                        


            

    def _initialize_callbacks(self):
        session = self.session
        st = session.triggers
        from chimerax.core.selection import SELECTION_CHANGED
        st.add_handler(SELECTION_CHANGED, self._selection_changed_cb)

    def _enable_if_single_peptide_selected_cb(self):
        enable=False
        isolde = getattr(self.session, 'isolde', None)
        if isolde is not None:
            m = isolde.selected_model
            if m is not None and not m.deleted:
                sel_res = m.residues[m.residues.selected]
                if len(sel_res)==1:
                    from chimerax.atomic import Residue
                    r = sel_res[0]
                    if r.polymer_type==Residue.PT_AMINO:
                        enable=True
        for key in self.enable_if_single_peptide_selected:
            self.set_enabled(key, enable)

    def _enable_sim_start_button_if_necessary(self):
        enable = False
        isolde = getattr(self.session, 'isolde', None)
        if isolde is not None and not isolde.simulation_running:
            m = isolde.selected_model
            if m is not None and not m.deleted:
                import numpy
                if numpy.any(m.atoms.selected):
                    enable = True
        self.set_enabled('Start simulation', enable)

    def _selection_changed_cb(self, *_):
        self._enable_if_single_peptide_selected_cb()
        self._enable_sim_start_button_if_necessary()
        self._update_rotamer_buttons()
        self._update_map_buttons()
        self._update_stepper_buttons()
    
    def _isolde_close_cb(self, *_):
        for triggerset, handlers in zip(self._handlers.items()):
            for h in handlers:
                triggerset.remove_handler(h)
        self._set_button_starting_states()


def _selected_model_has_maps(session):
    isolde = getattr(session, 'isolde', None)
    if isolde is None:
        return False
    m = session.isolde.selected_model
    if m is None:
        return False
    from chimerax.clipper import get_map_mgr
    mgr = get_map_mgr(m)
    if mgr is None:
        return False
    return (len(mgr.all_maps)>0)



