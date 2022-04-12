
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
    elif name == 'stop-revert':
        run(session, 'isolde sim stop discardTo checkpoint')
    elif name == 'stop-discard':
        run(session, 'isolde sim stop discardTo start')
    elif name == 'flip peptide':
        run(session, 'isolde pepflip sel')
    elif name == 'flip cis-trans':
        run(session, 'isolde cisflip sel')

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
        'Stop (revert)': ('ISOLDE', 'Control', 'stop-revert', 'Stop (revert)'),
        'Stop (discard)': ('ISOLDE', 'Control', 'stop-discard', 'Stop (discard)'),

        'flip peptide': ('ISOLDE', 'Peptide bond', 'flip peptide', 'Flip peptide'),
        'flip cis-trans': ('ISOLDE', 'Peptide bond', 'flip cis-trans', 'Flip cis<->trans'),

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
        session.triggers.add_handler('new frame', self._set_button_starting_states)
        self._initialize_callbacks()

    def _set_button_starting_states(self, *_):
        for key, (tab, section, name, display_name) in self.all_buttons.items():
            if key != 'Start ISOLDE':
                self.set_enabled(key, False)
            else:
                self.set_enabled(key, True)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def isolde_started(self):
        isolde = self.session.isolde
        self.set_enabled('Start ISOLDE', False)
        isolde.triggers.add_handler('isolde closed', self._isolde_close_cb)
        isolde.triggers.add_handler('simulation started', self._sim_start_cb)
        isolde.triggers.add_handler('simulation paused', self._sim_pause_cb)
        isolde.triggers.add_handler('simulation resumed', self._sim_resume_cb)
        isolde.triggers.add_handler('simulation terminated', self._sim_end_cb)

    def _sim_start_cb(self, *_):
        self.set_enabled('Start simulation', False)
        self.set_enabled('Pause simulation', True)
        self.set_enabled('Resume simulation', False)
        self.set_enabled('Store checkpoint', True)
        self.set_enabled('Revert to checkpoint', True)
        self.set_enabled('Stop (keep)', True)
        self.set_enabled('Stop (revert)', True)
        self.set_enabled('Stop (discard)', True)


    def _sim_pause_cb(self, *_):
        self.set_enabled('Pause simulation', False)
        self.set_enabled('Resume simulation', True)
    
    def _sim_resume_cb(self, *_):
        self.set_enabled('Resume simulation', False)
        self.set_enabled('Pause simulation', True)

    def _sim_end_cb(self, *_):
        self.set_enabled('Pause simulation', False)
        self.set_enabled('Resume simulation', False)
        self._enable_sim_start_button_if_necessary()
        self.set_enabled('Store checkpoint', False)
        self.set_enabled('Revert to checkpoint', False)
        self.set_enabled('Stop (keep)', False)
        self.set_enabled('Stop (revert)', False)
        self.set_enabled('Stop (discard)', False)


    def _initialize_callbacks(self):
        session = self.session
        st = session.triggers
        st.add_handler('selection changed', self._selection_changed_cb)

    def _enable_if_single_peptide_selected_cb(self):
        enable=False
        isolde = getattr(self.session, 'isolde', None)
        if isolde is not None:
            m = isolde.selected_model
            if m is not None and not m.was_deleted:
                sel_res = m.atoms[m.atoms.selecteds].unique_residues
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
            if m is not None and not m.was_deleted:
                import numpy
                if numpy.any(m.atoms.selected):
                    enable = True
        self.set_enabled('Start simulation', enable)

    def _selection_changed_cb(self, *_):
        self._enable_if_single_peptide_selected_cb()
        self._enable_sim_start_button_if_necessary()
    
    def _isolde_close_cb(self, *_):
        for triggerset, handlers in zip(self._handlers.items()):
            for h in handlers:
                triggerset.remove_handler(h)
        self._set_button_starting_states()






