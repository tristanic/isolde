
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

class ToolbarButtonMgr:
    # (tab, section, name, display_name)
    all_buttons = {
        'Start ISOLDE': ('ISOLDE', 'Main', 'isolde start', 'Start ISOLDE'),
        'flip peptide': ('ISOLDE', 'Peptide bond', 'flip peptide', 'Flip peptide'),
        'flip cis-trans': ('ISOLDE', 'Peptide bond', 'flip cis-trans', 'cis<->trans'),
    }

    enable_if_single_peptide_selected=('flip peptide', 'flip cis-trans')

    def __init__(self, session):
        self.session = session
        tb = session.toolbar
        from collections import defaultdict
        self._handlers = defaultdict(list)
        session.triggers.add_handler('new frame', self._set_button_starting_states)
        self._initialize_callbacks()

    def _set_button_starting_states(self, *_):
        tb = self.session.toolbar
        for key, (tab, section, name, display_name) in self.all_buttons.items():
            if key != 'Start ISOLDE':
                tb.set_enabled(False, tab, section, display_name)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def isolde_started(self):
        isolde = self.session.isolde
        tab, section, name, display_name = self.all_buttons['Start ISOLDE']
        self.session.toolbar.set_enabled(False, tab, section, display_name)
        isolde.triggers.add_handler('isolde closed', self._isolde_close_cb)

    
    def _initialize_callbacks(self):
        session = self.session
        st = session.triggers
        for key in self.enable_if_single_peptide_selected:
            tab, section, name, display_name = self.all_buttons[key]
            cb = self._create_cb_enable_if_peptide_selected(session, tab, section, display_name)
            st.add_handler('selection changed', cb)
            cb()



    def _create_cb_enable_if_peptide_selected(self, session, tab, section, display_name):
        def button_cb(*_,session=session, tab=tab, section=section, display_name=display_name):
            enable=False
            isolde = getattr(session, 'isolde', None)
            if isolde is not None:
                m = session.isolde.selected_model
                if m is not None:
                    sel_res = m.atoms[m.atoms.selecteds].unique_residues
                    if len(sel_res)==1:
                        from chimerax.atomic import Residue
                        r = sel_res[0]
                        if r.polymer_type==Residue.PT_AMINO:
                            enable=True
            session.toolbar.set_enabled(enable, tab, section, display_name)
        return button_cb
    
    def _isolde_close_cb(self, *_):
        for triggerset, handlers in zip(self._handlers.items()):
            for h in handlers:
                triggerset.remove_handler(h)
        self._set_button_starting_states()






