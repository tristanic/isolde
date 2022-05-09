from chimerax.ui.gui import MainToolWindow

from chimerax.core.tools import ToolInstance
from Qt.QtCore import Qt

from chimerax.isolde.ui.ui_base import (
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem, QComboBox
)

from Qt.QtWidgets import (
    QGridLayout, QFrame, QLabel, QCheckBox, QPushButton, QMenu
)

from chimerax.isolde.ui.ui_base import QComboBox, DefaultHLayout, DefaultVLayout, DefaultSpacerItem


class Rama_ToolUI(ToolInstance):
    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        tw = self.tool_window=RamaMainWin(self)
        tw.manage(placement=None, allowed_areas=Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)


SELECTED_MODEL_CHANGED = 'selected model changed'
DISPLAY_MODE_CHANGED = 'display mode changed'


class RamaMainWin(MainToolWindow):
    from enum import IntFlag
    class Restrictions(IntFlag):
        ALL_PROTEIN = 0
        RESTRICTED = 1
        DISFAVORED_ONLY = 2
        OUTLIERS_ONLY = 4


    def __init__(self, tool_instance, **kw):
        super().__init__(tool_instance, **kw)
    
        from chimerax.isolde import session_extensions as sx
        mgr = self._rama_mgr = sx.get_ramachandran_mgr(self.session)
        cenum = self._case_enum = mgr.RamaCase

        from chimerax.core.triggerset import TriggerSet
        t = self.triggers = TriggerSet()
        for tname in (SELECTED_MODEL_CHANGED, DISPLAY_MODE_CHANGED):
            t.add_trigger(tname)

        self._display_mode = self.Restrictions.ALL_PROTEIN
        self._current_model = None
        self._handlers = []

        parent = self.ui_area
        main_layout = self.main_layout = DefaultVLayout()
        menu_layout = self.menu_layout = DefaultHLayout()
        main_layout.addLayout(menu_layout)
        parent.setLayout(main_layout)

        menu_layout.addWidget(QLabel('Current model: ', parent=parent))
        msb = self.model_select_button = QPushButton(parent)
        msm = self.model_select_menu = QMenu(msb)
        msm.aboutToShow.connect(self._populate_models_menu)
        msb.setText('Choose a model')
        msb.setMenu(msm)
        menu_layout.addWidget(msb)

        from chimerax.core.models import MODEL_ID_CHANGED, REMOVE_MODELS
        self._handlers.append(self.session.triggers.add_handler(MODEL_ID_CHANGED, self._session_models_changed_cb))
        self._handlers.append(self.session.triggers.add_handler(REMOVE_MODELS, self._session_models_changed_cb))
        self._handlers.append(self.session.triggers.add_handler('selection changed', self._selection_changed_cb))

        menu_layout.addWidget(QLabel('Showing: ', parent=parent))
        dmb = self.display_mode_menu_button = QPushButton('All', parent)
        dmm = self.display_mode_menu = QMenu(dmb)
        self._populate_display_modes_menu()
        dmb.setMenu(dmm)
        menu_layout.addWidget(dmb)

        menu_layout.addWidget(QLabel('Rama case(s): ', parent=parent))
        cmb = self.rama_case_menu_button = QPushButton('General', parent=parent)
        cmm = self.rama_case_menu = QMenu(cmb)
        cmb.setMenu(cmm)
        menu_layout.addWidget(cmb)



        menu_layout.addItem(DefaultSpacerItem())

        plot_layout = self.plot_layout = QGridLayout()

        self._plots = {}
        from .ramaplot_new import RamaPlot
        from math import floor
        for i, case in enumerate(reversed(cenum)):
            row = floor(i/3)
            column = i%3
            if case != cenum.NONE:
                p = self._plots[case] = RamaPlot(self.session, self, case)
                plot_layout.addWidget(p, row, column, 1, 1)
        main_layout.addLayout(plot_layout)
        
        self._populate_rama_cases_menu()
        self._show_one_plot(cenum.GENERAL)

        self._selection_changed_cb()


    def _show_one_plot(self, case):
        for c, p in self._plots.items():
            p.setVisible(c==case)
    
    def _show_all_plots(self):
        for c, p in self._plots.items():
            p.setVisible(True)



    def _populate_rama_cases_menu(self):
        from chimerax.isolde.molobject import RamaMgr
        case_dict = RamaMgr.RAMA_CASE_DETAILS
        menu = self.rama_case_menu
        for case in self._plots.keys():
            details = case_dict[case]
            def select_case(*_,c=case):
                self._show_one_plot(c)
            action = menu.addAction(details['name'])
            action.triggered.connect(select_case)
        a = menu.addAction('All')
        a.triggered.connect(self._show_all_plots)



    def _populate_models_menu(self):
        menu = self.model_select_menu
        menu.clear()
        from chimerax.atomic import AtomicStructure
        structures = [m for m in self.session.models.list() if type(m)==AtomicStructure]
        for s in structures:
            title = f'#{s.id_string}: {s.name}'
            action = menu.addAction(title)
            def set_structure(*_, structure=s):
                self.current_model = structure
            action.triggered.connect(set_structure)

    def _populate_display_modes_menu(self):
        labels = {
            'All protein': self.Restrictions.ALL_PROTEIN,
            'Custom selection': self.Restrictions.RESTRICTED,
            'Chain': self.Restrictions.RESTRICTED    
        }
        dmb = self.display_mode_menu_button
        dmm = self.display_mode_menu
        def show_all():
            label = 'All protein'
            dmb.setText(label)
            self.restrict_to = None
        a = dmm.addAction('All protein')
        a.triggered.connect(show_all)
        chain_menu = dmm.addMenu('Chain')
        def populate_chain_menu(*_):
            chains = []
            cm = self.current_model
            if cm is not None:
                from chimerax.atomic import Residue
                chains = cm.chains[cm.chains.polymer_types==Residue.PT_AMINO]
            chain_menu.clear()
            for c in chains:
                def restrict_to_chain(*_, chain=c):
                    self.restrict_to = cm.residues[cm.residues.chain_ids==c.chain_id]
                    dmb.setText(f'Chain {chain.chain_id}')
                ca = chain_menu.addAction(c.chain_id)
                ca.triggered.connect(restrict_to_chain)
        chain_menu.aboutToShow.connect(populate_chain_menu)
        def restrict_to_selection():
            from chimerax.atomic import selected_residues
            self.restrict_to = selected_residues(self.session)
            dmb.setText('Custom selection')
        a = self._restrict_to_selection_action = dmm.addAction('Current selection')
        a.triggered.connect(restrict_to_selection)
        
                


        

    @property
    def current_model(self):
        return self._current_model
    
    @current_model.setter
    def current_model(self, structure):
        self._current_model = structure
        if structure is None:
            self.model_select_button.setText('Choose a model')
        else:
            self.model_select_button.setText(f'#{structure.id_string}')
        self.triggers.activate_trigger(SELECTED_MODEL_CHANGED, structure)
        self.update()
    
    @property
    def restrict_to(self):
        return getattr(self, '_restrict_to', None)
    
    @restrict_to.setter
    def restrict_to(self, residues):
        if residues is not None:
            cm = self.current_model
            if cm is not None:
                residues = residues.intersect(cm.residues)
                from chimerax.atomic import Residue
                residues = residues[residues.polymer_types==Residue.PT_AMINO]
                if len(residues)==0:
                    from chimerax.core.errors import UserError
                    raise UserError('RamaPlot: No protein residues found in selection. Ignoring restrict_to call.')
                    return
        self._restrict_to = residues
        if residues is not None:
            self.display_mode |= self.Restrictions.RESTRICTED
        else:
            self.display_mode &= ~self.Restrictions.RESTRICTED

    @property
    def display_mode(self):
        return self._display_mode
    
    @display_mode.setter
    def display_mode(self, mode):
        if mode != self._display_mode:
            self._display_mode = mode
            self.triggers.activate_trigger(DISPLAY_MODE_CHANGED, mode)

    def update(self):
        pass


    def _session_models_changed_cb(self, trigger_name, *_):
        cm = self.current_model
        from chimerax.core.models import MODEL_ID_CHANGED, REMOVE_MODELS
        if trigger_name==REMOVE_MODELS:
            if cm is not None and cm not in self.session.models.list():
                self.current_model = None
        elif trigger_name==MODEL_ID_CHANGED:
            if cm is not None:
                # Just in case it was this model
                self.model_select_button.setText(f'#{cm.id_string}')

    def _selection_changed_cb(self, *_):
        cm = self.current_model
        sel_contains_protein = False
        if cm is not None:
            sel = cm.atoms[cm.atoms.selecteds].unique_residues
            from chimerax.atomic import Residue
            sel = sel[sel.polymer_types==Residue.PT_AMINO]
            if len(sel)>0:
                sel_contains_protein = True
        self._restrict_to_selection_action.setEnabled(sel_contains_protein)
            

    def cleanup(self):
        for h in self._handlers:
            h.remove()
        super().cleanup()



