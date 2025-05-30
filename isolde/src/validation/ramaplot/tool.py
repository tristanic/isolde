from chimerax.ui.gui import MainToolWindow

from chimerax.core.tools import ToolInstance
from Qt.QtCore import Qt

from chimerax.isolde.ui.ui_base import (
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem, QComboBox
)

from Qt.QtWidgets import (
    QGridLayout, QFrame, QLabel, QCheckBox, QPushButton, QMenu,
    QSizePolicy
)

from chimerax.isolde.ui.ui_base import QComboBox, DefaultHLayout, DefaultVLayout, DefaultSpacerItem


class Rama_ToolUI(ToolInstance):
    SESSION_ENDURING = True
    TOOL_CLOSED = 'tool closed'
    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        tw = self.tool_window=RamaMainWin(self)
        tw.manage(placement=None, allowed_areas=Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)

        from chimerax.core.triggerset import TriggerSet
        triggers = self.triggers=TriggerSet()
        triggers.add_trigger(self.TOOL_CLOSED)

    def delete(self):
        self.triggers.activate_trigger(self.TOOL_CLOSED, None)
        super().delete()

SELECTED_MODEL_CHANGED = 'selected model changed'
DISPLAY_MODE_CHANGED = 'display mode changed'


class RamaMainWin(MainToolWindow):

    from enum import IntFlag
    class Restrictions(IntFlag):
        ALL_PROTEIN = 0
        RESTRICTED = 1
        DISFAVORED_ONLY = 2
        OUTLIERS_ONLY = 4

    SINGLE_PLOT_SCATTER_SIZE = 10
    ALL_PLOT_SCATTER_SIZE = 6


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
        self.manage(placement=None)

        parent = self.ui_area
        self.base_scatter_size = self.SINGLE_PLOT_SCATTER_SIZE
        main_layout = self.main_layout = DefaultVLayout()
        menu_layout = self.menu_layout = DefaultHLayout()
        main_layout.addLayout(menu_layout)
        parent.setLayout(main_layout)

        menu_layout.addWidget(QLabel('Model: ', parent=parent))
        msb = self.model_select_button = QPushButton(parent)
        msm = self.model_select_menu = QMenu(msb)
        msm.aboutToShow.connect(self._populate_models_menu)
        msb.setText('Choose a model')
        msb.setMenu(msm)
        menu_layout.addWidget(msb)

        from chimerax.core.models import MODEL_ID_CHANGED, REMOVE_MODELS
        self._handlers.append(self.session.triggers.add_handler(MODEL_ID_CHANGED, self._session_models_changed_cb))
        self._handlers.append(self.session.triggers.add_handler(REMOVE_MODELS, self._session_models_changed_cb))
        
        from chimerax.core.selection import SELECTION_CHANGED
        self._handlers.append(self.session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed_cb))

        menu_layout.addWidget(QLabel('Display: ', parent=parent))
        dmb = self.display_mode_menu_button = QPushButton('All protein', parent)
        dmm = self.display_mode_menu = QMenu(dmb)
        self._populate_display_modes_menu()
        dmb.setMenu(dmm)
        menu_layout.addWidget(dmb)

        menu_layout.addWidget(QLabel('Case(s): ', parent=parent))
        cmb = self.rama_case_menu_button = QPushButton('General', parent=parent)
        cmm = self.rama_case_menu = QMenu(cmb)
        cmb.setMenu(cmm)
        menu_layout.addWidget(cmb)



        menu_layout.addStretch()

        plot_layout = self.plot_layout = QGridLayout()
        plot_layout.setContentsMargins(1,0,1,0)
        plot_layout.setSpacing(0)

        self._plots = {}
        from .ramaplot import RamaPlot
        from math import floor
        for i, case in enumerate(reversed(cenum)):
            row = floor(i/3)
            column = i%3
            if case != cenum.NONE:
                p = self._plots[case] = RamaPlot(self.session, self, case)
                plot_layout.addWidget(p, row, column, 1, 1)
        main_layout.addLayout(plot_layout)
        main_layout.addStretch()
        
        self._populate_rama_cases_menu()
        self._visible_plots=[]
        self._show_one_plot(cenum.GENERAL)

        self._selection_changed_cb()

    def add_callback(self, triggerset, trigger_name, callback):
        '''
        Add a handler to a `chimerax.core.triggerset.TriggerSet` instance. The handler
        will be automatically removed when the Ramachandran plot is closed. If you need
        to remove it earlier, use `remove_callback()`.
        '''
        self._handlers.append(triggerset.add_handler(trigger_name, callback))

    def _show_one_plot(self, case):
        self.base_scatter_size = self.SINGLE_PLOT_SCATTER_SIZE
        for c, p in self._plots.items():
            p.setVisible(c==case)
        case_dict = self._rama_mgr.RAMA_CASE_DETAILS
        details = case_dict[case]
        self.rama_case_menu_button.setText(details['short name'])
        self._visible_plots = [self._plots[case]]
        self.update_scatter()


    def _show_all_plots(self):
        self.base_scatter_size = self.ALL_PLOT_SCATTER_SIZE
        for c, p in self._plots.items():
            p.setVisible(True)
        self._visible_plots = list(self._plots.values())
        self.rama_case_menu_button.setText('All')
        self.update_scatter()

    def _populate_rama_cases_menu(self):
        case_dict = self._rama_mgr.RAMA_CASE_DETAILS
        cmb = self.rama_case_menu_button
        menu = self.rama_case_menu
        for case in self._plots.keys():
            details = case_dict[case]
            def select_case(*_,c=case):
                self._show_one_plot(c)
            action = menu.addAction(details['name'])
            action.triggered.connect(select_case)
        # Unfortunately the Qt window system makes it *really* hard to 
        # reliably adjust the window size to fit, particularly when 
        # redrawing is slow. Disabling the "all" option for now.
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
        a = dmm.addAction('All protein')
        a.triggered.connect(self.display_all_residues)
        chain_menu = dmm.addMenu('Chain')
        def populate_chain_menu(*_):
            chains = []
            cm = self.current_model
            if cm is not None:
                from chimerax.atomic import Residue
                chains = cm.chains[cm.chains.polymer_types==Residue.PT_AMINO]
            chain_menu.clear()
            if not len(chains):
                dummy = chain_menu.addAction('No protein chains found')
                dummy.setEnabled(False)
            for c in chains:
                def restrict_to_chain(*_, chain=c):
                    self.restrict_to_chain(chain.chain_id)
                ca = chain_menu.addAction(c.chain_id)
                ca.triggered.connect(restrict_to_chain)
        chain_menu.aboutToShow.connect(populate_chain_menu)
        


        def restrict_to_selection():
            from chimerax.atomic import selected_residues
            self.restrict_to_selection(selected_residues(self.session))
        a = self._restrict_to_selection_action = dmm.addAction('Current selection')
        a.triggered.connect(restrict_to_selection)
        
        options_menu = dmm.addMenu('Options')
        a1 = options_menu.addAction('Hide favoured')
        a1.setCheckable(True)
        a2 = options_menu.addAction('Outliers only')
        a2.setCheckable(True)
        a3 = self._show_markup_action = options_menu.addAction('Show markup on model')
        a3.setCheckable(True)
        def _hide_favored(flag):
            if flag:
                self.display_mode |= self.Restrictions.DISFAVORED_ONLY
            else:
                self.display_mode &= ~self.Restrictions.DISFAVORED_ONLY
                if a2.isChecked():
                    a2.toggle()
                    return
            self.update_scatter()
        a1.toggled.connect(_hide_favored)
        def _outliers_only(flag):
            if flag:
                self.display_mode |= self.Restrictions.OUTLIERS_ONLY
                if not a1.isChecked():
                    a1.toggle()
                    return
            else:
                self.display_mode &= ~self.Restrictions.OUTLIERS_ONLY
            self.update_scatter()
        a2.toggled.connect(_outliers_only)
        def _show_markup(flag):
            self.show_markup = flag
        a3.toggled.connect(_show_markup)


                

    @property
    def show_markup(self):
        if self.current_model is None:
            return False
        from chimerax.isolde.validation import RamaAnnotator
        for m in self.current_model.child_models():
            if isinstance(m, RamaAnnotator):
                return True
        return False
    
    @show_markup.setter
    def show_markup(self, show):
        if show != self.show_markup:
            self._show_or_hide_model_markup(show)
        from Qt.QtCore import QSignalBlocker
        with QSignalBlocker(self._show_markup_action):
            self._show_markup_action.setChecked(show)
        
    def _show_or_hide_model_markup(self, show):
        if self.current_model is not None:
            from chimerax.core.commands import run
            if show:
                cmd = 'rama'
            else:
                cmd = '~rama'
            run(self.session, f'{cmd} #{self.current_model.id_string}', log=False)


    @property
    def current_model(self):
        return self._current_model
    
    @current_model.setter
    def current_model(self, structure):
        if self._current_model == structure:
            return
        self._current_model = structure
        if structure is None:
            self.model_select_button.setText('Choose a model')
            residues = None
        else:
            from Qt.QtCore import QSignalBlocker
            with QSignalBlocker(self._show_markup_action):
                self._show_markup_action.setChecked(self.show_markup)
            self.model_select_button.setText(f'#{structure.id_string}')
            residues = structure.residues
            ch = getattr(self, '_model_changes_handler', None)
            if ch is not None:
                ch.remove()
            self._model_changes_handler = structure.triggers.add_handler('changes', self._model_changes_cb)
        for plot in self._plots.values():
            plot.set_target_residues(residues)
        self.update_scatter()
        self.triggers.activate_trigger(SELECTED_MODEL_CHANGED, structure)

    def _model_changes_cb(self, trigger_name, changes):
        for p in self._visible_plots:
            p._model_changed_cb(trigger_name, changes)

    def update_scatter(self):
        for p in self._visible_plots:
            p.update_scatter()

    def restrict_to_selection(self, residues, display_text = 'Custom selection'):
        dmb = self.display_mode_menu_button
        self.restrict_to = residues
        if residues is not None:
            dmb.setText(display_text)
    
    def restrict_to_chain(self, chain_id):
        from chimerax.core.errors import UserError
        m = self.current_model
        if m is None:
            return
        residues = m.residues[m.residues.chain_ids == chain_id]
        from chimerax.atomic import Residue
        residues = residues[residues.polymer_types==Residue.PT_AMINO]
        if not len(residues):
            raise UserError(f'Chain {chain_id} of model #{m.id_string} contains no protein!')
        dmb = self.display_mode_menu_button
        self.restrict_to = residues
        dmb.setText(f'Chain {chain_id}')
    
    def display_all_residues(self):
        self.restrict_to = None
        self.display_mode_menu_button.setText('All protein')


    @property
    def restrict_to(self):
        '''
        Get/set the list of residues to be shown on the Ramachandran plot. Note that setting 
        via this property is mostly for internal use and *does not* update the text on the 
        display menu button - instead you should use `restrict_to_selection()`, 
        `restrict_to_chain()` or `display_all_residues()`.
        '''
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
            for p in self._plots.values():
                p.set_target_residues(residues)
        else:
            self.display_mode &= ~self.Restrictions.RESTRICTED
            if self.current_model is not None:
                for p in self._plots.values():
                    p.set_target_residues(self.current_model.residues)
        self.update_scatter()

    @property
    def display_mode(self):
        return self._display_mode
    
    @display_mode.setter
    def display_mode(self, mode):
        if mode != self._display_mode:
            self._display_mode = mode
            self.triggers.activate_trigger(DISPLAY_MODE_CHANGED, mode)


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
        if cm is not None:
            sel = cm.residues[cm.residues.selected]
            ramas = self._rama_mgr.get_ramas(sel)
            ramas = ramas[ramas.valids]
            if not len(ramas):
                self._restrict_to_selection_action.setEnabled(False)
                self.update_scatter()
                return
            if len(self._visible_plots)==1:
                cenum = self._case_enum
                import numpy
                unique_cases = numpy.unique(ramas.cases)
                case = cenum.GENERAL
                if len(unique_cases) == 1:
                    case = cenum(unique_cases[0])
                if case != self._visible_plots[0].case:
                    self._show_one_plot(case)
            self.update_scatter()
            self._restrict_to_selection_action.setEnabled(True)


        sel_contains_protein = False
        if cm is not None:
            sel = cm.atoms[cm.atoms.selecteds].unique_residues
            from chimerax.atomic import Residue
            sel = sel[sel.polymer_types==Residue.PT_AMINO]
            if len(sel)>0:
                sel_contains_protein = True
                import numpy

        self.update_scatter()
        self._restrict_to_selection_action.setEnabled(sel_contains_protein)
            

    def cleanup(self):
        for h in self._handlers:
            h.remove()
        ch = getattr(self, '_model_changes_handler', None)
        if ch is not None:
            ch.remove()
        super().cleanup()



