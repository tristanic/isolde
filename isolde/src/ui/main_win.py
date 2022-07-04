from contextlib import contextmanager
from chimerax.ui.gui import MainToolWindow

from Qt.QtWidgets import (
    QFrame, QLabel,
    QPushButton, QMenu, QTabWidget, QWidget, 
    QToolBar
)
from Qt import QtCore
from Qt.QtGui import QPixmap, QFont, QIcon

from .ui_base import (
    DefaultVLayout, DefaultHLayout, DefaultSpacerItem, QComboBox,
    UI_Panel_Base
)
from Qt.QtCore import Qt

from .util import slot_disconnected

from . import icon_dir, DEFAULT_ICON_SIZE

class IsoldeMainWin(MainToolWindow):

    def __init__(self, tool_instance, **kw):
        super().__init__(tool_instance, **kw)
        if hasattr(self.session, 'isolde'):
            self.isolde = self.session.isolde
            self.isolde.gui = self
        else:
            from ..isolde import Isolde
            self.isolde = Isolde(self.session, gui=self)
        self._gui_panels = []
        sth = self._session_trigger_handlers = []
        ith = self._isolde_trigger_handlers = []

        smch = self._selected_model_changed_handler = self.isolde.triggers.add_handler(
            'selected model changed', self._selected_model_changed_cb
        )
        self._isolde_trigger_handlers.append(smch)

        parent = self.ui_area
        main_layout = DefaultVLayout()
        parent.setLayout(main_layout)

        self._prepare_top_frame(main_layout)

        tabw = self.main_tab_widget = QTabWidget(parent)
        tabw.setMinimumHeight(350)
        tabw.setElideMode(QtCore.Qt.ElideNone)
        tabw.setUsesScrollButtons(True)
        tabw.setDocumentMode(False)
        main_layout.addWidget(tabw)
        from .general_tab import GeneralTab
        self.general_tab = GeneralTab(self.session, self.isolde, self, tabw)
        from .restraints_tab import RestraintsTab
        self.restraints_tab = RestraintsTab(self.session, self.isolde, self, tabw)
        from .validation_tab import ValidationTab
        self.validate_tab = ValidationTab(self.session, self.isolde, self, tabw)
        from ..problem_regions.ui import ProblemAggregatorTab
        self.problems_tab = ProblemAggregatorTab(self.session, self.isolde, self, tabw)

        selection_frame = QWidget()
        main_layout.addWidget(selection_frame)
        self.selection_toolbar = SelectionPanel(self.session, self.isolde, self, selection_frame)

    def register_panel(self, panel):
        '''
        Register a GUI subpanel to be managed by this one. At a minimum, the panel *must* implement
        a `cleanup()` method to remove its callbacks when the GUI is destroyed.
        '''
        self._gui_panels.append(panel)

    def _prepare_top_frame(self, main_layout):
        import os
        session = self.session
        tf = self._top_frame = QFrame()
        layout = DefaultHLayout()
        tf.setLayout(layout)
        icon = os.path.join(icon_dir, 'logo_small.png')
        il = QLabel(tf)
        il.setFixedSize(QtCore.QSize(50, 65))
        il.setPixmap(QPixmap(icon))
        il.setScaledContents(True)
        layout.addWidget(il)

        inner = QFrame()
        li = DefaultVLayout()
        inner.setLayout(li)
        tw = QWidget()
        l1 = DefaultHLayout()
        tw.setLayout(l1)

        wol = QLabel(tw)
        wol.setText('Working on: ')
        l1.addWidget(wol)
        mmb = self.master_model_menu_button = QPushButton(tw)
        self._update_model_menu_button(self.isolde.selected_model)
        mmm = self.master_model_menu = QMenu()
        mmm.aboutToShow.connect(self._populate_available_models_menu)
        mmb.setMenu(mmm)
        mmb.setMinimumSize(QtCore.QSize(150,0))
        mmb.setToolTip('Select the model to work on in ISOLDE')
        l1.addWidget(mmb)
        l1.addItem(DefaultSpacerItem(200))
        l1.addWidget(QLabel('Experience level: ', parent=tw))
        from .ui_base import ExpertModeSelector
        elcb = self.expert_mode_combo_box = ExpertModeSelector(tw)
        l1.addWidget(elcb)

        self._isolde_trigger_handlers.append(session.isolde.triggers.add_handler(self.isolde.SIMULATION_STARTED, self._sim_start_cb))
        self._isolde_trigger_handlers.append(session.isolde.triggers.add_handler(self.isolde.SIMULATION_TERMINATED, self._sim_end_cb))

        li.addWidget(tw)

        bw = QWidget()
        l2 = DefaultHLayout()
        bw.setLayout(l2)

        status_frame = QWidget(tw)
        self.sim_status_indicator = SimStatusIndicator(session, self.isolde, self, status_frame)
        l2.addWidget(status_frame)

        l2.addItem(DefaultSpacerItem())

        tcb = self.tutorials_button = QPushButton('Tutorials', parent=bw)
        tcb.setToolTip("Start an interactive tutorial")
        from chimerax.core.commands import run
        tcb.clicked.connect(lambda *_:run(session, 'isolde tut'))
        l2.addWidget(tcb)

        hb = self.main_help_button = QPushButton(bw)
        hb.setText('Help')
        hb.clicked.connect(self.isolde.show_master_help_in_browser)
        l2.addWidget(hb)
        
        li.addWidget(bw)
        layout.addWidget(inner)

        main_layout.addWidget(tf)


    def _sim_start_cb(self, *_):
        self.master_model_menu_button.setEnabled(False)
    
    def _sim_end_cb(self, *_):
        self.master_model_menu_button.setEnabled(True)

    def _populate_available_models_menu(self):
        from chimerax.atomic import AtomicStructure
        mmm = self.master_model_menu
        mmm.clear()
        models = [m for m in self.session.models if type(m)==AtomicStructure]
        models = sorted(models, key=lambda m: m.id)
        for m in models:
            a = mmm.addAction(f'{m.id_string}: {m.name}')
            def _set_selected_model(_, model=m):
                self.isolde.selected_model = model
            a.triggered.connect(_set_selected_model)

    def cleanup(self):
        for h in self._session_trigger_handlers:
            h.remove()
        for h in self._isolde_trigger_handlers:
            h.remove()
        for panel in self._gui_panels:
            panel.cleanup()

    ###
    # ISOLDE event callbacks
    ###

    def _selected_model_changed_cb(self, trigger_name, m):
        self._update_model_menu_button(m)
    
    def _update_model_menu_button(self, m=None):
        mmb = self.master_model_menu_button
        text = m.id_string if m is not None else "None"
        if m is None:
            m = self.isolde.find_next_available_model()
            if m is not None:
                self.isolde.selected_model = m
                return
            mmb.setText('Choose a model')
        else:    
            mmb.setText(f'{m.id_string}')


class SimStatusIndicator(UI_Panel_Base):
    def __init__(self, session, isolde, gui, main_frame):
        super().__init__(session, isolde, gui, main_frame, sim_sensitive=True)
        font = QFont()
        font.setFamily('Carlito')
        font.setPointSize(14)
        font.setBold(True)
        layout = DefaultVLayout()
        self.main_frame.setLayout(layout)
        sl = self.status_label = QLabel(self.main_frame)
        sl.setFont(font)
        layout.addWidget(sl)
        ith = self._isolde_trigger_handlers
        ith.append(isolde.triggers.add_handler(
            'simulation resumed', self.sim_resume_cb
        ))
        ith.append(isolde.triggers.add_handler(
            'simulation paused', self.sim_pause_cb
        ))
    
    def sim_start_cb(self, *_):
        self.status_label.setText("<font color='red'>SIMULATION RUNNING</font>")
    
    def sim_pause_cb(self, *_):
        self.status_label.setText("<font color='blue'>SIMULATION PAUSED</font>")
    
    def sim_resume_cb(self, *_):
        return self.sim_start_cb()
    
    def sim_end_cb(self, *_):
        self.status_label.setText('')


class SelectionPanel(UI_Panel_Base):
    def __init__(self, session, isolde, gui, main_frame):
        import os
        super().__init__(session, isolde, gui, main_frame, sim_sensitive=False)
        ml = self.main_layout = DefaultVLayout()
        main_frame.setLayout(ml)
        ml.addWidget(QLabel('Selection Tools'))
        tb = self.toolbar = QToolBar()
        ml.addWidget(tb)
        tb.setMovable(False)
        tb.setIconSize(DEFAULT_ICON_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        tb.setStyleSheet('background: rgb(80,80,80); border:none; color: white;')
        en = self.extend_backward_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'extend_selection_n.png')),
            'Extend\nbackward')
        en.triggered.connect(self._extend_backward_along_chain)
        en.setToolTip('<span>Extend each contiguous protein or nucleic acid selection by one residue towards the start of the chain, stopping at chain breaks.</span>')
        ec = self.extend_forward_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'extend_selection_c.png')),
            'Extend\nforward')
        ec.triggered.connect(self._extend_forward_along_chain)
        ec.setToolTip('<span>Extend each contiguous protein or nucleic acid selection by one residue towards the end of the chain, stopping at chain breaks.</span>')
        sn = self.shrink_from_start_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'shrink_selection_n.png')),
            'Shrink\nfrom start')
        sn.triggered.connect(self._shrink_sel_from_start)
        sn.setToolTip('<span>Shrink each contiguous protein or nucleic acid selection by one residue away from the start of the chain.</span>')
        sc = self.shrink_from_end_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'shrink_selection_c.png')),
            'Shrink\nfrom end')
        sc.triggered.connect(self._shrink_sel_from_end)
        sc.setToolTip('<span>Shrink each contiguous protein or nucleic acid selection by one residue away from the end of the chain.</span>')
        tb.addSeparator()
        av = self.all_visible_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'visible_selected.png')),
            'All\nvisible'
        )
        av.triggered.connect(self.select_all_visible)
        av.setToolTip('<span>Select all residues with atoms currently visible in the working model.</span>')

        im = self.in_map_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'mask_selected.png')),
            'Within\nmap mask'
        )
        im.triggered.connect(self.select_all_in_map)
        im.setToolTip('<span>Select all residues from the working model with all heavy atoms inside the current map mask.</span>')
    
    def select_all_visible(self, *_):
        m = self.isolde.selected_model
        if m is None:
            return
        m.atoms.selected=False
        m.bonds.selected=False
        visibles = m.atoms[m.atoms.visibles]
        visibles.selected=True
        visibles.intra_bonds.selected=True

    def select_all_in_map(self, *_):
        import numpy
        m = self.isolde.selected_model
        if m is None:
            return
        m.atoms.selected = False
        m.bonds.selected = False
        from chimerax.clipper import get_map_mgr
        mgr = get_map_mgr(m)
        if mgr is None:
            return
        mask = mgr.zone_mgr.mask
        in_map = m.atoms[mask.interpolated_values(m.atoms.coords)>0]
        in_map.selected=True
        from chimerax.atomic import Residues
        all_in = []
        for r in in_map.unique_residues:
            if numpy.all(r.atoms[r.atoms.element_names!='H'].selected):
                all_in.append(r)
        in_map.selected = False
        atoms = Residues(all_in).atoms
        atoms.selected=True
        atoms.intra_bonds.selected=True
        

    def current_selection(self):
        m = self.isolde.selected_model
        if m is None:
            return None
        return m.residues[m.residues.selected]
    
    def _extend_forward_along_chain(self, *_):
        sel = self.current_selection()
        if sel is None or not len(sel):
            return
        from chimerax.isolde.selections import extend_selection_along_chains
        extend_selection_along_chains(sel, 1)
    
    def _extend_backward_along_chain(self, *_):
        sel = self.current_selection()
        if sel is None or not len(sel):
            return
        from chimerax.isolde.selections import extend_selection_along_chains
        extend_selection_along_chains(sel, -1)
    
    def _shrink_sel_from_end(self, *_):
        sel = self.current_selection()
        if sel is None or not len(sel):
            return
        from chimerax.isolde.selections import shrink_selection_by_one_res
        shrink_selection_by_one_res(sel, 1)

    def _shrink_sel_from_start(self, *_):
        sel = self.current_selection()
        if sel is None or not len(sel):
            return
        from chimerax.isolde.selections import shrink_selection_by_one_res
        shrink_selection_by_one_res(sel, -1)







            
    


    
    

















    
    

