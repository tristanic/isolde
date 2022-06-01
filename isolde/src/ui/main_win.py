from contextlib import contextmanager
from chimerax.ui.gui import MainToolWindow

from Qt.QtWidgets import (
    QFrame, QLabel,
    QPushButton, QTabWidget, QWidget,
)
from Qt import QtCore
from Qt.QtGui import QPixmap, QFont

from .ui_base import (
    DefaultVLayout, DefaultHLayout, DefaultSpacerItem, QComboBox,
    UI_Panel_Base
)

from .util import slot_disconnected

from . import icon_dir

class IsoldeMainWin(MainToolWindow):

    def __init__(self, tool_instance, **kw):
        super().__init__(tool_instance, **kw)
        if hasattr(self.session, 'isolde'):
            self.isolde = self.session.isolde
            self.isolde._gui = self
        else:
            from ..isolde import Isolde
            self.isolde = Isolde(self.session, gui=self)
        self.isolde.gui_mode = True
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
        from .ui_base import IsoldeTab
        from .general_tab import GeneralTab
        self.general_tab = GeneralTab(self.session, self.isolde, self, tabw, "General")
        from .validation_tab import ValidationTab
        self.validate_tab = ValidationTab(self.session, self.isolde, self, tabw, "Validate")
        from ..problem_regions.ui import ProblemAggregatorTab
        self.problems_tab = ProblemAggregatorTab(self.session, self.isolde, self, tabw)



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
        mmcb = self.master_model_combo_box = QComboBox(tw)
        mmcb.setMinimumSize(QtCore.QSize(150,0))
        mmcb.setToolTip('Select the model to work on in ISOLDE')
        l1.addWidget(mmcb)
        l1.addItem(DefaultSpacerItem(200))
        l1.addWidget(QLabel('Experience level: ', parent=tw))
        from .ui_base import ExpertModeSelector
        elcb = self.expert_mode_combo_box = ExpertModeSelector(tw)
        l1.addWidget(elcb)

        from chimerax.core.models import ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS
        for event_type in (ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS):
            self._session_trigger_handlers.append(
                session.triggers.add_handler(event_type, self._update_model_list_cb)
            )
        
        self._update_model_list_cb(None, None)
        self._isolde_trigger_handlers.append(session.isolde.triggers.add_handler(self.isolde.SIMULATION_STARTED, self._sim_start_cb))
        self._isolde_trigger_handlers.append(session.isolde.triggers.add_handler(self.isolde.SIMULATION_TERMINATED, self._sim_end_cb))

        mmcb.currentIndexChanged.connect(self._change_selected_model_cb)

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


    def _change_selected_model_cb(self, *_):
        mmcb = self.master_model_combo_box
        from ..util import block_managed_trigger_handler
        with block_managed_trigger_handler(self, '_selected_model_changed_handler'):
            m = mmcb.currentData()
            if m is not None:
                if mmcb.itemData(0) is None:
                    mmcb.removeItem(0)
                self.isolde.change_selected_model(m)

    def _sim_start_cb(self, *_):
        self.master_model_combo_box.setEnabled(False)
    
    def _sim_end_cb(self, *_):
        self.master_model_combo_box.setEnabled(True)

    @contextmanager
    def _block_update_model_list_cb(self):
        self._update_model_list_cb_blocked = True
        yield
        self._update_model_list_cb_blocked = False

    def _update_model_list_cb(self, _, models=None):
        from chimerax.atomic import AtomicStructure
        blocked = getattr(self, '_update_model_list_cb_blocked', False)
        if blocked:
            return
        mmcb = self.master_model_combo_box
        from chimerax.core.models import Model
        if isinstance(models, Model):
            models = [models]
        cm = self.isolde.selected_model
        if models is None or cm is None:
            structures_need_update = True
        else:
            structures_need_update = False
            for m in models:
                # Use explicit type equality rather than isInstance because some previews 
                # are AtomicStructure subclasses
                if type(m) == AtomicStructure:
                    structures_need_update = True
                    break
            
            if cm is not None and cm not in self.session.models.list():
                # Model has been deleted
                cm = None
                structures_need_update = True
        
        if structures_need_update:
            with slot_disconnected(mmcb.currentIndexChanged, self._change_selected_model_cb), self._block_update_model_list_cb():
                mmcb.clear()
                models = [m for m in self.session.models.list() if type(m) == AtomicStructure]
                models = sorted(models, key=lambda m: m.id)
                if cm is None:
                    # If any models are already initialised with Clipper, choose the first. Otherwise, leave it
                    # up to the user
                    from chimerax.clipper import get_symmetry_handler
                    handlers = [get_symmetry_handler(m) for m in models]
                    handlers = [h for h in handlers if h is not None]
                    mmcb.addItem('Choose a model...', None)

                for m in models:
                    id_str = f'{m.id_string}. {m.name}'
                    mmcb.addItem(id_str, m)
                if cm is not None:
                    mmcb.setCurrentIndex(mmcb.findData(cm))
                else:
                    if len(handlers):
                        cm = handlers[0].structure
                        mmcb.setCurrentIndex(mmcb.findData(cm))
                    else:
                        mmcb.setCurrentIndex(0)
                    self._change_selected_model_cb()

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
        if m is None:
            self._update_model_list_cb(None)
            return
        mmcb = self.master_model_combo_box
        with slot_disconnected(mmcb.currentIndexChanged, self._change_selected_model_cb):
            mmcb.setCurrentIndex(mmcb.findData(m))



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


        
            
    


    
    

















    
    

