from contextlib import contextmanager
from chimerax.ui.gui import MainToolWindow

import sys
if 'win' in sys.platform.lower():
    from .util import WinAutoResizeQComboBox as QComboBox
else:
    from Qt.QtWidgets import QComboBox

from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFrame, QLayout, 
    QLabel, QSizePolicy, QSpacerItem,
    QPushButton, QTabWidget, QWidget, QScrollArea
)
from Qt import QtCore
from Qt.QtGui import QPixmap

from .ui_base import DefaultVLayout, DefaultHLayout, DefaultSpacerItem

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

        self._session_trigger_handlers = []
        parent = self.ui_area

        main_layout = DefaultVLayout()
        parent.setLayout(main_layout)

        self._prepare_top_frame(main_layout)

        tabw = self.main_tab_widget = QTabWidget(parent)
        tabw.setElideMode(QtCore.Qt.ElideNone)
        tabw.setUsesScrollButtons(True)
        tabw.setDocumentMode(False)
        main_layout.addWidget(tabw)
        from .ui_base import IsoldeTab
        from .general_tab import GeneralTab
        self.general_tab = GeneralTab(self.session, self.isolde, self, tabw, "General")
        self.validate_tab = IsoldeTab(self, tabw, "Validate")
        self.problems_tab = IsoldeTab(self, tabw, "Problem Zones")



    def register_panel(self, panel):
        '''
        Register a GUI subpanel to be managed by this one. At a minimum, the panel *must* implement
        a `cleanup()` method to remove its callbacks when the GUI is destroyed.
        '''
        self._gui_panels.append(panel)

    def _prepare_top_frame(self, main_layout):
        import os
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

        from chimerax.core.models import ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS
        for event_type in (ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS):
            self._session_trigger_handlers.append(
                self.session.triggers.add_handler(event_type, self._update_model_list_cb)
            )
        self._update_model_list_cb(None, None)
        mmcb.currentIndexChanged.connect(self._change_selected_model_cb)

        li.addWidget(tw)

        bw = QWidget()
        l2 = DefaultHLayout()
        bw.setLayout(l2)


        l2.addItem(DefaultSpacerItem())

        tcb = self.tutorials_button = QPushButton('Tutorials', parent=bw)
        tcb.setToolTip("Start an interactive tutorial")
        from chimerax.core.commands import run
        tcb.clicked.connect(lambda *_:run(self.session, 'isolde tut'))
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
        with slot_disconnected(mmcb.currentIndexChanged, self._change_selected_model_cb):
            m = mmcb.currentData()
            self.isolde.change_selected_model(m)
    
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
        if models is None:
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

                for m in models:
                    id_str = f'{m.id_string}. {m.name}'
                    mmcb.addItem(id_str, m)
                if cm is not None:
                    mmcb.setCurrentIndex(mmcb.findData(cm))
                else:
                    mmcb.setCurrentIndex(0)
                    self._change_selected_model_cb()

    def cleanup(self):
        for h in self._session_trigger_handlers:
            h.remove()
        for panel in self._gui_panels:
            panel.cleanup()


class ExpertModeSelector(QComboBox):
    _expert_modes = ('default', 'advanced', 'developer')
    expertModeChanged = QtCore.Signal()
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for mode in self._expert_modes:
            self.addItem(mode, mode)
        self.currentIndexChanged.connect(self._index_changed_cb)
    
    def _index_changed_cb(self):
        self.expertModeChanged.emit(self.currentData())


def test_collapse_button(tab, duration=300):
    layout = DefaultVLayout()

    layout.addWidget(QLabel('Test'))
    from .collapse_button import CollapsibleArea
    a = CollapsibleArea(tab, f'Test button (duration={duration})', duration=duration)
    a.setContentLayout(layout)   
    tab.addWidget(a)
    return a



















    
    

