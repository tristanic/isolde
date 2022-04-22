from chimerax.ui.gui import MainToolWindow

from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFrame, QLayout, 
    QLabel, QSizePolicy, QComboBox, QSpacerItem,
    QPushButton, QTabWidget, QWidget, QScrollArea
)
from Qt import QtCore
from Qt.QtGui import QPixmap

from .ui_base import DefaultVLayout, DefaultHLayout, DefaultSpacerItem

from .util import slot_disconnected

import os
icon_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','icons'))


class IsoldeMainWin(MainToolWindow):
    def __init__(self, tool_instance, **kw):
        super().__init__(tool_instance, **kw)
        if hasattr(self.session, 'isolde'):
            self.isolde = self.session.isolde
        else:
            from ..isolde import Isolde
            self.isolde = Isolde(self.session)

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
        self.general_tab = GeneralTab(self.isolde, self.session, tabw, "General")
        self.validate_tab = IsoldeTab(tabw, "Validate")
        self.problems_tab = IsoldeTab(tabw, "Problem Zones")



    def _prepare_top_frame(self, main_layout):
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
        tcb = self.tutorials_combo_box = QComboBox(tw)
        tcb.setToolTip("Start an interactive tutorial")
        tcb.addItem('Tutorials')
        from ..tutorials import populate_tutorial_combo_box
        populate_tutorial_combo_box(tcb)
        tcb.currentIndexChanged.connect(self._show_tutorial_cb)


        l1.addWidget(tcb)
        l1.addItem(DefaultSpacerItem(200))
        li.addWidget(tw)

        bw = QWidget()
        l2 = DefaultHLayout()
        bw.setLayout(l2)

        wol = QLabel(bw)
        wol.setText('Working on: ')
        l2.addWidget(wol)
        mmcb = self.master_model_combo_box = QComboBox(bw)
        mmcb.setMinimumSize(QtCore.QSize(150,0))
        mmcb.setToolTip('Select the model to work on in ISOLDE')
        l2.addWidget(mmcb)

        from chimerax.core.models import ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS
        for event_type in (ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS):
            self._session_trigger_handlers.append(
                self.session.triggers.add_handler(event_type, self._update_model_list_cb)
            )

        l2.addItem(DefaultSpacerItem())

        hb = self.main_help_button = QPushButton(bw)
        hb.setText('Help')
        hb.clicked.connect(self.isolde.show_master_help_in_browser)
        l2.addWidget(hb)
        
        li.addWidget(bw)
        layout.addWidget(inner)

        main_layout.addWidget(tf)
    
    def _show_tutorial_cb(self, *_):
        tcb = self.tutorials_combo_box
        tpath = tcb.currentData()
        from ..tutorials import show_tutorial
        show_tutorial(self.session, tpath)
        with slot_disconnected(tcb.currentIndexChanged, self._show_tutorial_cb):
            tcb.setCurrentIndex(0)
    
    def _change_selected_model_cb(self, *_):
        mmcb = self.master_model_combo_box
        with slot_disconnected(mmcb.currentIndexChanged, self._change_selected_model_cb):
            m = mmcb.currentData()
            self.isolde.change_selected_model(m)
    
    def _update_model_list_cb(self, trigger_name, models):
        mmcb = self.master_model_combo_box
        from chimerax.core.models import Model
        if isinstance(models, Model):
            models = [models]
        from chimerax.atomic import AtomicStructure
        structures_need_update = False
        for m in models:
            # Use explicit type equality rather than isInstance because some previews 
            # are AtomicStructure subclasses
            if type(m) == AtomicStructure:
                structures_need_update = True
                break
        
        cm = self.isolde.selected_model
        if cm is not None and cm not in self.session.models.list():
            # Model has been deleted
            cm = None
            structures_need_update = True
        
        if structures_need_update:
            with slot_disconnected(mmcb.currentIndexChanged, self._change_selected_model_cb):
                mmcb.clear()
                models = [m for m in self.session.models.list() if type(m) == AtomicStructure]
                models = sorted(models, key=lambda m: m.id)

                for m in models:
                    id_str = f'{m.id_string}. {m.name}'
                    mmcb.addItem(id_str, m)
                if cm is not None:
                    mmcb.setCurrentIndex(mmcb.findData(cm))

    def cleanup(self):
        for h in self._session_trigger_handlers:
            self.session.triggers.remove_handler(h)







def test_collapse_button(tab, duration=300):
    layout = DefaultVLayout()

    layout.addWidget(QLabel('Test'))
    from .collapse_button import CollapsibleArea
    a = CollapsibleArea(tab, f'Test button (duration={duration})', duration=duration)
    a.setContentLayout(layout)   
    tab.addWidget(a)
    return a



















    
    

