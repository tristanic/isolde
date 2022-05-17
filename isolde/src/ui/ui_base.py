# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from Qt.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QSpacerItem, QSizePolicy,
    QScrollArea, QWidget, QFrame
    )
from Qt import QtCore
from Qt.QtGui import QColor

import sys
if 'win' in sys.platform.lower():
    from .util import WinAutoResizeQComboBox as QComboBox
else:
    from Qt.QtWidgets import QComboBox


import os
_base_path = os.path.dirname(os.path.abspath(__file__))

class DefaultVLayout(QVBoxLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setContentsMargins(0,2,0,0)
        self.setSpacing(3)

class DefaultHLayout(QHBoxLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setContentsMargins(0,0,0,0)
        self.setSpacing(3)

class DefaultSpacerItem(QSpacerItem):
    def __init__(self, width=100):
        super().__init__(width, 0, QSizePolicy.Expanding, QSizePolicy.Minimum)

class ExpertModeSelector(QComboBox):
    DEFAULT=0
    ADVANCED=1
    DEVELOPER=2
    _expert_modes = {
        DEFAULT: ['Default', [230,230,230]], 
        ADVANCED: ['Advanced', [215, 237, 255]], 
        DEVELOPER: ['Developer', [255, 215, 215]],
    }
    with open(os.path.join(_base_path, 'intermediate.qss'), 'rt') as intfile, open(os.path.join(_base_path, 'developer.qss'), 'rt') as devfile:
        stylesheets = {
            DEFAULT: "",
            ADVANCED: intfile.read(),
            DEVELOPER: devfile.read()
        }
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for val, (mode, color) in self._expert_modes.items():
            self.addItem(mode, val)
        self.currentIndexChanged.connect(self._update_stylesheet)
        self._update_stylesheet(0)
    
    def _update_stylesheet(self, index):
        self.setStyleSheet(self.stylesheets[index])
        for val, (mode, color) in self._expert_modes.items():
            if color is not None:
                self.setItemData(val, QColor(*color), QtCore.Qt.BackgroundRole)

class UI_Panel_Base:
    def __init__(self, session, isolde, gui, main_frame, sim_sensitive=True, expert_level = ExpertModeSelector.DEFAULT):
        self.gui = gui
        self.isolde = isolde
        self.session = session
        self.main_frame = main_frame
        self.expert_level = expert_level
        self._set_expert_level()
        self._chimerax_trigger_handlers = []
        self._isolde_trigger_handlers = []
        if sim_sensitive:
            self._isolde_trigger_handlers.extend( [
                isolde.triggers.add_handler('simulation started', self.sim_start_cb),
                isolde.triggers.add_handler('simulation terminated', self.sim_end_cb)
            ])
        gui.register_panel(self)

    def _set_expert_level(self):
        el = self.expert_level
        emcb = self.gui.expert_mode_combo_box
        if el > ExpertModeSelector.DEFAULT:
            self.main_frame.setStyleSheet(ExpertModeSelector.stylesheets[el])
            emcb.currentIndexChanged.connect(self._expert_level_changed_cb)
            self._expert_level_changed_cb()
    
    def _expert_level_changed_cb(self, *_):
        emcb = self.gui.expert_mode_combo_box
        el = emcb.currentData()
        display = (el >= self.expert_level)
        self.main_frame.setVisible(display)
        



    def sim_start_cb(self, trigger_name, data):
        '''
        Override in derived class if panel behaviour should change on
        simulation start.
        '''
        pass


    def sim_end_cb(self, trigger_name, data):
        '''
        Override in derived class if panel behaviour should change on
        simulation end.
        '''
        pass


    def chimerax_models_changed(self, selected_model):
        '''
        Will be called by ISOLDE when the set of models loaded in ChimeraX
        changes.
        '''
        pass
    

    def enable(self):
        self.main_frame.setEnabled(True)

    def disable(self):
        self.main_frame.setEnabled(False)

    @property
    def enabled(self):
        return self.main_frame.isEnabled()

    def set_enabled(self, flag):
        if flag:
            self.enable()
        else:
            self.disable()

    def cleanup(self):
        '''
        Called when the GUI is closed to clean up all callbacks.
        '''
        self._remove_trigger_handlers()

    def _remove_trigger_handlers(self):
        for h in self._chimerax_trigger_handlers:
            h.remove()
        for h in self._isolde_trigger_handlers:
            h.remove()



class IsoldeTab(QWidget):
    def __init__(self, gui, tab_widget, tab_name):
        super().__init__()
        self.gui = gui
        self.tab_widget = tab_widget
        tab_widget.addTab(self, tab_name)
        hl = DefaultHLayout()
        self.setLayout(hl)
        sa = self.scroll_area = QScrollArea(self)
        hl.addWidget(sa)
        sa.setWidgetResizable(True)
        mf = self.main_frame = QFrame(sa)
        ml = self.main_layout = DefaultVLayout()
        ml.setContentsMargins(1,1,3,1)
        ml.addStretch()
        mf.setLayout(ml)
        sa.setWidget(mf)
    
    def addWidget(self, widget):
        # Last position is the spacer, so we always want to go just before that.
        self.main_layout.insertWidget(self.main_layout.count()-1, widget)
    
