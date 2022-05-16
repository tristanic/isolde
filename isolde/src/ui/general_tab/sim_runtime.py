from ..util import slot_disconnected
from ...util import block_managed_trigger_handler

from ..collapse_button import CollapsibleArea
from ..ui_base import (
    UI_Panel_Base, 
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem,
    QComboBox
)

from Qt.QtWidgets import (
    QLabel, QSlider, QDoubleSpinBox
)

from Qt.QtCore import Qt

class SimRuntimePanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Simulation Runtime', **kwargs)
        srd = self.content = SimRuntimeDialog(session, isolde, gui, self)
        self.setContentLayout(srd.main_layout)

class SimRuntimeDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=False)
        
        self.container=collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        tl = DefaultHLayout()
        tl.addWidget(QLabel('Temperature: ', parent=mf))
        tsl = self.temperature_slider = QSlider(Qt.Orientation.Horizontal, parent=mf)
        tsl.setTickPosition(QSlider.TickPosition.NoTicks)
        tsl.setStyleSheet(_temperature_stylesheet)
        tsl.setMinimum(0)
        tsl.setMaximum(100)
        tsl.setSingleStep(0)
        tl.addWidget(tsl)
        tsb = self.temperature_spinbox = QDoubleSpinBox(mf)
        tl.addWidget(tsb)
        tl.addWidget(QLabel(' K', parent=mf))
        ml.addLayout(tl)


_temperature_stylesheet = '''
QSlider::groove:horizontal {
    height: 10px;
    background-color: qlineargradient(x1:0, y1:0, x2:1, y1:0, stop: 0 #b465da, stop:0.4 #cf6cc9, stop:0.5 #ee609c, stop:1.0 #ff5040)
    border-radius: 5px;
    position: absolute;
    left: 10px;
    right: 10px;
}
QSlider::handle:horizontal {
    width: 10px;
    background: #0b1707;
    border: 1px solid #46992b;
    margin: 0px -10px;
    border-radius: 5px;
}
QSlider::handle:horizontal:hover {
    background-color: #46992b;
}
'''

