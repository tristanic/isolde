import os

from ..util import slot_disconnected
from ...util import block_managed_trigger_handler

from ..collapse_button import CollapsibleArea
from ..ui_base import (
    UI_Panel_Base, 
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem,
    QComboBox
)

from .. import icon_dir

from Qt.QtWidgets import (
    QLabel, QSlider, QDoubleSpinBox, QCheckBox, QWidget
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
        tsl = self.temperature_slider = TemperatureSlider(isolde, Qt.Orientation.Horizontal, parent=mf)
        tl.addWidget(tsl)
        tsb = self.temperature_spinbox = TemperatureSpinBox(isolde, mf)
        tl.addWidget(tsb)
        tl.addWidget(QLabel(' K', parent=mf))
        ml.addLayout(tl)
        sw = self.smoothing_widget = SmoothingWidget(isolde, mf)
        ml.addWidget(sw)
    
    def cleanup(self):
        self.temperature_slider.cleanup()
        self.temperature_spinbox.cleanup()
        self.smoothing_widget.cleanup()
        super().cleanup()


class TemperatureSlider(QSlider):
    _stylesheet = '''
QSlider::groove:horizontal {
    height: 10px;
    background-color: qlineargradient(x1:0, y1:0, x2:1, y1:0, stop: 0 #1030ff, stop: 0.2 #b465da, stop:0.4 #cf6cc9, stop:0.5 #ee609c, stop:1.0 #ff5040);
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
    def __init__(self, isolde, orientation, parent=None):
        super().__init__(orientation, parent=parent)
        self.setStyleSheet(self._stylesheet)
        self.setMinimumHeight(20)
        self.setMinimum(0)
        self.setMaximum(200)
        self.setSingleStep(1)
        self.setTickPosition(QSlider.TickPosition.NoTicks)
        self.isolde = isolde

        self._exponent_multiplier = 0.031073 # gives T=500 at position 200
        self.valueChanged.connect(self._slider_changed_cb)
        self._temperature_changed_handler = isolde.sim_params.triggers.add_handler(
            isolde.sim_params.PARAMETER_CHANGED, self._parameter_changed_cb
        )
        # Initialize to ISOLDE's current value
        self._parameter_changed_cb('', ('temperature', isolde.sim_params.temperature))
    
    def _slider_changed_cb(self, value):
        if value==0:
            t = 0
        else:
            from math import exp, floor
            t = floor(exp(self._exponent_multiplier*value))
        with block_managed_trigger_handler(self, '_temperature_changed_handler'):
            self.isolde.sim_params.temperature = t

    def _parameter_changed_cb(self, trigger_name, data):
        param, value = data
        if param != 'temperature':
            return
        from openmm import unit
        value = value.value_in_unit(unit.kelvin)
        with slot_disconnected(self.valueChanged, self._slider_changed_cb):
            if value<1:
                self.setValue(0)
            else:
                from math import log, floor
                tick = floor(log(value)/self._exponent_multiplier)
                self.setValue(tick)

    def cleanup(self):
        self._tempreature_changed_handler.remove()    
        
        
class TemperatureSpinBox(QDoubleSpinBox):
    def __init__(self, isolde, parent):
        super().__init__(parent)
        self.isolde = isolde
        self.valueChanged.connect(self._value_changed_cb)
        self.setKeyboardTracking(False)
        self.setDecimals(0)
        self.setMaximum(500)
        self._param_changed_handler = isolde.sim_params.triggers.add_handler(
            isolde.sim_params.PARAMETER_CHANGED, 
            self._parameter_changed_cb
        )
        # Initialize to ISOLDE's current value
        self._parameter_changed_cb('', ('temperature', isolde.sim_params.temperature))
    
    def _value_changed_cb(self, value):
        with block_managed_trigger_handler(self, '_param_changed_handler'):
            self.isolde.sim_params.temperature = value
    
    def _parameter_changed_cb(self, trigger_name, data):
        param, val = data
        if param != 'temperature':
            return
        from openmm import unit
        val = val.value_in_unit(unit.kelvin)
        with slot_disconnected(self.valueChanged, self._value_changed_cb):
            self.setValue(val)
    
    def cleanup(self):
        self._param_changed_handler.remove()
        

class SmoothingWidget(QWidget):
    _stylesheet = f'''
QSlider::groove:horizontal {{
    border-image: url({os.path.join(icon_dir,'smoothing.png').replace(os.sep, '/')}) 0 0 0 0 stretch stretch;
    border: 1px solid black;
    position: absolute;
    left: 10px;
    right: 10px;
}}
QSlider::handle:horizontal {{
    width: 10px;
    background: #0b1707;
    border: 1px solid #46992b;
    margin: 0px -10px;
    border-radius: 5px;
}}
QSlider::handle:horizontal:hover {{
    background-color: #46992b;
}}
'''
    def __init__(self, isolde, parent=None):
        super().__init__(parent=parent)
        self.isolde = isolde
        ml = self.main_layout = DefaultHLayout()
        ml.addWidget(QLabel('Trajectory smoothing: ', parent=self))
        ssl = self.slider = QSlider(Qt.Orientation.Horizontal, parent=self)
        ssl.setStyleSheet(self._stylesheet)
        ssl.setMinimumHeight(20)
        ssl.setMinimum(1)
        ssl.setMaximum(100)
        ssl.setValue(round(isolde.sim_params.smoothing_alpha*100))
        ssl.valueChanged.connect(self._value_changed_cb)
        ml.addWidget(ssl)
        scb = self.enable_smoothing_checkbox = QCheckBox(parent=self)
        scb.setChecked(isolde.sim_params.trajectory_smoothing)
        scb.toggled.connect(self._smoothing_checkbox_cb)
        ml.addWidget(scb)
        self._param_changed_handler = isolde.sim_params.triggers.add_handler(
            isolde.sim_params.PARAMETER_CHANGED, self._parameter_changed_cb
        )
        self.setLayout(ml)



    
    def _value_changed_cb(self, value):
        with block_managed_trigger_handler(self, '_param_changed_handler'):
            self.isolde.sim_params.smoothing_alpha = value/100
    
    def _parameter_changed_cb(self, trigger_name, data):
        param, val = data
        if param == 'trajectory_smoothing':
            with slot_disconnected(self.enable_smoothing_checkbox, self._smoothing_checkbox_cb):
                self.enable_smoothing_checkbox.setChecked(val)
        elif param == 'smoothing_alpha':
            with slot_disconnected(self.slider, self._value_changed_cb):
                from chimerax.isolde.constants import defaults
                val = max(defaults.SMOOTHING_ALPHA_MIN, max(val, defaults.SMOOTHING_ALPHA_MAX))
                self.slider.setValue(round(val*100))

    def _smoothing_checkbox_cb(self, flag):
        with block_managed_trigger_handler(self, '_param_changed_handler'):
            self.isolde.sim_params.trajectory_smoothing = flag

    def cleanup(self):
        self._param_changed_handler.remove()




