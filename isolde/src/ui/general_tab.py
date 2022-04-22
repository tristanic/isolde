from .collapse_button import CollapsibleArea
from .ui_base import (
    UI_Panel_Base, IsoldeTab,
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem
)

from Qt.QtWidgets import (
    QComboBox, QLabel, QToolButton, QPushButton, QCheckBox,
    QDoubleSpinBox
)

from Qt.QtCore import Qt

def populate_general_tab(session, isolde, tab):
    parent = tab.scroll_area
    map_settings_panel = CollapsibleArea(parent, title='Map settings')
    msd = MapSettingsDialog(session, isolde, map_settings_panel.content_area)
    map_settings_panel.content = msd
    map_settings_panel.setContentLayout(msd.main_layout)
    
    tab.addWidget(map_settings_panel)

class GeneralTab(IsoldeTab):
    def __init__(self, session, isolde, tab_widget, tab_name):
        super().__init__(tab_widget, tab_name)
        self.session = session
        self.isolde = isolde
        msp = self.map_settings_panel = CollapsibleArea(parent=self.scroll_area, title='Map settings')
        msd = self.map_settings_dialog = MapSettingsDialog(session, isolde, msp.content_area)
        msp.setContentLayout(msd.main_layout)
        self.addWidget(msp)


class MapSettingsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, main_frame, sim_sensitive=False):
        super().__init__(session, isolde, main_frame, sim_sensitive=sim_sensitive)
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        msl = self.map_selector_layout = DefaultHLayout()
        msl.addWidget(QLabel('Map: ', parent=mf))
        mcb = self.map_selector_combo_box = QComboBox(mf)
        msl.addWidget(mcb)
        msl.addItem(DefaultSpacerItem())
        dmcb = self.is_difference_map_checkbox = QCheckBox(mf)
        dmcb.setText('Is a difference map')
        msl.addWidget(dmcb)
        mdcb = self.enable_mdff_checkbox = QCheckBox(mf)
        mdcb.setText('Enable MDFF')
        msl.addWidget(mdcb)
        ml.addLayout(msl)

        mwl = self.map_weight_layout = DefaultHLayout()
        mwl.addWidget(QLabel('Weight', parent=mf))
        mwsb = self.map_weight_spin_box = QDoubleSpinBox(mf)
        mwsb.setDecimals(3)
        mwsb.setMaximum(10000.0)
        mwsb.setSingleStep(0.1)
        mwl.addWidget(mwsb)
        units_label = QLabel(mf)
        units_label.setTextFormat(Qt.RichText)
        units_label.setText('<html><head/><body><p>x 1000 kJ mol<span style=" vertical-align:super;">-1</span> (map units)<span style=" vertical-align:super;">-1</span> Å<span style=" vertical-align:super;">3</span></p></body></html>')
        mwl.addWidget(units_label)
        mwl.addItem(DefaultSpacerItem())
        mwb = self.map_weight_set_button = QPushButton(mf)
        mwb.setText('Set')
        mwl.addWidget(mwb)
        ml.addLayout(mwl)











    def _maps_changed_cb(self, *_):
        pass

