from .collapse_button import CollapsibleArea
from .ui_base import (
    UI_Panel_Base, IsoldeTab,
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem
)

import sys
if 'win' in sys.platform.lower():
    from .util import WinAutoResizeQComboBox as QComboBox
else:
    from Qt.QtWidgets import QComboBox

from Qt.QtWidgets import (
    QLabel, QToolButton, QPushButton, QCheckBox,
    QDoubleSpinBox, QWidget
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
    def __init__(self, session, isolde, gui, tab_widget, tab_name):
        super().__init__(gui, tab_widget, tab_name)
        self.session = session
        self.isolde = isolde
        msp = self.map_settings_panel = CollapsibleArea(parent=self.scroll_area, title='Map settings')
        msd = self.map_settings_dialog = MapSettingsDialog(session, isolde, gui, msp)
        msp.setContentLayout(msd.main_layout)
        self.addWidget(msp)


class MapSettingsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        msl = self.map_selector_layout = DefaultHLayout()
        msl.addWidget(QLabel('Map: ', parent=mf))
        mcb = self.map_selector_combo_box = QComboBox(mf)
        msl.addWidget(mcb)
        mcb.currentIndexChanged.connect(self._map_chosen_cb)
        msl.addItem(DefaultSpacerItem())
        dmcb = self.is_difference_map_checkbox = QCheckBox(mf)
        dmcb.setText('Is a difference map')
        msl.addWidget(dmcb)
        mdcb = self.enable_mdff_checkbox = QCheckBox(mf)
        mdcb.setText('Enable MDFF')
        msl.addWidget(mdcb)
        ml.addLayout(msl)

        mwf = self.map_weight_frame = QWidget(mf)
        mwl = self.map_weight_layout = DefaultHLayout()
        mwl.addWidget(QLabel('Weight', parent=mf))
        mwsb = self.map_weight_spin_box = MapWeightSpinBox(mf)
        mwl.addWidget(mwsb)
        units_label = QLabel(mf)
        units_label.setTextFormat(Qt.RichText)
        units_label.setText('<html><head/><body><p>x 1000 kJ mol<span style=" vertical-align:super;">-1</span> (map units)<span style=" vertical-align:super;">-1</span> Å<span style=" vertical-align:super;">3</span></p></body></html>')
        mwl.addWidget(units_label)
        mwl.addItem(DefaultSpacerItem())
        mwb = self.map_weight_set_button = QPushButton(mf)
        mwb.setText('Set')
        mwl.addWidget(mwb)
        mwf.setLayout(mwl)
        ml.addWidget(mwf)

        ca = self.container
        ca.expanded.connect(self._expanded_cb)
        ca.collapsed.connect(self._collapsed_cb)
        if not ca.is_collapsed:
            self._expanded_cb()
        
    def _update_map_selector_combo_box(self, *_):
        mcb = self.map_selector_combo_box
        mcb.clear()
        sm = self.isolde.selected_model
        if sm is None:
            return
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(sm)
        if mgr is None:
            return
        for v in mgr.all_maps:
            mcb.addItem(v.name, v)
        
    def _collapsed_cb(self, *_):
        self._selected_model_changed_handler.remove()
        self._selected_model_changed_handler = None

    def _expanded_cb(self, *_):
        self._update_map_selector_combo_box()
        self._selected_model_changed_handler = self.isolde.triggers.add_handler(
            'selected model changed', self._selected_model_changed_cb)

    def _selected_model_changed_cb(self, *_):
        self._update_map_selector_combo_box()

    def _map_is_difference_map_cb(self, checked):
        ''' TODO: implement difference map setting in Clipper '''
        mcb = self.map_selector_combo_box
        map = mcb.currentData()
        text = {True: "", False: "not"}
        print(f'Map {map.id_string} is {text[checked]} a difference map.')
    
    def _map_chosen_cb(self, *_):
        mcb = self.map_selector_combo_box
        mwb = self.map_weight_set_button
        mwf = self.map_weight_frame
        mwsb = self.map_weight_spin_box
        mdcb = self.enable_mdff_checkbox
        dmcb = self.is_difference_map_checkbox #TODO: this currently does nothing

        mgr = None
        this_map = mcb.currentData()
        if this_map is None:
            dmcb.setEnabled(False)
        else:
            from chimerax.isolde.session_extensions import get_mdff_mgr
            mgr = get_mdff_mgr(self.isolde.selected_model, this_map)
        if mgr is None:
            mwf.setEnabled(False)
            return
        dmcb.setEnabled(True)
        mwf.setEnabled(True)
        mdcb.setChecked(mgr.enabled)
        mwsb.setValue(mgr.global_k)






    def cleanup(self):
        if self._selected_model_changed_handler is not None:
            self._selected_model_changed_handler.remove()
        super().cleanup()

    def _maps_changed_cb(self, *_):
        pass

class MapWeightSpinBox(QDoubleSpinBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setKeyboardTracking(False)
        self.setMinimumWidth(50)
        self.setMaximum(1e5)
        self.setMinimum(1e-6)
        self.valueChanged.connect(self.update_display_and_step)

    def _displayed_decimal_places_and_step(self, number, sig_figs=3):
        from math import log, ceil, floor
        if number <= 0:
            return 2, 0.1
        places = max(ceil(-log(number, 10)+sig_figs-1), 0)
        step = 10**(floor(log(number, 10)-1))
        return places, step

    def update_display_and_step(self, *_):
        v = self.value()
        dps, step = self._displayed_decimal_places_and_step(v)
        self.setDecimals(dps)
        self.setSingleStep(step)
