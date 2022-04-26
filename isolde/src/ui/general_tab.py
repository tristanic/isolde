from .util import slot_disconnected
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

from Qt.QtGui import QIcon

from Qt.QtCore import Qt

from . import icon_dir, DEFAULT_ICON_SIZE

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
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        msl = DefaultHLayout()
        msl.addWidget(QLabel('Map: ', parent=mf))
        mcb = self.map_selector_combo_box = QComboBox(mf)
        msl.addWidget(mcb)
        mcb.currentIndexChanged.connect(self._map_chosen_cb)
        msl.addItem(DefaultSpacerItem())
        dmcb = self.is_difference_map_checkbox = QCheckBox(mf)
        dmcb.stateChanged.connect(self._is_difference_map_cb)
        dmcb.setText('Is a difference map')
        dmcb.setToolTip('Difference maps are shown with positive and negative contours. Ordinary maps are shown with a single positive contour.')
        msl.addWidget(dmcb)
        mdcb = self.enable_mdff_checkbox = QCheckBox(mf)
        mdcb.setToolTip('Choose whether this map contributes to the fitting potential during simulations')
        mdcb.stateChanged.connect(self._enable_mdff_cb)
        mdcb.setText('Enable MDFF')
        msl.addWidget(mdcb)
        ml.addLayout(msl)

        mwf = self.map_weight_frame = QWidget(mf)
        mwl = DefaultHLayout()
        mwl.addWidget(QLabel('Weight', parent=mf))
        mwsb = self.map_weight_spin_box = MapWeightSpinBox(mf)
        mwsb.setToolTip('Adjust how strongly this map "pulls" on atoms')
        mwl.addWidget(mwsb)
        units_label = QLabel(mf)
        units_label.setTextFormat(Qt.RichText)
        units_label.setText('<html><head/><body><p>x 1000 kJ mol<span style=" vertical-align:super;">-1</span> (map units)<span style=" vertical-align:super;">-1</span> Å<span style=" vertical-align:super;">3</span></p></body></html>')
        mwl.addWidget(units_label)
        mwl.addItem(DefaultSpacerItem())
        mwb = self.map_weight_set_button = QPushButton(mf)
        mwb.setText('Set')
        mwb.setToolTip('Click to apply any changes made to the weight box on left')
        mwl.addWidget(mwb)
        mwf.setLayout(mwl)
        ml.addWidget(mwf)

        buttons = self.button_frame = QWidget(mf)
        bl = DefaultHLayout()
        b_solid = self.solid_map_button = QToolButton(buttons)
        b_solid.setIcon(QIcon(os.path.join(icon_dir, 'mapsurf.png')))
        b_solid.setIconSize(DEFAULT_ICON_SIZE)
        b_solid.setToolTip('Set this map to opaque surface representation')
        b_solid.clicked.connect(self._solid_map_button_cb)
        bl.addWidget(b_solid)

        b_transp = self.transparent_map_button = QToolButton(buttons)
        b_transp.setIcon(QIcon(os.path.join(icon_dir, 'icecube.png')))
        b_transp.setIconSize(DEFAULT_ICON_SIZE)
        b_transp.setToolTip('Set this map to transparent surface representation')
        b_transp.clicked.connect(self._transparent_map_button_cb)
        bl.addWidget(b_transp)

        b_mesh = self.mesh_map_button = QToolButton(buttons)
        b_mesh.setIcon(QIcon(os.path.join(icon_dir, 'mesh.png')))
        b_mesh.setIconSize(DEFAULT_ICON_SIZE)
        b_mesh.setToolTip('Set this map to triangle mesh representation')
        b_mesh.clicked.connect(self._mesh_map_button_cb)
        bl.addWidget(b_mesh)

        b_color = self.map_color_button = QToolButton(buttons)
        b_color.setIcon(QIcon(os.path.join(icon_dir, 'rainbow.png')))
        b_color.setIconSize(DEFAULT_ICON_SIZE)
        b_color.setToolTip('Choose new color(s) for this map')
        b_color.clicked.connect(self._map_color_button_cb)
        bl.addWidget(b_color)
        bl.addItem(DefaultSpacerItem())
        
        buttons.setLayout(bl)
        ml.addWidget(buttons)


        ca = self.container
        ca.expanded.connect(self._expanded_cb)
        ca.collapsed.connect(self._collapsed_cb)
        if not ca.is_collapsed:
            self._expanded_cb()

    @property
    def current_map(self):
        return self.map_selector_combo_box.currentData()

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

    def _is_difference_map_cb(self, checked):
        self.current_map.is_difference_map = checked

    def _collapsed_cb(self, *_):
        self._selected_model_changed_handler.remove()
        self._selected_model_changed_handler = None

    def _expanded_cb(self, *_):
        self._update_map_selector_combo_box()
        self._selected_model_changed_handler = self.isolde.triggers.add_handler(
            'selected model changed', self._selected_model_changed_cb)

    def _selected_model_changed_cb(self, *_):
        self._update_map_selector_combo_box()
    
    def _enable_mdff_cb(self, checked):
        from chimerax.isolde.session_extensions import get_mdff_mgr
        mgr = get_mdff_mgr(self.isolde.selected_model, self.current_map)
        mgr.enabled = checked
    
    def _map_chosen_cb(self, *_):
        mcb = self.map_selector_combo_box
        mwb = self.map_weight_set_button
        mwf = self.map_weight_frame
        mwsb = self.map_weight_spin_box
        mdcb = self.enable_mdff_checkbox
        bf = self.button_frame
        dmcb = self.is_difference_map_checkbox #TODO: this currently does nothing

        mgr = None
        this_map = self.current_map
        if this_map is None:
            dmcb.setEnabled(False)
            bf.setEnabled(False)
        else:
            dmcb.setEnabled(True)
            bf.setEnabled(True)
            with slot_disconnected(self.is_difference_map_checkbox.stateChanged, self._is_difference_map_cb):
                dmcb.setChecked(this_map.is_difference_map)
            from chimerax.isolde.session_extensions import get_mdff_mgr
            mgr = get_mdff_mgr(self.isolde.selected_model, this_map)
        if mgr is None:
            mwf.setEnabled(False)
            mdcb.setEnabled(False)
            return
        mwf.setEnabled(True)
        mdcb.setEnabled(True)
        mdcb.setChecked(mgr.enabled)
        mwsb.setValue(mgr.global_k)

    def _solid_map_button_cb(self, *_):
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [self.current_map], **map_style_settings[map_styles.solid_opaque])

    def _transparent_map_button_cb(self, *_):
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        v = self.current_map
        if v.is_difference_map:
            style = map_styles.solid_t40
        else:
            style = map_styles.solid_t20
        volumecommand.volume(self.session, [v], **map_style_settings[style])

    def _mesh_map_button_cb(self, *_):
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [self.current_map], **map_style_settings[map_styles.mesh_triangle])

    def _map_color_button_cb(self, *_):
        v = self.current_map
        from Qt.QtWidgets import QColorDialog
        cd = QColorDialog(self.map_color_button)
        sa = cd.ShowAlphaChannel
        colors = []
        if v.is_difference_map:
            colors.append(cd.getColor(title="Colour for negative contour", options=sa))
            colors.append(cd.getColor(title="Colour for positive contour", options=sa))
        else:
            colors.append(cd.getColor(options=sa))

        import numpy
        carg = [
            numpy.array([c.red(), c.green(), c.blue(), c.alpha()], dtype=numpy.double)/255
                for c in colors
        ]
        v.set_parameters(surface_colors=carg)


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
