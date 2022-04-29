from ..util import slot_disconnected
from ..ui_base import (
    DefaultSpacerItem,
    UI_Panel_Base,
    DefaultHLayout
)

from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QLabel, QDoubleSpinBox

class MaskAndSpotlightSettingsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Mask and spotlight settings', **kwargs)
        mssd = self.content = MaskAndSpotlightSettingsDialog(session, isolde, gui, self)
        self.setContentLayout(mssd.main_layout)


class MaskAndSpotlightSettingsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area)
        self.container = collapse_area

        params = self.isolde.params
        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()

        ml.addWidget(QLabel('Mask radius (Å): ', parent=mf))
        mrsb = self.mask_radius_spinbox = QDoubleSpinBox(mf)
        mrsb.setValue(params.map_mask_radius)
        mrsb.setDecimals(1)
        mrsb.setSingleStep(0.5)
        mrsb.setMinimum(0.5)
        mrsb.setMaximum(100)
        mrsb.valueChanged.connect(self._mask_radius_changed_cb)
        ml.addWidget(mrsb)

        ml.addItem(DefaultSpacerItem())

        ml.addWidget(QLabel('Spotlight radius (Å): ', parent=mf))
        slsb = self.spotlight_radius_spinbox = QDoubleSpinBox(mf)
        slsb.setValue(params.spotlight_radius)
        slsb.setDecimals(1)
        slsb.setSingleStep(1)
        slsb.setMinimum(5)
        slsb.setMaximum(1000)
        slsb.valueChanged.connect(self._spotlight_radius_spinbox_changed_cb)
        ml.addWidget(slsb)

        self._selected_model_changed_cb('', self.isolde.selected_model)

    def _selected_model_changed_cb(self, _, selected_model):
        h = getattr(self, '_spotlight_radius_change_handler', None)
        if h is not None:
            h.remove()
        if selected_model is not None:
            from chimerax.clipper import get_map_mgr
            mgr = get_map_mgr(selected_model)
            if mgr is not None:
                self._spotlight_radius_change_handler = mgr.triggers.add_handler('spotlight changed', self._spotlight_radius_changed_cb)
            else:
                self._spotlight_radius_change_handler = None

    def _mask_radius_changed_cb(self, value):
        self.isolde.params.map_mask_radius = value

    def _spotlight_radius_spinbox_changed_cb(self, value):
        self.session.isolde.params.spotlight_radius=value
        from chimerax.core.commands import run
        run(self.session, f'clipper spotlight radius {value:.1f}', log=False)

    def _spotlight_radius_changed_cb(self, trigger_name, params):
        center, radius = params
        slsb = self.spotlight_radius_spinbox
        with slot_disconnected(slsb.valueChanged, self._spotlight_radius_spinbox_changed_cb):
            slsb.setValue(radius)
            self.session.isolde.params.spotlight_radius = radius
    
    def cleanup(self):
        super().cleanup()
        h = getattr(self, '_spotlight_radius_change_handler', None)
        if h is not None:
            h.remove()