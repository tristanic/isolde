from ..ui_base import (
    DefaultSpacerItem,
    UI_Panel_Base, 
    DefaultHLayout, DefaultVLayout,
    ExpertModeSelector
)
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QLabel, QComboBox

class ComputationalPlatformPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Computational Platform', expert_level=ExpertModeSelector.ADVANCED, **kwargs)
        cpd = self.content = ComputationalPlatformDialog(session, isolde, gui, self)
        self.setContentLayout(cpd.main_layout)


class ComputationalPlatformDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area)
        self.container = collapse_area
        
        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()
        pcb = self.platform_combo_box = QComboBox(mf)
        self._populate_platform_combo_box()
        pcb.currentIndexChanged.connect(self._choose_platform_cb)
        ml.addWidget(pcb)
        ml.addItem(DefaultSpacerItem())

    
    def _populate_platform_combo_box(self):
        from chimerax.isolde.openmm.openmm_interface import get_available_platforms
        sim_params = self.session.isolde.sim_params
        platform_names = get_available_platforms()
        platform_names = list(sorted(platform_names, key=lambda p:sim_params.platforms.index(p)))
        pcb = self.platform_combo_box
        pcb.addItems(platform_names)

        if "CUDA" not in platform_names and "OpenCL" not in platform_names:
            self.session.logger.warning('WARNING: no OpenCL or compatible CUDA '
                'drivers detected! While it is theoretically possible to run '
                'ISOLDE using CPU only, in practice it is prohibitively slow. '
                'If you have a suitable GPU in your machine, please check that you '
                'have the recommended drivers from the manufacturer installed. '
                'The current required CUDA version is 11.2 - if installed, please '
                'make sure this is on your library path before starting ChimeraX.')
        
        # Set to the preferred or, failing that, fastest available platform
        if sim_params.platform in platform_names:
            pcb.setCurrentIndex(pcb.findText(sim_params.platform))
        else:
            for p in sim_params.platforms:
                if p in platform_names:
                    pcb.setCurrentIndex(pcb.findText(p))
                    sim_params.platform = p
                    break
    
    def _choose_platform_cb(self, *_):
        self.session.isolde.sim_params.platform = self.platform_combo_box.currentText()

