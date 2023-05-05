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
        # It is theoretically possible for an advanced user to install their own custom OpenMM Platform implementation(s).
        # Probably not advisable to make those the default, but we should definitely make them *available*
        unknown_platforms = [p for p in platform_names if p not in sim_params.platforms]
        if len(unknown_platforms):
            self.session.logger.info('ISOLDE detected that you have one or more non-standard OpenMM Platform(s) installed:\n'
                f'{",".join(unknown_platforms)}\n'
                'If you wish to use these for simulation purposes, you can do so by first choosing "Advanced" from the '
                'Experience Level drop-down menu at top right of the ISOLDE panel, then using the Computational Platform '
                'widget on the General tab.')
        platform_names = [p for p in platform_names if p not in unknown_platforms]
        platform_names = list(sorted(platform_names, key=lambda p:sim_params.platforms.index(p)))
        platform_names.extend(unknown_platforms)
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

