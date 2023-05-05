from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QGroupBox, QRadioButton, QWidget

class SimFidelityPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, "Simulation Fidelity/Speed", **kwargs)
        sfd = self.content = SimFidelityDialog(session, isolde, gui, self.content_area)
        self.setContentLayout(sfd.main_layout)

class SimFidelityDialog(UI_Panel_Base):
    from chimerax.isolde.openmm.sim_param_mgr import fidelity_modes


    _fidelity_settings = fidelity_modes.keys()

    _tooltips = {
        'Lowest/Fastest': '<span>Lowest fidelity, highest speed - useful for low-end GPUs when working in good density. Nonbonded interactions cutoff at 9 Å; vacuum electrostatics.</span>',
        'Medium/Medium': '<span>Moderate fidelity, high speed. Nonbonded interactions cutoff at 9 Å; implicit solvent electrostatics.</span>',
        'Highest/Slowest': '<span>Best quality, recommended when using a high-end GPU. Nonbonded interactions cutoff at 17 Å; implicit solvent electrostatics.</span>',
    }

    def __init__(self, session, isolde, gui, main_frame):
        super().__init__(session, isolde, gui, main_frame)
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        gb = QGroupBox(mf)
        il = DefaultHLayout()
        buttons = []
        for button_name in self._fidelity_settings:
            button = QRadioButton(button_name, parent=gb)
            button.setToolTip(self._tooltips[button_name])
            buttons.append(button)
            il.addWidget(button)
        gb.setLayout(il)
        ml.addWidget(gb)
        


        self._buttons = buttons
        for b in buttons:
            b.toggled.connect(self._button_click_cb)
        
        # For Nvidia GPUs, set fidelity level to highest. Otherwise set it to medium.
        from chimerax.isolde.benchmark import _opengl_info
        if 'nvidia' in _opengl_info(self.session).lower():
            buttons[-1].setChecked(True)
        else:
            buttons[1].setChecked(True)

    def sim_start_cb(self, trigger_name, data):
        self.main_frame.setEnabled(False)
    
    def sim_end_cb(self, trigger_name, data):
        self.main_frame.setEnabled(True)
    
    def _button_click_cb(self, *_):
        from chimerax.core.commands import run
        for b in self._buttons:
            if b.isChecked():
                run(self.session, f'isolde set simFidelityMode "{b.text()}"')
                return


