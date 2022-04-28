from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QGroupBox, QRadioButton, QWidget

class SimFidelityPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui):
        super().__init__(gui, parent, "Sim Fidelity/Speed")
        sfd = self.content = SimFidelityDialog(session, isolde, gui, self.content_area)
        self.setContentLayout(sfd.main_layout)

class SimFidelityDialog(UI_Panel_Base):

    _fidelity_settings = {
        'Quick': {
            'nonbonded_cutoff_distance': 0.9,
            'use_gbsa': False,
        },
        'Medium': {
            'nonbonded_cutoff_distance': 0.9,
            'use_gbsa': True,
            'gbsa_cutoff': 1.1,
        },
        'High': {
            'nonbonded_cutoff_distance': 1.7,
            'use_gbsa': True,
            'gbsa_cutoff': 2.0,
        }

    }

    _tooltips = {
        'Quick': 'Lowest fidelity, highest speed - useful for low-end GPUs when working in good density. Nonbonded interactions cutoff at 9 Å; vacuum electrostatics.',
        'Medium': 'Moderate fidelity, high speed. Nonbonded interactions cutoff at 9 Å; implicit solvent electrostatics. ',
        'High': 'Best quality, recommended when using a high-end GPU. Nonbonded interactions cutoff at 17 Å; implicit solvent electrostatics.',
    }

    def __init__(self, session, isolde, gui, main_frame):
        super().__init__(session, isolde, gui, main_frame)
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        
        gb = QGroupBox(mf)
        il = DefaultHLayout()
        buttons = []
        for button_name in self._fidelity_settings.keys():
            button = QRadioButton(button_name, parent=gb)
            button.setToolTip(self._tooltips[button_name])
            buttons.append(button)
            il.addWidget(button)
        gb.setLayout(il)
        ml.addWidget(gb)
        


        self._buttons = buttons
        for b in buttons:
            b.toggled.connect(self._button_click_cb)
        buttons[-1].setChecked(True)

    def sim_start_cb(self, trigger_name, data):
        self.main_frame.setEnabled(False)
    
    def sim_end_cb(self, trigger_name, data):
        self.main_frame.setEnabled(True)
    
    def _button_click_cb(self, *_):
        for b in self._buttons:
            if b.isChecked():
                params = self._fidelity_settings[b.text()]
                break
        sp = self.isolde.sim_params
        for key, value in params.items():
            sp.set_param(key, value)


