from .ui_base import UI_Panel_Base

class SimFidelityPanel(UI_Panel_Base):

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

    def __init__(self, session, isolde, gui, main_frame, buttons):
        super().__init__(session, isolde, gui, main_frame)
        self._buttons = buttons
        for b in buttons:
            b.clicked.connect(self._button_click_cb)

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
            setattr(sp, key, value)


