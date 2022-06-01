# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 28-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



from .ui_base import UI_Panel_Base

class UI_Live_Map_Controls(UI_Panel_Base):
    def __init__(self, session, isolde, gui, main_frame,
        live_recalc_checkbox, manual_recalc_button,
        update_sim_potential_button):
        super().__init__(session, isolde, gui, main_frame)
        self._live_checkbox = live_recalc_checkbox
        self._recalc_button = manual_recalc_button
        self._update_sim_button = update_sim_potential_button
        self.set_enabled(self._enable_check(isolde.selected_model))
        self._isolde_trigger_handlers.append(
            isolde.triggers.add_handler(isolde.SELECTED_MODEL_CHANGED,
                self._model_changed_cb)
        )

    def _enable_check(self, selected_model):
        if selected_model is None:
            return False
        from chimerax.clipper.maps import XmapHandler_Live
        for m in selected_model.all_models():
            if isinstance(m, XmapHandler_Live):
                return True
        return False

    def sim_start_cb(self, trigger_name, data):
        if self.enabled:
            self._update_sim_button.setEnabled(True)

    def sim_end_cb(self, trigger_name, data):
        self._update_sim_button.setEnabled(False)

    def selected_model_changed_cb(self, trigger_name, selected_model):
        self.set_enabled(self._enable_check(selected_model))

    def _model_changed_cb(self, trigger_name, model):
        self.set_enabled(self._enable_check(model))

    @property
    def live_recalc(self):
        return self._live_checkbox.checkState()

    @live_recalc.setter
    def live_recalc(self, flag):
        self._live_checkbox.setChecked(flag)
        self.isolde.set_live_xmap_recalc(flag)

    def enable(self):
        self.main_frame.setVisible(True)
        super().enable()

    def disable(self):
        self.main_frame.setVisible(False)
        super().disable()
