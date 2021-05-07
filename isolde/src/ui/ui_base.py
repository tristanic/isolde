# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



class UI_Panel_Base:
    def __init__(self, session, isolde, main_frame):
        self.isolde = isolde
        self.session = session
        self.main_frame = main_frame
        self._chimerax_trigger_handlers = []
        self._isolde_trigger_handlers = [
            isolde.triggers.add_handler('simulation started', self.sim_start_cb),
            isolde.triggers.add_handler('simulation terminated', self.sim_end_cb)
        ]

    def sim_start_cb(self, trigger_name, data):
        '''
        Override in derived class if panel behaviour should change on
        simulation start.
        '''
        pass


    def sim_end_cb(self, trigger_name, data):
        '''
        Override in derived class if panel behaviour should change on
        simulation end.
        '''
        pass

    def chimerax_models_changed(self, selected_model):
        '''
        Will be called by ISOLDE when the set of models loaded in ChimeraX
        changes.
        '''
        pass


    def enable(self):
        self.main_frame.setEnabled(True)

    def disable(self):
        self.main_frame.setEnabled(False)

    @property
    def enabled(self):
        return self.main_frame.isEnabled()

    def set_enabled(self, flag):
        if flag:
            self.enable()
        else:
            self.disable()

    def remove_trigger_handlers(self):
        for h in self._chimerax_trigger_handlers:
            self.session.triggers.remove_handler(h)
        for h in self._isolde_trigger_handlers:
            self.isolde.triggers.remove_handler(h)
