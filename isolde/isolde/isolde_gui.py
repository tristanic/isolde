


class IsoldeGui():
    
    def __init__(self, session):
        self.session = session
        from PyQt5 import QtWidgets
        
        self.mainwin = QtWidgets.QDockWidget(parent=session.ui.main_window)
        self.mainwin.setFloating(True) 
        from . import isoldewidget
        
        self.iw = isoldewidget.Ui_isolde_widget()
        self.iw.setupUi(self.mainwin)
        # Should load saved state here
        
        # Define frames specific to crystallograpy, EM or free mode
        self._xtal_frames = [
            self.iw._sim_basic_xtal_map_frame
            ]
        self._em_frames = [
            self.iw._sim_basic_em_map_frame
            ]
        self._free_frames = [
            ]
        
        self._sim_mode_frame_lists = [
            self._xtal_frames,
            self._em_frames,
            self._free_frames
            ]
        
        # Radio buttons to choose different selection modes
        self._selection_mode_buttons = [
            self.iw._sim_basic_by_selected_atoms_button,
            self.iw._sim_basic_by_chain_button,
            self.iw._sim_basic_whole_structure_button,
            self.iw._sim_basic_custom_selection_button
            ]        
        
        
        # Define intermediate and expert frames
        self._intermediate_frames = [
            self.iw._sim_platform_frame
            ]
        self._expert_frames = [
            self.iw._force_field_selection_frame,
            self.iw._sim_basic_custom_selection_button
            ]
        
        
        # Apply custom palettes to intermediate and expert frames
        from . import palettes
        self._pi = palettes.IntermediatePalette()
        self._pe = palettes.ExpertPalette()
        
        for f in self._intermediate_frames:
            f.setPalette(self._pi.palette)
            f.setAutoFillBackground(True)
        
        for f in self._expert_frames:
            f.setPalette(self._pe.palette)
            f.setAutoFillBackground(True)


    def _change_experience_level_or_sim_mode(self):
        exp_index = self.iw._experience_level_combo_box.currentIndex()
        mode_index = self.iw._sim_basic_mode_combo_box.currentIndex()
        # Need to consider both at once to ensure we don't show/hide
        # something we shouldn't. For the simulation mode, we need to
        # ensure all frames associated with the *other* two modes remain
        # hidden.
        import copy
        hide_sim_modes = copy.copy(self._sim_mode_frame_lists)
        hide_sim_modes.pop(mode_index)
        # Flatten to a single list for easy searching
        hide_sim_modes = [item for sublist in hide_sim_modes for item in sublist]
        for f in hide_sim_modes:
            f.hide()
        show_sim_modes = self._sim_mode_frame_lists[mode_index]
        if (exp_index == 0):
            # Easy. Just hide everything intermediate or expert, and
            # everything belonging to other sim modes
            for f in self._intermediate_frames:
                f.hide()
            for f in self._expert_frames:
                f.hide()
            for f in show_sim_modes:
                if f not in self._intermediate_frames and \
                   f not in self._expert_frames:
                    f.show()
        elif (exp_index == 1):
            for f in self._intermediate_frames:
                if f not in hide_sim_modes:
                    f.show()
            for f in self._expert_frames:
                f.hide()
            for f in show_sim_modes:
                if f not in self._expert_frames:
                    f.show()
        else:
            for f in self._intermediate_frames:
                if f not in hide_sim_modes:
                    f.show()
            for f in self._expert_frames:
                if f not in hide_sim_modes:
                    f.show()
            for f in show_sim_modes:
                f.show()

        return

        
    
