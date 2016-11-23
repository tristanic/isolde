# vim: set expandtab shiftwidth=4 softtabstop=4:

# ToolUI should inherit from ToolInstance if they will be
# registered with the tool state manager.
#
# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
from chimerax.core.tools import ToolInstance


class ISOLDE_ToolUI(ToolInstance):

    SESSION_ENDURING = False
    # if SESSION_ENDURING is True, tool instance not deleted at session closure

    def __init__(self, session, tool_name):
        ToolInstance.__init__(self, session, tool_name)
        self.display_name = "ISOLDE"
        from chimerax.core.ui.gui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        self.tool_window.manage(placement=None)
        parent = self.tool_window.ui_area
        pp = parent.parent()
        pp.resize(480,700) 

        from PyQt5 import QtWidgets
        from . import isoldewidget
        self.mainwin = QtWidgets.QFrame(parent=parent)
        self.iw = isoldewidget.Ui_isolde_widget()
        self.iw.setupUi(self.mainwin)
        
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.mainwin)
        layout.setStretchFactor(self.mainwin, 1)
        parent.setLayout(layout)
        self.tool_window.manage(placement=None)
        # Should load saved state here
        
        # Define frames specific to crystallograpy, EM or free mode
        self._xtal_frames = [
            self.iw._sim_basic_xtal_map_frame,
            ]
        self._em_frames = [
            self.iw._sim_basic_em_map_frame,
            ]
        self._free_frames = [
            ]
        
        self._sim_mode_frame_lists = [
            self._xtal_frames,
            self._em_frames,
            self._free_frames,
            ]
        
        # Radio buttons to choose different selection modes
        self._selection_mode_buttons = (
            self.iw._sim_basic_by_selected_atoms_button,
            self.iw._sim_basic_by_chain_button,
            self.iw._sim_basic_whole_structure_button,
            self.iw._sim_basic_custom_selection_button,
            )        
        
        
        # Define intermediate and expert frames
        self._intermediate_frames = [
            self.iw._sim_platform_frame,
            ]
        self._expert_frames = [
            self.iw._force_field_selection_frame,
            self.iw._sim_basic_custom_selection_button,
            ]
        
        # Any other frames/widgets that should be hidden at the start
        self._hidden_at_start = [
            self.iw._validate_rama_main_frame,
            ]
        
        for f in self._hidden_at_start:
            f.hide()
                
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

        from . import isolde
        self.isolde = isolde.Isolde(self)


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

        
    
