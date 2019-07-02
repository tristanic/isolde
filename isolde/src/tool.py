# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 14-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



# ToolUI should inherit from ToolInstance if they will be
# registered with the tool state manager.
#
# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
from chimerax.core.tools import ToolInstance

def _find_help():
    # import os, pathlib
    # fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'doc', 'index.html')
    # return pathlib.Path(os.path.abspath(fname)).as_uri()
    return 'help:isolde/index.html'

class ISOLDE_ToolUI(ToolInstance):

    SESSION_ENDURING = False
    # if SESSION_ENDURING is True, tool instance not deleted at session closure

    help = _find_help()



    def __init__(self, session, tool_name):
        ToolInstance.__init__(self, session, tool_name)
        self.display_name = "ISOLDE"
        from chimerax.core.ui.gui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        self.tool_window.manage(placement=None)
        parent = self.tool_window.ui_area
        pp = parent.parent().parent()
        pp.resize(540,850)

        import os
        basedir = os.path.dirname(os.path.abspath(__file__))
        from PyQt5 import uic
        uifile = os.path.join(basedir, 'ui', 'IsoldeFrame.ui')
        from PyQt5 import QtWidgets, QtGui
        #from . import isoldewidget
        from .resources import resources_rc
        mw = self.mainwin = QtWidgets.QFrame(parent=parent)
        iw = self.iw = uic.loadUi(uifile, mw)
        # iw = self.iw = isoldewidget.Ui_isolde_widget()
        #iw = w.Ui_isolde_widget()
        # iw.setupUi(self.mainwin)

        import os
        icon_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources/logo_small.png')
        isolde_icon = QtGui.QPixmap(icon_file)
        iw._isolde_icon.setPixmap(isolde_icon)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.mainwin)
        layout.setStretchFactor(self.mainwin, 1)
        layout.setContentsMargins(1,1,1,3)
        parent.setLayout(layout)
        self.tool_window.manage(placement=None)
        # Should load saved state here

        ###
        # Selection manipulation buttons
        ###

        # Grow selection at N-terminus
        self._sel_grow_n_terminal_buttons = (
            iw._sel_extend_N_button,
            )

        # Shrink selection at N-terminus
        self._sel_shrink_n_terminal_buttons = (
            iw._sel_shrink_N_button,
            )

        # Shrink selection at C-terminus
        self._sel_shrink_c_terminal_buttons = (
            iw._sel_shrink_C_button,
            )

        # Grow selection at C-terminus
        self._sel_grow_c_terminal_buttons = (
            iw._sel_extend_C_button,
            )



        # Define intermediate and expert frames
        self._intermediate_frames = [
            iw._sim_platform_frame,
            ]
        self._expert_frames = [
            iw._force_field_selection_frame,
            ]

        # Any other frames/widgets that should be hidden at the start
        self._hidden_at_start = [
            iw._validate_rama_main_frame,
            iw._validate_pep_main_frame,
            iw._validate_rota_main_frame,
            iw._sim_basic_xtal_init_open_button,
            iw._sim_basic_xtal_init_main_frame,
            iw._sim_basic_xtal_settings_live_recalc_checkbox,
            iw._sim_basic_xtal_map_settings_frame,
            iw._real_space_map_from_volume_frame,
            iw._live_map_control_frame,
            iw._sim_status_indicator,
            ]

        for f in self._hidden_at_start:
            f.hide()

        # Any frames/widgets that should be disabled at the start
        self._disabled_at_start = [
            #iw._sim_basic_xtal_map_settings_frame,
            iw._map_masking_frame,
            iw._rebuild_sel_res_pep_flip_button,
            iw._rebuild_sel_res_cis_trans_flip_button,
            iw._rebuild_cycle_rotamer_frame,
            #iw._rebuild_sel_res_last_rotamer_button,
            #iw._rebuild_sel_res_next_rotamer_button,
            iw._rebuild_sel_res_rot_commit_button,
            iw._rebuild_sel_res_rot_target_button,
            iw._rebuild_sel_res_rot_discard_button,
            iw._rebuild_sel_res_rot_release_button,
            #~ iw._rebuild_pos_restraint_one_atom_frame,
            #iw._rebuild_pin_atom_container,
            iw._rebuild_pin_atom_to_current_pos_button,
            iw._rebuild_pin_atom_to_pivot_button,
            iw._rebuild_pos_restraint_clear_button,
            iw._rebuild_2ry_struct_restr_container,
            iw._rebuild_2ry_struct_restr_clear_button,
            iw._rebuild_register_shift_container,
            iw._rebuild_dist_restraint_container,
            iw._rebuild_dist_restraint_apply_button,
            iw._rebuild_grow_shrink_sel_frame,
            iw._right_mouse_modes_frame,
            iw._live_map_recalc_button,
            iw._live_map_update_sim_button,
            ]
        for f in self._disabled_at_start:
            f.setEnabled(False)



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



    def delete(self):
        try:
            self.isolde._on_close()
        except Exception as e:
            print(e)
            import traceback
            import sys
            traceback.print_stack(file=sys.__stderr__)
        super().delete()

    def _change_experience_level(self):
        iw = self.iw
        exp_index = iw._experience_level_combo_box.currentIndex()
        if (exp_index == 0):
            # Easy. Just hide everything intermediate or expert, and
            # everything belonging to other sim modes
            for f in self._intermediate_frames:
                f.hide()
            for f in self._expert_frames:
                f.hide()
        elif (exp_index == 1):
            for f in self._intermediate_frames:
                f.show()
            for f in self._expert_frames:
                f.hide()
        else:
            for f in self._intermediate_frames:
                f.show()
            for f in self._expert_frames:
                f.show()

        return
