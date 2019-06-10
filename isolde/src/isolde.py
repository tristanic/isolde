# @Author: Tristan Croll <tic20>
# @Date:   25-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 10-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os, sys
import numpy
from math import inf, degrees, radians, pi
from enum import IntEnum

import PyQt5
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QFileDialog

_openmm_initialized = False
def initialize_openmm():
    # On linux need to set environment variable to find plugins.
    # Without this it gives an error saying there is no "CPU" platform.
    global _openmm_initialized
    if not _openmm_initialized:
        _openmm_initialized = True
        from sys import platform
        if platform in ['linux', 'darwin']:
            from os import environ, path
            from chimerax import app_lib_dir
            environ['OPENMM_PLUGIN_DIR'] = path.join(app_lib_dir, 'plugins')
initialize_openmm()
from simtk import unit

import chimerax

from chimerax import clipper

from chimerax.core import triggerset
from chimerax.core.models import Drawing, Model
from chimerax.map import Volume
from chimerax.atomic import AtomicStructure

from .eventhandler import EventHandler
from .constants import defaults, sim_outcomes, control
from .param_mgr import Param_Mgr, autodoc, param_properties
from .checkpoint import CheckPoint
#from .openmm import sim_interface
from .openmm.sim_param_mgr import SimParams

from PyQt5.QtWidgets import QMessageBox

OPENMM_LENGTH_UNIT =        defaults.OPENMM_LENGTH_UNIT
OPENMM_FORCE_UNIT =         defaults.OPENMM_FORCE_UNIT
OPENMM_SPRING_UNIT =        defaults.OPENMM_SPRING_UNIT
OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT
OPENMM_ENERGY_UNIT =        defaults.OPENMM_ENERGY_UNIT
OPENMM_ANGLE_UNIT =         defaults.OPENMM_ANGLE_UNIT
OPENMM_TEMPERATURE_UNIT =   defaults.OPENMM_TEMPERATURE_UNIT
CHIMERAX_LENGTH_UNIT =      defaults.CHIMERAX_LENGTH_UNIT
CHIMERAX_FORCE_UNIT =       defaults.CHIMERAX_FORCE_UNIT
CHIMERAX_SPRING_UNIT =      defaults.CHIMERAX_SPRING_UNIT

@param_properties
@autodoc
class IsoldeParams(Param_Mgr):
    '''
    Contains the basic settings for defining and displaying the atoms and maps
    involved in an interactive simulation. Changes to parameter values will
    affect all future simulations in the ISOLDE session. At present, permanent
    changes can only be made by editing the corresponding entry in `constants.py`.

    '''
    _default_params = {
            # Number of residues before and after each selected residue to add
            # to the mobile selection
        'num_selection_padding_residues':       (defaults.SELECTION_SEQUENCE_PADDING, None),
            # Residues coming within this distance of the primary mobile selection
            # will also become mobile
        'soft_shell_cutoff_distance':           (defaults.SOFT_SHELL_CUTOFF, None),
            # Residues coming within this distance of any mobile atom will be
            # incorporated as fixed atoms in the simulation
        'hard_shell_cutoff_distance':           (defaults.HARD_SHELL_CUTOFF, None),
            # If true, residues in the soft shell will have their backbone atoms
            # fixed in space
        'fix_soft_shell_backbone':              (defaults.FIX_SOFT_SHELL_BACKBONE, None),
            # Multiplier to reduce the bond radii for bonds between fixed_atoms
        'fixed_bond_radius_ratio':              (defaults.FIXED_BOND_RADIUS_RATIO, None),
            # Limit the drawing to only the atoms involved in the simulation
            # (highly recommended for performance)
        'hide_surroundings_during_sim':         (defaults.HIDE_SURROUNDINGS_DURING_SIM, None),
            # Update map masking at regular intervals, to account for large
            # atom movements
        'remask_maps_during_sim':               (defaults.REMASK_MAPS_DURING_SIM, None),
        'rounds_per_map_remask':                (defaults.ROUNDS_PER_MAP_REMASK, None),
            # Mask cutoff for a 2mFo-DFc crystallographic or standard EM map
        'standard_map_mask_cutoff':             (defaults.STANDARD_MAP_MASK_RADIUS, None),
            # Mask cutoff for a difference map (e.g. mFo-DFc)
        'difference_map_mask_cutoff':           (defaults.DIFFERENCE_MAP_MASK_RADIUS, None),
            # Use multiprocessing instead of threads (only available on Linux)
        'map_shannon_rate':                     (defaults.MAP_SHANNON_RATE, None),

    }

class Isolde():
    '''
    The master ISOLDE class, primarily designed to be used by the ISOLDE GUI as
    the front-end interface for starting, controlling and analysing interactive
    simulations. Should only be used as a session-level singleton.
    '''
    ####
    # Environment information
    ####
    _root_dir = os.path.dirname(os.path.abspath(__file__))
    _platform = sys.platform

    ####
    # Enums for menu options
    ####

    _force_groups = {
        'main':                         0,
        'position restraints':          1,
        'map forces':                   2,
        }

    # Master switch to set the level of control the user has over the simulation.
    class _experience_levels(IntEnum):
        beginner        = 0
        intermediate    = 1
        expert          = 2

    class _sim_selection_modes(IntEnum):
        from_picked_atoms       = 0
        chain                   = 1
        whole_model             = 2
        custom                  = 3
        script                  = 99



    '''
    Named triggers to simplify handling of changes on key ISOLDE events,
    using the same code as ChimeraX's session.triggers. The same
    functionality could be achieved without too much trouble using the
    PyQt signal/slot system.
    '''
    trigger_names = (
        'selected model changed', # Changed the master model selection
        'simulation started',
        'simulation terminated',
        'simulation paused',
        'simulation resumed',
        'isolde closed'
        )

    def __init__(self, gui=None):
        '''
        Initialises the ISOLDE object and adds it to the ChimeraX session as
        ``session.isolde``.

        Args:
            * gui (default=None):
                - Used by the ISOLDE BundleAPI to prepare the ISOLDE GUI
                  interface (i.e. when running from Tools/General/ISOLDE in the
                  ChimeraX GUI menu).
        '''
        self.session = session = gui.session

        self.triggers = triggerset.TriggerSet()
        for t in self.trigger_names:
            self.triggers.add_trigger(t)

        self._isolde_events = EventHandler(self)
        self._logging = False
        self._log = Logger('isolde.log')
        self.checkpoint_disabled_reasons = {}

        '''
        A :class:`IsoldeParams` instance containing basic parameters controlling
        how the interactive simulation selection is defined and displayed.
        '''
        self.params = IsoldeParams()

        '''
        Parameters controlling the initialisation and behaviour of the
        simulation.
        '''
        sp = self.sim_params = SimParams()

        from .openmm.forcefields import Forcefield_Mgr
        ffmgr = self._ff_mgr = Forcefield_Mgr()

        self._status = self.session.logger.status

        self._event_handler = EventHandler(self.session)

        # Available pre-defined colors
        from chimerax.core import colors
        import copy

        self._available_models = {}

        self._available_colors = copy.copy(colors.BuiltinColors)
        # Remove duplicates
        for key in self._available_colors:
            if ' ' in key:
                stripped_key = key.replace(' ',"")
                self._available_colors.pop(stripped_key, None)

        ####
        # Settings for handling of atomic coordinates
        ####

        # Currently chosen mode for selecting the mobile simulation
        self._sim_selection_mode = None
        # Selected model on which we are actually going to run a simulation
        self._selected_model = None

        ####
        # Settings for live tracking of structural quality
        ####
        # Generic widget object holding the Ramachandran plot
        self._rama_plot_window = None
        # Object holding Ramachandran plot information and controls
        self._rama_plot = None

        ####
        # Settings for handling of map visualisation
        ####
        # For updating map zoning to current coordinates
        self._map_rezone_counter = 0

        ####
        # Mouse modes
        ####
        from . import mousemodes

        # Object to initialise and hold ISOLDE-specific mouse modes,
        # and keep track of standard modes they've replaced. The pre-existing
        # mode will be reinstated when a mode is removed.
        # TODO: Replace with system that plays more nicely with ChimeraX built-in
        #       modes
        mm = self._mouse_modes = mousemodes.MouseModeRegistry(self.session, self)

        # Placeholder for mouse tugging object
        self._mouse_tugger = None

        ####
        # Variables for live information/control
        ####

        self._selected_rotamer = None


        ####
        # Haptic interfaces
        ####
        # TODO: Move out into separate class
        # Will be set to true if haptic devices are detected
        self._use_haptics = False
        # Number of connected haptic devices
        self._num_haptic_devices = None
        # Array of HapticTugger objects, one for each device
        self._haptic_devices = []
        # Atoms currently being tugged by haptic devices
        self._haptic_tug_atom = []
        # Highlight the nearest atom?
        self._haptic_highlight_nearest_atom = []
        # Current nearest atom
        self._haptic_current_nearest_atom = []
        # Nearest atom's normal color
        self._haptic_current_nearest_atom_color = []
        # Spring constant felt by atoms tugged by the haptic device
        self._haptic_force_constant = defaults.HAPTIC_SPRING_CONSTANT

        # Handler for shifting stretches in register.
        self._register_shifter = None

        ####
        # Settings for OpenMM
        ####

        # Placeholder for the openmm_interface.Sim_Handler object used to perform
        # simulation setup and management
        self._sim_manager = None

        # self.initialize_haptics()



        ####
        # Internal trigger handlers
        ####

        self.gui_mode = False

        self._ui_panels = []

        if gui is not None:
            from PyQt5.QtGui import QPixmap
            from PyQt5.QtWidgets import QSplashScreen
            from PyQt5.QtCore import Qt
            import os

            splash_pix = QPixmap(os.path.join(
                self._root_dir,'resources/isolde_splash_screen.jpg'))
            splash = self._splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
            splash.setMask(splash_pix.mask())
            splash.show()
            # Make sure the splash screen is actually shown
            from time import sleep
            for i in range(5):
                self.session.ui.processEvents()
                sleep(0.01)
            from PyQt5 import QtCore
            # Close the splash after 5 seconds
            QtCore.QTimer.singleShot(5000, splash.close)

            self._start_gui(gui)

        session.isolde = self
        ffmgr.background_load_ff(sp.forcefield)


    def _prepare_environment(self):
        session = self.session
        from chimerax.std_commands import cofr, camera
        cofr.cofr(session, 'centerOfView', show_pivot=True)
        camera.camera(session, 'ortho')
        self._mouse_modes.register_all_isolde_modes()

    @property
    def forcefield_mgr(self):
        return self._ff_mgr

    @property
    def ignored_residues(self):
        '''
        A :class:`chimerax.Residues` object defining residues to be ignored for
        simulation purposes. Residues in this list will be retained in the model
        and used for map calculations, but will not be included in simulations.
        In most cases, this should only be populated with residues that have not
        been parameterised for MD.
        '''
        m = self.selected_model
        if m is None:
            return None
        if not hasattr(m, '_isolde_ignored_residues') or m._isolde_ignored_residues is None:
            from chimerax.atomic import Residues
            m._isolde_ignored_residues = Residues()
        return m._isolde_ignored_residues

    @ignored_residues.setter
    def ignored_residues(self, residues):
        self.selected_model._isolde_ignored_residues = residues

    @property
    def sim_handler(self):
        '''
        Returns the :class:`isolde.openmm.Sim_Handler` instance controlling
        the currently running simulation, or None if no simulation is running.
        Read only.
        '''
        if self.sim_manager is None:
            return None
        if not hasattr(self.sim_manager, 'sim_handler'):
            return None
        return self.sim_manager.sim_handler

    @property
    def simulation_running(self):
        '''
        Returns True if a simulation is running, otherwise False. Read only.
        '''
        if self.sim_handler is None:
            return False
        return self.sim_handler.sim_running

    @property
    def sim_manager(self):
        '''
        Returns the :class:`isolde.openmm.Sim_Manager` instance providing
        high-level control over the current simulation. If no simulation is
        currently defined, returns None. Read only.
        '''
        return self._sim_manager

    @property
    def simulation_mode(self):
        '''
        Returns 'equil' or 'min' depending on the state of the simulation,
        or None if no simulation is defined. Read only. To switch modes, use
        :func:`isolde.minimize` and :func@`isolde.equilibrate`.
        '''
        sh = self.sim_handler
        if sh is None:
            return None
        if not sh.minimize:
            return 'equil'
        return 'min'

    @property
    def selected_model(self):
        '''
        Returns the atomic model on which ISOLDE is currently operating.
        Can be set.
        '''
        return self._selected_model

    @selected_model.setter
    def selected_model(self, model):
        from chimerax.atomic import AtomicStructure
        if not isinstance(model, AtomicStructure):
            raise TypeError('Selection must be a single AtomicStructure model!')
        self.change_selected_model(model)


    ###################################################################
    # GUI related functions
    ###################################################################

    def _start_gui(self, gui):
        '''
        Connect and initialise ISOLDE widget
        '''

        self.gui = gui
        iw = self.iw = gui.iw

        # FIXME Remove 'Custom restraints' tab from the gui while I decide
        # whether to keep it
        iw._sim_tab_widget.removeTab(1)

        self.gui_mode = True

        # Prepare display and register ISOLDE-specific mouse modes
        self._prepare_environment()



        # Function to remove all event handlers, mouse modes etc. when
        # ISOLDE is closed, and return ChimeraX to its standard state.

        self._event_handler.add_event_handler('close on app quit',
                                      'app quit',
                                      self._on_close)


        # Any values in the Qt Designer .ui file are placeholders only.
        # Combo boxes need to be repopulated, and variables need to be
        # set to match any saved state.
        self._populate_menus_and_update_params()

        # Make sure everything in the widget actually does something.
        self._connect_functions()

        ####
        # Add handlers for GUI events, and run each callback once to
        # initialise to current conditions
        ####

        from chimerax.core.models import ADD_MODELS, MODEL_ID_CHANGED, REMOVE_MODELS
        eh = self._event_handler
        eh.add_event_handler('update_menu_on_selection',
                             'selection changed', self._selection_changed)
        eh.add_event_handler('update_menu_on_model_add',
                            ADD_MODELS, self._update_model_list)
        eh.add_event_handler('update_menu_on_model_id_change',
                            MODEL_ID_CHANGED, self._update_model_list)
        eh.add_event_handler('update_menu_on_model_remove',
                             REMOVE_MODELS, self._update_model_list)
        self._selection_changed()
        self._update_model_list(None, None, force=True)

        # Work out menu state based on current ChimeraX session
        self._update_sim_control_button_states()

        st = iw._sim_tab_widget
        st.setCurrentIndex(0)


    def _populate_menus_and_update_params(self):
        iw = self.iw
        params = self.params
        sim_params = self.sim_params
        # Clear experience mode placeholders from QT Designer and repopulate
        cb = iw._experience_level_combo_box
        cb.clear()
        for lvl in self._experience_levels:
            cb.addItem(lvl.name, lvl)

        iw._sim_temp_spin_box.setValue(self.sim_params.temperature.value_in_unit(OPENMM_TEMPERATURE_UNIT))

        # Populate force field combo box with available forcefields
        cb = iw._sim_force_field_combo_box
        cb.clear()
        cb.addItems(self.forcefield_mgr.available_forcefields)
        cb.setCurrentIndex(cb.findText(self.sim_params.forcefield))

        # Populate OpenMM platform combo box with available platforms
        cb = iw._sim_platform_combo_box
        cb.clear()

        from .openmm.openmm_interface import get_available_platforms
        platform_names = get_available_platforms()
        cb.addItems(platform_names)

        # Set to the preferred or, failing that, the fastest available platform
        if sim_params.platform in platform_names:
            cb.setCurrentIndex(cb.findText(sim_params.platform))
        else:
            for p in sim_params.platforms:
                if p in platform_names:
                    cb.setCurrentIndex(cb.findText(p))
                    sim_params.platform = p
                    break


        # Basic settings for defining the mobile selection
        iw._sim_basic_mobile_b_and_a_spinbox.setValue(params.num_selection_padding_residues)
        iw._sim_basic_mobile_sel_within_spinbox.setValue(params.soft_shell_cutoff_distance)

        iw._sim_basic_xtal_settings_mask_radius_spinbox.setValue(params.standard_map_mask_cutoff)
        # iw._sim_basic_xtal_map_weight_spin_box.setValue(params.difference_map_mask_cutoff)
        ####
        # Rebuild tab
        ####

        ## Info for a single selected residue
        iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        iw._rebuild_sel_res_pep_info.setText('')
        iw._rebuild_sel_res_rot_info.setText('')

        iw._rebuild_pos_restraint_spring_constant.setValue(
            self.sim_params.position_restraint_spring_constant.value_in_unit(CHIMERAX_SPRING_UNIT))
        iw._rebuild_dist_restraint_spring_constant.setValue(
            self.sim_params.distance_restraint_spring_constant.value_in_unit(CHIMERAX_SPRING_UNIT))
        ####
        # Validate tab
        ####

        self._update_model_list(None, None, force=True)

    def _connect_functions(self):
        '''
        Connect PyQt events from the ISOLDE gui widget to functions.
        '''
        iw = self.iw
        ####
        # Master switches
        ####

        iw._experience_level_combo_box.currentIndexChanged.connect(
            self.gui._change_experience_level
            )
        self.gui._change_experience_level()


        iw._demo_load_button.clicked.connect(self.load_demo_data)

        ####
        # Help
        ####

        iw._main_help_button.clicked.connect(
            self.show_master_help_in_browser
        )

        ####
        # Right mouse button modes
        ####
        iw._right_mouse_tug_atom_button.clicked.connect(
            self._set_right_mouse_mode_tug_atom
        )
        iw._right_mouse_tug_residue_button.clicked.connect(
            self._set_right_mouse_mode_tug_residue
        )
        iw._right_mouse_tug_selection_button.clicked.connect(
            self._set_right_mouse_mode_tug_selection
        )
        ####
        # Simulation global parameters (can only be set before running)
        ####
        iw._sim_force_field_combo_box.currentIndexChanged.connect(
            self._change_force_field
            )
        iw._master_model_combo_box.currentIndexChanged.connect(
            self._change_selected_model
            )
        iw._sim_basic_mobile_sel_within_spinbox.valueChanged.connect(
            self._change_soft_shell_cutoff_from_sel_menu
            )
        iw._sim_basic_mobile_b_and_a_spinbox.valueChanged.connect(
            self._change_b_and_a_padding
            )
        iw._sim_platform_combo_box.currentIndexChanged.connect(
            self._change_sim_platform
            )
        iw._sim_basic_xtal_settings_mask_radius_spinbox.valueChanged.connect(
            self._change_mask_radius
        )
        iw._sim_basic_xtal_settings_spotlight_radius_spinbox.valueChanged.connect(
            self._change_spotlight_radius
        )
        # Run all connected functions once to initialise
        self._change_force_field()
        self._change_selected_model(force=True)
        self._change_b_and_a_padding()
        self._change_sim_platform()


        iw._save_model_button.clicked.connect(
            self._save_cif_file
            )

        ####
        # Buttons to grow/shrink a continuous selection
        ####
        for b in self.gui._sel_grow_n_terminal_buttons:
            b.clicked.connect(self._extend_selection_by_one_res_N)

        for b in self.gui._sel_shrink_n_terminal_buttons:
            b.clicked.connect(self._shrink_selection_by_one_res_N)

        for b in self.gui._sel_shrink_c_terminal_buttons:
            b.clicked.connect(self._shrink_selection_by_one_res_C)

        for b in self.gui._sel_grow_c_terminal_buttons:
            b.clicked.connect(self._extend_selection_by_one_res_C)

        ####
        # Real space map parameters (can only be set before starting simulation)
        ####

        iw._real_space_map_from_volume_show_button.clicked.connect(
            self._show_real_space_map_dialog
        )
        iw._real_space_map_from_volume_done_button.clicked.connect(
            self._hide_real_space_map_dialog
        )
        iw._real_space_map_from_volume_button.clicked.connect(
            self._add_real_space_map_from_gui
        )

        ####
        # Trajectory smoothing
        ####

        iw._trajectory_smooth_button.clicked.connect(
            self._change_smoothing_state_from_gui
        )
        iw._smoothing_amount_dial.valueChanged.connect(
            self._change_smoothing_amount_from_gui
        )


        ####
        # Xtal map parameters (can only be set before starting simulation)
        ####
        iw._sim_basic_xtal_init_open_button.clicked.connect(
            self._show_xtal_init_frame
        )
        iw._sim_basic_xtal_init_done_button.clicked.connect(
            self._hide_xtal_init_frame
        )
        iw._sim_basic_xtal_init_reflections_file_button.clicked.connect(
            self._choose_mtz_file
        )
        iw._sim_basic_xtal_init_go_button.clicked.connect(
            self._initialize_xtal_structure
        )
        iw._sim_basic_xtal_map_weight_spin_box.valueChanged.connect(
            self._update_map_weight_box_settings
        )

        # Load custom forcefield file(s)
        iw._sim_basic_ff_file_load_button.clicked.connect(
            self._add_ff_files_gui
        )

        # Live maps direct from structure factors
        iw._sim_basic_xtal_init_exp_data_button.clicked.connect(
            self._choose_reflections_file
        )

        iw._sim_basic_xtal_map_settings_show_button.clicked.connect(
            self._toggle_xtal_map_dialog
            )
        iw._sim_basic_xtal_settings_live_recalc_checkbox.stateChanged.connect(
            self._set_live_xmap_recalc
        )
        iw._sim_basic_xtal_settings_map_combo_box.currentIndexChanged.connect(
            self._populate_xtal_map_params
            )
        iw._sim_basic_xtal_settings_enable_mdff_checkbox.clicked.connect(
            self._enable_or_disable_mdff_potential
        )
        iw._sim_basic_xtal_settings_set_button.clicked.connect(
            self._apply_xtal_map_params
            )

        # Visualisation tools
        iw._map_settings_solid_surface_button.clicked.connect(
            self._set_map_to_solid_surface
        )
        iw._map_settings_transparent_surface_button.clicked.connect(
            self._set_map_to_transparent_surface
        )
        iw._map_settings_mesh_button.clicked.connect(
            self._set_map_to_mesh
        )
        iw._map_settings_color_button.clicked.connect(
            self._choose_map_color
        )

        iw._vis_step_mask_forward_button.clicked.connect(
            self._xtal_step_forward
            )
        iw._vis_step_mask_backward_button.clicked.connect(
            self._xtal_step_backward
            )
        iw._vis_mask_to_sel_button.clicked.connect(
            self._xtal_mask_to_selection
            )
        iw._vis_spotlight_mode_button.clicked.connect(
            self._xtal_enable_live_scrolling
            )


        ####
        # Rebuild tab
        ####
        iw._rebuild_sel_res_cis_trans_flip_button.clicked.connect(
            self._flip_cis_trans
            )
        iw._rebuild_sel_res_pep_flip_button.clicked.connect(
            self._flip_peptide_bond
            )
        iw._rebuild_sel_res_last_rotamer_button.clicked.connect(
            self._prev_rotamer
            )
        iw._rebuild_sel_res_next_rotamer_button.clicked.connect(
            self._next_rotamer
            )
        iw._rebuild_sel_res_rot_commit_button.clicked.connect(
            self._commit_rotamer
            )
        iw._rebuild_sel_res_rot_target_button.clicked.connect(
            self._set_rotamer_target
            )
        iw._rebuild_sel_res_rot_discard_button.clicked.connect(
            self._clear_rotamer
            )
        iw._rebuild_sel_res_rot_release_button.clicked.connect(
            self._release_rotamer
            )
        iw._rebuild_sel_res_rot_backrub_button.clicked.connect(
            self._backrub_rotamer
        )

        iw._rebuild_restrain_helix_button.clicked.connect(
            self._restrain_selection_as_alpha_helix
            )
        iw._rebuild_restrain_anti_beta_button.clicked.connect(
            self._restrain_selection_as_antiparallel_beta
            )
        iw._rebuild_restrain_par_beta_button.clicked.connect(
            self._restrain_selection_as_parallel_beta
            )
        iw._rebuild_2ry_struct_restr_clear_button.clicked.connect(
            self.clear_secondary_structure_restraints_for_selection
            )
        iw._rebuild_register_shift_reduce_button.clicked.connect(
            self._decrement_register_shift
            )
        iw._rebuild_register_shift_increase_button.clicked.connect(
            self._increment_register_shift
            )
        iw._rebuild_register_shift_go_button.clicked.connect(
            self._apply_register_shift
            )
        iw._rebuild_register_shift_release_button.clicked.connect(
            self._release_register_shifter
            )
        iw._rebuild_pin_atom_to_current_pos_button.clicked.connect(
            self._restrain_selected_atom_to_current_xyz
            )
        iw._rebuild_pin_atom_to_pivot_button.clicked.connect(
            self._restrain_selected_atom_to_pivot_xyz
            )
        iw._rebuild_pos_restraint_clear_button.clicked.connect(
            self.release_xyz_restraints_on_selected_atoms
            )
        iw._rebuild_dist_restraint_apply_button.clicked.connect(
            self._add_distance_restraint_between_selected_atoms
        )
        iw._rebuild_dist_restraint_set_target_to_current_distance_button.clicked.connect(
            self._set_distance_restraint_target_to_current_distance
        )
        iw._rebuild_remove_distance_restraint_button.clicked.connect(
            self._release_distance_restraint_between_selected_atoms
        )
        iw._rebuild_remove_all_distance_restraints_button.clicked.connect(
            self._release_all_distance_restraints_for_selected_atoms
        )

        ####
        # Validation tab
        ####

        iw._validate_rama_show_button.clicked.connect(
            self._show_rama_plot
            )
        iw._validate_rama_hide_button.clicked.connect(
            self._hide_rama_plot
            )

        iw._validate_pep_show_button.clicked.connect(
            self._show_peptide_validation_frame
            )
        iw._validate_pep_hide_button.clicked.connect(
            self._hide_peptide_validation_frame
            )
        iw._validate_pep_update_button.clicked.connect(
            self._update_iffy_peptide_lists
            )
        iw._validate_pep_iffy_table.itemClicked.connect(
            self._show_selected_iffy_peptide
            )

        iw._validate_rota_show_button.clicked.connect(
            self._show_rota_validation_frame
            )
        iw._validate_rota_hide_button.clicked.connect(
            self._hide_rota_validation_frame
            )
        iw._validate_rota_update_button.clicked.connect(
            self._update_iffy_rota_list
            )
        iw._validate_rota_table.itemClicked.connect(
            self._show_selected_iffy_rota
        )

        from .validation.clashes import Clash_Table_Mgr
        self._clash_mgr = Clash_Table_Mgr(self)

        ####
        # Simulation control functions
        ####

        iw._sim_temp_spin_box.valueChanged.connect(
            self._update_sim_temperature
            )
        iw._sim_go_button.clicked.connect(
            self._start_sim_or_toggle_pause
            )
        iw._sim_save_checkpoint_button.clicked.connect(
            self.checkpoint)
        iw._sim_revert_to_checkpoint_button.clicked.connect(
            self.revert_to_checkpoint)
        iw._sim_commit_button.clicked.connect(
            self.commit_sim
            )
        iw._sim_stop_and_revert_to_checkpoint_button.clicked.connect(
            self._stop_sim_and_revert_to_checkpoint
            )
        iw._sim_stop_and_discard_button.clicked.connect(
            self._discard_sim
            )

        iw._sim_min_button.clicked.connect(
            self.minimize
            )
        iw._sim_equil_button.clicked.connect(
            self.equilibrate
        )

    def _disable_chimerax_mouse_mode_panel(self, *_):
        self._set_chimerax_mouse_mode_panel_enabled(False)

    def _enable_chimerax_mouse_mode_panel(self, *_):
        self._set_chimerax_mouse_mode_panel_enabled(True)

    def _set_chimerax_mouse_mode_panel_enabled(self, state):
        from chimerax.mouse_modes.tool import MouseModePanel
        mm = self.session.tools.find_by_class(MouseModePanel)[0]
        mm.buttons.setEnabled(state)

    # def initialize_haptics(self):
    #     '''
    #     If the HapticHandler plugin is installed, start it and create
    #     a HapticTugger object for each connected haptic device.
    #     '''
    #     if hasattr(self.session, 'HapticHandler'):
    #         self._status('Initialising haptic interface(s)')
    #         hh = self.session.HapticHandler
    #         hh.startHaptics()
    #         n = self._num_haptic_devices = hh.getNumDevices()
    #         if n == 0:
    #             print('No haptic devices found. Stopping haptic driver.')
    #             hh.stopHaptics()
    #             self._use_haptics = False
    #             return
    #
    #         d = self._haptic_devices = [None] * n
    #         self._haptic_tug_atom = [None] * n
    #         self._haptic_tug_index = [-1] * n
    #         self._haptic_highlight_nearest_atom = [True] * n
    #         self._haptic_current_nearest_atom = [None] * n
    #         self._haptic_current_nearest_atom_color = [None] * n
    #         self._haptic_current_nearest_atom_draw_mode = [None] * n
    #         from . import tugging
    #         for i in range(n):
    #             d[i] = tugging.HapticTugger(self.session, i)
    #         self._use_haptics = True
    #         self._status('')
    #         self._isolde_events.add_event_handler('Enable haptics during simulation',
    #                                                'simulation started',
    #                                                self._start_sim_haptics)
    #         self._isolde_events.add_event_handler('Disable haptics on sim end',
    #                                               'simulation terminated',
    #                                               self._stop_sim_haptics)
    #     else:
    #         self._use_haptics = False

    def _update_sim_temperature(self):
        t = self.iw._sim_temp_spin_box.value()
        if self.simulation_running:
            self.sim_handler.temperature = t
        self.sim_params.temperature = t

    ##############################################################
    # Menu control functions to run on key events
    ##############################################################
    # TODO: This function has become somewhat monolithic and has expanded well
    #       beyond its original remit. Needs some rethinking.
    def _update_model_list(self, trigger_name, models, force=False):
        from chimerax.core.models import Model
        if isinstance(models, Model):
            models = [models]
        from chimerax.atomic import AtomicStructure
        from chimerax.map import Volume
        from chimerax.clipper.symmetry import Symmetry_Manager
        from chimerax.clipper.maps import XmapHandler_Live
        from .restraints import MDFF_Mgr

        if force:
            structures_need_update=True
            volumes_need_update=True
        else:
            needs_update = False
            structures_need_update = False
            volumes_need_update = False

            for m in models:
                if type(m) == AtomicStructure:
                    structures_need_update = True
                if isinstance(m, Volume) or isinstance(m, MDFF_Mgr):
                    volumes_need_update = True
                if structures_need_update and volumes_need_update:
                    break
        mmcb = self.iw._master_model_combo_box
        current_model = mmcb.currentData()

        if structures_need_update:
            mmcb = self.iw._master_model_combo_box
            current_model = mmcb.currentData()
            mmcb.clear()

            _models = self.session.models.list()
            models = sorted(_models, key = lambda m: m.id)
            mtd = {
                AtomicStructure: [],
                Volume: []
            }

            for m in models:
                for mtype in mtd.keys():
                    if type(m) == mtype:
                        mtd[mtype].append(m)


            valid_models = mtd[AtomicStructure]
            valid_models = sorted(valid_models, key=lambda m: m.id)

            self._available_models = {}

            for m in valid_models:
                id_str = '{}. {}'.format(m.id_string, m.name)
                mmcb.addItem(id_str, m)
                self._available_models[id_str] = m

        if volumes_need_update:
            self._populate_available_volumes_combo_box()

        self._change_selected_model(model=current_model)
        self._initialize_maps(current_model)
        for p in self._ui_panels:
            p.chimerax_models_changed(self.selected_model)
        self._update_sim_control_button_states()

    def _selection_changed(self, *_):
        if not self.gui_mode:
            return
        from chimerax.atomic import selected_atoms
        from .util import is_continuous_protein_chain
        sel = selected_atoms(self.session)
        selres = sel.unique_residues
        if self.selected_model is not None: # self.simulation_running:
            natoms = len(sel)
            if natoms == 1:
                self._enable_atom_position_restraints_frame()
            else:
                self._disable_atom_position_restraints_frame()
            if natoms:
                self._enable_position_restraints_clear_button()
                self._enable_secondary_structure_restraints_clear_button()
                self._enable_distance_restraints_frame()
            else:
                self._disable_position_restraints_clear_button()
                self._disable_secondary_structure_restraints_clear_button()
                self._disable_distance_restraints_frame()

            if natoms == 2:
                self._enable_distance_restraint_apply_button()
                self._enable_distance_restraint_set_distance_to_current_button()
                self._enable_distance_restraint_remove_single_button()
            else:
                self._disable_distance_restraint_apply_button()
                self._disable_distance_restraint_set_distance_to_current_button()
                self._disable_distance_restraint_remove_single_button()

            if len(selres) == 1:
                r = selres[0]
                self._enable_rebuild_residue_frame(r)
            else:
                self._clear_rotamer()
                self._disable_rebuild_residue_frame()
            if is_continuous_protein_chain(sel):
                self._enable_secondary_structure_restraints_frame()
                self._enable_register_shift_frame()
            else:
                self._disable_secondary_structure_restraints_frame()
                self._disable_register_shift_frame()
            if is_continuous_protein_chain(sel, allow_single=True):
                self._enable_selection_extend_frame()
            else:
                self._disable_selection_extend_frame()

            # A running simulation takes precedence for memory control
            #return
        iw = self.iw
        if not self.simulation_running:
            self._disable_register_shift_frame()
            self._disable_peptide_bond_manipulation_frame()
            flag = not(self.session.selection.empty())
            iw._sim_go_button.setEnabled(flag)
        else:
            self._enable_peptide_bond_manipulation_frame()
            iw._sim_go_button.setEnabled(True)


    def _update_sim_control_button_states(self, *_):
        # Set enabled/disabled states of main simulation control panel
        # based on whether a simulation is currently running
        if self.gui_mode:
            running = self.simulation_running
            iw = self.iw
            paused = self.sim_paused
            go_button = iw._sim_go_button
            mdff_cb = iw._sim_basic_xtal_settings_enable_mdff_checkbox
            mdff_cb.setEnabled(not running)
            if paused and not running:
                go_button.setChecked(False)
            elif paused:
                go_button.setChecked(False)
                go_button.setToolTip('Resume')
            elif running:
                go_button.setChecked(True)
                go_button.setToolTip('Pause')
            if not running:
                go_button.setChecked(False)
                go_button.setToolTip('Start a simulation')

            iw._map_masking_frame.setDisabled(
                running
                    or not self.selected_model_has_maps
                    or self.selected_model is None)
            iw._right_mouse_modes_frame.setEnabled(running)
            iw._sim_save_checkpoint_button.setEnabled(running)
            iw._sim_revert_to_checkpoint_button.setEnabled(running)
            iw._sim_commit_button.setEnabled(running)
            iw._sim_stop_and_revert_to_checkpoint_button.setEnabled(running)
            iw._sim_stop_and_discard_button.setEnabled(running)
            if self.simulation_mode == 'equil':
                iw._sim_equil_button.setChecked(True)
            else:
                iw._sim_min_button.setChecked(True)

        # Update the status of the Go button
        self._selection_changed()

    def _model_changes_cb(self, trigger_name, changes):
        if changes is not None:
            changes = changes[1]
            added = len(changes.created_atoms())
            deleted = changes.num_deleted_atoms()
            if added or deleted:
                self._update_iffy_rota_list()
                self._update_iffy_peptide_lists()

    ####
    # General
    ####

    def _save_cif_file(self, *_):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        caption = 'Save structure as...'
        filetypes = 'mmCIF files (*.cif)'
        filename, _ = QFileDialog.getSaveFileName(None, caption, '', filetypes, options = options)
        if filename:
            self.save_cif_file(self._selected_model, filename)

    def save_cif_file(self, model, filename):
        '''
        Save a model in mmCIF format.

        Args:
            * model:
                - an :class:`AtomicStructure`
            * filename:
                - a text string providing a valid filename (must end with .cif).
        '''
        from chimerax.core.commands import save
        save.save(self.session, filename, [model])

    ####
    # Xtal
    ####
    def _show_xtal_init_frame(self, *_):
        self.iw._sim_basic_xtal_init_main_frame.show()
        self.iw._sim_basic_xtal_init_open_button.setEnabled(False)

    def _hide_xtal_init_frame(self, *_):
        self.iw._sim_basic_xtal_init_main_frame.hide()
        self.iw._sim_basic_xtal_init_open_button.setEnabled(True)

    def _check_for_valid_xtal_init(self, *_):
        sm = self.selected_model
        valid = True
        if sm is not None:
            if not os.path.isfile(self.iw._sim_basic_xtal_init_reflections_file_name.text()):
                valid = False
            if valid:
                from chimerax.clipper.symmetry import Symmetry_Manager
                for m in sm.all_models():
                    if isinstance(m, Symmetry_Manager):
                        valid = False
                        break
        else:
            valid = False
        if valid:
            self.iw._sim_basic_xtal_init_go_button.setEnabled(True)
        else:
            self.iw._sim_basic_xtal_init_go_button.setEnabled(False)
            self.iw._sim_basic_xtal_init_reflections_file_name.setText('')

    def _add_ff_files_gui(self, *_):
        files = self._choose_ff_files(self)
        if files is not None and len(files):
            ff = self.forcefield_mgr[self.sim_params.forcefield]
            self.add_ffxml_files(ff, files)

    def _choose_ff_files(self, *_):
        options = QFileDialog.Options()
        caption = 'Choose one or more ffXML files'
        filetypes = 'ffXML files (*.xml)'
        dlg = QFileDialog(caption=caption)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilter(filetypes)
        dlg.setFileMode(QFileDialog.ExistingFiles)
        import os
        dlg.setDirectory(os.getcwd())
        if dlg.exec():
            return dlg.selectedFiles()

    def add_ffxml_files(self, forcefield, file_list):
        log = self.session.logger
        for f in file_list:
            try:
                forcefield.loadFile(f)
            except ValueError as e:
                log.warning('Failed to add {}: {}'.format(f, str(e)))
            except:
                raise

    def _choose_mtz_file(self, *_):
        options = QFileDialog.Options()
        caption = 'Choose a file containing map structure factors'
        filetypes = 'MTZ files (*.mtz)'
        dlg = QFileDialog(caption=caption)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilter(filetypes)
        dlg.setFileMode(QFileDialog.ExistingFile)
        import os
        dlg.setDirectory(os.getcwd())

        #filename, _ = QFileDialog.getOpenFileName(None, caption, '', filetypes, options = options)
        if dlg.exec():
            return dlg.selectedFiles()[0]
        # self.iw._sim_basic_xtal_init_reflections_file_name.setText(filename)
        # self._check_for_valid_xtal_init()

    def _initialize_xtal_structure(self, *_):
        fname = self.iw._sim_basic_xtal_init_reflections_file_name.text()
        if not os.path.isfile(fname):
            from .dialog import generic_warning
            errstring = 'Please select a valid MTZ file!'
            generic_warning(errstring)
        m = self.selected_model
        spotlight_radius = self.iw._sim_basic_xtal_settings_spotlight_radius_spinbox.value()
        from chimerax.clipper import symmetry
        sym_handler = symmetry.Symmetry_Manager(m, fname,
            map_oversampling=self.params.map_shannon_rate,
            spotlight_radius=spotlight_radius)
        self.iw._sim_basic_xtal_init_reflections_file_name.setText('')
        self.iw._sim_basic_xtal_init_go_button.setEnabled(False)
        self._update_sim_control_button_states()

    def _choose_reflections_file(self, *_):
        # options = QFileDialog.Options()
        # caption = ('Choose a file containing amplitudes/phases or F/sigF/free flags')
        # filetypes = 'MTZ files (*.mtz)'
        # filename, _ = QFileDialog.getOpenFileName(None, caption, filetypes, filetypes, options=options)
        filename = self._choose_mtz_file()
        #import os
        #if os.path.isfile(filename):
        if filename is not None:
            self.add_xtal_data(filename)

    def add_xtal_data(self, filename, model=None):
        '''
        Add a set of maps derived from a crystal dataset (in MTZ format) to the
        target model. Cell dimensions and space group must match. If the file
        contains experimental measurements (F/sigF), a set of "live" maps will
        be generated and will be recalculated whenever the model coordinates
        change. If it contains precalculated map data (F/phi), one static map
        will be generated for each F/phi pair found.

        Args:
            * filename:
                - a file in .mtz format
            * model (default None):
                - A :class:`chimerax.AtomicStructure` or None. If None, the
                  model currently selected by ISOLDE will be used as the target.
        '''
        if model is None:
            model = self.selected_model
        if model is None:
            from .dialog import generic_warning
            generic_warning("You must have the corresponding model loaded "
                "before loading reflection data!")
        m = model
        from chimerax.clipper.symmetry import get_map_mgr
        map_mgr = get_map_mgr(m, create=True, auto_add_to_session=True)
        xmapset = map_mgr.add_xmapset_from_mtz(filename, self.params.map_shannon_rate)
        if xmapset.live_xmap_mgr is not None:
            xmapset.live_update = self.iw._sim_basic_xtal_settings_live_recalc_checkbox.checkState()
            # 2mFo-DFc and mFo-DFc maps are created automatically, but should
            # not be used as MDFF potentials. For that, we need a specific map
            # that we know excludes the free reflections.
            # self._add_mdff_potential_to_live_xmapset(xmapset)
        map_mgr.spotlight_radius = self.iw._sim_basic_xtal_settings_spotlight_radius_spinbox.value()
        self._change_selected_model(model=m, force=True)

    def _add_mdff_potential_to_live_xmapset(self, xmapset):
        from chimerax.clipper.maps.xmapset import map_potential_recommended_bsharp
        mdff_b = map_potential_recommended_bsharp(xmapset.resolution)
        mdff_p = xmapset.add_live_xmap('MDFF potential', b_sharp=mdff_b,
            is_difference_map=False,
            exclude_free_reflections=True,
            fill_with_fcalc = True,
            exclude_missing_reflections=True,
            display=False)
        return mdff_p


    def add_real_space_map(self,
            existing_volume=None, filename=None, to_model=None):
        '''
        Associate a real-space map with the given model. The input map will be
        "slaved" to the target model to behave similarly to a crystallographic
        map. Exactly one of `existing_volume` or `filename` should be given.

        Args:
            * existing_volume (default None):
                - A :class:`chimerax.Volume` object or None. If not None, the
                  chosen :class:`Volume` will be removed from the model tree,
                  and replaced with a :class:`clipper.maps.NXmapHandler` object.
            * filename (default None):
                - A map file of any format recognised by ChimeraX.
            * to_model (default None):
                - A :class:`chimerax.AtomicStructure` or None. If None, the
                  model currently selected by ISOLDE will be used as the target.
        '''
        if to_model is None:
            to_model = self.selected_model
        if to_model is None:
            from .dialog import generic_warning
            generic_warning("You must have the corresponding model loaded "
                "before loading reflection data!")
        m = to_model
        if existing_volume is None and filename is None:
            raise TypeError('Either existing_volume or filename must be provided!')
        if existing_volume and filename:
            raise TypeError('Cannot provide both existing_volume and filename!')
        from chimerax.clipper.symmetry import get_map_mgr
        map_mgr = get_map_mgr(m, create=True, auto_add_to_session=True)
        if existing_volume:
            nxmap = map_mgr.nxmapset.add_nxmap_handler_from_volume(existing_volume)
        else:
            nxmap = map_mgr.nxmapset.add_nxmap_handler_from_file(filename)
        self.change_selected_model(self.selected_model)


    @property
    def selected_model_has_maps(self):
        '''
        Returns True if the selected model has at least one map associated with
        it, otherwise False. Read only.
        '''
        m = self.selected_model
        if m is None:
            return False
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(m)
        if mgr is None:
            # No maps associated with this model.
            return False
        if not len(mgr.all_maps):
            return False
        return True

    def _initialize_maps(self, model):
        '''
        Set up all volumetric maps associated with a model, ready for simulation.
        '''
        from chimerax.map import Volume
        from chimerax.clipper.symmetry import get_map_mgr
        if model is None:
            return False
        mgr = get_map_mgr(model)
        if mgr is None:
            # No maps associated with this model.
            return False

        recalc_cb = self.iw._sim_basic_xtal_settings_live_recalc_checkbox

        from chimerax.clipper.maps import XmapHandler_Live
        from .session_extensions import get_mdff_mgr
        mdff_mgrs = []
        live_maps = False
        for xmapset in mgr.xmapsets:
            if xmapset.live_xmap_mgr is not None:
                # Need to make absolutely sure the MDFF potential map is created.
                # Handle this case specially.
                mdff_ps = [v for v in xmapset if isinstance(v, XmapHandler_Live) and 'MDFF potential' in v.name]
                if not len(mdff_ps):
                    mdff_ps = [self._add_mdff_potential_to_live_xmapset(xmapset)]
                for mdff_p in mdff_ps:
                    mdff_mgrs.append(get_mdff_mgr(model, mdff_p, create=True))

        for v in mgr.all_maps:
            # Only XmapHandler_Live objects specifically created to exclude the
            # free set are allowed to act as MDFF potentials
            if isinstance(v, XmapHandler_Live):
                live_maps = True
                continue
            # Difference maps are excluded from MDFF by default
            if v.is_difference_map:
                continue
            mdff_mgrs.append(get_mdff_mgr(model, v, create=True))
        recalc_cb.setVisible(live_maps)
        if len(mdff_mgrs):
            return True
        return False


    def _populate_available_volumes_combo_box(self, *_):
        '''
        Only true Volume instances (not subclasses) will be considered
        '''
        from chimerax.map import Volume
        shortlist = self.session.models.list(type = Volume)
        cb = self.iw._real_space_map_volume_combo_box
        cb.clear()
        for v in shortlist:
            if type(v) == Volume:
                label = '{}  {}'.format(v.id_string, v.name)
                cb.addItem(label, v)

    def _show_real_space_map_dialog(self, *_):
        self._populate_available_volumes_combo_box()
        self.iw._real_space_map_from_volume_frame.show()

    def _hide_real_space_map_dialog(self, *_):
        self._populate_available_volumes_combo_box()
        self.iw._real_space_map_from_volume_frame.hide()

    def _add_real_space_map_from_gui(self, *_):
        cb = self.iw._real_space_map_volume_combo_box
        i = cb.currentIndex()
        if i == -1:
            return
        v = cb.itemData(i)
        self.add_real_space_map(existing_volume=v, to_model=self.selected_model)

    def _toggle_xtal_map_dialog(self, *_):
        button = self.iw._sim_basic_xtal_map_settings_show_button
        frame = self.iw._sim_basic_xtal_map_settings_frame
        show_text = 'Show map settings dialogue'
        hide_text = 'Hide map settings dialogue'
        if button.text() == show_text:
            if (not hasattr(self, '_xtal_dialog_model_change_handler') or
                    self._xtal_dialog_model_change_handler is None):
                self._xtal_dialog_model_change_handler = self.triggers.add_handler(
                    'selected model changed',
                    self._xtal_dialog_model_changed_cb
                )
            self._populate_xtal_map_combo_box()
            frame.show()
            self._populate_xtal_map_params()
            button.setText(hide_text)
        else:
            frame.hide()
            button.setText(show_text)
            if self._xtal_dialog_model_change_handler is not None:
                self.triggers.remove_handler(self._xtal_dialog_model_change_handler)
                self._xtal_dialog_model_change_handler = None

    def _xtal_dialog_model_changed_cb(self, trigger_name, changes):
        self._populate_xtal_map_combo_box()
        #self._populate_xtal_map_params()

    def _set_live_xmap_recalc(self, state):
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(self.selected_model)
        if mgr is not None:
            for xs in mgr.xmapsets:
                xs.live_update = state

    def _populate_xtal_map_combo_box(self, *_):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        cb.clear()
        sm = self.selected_model
        if sm is None:
            return
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(sm)
        if mgr is None:
            return
        for v in mgr.all_maps:
            cb.addItem(v.name, v)

    def _displayed_decimal_places_and_step(self, number, sig_figs=3):
        from math import log, ceil, floor
        if number <= 0:
            return 2, 0.1
        places = max(ceil(-log(number, 10)+sig_figs-1), 0)
        step = 10**(floor(log(number, 10)-1))
        return places, step

    def _update_map_weight_box_settings(self, *_):
        sb = self.iw._sim_basic_xtal_map_weight_spin_box
        v = sb.value()
        dps, step = self._displayed_decimal_places_and_step(v)
        sb.setDecimals(dps)
        sb.setSingleStep(step)


    def _populate_xtal_map_params(self, *_):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        gb = self.iw._sim_basic_xtal_settings_set_button
        wf = self.iw._sim_basic_xtal_map_weight_frame
        eb = self.iw._sim_basic_xtal_settings_enable_mdff_checkbox
        this_map = cb.currentData()
        if cb.currentIndex() == -1 or this_map is None:
            gb.setEnabled(False)
            wf.setEnabled(False)
            eb.setEnabled(False)
            return
        gb.setEnabled(True)
        from .session_extensions import get_mdff_mgr
        mgr = get_mdff_mgr(self.selected_model, this_map)
        if mgr is not None:
            wf.setEnabled(True)
            eb.setEnabled(True)
            eb.setChecked(mgr.enabled)
            gk = mgr.global_k
            sb = self.iw._sim_basic_xtal_map_weight_spin_box
            sb.setValue(gk)
        else:
            wf.setEnabled(False)
            eb.setEnabled(False)
            eb.setChecked(False)

    def _enable_or_disable_mdff_potential(self, flag):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        this_map = cb.currentData()
        if this_map is None:
            return
        from .session_extensions import get_mdff_mgr
        mgr = get_mdff_mgr(self.selected_model, this_map)
        mgr.enabled = flag

    def _apply_xtal_map_params(self, *_):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        this_map = cb.currentData()
        if this_map is None:
            return
        from .session_extensions import get_mdff_mgr
        mgr = get_mdff_mgr(self.selected_model, this_map)
        mgr.global_k = self.iw._sim_basic_xtal_map_weight_spin_box.value()


    # Update button states after a simulation has finished
    def _update_menu_after_sim(self):
        self._update_sim_control_button_states()

    ####
    # Right mouse button modes
    ####

    def _set_right_mouse_tug_mode(self, mode):
        from .tugging import TugAtomsMode
        sm = self.sim_manager
        sc = sm.sim_construct
        t = TugAtomsMode(self.session, sm.tuggable_atoms_mgr, sc.mobile_atoms,
                spring_constant = self.sim_params.mouse_tug_spring_constant,
                mode=mode)
        mm = self.session.ui.mouse_modes
        mm.bind_mouse_mode('right', [], t)

    def _set_right_mouse_mode_tug_atom(self, *_):
        self._set_right_mouse_tug_mode('atom')

    def _set_right_mouse_mode_tug_residue(self, *_):
        self._set_right_mouse_tug_mode('residue')

    def _set_right_mouse_mode_tug_selection(self, *_):
        self._set_right_mouse_tug_mode('selection')



    ####
    # Rebuild tab
    ####
    def _enable_rebuild_residue_frame(self, res):
        from . import session_extensions
        pdmgr = session_extensions.get_proper_dihedral_mgr(self.session)
        # Peptide cis/trans flips
        self._rebuild_residue = res
        omega = self._rebuild_res_omega = pdmgr.get_dihedral(res, 'omega')
        if omega is None:
            # This is a terminal residue without an omega dihedral
            self.iw._rebuild_sel_res_pep_info.setText('N/A')
            self.iw._rebuild_sel_res_cis_trans_flip_button.setDisabled(True)
            self.iw._rebuild_sel_res_pep_flip_button.setDisabled(True)
            # return
        else:
            self.iw._rebuild_sel_res_cis_trans_flip_button.setEnabled(True)
            self.iw._rebuild_sel_res_pep_flip_button.setEnabled(True)

        # Rotamer manipulations
        rot_text = self.iw._rebuild_sel_res_rot_info
        rot_m = session_extensions.get_rotamer_mgr(self.session)
        rota = rot_m.get_rotamer(res)
        if rota != self._selected_rotamer:
            self._clear_rotamer()
            self._selected_rotamer = rota
        if rota is None:
            # This residue has no rotamers
            rot_text.setText('')
            self._set_rotamer_buttons_enabled(False)
        else:
            self._set_rotamer_buttons_enabled(True)

        self.iw._rebuild_sel_residue_frame.setEnabled(True)
        chain_id, resnum, resname = (res.chain_id, res.number, res.name)
        self.iw._rebuild_sel_residue_info.setText(str(chain_id) + ' ' + str(resnum) + ' ' + resname)

        self._steps_per_sel_res_update = 10
        self._sel_res_update_counter = 0
        self._update_selected_residue_info_live(None, None)
        if not hasattr(self, '_res_info_update_handler') or self._res_info_update_handler is None:
            self._res_info_update_handler = self.selected_model.triggers.add_handler(
                'changes', self._update_selected_residue_info_live
            )

    def _disable_rebuild_residue_frame(self):
        if hasattr(self, '_res_info_update_handler') and self._res_info_update_handler is not None:
            self.selected_model.triggers.remove_handler(self._res_info_update_handler)
            self._res_info_update_handler = None
        self.iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        self.iw._rebuild_sel_res_pep_info.setText('')
        self.iw._rebuild_sel_res_rot_info.setText('')
        self.iw._rebuild_sel_res_rot_target_button.setText('Set target')
        self.iw._rebuild_sel_residue_frame.setDisabled(True)

    def _enable_peptide_bond_manipulation_frame(self):
        self.iw._rebuild_sel_res_pep_frame.setEnabled(True)

    def _disable_peptide_bond_manipulation_frame(self):
        self.iw._rebuild_sel_res_pep_frame.setEnabled(False)

    def _enable_atom_position_restraints_frame(self):
        self.iw._rebuild_pin_atom_to_current_pos_button.setEnabled(True)
        self.iw._rebuild_pin_atom_to_pivot_button.setEnabled(True)

        #~ self.iw._rebuild_pin_atom_container.setEnabled(True)

    def _disable_atom_position_restraints_frame(self):
        self.iw._rebuild_pin_atom_to_current_pos_button.setEnabled(False)
        self.iw._rebuild_pin_atom_to_pivot_button.setEnabled(False)
        #self.iw._rebuild_pin_atom_container.setEnabled(False)

    def _enable_position_restraints_clear_button(self):
        self.iw._rebuild_pos_restraint_clear_button.setEnabled(True)

    def _disable_position_restraints_clear_button(self):
        self.iw._rebuild_pos_restraint_clear_button.setEnabled(False)

    def _enable_secondary_structure_restraints_frame(self, *_):
        self.iw._rebuild_2ry_struct_restr_container.setEnabled(True)

    def _disable_secondary_structure_restraints_frame(self, *_):
        self.iw._rebuild_2ry_struct_restr_container.setEnabled(False)

    def _enable_secondary_structure_restraints_clear_button(self, *_):
        self.iw._rebuild_2ry_struct_restr_clear_button.setEnabled(True)

    def _disable_secondary_structure_restraints_clear_button(self, *_):
        self.iw._rebuild_2ry_struct_restr_clear_button.setEnabled(False)

    def _enable_register_shift_frame(self, *_):
        self.iw._rebuild_register_shift_container.setEnabled(True)

    def _enable_selection_extend_frame(self, *_):
        self.iw._rebuild_grow_shrink_sel_frame.setEnabled(True)

    def _disable_register_shift_frame(self, *_):
        self.iw._rebuild_register_shift_container.setEnabled(False)

    def _disable_selection_extend_frame(self, *_):
        self.iw._rebuild_grow_shrink_sel_frame.setEnabled(False)

    def _extend_selection_by_one_res_N(self, *_):
        self._extend_selection_by_one_res(-1)

    def _shrink_selection_by_one_res_N(self, *_):
        self._shrink_selection_by_one_res(1)

    def _extend_selection_by_one_res_C(self, *_):
        self._extend_selection_by_one_res(1)

    def _shrink_selection_by_one_res_C(self, *_):
        self._shrink_selection_by_one_res(-1)

    def _enable_distance_restraints_frame(self, *_):
        self.iw._rebuild_dist_restraint_container.setEnabled(True)

    def _disable_distance_restraints_frame(self, *_):
        self.iw._rebuild_dist_restraint_container.setEnabled(False)

    def _enable_distance_restraint_apply_button(self, *_):
        self.iw._rebuild_dist_restraint_apply_button.setEnabled(True)

    def _disable_distance_restraint_apply_button(self, *_):
        self.iw._rebuild_dist_restraint_apply_button.setEnabled(False)

    def _enable_distance_restraint_set_distance_to_current_button(self, *_):
        self.iw._rebuild_dist_restraint_set_target_to_current_distance_button.setEnabled(True)

    def _disable_distance_restraint_set_distance_to_current_button(self, *_):
        self.iw._rebuild_dist_restraint_set_target_to_current_distance_button.setEnabled(False)

    def _enable_distance_restraint_remove_single_button(self, *_):
        self.iw._rebuild_remove_distance_restraint_button.setEnabled(True)

    def _disable_distance_restraint_remove_single_button(self, *_):
        self.iw._rebuild_remove_distance_restraint_button.setEnabled(False)


    def _extend_selection_by_one_res(self, direction):
        '''
        Extends a selection by one residue in the given direction, stopping
        when it hits a chain break or the end of a chain
        Args:
            direction:
                -1 or 1
        '''
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        residues = sel.unique_residues
        # expand the selection to cover all atoms in the existing residues
        residues.atoms.selected = True
        m = sel.unique_structures[0]
        polymers = m.polymers(
            missing_structure_treatment = m.PMS_NEVER_CONNECTS)
        for p in polymers:
            p = p[0]
            indices = p.indices(residues)
            if indices[0] != -1:
                break
        indices = numpy.sort(indices)
        residues = p[indices]
        first = residues[0]
        last = residues[-1]
        if direction == -1:
            i0 = indices[0]
            if i0 > 0:
                first = p[i0-1]
                first.atoms.selected = True
        else:
            iend = indices[-1]
            if iend < len(p) - 1:
                last = p[iend+1]
                last.atoms.selected = True

    def _shrink_selection_by_one_res(self, direction):
        '''
        Shrinks the current selection by one residue from one end. If
        direction == 1 the first residue will be removed from the selection,
        otherwise the last will be removed. The selection will never be
        shrunk to less than a single residue.
        '''
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        residues = sel.unique_residues
        if len(residues) > 1:
            if direction == 1:
                r = residues[0]
            elif direction == -1:
                r = residues[-1]
            else:
                raise TypeError('Direction must be either 1 or -1!')
            r.atoms.selected = False
            r.atoms.intra_bonds.selected = False

    def _restrain_selection_as_alpha_helix(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self.restrain_secondary_structure(sel, 'Helix')

    def _restrain_selection_as_antiparallel_beta(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self.restrain_secondary_structure(sel, 'Antiparallel Beta')

    def _restrain_selection_as_parallel_beta(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self.restrain_secondary_structure(sel, 'Parallel Beta')


    def restrain_secondary_structure(self, atoms, target):
        '''
        Restrain all amino acid residues in a selection to a target secondary
        structure. Secondary structure restraints consist of dihedral restraints
        on phi and psi, and distance restraints on (CA(n)-CA(n+2)) and
        (O(n)-O(n+4)).

        Args:
            * atoms:
                - a :class:`chimerax.Atoms` instance defining the selection. Any
                  residue with at least one atom in the selection will be
                  restrained. No distance restraints involving residues outside
                  the selection will be applied.
            * target:
                - a string selected from the following:
                    . 'Helix'
                    . 'Parallel Beta'
                    . 'Antiparallel Beta'

        Target distances/angles for each secondary structure definition are
        stored in `constants.py`.
        '''
        dihed_k = self.sim_params.phi_psi_spring_constant.value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
        dist_k = self.sim_params.distance_restraint_spring_constant.value_in_unit(OPENMM_SPRING_UNIT)
        residues = atoms.unique_residues
        m = self.selected_model
        from .restraints.constants import ss_restraints
        restraint_params = ss_restraints[target]
        from . import session_extensions as sx
        dr_m = sx.get_distance_restraint_mgr(m)
        o_to_n_plus_four, ca_to_ca_plus_two = dr_m.add_ss_restraints(residues)
        o_to_n_plus_four.targets = restraint_params.O_TO_N_PLUS_FOUR_DISTANCE
        ca_to_ca_plus_two.targets = restraint_params.CA_TO_CA_PLUS_TWO_DISTANCE
        o_to_n_plus_four.spring_constants = dist_k
        ca_to_ca_plus_two.spring_constants = dist_k
        o_to_n_plus_four.enableds = True
        ca_to_ca_plus_two.enableds = True

        pdr_m = sx.get_proper_dihedral_restraint_mgr(m)
        phi = pdr_m.add_restraints_by_residues_and_name(residues, 'phi')
        phi.targets = restraint_params.PHI_ANGLE
        phi.spring_constants = dihed_k
        phi.cutoffs = restraint_params.CUTOFF_ANGLE
        phi.enableds = True

        psi = pdr_m.add_restraints_by_residues_and_name(residues, 'psi')
        psi.targets = restraint_params.PSI_ANGLE
        psi.spring_constants = dihed_k
        psi.cutoffs = restraint_params.CUTOFF_ANGLE
        psi.enableds = True

    def _increment_register_shift(self, *_):
        self.iw._rebuild_register_shift_nres_spinbox.stepUp()

    def _decrement_register_shift(self, *_):
        self.iw._rebuild_register_shift_nres_spinbox.stepDown()

    def _apply_register_shift(self, *_):
        from chimerax.atomic import selected_atoms
        from .manipulations import Protein_Register_Shifter
        nshift = self.iw._rebuild_register_shift_nres_spinbox.value()
        sel = selected_atoms(self.session)
        rs = self._register_shifter = Protein_Register_Shifter(self.session, self, sel)
        rs.shift_register(nshift)
        self.iw._rebuild_register_shift_release_button.setEnabled(False)
        self.iw._rebuild_register_shift_go_button.setEnabled(False)
        self.add_checkpoint_block(rs, 'Register shift in progress')
        rs.triggers.add_handler('register shift finished', self._register_shift_finished)
        rs.triggers.add_handler('register shift released', self._register_shift_released_cb)

    def _register_shift_finished(self, *_):
        self.iw._rebuild_register_shift_release_button.setEnabled(True)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _release_register_shifter(self, *_):
        rs = self._register_shifter
        if rs is not None:
            rs.release_all()

    def _register_shift_released_cb(self, *_):
        rs = self._register_shifter
        if rs is not None:
            self.remove_checkpoint_block(rs)
        self._register_shifter = None
        self.iw._rebuild_register_shift_go_button.setEnabled(True)
        self.iw._rebuild_register_shift_release_button.setEnabled(False)

    def _restrain_selected_atom_to_current_xyz(self, *_):
        from chimerax.atomic import selected_atoms
        atom = selected_atoms(self.session)[0]
        # Displayed value in kJ mol^-1, A^2, OpenMM in kJ mol^-1 nm^2
        k = self.iw._rebuild_pos_restraint_spring_constant.value()*100
        self.restrain_atom_to_xyz(atom, atom.coord, k)

    def _restrain_selected_atom_to_pivot_xyz(self, *_):
        from chimerax.atomic import selected_atoms
        m = self.selected_model
        atom = selected_atoms(self.session)[0]
        # Displayed value in kJ mol^-1, A^2, OpenMM in kJ mol^-1 nm^2
        k = self.iw._rebuild_pos_restraint_spring_constant.value()*100
        cofr = self.session.view.center_of_rotation
        # Transform target from scene to model coordinates
        target = m.position.inverse(is_orthonormal=True)*cofr
        self.restrain_atom_to_xyz(atom, target, k)

    def restrain_atom_to_xyz(self, atom, target, spring_constant):
        '''
        Restrain the given atom to a (x,y,z) position.

        Args:
            * atom:
                - A :class:`chimerax.Atom` instance pointing to an atom in
                  the currently selected model. By default hydrogen atoms are
                  not restrainable.
            * target:
                - An iterable giving the target (x,y,z) position in Angstroms
            * spring_constant:
                - The desired spring constant in kJ mol-1 nm-2
        '''
        from . import session_extensions as sx
        pr_m = sx.get_position_restraint_mgr(self.selected_model)
        pr = pr_m.add_restraint(atom)
        pr.target = target
        pr.spring_constant = spring_constant
        pr.enabled = True

    def release_xyz_restraints_on_selected_atoms(self, *_, sel = None):
        '''
        Release current position restraints on a set of atoms.

        Args:
            * sel (default: None):
                - A :class:`chimerax.Atoms` instance giving the atoms to
                  release, or None. If None, all currently-selected atoms from
                  the current model will be released.
        '''
        if sel is None:
            from chimerax.atomic import selected_atoms
            sel = selected_atoms(self.session)
        from . import session_extensions as sx
        pr_m = sx.get_position_restraint_mgr(self.selected_model)
        prs = pr_m.get_restraints(sel)
        prs.enableds = False

    def _set_distance_restraint_target_to_current_distance(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        if len(sel) != 2:
            raise TypeError('A distance restraint must involve exactly two atoms!')
        coords = sel.coords
        d = numpy.linalg.norm(coords[1]-coords[0])
        self.iw._rebuild_dist_restraint_target_distance.setValue(d)

    def _add_distance_restraint_between_selected_atoms(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        if len(sel) != 2:
            raise TypeError('A distance restraint must involve exactly two atoms!')
        target = self.iw._rebuild_dist_restraint_target_distance.value()
        k = self.iw._rebuild_dist_restraint_spring_constant.value()
        self.add_distance_restraint(*sel, target, k * 100)

    def add_distance_restraint(self, atom1, atom2, target, spring_constant):
        '''
        Adds a distance restraint between two atoms with the given target
        distance and spring constant. If the restraint already exists, just
        changes the target/spring constant and enables it.

        Args:
            * atom1, atom2:
                - :class:`chimerax.Atom` instances pointing to the pair of
                  atoms to be restrained.
            * target:
                - The target distance in Angstroms
            * spring_constant:
                - The restraint spring constant in kJ mol-1 nm-2
        '''
        from . import session_extensions as sx
        dr_m = sx.get_distance_restraint_mgr(self.selected_model)
        dr = dr_m.add_restraint(atom1, atom2)
        dr.target = target
        dr.spring_constant = spring_constant
        dr.enabled = True

    def _release_distance_restraint_between_selected_atoms(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        from . import session_extensions as sx
        dr_m = sx.get_distance_restraint_mgr(self.selected_model)
        drs = dr_m.intra_restraints(sel)
        if len(drs):
            drs.enableds=False

    def _release_all_distance_restraints_for_selected_atoms(self, *_):
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(self.session)
        from . import session_extensions as sx
        dr_m = sx.get_distance_restraint_mgr(self.selected_model)
        drs = dr_m.atoms_restraints(sel)
        if len(drs):
            drs.enableds=False

    def _set_rotamer_buttons_enabled(self, switch):
        iw = self.iw
        iw._rebuild_cycle_rotamer_frame.setEnabled(switch)
        iw._rebuild_sel_res_last_rotamer_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_commit_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_target_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_discard_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_release_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_backrub_button.setEnabled(switch)

    def _update_rotamer_preview_text(self, name, freq):
        f = self.iw._rebuild_sel_res_rot_preview_info
        f.setText('<html><body><p>{}<br/>{:.3f}</p></body></html>'.format(
            name, freq
        ))

    def _clear_rotamer_preview_text(self):
        self._update_rotamer_preview_text('name', 0.0)

    def _next_rotamer(self, *_):
        rot = self._selected_rotamer
        from . import session_extensions
        rrm = session_extensions.get_rotamer_restraint_mgr(self.selected_model)
        target_def = rrm.next_preview(rot)
        self._update_rotamer_preview_text(target_def['Name'], target_def['Frequency'])

    def _prev_rotamer(self, *_):
        rot = self._selected_rotamer
        from . import session_extensions
        rrm = session_extensions.get_rotamer_restraint_mgr(self.selected_model)
        target_def = rrm.prev_preview(rot)
        self._update_rotamer_preview_text(target_def['Name'], target_def['Frequency'])

    def _commit_rotamer(self, *_):
        rot = self._selected_rotamer
        from . import session_extensions
        rrm = session_extensions.get_rotamer_restraint_mgr(self.selected_model)
        if not rrm.valid_preview(rot):
            self.session.logger.warning('No rotamer preview selected! Ignoring command.')
            return
        rrm.commit_preview(rot)
        if self.simulation_running:
            self.sim_handler.push_coords_to_sim()

    def _set_rotamer_target(self, *_):
        rot = self._selected_rotamer
        from . import session_extensions
        rrm = session_extensions.get_rotamer_restraint_mgr(self.selected_model)
        if not rrm.valid_preview(rot):
            self.session.logger.warning('No rotamer preview selected! Ignoring command.')
            return
        rrm.set_targets(rot)
        rr = rrm.get_restraint(rot)
        rr.set_spring_constant(self.sim_params.rotamer_spring_constant.value_in_unit(OPENMM_RADIAL_SPRING_UNIT))
        rr.enabled = True

    def _clear_rotamer(self, *_):
        from . import session_extensions
        rrm = session_extensions.get_rotamer_restraint_mgr(self.selected_model)
        rrm.remove_preview()
        self._target_rotamer = None
        self._clear_rotamer_preview_text()

    def _release_rotamer(self, *_):
        rot = self._selected_rotamer
        if rot is not None:
            self.release_rotamer(rot)

    def _backrub_rotamer(self, *_):
        self._clear_rotamer()
        self._release_rotamer()
        if self.simulation_running and not self.sim_manager.pause:
            # We need to pause the simulation and the backrub until after the
            # coordinates update
            self.sim_manager.pause=True
            self.sim_handler.triggers.add_handler('sim paused', self._apply_backrub)
        else:
            # Just do it immediately
            self._apply_backrub(None, None)

    def _apply_backrub(self, trigger_name, data):
        res = self._selected_rotamer.residue
        from.refine.backrub_rotamer import apply_backrub
        mdff_mgr = self.iw._rebuild_sel_res_rot_backrub_map_combo_box.currentData()
        if mdff_mgr is None:
            self._populate_rot_mdff_target_combo_box()
            mdff_mgr = self.iw._rebuild_sel_res_rot_backrub_map_combo_box.currentData()
        if mdff_mgr is not None:
            apply_backrub(self, mdff_mgr, res)
        if trigger_name == 'sim paused':
            # Function was called after automatically pausing the simulation.
            # Resume now
            self.session.triggers.add_handler('frame drawn', self._sim_delayed_resume_cb)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _populate_rot_mdff_target_combo_box(self, *_):
        from chimerax.map import Volume
        from .session_extensions import get_mdff_mgr
        m = self.selected_model
        cb = self.iw._rebuild_sel_res_rot_backrub_map_combo_box
        current_mgr = cb.currentData()
        cb.clear()
        if m is None:
            return
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(m)
        if mgr is None:
            return
        for v in mgr.all_maps:
            mdff_mgr = get_mdff_mgr(m, v)
            if mdff_mgr is not None:
                label = '{} {}'.format(v.id_string, v.name)
                cb.addItem(label, mdff_mgr)
        index = cb.findData(current_mgr)
        if index != -1:
            cb.setCurrentIndex(index)

    def release_rotamers(self, rotamers_or_residues):
        '''
        Release all chi dihedral restraints on a set of rotamers.

        Args:
            * rotamers_or_residues:
                - Either a :class:`Rotamers` or a :class:`chimerax.Residues`
                  instance
        '''
        from . import session_extensions as sx
        rrm = sx.get_rotamer_restraint_mgr(self.selected_model)
        rrs = rrm.get_restraints(rotamers_or_residues)
        rrs.enableds = False

    def release_rotamer(self, rotamer_or_residue):
        '''
        Release all chi dihedral restraints on a single rotamer.

        Args:
            * rotamer_or_residue:
                - A :class:`Rotamer` or a :class:`chimerax.Residue` instance
        '''
        from . import session_extensions as sx
        rrm = sx.get_rotamer_restraint_mgr(self.selected_model)
        rr = rrm.get_restraint(rotamer_or_residue)
        if rr is not None:
            rr.enabled = False

    def _update_selected_residue_info_live(self, trigger_name, changes):
        if changes is not None:
            changes = changes[1]
            if not 'coord changed' in changes.atom_reasons():
                return
        from math import degrees
        res = self._rebuild_residue
        if True:
            # Get the residue's omega value
            omega = self._rebuild_res_omega
            if omega is not None:
                oval = degrees(omega.angle)
                if oval <= -150 or oval >= 150:
                    pep_type = 'trans'
                elif oval >= -30 and oval <= 30:
                    pep_type = 'cis'
                else:
                    pep_type = 'twisted'
                self.iw._rebuild_sel_res_pep_info.setText(pep_type)
            rot = self._selected_rotamer
            rot_text = self.iw._rebuild_sel_res_rot_info
            from .session_extensions import get_rotamer_mgr
            rot_m = get_rotamer_mgr(self.session)

            if rot is not None:
                t_info = rot.nearest_target
                r_desc = "{}, f={}, Z=({})".format(
                    t_info['Name'],
                    t_info['Frequency'],
                    ','.join(['{:0.2f}'.format(z) for z in t_info['Z scores']])
                    )
                score = rot.score
                from . import session_extensions as sx
                rot_m = sx.get_rotamer_mgr(self.session)
                cutoffs = rot_m.cutoffs
                if score > cutoffs[0]:
                    desc_color = 'green'
                elif score > cutoffs[1]:
                    desc_color = 'orange'
                else:
                    desc_color = 'red'
                rot_text.setText("<font color='{}'>{}</font>".format(desc_color, r_desc))


    def _flip_peptide_bond(self, *_):
        res = self._rebuild_residue
        self.flip_peptide_bond(res)

    def _flip_cis_trans(self, *_):
        res = self._rebuild_residue
        self.flip_peptide_omega(res)

    ####
    # Validation tab
    ####
    def _prepare_ramachandran_plot(self):
        '''
        Prepare an empty MatPlotLib figure to put the Ramachandran plots in.
        '''
        iw = self.iw
        iw._validate_rama_go_button.setEnabled(False)
        container = self._rama_plot_window = iw._validate_rama_plot_layout
        from .validation.ramaplot import RamaPlot
        mode_menu = iw._validate_rama_sel_combo_box
        case_menu = iw._validate_rama_case_combo_box
        restrict_button = iw._validate_rama_go_button
        self._rama_plot = RamaPlot(self.session, self, container, mode_menu,
            case_menu, restrict_button)

    def _show_rama_plot(self, *_):
        self.iw._validate_rama_stub_frame.hide()
        self.iw._validate_rama_main_frame.show()
        if self._rama_plot is None:
            # Create the basic MatPlotLib canvas for the Ramachandran plot
            self._prepare_ramachandran_plot()

    def _hide_rama_plot(self, *_):
        self.iw._validate_rama_main_frame.hide()
        self.iw._validate_rama_stub_frame.show()

    def _show_peptide_validation_frame(self, *_):
        self.iw._validate_pep_stub_frame.hide()
        self.iw._validate_pep_main_frame.show()
        self._update_iffy_peptide_lists()

    def _hide_peptide_validation_frame(self, *_):
        self.iw._validate_pep_main_frame.hide()
        self.iw._validate_pep_stub_frame.show()

    def _update_iffy_peptide_lists(self, *_):
        from .session_extensions import get_proper_dihedral_mgr
        pd_mgr = get_proper_dihedral_mgr(self.session)
        model = self.selected_model
        table = self.iw._validate_pep_iffy_table
        table.setRowCount(0)
        if model is None:
            return
        if self.simulation_running:
            residues = self.sim_manager.sim_construct.mobile_residues
        else:
            residues = model.residues
        omegas = pd_mgr.get_dihedrals(residues, 'omega')
        abs_angles = numpy.abs(omegas.angles)
        from math import pi
        from .constants import defaults
        cc = defaults.CIS_PEPTIDE_BOND_CUTOFF
        tc = defaults.TWISTED_PEPTIDE_BOND_DELTA
        cis_mask = abs_angles < cc
        twisted_mask = numpy.logical_and(abs_angles >= cc, abs_angles < pi-tc)

        iffy_mask = numpy.logical_or(cis_mask, twisted_mask)
        iffy = omegas[iffy_mask]
        angles = numpy.degrees(iffy.angles)
        cis_mask = cis_mask[iffy_mask]

        table.setRowCount(len(iffy))
        from PyQt5.QtWidgets import QListWidgetItem
        from PyQt5.Qt import QColor, QBrush
        from PyQt5.QtCore import Qt
        cis_nonpro_color = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        cis_pro_color = QBrush(QColor(100,255,100), Qt.SolidPattern)
        twisted_color = QBrush(QColor(240, 200, 160), Qt.SolidPattern)
        from PyQt5.QtWidgets import QTableWidgetItem
        for i, (omega, angle, cis) in enumerate(zip(iffy, angles, cis_mask)):
            res1, res2 = omega.atoms.unique_residues
            if cis:
                conf_text = 'cis'
            else:
                conf_text = 'twisted'
            data = (
                res1.chain_id,
                '{}-{}'.format(res1.number, res2.number),
                '{}-{}'.format(res1.name, res2.name),
                '{:.0f}° ({})'.format(angle, conf_text)
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.data = res2
                if cis:
                    if res2.name == 'PRO':
                        color = cis_pro_color
                    else:
                        color = cis_nonpro_color
                else:
                    color = twisted_color
                item.setBackground(color)
                table.setItem(i, j, item)

    def _show_selected_iffy_peptide(self, item):
        res = item.data
        from . import view
        view.focus_on_selection(self.session, res.atoms)
        self.session.selection.clear()
        res.atoms.selected = True

    def _show_rota_validation_frame(self, *_):
        self.iw._validate_rota_stub_frame.hide()
        self.iw._validate_rota_main_frame.show()
        self._update_iffy_rota_list()

    def _hide_rota_validation_frame(self, *_):
        self.iw._validate_rota_stub_frame.show()
        self.iw._validate_rota_main_frame.hide()

    def _update_iffy_rota_list(self, *_):
        from .session_extensions import get_rotamer_mgr
        rota_m = get_rotamer_mgr(self.session)
        table = self.iw._validate_rota_table
        table.setRowCount(0)
        if self.selected_model is None:
            return
        if self.simulation_running:
            residues = self.sim_manager.sim_construct.mobile_residues
        else:
            residues = self.selected_model.residues
        rotas = rota_m.get_rotamers(residues)
        iffy, scores = rota_m.non_favored_rotamers(rotas)
        order = numpy.argsort(scores)
        outlier_cutoff = rota_m.cutoffs[1]
        from PyQt5.Qt import QColor, QBrush
        from PyQt5.QtCore import Qt
        badColor = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        table.setRowCount(len(iffy))
        from PyQt5.QtWidgets import QTableWidgetItem
        for i, index in enumerate(order):
            r = iffy[index]
            score = scores[index]
            res = r.residue
            data = (
                res.chain_id,
                str(res.number),
                res.name,
                '{:.4f}'.format(score*100)
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.data = res
                if score < outlier_cutoff:
                    item.setBackground(badColor)
                table.setItem(i, j, item)
        table.resizeColumnsToContents()

    def _show_selected_iffy_rota(self, item):
        res = item.data
        from . import view
        view.focus_on_selection(self.session, res.atoms)
        self.session.selection.clear()
        res.atoms.selected = True

    ##############################################################
    # Simulation global settings functions
    ##############################################################

    def _change_force_field(self):
        ff_key = self.iw._sim_force_field_combo_box.currentText()
        self.sim_params.forcefield = ff_key

    def change_selected_model(self, model):
        '''
        Change the model upon which ISOLDE is focused (that is, upon which
        simulations will be run). This performs some basic preparations,
        including preparing any maps associated with the model for simulation.
        While ISOLDE is running, only atoms from the selected model will be
        selectable using the mouse.

        Args:
            * model:
                - a :class:`chimerax.AtomicStructure` instance.
        '''
        if self.gui_mode:
            self._change_selected_model(self, model=model, force=True)
        else:
            self._selected_model = m
            self.session.selection.clear()
            self._selected_model.selected = True
            self._initialize_maps(m)

    def _change_selected_model(self, *_, model = None, force = False):
        if not hasattr(self, '_model_changes_handler'):
            self._model_changes_handler = None
        sm = self._selected_model
        if sm is not None:
            if self._model_changes_handler is not None:
                sm.triggers.remove_handler(self._model_changes_handler)
                self._model_changes_handler = None

        if len(self._available_models) == 0:
            self._selected_model = None
            self._update_iffy_rota_list()
            self._update_iffy_peptide_lists()
            return
        if self.simulation_running:
            return
        iw = self.iw
        mmcb = iw._master_model_combo_box
        if not mmcb.count():
            self._selected_model = None
            self._update_iffy_rota_list()
            self._update_iffy_peptide_lists()
            return
        old_index = mmcb.currentIndex()
        if model is not None:
            # Find and select the model in the master combo box, which
            # will automatically call this function again with model = None
            index = iw._master_model_combo_box.findData(model)
            if index == -1 and mmcb.count():
                index = 0
            mmcb.setCurrentIndex(index)
            # If the new index happens to be the same as the old one, then the
            # currentIndexChanged trigger won't fire.
            if index != old_index:
                return
        if mmcb.currentIndex() == -1:
            mmcb.setCurrentIndex(0)
        m = mmcb.currentData()
        if force or (self._selected_model != m and m is not None):
            from chimerax.clipper.symmetry import get_symmetry_handler
            sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
            self._selected_model = m
            self._model_changes_handler = m.triggers.add_handler('changes',
                self._model_changes_cb)
            self.session.selection.clear()
            m.selected = True
            has_maps = self._initialize_maps(m)
            if has_maps:
                iw._map_masking_frame.setEnabled(True)
            else:
                iw._map_masking_frame.setEnabled(False)

            # Load/create validation managers
            from . import session_extensions as sx
            sx.get_rota_annotator(m)
            sx.get_rama_annotator(m)
            self._populate_rot_mdff_target_combo_box()
            self.triggers.activate_trigger('selected model changed', data=m)
            self._update_iffy_rota_list()
            self._update_iffy_peptide_lists()
        self._status('')

    def _change_b_and_a_padding(self, *_):
        self.params.num_selection_padding_residues = self.iw._sim_basic_mobile_b_and_a_spinbox.value()

    def _change_soft_shell_cutoff_from_sel_menu(self, *_):
        iw = self.iw
        val = iw._sim_basic_mobile_sel_within_spinbox.value()
        self.params.soft_shell_cutoff_distance = val

    def _change_sim_platform(self, *_):
        self.sim_platform = self.iw._sim_platform_combo_box.currentText()

    def _change_mask_radius(self, *_):
        self.params.standard_map_mask_cutoff = self.iw._sim_basic_xtal_settings_mask_radius_spinbox.value()

    def _change_spotlight_radius(self, *_):
        from chimerax.core.commands import run
        radius = self.iw._sim_basic_xtal_settings_spotlight_radius_spinbox.value()
        run(self.session, 'clipper spotlight radius {:.2f}'.format(radius))


    ##############################################################
    # Visualisation functions
    ##############################################################

    def _xtal_step_forward(self, *_):
        m = self.selected_model
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        focus = self.iw._vis_focus_on_sel_checkbox.checkState()
        m.atoms.selected = False
        sel = sh.stepper.step_forward()
        self._xtal_mask_to_atoms(sel, focus=False)
        if focus:
            from .view import focus_on_selection
            focus_on_selection(self.session, sel[0].residue.atoms)
        sel.selected = True


    def _xtal_step_backward(self, *_):
        m = self.selected_model
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        focus = self.iw._vis_focus_on_sel_checkbox.checkState()
        m.atoms.selected = False
        sel = sh.stepper.step_backward()
        self._xtal_mask_to_atoms(sel, focus=False)
        if focus:
            from .view import focus_on_selection
            focus_on_selection(self.session, sel[-1].residue.atoms)
        sel.selected = True

    def _xtal_mask_to_selection(self, *_):
        atoms = self.selected_model.atoms
        sel = atoms[atoms.selecteds]
        focus = self.iw._vis_focus_on_sel_checkbox.checkState()
        self._xtal_mask_to_atoms(sel, focus)

    def _xtal_mask_to_atoms(self, atoms, focus):
        m = self.selected_model
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        cutoff = self.params.standard_map_mask_cutoff
        context = self.params.soft_shell_cutoff_distance
        sh.isolate_and_cover_selection(
            atoms, 0, context, cutoff, focus=focus)

    def _xtal_enable_live_scrolling(self, *_):
        m = self.selected_model
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        sh.spotlight_mode = True

    def _set_map_to_solid_surface(self, *_):
        v = self.iw._sim_basic_xtal_settings_map_combo_box.currentData()
        if v is None:
            return
        from .visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [v], **map_style_settings[map_styles.solid_opaque])

    def _set_map_to_transparent_surface(self, *_):
        v = self.iw._sim_basic_xtal_settings_map_combo_box.currentData()
        if v is None:
            return
        from .visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        if v.is_difference_map:
            style = map_styles.solid_t40
        else:
            style = map_styles.solid_t20
        volumecommand.volume(self.session, [v], **map_style_settings[style])

    def _set_map_to_mesh(self, *_):
        v = self.iw._sim_basic_xtal_settings_map_combo_box.currentData()
        if v is None:
            return
        from .visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [v], **map_style_settings[map_styles.mesh_triangle])

    def _choose_map_color(self, *_):
        v = self.iw._sim_basic_xtal_settings_map_combo_box.currentData()
        from PyQt5.QtWidgets import QColorDialog
        cd = QColorDialog(self.iw._map_settings_color_button)
        sa = cd.ShowAlphaChannel
        colors = []
        if v.is_difference_map:
            colors.append(cd.getColor(title="Colour for negative contour", options=sa))
            colors.append(cd.getColor(title="Colour for positive contour", options=sa))
        else:
            colors.append(cd.getColor(options=sa))

        import numpy
        carg = [
            numpy.array([c.red(), c.green(), c.blue(), c.alpha()], dtype=numpy.double)/255
                for c in colors
        ]
        v.set_parameters(surface_colors=carg)

    def _change_smoothing_state_from_gui(self, *_):
        flag = self.iw._trajectory_smooth_button.isChecked()
        self.sim_params.trajectory_smoothing = flag
        if self.simulation_running:
            self.sim_handler.smoothing = flag

    def _change_smoothing_amount_from_gui(self, *_):
        sval = self.iw._smoothing_amount_dial.value()
        alpha = 10**-(sval/100)
        mina, maxa = defaults.SMOOTHING_ALPHA_MIN, defaults.SMOOTHING_ALPHA_MAX
        if alpha < mina:
            alpha = mina
        elif alpha > maxa:
            alpha = maxa
        self.sim_params.smoothing_alpha = alpha
        if self.simulation_running:
            self.sim_handler.smoothing_alpha = alpha


    ##############################################################
    # Simulation prep
    ##############################################################

    def reset_sim_params(self):
        '''
        Reset all the simulation parameters back to the defaults found in
        `constants.py`.
        '''
        from .openmm.sim_param_mgr import SimParams
        self.sim_params = SimParams()

    def _start_sim_or_toggle_pause(self, *_):
        if not self.simulation_running:
            self.start_sim()
        else:
            self.pause_sim_toggle()

    def start_sim(self):
        '''
        Start an interactive simulation based around the currently-selected
        atoms.
        '''
        if self.simulation_running:
            raise TypeError('Simulation already running!')
        self.sim_params.platform = self.iw._sim_platform_combo_box.currentText()
        from .openmm.openmm_interface import Sim_Manager
        main_sel = self._last_main_sel = self._get_main_sim_selection()
        try:
            sm = self._sim_manager = Sim_Manager(self, self.selected_model, main_sel,
                self.params, self.sim_params, excluded_residues = self.ignored_residues)
        except ValueError as e:
            err_text = str(e)
            if err_text.startswith("No template found") or \
               err_text.startswith("User-supplied template"):
                # These errors are already handled
                return
            raise
        except:
            raise
        sm.start_sim()
        self._sim_start_cb()

    def _start_sim_haptics(self, *_):
        self._event_handler.add_event_handler('sim haptic update',
                                              'new frame',
                                              self._update_haptics)

    def _stop_sim_haptics(self, *_):
        self._event_handler.remove_event_handler('sim haptic update')

    #TODO: Remove all haptics code from isolde.py
    # def _update_haptics(self, *_):
    #     hh = self.session.HapticHandler
    #     si = self._sim_interface
    #     from . import picking
    #     for d in self._haptic_devices:
    #         i = d.index
    #         pointer_pos = hh.getPosition(i, scene_coords = True)
    #         if self._haptic_highlight_nearest_atom[i] and not d.tugging:
    #             cur_a = self._haptic_current_nearest_atom[i]
    #             new_a = self._haptic_current_nearest_atom[i] = picking.pick_closest_to_point(
    #                     self.session, pointer_pos, self.sim_manager.sim_construct.mobile_atoms, 3.0,
    #                     displayed_only = True, hydrogens = self.tug_hydrogens)
    #             if new_a != cur_a:
    #                 if cur_a is not None and not cur_a.deleted:
    #                     # set the previously picked atom back to standard
    #                     cur_a.draw_mode = self._haptic_current_nearest_atom_draw_mode[i]
    #                     cur_a.color = self._haptic_current_nearest_atom_color[i]
    #                 if new_a is not None:
    #                     self._haptic_current_nearest_atom_draw_mode[i] = new_a.draw_mode
    #                     new_a.draw_mode = 1
    #                     self._haptic_current_nearest_atom_color[i] = new_a.color
    #                     new_a.color = [0, 255, 255, 255] # bright cyan
    #
    #
    #         b = hh.getButtonStates(i)
    #         # Button 0 controls tugging
    #         if b[0]:
    #             if not d.tugging:
    #                 a = self._haptic_tug_atom[i] = cur_a
    #                 if a is not None:
    #                     d.start_tugging(a)
    #                     target = hh.getPosition(i, scene_coords = True)
    #                     si.tug_atom_to(a, target)
    #
    #             else:
    #                 d.update_target()
    #                 a = self._haptic_tug_atom[i]
    #                 si.tug_atom_to(a, hh.getPosition(i, scene_coords = True))
    #         else:
    #             if d.tugging:
    #                 d.stop_tugging()
    #                 a = self._haptic_tug_atom[i]
    #                 si.release_tugged_atom(a)
    #                 self._haptic_tug_atom[i] = None

    def _handle_bad_template(self, residue):
        '''
        Called if OpenMM encounters a residue it doesn't recognise while
        attempting to prepare a simulation.
        '''
        self.selected_model.atoms.selected = False
        residue.atoms.selected = True
        from .view import focus_on_selection
        focus_on_selection(self.session, residue.atoms)
        from .dialog import failed_template_warning
        choice = failed_template_warning(residue)
        if choice == 'addh':
            print('Adding hydrogens')
            from chimerax.atomic import AtomicStructures
            from chimerax.atomic.addh import cmd
            cmd.cmd_addh(self.session, AtomicStructures([self.selected_model]), hbond=True)
            self._sim_end_cb()
            self.selected_model.atoms.selected = False
            self._last_main_sel.selected = True
            self.start_sim()
        elif choice == 'exclude':
            print('Excluding residue')
            from chimerax.atomic import Residues
            self.ignored_residues = self.ignored_residues.merge(Residues([residue]))
            self._sim_end_cb()
            self.selected_model.atoms.selected = False
            self._last_main_sel.selected=True
            self.start_sim()
        else:
            print('Doing nothing')
            self._sim_end_cb()

    def _sim_start_cb(self, *_):
        '''
        Register all event handlers etc. that have to be running during the
        simulation.
        '''
        self._update_sim_control_button_states()
        self._set_right_mouse_mode_tug_atom()
        self.triggers.activate_trigger('simulation started', None)
        sh = self.sim_handler
        sh.triggers.add_handler('sim terminated', self._sim_end_cb)
        sh.triggers.add_handler('sim paused', self._update_sim_control_button_states)
        sh.triggers.add_handler('sim resumed', self._update_sim_control_button_states)

    def _sim_end_cb(self, *_):
        self._update_menu_after_sim()
        for d in self._haptic_devices:
            d.cleanup()
        from chimerax.mouse_modes import TranslateMouseMode
        self.session.ui.mouse_modes.bind_mouse_mode('right', [], TranslateMouseMode(self.session))
        self.triggers.activate_trigger('simulation terminated', None)

    def _get_main_sim_selection(self):
        '''
        Get the primary selection that will be the focus of the simulation. In
        most cases, this involves taking the existing selected atoms, expanding
        the selection to cover whole residues, and further expanding the
        selection by the desired padding up and down the chain (stopping at
        chain breaks and ends).
        '''
        sm = self.selected_model
        all_atoms = sm.atoms
        pad = self.params.num_selection_padding_residues
        selatoms = sm.atoms[sm.atoms.selected]
        from chimerax.atomic import selected_atoms
        selected_atoms(self.session).selected=False
        selatoms.selected = True
        if not len(selatoms):
            raise TypeError('You must select at least one atom from the current '
                'working model prior to starting a simulation!')
        all_res = sm.residues
        sel_res = selatoms.unique_residues
        sel_res_indices = all_res.indices(sel_res)
        from chimerax.atomic import Structure

        allfrags = sm.polymers(missing_structure_treatment=Structure.PMS_NEVER_CONNECTS)

        sel_frags = []
        sel_frag_res_indices = []
        for frag in allfrags:
            frag = frag[0]
            if not frag.atoms.num_selected:
                continue
            sel_frags.append(frag)
            sel_frag_res_indices.append(all_res.indices(frag))

        for frag, frag_indices in zip(sel_frags, sel_frag_res_indices):
            frag_nres = len(frag_indices)
            sel_mask = numpy.isin(frag_indices, sel_res_indices, assume_unique=True)
            sel_pol_indices = numpy.where(sel_mask)[0]
            for i in sel_pol_indices:
                lb = i-pad
                ub = i+pad+1
                if lb<0:
                    lb = 0
                if ub > frag_nres:
                    ub = frag_nres
                sel_mask[lb:ub] = True
            frag[sel_mask].atoms.selected = True
        return selected_atoms(self.session)

    ##############################################################
    # Simulation on-the-fly control functions
    ##############################################################

    def _sim_delayed_resume_cb(self, *_):
        if self.simulation_running and self.sim_paused:
            self.sim_manager.pause=False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def pause_sim_toggle(self):
        '''
        Pause/resume the simulation
        '''
        if self.simulation_running:
            self.sim_manager.toggle_pause()

    @property
    def sim_paused(self):
        '''
        Returns True if the simulation is paused or no simulation is running,
        otherwise False.
        '''
        if self.simulation_running:
            return self.sim_manager.pause
        return True

    def _sim_pause_cb(self, *_):
        self._status('Simulation paused')
        go_button = self.iw._sim_go_button
        go_button.setChecked(False)
        go_button.setToolTip('Resume')
        self.triggers.activate_trigger('simulation paused', None)

    def _sim_resume_cb(self, *_):
        self._status('Simulation running')
        go_button = self.iw._sim_go_button
        go_button.setChecked(True)
        go_button.setToolTip('Pause')
        self.triggers.activate_trigger('simulation resumed', None)

    def _stop_sim_and_revert_to_checkpoint(self, *_):
        self.discard_sim(revert_to='checkpoint')

    def _discard_sim(self, *_):
        self.discard_sim(revert_to='start')

    def discard_sim(self, revert_to='checkpoint', warn=True):
        '''
        Stop the simulation and revert to either the starting state or
        the last saved checkpoint.

        Args:
            * revert_to (default: 'checkpoint'):
                - Either 'checkpoint' or 'start'
            * warn (default: True):
                - if warn is `True` and checkpoint=='start', a warning dialog
                  will be raised allowing the user to back out before discarding
                  the simulation coordinates.
        '''
        if not revert_to in ('checkpoint', 'start'):
            raise TypeError('Unrecognised option! Argument should be '\
                +'either "checkpoint" or "start".')
        if warn and revert_to == 'start':
            msg = 'All changes since you started this simulation will be '\
                +'lost! Are you sure you want to continue?'
            if self.gui_mode:
                from .dialog import choice_warning
                ok = choice_warning(msg)
            else:
                result = input(msg + ' (y or n)')
                while result.lower() not in ('y', 'n'):
                    result = input ('Please enter either "y" or "n"')
                ok = result.lower()=='y'
            if not ok:
                return
        self._release_register_shifter()
        self.sim_manager.stop_sim(revert=revert_to)

    def commit_sim(self):
        '''
        Stop the simulation and keep the current coordinates.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._release_register_shifter()
        self.sim_manager.stop_sim()

    def minimize(self):
        '''
        Switch the simulation to energy minimisation mode
        '''
        if self.simulation_running:
            self.sim_handler.minimize = True
        self.simulation_type = 'min'
        self._update_sim_control_button_states()

    def equilibrate(self):
        '''
        Switch the imulation to equilibration mode
        '''
        if self.simulation_running:
            self.sim_handler.minimize = False
        self.simulation_type = 'equil'
        self._update_sim_control_button_states()

    def tug_atom_to(self, atom, target, spring_constant = None):
        '''
        Tug an atom towards the given target coordinates in the context
        of a running simulation. If no spring constant is given, a
        default value will be used. NOTE: the tugging effect will
        remain in place until updated with a new call to tug_atom_to()
        or removed with stop_tugging(atom).

        Args:
            * atom:
                - The atom to be tugged. Must be a heavy atom that is
                  mobile in the current simulation
            * target:
                - An iterable giving the (x,y,z) target in Angstroms
            * spring constant (default None):
                - An optional spring constant (in kJ mol-1 A-2). If none is
                  provided, the default in `sim_params.mouse_tug_spring_constant`
                  will be used.
        '''
        if not self.simulation_running:
            return
        from . import session_extensions as sx
        tugm = sx.get_tuggable_atoms_mgr(self.selected_model)
        t_atom = tugm.get_tuggable(atom)
        t_atom.target = target
        if spring_constant is None:
            spring_constant = self.sim_params.mouse_tug_spring_constant
        t_atom.spring_constant = spring_constant
        t_atom.enabled = True

    def stop_tugging(self, atom_or_atoms):
        '''
        Release one or more atoms from the tugging force.

        Args:
            * atom_or_atoms:
                - A :class:`chimerax.Atom` or :class:`Atoms` instance
        '''
        if not self.simulation_running:
            return
        from chimerax.atomic import Atom, Atoms
        if isinstance(atom_or_atoms, Atom):
            atoms = Atoms([atom])
        elif isinstance(atom_or_atoms, Atoms):
            atoms = atom_or_atoms
        else:
            raise TypeError('Not an Atom or Atoms instance!')
        from . import session_extensions as sx
        tugm = sx.get_tuggable_atoms_mgr(self.selected_model)
        t_atoms = tugm.get_tuggables(atoms)
        t_atoms.enableds = False

    def add_checkpoint_block(self, obj, reason):
        '''
        Some processes are incompatible with checkpointing. We need to
        know when these are running and disable checkpointing until
        they're done.

        Args:
            * obj:
                - The Python object blocking checkpointing
            * reason:
                - A text string.
        '''
        self.checkpoint_disabled_reasons[obj] = reason
        self.iw._sim_save_checkpoint_button.setEnabled(False)
        self.iw._sim_revert_to_checkpoint_button.setEnabled(False)

    def remove_checkpoint_block(self, obj):
        '''
        Release a block on checkpointing. When all blocks are removed,
        checkpointing will be re-enabled.

        Args:
            * obj:
                - the Python object no longer blocking checkpointing
        '''
        r = self.checkpoint_disabled_reasons
        r.pop(obj)
        if len(r) ==0:
            self.iw._sim_save_checkpoint_button.setEnabled(True)
            self.iw._sim_revert_to_checkpoint_button.setEnabled(True)

    @property
    def can_checkpoint(self):
        '''Is checkpoint save/revert currently allowed?'''
        return self.simulation_running and not len(self.checkpoint_disabled_reasons)

    def checkpoint(self, *_):
        '''
        Save the current state of the simulation (coordinates and all
        restraints).
        '''
        if self.can_checkpoint:
            self.sim_manager.checkpoint()
        else:
            err_str = 'Checkpointing is currently disabled by the '\
                +'following scripts and will be re-enabled when they '\
                +'terminate: \n{}'.format(
                    '\n'.join([r for r in self.checkpoint_disabled_reasons.values()]))
            self.session.logger.warning(err_str)

    def revert_to_checkpoint(self, *_):
        '''
        Return the simulation to the last saved checkpoint. Note: a checkpoint
        is automatically saved at the beginning of each simulation.
        '''
        if self.can_checkpoint:
            self.sim_manager.revert_to_checkpoint()
        else:
            err_str = 'Checkpointing is currently disabled by the '\
                +'following scripts and will be re-enabled when they '\
                +'terminate: \n{}'.format(
                    '\n'.join([r for r in self.checkpoint_disabled_reasons.values()]))
            self._log.warning(err_str)


    @property
    def smoothing(self):
        '''
        Turn smoothing of the simulation trajectory visualisation on or off.
        '''
        return self.sim_params.trajectory_smoothing

    @smoothing.setter
    def smoothing(self, flag):
        if self.gui_mode:
            # then just let the GUI controls handle it
            self.iw._smoothing_button.setChecked(flag)
        self.sim_params.trajectory_smoothing = flag
        if self.simulation_running:
            self.sim_handler.smoothing = flag

    @property
    def smoothing_alpha(self):
        '''
        Get/set the value of the constant defining how heavily smoothed the
        trajectory visualisation is. Value must be greater than zero and less
        than or equal to one. A value of 1 indicates no smoothing, while a
        value of 0.01 indicates that a given frame contributes only 1% to the
        displayed moving average.
        '''
        return self.sim_params.smoothing_alpha

    @smoothing_alpha.setter
    def smoothing_alpha(self, alpha):
        if self.gui_mode:
            # then just let the GUI controls handle it
            slider = self.iw._smoothing_amount_dial
            from math import log
            la = log(alpha, 10)
            slider.setValue(int(-100*la))
        else:
            if alpha < defaults.SMOOTHING_ALPHA_MIN:
                alpha = defaults.SMOOTHING_ALPHA_MIN
            elif alpha > defaults.SMOOTHING_ALPHA_MAX:
                alpha = defaults.SMOOTHING_ALPHA_MAX
            self.sim_params.smoothing_alpha = alpha
            if self.simulation_running:
                self.sim_handler.smoothing_alpha = alpha


    ####
    # Restraint controls
    ####

    def flip_peptide_bond(self, res):
        '''
        Flip the peptide bond N-terminal to the given residue by the addition
        of temporary phi/psi restraints. Requires a simulation to be running.

        Args:
            * res:
                - A ChimeraX :class:`Residue` instance
        '''
        from .manipulations.peptide_flip import Peptide_Bond_Flipper
        pf = Peptide_Bond_Flipper(self, res)

    def clear_secondary_structure_restraints_for_selection(self, *_, atoms = None, residues = None):
        '''
        Clear all secondary structure restraints in the selection. If
        no atoms or residues are provided, restraints will be cleared
        for any atoms selected in the main window.

        Args:
            * atoms (default: None):
                - A :class:`chimerax.Atoms` instance or None
            * residues (default: None):
                - A :class:`chimerax.Residues` instance or None

        `atoms` and `residues` arguments should not both be given in the same
        call. If neither is given, then all currently-selected atoms in the
        model will be released.
        '''
        from chimerax.atomic import selected_atoms
        sel = None
        if atoms is not None:
            if residues is not None:
                raise TypeError('Cannot provide both atoms and residues!')
            sel = atoms
        elif residues is None:
            sel = selected_atoms(self.session)
        if sel is not None:
            residues = sel.unique_residues

        atoms = residues.atoms
        from . import session_extensions as sx
        m = self.selected_model
        dr_m = sx.get_distance_restraint_mgr(m)
        all_ss_restraints = dr_m.get_ss_restraints(m.residues)
        on_restraints = dr_m.atoms_restraints(atoms[numpy.in1d(atoms.names, ('O', 'N'))])
        on_restraints[all_ss_restraints[0].indices(on_restraints)!=-1].enableds = False
        ca_restraints = dr_m.atoms_restraints(atoms[atoms.names == 'CA'])
        ca_restraints[all_ss_restraints[1].indices(ca_restraints)!=-1].enableds = False

        pdr_m = sx.get_proper_dihedral_restraint_mgr(m)
        phi = pdr_m.get_restraints_by_residues_and_name(residues, 'phi')
        psi = pdr_m.get_restraints_by_residues_and_name(residues, 'psi')
        phi.enableds = False
        psi.enableds = False

    def flip_peptide_omega(self, res):
        '''
        Flip the peptide bond N-terminal to this residue from cis to
        trans or vice versa. Only usable when a simulation is running.

        Args:
            * res:
                - A :class:`chimerax.Residue` instance pointing to a
                  non-N-terminal amino acid residue.
        '''
        from . import session_extensions as sx
        pdr_m = sx.get_proper_dihedral_restraint_mgr(self.selected_model)
        omega = pdr_m.get_restraint_by_residue_and_name(res, 'omega')
        if omega is None:
            raise TypeError('This residue has no N-terminal peptide bond!')
        if omega.sim_index == -1:
            raise TypeError('Bond must be mobile in a running simulation!')
        current_angle = omega.dihedral.angle
        if abs(current_angle) <= defaults.CIS_PEPTIDE_BOND_CUTOFF:
            from math import pi
            omega.target = pi
        else:
            omega.target = 0
        omega.enabled = True

    #############################################
    # Final cleanup
    #############################################
    def _on_close(self, *_):
        self.gui_mode = False
        self.session.logger.status('Closing ISOLDE and cleaning up')

        for p in self._ui_panels:
            p.remove_trigger_handlers()

        if self.simulation_running:
            self.sim_manager.stop_sim()

        # Remove all registered event handlers
        self._event_handler.remove_all_handlers()
        self._isolde_events.remove_all_handlers()

        # Revert mouse modes
        # self._set_chimerax_mouse_mode_panel_enabled(True)
        self._mouse_modes.remove_all_modes()
        self.triggers.activate_trigger('isolde closed', None)




    ##############################################
    # General housekeeping
    ##############################################


    ##############################################
    # Help
    ##############################################

    def show_master_help_in_browser(self, *_):
        '''
        Open the ISOLDE help documentation.
        '''
        from chimerax.help_viewer import show_url
        import pathlib
        fname = os.path.join(self._root_dir, 'doc', 'index.html')
        show_url(self.session, pathlib.Path(os.path.abspath(fname)).as_uri())

    ##############################################
    # Demo
    ##############################################

    def load_model_and_data(self, model_file, mtzfile):
        '''
        Load a model and MTZ file ready for simulation

        Args:
            * model_file:
                - a PDB or mmCIF file
            * mtzfile:
                - a MTZ file containing experimental F/sigF and/or precalculated
                  F/phi columns.
        '''
        from chimerax.core.commands import open
        model = open.open(self.session, model_file)[0]
        from chimerax.std_commands import color
        color.color(self.session, model, color='bychain', target='ac')
        color.color(self.session, model, color='byhetero', target='a')
        from chimerax.clipper import symmetry
        sym_handler = symmetry.Symmetry_Manager(model,
            #mtzfile=os.path.join(data_dir, 'before_maps.mtz'),
            mtzfile=mtzfile,
            map_oversampling = self.params.map_shannon_rate)
        standard_map = sym_handler.xmapset['2mFo-DFc']
        diff_map = sym_handler.xmapset['mFo-DFc']
        from . import visualisation as v
        from chimerax.map import volumecommand
        styleargs = v.map_style_settings[v.map_styles.solid_t40]
        volumecommand.volume(self.session, [diff_map], **styleargs)
        sd = standard_map.mean_sd_rms()[1]
        standard_map.set_parameters(surface_levels = (2.5*sd,))
        diff_map.set_parameters(surface_levels = (-0.8*sd, 0.8*sd))
        from chimerax.std_commands import set
        from chimerax.core.colors import Color
        set.set(self.session, bg_color=Color([255,255,255,255]))

        self._change_selected_model(model=model, force=True)
        model.atoms[model.atoms.idatm_types != 'HC'].displays = True
        from . import view
        view.focus_on_selection(self.session, model.atoms)


    def load_demo_data(self):
        '''
        Load a small protein model with crystallographic maps to explore.
        '''
        from chimerax.core.commands import open
        data_dir = os.path.join(self._root_dir, 'demo_data', '3io0')
        before_struct = open.open(self.session, os.path.join(data_dir, 'before.pdb'))[0]
        from chimerax.std_commands import color
        color.color(self.session, before_struct, color='bychain', target='ac')
        color.color(self.session, before_struct, color='byhetero', target='a')
        from chimerax.clipper import symmetry
        sym_handler = symmetry.get_symmetry_handler(before_struct, create=True,
            auto_add_to_session=True)
        sym_handler.spotlight_radius = self.iw._sim_basic_xtal_settings_spotlight_radius_spinbox.value()
        xmapset = sym_handler.map_mgr.add_xmapset_from_mtz(os.path.join(data_dir, '3io0-sf.mtz'),
            oversampling_rate = self.params.map_shannon_rate)
        from chimerax.clipper.maps.xmapset import (map_potential_recommended_bsharp,
            viewing_recommended_bsharp)
        mdff_b = map_potential_recommended_bsharp(xmapset.resolution)
        view_b = viewing_recommended_bsharp(xmapset.resolution)
        standard_map = xmapset['(LIVE) 2mFo-DFc']
        sharp_map = xmapset['(LIVE) 2mFo-DFc_sharp_{:.0f}'.format(view_b)]
        diff_map = xmapset['(LIVE) mFo-DFc']
        from . import visualisation as v
        styleargs= v.map_style_settings[v.map_styles.solid_t20]
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [sharp_map], **styleargs)
        styleargs = v.map_style_settings[v.map_styles.solid_t40]
        volumecommand.volume(self.session, [diff_map], **styleargs)
        # standard_map.set_parameters(surface_levels = (2.5*standard_map.sigma,))
        # sharp_map.set_parameters(surface_levels = (3.0*sharp_map.sigma,))
        diff_map.set_parameters(surface_levels = (-4.0*diff_map.sigma, 4.0*diff_map.sigma))
        from chimerax.std_commands import set
        from chimerax.core.colors import Color
        set.set(self.session, bg_color=Color([255,255,255,255]))

        self._change_selected_model(model=before_struct, force=True)
        before_struct.atoms[before_struct.atoms.idatm_types != 'HC'].displays = True
        from . import view
        view.focus_on_selection(self.session, before_struct.atoms)


class Logger:
    def __init__(self, filename = None):
        self.filename = filename
        self._log_file = None
    def __call__(self, message, close = False):
        if self.filename is None:
            return    # No logging
        f = self._log_file
        if f is None:
            self._log_file = f = open(self.filename,'w')
            self._log_counter = 0
        f.write(message)
        f.write(' %d' % self._log_counter)
        f.write("\n")
        f.flush()
        self._log_counter += 1
        if close:
            f.close()
            self._log_file = None
