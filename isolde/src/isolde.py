# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

# ISOLDE: Interactive Structure Optimisation by Local Direct Exploration
# Copyright: 2016
# Author:    Tristan Croll
#            Cambridge Institute for Medical Research
#            University of Cambridge
import os
import numpy
from math import inf, degrees, radians, pi

import PyQt5
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QFileDialog

from simtk import unit

import chimerax

from chimerax import clipper

from chimerax.core import triggerset

from . import rotamers
from .eventhandler import EventHandler
from .constants import defaults, sim_outcomes, control
from .param_mgr import Param_Mgr, autodoc, param_properties
from .openmm.sim_interface import SimParams

OPENMM_LENGTH_UNIT = defaults.OPENMM_LENGTH_UNIT
OPENMM_FORCE_UNIT = defaults.OPENMM_FORCE_UNIT
OPENMM_SPRING_UNIT = defaults.OPENMM_SPRING_UNIT
OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT
OPENMM_ENERGY_UNIT = defaults.OPENMM_ENERGY_UNIT
OPENMM_ANGLE_UNIT = defaults.OPENMM_ANGLE_UNIT
CHIMERAX_LENGTH_UNIT        = defaults.CHIMERAX_LENGTH_UNIT
CHIMERAX_FORCE_UNIT         = defaults.CHIMERAX_FORCE_UNIT
CHIMERAX_SPRING_UNIT        = defaults.CHIMERAX_SPRING_UNIT

@param_properties
@autodoc
class IsoldeParams(Param_Mgr):
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
            # Provide live updates of the Ramachandran status of residues,
            # mapped to the colour of the C-alpha atoms?
        'track_ramachandran_status':            (defaults.TRACK_RAMACHANDRAN_STATUS, None),
            # Number of simulation updates to pass before updating Ramachandran
            # status
        'rounds_per_rama_update':               (defaults.ROUNDS_PER_RAMA_UPDATE, None),
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

    }

class Isolde():

    ####
    # Environment information
    ####
    import sys, os
    _root_dir = os.path.dirname(os.path.abspath(__file__))
    _platform = sys.platform

    ####
    # Enums for menu options
    ####
    from enum import IntEnum

    # Different simulation modes to set map, symmetry etc. parameters.
    class _sim_modes(IntEnum):
        xtal    = 0
        em      = 1
        free    = 2

    _force_groups = {
        'main':                         0,
        'position restraints':          1,
        }

    _human_readable_sim_modes = {
        _sim_modes.xtal:  "Crystallography mode",
        _sim_modes.em: "Single-particle EM mode",
        _sim_modes.free: "Free mode (no maps)"
        }

    # Master switch to set the level of control the user has over the simulation.
    class _experience_levels(IntEnum):
        beginner        = 0
        intermediate    = 1
        expert          = 2

    # Different modes for defining the mobile selection. If you want a
    # menu button associated with a mode, add it to the list in
    # _selection_mode_buttons
    class _sim_selection_modes(IntEnum):
        from_picked_atoms       = 0
        chain                   = 1
        whole_model             = 2
        custom                  = 3
        script                  = 99

    class _map_styles(IntEnum):
        mesh_square             = 0
        mesh_triangle           = 1
        solid_t20               = 2
        solid_t40               = 3
        solid_t60               = 4
        solid_t80               = 5
        solid_opaque            = 6

    _human_readable_map_display_styles = {
        _map_styles.mesh_square: "Mesh (squares)",
        _map_styles.mesh_triangle: "Mesh (triangles)",
        _map_styles.solid_t20: "Solid (20% opacity)",
        _map_styles.solid_t40: "Solid (40% opacity)",
        _map_styles.solid_t60: "Solid (60% opacity)",
        _map_styles.solid_t80: "Solid (80% opacity)",
        _map_styles.solid_opaque: "Solid (opaque)"
        }

    # array of settings to apply to the map depending on the chosen
    # representation. Order is: [style,
    _map_style_settings = {
        _map_styles.mesh_square: {'style': 'mesh', 'square_mesh': True, 'transparency':0},
        _map_styles.mesh_triangle: {'style': 'mesh', 'square_mesh': False, 'transparency':0},
        _map_styles.solid_t20: {'style': 'surface', 'transparency': 0.8},
        _map_styles.solid_t40: {'style': 'surface', 'transparency': 0.6},
        _map_styles.solid_t60: {'style': 'surface', 'transparency': 0.4},
        _map_styles.solid_t80: {'style': 'surface', 'transparency': 0.2},
        _map_styles.solid_opaque: {'style': 'surface', 'transparency': 0.0}
        }


    '''
    Named triggers to simplify handling of changes on key ISOLDE events,
    using the same code as ChimeraX's session.triggers. The same
    functionality could be achieved without too much trouble using the
    PyQt signal/slot system.
    '''
    trigger_names = (
        'simulation started',    # Successful initiation of a simulation
        'simulation terminated', # Simulation finished, crashed or failed to start
        'completed simulation step',
        'simulation paused',
        'simulation resumed',
        'selected model changed', # Changed the master model selection
        'position restraint added',
        'position restraint removed',
        )


    def __init__(self, gui):
        self.triggers = triggerset.TriggerSet()
        for t in self.trigger_names:
            self.triggers.add_trigger(t)

        self._isolde_events = EventHandler(self)
        self._logging = False
        self._log = Logger('isolde.log')

        self.session = gui.session

        self.params = IsoldeParams()
        self.sim_params = SimParams()

        self._status = self.session.logger.status

        #~ from .eventhandler import EventHandler
        self._event_handler = EventHandler(self.session)

        initialize_openmm()

        # Available pre-defined colors
        from chimerax.core import colors
        import copy
        self._available_colors = copy.copy(colors.BuiltinColors)
        # Remove duplicates
        for key in self._available_colors:
            if ' ' in key:
                stripped_key = key.replace(' ',"")
                self._available_colors.pop(stripped_key, None)

        # Model object to hold annotations (arrows, etc.)
        from chimerax.core.models import Model
        self._annotations = Model('ISOLDE annotations', self.session)
        self.session.models.add([self._annotations])

        ####
        # Settings for handling of atomic coordinates
        ####

        # Currently chosen mode for selecting the mobile simulation
        self._sim_selection_mode = None
        # Dict containing list of all currently loaded atomic models
        self._available_models = {}
        # Selected model on which we are actually going to run a simulation
        self._selected_model = None
        # Atoms within the model that are the primary focus of the simulation.
        # Should be whole residues
        self._selected_atoms = None
        # Extra mobile shell of surrounding atoms to provide a soft buffer to
        # the simulation. Whole residues only.
        self._soft_shell_atoms = None
        # Shell of fixed atoms surrounding all mobile atoms to maintain
        # the context of the simulation. Whole residues only.
        self._hard_shell_atoms = None
        # Construct containing all mobile atoms in the simulation
        self._total_mobile = None
        # Indices of mobile atoms in the total simulation construct
        self._total_mobile_indices = None
        # Construct containing all atoms that will actually be simulated
        self._total_sim_construct = None
        # List of all bonds in _total_sim_construct
        self._total_sim_bonds = None
        # All other atoms in the current model
        self._surroundings = None
        # Cache of atom display parameters prior to starting the simulation, so we
        # can revert once we're done.
        self._original_atom_colors = None
        self._original_atom_draw_modes = None
        self._original_bond_radii = None
        self._original_atom_radii = None
        self._original_display_state = None
        # Are the surroundings currently hidden?
        self._surroundings_hidden = False


        ####
        # Settings for live tracking of structural quality
        ####
        # TODO: Handle Ramachandran plot updating through trigger system
        # Load in Ramachandran maps
        from . import validation
        # object containing all the Ramachandran contours and lookup functions
        self._status('Preparing Ramachandran contours...')
        self.rama_validator = validation.RamaValidator()
        self._status('')
        # object that handles checking and annotation of peptide bond geometry
        self.omega_validator = validation.OmegaValidator(self._annotations)
        # Generic widget object holding the Ramachandran plot
        self._rama_plot_window = None
        # Object holding Ramachandran plot information and controls
        self._rama_plot = None
        # Is the Ramachandran plot running live?
        self._update_rama_plot = False
        # Object holding the protein phi, psi and omega dihedrals for the
        # currently selected model.
        self.backbone_dihedrals = None
        # Object holding only the backbone dihedrals that are mobile in
        # the current simulation
        self._mobile_backbone_dihedrals = None


        # Internal counter for Ramachandran update
        self._rama_counter = 0

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
        self._mouse_modes = mousemodes.MouseModeRegistry(self.session, self)
        # Placeholder for mouse tugging object
        self._mouse_tugger = None
        # Are we currently tugging an atom?
        self._currently_tugging = False

        ####
        # Haptic interfaces
        ####
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

        ####
        # Settings for handling of maps
        ####
        # List of currently available volumetric data sets
        self._available_volumes = {}
        # Master list of maps and their parameters
        self.master_map_list = {}
        # Are we adding a new map to the simulation list?
        self._add_new_map = True
        # Default coupling constants
        self.default_standard_map_k = defaults.STANDARD_MAP_K
        self.default_difference_map_k = defaults.DIFFERENCE_MAP_K

        ####
        # Restraints settings
        ####
        self.restrain_peptide_bonds = True
        self.peptide_bond_restraints_k = defaults.PEPTIDE_SPRING_CONSTANT
        self.secondary_structure_restraints_k = defaults.PHI_PSI_SPRING_CONSTANT

        # The difference between the dihedral angle and target beyond which
        # restraints begin to be applied. Below this angle there is no
        # extra biasing force. This cutoff can be changed on a per-dihedral
        # basis.
        self.default_dihedral_restraint_cutoff_angle = defaults.DIHEDRAL_RESTRAINT_CUTOFF

        ##
        # Distance retraint objects
        ##
        # (CA_n - CA_n+2) atom pairs (mainly for stretching beta strands)
        self.ca_to_ca_plus_two = None
        # (O_n - N_n+4) atom pairs (for alpha helix H-bonds)
        self.o_to_n_plus_four = None
        self.distance_restraints_k = defaults.DISTANCE_RESTRAINT_SPRING_CONSTANT

        ##
        # Position restraint objects
        ##

        # A Position_Restraints object defining all restrainable atoms
        self.position_restraints = None
        self.position_restraints_default_k = defaults.POSITION_RESTRAINT_SPRING_CONSTANT


        # A {Residue: Rotamer} dict encompassing all mobile rotameric residues
        self.rotamers = None
        self.rotamer_restraints_k = defaults.ROTAMER_SPRING_CONSTANT
        self.rotamer_restraint_cutoff_angle = defaults.ROTAMER_RESTRAINT_CUTOFF

        # Range of dihedral values which will be interpreted as a cis peptide
        # bond (-30 to 30 degrees). If restrain_peptide_bonds is True, anything
        # outside of this range at the start of the simulation will be forced
        # to trans.
        cis_offset = defaults.CIS_PEPTIDE_BOND_CUTOFF
        self.cis_peptide_bond_range = (-cis_offset, cis_offset)

        # Handler for shifting stretches in register.
        self._register_shifter = None


        ####
        # Settings for OpenMM
        ####

        # The actual simulation object
        self.sim = None
        # Placeholder for the sim_interface.SimHandler object used to perform
        # simulation setup and management
        self._sim_handler = None
        from . import sim_interface
        # List of forcefields available to the MD package
        self._available_ffs = sim_interface.available_forcefields()
        # Variables holding current forcefield choices
        self._sim_main_ff = None
        self._sim_implicit_solvent_ff = None
        self._sim_water_ff = None

        from simtk import unit
        # Simulation topology
        self._topology = None
        # OpenMM system containing the simulation
        self._system = None
        # Computational platform to run the simulation on
        self.sim_platform = None
        # Number of steps to run in before updating coordinates in ChimeraX
        self.sim_steps_per_update = defaults.SIM_STEPS_PER_GUI_UPDATE
        # Number of steps per GUI update in minimization mode
        self.min_steps_per_update = defaults.MIN_STEPS_PER_GUI_UPDATE
        # If using the VariableLangevinIntegrator, we define a tolerance
        self._integrator_tolerance = defaults.OPENMM_VAR_INTEGRATOR_TOL
        # ... otherwise, we simply set the time per step
        self._sim_time_step = defaults.OPENMM_FIXED_INTEGRATOR_TS
        # Type of integrator to use. Should give the choice in the expert level
        # of the menu. Variable is more stable, but simulated time per gui update
        # is harder to determine
        self._integrator_type = defaults.OPENMM_INTEGRATOR_TYPE
        # Constraints (e.g. rigid bonds) need their own tolerance
        self._constraint_tolerance = defaults.OPENMM_CONSTRAINT_TOL
        # Friction term for coupling to heat bath. Using a relatively high
        # value helps keep the simulation under control in response to
        # very large forces.
        self._friction = defaults.OPENMM_FRICTION
        # Limit on the net force on a single atom
        self._max_allowable_force = defaults.MAX_ALLOWABLE_FORCE # kJ mol-1 nm-1
        # For dynamics it's more efficient to just check that atoms aren't
        # moving too quickly.
        self._max_atom_movement_per_step = defaults.MAX_ATOM_MOVEMENT_PER_STEP * OPENMM_LENGTH_UNIT
        # We need to store the last measured maximum force to determine
        # when minimisation has converged.
        self._last_max_force = inf
        # Flag for unstable simulation
        self._sim_is_unstable = False
        # Counter for the number of simulation rounds that have been unstable
        self._unstable_min_rounds = 0
        # Maximum number of rounds to attempt minimisation of an unstable
        # simulation
        self.max_unstable_rounds = defaults.MAX_UNSTABLE_ROUNDS

        # Are we currently tugging on an atom?
        self._currently_tugging = False
        # Placeholder for tugging forces
        self._tugging_force = None
        # Force constant for mouse/haptic tugging. Need to make this user-adjustable
        self.tug_force_constant = defaults.MOUSE_TUG_SPRING_CONSTANT
        # Upper limit on the strength of the tugging force
        self.tug_max_force = defaults.MAX_TUG_FORCE # kJ/mol/nm



        # Holds the current simulation mode, to be chosen from the GUI
        # drop-down menu or loaded from saved settings.
        self.sim_mode = None
        # Do we have a simulation running right now?
        self._simulation_running = False
        # If running, is the simulation in startup mode?
        self._sim_startup = True
        # Maximum number of rounds of minimisation to run on startup
        self._sim_startup_rounds = defaults.SIM_STARTUP_ROUNDS
        # Counter for how many rounds we've done on startup
        self._sim_startup_counter = 0

        # Simulation temperature in Kelvin
        self.simulation_temperature = defaults.TEMPERATURE
        # Flag to update the temperature of a running simulation
        self._temperature_changed = False

        # If a simulation is running, is it paused?
        self._sim_paused = False

        # Are we equilibrating or minimising?
        self.simulation_type = 'equil'

        # Current positions of all particles in the simulation
        self._particle_positions = None
        # Saved particle positions in case we want to discard the simulation
        self._saved_positions = None

        # To ensure we do only one simulation round per graphics redraw
        self._last_frame_number = None

        self.tug_hydrogens = False
        self.hydrogens_feel_maps = False


        self.initialize_haptics()

        ####
        # Internal trigger handlers
        ####


        # During simulations, ISOLDE needs to reserve the right mouse
        # button for tugging atoms. Therefore, the ChimeraX right mouse
        # mode selection panel needs to be disabled for the duration of
        # each simulation.
        ie = self._isolde_events
        ie.add_event_handler('disable chimerax mouse mode panel',
                                              'simulation started',
                                              self._disable_chimerax_mouse_mode_panel)
        ie.add_event_handler('enable chimerax mouse mode panel',
                                              'simulation terminated',
                                              self._enable_chimerax_mouse_mode_panel)
        ie.add_event_handler('on sim start',
                                              'simulation started',
                                              self._sim_start_cb)
        ie.add_event_handler('cleanup after sim',
                                              'simulation terminated',
                                              self._sim_end_cb)
        ie.add_event_handler('sim pause', 'simulation paused', self._sim_pause_cb)
        ie.add_event_handler('sim resume', 'simulation resumed', self._sim_resume_cb)

        self.gui_mode = False

        from PyQt5.QtGui import QPixmap
        from PyQt5.QtWidgets import QSplashScreen
        from PyQt5.QtCore import Qt
        import os

        splash_pix = QPixmap(os.path.join(
            self._root_dir,'resources/isolde_splash_screen.jpg'))
        splash = self._splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.showMessage('\n\n\n\nAnalysing your structure. Please be patient...',
                alignment = Qt.AlignRight | Qt.AlignVCenter)
        splash.show()
        # Make sure the splash screen is actually shown
        for i in range(5):
            self.session.ui.processEvents()
        self._splash_destroy_countdown = 100
        self._splash_handler = self.session.triggers.add_handler('new frame',
            self._splash_destroy_cb)

        self.start_gui(gui)




    @property
    def selected_model(self):
        return self._selected_model

    @selected_model.setter
    def selected_model(self, model):
        if not isinstance(model, chimerax.core.atomic.AtomicStructure):
            raise TypeError('Selection must be a single AtomicStructure model!')
        self._change_selected_model(model = model)

    @property
    def fixed_atoms(self):
        return self._hard_shell_atoms
    @property
    def mobile_atoms(self):
        return self._total_mobile
    @property
    def all_simulated_atoms(self):
        return self._total_sim_construct



    ###################################################################
    # GUI related functions
    ###################################################################

    def start_gui(self, gui):
        ####
        # Connect and initialise ISOLDE widget
        ####

        self.gui = gui
        self.iw = gui.iw
        self.gui_mode = True

         # Register ISOLDE-specific mouse modes
        self._mouse_modes.register_all_isolde_modes()



        # Function to remove all event handlers, mouse modes etc. when
        # ISOLDE is closed, and return ChimeraX to its standard state.

        #self.gui.tool_window.ui_area.parentWidget().destroyed.connect(self._on_close)
        #self.gui.tool_window.destroyed.connect(self._on_close)
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


        self._event_handler.add_event_handler('update_menu_on_selection',
                                              'selection changed',
                                              self._selection_changed)
        self._event_handler.add_event_handler('update_menu_on_model_add',
                                              'add models',
                                              self._update_model_list)
        self._event_handler.add_event_handler('update_menu_on_model_remove',
                                              'remove models',
                                              self._update_model_list)
        self._selection_changed()
        self._update_model_list()

        # Work out menu state based on current ChimeraX session
        self._update_sim_control_button_states()




    def _populate_menus_and_update_params(self):
        iw = self.iw
        # Clear experience mode placeholders from QT Designer and repopulate
        cb = iw._experience_level_combo_box
        cb.clear()
        for lvl in self._experience_levels:
            cb.addItem(lvl.name, lvl)

        # Clear simulation mode placeholders from QT Designer and repopulate
        cb = iw._sim_basic_mode_combo_box
        cb.clear()
        for mode in self._sim_modes:
            text = self._human_readable_sim_modes[mode]
            cb.addItem(text, mode)

        iw._sim_temp_spin_box.setProperty('value', self.simulation_temperature)

        # Populate force field combo box with available forcefields
        cb = iw._sim_force_field_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.main_file_descriptions)

        # Populate water model combo box with available models
        cb = iw._sim_water_model_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.explicit_water_descriptions)

        # Populate OpenMM platform combo box with available platforms
        cb = iw._sim_platform_combo_box
        cb.clear()
        from . import sim_interface as si
        platform_names = si.get_available_platforms()
        cb.addItems(platform_names)

        # Set to the fastest available platform
        if 'CUDA' in platform_names:
            cb.setCurrentIndex(platform_names.index('CUDA'))
        elif 'OpenCL' in platform_names:
            cb.setCurrentIndex(platform_names.index('OpenCL'))
        elif 'CPU' in platform_names:
            cb.setCurrentIndex(platform_names.index('CPU'))

        # The last entry in the EM map chooser combo box should always be "Add map"
        cb = iw._em_map_chooser_combo_box
        cb.clear()
        cb.addItem('Add map')
        cb.setCurrentText('Add map')

        # Map display style options
        cb = iw._em_map_style_combo_box
        cb.clear()
        for mode in self._map_styles:
            text = self._human_readable_map_display_styles[mode]
            cb.addItem(text, mode)
        cb.setCurrentIndex(-1)

        cb = iw._em_map_color_combo_box
        cb.clear()
        for key, cval in self._available_colors.items():
            cb.addItem(key, cval)
        cb.setCurrentIndex(-1)

        cb = iw._em_map_contour_units_combo_box
        cb.clear()
        cb.addItem('sigma')
        cb.addItem('map units')
        cb.setCurrentIndex(0)

        ####
        # Rebuild tab
        ####

        ## Info for a single selected residue
        iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        iw._rebuild_sel_res_pep_info.setText('')
        iw._rebuild_sel_res_rot_info.setText('')

        from . import dihedrals
        phipsi = dihedrals.Backbone_Dihedrals.standard_phi_psi_angles
        iw._rebuild_2ry_struct_restr_chooser_combo_box.clear()
        for key, pair in phipsi.items():
            iw._rebuild_2ry_struct_restr_chooser_combo_box.addItem(key, pair)

        ####
        # Validate tab
        ####

        # Populate the Ramachandran plot case selector with available
        # cases
        cb = iw._validate_rama_case_combo_box
        cb.clear()
        from . import validation
        # First two keys are N- and C-terminal residues, which we don't plot
        keys = validation.RAMA_CASES[2:]
        for key in reversed(keys):
            cb.addItem(validation.RAMA_CASE_DETAILS[key]['name'], key)


    def _prepare_ramachandran_plot(self):
        '''
        Prepare an empty MatPlotLib figure to put the Ramachandran plots in.
        '''
        from . import validation
        iw = self.iw
        container = self._rama_plot_window = iw._validate_rama_plot_layout
        self._rama_plot = validation.RamaPlot(self.session, container, self.rama_validator)




    def _connect_functions(self):
        '''
        Connect PyQt events from the ISOLDE gui widget to functions.
        '''
        iw = self.iw
        ####
        # Master switches
        ####

        iw._experience_level_combo_box.currentIndexChanged.connect(
            self.gui._change_experience_level_or_sim_mode
            )
        ## Initialise to current level
        #self._change_experience_level_or_sim_mode()

        iw._sim_basic_mode_combo_box.currentIndexChanged.connect(
            self._change_sim_mode
            )
        # Initialise to selected mode.
        self._change_sim_mode()
        # automatically runs self.gui._change_experience_level_or_sim_mode()

        for button in self.gui._selection_mode_buttons:
            button.clicked.connect(self._change_sim_selection_mode)

        self._change_sim_selection_mode()


        ####
        # Simulation global parameters (can only be set before running)
        ####
        iw._sim_force_field_combo_box.currentIndexChanged.connect(
            self._change_force_field
            )
        iw._sim_water_model_combo_box.currentIndexChanged.connect(
            self._change_water_model
            )
        iw._master_model_combo_box.currentIndexChanged.connect(
            self._change_selected_model
            )
        iw._sim_basic_mobile_chains_list_box.itemSelectionChanged.connect(
            self._change_selected_chains
            )
        iw._sim_basic_mobile_sel_within_spinbox.valueChanged.connect(
            self._change_soft_shell_cutoff
            )
        iw._sim_basic_mobile_b_and_a_spinbox.valueChanged.connect(
            self._change_b_and_a_padding
            )
        iw._sim_basic_mobile_sel_backbone_checkbox.stateChanged.connect(
            self._change_soft_shell_fix_backbone
            )
        iw._sim_platform_combo_box.currentIndexChanged.connect(
            self._change_sim_platform
            )
        iw._sim_basic_whole_structure_button.clicked.connect(
            self._select_whole_model
            )

        # Run all connected functions once to initialise
        self._change_force_field()
        self._change_water_model()
        self._change_selected_model()
        self._change_selected_chains()
        self._change_soft_shell_cutoff()
        self._change_b_and_a_padding()
        self._change_soft_shell_fix_backbone()
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
        # Xtal map parameters (can only be set before starting simulation)
        ####
        iw._sim_basic_xtal_init_open_button.clicked.connect(
            self._show_xtal_init_frame
            )
        iw._sim_basic_xtal_init_done_button.clicked.connect(
            self._hide_xtal_init_frame
            )
        iw._sim_basic_xtal_init_model_combo_box.currentIndexChanged.connect(
            self._check_for_valid_xtal_init
            )
        iw._sim_basic_xtal_init_reflections_file_button.clicked.connect(
            self._choose_mtz_file
            )
        iw._sim_basic_xtal_init_go_button.clicked.connect(
            self._initialize_xtal_structure
            )
        iw._sim_basic_xtal_map_settings_show_button.clicked.connect(
            self._toggle_xtal_map_dialog
            )
        iw._sim_basic_xtal_settings_map_combo_box.currentIndexChanged.connect(
            self._populate_xtal_map_params
            )
        iw._sim_basic_xtal_settings_set_button.clicked.connect(
            self._apply_xtal_map_params
            )

        ####
        # EM map parameters (can only be set before starting simulation)
        ####
        iw._sim_basic_em_map_button.clicked.connect(
            self._show_em_map_chooser
            )
        iw._em_map_done_button.clicked.connect(
            self._hide_em_map_chooser
            )
        iw._em_map_set_button.clicked.connect(
            self._add_or_change_em_map_from_gui
            )
        iw._em_map_remove_button.clicked.connect(
            self._remove_em_map_from_gui
            )
        iw._em_map_chooser_combo_box.currentIndexChanged.connect(
            self._show_em_map_in_menu_or_add_new
            )
        iw._em_map_style_combo_box.currentIndexChanged.connect(
            self._change_display_of_selected_map
            )
        iw._em_map_color_combo_box.currentIndexChanged.connect(
            self._change_display_of_selected_map
            )
        iw._em_map_contour_spin_box.valueChanged.connect(
            self._change_contour_level_of_selected_map
            )
        iw._em_map_contour_units_combo_box.currentIndexChanged.connect(
            self._change_contour_units_of_selected_map
            )
         # We want to start with the EM map chooser hidden
        self._hide_em_map_chooser()

        ####
        # Restraints tab
        ####
        iw._restraints_pep_go_button.clicked.connect(
            self._change_peptide_bond_restraints
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

        iw._rebuild_2ry_struct_restr_show_button.clicked.connect(
            self._toggle_secondary_structure_dialog
            )
        iw._rebuild_2ry_struct_restr_chooser_go_button.clicked.connect(
            self._apply_selected_secondary_structure_restraints
            )
        iw._rebuild_2ry_struct_restr_clear_button.clicked.connect(
            self.clear_secondary_structure_restraints_for_selection
            )

        iw._rebuild_register_shift_dialog_toggle_button.clicked.connect(
            self._toggle_register_shift_dialog
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


        iw._rebuild_pos_restraint_dialog_toggle_button.clicked.connect(
            self._toggle_position_restraints_dialog
            )
        iw._rebuild_pos_restraint_go_button.clicked.connect(
            self._restrain_selected_atom_to_xyz
            )
        iw._rebuild_pos_restraint_clear_button.clicked.connect(
            self.release_xyz_restraints_on_selected_atoms
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
        iw._validate_rama_case_combo_box.currentIndexChanged.connect(
            self._change_rama_case
            )
        iw._validate_rama_go_button.clicked.connect(
            self._rama_static_plot
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
        iw._validate_pep_cis_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )
        iw._validate_pep_twisted_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )


        ####
        # Simulation control functions
        ####

        iw._sim_temp_spin_box.valueChanged.connect(
            self._update_sim_temperature
            )
        iw._sim_go_button.clicked.connect(
            self.start_sim
            )
        iw._sim_pause_button.clicked.connect(
            self.pause_sim_toggle
            )
        iw._sim_commit_button.clicked.connect(
            self.commit_sim
            )
        iw._sim_discard_button.clicked.connect(
            self.discard_sim
            )
        iw._sim_min_button.clicked.connect(
            self.minimize
            )
        iw._sim_equil_button.clicked.connect(
            self.equilibrate
        )
        iw._sim_hide_surroundings_toggle.stateChanged.connect(
            self._set_hide_surroundings
        )

    def _disable_chimerax_mouse_mode_panel(self, *_):
        self._set_chimerax_mouse_mode_panel_enabled(False)

    def _enable_chimerax_mouse_mode_panel(self, *_):
        self._set_chimerax_mouse_mode_panel_enabled(True)

    def _set_chimerax_mouse_mode_panel_enabled(self, state):
        from chimerax.mouse_modes.tool import MouseModePanel
        mm = self.session.tools.find_by_class(MouseModePanel)[0]
        mm.buttons.setEnabled(state)


    def initialize_haptics(self):
        '''
        If the HapticHandler plugin is installed, start it and create
        a HapticTugger object for each connected haptic device.
        '''
        if hasattr(self.session, 'HapticHandler'):
            self._status('Initialising haptic interface(s)')
            hh = self.session.HapticHandler
            hh.startHaptics()
            n = self._num_haptic_devices = hh.getNumDevices()
            d = self._haptic_devices = [None] * n
            self._haptic_tug_atom = [None] * n
            self._haptic_tug_index = [-1] * n
            self._haptic_highlight_nearest_atom = [True] * n
            self._haptic_current_nearest_atom = [None] * n
            self._haptic_current_nearest_atom_color = [None] * n
            self._haptic_current_nearest_atom_draw_mode = [None] * n
            from . import haptics
            for i in range(n):
                d[i] = haptics.HapticTugger(self.session, i, self._annotations)
            self._use_haptics = True
            self._status('')
        else:
            self._use_haptics = False

    def _update_sim_temperature(self):
        t = self.iw._sim_temp_spin_box.value()
        self.simulation_temperature = t
        # So we know to update the temperature in any running simulation
        self._temperature_changed = True

    ##############################################################
    # Menu control functions to run on key events
    ##############################################################

    def _update_model_list(self, *_):
        self.iw._master_model_combo_box.clear()
        self.iw._em_map_model_combo_box.clear()
        self.iw._sim_basic_xtal_init_model_combo_box.clear()
        models = self.session.models.list()
        atomic_model_list = []
        atomic_model_name_list = []
        potential_xtal_list = []
        potential_xtal_name_list = []
        volume_model_list = []
        volume_model_name_list = []
        sorted_models = sorted(models, key=lambda m: m.id)
        if len(sorted_models) != 0:
            # Find atomic and volumetric models and sort them into the
            # appropriate lists
            for i, m in enumerate(sorted_models):
                if m.atomspec_has_atoms():
                    # Ignore all symmetry-equivalent copies
                    if m.parent.name == 'symmetry equivalents':
                        continue
                    id_str = m.id_string() + ' ' + m.name
                    if type(m.parent) == clipper.CrystalStructure:
                        if self.sim_mode == self._sim_modes.em:
                            continue
                    elif self.sim_mode == self._sim_modes.xtal:
                        # Add it to the list of potential models to turn into crystal structures
                        potential_xtal_name_list.append(id_str)
                        potential_xtal_list.append(m)
                        continue
                    self._available_models[id_str] = m
                    atomic_model_name_list.append(id_str)
                    atomic_model_list.append(m)
            for i, m in enumerate(sorted_models):
                if hasattr(m, 'grid_data'):
                    # This menu is only for real-space maps. Ignore clipper maps
                    if type(m) == clipper.crystal.XmapHandler:
                        continue
                    id_str = m.id_string() + ' ' + m.name
                    self._available_volumes[id_str] = m
                    volume_model_name_list.append(id_str)
                    volume_model_list.append(m)
                else:
                    # This is a model type we don't currently handle. Ignore.
                    continue
        for l, m in zip(atomic_model_name_list, atomic_model_list):
            self.iw._master_model_combo_box.addItem(l, m)
        for l, m in zip(potential_xtal_name_list, potential_xtal_list):
            self.iw._sim_basic_xtal_init_model_combo_box.addItem(l, m)
        for l, m in zip(volume_model_name_list, volume_model_list):
            self.iw._em_map_model_combo_box.addItem(l, m)

    def _update_chain_list(self):
        m = self._selected_model
        chains = m.chains.chain_ids
        lb = iw = self.iw._sim_basic_mobile_chains_list_box
        lb.clear()
        lb.addItems(chains)


    def _selection_changed(self, *_):
        from chimerax.core.atomic import selected_atoms
        from .util import is_continuous_protein_chain
        sel = selected_atoms(self.session)
        selres = sel.unique_residues
        if self._simulation_running:
            natoms = len(sel)
            if natoms == 1:
                self._enable_atom_position_restraints_frame()
            else:
                self._disable_atom_position_restraints_frame()
            if natoms:
                self._enable_position_restraints_clear_button()
            else:
                self._disable_position_restraints_clear_button()

            if len(selres) == 1:
                self._enable_rebuild_residue_frame(selres[0])
            else:
                self._disable_rebuild_residue_frame()
            if is_continuous_protein_chain(sel):
                self._enable_secondary_structure_restraints_frame()
                self._enable_register_shift_frame()
            else:
                self._disable_secondary_structure_restraints_frame()
                self._disable_register_shift_frame()

            # A running simulation takes precedence for memory control
            return
        flag = not(self.session.selection.empty())
        iw = self.iw
        iw._sim_basic_mobile_selection_frame.setEnabled(flag)
        iw._sim_go_button.setEnabled(flag)


    def _update_sim_control_button_states(self):
        # Set enabled/disabled states of main simulation control panel
        # based on whether a simulation is currently running
        flag = self._simulation_running
        iw = self.iw
        iw._sim_go_button.setDisabled(flag)
        if self._sim_paused and not flag:
            self._sim_paused = False
            iw._sim_pause_button.setText('Pause')
        iw._sim_pause_button.setEnabled(flag)
        iw._sim_commit_button.setEnabled(flag)
        iw._sim_discard_button.setEnabled(flag)
        # Change colour of minimisation and equilibration buttons according
        # to current choice
        if self.simulation_type == 'equil':
            iw._sim_equil_button.setStyleSheet('background-color: green')
            iw._sim_min_button.setStyleSheet('background-color: red')
        else:
            iw._sim_equil_button.setStyleSheet('background-color: red')
            iw._sim_min_button.setStyleSheet('background-color: green')

        # Undo/redo will only lead to trouble while a simulation is running
        iw._undo_button.setDisabled(flag)
        iw._redo_button.setDisabled(flag)

        # Update the status of the Go button
        self._selection_changed()

    def _change_sim_selection_mode(self):
        iw = self.iw
        iw._sim_basic_mobile_selection_frame.hide()
        iw._sim_basic_mobile_by_chain_frame.hide()
        iw._sim_basic_mobile_custom_frame.hide()

        for i, b in enumerate(self.gui._selection_mode_buttons):
            if b.isChecked():
                break

        if i == 0:
            self._sim_selection_mode = self._sim_selection_modes.from_picked_atoms
            iw._sim_basic_mobile_selection_frame.show()
        elif i == 1:
            self._sim_selection_mode = self._sim_selection_modes.chain
            iw._sim_basic_mobile_by_chain_frame.show()
            self._change_selected_chains()
        elif i == 2:
            self._sim_selection_mode = self._sim_selection_modes.whole_model
        elif i == 3:
            self._sim_selection_mode = self._sim_selection_modes.custom
            iw._sim_basic_mobile_custom_frame.show()
        else:
            raise Exception('No or unrecognised mode selected!')

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
        cb = self.iw._sim_basic_xtal_init_model_combo_box
        if cb.currentData() is not None:
            if os.path.isfile(self.iw._sim_basic_xtal_init_reflections_file_name.text()):
                self.iw._sim_basic_xtal_init_go_button.setEnabled(True)
                return
        self.iw._sim_basic_xtal_init_go_button.setEnabled(False)
        self.iw._sim_basic_xtal_init_reflections_file_name.setText('')

    def _choose_mtz_file(self, *_):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        caption = 'Choose a file containing map structure factors'
        filetypes = 'MTZ files (*.mtz)'
        filename, _ = QFileDialog.getOpenFileName(None, caption, filetypes, filetypes, options = options)
        self.iw._sim_basic_xtal_init_reflections_file_name.setText(filename)
        self._check_for_valid_xtal_init()


    def _initialize_xtal_structure(self, *_):
        cb = self.iw._sim_basic_xtal_init_model_combo_box
        fname = self.iw._sim_basic_xtal_init_reflections_file_name.text()
        if not cb.count:
            errstring = 'No atomic structures are available that are not \
                already part of an existing crystal structure. Please load \
                one first.'
            _generic_waring(errstring)
        if not os.path.isfile(fname):
            errstring = 'Please select a valid MTZ file!'
            _generic_warning(errstring)
        m = cb.currentData()
        clipper.CrystalStructure(self.session, m, fname)
        self.iw._sim_basic_xtal_init_reflections_file_name.setText('')
        self.iw._sim_basic_xtal_init_go_button.setEnabled(False)
        self._change_selected_model(force = True)

    def _initialize_xtal_maps(self, model):
        '''
        Set up all crystallographic maps associated with a model, ready
        for simulation.
        '''
        self.master_map_list.clear()
        xmaps = model.xmaps.child_models()
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        cb.clear()
        for xmap in xmaps:
            name = xmap.name
            is_difference_map = xmap.is_difference_map
            if is_difference_map:
                cutoff = self.params.difference_map_mask_cutoff
                coupling_constant = self.default_difference_map_k
            else:
                cutoff = self.params.standard_map_mask_cutoff
                coupling_constant = self.default_standard_map_k
            new_map = self.add_map(name, xmap, cutoff, coupling_constant,
                is_difference_map = is_difference_map, mask = True,
                crop = False)
            cb.addItem(name, new_map)

    def _toggle_xtal_map_dialog(self, *_):
        button = self.iw._sim_basic_xtal_map_settings_show_button
        frame = self.iw._sim_basic_xtal_map_settings_frame
        show_text = 'Show map settings dialogue'
        hide_text = 'Hide map settings dialogue'
        if button.text() == show_text:
            frame.show()
            button.setText(hide_text)
        else:
            frame.hide()
            button.setText(show_text)

    def _populate_xtal_map_params(self, *_):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        tb = self.iw._sim_basic_xtal_settings_map_name
        this_map = cb.currentData()
        if cb.currentIndex() == -1 or this_map is None:
            tb.setText('No maps loaded!')
            return
        self.iw._sim_basic_xtal_map_cutoff_spin_box.setValue(
            this_map.get_mask_cutoff())
        self.iw._sim_basic_xtal_map_weight_spin_box.setValue(
            this_map.get_coupling_constant())
        self.iw._sim_basic_xtal_settings_map_masked_checkbox.setCheckState(
            this_map.get_mask_vis())

    def _apply_xtal_map_params(self, *_):
        cb = self.iw._sim_basic_xtal_settings_map_combo_box
        this_map = cb.currentData()
        this_map.set_mask_cutoff(
            self.iw._sim_basic_xtal_map_cutoff_spin_box.value())
        this_map.set_coupling_constant(
            self.iw._sim_basic_xtal_map_weight_spin_box.value())
        this_map.set_mask_vis(
            self.iw._sim_basic_xtal_settings_map_masked_checkbox.checkState())


    ####
    # EM
    ####

    def _show_em_map_chooser(self, *_):
        self.iw._em_map_chooser_frame.show()
        self.iw._sim_basic_em_map_button.setEnabled(False)

    def _hide_em_map_chooser(self, *_):
        self.iw._em_map_chooser_frame.hide()
        self.iw._sim_basic_em_map_button.setEnabled(True)

    def _show_em_map_in_menu_or_add_new(self, *_):
        if self.sim_mode != self._sim_modes.em:
            return
        iw = self.iw
        seltext = iw._em_map_chooser_combo_box.currentText()
        if seltext == 'Add map':
            self._add_new_map = True
            iw._em_map_name_field.setText('')
            iw._em_map_model_combo_box.setCurrentIndex(-1)
            iw._em_map_name_field.setEnabled(True)
            iw._em_map_style_combo_box.setCurrentIndex(-1)
            iw._em_map_color_combo_box.setCurrentIndex(-1)
        elif len(seltext):
            self._add_new_map = False
            current_map = self.master_map_list[seltext]
            name, vol, cutoff, coupling, style, color, contour, contour_units, \
                mask_vis, is_per_atom, per_atom_k = current_map.get_map_parameters()
            iw._em_map_name_field.setText(name)
            id_str = vol.id_string() + ' ' + vol.name
            iw._em_map_model_combo_box.setCurrentText(id_str)
            iw._em_map_cutoff_spin_box.setValue(cutoff)
            iw._em_map_coupling_spin_box.setValue(coupling)
            if style is not None:
                iw._em_map_style_combo_box.setCurrentText(style)
            else:
                iw._em_map_style_combo_box.setCurrentIndex(-1)
            if color is not None:
                iw._em_map_color_combo_box.setCurrentText(color)
            else:
                iw._em_map_color_combo_box.setCurrentIndex(-1)
            if contour is None:
                contour = vol.surface_levels[0]
                if iw._em_map_contour_units_combo_box.currentText() == 'sigma':
                    sigma = vol.mean_sd_rms()[1]
                    contour = contour/sigma
                iw._em_map_contour_spin_box.setValue(contour)

            if contour_units is not None:
                iw._em_map_contour_units_combo_box.setCurrentText(contour_units)

            iw._em_map_masked_checkbox.setCheckState(mask_vis)
            # Map name is our lookup key, so can't change it after it's been added
            iw._em_map_name_field.setEnabled(False)
        else:
            # Do nothing. List is empty.
            return

    def _update_master_map_list_combo_box(self):
        cb = self.iw._em_map_chooser_combo_box
        cb.clear()
        keylist = sorted([key for key in self.master_map_list])
        cb.addItems(keylist)
        cb.addItem('Add map')


    # Update button states after a simulation has finished
    def _update_menu_after_sim(self):
        self._update_sim_control_button_states()





    ####
    # Rebuild tab
    ####
    def _enable_rebuild_residue_frame(self, res):
        if not self._simulation_running:
            self._disable_rebuild_residue_frame()
            return
        if -1 in self._total_mobile.indices(res.atoms):
            self._disable_rebuild_residue_frame()
            return

        # Peptide cis/trans flips
        self._get_rotamer_list_for_selected_residue(res)
        bd = self._mobile_backbone_dihedrals
        self._rebuild_residue = res
        omega = bd.omega.by_residue(res)
        if omega is None:
            # This is a terminal residue without an omega dihedral
            self.iw._rebuild_sel_res_pep_info.setText('N/A')
            self.iw._rebuild_sel_res_cis_trans_flip_button.setDisabled(True)
            self.iw._rebuild_sel_res_pep_flip_button.setDisabled(True)
            self._rebuild_res_update_omega = False
            return
        self._rebuild_res_update_omega = True
        self._rebuild_res_omega = omega
        self.iw._rebuild_sel_res_cis_trans_flip_button.setEnabled(True)
        self.iw._rebuild_sel_res_pep_flip_button.setEnabled(True)

        try:
            rot = self._selected_rotamer = self.rotamers[res]
            if rot != self._selected_rotamer:
                self.iw._rebuild_sel_res_rot_info.setText('')
            self._selected_rotamer = rot
            self._set_rotamer_buttons_enabled(True)
        except KeyError:
            # This residue has no rotamers
            self._selected_rotamer = None
            self._set_rotamer_buttons_enabled(False)

        self.iw._rebuild_sel_residue_frame.setEnabled(True)
        chain_id, resnum, resname = (res.chain_id, res.number, res.name)
        self.iw._rebuild_sel_residue_info.setText(str(chain_id) + ' ' + str(resnum) + ' ' + resname)

        self._steps_per_sel_res_update = 10
        self._sel_res_update_counter = 0
        if 'update_selected_residue_info' not in self._event_handler.list_event_handlers():
            self._event_handler.add_event_handler('update_selected_residue_info',
                    'shape changed', self._update_selected_residue_info_live)

    def _enable_atom_position_restraints_frame(self):
        self.iw._rebuild_pos_restraint_one_atom_frame.setEnabled(True)

    def _disable_atom_position_restraints_frame(self):
        self.iw._rebuild_pos_restraint_one_atom_frame.setEnabled(False)

    def _enable_position_restraints_clear_button(self):
        self.iw._rebuild_pos_restraint_clear_button.setEnabled(True)

    def _disable_position_restraints_clear_button(self):
        self.iw._rebuild_pos_restraint_clear_button.setEnabled(False)


    def _enable_secondary_structure_restraints_frame(self, *_):
        self.iw._rebuild_2ry_struct_restr_container.setEnabled(True)
        self.iw._rebuild_2ry_struct_restr_top_label.setText('')

    def _enable_register_shift_frame(self, *_):
        self.iw._rebuild_register_shift_container.setEnabled(True)

    def _disable_secondary_structure_restraints_frame(self, *_):
        self.iw._rebuild_2ry_struct_restr_container.setEnabled(False)
        t = 'Select a protein atom or residue to begin'
        self.iw._rebuild_2ry_struct_restr_top_label.setText(t)
        t = 'Invalid selection'
        self.iw._rebuild_2ry_struct_restr_sel_text.setText(t)

    def _disable_register_shift_frame(self, *_):
        self.iw._rebuild_register_shift_container.setEnabled(False)


    def _toggle_secondary_structure_dialog(self, *_):
        button = self.iw._rebuild_2ry_struct_restr_show_button
        frame = self.iw._rebuild_2ry_struct_restr_container
        show_text = 'Show secondary structure dialogue'
        hide_text = 'Hide secondary structure dialogue'
        if button.text() == show_text:
            frame.show()
            button.setText(hide_text)
        else:
            frame.hide()
            button.setText(show_text)

    def _toggle_register_shift_dialog(self, *_):
        button = self.iw._rebuild_register_shift_dialog_toggle_button
        frame = self.iw._rebuild_register_shift_container
        show_text = 'Show register shift dialogue'
        hide_text = 'Hide register shift dialogue'
        if button.text() == show_text:
            frame.show()
            button.setText(hide_text)
        else:
            frame.hide()
            button.setText(show_text)

    def _toggle_position_restraints_dialog(self, *_):
        button = self.iw._rebuild_pos_restraint_dialog_toggle_button
        frame = self.iw._rebuild_pos_restraint_one_atom_frame
        show_text = 'Show position restraints dialogue'
        hide_text = 'Hide position restraints dialogue'
        if button.text() == show_text:
            frame.show()
            button.setText(hide_text)
        else:
            frame.hide()
            button.setText(show_text)



    def _extend_selection_by_one_res_N(self, *_):
        self._extend_selection_by_one_res(-1)

    def _shrink_selection_by_one_res_N(self, *_):
        self._shrink_selection_by_one_res(1)

    def _extend_selection_by_one_res_C(self, *_):
        self._extend_selection_by_one_res(1)

    def _shrink_selection_by_one_res_C(self, *_):
        self._shrink_selection_by_one_res(-1)

    def _extend_selection_by_one_res(self, direction):
        '''
        Extends a selection by one residue in the given direction, stopping
        at the ends of the chain
        Args:
            direction:
                -1 or 1
        '''
        from chimerax.core.atomic import selected_atoms
        sel = selected_atoms(self.session)
        residues = sel.unique_residues
        # expand the selection to cover all atoms in the existing residues
        residues.atoms.selected = True
        m = sel.unique_structures[0]
        polymers = m.polymers(
            missing_structure_treatment = m.PMS_NEVER_CONNECTS)
        for p in polymers:
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
        from chimerax.core.atomic import selected_atoms
        sel = selected_atoms(self.session)
        residues = sel.unique_residues
        if len(residues) > 1:
            if direction == 1:
                residues[0].atoms.selected = False
            elif direction == -1:
                residues[-1].atoms.selected = False
            else:
                raise TypeError('Direction must be either 1 or -1!')


        #~ seltext = 'Chain {}: residues {} to {}'.format(
            #~ p.unique_chain_ids[0], first.number, last.number)
        #~ self.iw._rebuild_2ry_struct_restr_sel_text.setText(seltext)


    def _apply_selected_secondary_structure_restraints(self, *_):
        from chimerax.core.atomic import selected_atoms
        sh = self._sim_handler
        sc = self._total_sim_construct
        dihed_k = self.secondary_structure_restraints_k
        sel = selected_atoms(self.session)
        residues = sel.unique_residues
        sel = residues.atoms
        cb = self.iw._rebuild_2ry_struct_restr_chooser_combo_box
        structure_type = cb.currentText()
        target_phi, target_psi = cb.currentData()
        phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
        phi.targets = target_phi
        psi.targets = target_psi
        phi.spring_constants = dihed_k
        psi.spring_constants = dihed_k
        self.apply_dihedral_restraints(phi)
        self.apply_dihedral_restraints(psi)


        dist_k = self.distance_restraints_k
        rca = self._sim_distance_restraint_dict['ca_to_ca_plus_two']
        ron = self._sim_distance_restraint_dict['o_to_n_plus_four']
        if structure_type == 'helix':
            ca_distance = rca.HELIX_DISTANCE
            on_distance = ron.HELIX_DISTANCE
        elif structure_type in ('parallel beta', 'antiparallel beta'):
            ca_distance = rca.STRAND_DISTANCE
            on_distance = ron.STRAND_DISTANCE
        else:
            raise TypeError('Unknown secondary structure type: {}!'.format(structure_type))
        for r in residues:
            cad_candidates = rca[r]
            #print('Residue number: {}, cad residue numbers: {}'.format(r.number, cad.atoms.residues.numbers))
            for cad in cad_candidates:
                if cad is not None:
                    if -1 not in sel.indices(cad.atoms):
                        print('setting CA-CA distance to {}'.format(ca_distance))
                        cad.target_distance = ca_distance
                        cad.spring_constant = dist_k
                    else:
                        cad.spring_constant = 0
                    self.apply_distance_restraint(cad)
            ond_candidates = ron[r]
            for ond in ond_candidates:
                if ond is not None:
                    if -1 not in sel.indices(ond.atoms):
                        ond.target_distance = on_distance
                        print('setting O-N distance to {}'.format(on_distance))
                        ond.spring_constant = dist_k
                    else:
                        ond.spring_constant = 0
                    self.apply_distance_restraint(ond)

    def clear_secondary_structure_restraints_for_selection(self, *_, atoms = None, residues = None):
        '''
        Clear all secondary structure restraints selection. If no atoms
        or residues are provided, restraints will be cleared for any
        atoms selected in the main window.
        '''
        from chimerax.core.atomic import selected_atoms
        si = self._sim_interface
        sc = self._total_sim_construct
        sel = None
        if atoms is not None:
            if residues is not None:
                raise TypeError('Cannot provide both atoms and residues!')
            sel = atoms
        elif residues is None:
            sel = selected_atoms(self.session)
        if sel is not None:
            residues = sel.unique_residues
        phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
        for dlist in (phi, psi):
            if dlist is not None:
                dlist.targets = 0
                dlist.spring_constants = 0
                self.apply_dihedral_restraints(dlist)

        for r in residues:
            for key, rlist in self._sim_distance_restraint_dict.items():
                drs = rlist[r]
                for dr in drs:
                    if dr is not None:
                        dr.spring_constant = 0
                        self.apply_distance_restraint(dr)


    def apply_distance_restraints(self, distance_restraints):
        self._sim_interface.update_distance_restraints(distance_restraints)

    def apply_distance_restraint(self, distance_restraint):
        self._sim_interface.update_distance_restraint(distance_restraint)


    def _increment_register_shift(self, *_):
        self.iw._rebuild_register_shift_nres_spinbox.stepUp()

    def _decrement_register_shift(self, *_):
        self.iw._rebuild_register_shift_nres_spinbox.stepDown()

    def _apply_register_shift(self, *_):
        from chimerax.core.atomic import selected_atoms
        from .register_shift import ProteinRegisterShifter
        nshift = self.iw._rebuild_register_shift_nres_spinbox.value()
        sel = selected_atoms(self.session)
        rs = self._register_shifter = ProteinRegisterShifter(self.session, self, sel)
        rs.shift_register(nshift)
        self.iw._rebuild_register_shift_release_button.setEnabled(False)
        self.iw._rebuild_register_shift_go_button.setEnabled(False)
        self._isolde_events.add_event_handler('register shift finish check',
                                              'completed simulation step',
                                              self._check_if_register_shift_finished)

        #~ self._register_finished_check_handler = self.triggers.add_handler(
            #~ 'completed simulation step', self._check_if_register_shift_finished)

    def _check_if_register_shift_finished(self, *_):
        if self._register_shifter.finished:
            self.iw._rebuild_register_shift_release_button.setEnabled(True)
            self._isolde_events.remove_event_handler('register shift finish check')
            #~ self.triggers.remove_handler(self._register_finished_check_handler)
            #~ self._register_finished_check_handler = None

    def _release_register_shifter(self, *_):
        if self._register_shifter is not None:
            self._register_shifter.release_all()
        self._register_shifter = None
        self.iw._rebuild_register_shift_go_button.setEnabled(True)
        self.iw._rebuild_register_shift_release_button.setEnabled(False)

    def _restrain_selected_atom_to_xyz(self, *_):
        from chimerax.core.atomic import selected_atoms
        atom = selected_atoms(self.session)[0]
        choice = self.iw._rebuild_pos_restraint_combo_box.currentIndex()
        if choice == 0:
            target = atom.coord
        else:
            target = self.session.view.center_of_rotation
        spring_constant = self.iw._rebuild_pos_restraint_spring_constant.value()
        pr = self.position_restraints[atom]
        pr.target = target
        pr.spring_constant = spring_constant
        self._sim_interface.update_position_restraint(pr)


    def release_xyz_restraints_on_selected_atoms(self, *_, sel = None):
        if sel is None:
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        restraints = self.position_restraints.in_selection(sel)
        restraints.release()
        self._sim_interface.update_position_restraints(restraints)

    def apply_xyz_restraints_to_selected_atoms(self, *_, sel = None):
        if sel is None:
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        restraints = self.position_restraints.in_selection(sel)
        self._sim_interface.update_position_restraints(restraints)
         
        
    

    def _set_rotamer_buttons_enabled(self, switch):
        self.iw._rebuild_sel_res_last_rotamer_button.setEnabled(switch)
        self.iw._rebuild_sel_res_next_rotamer_button.setEnabled(switch)
        self.iw._rebuild_sel_res_last_rotamer_button.setEnabled(switch)
        self.iw._rebuild_sel_res_rot_commit_button.setEnabled(switch)
        self.iw._rebuild_sel_res_rot_target_button.setEnabled(switch)
        self.iw._rebuild_sel_res_rot_discard_button.setEnabled(switch)



    def _next_rotamer(self, *_):
        r = self._selected_rotamer
        thisr = self._target_rotamer = r.next_rotamer(preview = True)
        self.iw._rebuild_sel_res_rot_info.setText('{} ({:.3f})'\
            .format(thisr.name, thisr.relative_abundance(self._rebuild_residue)))

    def _prev_rotamer(self, *_):
        r = self._selected_rotamer
        thisr = self._target_rotamer = r.previous_rotamer(preview = True)
        self.iw._rebuild_sel_res_rot_info.setText('{} ({:.3f})'\
            .format(thisr.name, thisr.relative_abundance(self._rebuild_residue)))

    def _commit_rotamer(self, *_):
        self._selected_rotamer.commit_current_preview()

    def _set_rotamer_target(self, *_):
        rot = self._selected_rotamer
        target = rot.target = self._target_rotamer.angles
        rot.restrained = True
        self._apply_rotamer_target_to_sim(rot)
        self._clear_rotamer()

    def _apply_rotamer_target_to_sim(self, rotamer):
        self._sim_interface.update_rotamer_target(rotamer)

    def _clear_rotamer(self, *_):
        if self._selected_rotamer is not None:
            self._selected_rotamer.cleanup()

    def release_rotamers(self, residues):
        for r in residues:
            self.release_rotamer(r)

    def release_rotamer(self, residue):
        try:
            rot = self.rotamers[residue]
        except KeyError:
            return
        rot.restrained = False
        self._apply_rotamer_target_to_sim(rot)


    def _disable_rebuild_residue_frame(self):
        if 'update_selected_residue_info' in self._event_handler.list_event_handlers():
            self._event_handler.remove_event_handler('update_selected_residue_info')
        self.iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        self.iw._rebuild_sel_res_pep_info.setText('')
        self.iw._rebuild_sel_res_rot_info.setText('')
        self.iw._rebuild_sel_res_rot_target_button.setText('Set target')
        self.iw._rebuild_sel_residue_frame.setDisabled(True)

    def _get_rotamer_list_for_selected_residue(self, res):
        pass


    def _update_selected_residue_info_live(self, *_):
        from math import degrees
        res = self._rebuild_residue
        c = self._sel_res_update_counter
        s = self._steps_per_sel_res_update
        c = (c + 1) % s
        if c == 0:
            # Get the residue's omega value
            if self._rebuild_res_update_omega:
                omega = self._rebuild_res_omega
                oval = degrees(omega.value)
                if oval <= -150 or oval >= 150:
                    pep_type = 'trans'
                elif oval >= -30 and oval <= 30:
                    pep_type = 'cis'
                else:
                    pep_type = 'twisted'
                self.iw._rebuild_sel_res_pep_info.setText(pep_type)
        self._sel_res_update_counter = c

    def _flip_peptide_bond(self, *_):
        res = self._rebuild_residue
        self.flip_peptide_bond(res)

    def flip_peptide_bond(self, res):
        '''
        A bit tricky. This involves flipping phi for this residue and
        psi for the preceding residue. Ideally, we don't want to leave
        them restrained once the flip is complete.
        '''

        bd = self._mobile_backbone_dihedrals
        phi = bd.phi.by_residue(res)
        prev_c = phi.atoms.filter(phi.atoms.names == 'C')[0]
        prev_r = prev_c.residue
        psi = bd.psi.by_residue(prev_r)

        targets = []
        for d in (phi, psi):
            v = d.value
            if v < 0:
                d.target = v+pi
            else:
                d.target = v-pi
            d.spring_constant = defaults.PHI_PSI_SPRING_CONSTANT

        self.apply_dihedral_restraint(phi)
        self.apply_dihedral_restraint(psi)

        self._pep_flip_timeout_counter = 0
        self._pep_flip_dihedrals = (phi, psi)
        self.iw._rebuild_sel_res_pep_flip_button.setEnabled(False)
        self._isolde_events.add_event_handler('pep flip timeout',
                                              'completed simulation step',
                                              self._check_pep_flip)

    def apply_dihedral_restraint(self, dihedral):
        '''
        Restrain a dihedral to a desired angle with a spring constant.
        @param dihedral_type:
            One of 'phi', 'psi' or 'omega'. Sidechain torsions are handled
            collectively as rotamers
        @param dihedral:
            The dihedral object. The target and spring constant will be
            taken from its properties
        '''
        self._sim_interface.update_dihedral_restraint(dihedral)

    def apply_dihedral_restraints(self, dihedrals):
        '''
        Apply restraints for a set of dihedrals at once. All dihedrals
        must be of the same type.
        '''
        all_names = dihedrals.names
        unique_names = numpy.unique(all_names)
        if len(unique_names) == 1:
            self._sim_interface.update_dihedral_restraints(dihedrals)
        else:
            for name in unique_names:
                indices = numpy.where(all_names == name)[0]
                self._sim_interface.update_dihedral_restraints(dihedrals[indices])


    def release_dihedral_restraint(self, dihedral):
        dihedral.target = 0
        dihedral.spring_constant = 0
        self.apply_dihedral_restraint(dihedral)

    def _check_pep_flip(self, *_):
        done = False
        if self._pep_flip_timeout_counter >= 20:
            print('Unable to flip peptide. Giving up.')
            done = True
        dihedrals = self._pep_flip_dihedrals
        if not done:
            done = True
            for d in dihedrals:
                if abs(d.value - d.target) > pi/6:
                    done = False
        if not done:
            self._pep_flip_timeout_counter += 1
            return
        else:
            self.release_dihedral_restraint(dihedrals[0])
            self.release_dihedral_restraint(dihedrals[1])
            self._isolde_events.remove_event_handler('pep flip timeout')
            self.iw._rebuild_sel_res_pep_flip_button.setEnabled(True)



    def _flip_cis_trans(self, *_):
        res = self._rebuild_residue
        self.flip_peptide_omega(res)

    def flip_peptide_omega(self, res):
        bd = self._mobile_backbone_dihedrals
        omega = bd.omega.by_residue(res)
        if omega is None:
            import warnings
            warnings.warn('Could not find residue or residue has no omega dihedral.')
            return
        oval = omega.value
        if abs(oval) <= defaults.CIS_PEPTIDE_BOND_CUTOFF:
            # cis, flip to trans
            target = pi
        else:
            target = 0
        omega.target = target
        self.apply_dihedral_restraint(omega)



    ####
    # Validation tab
    ####
    def _show_rama_plot(self, *_):
        self.iw._validate_rama_stub_frame.hide()
        self.iw._validate_rama_main_frame.show()
        if self._rama_plot is None:
            # Create the basic MatPlotLib canvas for the Ramachandran plot
            self._prepare_ramachandran_plot()
        if self._simulation_running and self.params.track_ramachandran_status:
            self._rama_go_live()
        else:
            self._rama_static_plot()

    def _hide_rama_plot(self, *_):
        self.iw._validate_rama_main_frame.hide()
        self.iw._validate_rama_stub_frame.show()
        self._rama_go_static()


    def _change_rama_case(self, *_):
        case_key = self.iw._validate_rama_case_combo_box.currentData()
        self._rama_plot.change_case(case_key)

    def _rama_go_live(self, *_):
        self._update_rama_plot = True
        self.iw._validate_rama_sel_combo_box.setDisabled(True)
        self.iw._validate_rama_go_button.setDisabled(True)

    def _rama_go_static(self, *_):
        self._update_rama_plot = False
        self.iw._validate_rama_sel_combo_box.setEnabled(True)
        self.iw._validate_rama_go_button.setEnabled(True)

    def _rama_static_plot(self, *_):
        model = self._selected_model
        whole_model = bool(self.iw._validate_rama_sel_combo_box.currentIndex())
        if whole_model:
            sel = model.atoms
            bd = self.backbone_dihedrals
            self._rama_plot.update_scatter(bd, force_update = True)

        else:
            sel = model.atoms.filter(model.atoms.selected)
            residues = sel.unique_residues
            if len(residues):
                phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
                from . import dihedrals
                bd = dihedrals.Backbone_Dihedrals(self.session, phi=phi, psi=psi, omega=omega)
                self._rama_plot.update_scatter(bd, force_update = True)
            else:
                self._rama_plot.update_scatter(force_update = True)

    def _show_peptide_validation_frame(self, *_):
        self.iw._validate_pep_stub_frame.hide()
        self.iw._validate_pep_main_frame.show()

    def _hide_peptide_validation_frame(self, *_):
        self.iw._validate_pep_main_frame.hide()
        self.iw._validate_pep_stub_frame.show()

    def _update_iffy_peptide_lists(self, *_):
        ov = self.omega_validator
        model = self._selected_model
        clist = self.iw._validate_pep_cis_list
        tlist = self.iw._validate_pep_twisted_list
        clist.clear()
        tlist.clear()
        if model != ov.current_model:
            sel = model.atoms
            from . import dihedrals
            bd = self.backbone_dihedrals
            ov.load_structure(model, bd.omega)
        cis, twisted = ov.find_outliers()
        ov.draw_outliers(cis, twisted)
        from PyQt5.QtWidgets import QListWidgetItem
        from PyQt5.Qt import QColor, QBrush
        from PyQt5.QtCore import Qt
        badColor = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        for c in cis:
            pre, r = c.residues.unique()
            label = r.chain_id + ' ' \
                    + str(pre.number) + ' - ' + str(r.number) + '\t' \
                    + pre.name + ' - ' + r.name
            list_item = QListWidgetItem(label)
            list_item.data = r
            if r.name != 'PRO':
                list_item.setBackground(badColor)
            clist.addItem(list_item)
        for t in twisted:
            pre, r = t.residues.unique()
            label = r.chain_id + ' ' \
                    + str(pre.number) + ' - ' + str(r.number) + '\t' \
                    + pre.name + ' - ' + r.name
            list_item = QListWidgetItem(label)
            list_item.data = r
            list_item.setBackground(badColor)
            tlist.addItem(list_item)

    def _show_selected_iffy_peptide(self, item):
        res = item.data
        from . import view
        view.focus_on_selection(self.session, self.session.main_view, res.atoms)
        self.session.selection.clear()
        res.atoms.selected = True


    ##############################################################
    # Simulation global settings functions
    ##############################################################
    def _change_sim_mode(self, *_):
        sm = self.iw._sim_basic_mode_combo_box.currentData()
        self.sim_mode = sm
        self.gui._change_experience_level_or_sim_mode()
        self._update_model_list()

    def _change_force_field(self):
        ffindex = self.iw._sim_force_field_combo_box.currentIndex()
        ffs = self._available_ffs
        self._sim_main_ff = ffs.main_files[ffindex]
        self._sim_implicit_solvent_ff = ffs.implicit_solvent_files[ffindex]

    def _change_water_model(self):
        ffindex = self.iw._sim_water_model_combo_box.currentIndex()
        self._sim_water_ff = self._available_ffs.explicit_water_files[ffindex]

    def _change_selected_model(self, *_, model = None, force = False):
        if len(self._available_models) == 0:
            return
        if self._simulation_running:
            return
        iw = self.iw
        if model is not None:
            # Find and select the model in the master combo box, which
            # will automatically call this function again with model = None
            index = iw._master_model_combo_box.findData(model)
            iw._master_model_combo_box.setCurrentIndex(index)
            return
        m = iw._master_model_combo_box.currentData()
        if force or (self._selected_model != m and m is not None):
            self._status('Analysing model and preparing restraints. Please be patient.')
            from . import util
            util.add_disulfides_from_model_metadata(m)
            self._selected_model = m
            self.session.selection.clear()
            self._selected_model.selected = True
            self._prepare_all_interactive_restraints(m)
            self._update_chain_list()
            if isinstance(m.parent, clipper.CrystalStructure):
                self._initialize_xtal_maps(m.parent)
            self.triggers.activate_trigger('selected model changed', data=m)
        self._status('')

    def _prepare_all_interactive_restraints(self, model):
        '''
        Any restraints you may want to turn on and off interactively must
        be defined at the beginning of each simulation. Otherwise, adding
        a new restraint would require a costly reinitialisation of the
        simulation context. For the most commonly-used restraints it's
        therefore best to pre-define the necessary objects for the entire
        model at the point the model is selected, and pull out the
        necessary pieces from these each time a simulation is started.
        '''
        m = model
        from . import dihedrals, rotamers
        from . import backbone_restraints as br
        from . import position_restraints as pr
        # Torsion restraints
        self.backbone_dihedrals = dihedrals.Backbone_Dihedrals(self.session, m)
        self.rotamers = rotamers.all_rotamers_in_selection(self.session, m.atoms)
        # Distance restraints
        self.ca_to_ca_plus_two = br.CA_to_CA_plus_Two(self.session, m)
        self.o_to_n_plus_four = br.O_to_N_plus_Four(self.session, m)
        # Positional restraints (one per heavy atom)
        self.position_restraints = pr.Atom_Position_Restraints(
                self.session, self.selected_model,
                m.atoms, triggers = self.triggers, create_target=True)


    def _select_whole_model(self, *_):
        if self._selected_model:
            self._selected_model.selected = True
            self._selected_atoms = self._selected_model.atoms


    def _change_selected_chains(self,*_):
        if len(self._available_models) == 0:
            return
        if self._simulation_running:
            return
        lb = self.iw._sim_basic_mobile_chains_list_box
        m = self._selected_model
        all_chains = m.chains
        all_chain_ids = list(all_chains.chain_ids)
        lb_sels = lb.selectedItems()
        self.session.selection.clear()
        for s in lb_sels:
            chain_text = s.text()
            chain_index = all_chain_ids.index(chain_text)
            all_chains[chain_index].existing_residues.atoms.selected = True
        from chimerax.core.atomic import selected_atoms
        self._selected_atoms = selected_atoms(self.session)

    def _change_b_and_a_padding(self, *_):
        self.params.num_selection_padding_residues = self.iw._sim_basic_mobile_b_and_a_spinbox.value()

    def _change_soft_shell_cutoff(self, *_):
        self.params.soft_shell_cutoff_distance = self.iw._sim_basic_mobile_sel_within_spinbox.value()

    def _change_soft_shell_fix_backbone(self, *_):
        self.params.fix_soft_shell_backbone = not self.iw._sim_basic_mobile_sel_backbone_checkbox.checkState()

    def _change_sim_platform(self, *_):
        self.sim_platform = self.iw._sim_platform_combo_box.currentText()

    def _add_or_change_em_map_from_gui(self, *_):
        iw = self.iw
        name = iw._em_map_name_field.text()
        m_id = iw._em_map_model_combo_box.currentText()
        model = self._available_volumes[m_id]
        cutoff = iw._em_map_cutoff_spin_box.value()
        coupling_constant = iw._em_map_coupling_spin_box.value()
        style = iw._em_map_style_combo_box.currentText()
        color = iw._em_map_color_combo_box.currentText()
        contour = iw._em_map_contour_spin_box.value()
        contour_units = iw._em_map_contour_units_combo_box.currentText()
        mask = iw._em_map_masked_checkbox.checkState()
        if self._add_new_map:
            self.add_map(name, model, cutoff, coupling_constant, style, color, contour, contour_units, mask)
        else:
            m = self.master_map_list[name]
            m.change_map_parameters(model, cutoff, coupling_constant, style, color, contour, contour_units, mask)

    def _remove_em_map_from_gui(self, *_):
        name = self.iw._em_map_name_field.text()
        self.remove_map(name)

    def _change_display_of_selected_map(self, *_):
        iw = self.iw
        m_id = iw._em_map_model_combo_box.currentText()
        if m_id == "":
            return
        model = self._available_volumes[m_id]
        styleargs = {}
        style = iw._em_map_style_combo_box.currentData()
        if style is not None:
            styleargs = self._map_style_settings[style]
        map_color = iw._em_map_color_combo_box.currentData()
        if map_color is not None:
            styleargs['color'] = [map_color]
        if len(styleargs) != 0:
            from chimerax.core.map import volumecommand
            volumecommand.volume(self.session, [model], **styleargs)

    def _change_contour_level_of_selected_map(self, *_):
        iw = self.iw
        m_id = iw._em_map_model_combo_box.currentText()
        model = self._available_volumes[m_id]
        contour_val = iw._em_map_contour_spin_box.value()
        contour_units = iw._em_map_contour_units_combo_box.currentText()
        if contour_units == 'sigma':
            map_sigma = model.mean_sd_rms()[1]
            contour_val = contour_val * map_sigma
        from chimerax.core.map import volumecommand
        volumecommand.volume(self.session,[model], level=[[contour_val]], cap_faces = False)


    def _change_contour_units_of_selected_map(self, *_):
        iw = self.iw
        sb = iw._em_map_contour_spin_box
        cb = iw._em_map_contour_units_combo_box
        m_id = iw._em_map_model_combo_box.currentText()
        model = self._available_volumes[m_id]
        contour_units = cb.currentText()
        contour_val = sb.value()
        map_sigma = model.mean_sd_rms()[1]
        if contour_units == 'sigma':
            contour_val = contour_val / map_sigma
        else:
            contour_val = contour_val * map_sigma
        sb.setValue(contour_val)

    def _set_hide_surroundings(self,*_):
        self.params.hide_surroundings_during_sim = self.iw._sim_hide_surroundings_toggle.checkState()



    def add_map(self, name, vol, cutoff, coupling_constant,
        is_difference_map = False, style = None, color = None,
        contour = None, contour_units = None, mask = True, crop = True):
        if name in self.master_map_list:
            for key in self.master_map_list:
                print(key)
            raise Exception('Each map must have a unique name!')
        # Check if this model is a unique volumetric map
        if len(vol.models()) !=1 or not hasattr(vol, 'grid_data'):
            raise Exception('vol must be a single volumetric map object')

        from .volumetric import IsoldeMap
        new_map = IsoldeMap(self.session, name, vol, cutoff, coupling_constant,
        is_difference_map = is_difference_map, style = style,
        color = color, contour = contour, contour_units = contour_units,
        mask = mask, crop = crop)
        self.master_map_list[name] = new_map
        self._update_master_map_list_combo_box()
        return new_map


    def remove_map(self, name):
        result = self.master_map_list.pop(name, 'Not present')
        if result == 'Not present':
            print(name + ' is an unrecognised key.')
        self._update_master_map_list_combo_box()

    ##############################################################
    # Interactive restraints functions
    ##############################################################

    def _change_peptide_bond_restraints(self, *_):
        '''
        Menu-driven function to turn peptide bond restraints on/off
        for some subset of residues.
        '''
        enable = bool(self.iw._restraints_pep_on_off_combo_box.currentIndex())
        selected = bool(self.iw._restraints_pep_all_sel_combo_box.currentIndex())

        bd = self._mobile_backbone_dihedrals

        if selected:
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        else:
            sel = self._total_mobile

        res = sel.unique_residues

        if enable:
            k = self.peptide_bond_restraints_k
        else:
            k = 0

        omegas = bd.omega.by_residues(res)
        #omegas = []

        #for r in res:
            #rindex = bd.residues.index(r)
            #if rindex != -1:
                #omegas.append(bd.omega[rindex])

        if enable:
            self.apply_peptide_bond_restraints(omegas)
        else:
            self.remove_peptide_bond_restraints(omegas)






    ##############################################################
    # Simulation prep
    ##############################################################

    def reset_sim_params(self):
        '''
        Reset all the simulation parameters back to their defaults.
        '''
        from .openmm.sim_interface import SimParams
        self.sim_params = SimParams()


    def start_sim(self):
        try:
            #~ self._start_sim()
            self._start_threaded_sim()
        except Exception as e:
            self.triggers.activate_trigger('simulation terminated',
                                        (sim_outcomes.FAIL_TO_START, e))
            raise


    def _start_threaded_sim(self):
        session = self.session
        log = session.logger.info
        sel_model = self.selected_model

        if self._simulation_running:
            print('You already have a simulation running!')
            return
        if self._logging:
            self._log('Initialising simulation')
        log('Simulation should start now')

        sm = self._sim_modes
        if self.sim_mode in [sm.xtal, sm.em]:
            if not len(self.master_map_list):
                self._no_maps_warning()
                return

        self._status('Defining simulation selection...')

        # Define final mobile selection
        main_sel = self._get_main_sim_selection()

        # Define "soft shell" of mobile atoms surrounding main selection
        soft_shell = self._soft_shell_atoms = self.get_shell_of_residues(
            main_sel, self.params.soft_shell_cutoff_distance)


        # Define fixed selection (whole residues with atoms coming within
        # a cutoff of the mobile selection
        total_mobile = self._total_mobile = main_sel.merge(soft_shell)
        hard_shell = self._hard_shell_atoms = self.get_shell_of_residues(
            total_mobile, self.params.hard_shell_cutoff_distance)

        sc = self._total_sim_construct = total_mobile.merge(hard_shell)

        # We need to make sure that the atoms in our important selections are
        # sorted in the same order as in the master model to avoid Bad Things(TM)
        # happening.
        all_a = sel_model.atoms
        tm_i = all_a.indices(total_mobile)
        sc_i = all_a.indices(sc)
        tm_i.sort()
        sc_i.sort()
        total_mobile = self._total_mobile = all_a[tm_i]
        sc = self._total_sim_construct = all_a[sc_i]

        self._total_mobile_indices = sc.indices(total_mobile)

        sb = self._total_sim_bonds = sc.intra_bonds
        surr = self._surroundings = all_a.subtract(sc)

        if self.sim_mode == sm.xtal:
            sel_model.parent.isolate_and_cover_selection(
                total_mobile, include_surrounding_residues = 0,
                show_context = self.params.hard_shell_cutoff_distance,
                mask_radius = 4, extra_padding = 10,
                hide_surrounds = self.params.hide_surroundings_during_sim, focus = False)

        else:
            surr.hides |= control.HIDE_ISOLDE
            self._surroundings_hidden = True

        # Cache all the colors so we can revert at the end of the simulation
        self._original_atom_colors = sc.colors
        self._original_atom_draw_modes = sc.draw_modes
        self._original_bond_radii = sb.radii
        self._original_atom_radii = sc.radii
        self._original_display_state = self._selected_model.atoms.displays

        hsb = hard_shell.intra_bonds
        hsb.radii = 0.1
        hard_shell.radii = 0.1
        hard_shell.draw_modes = hard_shell.BALL_STYLE

        fixed_flags = numpy.zeros(len(sc), numpy.bool)
        fixed_flags[sc.indices(hard_shell)] = True

        if self.params.fix_soft_shell_backbone:
            indices = sc.indices(soft_shell)
            names = soft_shell.names
            for i, n in zip(indices, names):
                if n in ['N','C','O','H','H1','H2','H3']:
                    fixed[i] = True


        # Collect backbone dihedral atoms and prepare the Ramachandran
        # validator

        self._status('Organising dihedrals...')

        from . import dihedrals
        all_bd = self.backbone_dihedrals
        sim_phi, sim_psi, sim_omega = all_bd.by_residues(total_mobile.unique_residues)
        bd = self._mobile_backbone_dihedrals \
            = dihedrals.Backbone_Dihedrals(self.session, phi = sim_phi, psi = sim_psi, omega = sim_omega)

        if self.params.track_ramachandran_status:
            bd.CAs.draw_modes = 1

            self.omega_validator.load_structure(sel_model, bd.omega)


        distance_restraints = self._sim_distance_restraint_dict = {
            'ca_to_ca_plus_two':    self.ca_to_ca_plus_two.in_selection(sc),
            'o_to_n_plus_four':     self.o_to_n_plus_four.in_selection(sc),
            }

        position_restraints = self._sim_pos_restr =\
            self.position_restraints.in_selection(total_mobile)


        tuggable_atoms = total_mobile[total_mobile.element_names != 'H']

        from .openmm.sim_interface import ChimeraXSimInterface
        sp = self.sim_params
        si = self._sim_interface = ChimeraXSimInterface(self.session, self)
        si.start_sim_thread(sp, sc, tuggable_atoms, fixed_flags, bd, self.rotamers,
                            distance_restraints, position_restraints, self.master_map_list)





        self._status('Simulation running')


    def _sim_start_cb(self, *_):
        '''
        Register all event handlers etc. that have to be running during the
        simulation.
        '''
        self._simulation_running = True
        self._update_sim_control_button_states()


        if self.params.track_ramachandran_status and self._rama_plot is not None:
            self._rama_go_live()

        if self.params.remask_maps_during_sim:
            self._isolde_events.add_event_handler('rezone maps during sim',
                                                  'completed simulation step',
                                                  self._rezone_maps)

    def _sim_end_cb(self, name, outcome):
        pass
        if outcome[0] == sim_outcomes.UNSTABLE:
            print('''Unable to minimise excessive force. Giving up. Try starting
                a simulation on a larger selection to give the surroundings
                more room to settle.''')
            '''
            TODO: Pop up a dialog box, giving the user the opportunity to save
            or discard the existing coordinates, and/or start a simulation
            with a wider selection to give the minimiser more room to work.
            '''
        elif outcome[0] == sim_outcomes.FAIL_TO_START:
            '''
            TODO: Pop up a dialog box giving the details of the exception, then
            clean up.
            '''
            print('Failed to start')
            pass
        si = self._sim_interface
        # Otherwise just clean up
        self._simulation_running = False
        self._event_handler.remove_event_handler(
            'do_sim_steps_on_gui_update', error_on_missing = False)
        self._isolde_events.remove_event_handler(
            'rezone maps during sim', error_on_missing = False)

            #~ self.triggers.remove_handler(self._map_rezone_handler)
            #~ self._map_rezone_handler = None

        self._disable_rebuild_residue_frame()

        self.sim = None
        self._system = None
        self._update_menu_after_sim()
        #mouse_mode_names = self._mouse_modes.get_names()
        #for n in mouse_mode_names:
            #self._mouse_modes.remove_mode(n)
        self._mouse_tugger = None
        for d in self._haptic_devices:
            d.cleanup()
        self._selected_model.atoms.displays = self._original_display_state
        self._surroundings.hides &= ~control.HIDE_ISOLDE
        self._surroundings_hidden = False
        self._total_sim_construct.colors = self._original_atom_colors
        self._total_sim_construct.draw_modes = self._original_atom_draw_modes
        self._total_sim_construct.radii = self._original_atom_radii
        self._total_sim_bonds.radii = self._original_bond_radii
        self._total_sim_construct = None
        self._surroundings = None
        if self.params.track_ramachandran_status:
            # Update one last time
            if self._update_rama_plot:
                self._rama_plot.update_scatter(self._mobile_backbone_dihedrals)
                self.update_omega_check()
        self._rama_go_static()
        self.iw._rebuild_sel_residue_frame.setDisabled(True)
        self.omega_validator.current_model = None
        self._last_max_force = inf
        self._unstable_min_rounds = 0
        self._sim_is_unstable = False
        if self.sim_mode == self._sim_modes['xtal']:
            cs = self._selected_model.parent
            cs.xmaps.live_scrolling = True
            cs.live_atomic_symmetry = True
        self._status('')


    def _get_main_sim_selection(self):
        '''
        Get the primary selection that will be the focus of the simulation. In
        most cases, this involves taking the existing selected atoms, expanding
        the selection to cover whole residues, and further expanding the
        selection by the desired padding up and down the chain (stopping at
        chain breaks and ends).
        '''

        mode = self._sim_selection_mode
        modes = self._sim_selection_modes

        if mode == modes.chain or mode == modes.whole_model:
            # Then everything is easy. The selection is already defined
            pass
        elif mode == modes.from_picked_atoms:
            # A bit more complex. Have to work through the model to find
            # the picked atoms (making sure only one model is selected!),
            # then work back and forward from each picked atom to expand
            # the selection by the specified number of residues.
            pad = self.params.num_selection_padding_residues
            from chimerax.core.atomic import selected_atoms
            selatoms = selected_atoms(self.session)
            us = selatoms.unique_structures
            if len(us) != 1:
                print(len(us))
                for m in us:
                    print(m.category)
                raise Exception('Selected atoms must all be in the same model!')
            sm = self._selected_model = us[0]
            # selatoms_by_chain = selatoms.by_chain
            # selchains = [row[1] for row in selatoms_by_chain]
            all_atoms = sm.atoms
            all_res = sm.residues
            sel_res = selatoms.unique_residues
            sel_res_indices = all_res.indices(sel_res)
            from chimerax.core.atomic.structure import Structure

            allfrags = sm.polymers(missing_structure_treatment=Structure.PMS_NEVER_CONNECTS)

            sel_frags = []
            sel_frag_res_indices = []
            for frag in allfrags:
                if not frag.atoms.num_selected:
                    continue
                sel_frags.append(frag)
                sel_frags.append(all_res.indices(frag))


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

        elif mode == modes.custom:
            # relatively simple. Just need to apply the in-built selection
            # text parser. To be completed.
            pass


    # Get a shell of whole residues within a user-defined cut-off surrounding
    # an existing set of selected atoms. Expects the existing_sel set to be
    # whole residues, and all within the same model.

    def get_shell_of_residues(self, existing_sel, dist_cutoff):
        from chimerax.core.geometry import find_close_points
        from chimerax.core.atomic import selected_atoms, Atoms, concatenate
        selatoms = existing_sel
        us = selatoms.unique_structures
        if len(us) !=1:
            raise Exception('selection should contain atoms from a single molecule!')
        allatoms = us[0].atoms
        unselected_atoms = allatoms.subtract(selatoms)
        selcoords = selatoms.coords
        unselcoords = unselected_atoms.coords
        ignore, shell_indices = find_close_points(selcoords, unselcoords, dist_cutoff)
        shell_atoms = unselected_atoms[shell_indices].unique_residues.atoms
        return shell_atoms







    ##############################################################
    # Simulation on-the-fly control functions
    ##############################################################



    def pause_sim_toggle(self):
        print('This function should toggle pause/resume of the sim')
        if self._simulation_running:
            self._sim_interface.toggle_pause()

    def _sim_pause_cb(self, *_):
        self._sim_paused = True
        self._status('Simulation paused')
        self.iw._sim_pause_button.setText('Resume')

    def _sim_resume_cb(self, *_):
        self._sim_paused = False
        self._status('Simulation running')
        self.iw._sim_pause_button.setText('Pause')


    def discard_sim(self):
        print("""This function should stop the simulation and revert to
                 the original coordinates""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.stop_sim(sim_outcomes.DISCARD, None)

    def commit_sim(self):
        print("""This function should stop the simulation and write the
                 coordinates to the target molecule""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.stop_sim(sim_outcomes.COMMIT, None)

    def minimize(self):
        print('Minimisation mode')
        self.simulation_type = 'min'
        self._update_sim_control_button_states()

    def equilibrate(self):
        print('Equilibration mode')
        self.simulation_type = 'equil'
        self._update_sim_control_button_states()


    ####
    # Restraint controls
    ####
    def apply_peptide_bond_restraints(self, dihedrals, target = None, units = 'deg'):
        '''
        Apply restraints to a list of omega dihedrals. If target is None,
        then each dihedral will be automatically restrained to trans or
        cis according to its starting configuration. Otherwise, the target
        may be a value (interpreted as degrees if units = 'deg' or radians
        if units = 'rad'), or either 'trans' or 'cis'.
        '''
        from numbers import Number
        if units not in ('deg','rad'):
            raise Exception('Units should be either "deg" or "rad"')

        cr = self.cis_peptide_bond_range
        sh = self._sim_handler
        sc = self._total_sim_construct
        k = self.peptide_bond_restraints_k
        sim = self.sim
        context = None
        if sim is not None and hasattr(sim, 'context'):
            context = sim.context

        from math import pi, radians
        if target == 'trans':
            t = pi
        elif target == 'cis':
            t = 0.0
        elif isinstance(target, Number):
            if units == 'deg':
                t = radians(target)
            else:
                t = target
        elif target is not None:
            raise Exception('Target must be either a number, "trans", "cis", or None')
        else:
            dvals = dihedrals.values
            targets = (numpy.logical_or(dvals > cr[1], dvals < cr[0])).astype(float) * pi

        # Get all atom indices in one go because the lookup is expensive
        #indices = numpy.reshape(sc.indices(dihedrals.atoms),[len(dihedrals),4])


        for i, d in enumerate(dihedrals):
            if target is None:
                t = targets[i]
            #sh.set_dihedral_restraint(d.sim_index, indices[i], t, k)
            sh.update_dihedral_restraint(d.sim_index, target=t, k=k)

    def remove_peptide_bond_restraints(self, dihedrals):
        '''
        Remove restraints on a list of peptide bonds (actually, just set
        their spring constants to zero). Simulation must already be running.
        '''
        sc = self._total_sim_construct
        sh = self._sim_handler
        sh.set_dihedral_restraints(dihedrals, 0, 0)


    #############################################
    # Main simulation functions to be run once per GUI update
    #############################################

    def do_sim_steps(self,*_):
        if self._logging:
            self._log('Running ' + str(self.sim_steps_per_update) + ' steps')

        v = self.session.main_view
#        if v.frame_number == self._last_frame_number:
#            return # Make sure we draw a frame before doing another MD calculation
        sh = self._sim_handler

        s = self.sim
        c = s.context
        integrator = c.getIntegrator()
        steps = self.sim_steps_per_update
        minsteps = self.min_steps_per_update
        mode = self.simulation_type
        startup = self._sim_startup
        self._sim_startup_counter += 1
        s_max_count = self._sim_startup_rounds
        if self._sim_startup and self._sim_startup_counter >= s_max_count:
            print('Maximum number of startup minimisation rounds reached. \
                    starting dynamics anyway...')
            startup = self._sim_startup = False
        pos = self._particle_positions
        sc = self._total_sim_construct
        surr = self._surroundings

        # Check if we need to hide or show the surroundings
        if self.params.hide_surroundings_during_sim:
            if not self._surroundings_hidden:
                surr.displays = False
                self._surroundings_hidden = True
        elif self._surroundings_hidden:
            surr.displays = True
            self._surroundings_hidden = False

        sh.update_restraints_in_context(c)
        if startup or mode == 'min':
            newpos, max_force, max_index = self._get_positions_and_max_force()
        else:
            newpos, fast_indices = self._get_and_check_positions()
        sc.coords = self._particle_positions = newpos
        if startup:
            print('Startup round {} max force: {:0.0f} kJ/mol/nm'
                    .format(self._sim_startup_counter, max_force))
            if abs(max_force - self._last_max_force) < 1.0 and max_force < self._max_allowable_force:
                print('Minimisation converged. Starting dynamics.')
                startup = self._sim_startup = False
            self._last_max_force = max_force
        elif mode == 'equil' and fast_indices is not None:
            self._sim_is_unstable = True
            if mode == 'equil':
                # revert to the coordinates before this simulation step
                #c.setPositions(pos/10)
                self.simulation_type = 'min'
                return
        if self._sim_is_unstable:
            if self._unstable_min_rounds >= self.max_unstable_rounds:
                # We have a problem we can't fix. Terminate the simulation
                self.triggers.activate_trigger('simulation terminated', (sim_outcomes.UNSTABLE, None))
                return
            bad_atom = self._total_mobile[max_index]
            bad_res = bad_atom.residue
            outstr = '''
            Simulation is unstable! Atom {} from residue {} of chain {}
            is experiencing a net force of {:0.0f} kJ/mol/nm.
            '''.format(bad_atom.name, bad_res.number, bad_atom.chain_id,
                        max_force)
            print(outstr)
            self._unstable_min_rounds += 1
            if max_force < self._max_allowable_force:
                # We're back to stability. We can go back to equilibrating
                self._sim_is_unstable = False
                self.simulation_type = 'equil'
                self._unstable_min_rounds = 0



        # Mouse interaction
        mtug = self._mouse_tugger
        t_force = self._tugging_force
        t_k = self.tug_force_constant
        cur_tug = self._currently_tugging
        tugging, tug_atom, xyz0 = mtug.status
        # If tugging is true, we need to update the parameters for the picked atom
        if tugging:
            xyz0 = xyz0 / 10 # OpenMM coordinates are in nanometres
            params = [t_k, *xyz0]
            tug_index = self._total_sim_construct.index(tug_atom)
            if not cur_tug:
                self._last_tugged_index = tug_index
                self._currently_tugging = True
            sh.set_custom_external_force_particle_params('tug', tug_index, params)
            sh.update_force_in_context('tug', c)
        # If cur_tug is true and tugging is false, we need to disconnect the atom from the tugging force
        elif cur_tug:
            sh.set_custom_external_force_particle_params('tug', self._last_tugged_index, [0,0,0,0])
            sh.update_force_in_context('tug', c)
            self._currently_tugging = False
            self._last_tugged_index = None

        # If both cur_tug and tugging are false, do nothing



        # Haptic interaction
        if self._use_haptics:
            hh = self.session.HapticHandler
            from . import picking
            for d in self._haptic_devices:
                i = d.index
                pointer_pos = hh.getPosition(i, scene_coords = True)
                if self._haptic_highlight_nearest_atom[i] and not d.tugging:
                    a = self._haptic_current_nearest_atom[i]
                    if a is not None and not a.deleted:
                        # set the previously picked atom back to standard
                        a.draw_mode = self._haptic_current_nearest_atom_draw_mode[i]
                        a.color = self._haptic_current_nearest_atom_color[i]
                    a = self._haptic_current_nearest_atom[i] = picking.pick_closest_to_point(
                            self.session, pointer_pos, self._total_mobile,3.0)
                    if a is not None:
                        self._haptic_current_nearest_atom_draw_mode[i] = a.draw_mode
                        a.draw_mode = 1
                        self._haptic_current_nearest_atom_color[i] = a.color
                        a.color = [0, 255, 255, 255] # bright cyan

                b = hh.getButtonStates(i)
                # Button 0 controls tugging
                if b[0]:
                    if not d.tugging:
                        a = self._haptic_tug_atom[i] = picking.pick_closest_to_point(
                            self.session, pointer_pos, self._total_mobile,
                            3.0, displayed_only = True,
                            tug_hydrogens = self.tug_hydrogens)
                        if a is not None:
                            tug_index = self._haptic_tug_index[i] = self._total_sim_construct.index(a)
                            # Start the haptic device feedback loop
                            d.start_tugging(a)
                            params = [self._haptic_force_constant, *(hh.getPosition(i, scene_coords = True)/10)]
                            sh.set_custom_external_force_particle_params('tug', tug_index, params)
                            sh.update_force_in_context('tug', c)

                    else:
                        d.update_target()
                        # Update the target for the tugged atom in the simulation
                        tug_index = self._haptic_tug_index[i]
                        params = [self._haptic_force_constant, *(hh.getPosition(i, scene_coords = True)/10)]
                        sh.set_custom_external_force_particle_params('tug', tug_index, params)
                        sh.update_force_in_context('tug', c)
                else:
                    if d.tugging:
                        d.stop_tugging()
                        tug_index = self._haptic_tug_index[i]
                        sh.set_custom_external_force_particle_params('tug', tug_index, [0,0,0,0])
                        self._haptic_tug_atom[i] = None




        if self._temperature_changed:
            integrator.setTemperature(self.simulation_temperature)
            c.setVelocitiesToTemperature(self.simulation_temperature)
            self._temperature_changed = False


        if startup:
            #s.minimizeEnergy(maxIterations = steps)
            s.minimizeEnergy(maxIterations = 1000)
            self._sim_startup_counter += 1
        elif self._sim_is_unstable:
            # Run a few timesteps to jiggle out of local minima
            s.step(5)
            s.minimizeEnergy()
            c.setVelocitiesToTemperature(self.simulation_temperature)
        elif mode == 'min':
            s.minimizeEnergy(maxIterations = minsteps)
            c.setVelocitiesToTemperature(self.simulation_temperature)
        elif mode == 'equil':
            s.step(steps)
        else:
            raise Exception('Unrecognised simulation mode!')


        from simtk import unit
        self._last_frame_number = v.frame_number
        if self._logging:
            self._log('Ran ' + str(self.sim_steps_per_update) + ' steps')
        if self.params.track_ramachandran_status:
            self.update_ramachandran()
            self.update_omega_check()

        self.triggers.activate_trigger('completed simulation step', data=None)


    def _rezone_maps(self, *_):
        self._map_rezone_counter += 1
        self._map_rezone_counter %= self.params.rounds_per_map_remask
        if self._map_rezone_counter == 0:
            self.rezone_maps()

    def rezone_maps(self):
        from chimerax.core.commands.sop import surface_zone
        for key, m in self.master_map_list.items():
            v = m.get_source_map()
            cutoff = m.get_mask_cutoff()
            surface_zone(self.session, v.surface_drawings,
                near_atoms = self._total_mobile, range = cutoff)


    def update_ramachandran(self, *_):
        self._rama_counter = (self._rama_counter + 1) % self.params.rounds_per_rama_update
        if self._rama_counter == 0:
            rv = self.rama_validator
            bd = self._mobile_backbone_dihedrals
            rv.get_scores(bd, update_colors = True)
            bd.CAs.colors = bd.rama_colors
            if self._update_rama_plot:
                self._rama_plot.update_scatter(bd)

    def update_omega_check(self, *_):
        rc = self._rama_counter
        ov = self.omega_validator
        # Stagger Ramachandran and omega validation to reduce jerkiness
        if self._rama_counter == self.params.rounds_per_rama_update//2:
            cis, twisted = ov.find_outliers()
            ov.draw_outliers(cis, twisted)
        ov.update_coords()


    def _get_and_check_positions(self):
        from simtk.unit import angstrom
        c = self.sim.context
        state = c.getState(getPositions = True)
        indices = self._total_mobile_indices
        old_pos = self._particle_positions
        pos = state.getPositions(asNumpy = True).value_in_unit(CHIMERAX_LENGTH_UNIT)
        delta = pos[indices] - old_pos[indices]
        distances = numpy.linalg.norm(delta, axis=1)
        max_distance = distances.max()
        max_allowed = self._max_atom_movement_per_step.value_in_unit(CHIMERAX_LENGTH_UNIT)
        if max_distance > max_allowed:
            fast_indices = numpy.where(distances > max_allowed)[0]
        else:
            fast_indices = None
        return pos, fast_indices


    def _get_positions_and_max_force (self, save_forces = False):
        c = self.sim.context
        from simtk.unit import kilojoule_per_mole, nanometer, angstrom
        state = c.getState(getForces = True, getPositions = True)
        # We only want to consider forces on the mobile atoms to decide
        # if the simulation is unstable
        forces = (state.getForces(asNumpy = True) \
            /(kilojoule_per_mole/nanometer))[self._total_mobile_indices]
        magnitudes = numpy.linalg.norm(forces, axis=1)
        #~ rstate = c.getState(Unab
            #~ getForces=True,
            #~ groups = {self._force_groups['position restraints']})
        #~ f = (rstate.getForces(asNumpy=True) \
            #~ / (kilojoule_per_mole/nanometer))[self._sim_pos_restr_indices_in_sim]
        #self._pos_restraint_forces = numpy.linalg.norm(f, axis = 1)
        if save_forces:
            self.forces = forces
            self.starting_positions = state.getPositions(asNumpy=True) / angstrom
            for i, f in enumerate(self._system.getForces()):
                print(f, c.getState(getEnergy=True, groups = {i}).getPotentialEnergy())
        max_mag = magnitudes.max()
        # Only look up the index if the maximum force is too high
        if max_mag > self._max_allowable_force:
            max_index = numpy.where(magnitudes == max_mag)[0][0]
        else:
            max_index = -1
        pos = state.getPositions(asNumpy = True)/angstrom
        return pos, max_mag, max_index

    #############
    # Commands for script/command-line control
    #############

    def set_sim_selection_mode(self, mode):
        try:
            self._sim_selection_mode = self._sim_selection_modes[mode]
        except KeyError:
            e = "mode must be one of 'from_picked_atoms', 'chain', 'whole_model', \
                    'custom', or 'script'"
            raise KeyError(e)

    def set_sim_mode(self, mode):
        try:
            self.sim_mode = self._sim_modes[mode]
        except KeyError:
            e = "Mode must be one of 'xtal', 'em' or 'free'"
            raise Exception(e)

    ##############################################
    # Warning dialogs
    ##############################################
    def _no_maps_warning(self):
        '''
        If the user has chosen crystallography or EM mode and has not
        set up any maps.
        '''
        from PyQt5.QtWidgets import QMessageBox
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        errstring = 'You have selected {} but have not selected any maps. '\
            .format(self._human_readable_sim_modes[self.sim_mode])\
            + 'Click Ok to switch to Free mode and run without maps, or '\
            + 'Cancel to stop and choose map(s).'
        msg.setText(errstring)
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        choice = msg.exec()
        if choice == QMessageBox.Ok:
            self.iw._sim_basic_mode_combo_box.setCurrentIndex(2)
            self.start_sim()
        else:
            return

    #############################################
    # Final cleanup
    #############################################
    def _on_close(self, *_):
        self.session.logger.status('Closing ISOLDE and cleaning up')

        # Remove all registered event handlers
        self._event_handler.remove_all_handlers()
        self._isolde_events.remove_all_handlers()

        # Revert mouse modes
        self._set_chimerax_mouse_mode_panel_enabled(True)
        self._mouse_modes.remove_all_modes()

    ##############################################
    # General housekeeping
    ##############################################

    def _splash_destroy_cb(self, *_):
        self._splash_destroy_countdown -= 1
        if self._splash_destroy_countdown <= 0:
            self.session.triggers.remove_handler(self._splash_handler)
            self._splash_handler = None
            self._splash.close()


def _generic_warning(message):
    msg = QMessageBox()
    ms.setIcon(QMessageBox.Warning)
    msg.setText(message)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec()


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