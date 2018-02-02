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
from chimerax.core.map import Volume
from chimerax.core.atomic import AtomicStructure

from .isolde_model import IsoldeCrystalModel, IsoldeEMModel, IsoldeFreeModel
from . import rotamers, dihedrals
from .eventhandler import EventHandler
from .constants import defaults, sim_outcomes, control
from .param_mgr import Param_Mgr, autodoc, param_properties
from .checkpoint import CheckPoint
from .openmm import sim_interface
from .openmm.sim_interface import SimParams

from .validation import Validation_Mgr

from .validation_new import validation_interface #TODO: Remove

from PyQt5.QtWidgets import QMessageBox

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
        'use_multiprocessing':                  (defaults.USE_FORK_INSTEAD_OF_THREADS, None),
            # The Shannon rate for oversampling maps from FFTs
        'map_shannon_rate':                     (defaults.MAP_SHANNON_RATE, None),

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
        _sim_modes.xtal:  "Map Fitting mode",
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
        _map_styles.solid_t20: {'style': 'surface', 'transparency': 0.2},
        _map_styles.solid_t40: {'style': 'surface', 'transparency': 0.4},
        _map_styles.solid_t60: {'style': 'surface', 'transparency': 0.6},
        _map_styles.solid_t80: {'style': 'surface', 'transparency': 0.8},
        _map_styles.solid_opaque: {'style': 'surface', 'transparency': 1.0}
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
        self.session = session = gui.session

        # Find or create the validation managers
        from . import session_extensions
        self._validation_mgr = Validation_Mgr(session)

        self.triggers = triggerset.TriggerSet()
        for t in self.trigger_names:
            self.triggers.add_trigger(t)

        self._isolde_events = EventHandler(self)
        self._logging = False
        self._log = Logger('isolde.log')

        self._sim_interface = None

        self._can_checkpoint = True
        self.checkpoint_disabled_reasons = {}
        self._last_checkpoint = None

        self.params = IsoldeParams()
        self.sim_params = SimParams()

        self._status = self.session.logger.status

        self._event_handler = EventHandler(self.session)

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
        # TODO: There needs to be a unique annotation drawing for each model
        # considered by ISOLDE
        #~ from chimerax.core.models import Model
        #~ self._annotations = Model('ISOLDE annotations', self.session)
        #~ self.session.models.add([self._annotations])

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
        self._selected_rotamer = None
        self._target_rotamer = None
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
        #~ self.omega_validator = validation.OmegaValidator(self._annotations)
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
        # A single Dihedrals array containing all proper dihedrals that
        # we want to be able to draw annotations for.
        #~ self._all_annotated_dihedrals = dihedrals.Dihedrals(drawing=self._annotations, session=self.session)


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
        splash.show()
        # Make sure the splash screen is actually shown
        for i in range(5):
            self.session.ui.processEvents()
        from PyQt5 import QtCore
        # Close the splash after 2 seconds
        QtCore.QTimer.singleShot(2000, splash.close)

        self.start_gui(gui)


    @property
    def can_checkpoint(self):
        '''Is checkpoint save/revert currently allowed?'''
        return self._can_checkpoint

    @can_checkpoint.setter
    def can_checkpoint(self, flag):
        self._can_checkpoint = flag


    @property
    def sim_interface(self):
        if self._sim_interface is None:
            raise TypeError('No simulation is currently running!')
        return self._sim_interface

    @property
    def simulation_running(self):
        return self._sim_interface is not None

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

        # FIXME Remove 'Custom restraints' tab from the gui while I decide
        # whether to keep it
        self.iw._sim_tab_widget.removeTab(1)

        self.gui_mode = True

         # Register ISOLDE-specific mouse modes
        self._mouse_modes.register_all_isolde_modes()



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

        eh = self._event_handler
        eh.add_event_handler('update_menu_on_selection',
                             'selection changed', self._selection_changed)
        eh.add_event_handler('update_menu_on_model_add',
                            'add models', self._update_model_list)
        eh.add_event_handler('update_menu_on_model_remove',
                             'remove models', self._update_model_list)
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
            if mode == self._sim_modes.em:
                #FIXME EM Mode is under development. Hide for now.
                continue
            text = self._human_readable_sim_modes[mode]
            cb.addItem(text, mode)

        iw._sim_temp_spin_box.setProperty('value',
            defaults.TEMPERATURE)

        # Populate force field combo box with available forcefields
        cb = iw._sim_force_field_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.main_file_descriptions)

        # Populate water model combo box with available models
        cb = iw._sim_water_model_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.explicit_water_descriptions)

        # Populate OpenMM platform combo box with available platforms
        # The current threaded implementation only works for the CPU
        # platform on Mac systems.
        cb = iw._sim_platform_combo_box
        cb.clear()

        platform_names = sim_interface.get_available_platforms()
        cb.addItems(platform_names)

        # Set to the fastest available platform
        if 'CUDA' in platform_names:
            cb.setCurrentIndex(platform_names.index('CUDA'))
        elif 'OpenCL' in platform_names:
            cb.setCurrentIndex(platform_names.index('OpenCL'))
        elif 'CPU' in platform_names:
            cb.setCurrentIndex(platform_names.index('CPU'))

        iw._sim_basic_mobile_b_and_a_spinbox.setProperty('value',
            defaults.SELECTION_SEQUENCE_PADDING)
        iw._sim_basic_mobile_sel_within_spinbox.setProperty('value',
            defaults.SOFT_SHELL_CUTOFF)
        iw._sim_basic_mobile_chains_within_spinbox.setProperty('value',
            defaults.SOFT_SHELL_CUTOFF)

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

        iw._sim_basic_xtal_map_cutoff_spin_box.setProperty('value',
            defaults.STANDARD_MAP_MASK_RADIUS)
        iw._sim_basic_xtal_map_weight_spin_box.setProperty('value',
            defaults.STANDARD_MAP_K)
        iw._em_map_cutoff_spin_box.setProperty('value',
            defaults.STANDARD_MAP_MASK_RADIUS)
        iw._em_map_coupling_spin_box.setProperty('value',
            defaults.STANDARD_MAP_K)
        ####
        # Rebuild tab
        ####

        ## Info for a single selected residue
        iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        iw._rebuild_sel_res_pep_info.setText('')
        iw._rebuild_sel_res_rot_info.setText('')

        iw._rebuild_pos_restraint_spring_constant.setProperty('value',
            self.sim_params.position_restraint_spring_constant.value_in_unit(CHIMERAX_SPRING_UNIT))

        ####
        # Validate tab
        ####

        # Populate the Ramachandran plot case selector with available
        # cases
        cb = iw._validate_rama_case_combo_box
        cb.clear()
        rm = self._validation_mgr.rama_mgr
        #~ from . import validation
        # First key is the null (N/A) case, which doesn't get plotted
        keys = rm.RAMA_CASES[1:]
        for key in reversed(keys):
            cb.addItem(rm.RAMA_CASE_DETAILS[key]['name'], key)


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

        iw._demo_load_button.clicked.connect(self.load_demo_data)

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
            self._change_soft_shell_cutoff_from_sel_menu
            )
        iw._sim_basic_mobile_b_and_a_spinbox.valueChanged.connect(
            self._change_b_and_a_padding
            )
        iw._sim_basic_mobile_sel_backbone_checkbox.stateChanged.connect(
            self._change_soft_shell_fix_backbone_from_sel_menu
            )
        iw._sim_basic_mobile_chains_within_spinbox.valueChanged.connect(
            self._change_soft_shell_cutoff_from_chains_menu
            )
        iw._sim_basic_mobile_chain_backbone_checkbox.stateChanged.connect(
            self._change_soft_shell_fix_backbone_from_chains_menu
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

        # Visualisation tools
        iw._sim_basic_xtal_step_forward_button.clicked.connect(
            self._xtal_step_forward
            )
        iw._sim_basic_xtal_mask_to_selection_button.clicked.connect(
            self._xtal_mask_to_selection
            )
        iw._sim_basic_xtal_live_scrolling_button.clicked.connect(
            self._xtal_enable_live_scrolling
            )
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
        iw._rebuild_sel_res_rot_release_button.clicked.connect(
            self._release_rotamer
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
            if n == 0:
                print('No haptic devices found. Stopping haptic driver.')
                hh.stopHaptics()
                self._use_haptics = False
                return

            d = self._haptic_devices = [None] * n
            self._haptic_tug_atom = [None] * n
            self._haptic_tug_index = [-1] * n
            self._haptic_highlight_nearest_atom = [True] * n
            self._haptic_current_nearest_atom = [None] * n
            self._haptic_current_nearest_atom_color = [None] * n
            self._haptic_current_nearest_atom_draw_mode = [None] * n
            from . import tugging
            for i in range(n):
                d[i] = tugging.HapticTugger(self.session, i)
            self._use_haptics = True
            self._status('')
            self._isolde_events.add_event_handler('Enable haptics during simulation',
                                                   'simulation started',
                                                   self._start_sim_haptics)
            self._isolde_events.add_event_handler('Disable haptics on sim end',
                                                  'simulation terminated',
                                                  self._stop_sim_haptics)
        else:
            self._use_haptics = False

    def _update_sim_temperature(self):
        t = self.iw._sim_temp_spin_box.value()
        self.sim_params.temperature = t
        if self.simulation_running:
            self._sim_interface.temperature = t

    ##############################################################
    # Menu control functions to run on key events
    ##############################################################
    def _update_model_list(self, *_):
        sim_mode = self.sim_mode
        modes = self._sim_modes
        mmcb = self.iw._master_model_combo_box
        mmcb.clear()
        emcb = self.iw._em_map_model_combo_box
        emcb.clear()
        xmcb = self.iw._sim_basic_xtal_init_model_combo_box
        xmcb.clear()
        potential_model_combo_box = {
            modes.xtal:    xmcb,
            modes.em:      None,
            modes.free:    None,
        }

        _models = self.session.models.list()
        #Consider top-level models only
        models = []
        for m in _models:
            if len(m.id) == 1:
                models.append(m)
        models = sorted(models, key = lambda m: m.id)
        mtd = {
            AtomicStructure: [],
            clipper.CrystalStructure: [],
            IsoldeCrystalModel: [],
            IsoldeEMModel: [],
            IsoldeFreeModel: [],
            Volume: []
        }

        for m in models:
            for mtype in mtd.keys():
                if isinstance(m, mtype):
                    mtd[mtype].append(m)

        if sim_mode == modes.xtal:
            valid_models = mtd[clipper.CrystalStructure] + mtd[IsoldeCrystalModel]
            potential_models = mtd[AtomicStructure] + mtd[IsoldeFreeModel]
        elif sim_mode == modes.free:
            valid_models = mtd[AtomicStructure] + mtd[IsoldeFreeModel]
            potential_models = []
        elif sim_mode == modes.em:
            valid_models = mtd[IsoldeEMModel]
            potential_models = mtd[AtomicStructure] + mtd[IsoldeFreeModel]

        valid_models = sorted(valid_models, key=lambda m: m.id)
        potential_models = sorted(potential_models, key=lambda m: m.id)

        for m in valid_models:
            id_str = '{}. {}'.format(m.id_string(), m.name)
            mmcb.addItem(id_str, _get_atomic_model(m))
            self._available_models[id_str] = _get_atomic_model(m)

        for m in potential_models:
            xmcb.addItem('{}. {}'.format(m.id_string(), m.name), _get_atomic_model(m))

        pmcb = potential_model_combo_box[self.sim_mode]
        if pmcb is not None:
            for m in potential_models:
                pmcb.addItem('{}. {}'.format(m.id_string(), m.name), _get_atomic_model(m))

        if sim_mode == modes.em:
            for m in mtd[Volume]:
                emcb.addItem('{}. {}'.format(m.id_string(), m.name), m)



    def _update_model_list_old(self, *_):
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
                r = selres[0]
                self._enable_rebuild_residue_frame(r)
            else:
                self._clear_rotamer()
                self._disable_rebuild_residue_frame()
            if is_continuous_protein_chain(sel):
                self._enable_secondary_structure_restraints_frame()
                self._enable_register_shift_frame()
                self._enable_selection_extend_frame()
            else:
                self._disable_secondary_structure_restraints_frame()
                self._disable_register_shift_frame()
                self._disable_selection_extend_frame()

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
        paused = self._sim_paused
        go_button = iw._sim_go_button
        if paused and not flag:
            self._sim_paused = False
            go_button.setChecked(False)
        elif paused:
            go_button.setToolTip('Resume')
        elif flag:
            go_button.setToolTip('Pause')
        if not flag:
            iw._sim_go_button.setToolTip('Start a simulation')

        iw._sim_save_checkpoint_button.setEnabled(flag)
        iw._sim_revert_to_checkpoint_button.setEnabled(flag)
        iw._sim_commit_button.setEnabled(flag)
        iw._sim_stop_and_revert_to_checkpoint_button.setEnabled(flag)
        iw._sim_stop_and_discard_button.setEnabled(flag)
        if self.simulation_type == 'equil':
            iw._sim_equil_button.setChecked(True)
        else:
            iw._sim_min_button.setChecked(True)

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
            _generic_warning(errstring)
        if not os.path.isfile(fname):
            errstring = 'Please select a valid MTZ file!'
            _generic_warning(errstring)
        m = cb.currentData()
        cs = clipper.CrystalStructure(self.session, m, mtzfile=fname,
                            map_oversampling=self.params.map_shannon_rate)
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
        rot_text = self.iw._rebuild_sel_res_rot_info

        try:
            rot = self.rotamers[res]
            if rot != self._selected_rotamer:
                self._clear_rotamer()
            self._selected_rotamer = rot
            t = rot.current_target_rotamer
            if t is None:
                rot_text.setText('')
            else:
                rot_text.setText('{} ({:.3f})'\
                    .format(t.name, t.relative_abundance(res)))
            self._selected_rotamer = rot
            self._set_rotamer_buttons_enabled(True)
        except KeyError:
            # This residue has no rotamers
            rot_text.setText('')
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
        self.iw._rebuild_2ry_struct_restr_clear_button.setEnabled(True)

    def _enable_register_shift_frame(self, *_):
        self.iw._rebuild_register_shift_container.setEnabled(True)

    def _enable_selection_extend_frame(self, *_):
        self.iw._rebuild_grow_shrink_sel_frame.setEnabled(True)

    def _disable_secondary_structure_restraints_frame(self, *_):
        self.iw._rebuild_2ry_struct_restr_container.setEnabled(False)
        self.iw._rebuild_2ry_struct_restr_clear_button.setEnabled(False)

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

    def _extend_selection_by_one_res(self, direction):
        '''
        Extends a selection by one residue in the given direction, stopping
        when it hits a chain break or the end of a chain
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

    def _restrain_selection_as_alpha_helix(self, *_):
        from chimerax.core.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self._restrain_secondary_structure(sel, 'helix')

    def _restrain_selection_as_antiparallel_beta(self, *_):
        from chimerax.core.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self._restrain_secondary_structure(sel, 'antiparallel beta')

    def _restrain_selection_as_parallel_beta(self, *_):
        from chimerax.core.atomic import selected_atoms
        sel = selected_atoms(self.session)
        self._restrain_secondary_structure(sel, 'parallel beta')


    def _restrain_secondary_structure(self, atoms, target):
        sh = self._sim_handler
        sc = self._total_sim_construct
        dihed_k = self.secondary_structure_restraints_k
        sel = atoms
        residues = sel.unique_residues
        sel = residues.atoms
        structure_type = target
        from .dihedrals import Backbone_Dihedrals
        target_phi, target_psi = Backbone_Dihedrals.standard_phi_psi_angles[structure_type]
        phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
        phi.targets = target_phi
        psi.targets = target_psi
        phi.spring_constants = dihed_k
        psi.spring_constants = dihed_k
        self.apply_dihedral_restraints(phi)
        self.apply_dihedral_restraints(psi)
        self._update_dihedral_restraints_drawing()


        dist_k = self.sim_params.distance_restraint_spring_constant
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
            for cad in cad_candidates:
                if cad is not None:
                    if -1 not in sel.indices(cad.atoms):
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
                        ond.spring_constant = dist_k
                    else:
                        ond.spring_constant = 0
                    self.apply_distance_restraint(ond)



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
        self.add_checkpoint_block(rs, 'Register shift in progress')
        self._isolde_events.add_event_handler('register shift finish check',
                                              'completed simulation step',
                                              self._check_if_register_shift_finished)


    def _check_if_register_shift_finished(self, *_):
        if self._register_shifter.finished:
            self.iw._rebuild_register_shift_release_button.setEnabled(True)
            self._isolde_events.remove_event_handler('register shift finish check')

    def _release_register_shifter(self, *_):
        rs = self._register_shifter
        if rs is not None:
            rs.release_all()
            self.remove_checkpoint_block(rs)
        self._register_shifter = None
        self.iw._rebuild_register_shift_go_button.setEnabled(True)
        self.iw._rebuild_register_shift_release_button.setEnabled(False)

    def _restrain_selected_atom_to_current_xyz(self, *_):
        from chimerax.core.atomic import selected_atoms
        atom = selected_atoms(self.session)[0]
        k = self.iw._rebuild_pos_restraint_spring_constant.value()
        self.restrain_atom_to_xyz(atom, atom.coord, k)

    def _restrain_selected_atom_to_pivot_xyz(self, *_):
        from chimerax.core.atomic import selected_atoms
        atom = selected_atoms(self.session)[0]
        k = self.iw._rebuild_pos_restraint_spring_constant.value()
        self.restrain_atom_to_xyz(atom, self.session.view.center_of_rotation, k)


    def restrain_atom_to_xyz(self, atom, target, spring_constant):
        pr = self.position_restraints[atom]
        pr.target = target
        pr.spring_constant = spring_constant
        if self.simulation_running:
            self._sim_interface.update_position_restraint(pr)


    def release_xyz_restraints_on_selected_atoms(self, *_, sel = None):
        '''
        Release current position restraints on a set of atoms. If a
        simulation is currently running, it will be automatically
        notified of the change.
        '''
        if sel is None:
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        restraints = self.position_restraints.in_selection(sel)
        restraints.release()
        if self.simulation_running:
            self._sim_interface.update_position_restraints(restraints)

    def apply_xyz_restraints_to_selected_atoms(self, *_, sel = None):
        '''
        After setting targets and/or spring constants for a set of
        position restraints, run this to apply them to the currently
        running simulation. Note: If no simulation is running this step
        is unnecessary. The restraints will be automatically applied to
        all future simulations.
        '''
        if not self.simulation_running:
            return
        if sel is None:
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        restraints = self.position_restraints.in_selection(sel)
        self._sim_interface.update_position_restraints(restraints)

    def _set_rotamer_buttons_enabled(self, switch):
        iw = self.iw
        #iw._rebuild_sel_res_last_rotamer_button.setEnabled(switch)
        #iw._rebuild_sel_res_next_rotamer_button.setEnabled(switch)
        iw._rebuild_cycle_rotamer_frame.setEnabled(switch)
        iw._rebuild_sel_res_last_rotamer_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_commit_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_target_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_discard_button.setEnabled(switch)
        iw._rebuild_sel_res_rot_release_button.setEnabled(switch)

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
        rot = self._selected_rotamer
        rot.commit_current_preview()
        if self.simulation_running:
            self._sim_interface.update_coords(rot.residue.atoms)
        self._clear_rotamer()

    def _set_rotamer_target(self, *_):
        rot = self._selected_rotamer
        tr = self._target_rotamer
        if tr is None:
            print('No target set!')
            return
        target = rot.target = tr.angles
        rot.spring_constant = \
            self.sim_params.rotamer_spring_constant
        rot.restrained = True
        if self.simulation_running:
            self._apply_rotamer_target_to_sim(rot)
        self._clear_rotamer()
        self._update_dihedral_restraints_drawing()

    def _apply_rotamer_target_to_sim(self, rotamer):
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.update_rotamer_target(rotamer)
        self._all_annotated_dihedrals.update_needed = True

    def _clear_rotamer(self, *_):
        if self._selected_rotamer is not None:
            self._selected_rotamer.cleanup()
        self._target_rotamer = None

    def _release_rotamer(self, *_):
        rot = self._selected_rotamer
        rot.restrained = False
        rot.cleanup()
        if self.simulation_running:
            self._apply_rotamer_target_to_sim(rot)
        self._update_dihedral_restraints_drawing()

    def release_rotamers(self, residues):
        for r in residues:
            self.release_rotamer(r)

    def release_rotamer(self, residue):
        try:
            rot = self.rotamers[residue]
        except KeyError:
            return
        rot.restrained = False
        if self.simulation_running:
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
        #TODO
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

    def _flip_cis_trans(self, *_):
        res = self._rebuild_residue
        self.flip_peptide_omega(res)



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
        res = self._total_mobile.unique_residues
        self._rama_plot.set_target_residues(res)
        self._update_rama_plot = True
        self.iw._validate_rama_sel_combo_box.setDisabled(True)
        self.iw._validate_rama_go_button.setDisabled(True)

    def _rama_go_static(self, *_):
        self._update_rama_plot = False
        self.iw._validate_rama_sel_combo_box.setEnabled(True)
        self.iw._validate_rama_go_button.setEnabled(True)

    def _rama_static_plot(self, *_):
        model = self._selected_model
        rplot = self._rama_plot
        whole_model = bool(self.iw._validate_rama_sel_combo_box.currentIndex())
        if whole_model:
            res = model.residues
            rplot.update_scatter(residues = res)
            # sel = model.atoms
            # bd = self.backbone_dihedrals
            # self._rama_plot.update_scatter(bd, force_update = True)

        else:
            sel = model.atoms.filter(model.atoms.selected)
            residues = sel.unique_residues
            rplot.set_target_residues(residues)
            rplot.update_scatter()
            # if not len(residues):
            #
            #     rplot.update_scatter(res)
            # else:
            #     phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
            #     from . import dihedrals
            #     bd = dihedrals.Backbone_Dihedrals(self.session, phi=phi, psi=psi, omega=omega)
            #     self._rama_plot.update_scatter(bd, force_update = True)
            # else:
            #     self._rama_plot.update_scatter(force_update = True)

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
            #self._prepare_all_interactive_restraints(m)
            self._update_chain_list()
            if isinstance(m.parent, clipper.CrystalStructure):
                self._initialize_xtal_maps(m.parent)
                iw._sim_basic_xtal_stepper_frame.setEnabled(True)
            else:
                iw._sim_basic_xtal_stepper_frame.setEnabled(False)
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
        from . import backbone_restraints as br
        from . import position_restraints as pr
        # Torsion restraints
        bd = self.backbone_dihedrals = dihedrals.Backbone_Dihedrals(self.session, m)
        ad = self._all_annotated_dihedrals
        ad.append(bd.phi)
        ad.append(bd.psi)
        self.rotamers = rotamers.all_rotamers_in_selection(self.session, m.atoms)
        for r in self.rotamers.values():
            try:
                ad.append(r.dihedrals)
            except:
                raise Exception('{}: {}'.format(r.residue.name, r.dihedrals))
        # Distance restraints
        self.ca_to_ca_plus_two = br.CA_to_CA_plus_Two(self.session, m)
        self.o_to_n_plus_four = br.O_to_N_plus_Four(self.session, m)
        # Positional restraints (one per heavy atom)
        if self.position_restraints is not None:
            self.position_restraints.cleanup()
        self.position_restraints = pr.Atom_Position_Restraints(
                self.session, m,
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

    def _change_soft_shell_cutoff_from_sel_menu(self, *_):
        iw = self.iw
        val = iw._sim_basic_mobile_sel_within_spinbox.value()
        sb2 = iw._sim_basic_mobile_chains_within_spinbox
        if sb2.value() != val:
            sb2.setValue(val)
        self.params.soft_shell_cutoff_distance = val

    def _change_soft_shell_cutoff_from_chains_menu(self, *_):
        iw = self.iw
        val = iw._sim_basic_mobile_chains_within_spinbox.value()
        sb2 = iw._sim_basic_mobile_sel_within_spinbox
        if sb2.value() != val:
            sb2.setValue(val)
        self.params.soft_shell_cutoff_distance = val

    def _change_soft_shell_fix_backbone_from_sel_menu(self, *_):
        mobile = self.iw._sim_basic_mobile_sel_backbone_checkbox.checkState()
        self.params.fix_soft_shell_backbone = not mobile
        cb2 = self.iw._sim_basic_mobile_chain_backbone_checkbox
        if cb2.checkState() != mobile:
            cb2.setChecked(mobile)

    def _change_soft_shell_fix_backbone_from_chains_menu(self, *_):
        mobile = self.iw._sim_basic_mobile_chain_backbone_checkbox.checkState()
        self.params.fix_soft_shell_backbone = not mobile
        cb2 = self.iw._sim_basic_mobile_sel_backbone_checkbox
        if cb2.checkState() != mobile:
            cb2.setChecked(mobile)


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
    # Visualisation functions
    ##############################################################

    def _xtal_step_forward(self, *_):
        m = self.selected_model
        cs = m.parent
        focus = self.iw._sim_basic_xtal_step_focus_checkbox.checkState()
        m.atoms.selected = False
        sel = cs.stepper.step_forward()
        sel.selected = True
        self._xtal_mask_to_atoms(sel, focus)

    def _xtal_mask_to_selection(self, *_):
        atoms = self.selected_model.atoms
        sel = atoms[atoms.selecteds]
        self._xtal_mask_to_atoms(sel, False)

    def _xtal_mask_to_atoms(self, atoms, focus):
        m = self.selected_model
        cs = m.parent
        cutoff = self.params.standard_map_mask_cutoff
        context = self.params.soft_shell_cutoff_distance
        cs.isolate_and_cover_selection(
            atoms, 0, context, cutoff, focus=focus)


    def _xtal_enable_live_scrolling(self, *_):
        m = self.selected_model
        cs = m.parent
        cs.xmaps.live_scrolling = True
        cs.live_atomic_symmetry = True

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

    def _start_sim_or_toggle_pause(self, *_):
        if not self._simulation_running:
            self.start_sim()
        else:
            self.pause_sim_toggle()

    def start_sim(self):
        self.sim_params.platform = self.iw._sim_platform_combo_box.currentText()
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
        all_a = sel_model.residues.atoms
        tm_i = all_a.indices(total_mobile)
        sc_i = all_a.indices(sc)
        tm_i.sort()
        sc_i.sort()
        total_mobile = self._total_mobile = all_a[tm_i]
        sc = self._total_sim_construct = all_a[sc_i]

        self._total_mobile_indices = sc.indices(total_mobile)

        sb = self._total_sim_bonds = sc.intra_bonds
        surr = self._surroundings = all_a.subtract(sc)

        # Cache all the colors so we can revert at the end of the simulation
        self._original_atom_colors = sc.colors
        self._original_atom_draw_modes = sc.draw_modes
        self._original_bond_radii = sb.radii
        self._original_atom_radii = sc.radii
        self._original_display_state = sel_model.atoms.displays
        self._original_ribbon_state = sel_model.residues.ribbon_displays

        if self.sim_mode == sm.xtal:
            sel_model.parent.isolate_and_cover_selection(
                total_mobile, include_surrounding_residues = 0,
                show_context = self.params.hard_shell_cutoff_distance,
                mask_radius = 4, extra_padding = 10,
                hide_surrounds = self.params.hide_surroundings_during_sim, focus = False)

        # else:
        surr.hides |= control.HIDE_ISOLDE
        surr.displays = False
        self._surroundings_hidden = True

        sel_model.residues.ribbon_displays = False

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
                    fixed_flags[i] = True


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

        from . import rotamers
        mobile_residues = total_mobile.unique_residues
        mobile_rotamers = self._mobile_rotamers = rotamers.Rotamers(
                        self.session, mobile_residues, self.rotamers)

        from .openmm.sim_interface import ChimeraXSimInterface
        sp = self.sim_params
        si = self._sim_interface = ChimeraXSimInterface(self.session, self)
        si.start_sim_thread(sp, sc, tuggable_atoms, fixed_flags, bd, mobile_rotamers,
                            distance_restraints, position_restraints, self.master_map_list)
        self._last_checkpoint = si.starting_checkpoint

        #~ if self.params.track_rotamer_status:
            #~ vi = self.live_validation_interface
            #~ vi.start_validation_threads(mobile_rotamers)


        self._status('Simulation running')

    def _start_sim_haptics(self, *_):
        self._event_handler.add_event_handler('sim haptic update',
                                              'new frame',
                                              self._update_haptics)

    def _stop_sim_haptics(self, *_):
        self._event_handler.remove_event_handler('sim haptic update')

    def _update_haptics(self, *_):
        hh = self.session.HapticHandler
        si = self._sim_interface
        from . import picking
        for d in self._haptic_devices:
            i = d.index
            pointer_pos = hh.getPosition(i, scene_coords = True)
            if self._haptic_highlight_nearest_atom[i] and not d.tugging:
                cur_a = self._haptic_current_nearest_atom[i]
                new_a = self._haptic_current_nearest_atom[i] = picking.pick_closest_to_point(
                        self.session, pointer_pos, self._total_mobile, 3.0,
                        displayed_only = True, tug_hydrogens = self.tug_hydrogens)
                if new_a != cur_a:
                    if cur_a is not None and not cur_a.deleted:
                        # set the previously picked atom back to standard
                        cur_a.draw_mode = self._haptic_current_nearest_atom_draw_mode[i]
                        cur_a.color = self._haptic_current_nearest_atom_color[i]
                    if new_a is not None:
                        self._haptic_current_nearest_atom_draw_mode[i] = new_a.draw_mode
                        new_a.draw_mode = 1
                        self._haptic_current_nearest_atom_color[i] = new_a.color
                        new_a.color = [0, 255, 255, 255] # bright cyan


            b = hh.getButtonStates(i)
            # Button 0 controls tugging
            if b[0]:
                if not d.tugging:
                    a = self._haptic_tug_atom[i] = cur_a
                    if a is not None:
                        d.start_tugging(a)
                        target = hh.getPosition(i, scene_coords = True)
                        si.tug_atom_to(a, target)

                else:
                    d.update_target()
                    a = self._haptic_tug_atom[i]
                    si.tug_atom_to(a, hh.getPosition(i, scene_coords = True))
            else:
                if d.tugging:
                    d.stop_tugging()
                    a = self._haptic_tug_atom[i]
                    si.release_tugged_atom(a)
                    self._haptic_tug_atom[i] = None



    def _sim_start_cb(self, *_):
        '''
        Register all event handlers etc. that have to be running during the
        simulation.
        '''
        self._simulation_running = True
        self._update_sim_control_button_states()

        vm = self._validation_mgr
        r = self._total_mobile.unique_residues
        vm.start_tracking(r)
        self._isolde_events.add_event_handler('live validation during sim',
                                              'completed simulation step',
                                              vm.update)

        # if self.params.track_ramachandran_status:
        #     rm = self._validation_mgr.rama_mgr
        #     tm = self._total_mobile
        #     r = self._mobile_res = tm.unique_residues
        #     rcas = self._mobile_CA_atoms = r.atoms[r.atoms.names=='CA']
        #     default_color = rm.color_scale[-1]
        #     rcas.colors = default_color
        #     if self._rama_plot is not None:
        #         self._rama_go_live()

        self._isolde_events.add_event_handler('rezone maps during sim',
                                              'completed simulation step',
                                              self._rezone_maps_if_required)




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
            self.revert_to_checkpoint()

        elif outcome[0] == sim_outcomes.TIMEOUT:
            raise outcome[1]

        elif outcome[0] == sim_outcomes.FAIL_TO_START:
            '''
            TODO: Pop up a dialog box giving the details of the exception, then
            clean up.
            '''
            print('Failed to start')
            pass
        # Otherwise just clean up
        si = self._sim_interface
        xeh = self._event_handler
        ieh = self._isolde_events
        self._simulation_running = False
        ieh.remove_event_handler(
            'rezone maps during sim', error_on_missing = False)
        ieh.remove_event_handler(
            'live validation during sim', error_on_missing = False)


        #~ if self.params.track_rotamer_status:
            #~ self.live_validation_interface.stop_validation_threads()
        self._update_dihedral_restraints_drawing()

        self._disable_rebuild_residue_frame()

        self._update_menu_after_sim()
        self._mouse_tugger = None
        for d in self._haptic_devices:
            d.cleanup()
        sel_model = self._selected_model
        sc = self._total_sim_construct
        sb = self._total_sim_bonds
        sel_model.atoms.displays = self._original_display_state
        sel_model.residues.ribbon_displays = self._original_ribbon_state
        self._surroundings.hides &= ~control.HIDE_ISOLDE
        self._surroundings_hidden = False
        sc.colors = self._original_atom_colors
        sc.draw_modes = self._original_atom_draw_modes
        sc.radii = self._original_atom_radii
        sb.radii = self._original_bond_radii
        self._total_sim_construct = None

        self._surroundings = None
        if self.params.track_ramachandran_status:
            # Update one last time
            if self._update_rama_plot:
                self._rama_plot.update_scatter(self._mobile_backbone_dihedrals)
                self.update_omega_check()
        self._rama_go_static()
        self.iw._rebuild_sel_residue_frame.setDisabled(True)
        self.omega_validator.load_structure(self.selected_model, self.backbone_dihedrals.omega)
        self.update_omega_check(force=True)

        self._sim_is_unstable = False
        if self.sim_mode == self._sim_modes['xtal']:
            cs = self._selected_model.parent
            cs.xmaps.live_scrolling = True
            cs.live_atomic_symmetry = True
        self._status('')
        self._sim_interface = None
        self.equilibrate()
        self.iw._sim_go_button.setChecked(False)

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
        sm = self._selected_model
        all_atoms = sm.atoms
        if mode == modes.chain or mode == modes.whole_model:
            # Then everything is easy. The selection is already defined
            # FIXME: Amend this pipeline to allow better command-line
            # control
            return self._selected_atoms
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
            if sm != us[0]:
                raise Exception('Selection must be in the model currently chosen for ISOLDE!')
            # selatoms_by_chain = selatoms.by_chain
            # selchains = [row[1] for row in selatoms_by_chain]
            all_res = sm.residues
            sel_res = selatoms.unique_residues
            sel_res_indices = all_res.indices(sel_res)
            from chimerax.core.atomic.structure import Structure

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
        if self.simulation_running:
            self._sim_interface.toggle_pause()

    def _sim_pause_cb(self, *_):
        self._sim_paused = True
        self._status('Simulation paused')
        go_button = self.iw._sim_go_button
        go_button.setChecked(False)
        go_button.setToolTip('Resume')

    def _sim_resume_cb(self, *_):
        self._sim_paused = False
        self._status('Simulation running')
        go_button = self.iw._sim_go_button
        go_button.setChecked(True)
        go_button.setToolTip('Pause')

    def _stop_sim_and_revert_to_checkpoint(self, *_):
        self.discard_sim('checkpoint')

    def _discard_sim(self, *_):
        self.discard_sim('start')



    def discard_sim(self, revert_to='checkpoint'):
        '''
        Stop the simulation and revert to either the starting state or
        the last saved checkpoint.
        Args:
            revert_to (default: 'checkpoint'):
                Either 'checkpoint' or 'start'
        '''

        if not self.simulation_running:
            print('No simulation running!')
            return
        self._release_register_shifter()
        if revert_to == 'start':
            msg = 'All changes since you started this simulation will be '\
                +'lost! Are you sure you want to continue?'
            ok = _choice_warning(msg)
            if ok:
                self._sim_interface.stop_sim(sim_outcomes.DISCARD, None)
            return
        elif revert_to == 'checkpoint':
            print('Stopping and reverting to the last saved checkpoint.')
            self.revert_to_checkpoint()
            self._sim_interface.stop_sim(sim_outcomes.COMMIT, None)
        else:
            raise TypeError('Unrecognised option! Argument should be '\
                +'either "checkpoint" or "start".')


    def commit_sim(self):
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._release_register_shifter()
        self._sim_interface.stop_sim(sim_outcomes.COMMIT, None)

    def minimize(self):
        print('Minimisation mode')
        if self.simulation_running:
            self._sim_interface.sim_mode = 'min'
        self.simulation_type = 'min'
        self._update_sim_control_button_states()

    def equilibrate(self):
        print('Equilibration mode')
        if self.simulation_running:
            self._sim_interface.sim_mode = 'equil'
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
            atom:
                The atom to be tugged. Must be a heavy atom that is
                mobile in the current simulation
            target:
                An (x,y,z) Numpy array (in Angstroms)
            spring constant (default None):
                An optional spring constant (in kJ mol-1 A-2)
        '''
        if not self.simulation_running:
            return
        self._sim_interface.tug_atom_to(atom, target, spring_constant)

    def stop_tugging(self, atom):
        if not self.simulation_running:
            return
        self._sim_interface.release_tugged_atom(atom)

    def add_checkpoint_block(self, obj, reason):
        '''
        Some processes are incompatible with checkpointing. We need to
        know when these are running and disable checkpointing until
        they're done.
        '''
        self.can_checkpoint = False
        self.checkpoint_disabled_reasons[obj] = reason
        self.iw._sim_save_checkpoint_button.setEnabled(False)
        self.iw._sim_revert_to_checkpoint_button.setEnabled(False)

    def remove_checkpoint_block(self, obj):
        '''
        Release a block on checkpointing. When all blocks are removed,
        checkpointing will be re-enabled.
        '''
        r = self.checkpoint_disabled_reasons
        r.pop(obj)
        if len(r) ==0:
            self.can_checkpoint = True
            self.iw._sim_save_checkpoint_button.setEnabled(True)
            self.iw._sim_revert_to_checkpoint_button.setEnabled(True)



    def checkpoint(self, *_):
        if self.can_checkpoint:
            self._last_checkpoint = CheckPoint(self)
        else:
            err_str = 'Checkpointing is currently disabled by the '\
                +'following scripts and will be re-enabled when they '\
                +'terminate: \n{}'.format(
                    '\n'.join([r for r in self.checkpoint_disabled_reasons.values()]))
            raise TypeError(err_str)

    def revert_to_checkpoint(self, *_):
        if self.can_checkpoint:
            if self._last_checkpoint is None:
                raise TypeError('No saved checkpoint available!')
            self._last_checkpoint.revert()
            self._update_dihedral_restraints_drawing()
        else:
            err_str = 'Checkpointing is currently disabled by the '\
                +'following scripts and will be re-enabled when they '\
                +'terminate: \n{}'.format(
                    '\n'.join([r for r in self.checkpoint_disabled_reasons.values()]))
            raise TypeError(err_str)


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
        sc = self._total_sim_construct
        k = self.sim_params.peptide_bond_spring_constant
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
            t = (numpy.logical_or(dvals > cr[1], dvals < cr[0])).astype(float) * pi

        dihedrals.targets = t
        dihedrals.spring_constants = k
        if self.simulation_running:
            self._sim_interface.update_dihedral_restraints(dihedrals)

    def remove_peptide_bond_restraints(self, dihedrals):
        '''
        Remove restraints on a list of peptide bonds (actually, just set
        their spring constants to zero). Simulation must already be running.
        '''
        if self.simulation_running:
            self._sim_interface.update_dihedral_restraints(dihedrals)

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
            d.spring_constant = self.sim_params.phi_psi_spring_constant

        self.apply_dihedral_restraint(phi)
        self.apply_dihedral_restraint(psi)
        self._update_dihedral_restraints_drawing()

        self._pep_flip_timeout_counter = 0
        self._pep_flip_polish_counter = 0
        self._pep_flip_dihedrals = (phi, psi)
        self.iw._rebuild_sel_res_pep_flip_button.setEnabled(False)
        self._isolde_events.add_event_handler('pep flip timeout',
                                              'completed simulation step',
                                              self._check_pep_flip)

    def _check_pep_flip(self, *_):
        if self._pep_flip_timeout_counter * self.sim_params.sim_steps_per_gui_update < 100:
            self._pep_flip_timeout_counter += 1
            # Need to give it some time to settle first
            return
        done = False
        if self._pep_flip_timeout_counter * self.sim_params.sim_steps_per_gui_update >= 1000:
            print('Unable to flip peptide. Giving up.')
            done = True
        dihedrals = self._pep_flip_dihedrals
        if not done:
            done = True
            for d in dihedrals:
                diff = abs(d.value-d.target)
                if diff > pi:
                    diff -= 2*pi
                if abs(diff)*OPENMM_ANGLE_UNIT > self.sim_params.dihedral_restraint_cutoff_angle:
                    done = False
        if not done:
            self._pep_flip_timeout_counter += 1
            return
        else:
            self._isolde_events.remove_event_handler('pep flip timeout')
            self._isolde_events.add_event_handler('pep flip polish',
                                            'completed simulation step',
                                            self._polish_pep_flip)

    def _polish_pep_flip(self, *_):
        if self._pep_flip_polish_counter == 0:
            self._sim_interface.sim_mode = 'min'
            self._pep_flip_polish_counter += 1
        else:
            dihedrals = self._pep_flip_dihedrals
            self.release_dihedral_restraint(dihedrals[0])
            self.release_dihedral_restraint(dihedrals[1])
            self._sim_interface.sim_mode = 'equil'
            self.iw._rebuild_sel_res_pep_flip_button.setEnabled(True)
            self._isolde_events.remove_event_handler('pep flip polish')

    def apply_dihedral_restraint(self, dihedral):
        '''
        Apply the current restraint parameters for a dihedral to the
        currently-running simulation.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.update_dihedral_restraint(dihedral)
        self._all_annotated_dihedrals.update_needed = True

    def apply_dihedral_restraints(self, dihedrals):
        '''
        Apply restraints for a set of dihedrals at once.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
        all_names = dihedrals.names
        unique_names = numpy.unique(all_names)
        if len(unique_names) == 1:
            self._sim_interface.update_dihedral_restraints(dihedrals)
        else:
            for name in unique_names:
                indices = numpy.where(all_names == name)[0]
                self._sim_interface.update_dihedral_restraints(dihedrals[indices])
        self._all_annotated_dihedrals.update_needed = True

    def release_dihedral_restraint(self, dihedral):
        dihedral.target = 0
        dihedral.spring_constant = 0
        if self.simulation_running:
            self.apply_dihedral_restraint(dihedral)
        else:
            self._update_dihedral_restraints_drawing()

    def clear_secondary_structure_restraints_for_selection(self, *_, atoms = None, residues = None):
        '''
        Clear all secondary structure restraints in the selection. If
        no atoms or residues are provided, restraints will be cleared
        for any atoms selected in the main window.
        '''
        from chimerax.core.atomic import selected_atoms
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
                if self.simulation_running:
                    self.apply_dihedral_restraints(dlist)
        self._update_dihedral_restraints_drawing()

        for r in residues:
            for key, rlist in self._sim_distance_restraint_dict.items():
                drs = rlist[r]
                for dr in drs:
                    if dr is not None:
                        dr.spring_constant = 0
                        if self.simulation_running:
                            self.apply_distance_restraint(dr)


    def apply_distance_restraints(self, distance_restraints):
        '''
        After changing the target distances and/or spring constants of a
        set of distance restraints, call this function to apply them to
        a currently running simulation.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.update_distance_restraints(distance_restraints)

    def apply_distance_restraint(self, distance_restraint):
        '''
        After changing the target distance and/or spring constant of a
        distance restraint, call this function to apply the changes to
        a currently running simulation.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
        self._sim_interface.update_distance_restraint(distance_restraint)

    def flip_peptide_omega(self, res):
        '''
        Flip the peptide bond N-terminal to this residue from cis to
        trans or vice versa. Only usable when a simulation is running.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
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



    #############################################
    # Main simulation functions to be run once per coordinate update
    #############################################


    def _rezone_maps_if_required(self, *_):
        self._map_rezone_counter += 1
        self._map_rezone_counter %= self.params.rounds_per_map_remask
        if self._map_rezone_counter == 0 and self.params.remask_maps_during_sim:
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
            mgr = self._validation_mgr.rama_mgr
            #rv = self.rama_validator
            #bd = self._mobile_backbone_dihedrals
            #rv.get_scores(bd, update_colors = True)
            bd.CAs.colors = bd.rama_colors
            if self._update_rama_plot:
                self._rama_plot.update_scatter(bd)

    def update_omega_check(self, *_, force = False):
        rc = self._rama_counter
        ov = self.omega_validator
        # Stagger Ramachandran and omega validation to reduce jerkiness
        if force or self._rama_counter == self.params.rounds_per_rama_update//2:
            cis, twisted = ov.find_outliers()
            ov.draw_outliers(cis, twisted)
        ov.update_coords()

    def _update_dihedral_restraints_drawing(self):
        d = self._all_annotated_dihedrals
        d.update_graphics(update_needed = True)



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


    ##############################################
    # Demo
    ##############################################

    def load_demo_data(self):
        from chimerax.core.commands import open
        data_dir = os.path.join(self._root_dir, 'demo_data', '3io0')
        before_struct = open.open(self.session, os.path.join(data_dir, 'before.pdb'))[0]
        from chimerax.core.commands import color
        color.color(self.session, before_struct, color='bychain', target='ac')
        color.color(self.session, before_struct, color='byhetero', target='a')
        before_cs = clipper.CrystalStructure(self.session, before_struct,
            mtzfile=os.path.join(data_dir, 'before_maps.mtz'),
            map_oversampling = self.params.map_shannon_rate)
        #~ from chimerax.core.commands import lighting
        #~ lighting.lighting(self.session, depth_cue=True)
        sharp_map = before_cs.xmaps['2FOFCWT_sharp, PH2FOFCWT_sharp']
        sd = sharp_map.mean_sd_rms()[1]
        styleargs= self._map_style_settings[self._map_styles.solid_t60]
        from chimerax.core.map import volumecommand
        volumecommand.volume(self.session, [sharp_map], **styleargs)
        sharp_map.set_parameters(surface_levels = (2.5*sd,))

        from chimerax.clipper import crystal
        crystal.set_to_default_cartoon(self.session)
        self._change_selected_model(force=True)
        from . import view
        view.focus_on_selection(self.session, self.session.main_view, before_struct.atoms)


def _generic_warning(message):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText(message)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec()

def _choice_warning(message):
    '''
    Pop up a warning dialog box with the given message, and return True
    if the user wants to go ahead.
    '''
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText(message)
    msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
    reply = msg.exec()
    if reply == QMessageBox.Ok:
        return True
    return False

def _get_atomic_model(m):
    if hasattr(m, 'master_model'):
        am = m.master_model
    else:
        am = m
    return am



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
