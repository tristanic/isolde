# vim: set expandtab shiftwidth=4 softtabstop=4:

# ISOLDE: Interactive Structure Optimisation by Local Direct Exploration
# Copyright: 2016
# Author:    Tristan Croll
#            Cambridge Institute for Medical Research
#            University of Cambridge

from PyQt5.QtCore import QObject, pyqtSignal
import chimerax
from math import inf

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
    
    # Possible outcomes of simulation
    class sim_outcome(IntEnum):
        COMMIT          = 0
        DISCARD         = 1
        FAIL_TO_START   = 2
        UNSTABLE        = 3
    
    '''
    Named triggers to simplify handling of changes on key ISOLDE events,
    using the same code as ChimeraX's session.triggers. The same
    functionality could be achieved without too much trouble using the
    PyQt signal/slot system.
    '''    
    from chimerax.core import triggerset
    triggers = triggerset.TriggerSet()
    trigger_names = [
        'simulation started',    # Successful initiation of a simulation
        'simulation terminated', # Simulation finished, crashed or failed to start
        'selected model changed', # Changed the master model selection
        ]
    for t in trigger_names:
        triggers.add_trigger(t)
    
    
    
    def __init__(self, gui):
        self._logging = False
        self._log = Logger('isolde.log')
        
        self.session = gui.session
        
        self._status = self.session.logger.status

        from .eventhandler import EventHandler
        self._event_handler = EventHandler(self.session)
        
        initialize_openmm()
 
        # Available pre-defined colors
        from chimerax.core import colors
        self._available_colors = colors.BuiltinColors
        # Remove duplicates
        for key in self._available_colors:
            stripped_key = key.replace(" ","")
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
        # Number of residues before and after each selected residue to add
        # to the mobile selection
        self.b_and_a_padding = 5
        # Extra mobile shell of surrounding atoms to provide a soft buffer to
        # the simulation. Whole residues only.
        self._soft_shell_atoms = None
        # User-definable distance cutoff to define the soft shell
        self.soft_shell_cutoff = 5      # Angstroms
        # Do we fix the backbone atoms of the soft shell?
        self.fix_soft_shell_backbone = False
        # Shell of fixed atoms surrounding all mobile atoms to maintain 
        # the context of the simulation. Whole residues only.
        self._hard_shell_atoms = None
        # User-definable distance cutoff to define the hard shell
        self.hard_shell_cutoff = 8      # Angstroms
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
        # Do we want to hide surroundings during the simulation?
        self.hide_surroundings_during_sim = True
        
        ####
        # Settings for live tracking of structural quality
        ####
        
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
        
        
        # Do we track and display residues' status on the Ramachandran plot?
        self.track_rama = True
        # How often do we want to update our Ramachandran statistics?
        self.steps_per_rama_update = 5
        # Internal counter for Ramachandran update
        self._rama_counter = 0
        
        
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
        self._haptic_force_constant = 2500
        
        ####
        # Settings for handling of maps
        ####
        # List of currently available volumetric data sets
        self._available_volumes = {}
        # Master list of maps and their parameters
        self.master_map_list = {}
        # Are we adding a new map to the simulation list?
        self._add_new_map = True
        
        
        ####
        # Restraints settings
        ####
        self.restrain_peptide_bonds = True
        self.peptide_bond_restraints_k = 200
        # Range of dihedral values which will be interpreted as a cis peptide
        # bond (-30 to 30 degrees). If restrain_peptide_bonds is True, anything
        # outside of this range at the start of the simulation will be forced
        # to trans.
        from math import pi 
        self.cis_peptide_bond_range = (-pi/6, pi/6)
        
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
        self.sim_steps_per_update = 50
        # Number of steps per GUI update in minimization mode
        self.min_steps_per_update = 100
        # If using the VariableLangevinIntegrator, we define a tolerance
        self._integrator_tolerance = 0.0001
        # ... otherwise, we simply set the time per step
        self._sim_time_step = 1.0*unit.femtoseconds
        # Type of integrator to use. Should give the choice in the expert level
        # of the menu. Variable is more stable, but simulated time per gui update
        # is harder to determine
        self._integrator_type = 'fixed'
        # Constraints (e.g. rigid bonds) need their own tolerance
        self._constraint_tolerance = 0.001
        # Friction term for coupling to heat bath. Using a relatively high
        # value helps keep the simulation under control in response to
        # very large forces.
        self._friction = 5.0/unit.picoseconds
        # Limit on the net force on a single atom to detect instability and
        # force a minimisation
        self._max_allowable_force = 10000.0 # kJ mol-1 nm-1
        # We need to store the last measured maximum force to determine
        # when minimisation has converged.
        self._last_max_force = inf
        # Flag for unstable simulation
        self._sim_is_unstable = False
        # Counter for the number of simulation rounds that have been unstable
        self._unstable_min_rounds = 0
        # Maximum number of rounds to attempt minimisation of an unstable
        # simulation
        self.max_unstable_rounds = 20
        
        # Are we currently tugging on an atom?
        self._currently_tugging = False
        # Placeholder for tugging forces
        self._tugging_force = None
        # Force constant for mouse/haptic tugging. Need to make this user-adjustable
        self.tug_force_constant = 10000 # kJ/mol/nm
        
        
        
        # Holds the current simulation mode, to be chosen from the GUI
        # drop-down menu or loaded from saved settings.
        self.sim_mode = None
        # Do we have a simulation running right now?
        self._simulation_running = False
        # If running, is the simulation in startup mode?
        self._sim_startup = True
        # Maximum number of rounds of minimisation to run on startup
        self._sim_startup_rounds = 10
        # Counter for how many rounds we've done on startup
        self._sim_startup_counter = 0
        
        # Simulation temperature in Kelvin
        self.simulation_temperature = 100.0
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
        
        self.initialize_haptics()
        
        ####
        # Internal trigger handlers
        ####
        self.triggers.add_handler('simulation terminated', self._cleanup_after_sim)
        
        self.gui_mode = False
        
        from PyQt5.QtGui import QPixmap
        from PyQt5.QtWidgets import QSplashScreen
        from PyQt5.QtCore import Qt
        import os
        
        splash_pix = QPixmap(os.path.join(
            self._root_dir,'resources/isolde_splash_screen.jpg'))
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.showMessage('\n\n\n\nAnalysing your structure. Please be patient...', 
                alignment = Qt.AlignRight | Qt.AlignVCenter)
        splash.show()
        # Make sure the splash screen is actually shown
        for i in range(5):
            self.session.ui.processEvents()
        from time import sleep
        sleep(1)
        self.start_gui(gui)
        
        splash.destroy()
    
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
        
        # ISOLDE will need to define a number of custom mouse modes to
        # give a familiar environment to users of existing packages. 
        # These will not work well with the standard ChimeraX right mouse
        # button selection panel, so best if we disable it for the
        # duration of the ISOLDE session. We'll put them back the way
        # we found them when the ISOLDE gui is closed.
        self._set_chimerax_mouse_mode_panel_enabled(False)
        # Register ISOLDE-specific mouse modes
        self._mouse_modes.register_all_isolde_modes()
        
        
        
        # Function to remove all event handlers, mouse modes etc. when
        # ISOLDE is closed, and return ChimeraX to its standard state.
        self.gui.tool_window.ui_area.destroyed.connect(self._on_close)
        
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
        iw._rebuild_sel_res_pep_flip_button.clicked.connect(
            self._flip_peptide_bond
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
            self.session.view.center_of_rotation_method = 'fixed'
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
        models = self.session.models.list()
        atomic_model_list = []
        atomic_model_name_list = []
        volume_model_list = []
        volume_model_name_list = []
        sorted_models = sorted(models, key=lambda m: m.id)
        if len(sorted_models) != 0:
            # Find atomic and volumetric models and sort them into the
            # appropriate lists
            for i, m in enumerate(sorted_models):
                if m.atomspec_has_atoms():
                    id_str = m.id_string() + ' ' + m.name
                    self._available_models[id_str] = m
                    atomic_model_name_list.append(id_str)
                    atomic_model_list.append(m)
                elif hasattr(m, 'grid_data'):
                    id_str = m.id_string() + ' ' + m.name
                    self._available_volumes[id_str] = m
                    volume_model_name_list.append(id_str)
                    volume_model_list.append(m)
                else:
                    # This is a model type we don't currently handle. Ignore.
                    continue
        for l, m in zip(atomic_model_name_list, atomic_model_list):
            self.iw._master_model_combo_box.addItem(l, m)
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
        sel = selected_atoms(self.session)
        selres = sel.unique_residues
        if self._simulation_running:
            if len(selres) == 1:
                self._enable_rebuild_residue_frame(selres[0])
            else:
                self._disable_rebuild_residue_frame()
            # A running simulation takes precedence for memory control
            return
        
        if self.session.selection.empty():
            flag = False
        else:
            flag = True
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
        
    def _show_em_map_chooser(self, *_):
        self.iw._em_map_chooser_frame.show()
        self.iw._sim_basic_em_map_button.setEnabled(False)

    def _hide_em_map_chooser(self, *_):
        self.iw._em_map_chooser_frame.hide()
        self.iw._sim_basic_em_map_button.setEnabled(True)
    
    def _show_em_map_in_menu_or_add_new(self, *_):
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
            iw._em_map_model_combo_box.setCurrentText(vol.id_string())
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
        self._get_rotamer_list_for_selected_residue()
        bd = self._mobile_backbone_dihedrals
        self._rebuild_residue = res
        omega = bd.omega.by_residue(res)
        if omega is None:
            # This is a terminal residue without an omega dihedral
            self.iw._rebuild_sel_res_pep_info.setText('N/A')
            self.iw._rebuild_sel_res_pep_flip_button.setDisabled(True)
            self._rebuild_res_update_omega = False
            return
        self._rebuild_res_update_omega = True
        self._rebuild_res_omega = omega
        self.iw._rebuild_sel_res_pep_flip_button.setEnabled(True)
        self.iw._rebuild_sel_residue_frame.setEnabled(True)
        chain_id, resnum, resname = (res.chain_id, res.number, res.name)
        self.iw._rebuild_sel_residue_info.setText(str(chain_id) + ' ' + str(resnum) + ' ' + resname)
        
        self._steps_per_sel_res_update = 10
        self._sel_res_update_counter = 0
        if 'update_selected_residue_info' not in self._event_handler.list_event_handlers():
            self._event_handler.add_event_handler('update_selected_residue_info',
                    'atomic changes', self._update_selected_residue_info_live)
                    
                    
    def _disable_rebuild_residue_frame(self):
        if 'update_selected_residue_info' in self._event_handler.list_event_handlers():
            self._event_handler.remove_event_handler('update_selected_residue_info')
        self.iw._rebuild_sel_residue_info.setText('(Select a mobile residue)')
        self.iw._rebuild_sel_res_pep_info.setText('')
        self.iw._rebuild_sel_res_rot_info.setText('')
        self.iw._rebuild_sel_res_rot_target_button.setText('Set target')
        self.iw._rebuild_sel_res_rot_combo_box.clear()
        self.iw._rebuild_sel_res_rot_combo_box.addItem('Pause sim to adjust')
        self.iw._rebuild_sel_residue_frame.setDisabled(True)
        
    def _get_rotamer_list_for_selected_residue(self, *_):
        # stub
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
        self.flip_peptide_omega(res)
    
    def flip_peptide_omega(self, res):
        from math import degrees
        bd = self._mobile_backbone_dihedrals
        omega = bd.omega.by_residue(res)
        if omega is None:
            import warnings
            warnings.warn('Could not find residue or residue has no omega dihedral.')
            return
        oval = degrees(omega.value)
        if oval >= -30 and oval <= 30:
            # cis, flip to trans
            target = 180
        else:
            target = 0
        from . import dihedrals
        self.apply_peptide_bond_restraints(dihedrals.Dihedrals([omega]), target = target, units = 'deg')


    
    ####
    # Validation tab
    ####
    def _show_rama_plot(self, *_):
        self.iw._validate_rama_stub_frame.hide()
        self.iw._validate_rama_main_frame.show()
        if self._rama_plot is None:
            # Create the basic MatPlotLib canvas for the Ramachandran plot
            self._prepare_ramachandran_plot()
        if self._simulation_running and self.track_rama:
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
            
        else:
            sel = model.atoms.filter(model.atoms.selected)
            residues = sel.unique_residues
            phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
            from . import dihedrals
            bd = dihedrals.Backbone_Dihedrals(self.session, phi=phi, psi=psi, omega=omega)
        #phi_vals, psi_vals, omega_vals = bd.all_vals
        #self.rama_validator.load_structure(bd.residues, bd.resnames, phi_vals, psi_vals, omega_vals)
        self._rama_plot.update_scatter(bd, force_update = True)
    
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
    
    def _change_force_field(self):
        ffindex = self.iw._sim_force_field_combo_box.currentIndex()
        ffs = self._available_ffs
        self._sim_main_ff = ffs.main_files[ffindex]
        self._sim_implicit_solvent_ff = ffs.implicit_solvent_files[ffindex]
    
    def _change_water_model(self):
        ffindex = self.iw._sim_water_model_combo_box.currentIndex()
        self._sim_water_ff = self._available_ffs.explicit_water_files[ffindex]
    
    def _change_selected_model(self, *_, model = None):
        self._status('Finding backbone dihedrals. Please be patient.')
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
        if self._selected_model != m and m is not None:
            self._selected_model = m
            self.session.selection.clear()
            self._selected_model.selected = True
            from . import dihedrals
            self.backbone_dihedrals = dihedrals.Backbone_Dihedrals(self.session, m)            
            self._update_chain_list()
            self.triggers.activate_trigger('selected model changed', data=m)
        self._status('')
    
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
        self.b_and_a_padding = self.iw._sim_basic_mobile_b_and_a_spinbox.value()
        
    def _change_soft_shell_cutoff(self, *_):
        self.soft_shell_cutoff = self.iw._sim_basic_mobile_sel_within_spinbox.value()
    
    def _change_soft_shell_fix_backbone(self, *_):
        self.fix_soft_shell_backbone = not self.iw._sim_basic_mobile_sel_backbone_checkbox.checkState()
    
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
        self.hide_surroundings_during_sim = self.iw._sim_hide_surroundings_toggle.checkState() 



    def add_map(self, name, vol, cutoff, coupling_constant, style = None, 
                color = None, contour = None, contour_units = None, mask = True):
        if name in self.master_map_list:
            for key in self.master_map_list:
                print(key)
            raise Exception('Each map must have a unique name!')
        # Check if this model is a unique volumetric map
        if len(vol.models()) !=1 or not hasattr(vol, 'grid_data'):
            raise Exception('vol must be a single volumetric map object')
        
        from .volumetric import IsoldeMap
        new_map = IsoldeMap(self.session, name, vol, cutoff, coupling_constant, style, color, contour, contour_units, mask)
        self.master_map_list[name] = new_map
        self._update_master_map_list_combo_box()

        
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
    
    def start_sim(self):
        try:
            self._start_sim()
        except Exception as e:
            self.triggers.activate_trigger('simulation terminated', 
                                        (self.sim_outcome.FAIL_TO_START, e))
    
    
    def _start_sim(self):
        log = self.session.logger.info
        from time import time
        last_time = time()
        if self._logging:
            self._log('Initialising simulation')
            
        print('Simulation should start now')
        if self._simulation_running:
            print('You already have a simulation running!')
            return
        
        sm = self._sim_modes
        if self.sim_mode in [sm.xtal, sm.em]:
            if not len(self.master_map_list):
                self._no_maps_warning()
                return
                
        self._status('Defining simulation selection...')
        
        # Define final mobile selection
        self._get_final_sim_selection()
        
        # Define "soft shell" of mobile atoms surrounding main selection
        self._soft_shell_atoms = self.get_shell_of_residues(
            self._selected_atoms,
            self.soft_shell_cutoff
            )
        
        
        # Define fixed selection (whole residues with atoms coming within
        # a cutoff of the mobile selection
        total_mobile = self._total_mobile = self._selected_atoms.merge(self._soft_shell_atoms)
        self._hard_shell_atoms = self.get_shell_of_residues(
            total_mobile,
            self.hard_shell_cutoff
            )
        
        sc = self._total_sim_construct = total_mobile.merge(self._hard_shell_atoms)
        #sb = self._total_sim_bonds = sc.inter_bonds
        
        self._total_mobile_indices = sc.indices(total_mobile)
            
        sc.selected = True
        from chimerax.core.atomic import selected_atoms, selected_bonds
        sc = self._total_sim_construct = selected_atoms(self.session)
        sb = self._total_sim_bonds = selected_bonds(self.session)
        import numpy
        surr = self._surroundings = self._selected_model.atoms.filter(numpy.invert(
            self._selected_model.atoms.selected))
        
        
        # Cache all the colors so we can revert at the end of the simulation
        self._original_atom_colors = sc.colors
        self._original_atom_draw_modes = sc.draw_modes
        self._original_bond_radii = sb.radii
        self._original_atom_radii = sc.radii
        self._original_display_state = self._selected_model.atoms.displays
        
        sc.selected = False
        self._hard_shell_atoms.selected = True
        hsb = selected_bonds(self.session)
        hsb.radii = 0.1
        self._hard_shell_atoms.radii = 0.1
        self._hard_shell_atoms.draw_modes = 3
        
        sc.selected = True
        
        # Collect backbone dihedral atoms and prepare the Ramachandran
        # validator
        log('Initial phase took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()
        
        self._status('Organising dihedrals...')
        
        from . import dihedrals
        all_bd = self.backbone_dihedrals
        sim_phi, sim_psi, sim_omega = all_bd.by_residues(total_mobile.unique_residues)
        bd = self._mobile_backbone_dihedrals \
            = dihedrals.Backbone_Dihedrals(self.session, phi = sim_phi, psi = sim_psi, omega = sim_omega)
        log('Preparing backbone dihedrals took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()
        if self.track_rama:
            bd.CAs.draw_modes = 1

            self.omega_validator.load_structure(self._selected_model, bd.omega)
            log('Preparing peptide bond validation took {0:0.4f} seconds'.format(time() - last_time))
            last_time = time()

        from . import sim_interface as si
        sh = self._sim_handler = si.SimHandler(self.session)        
        
        self._status('Preparing restraints...')            
        # Register each dihedral with the CustomTorsionForce, so we can
        # restrain them at will during the 
        for dlist in (bd.phi, bd.psi, bd.omega):
            atoms = dlist.atoms
            indices = numpy.reshape(sc.indices(atoms),[len(dlist),4])
            for i, d in enumerate(dlist):
                d_indices = indices[i]
                d.sim_index = sh.initialize_dihedral_restraint(d, d_indices)
        log('Defining dihedral restraint forces took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()

        
        # If we wish to apply peptide bond restraints, initialise them here
        if self.restrain_peptide_bonds:
            self.apply_peptide_bond_restraints(bd.omega, target = None)
        log('Adding peptide bond restraints took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()

        
        if self._logging:
            self._log('Generating topology')
            
        # Setup pulling force, to be used by mouse and/or haptic tugging
        e = 'k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        per_particle_parameters = ['k','x0','y0','z0']
        per_particle_defaults = [0,0,0,0]
        global_parameters = None
        global_defaults = None
        tug_force = sh.initialize_tugging_force(e, global_parameters, global_defaults, per_particle_parameters)
        
        sh.register_custom_external_force('tug', tug_force, global_parameters,
                per_particle_parameters, per_particle_defaults)
        
        
        self._status('Preparing maps...')
        # Crop down maps and convert to potential fields
        if self.sim_mode in [sm.xtal, sm.em]:
            for mkey in self.master_map_list:
                m = self.master_map_list[mkey]
                vol = m.mask_volume_to_selection(
                    total_mobile, invert = False)
                c3d = sh.continuous3D_from_volume(vol)
                m.set_c3d_function(c3d)
                f = sh.map_potential_force_field(c3d, m.get_coupling_constant())
                m.set_potential_function(f)
                # Register the map with the SimHandler
                sh.register_map(m)
                
                # If required, mask the map visualisation down to the mobile selection
                do_mask = m.get_mask_vis()
                if do_mask:
                    v = m.get_source_map()
                    cutoff = m.get_mask_cutoff()
                    from chimerax.core.commands import sop
                    sop.sop_zone(self.session, v.surface_drawings, near_atoms = total_mobile, range = cutoff)
        
        log('Preparing maps took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()

        self._status('Preparing simulation topology...')
        # Generate topology
        self._topology, self._particle_positions = sh.openmm_topology_and_external_forces(
            sc, sb, tug_hydrogens = False, hydrogens_feel_maps = False,
            logging = self._logging, log = self._log)
        
        log('Preparing topology took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()
        
        self._status('Initialising OpenMM simulation...')
        if self._logging:
            self._log('Generating forcefield')
        
        forcefield_list = [self._sim_main_ff,
                            self._sim_implicit_solvent_ff,
                            self._sim_water_ff]
        
        self._ff = si.define_forcefield(forcefield_list)
        
        if self._logging:
            self._log('Preparing system')
        
        # Define simulation System
        sys = self._system = si.create_openmm_system(self._topology, self._ff)
        
        log('Creating system took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()
        
        # Register extra forces
        for key, f in sh.get_all_custom_external_forces().items():
            sys.addForce(f[0])
        
        sm = self._sim_modes
        if self.sim_mode in [sm.xtal, sm.em]:
            for key,m in self.master_map_list.items():
                sys.addForce(m.get_potential_function())
        
        sys.addForce(sh._dihedral_restraint_force)

        if self._logging:
            self._log('Choosing integrator')
        
        integrator = si.integrator(self._integrator_type,
                                    self.simulation_temperature,
                                    self._friction,
                                    self._integrator_tolerance,
                                    self._sim_time_step)
        
        if self._logging:
            self._log('Setting platform to ' + self.sim_platform)
        platform = si.platform(self.sim_platform)

        if self._logging:
            self._log('Generating simulation')
                                    
        self.sim = si.create_sim(self._topology, self._system, integrator, platform)
        
        log('Finalising sim platform took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()
            
        # Apply fixed atoms to System
        if self._logging:
            self._log('Applying fixed atoms')
        
        fixed = numpy.array(len(sc)*[False],numpy.bool)
        fixed[sc.indices(self._hard_shell_atoms)] = True
        if self.fix_soft_shell_backbone:
            indices = sc.indices(self._soft_shell_atoms)
            names = self._soft_shell_atoms.names
            for i, n in zip(indices, names):
                if n in ['N','C','O','H','H1','H2','H3']:
                    fixed[i] = True        
        for i, f in enumerate(fixed):
            if f:
                self.sim.system.setParticleMass(i, 0)
        if True in fixed:
            self.sim.context.reinitialize()
            
        log('Applying fixed atoms took {0:0.4f} seconds'.format(time() - last_time))
        last_time = time()

        
        # Save the current positions in case of reversion
        self._saved_positions = self._particle_positions
        # Go
        c = self.sim.context
        from simtk import unit
        c.setPositions(self._particle_positions/10) # OpenMM uses nanometers
        c.setVelocitiesToTemperature(self.simulation_temperature)
        
        self._sim_startup = True
        self._sim_startup_counter = 0

        if self._logging:
            self._log('Starting sim')
        
        # Register simulation-specific mouse modes
        from . import mousemodes
        mt = self._mouse_tugger = mousemodes.TugAtomsMode(self.session, total_mobile, self._annotations)
        self._mouse_modes.register_mode(mt.name, mt, 'right', [])
        
        self._event_handler.add_event_handler('do_sim_steps_on_gui_update',
                                              'new frame',
                                              self.do_sim_steps)
        self._simulation_running = True
        self._update_sim_control_button_states()

        if self.track_rama and self._rama_plot is not None:
            self._rama_go_live()
        
        log('Everything else took {0:0.4f} seconds'.format(time() - last_time))
        self._status('Simulation running')
        
    # Get the mobile selection. The method will vary depending on
    # the selection mode chosen
    def _get_final_sim_selection(self):
                        
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
            pad = self.b_and_a_padding
            from chimerax.core.atomic import selected_atoms
            import numpy
            selatoms = selected_atoms(self.session)
            us = selatoms.unique_structures
            if len(us) != 1:
                print(len(us))
                for m in us:
                    print(m.category)
                raise Exception('Selected atoms must all be in the same model!')
            self._selected_model = us[0]    
            selatoms_by_chain = selatoms.by_chain
            selchains = [row[1] for row in selatoms_by_chain]
            allatoms = self._selected_model.atoms
            
            allchains = self._selected_model.chains
            allchainids = list(allchains.chain_ids)
            numchains = len(allchains)
            
            chain_mask = numpy.zeros(numchains,dtype=bool)
            for c in selchains:
                i = allchainids.index(c)
                chain_mask[i] = True
                
            # Throw out the chains containing no selected atoms
            from itertools import compress
            allchains = list(compress(allchains, chain_mask))
            
            for selchain, wholechain in zip(selatoms_by_chain, allchains):
                selatoms = selchain[2]
                sel_resnums_in_chain = selatoms.residues.unique().numbers
                all_residues_in_chain = wholechain.existing_residues
                max_resid_in_chain = all_residues_in_chain.numbers.max()
                min_resid_in_chain = all_residues_in_chain.numbers.min()
                resid_in_range = numpy.zeros(max_resid_in_chain+1,dtype=bool)
                for r in sel_resnums_in_chain:
                    minr = r-pad
                    maxr = r+pad+1
                    if maxr > max_resid_in_chain:
                        maxr = max_resid_in_chain
                    if minr < 0:
                        minr = 0
                    resid_in_range[minr:maxr] = [True] * (maxr-minr)
                all_residues_in_chain.filter(resid_in_range[all_residues_in_chain.numbers]).atoms.selected = True
                    
            
            self._selected_atoms = selected_atoms(self.session)
                       
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
        shell_atoms = unselected_atoms[shell_indices].residues.atoms
        return shell_atoms

        
        
        
        
        
        
    ##############################################################
    # Simulation on-the-fly control functions
    ##############################################################
            
    def pause_sim_toggle(self):
        print('This function should toggle pause/resume of the sim')
        if self._simulation_running:
            if not self._sim_paused:
                self._event_handler.remove_event_handler('do_sim_steps_on_gui_update')
                self._sim_paused = True
                self._status('Simulation paused')
                self.iw._sim_pause_button.setText('Resume')
            else:
                self._event_handler.add_event_handler('do_sim_steps_on_gui_update',
                                      'new frame',
                                      self.do_sim_steps)
                self._sim_paused = False
                self._status('Simulation running')
                self.iw._sim_pause_button.setText('Pause')    
    
    def discard_sim(self):
        print("""This function should stop the simulation and revert to
                 the original coordinates""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        if self._saved_positions is not None:
            self._total_sim_construct.coords = self._saved_positions
        if self.track_rama:
            self.update_omega_check()
        self.triggers.activate_trigger('simulation terminated', (self.sim_outcome.DISCARD, None))
    
    def commit_sim(self):
        print("""This function should stop the simulation and write the
                 coordinates to the target molecule""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        self.triggers.activate_trigger('simulation terminated', (self.sim_outcome.COMMIT, None))
        
    def minimize(self):
        print('Minimisation mode')
        self.simulation_type = 'min'
        self._update_sim_control_button_states()
    
    def equilibrate(self):
        print('Equilibration mode')
        self.simulation_type = 'equil'
        self._update_sim_control_button_states()
    
    def _cleanup_after_sim(self, name, outcome):
        outcomes = self.sim_outcome
        if outcome[0] == outcomes.UNSTABLE:
            print('''Unable to minimise excessive force. Giving up. Try starting 
                a simulation on a larger selection to give the surroundings 
                more room to settle.''')
            '''
            TODO: Pop up a dialog box, giving the user the opportunity to save
            or discard the existing coordinates, and/or start a simulation
            with a wider selection to give the minimiser more room to work.
            '''
            pass
        elif outcome[0] == outcomes.FAIL_TO_START:
            '''
            TODO: Pop up a dialog box giving the details of the exception, then
            clean up.
            '''
            pass
        # Otherwise just clean up            
        self._simulation_running = False
        if 'do_sim_steps_on_gui_update' in self._event_handler.list_event_handlers():
            self._event_handler.remove_event_handler('do_sim_steps_on_gui_update')
        self.sim = None
        self._system = None
        self._update_menu_after_sim()
        self._mouse_modes.remove_mode(self._mouse_tugger.name)
        #mouse_mode_names = self._mouse_modes.get_names()
        #for n in mouse_mode_names:
            #self._mouse_modes.remove_mode(n)
        self._mouse_tugger = None
        for d in self._haptic_devices:
            d.cleanup()
        self._selected_model.atoms.displays = self._original_display_state
        self._surroundings_hidden = False
        self._total_sim_construct.colors = self._original_atom_colors
        self._total_sim_construct.draw_modes = self._original_atom_draw_modes
        self._total_sim_construct.radii = self._original_atom_radii
        self._total_sim_bonds.radii = self._original_bond_radii
        self._total_sim_construct = None
        self._surroundings = None
        if self.track_rama:
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
        self._status('')
    
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
        
        import numpy
        # Get all atom indices in one go because the lookup is expensive
        indices = numpy.reshape(sc.indices(dihedrals.atoms),[len(dihedrals),4])
        
        for i, d in enumerate(dihedrals):
            if target is None:
                dval = d.value
                if dval < cr[1] and dval > cr[0]:
                    t = 0
                else:
                    t = pi
            sh.set_dihedral_restraint(context, sc, d, indices[i], t, k)
                    
    def remove_peptide_bond_restraints(self, dihedrals):
        '''
        Remove restraints on a list of peptide bonds (actually, just set
        their spring constants to zero). Simulation must already be running.
        '''
        import numpy
        sc = self._total_sim_construct
        sh = self._sim_handler
        context = self.sim.context
        # Get all atom indices in one go because the lookup is expensive
        indices = numpy.reshape(sc.indices(dihedrals.atoms),[len(dihedrals),4])
        for i, d in enumerate(dihedrals):
            if d is not None:
                sh.set_dihedral_restraint(context, sc, d, indices[i], 0, 0)        


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
        self._sim_startup_counter
        s_max_count = self._sim_startup_rounds
        if self._sim_startup and self._sim_startup_counter >= s_max_count:
            print('Maximum number of startup minimisation rounds reached. \
                    starting dynamics anyway...')
            startup = self._sim_startup = False
        pos = self._particle_positions
        sc = self._total_sim_construct
        surr = self._surroundings
        
        # Check if we need to hide or show the surroundings
        if self.hide_surroundings_during_sim:
            if not self._surroundings_hidden:
                surr.displays = False
                self._surroundings_hidden = True
        elif self._surroundings_hidden:
            surr.displays = True
            self._surroundings_hidden = False
        
        newpos, max_force, max_index = self._get_positions_and_max_force()
        if startup:
            print('Startup round {} max force: {:0.0f} kJ/mol/nm'
                    .format(self._sim_startup_counter, max_force))
            if abs(max_force - self._last_max_force) < 1.0 and max_force < self._max_allowable_force:
                print('Minimisation converged. Starting dynamics.')
                startup = self._sim_startup = False
            self._last_max_force = max_force
        elif max_force > self._max_allowable_force and not self._sim_is_unstable:
            self._sim_is_unstable = True
            self._oldmode = mode
            if mode == 'equil':
                # revert to the coordinates before this simulation step
                #c.setPositions(pos/10)
                self.simulation_type = 'min'
                return
        if self._sim_is_unstable:
            if self._unstable_min_rounds >= self.max_unstable_rounds:
                # We have a problem we can't fix. Terminate the simulation
                self.triggers.activate_trigger('simulation terminated', (self.sim_outcome.UNSTABLE, None))
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
                self.simulation_type = self._oldmode
                self._unstable_min_rounds = 0
            

        
        # Mouse interaction
        mtug = self._mouse_tugger
        t_force = self._tugging_force
        t_k = self.tug_force_constant
        cur_tug = self._currently_tugging
        tugging, tug_atom, xyz0 = mtug.get_status()
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
                    if a is not None:
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
                        a = self._haptic_tug_atom[i] = picking.pick_closest_to_point(self.session, 
                            pointer_pos, self._total_mobile, 3.0)
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
            s.minimizeEnergy(maxIterations = steps)
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
        self._particle_positions = newpos
        sc.coords = self._particle_positions
        self._last_frame_number = v.frame_number
        if self._logging:
            self._log('Ran ' + str(self.sim_steps_per_update) + ' steps')
        if self.track_rama:
            self.update_ramachandran()
            self.update_omega_check()
        
        
    def update_ramachandran(self):
        self._rama_counter = (self._rama_counter + 1) % self.steps_per_rama_update
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
        if self._rama_counter == 0:
            cis, twisted = ov.find_outliers()
            ov.draw_outliers(cis, twisted)
        ov.update_coords()
        
    def _get_positions_and_max_force (self):
        import numpy
        c = self.sim.context
        from simtk.unit import kilojoule_per_mole, nanometer, angstrom
        state = c.getState(getForces = True, getPositions = True)
        # We only want to consider forces on the mobile atoms to decide
        # if the simulation is unstable
        forces = (state.getForces(asNumpy = True) \
            /(kilojoule_per_mole/nanometer))[self._total_mobile_indices]
        forcesx = forces[:,0]
        forcesy = forces[:,1]
        forcesz = forces[:,2]
        magnitudes =numpy.sqrt(forcesx*forcesx + forcesy*forcesy + forcesz*forcesz)
        max_mag = max(magnitudes)
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
    def _on_close(self):
        self.session.logger.status('Closing ISOLDE and cleaning up')
    
        # Remove all registered event handlers
        eh_keys = list(self._event_handler.list_event_handlers())
        for e in eh_keys:
            self._event_handler.remove_event_handler(e)
        
        # Revert mouse modes
        self._set_chimerax_mouse_mode_panel_enabled(True)
        self._mouse_modes.remove_all_modes()


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
    
    

