
# ISOLDE: Interactive Structure Optimisation by Local Direct Exploration
# Copyright: 2016
# Author:    Tristan Croll
#            Cambridge Institute for Medical Research
#            University of Cambridge



class Isolde():
    

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

    
    def __init__(self, session):
        self._logging = False
        self._log = Logger('isolde.log')

        self.session = session

        from .eventhandler import EventHandler
        self._event_handler = EventHandler(self.session)
        
        initialize_openmm()
        from . import sim_interface
 
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
        self._b_and_a_padding = 5
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
        # Construct containing all atoms that will actually be simulated
        self._total_sim_construct = None
        # List of all bonds in _total_sim_construct
        self._total_sim_bonds = None
        
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
        # Settings for OpenMM
        ####
        
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
        self.min_steps_per_update = 50
        # If using the VariableLangevinIntegrator, we define a tolerance
        self._integrator_tolerance = 0.0001
        # ... otherwise, we simply set the time per step
        self._sim_time_step = 1.0*unit.femtoseconds
        # Type of integrator to use. Should give the choice in the expert level
        # of the menu. Variable is more stable, but simulated time per gui update
        # is harder to determine
        self._integrator_type = 'variable'
        # Constraints (e.g. rigid bonds) need their own tolerance
        self._constraint_tolerance = 0.001
        # Friction term for coupling to heat bath
        self._friction = 1.0/unit.picoseconds
        # Limit on the net force on a single atom to detect instability and
        # force a minimisation
        self._max_allowable_force = 20000.0 # kJ mol-1 nm-1
        # Flag for unstable simulation
        self._sim_is_unstable = False
        
        
        
        
        # Holds the current simulation mode, to be chosen from the GUI
        # drop-down menu or loaded from saved settings.
        self.sim_mode = None
        # Do we have a simulation running right now?
        self._simulation_running = False
        # If running, is the simulation in startup mode?
        self._sim_startup = True
        # Maximum number of rounds of minimisation to run on startup
        self._sim_startup_rounds = 50
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
    
    ###################################################################
    # GUI related functions
    ###################################################################
    
    def start_gui(self):
        ####
        # Connect and initialise ISOLDE widget
        ####
        
        from . import isolde_gui
        self.gui = isolde_gui.IsoldeGui(self.session)
        self.iw = self.gui.iw
        self.gui.mainwin.show()
        
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



    def _connect_functions(self):
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
        iw._sim_basic_whole_model_combo_box.currentIndexChanged.connect(
            self._change_selected_model
            )
        iw._sim_basic_by_chain_model_combo_box.currentIndexChanged.connect(
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
         # We want to start with the EM map chooser hidden
        self._hide_em_map_chooser()
        
                
        
        
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
        
    
        
    def _update_sim_temperature(self):
        t = self.iw._sim_temp_spin_box.value()
        self.simulation_temperature = t
        # So we know to update the temperature in any running simulation
        self._temperature_changed = True
        
    
    ##############################################################
    # Menu control functions to run on key events
    ##############################################################

    def _update_model_list(self, *_):
        self.iw._sim_basic_whole_model_combo_box.clear()
        self.iw._sim_basic_by_chain_model_combo_box.clear()
        self.iw._em_map_model_combo_box.clear()
        models = self.session.models.list()
        atomic_model_list = []
        volume_model_list = []
        sorted_models = sorted(models, key=lambda m: m.id)
        if len(sorted_models) != 0:
            # Find atomic and volumetric models and sort them into the
            # appropriate lists
            for i, m in enumerate(sorted_models):
                if m.atomspec_has_atoms():
                    id_str = m.id_string()
                    self._available_models[id_str] = m
                    atomic_model_list.append(id_str)
                elif hasattr(m, 'grid_data'):
                    id_str = m.id_string()
                    self._available_volumes[id_str] = m
                    volume_model_list.append(id_str)
                else:
                    # This is a model type we don't currently handle. Ignore.
                    continue
        self.iw._sim_basic_whole_model_combo_box.addItems(atomic_model_list)
        self.iw._sim_basic_by_chain_model_combo_box.addItems(atomic_model_list)
        self.iw._em_map_model_combo_box.addItems(volume_model_list)

    def _update_chain_list(self):
        m = self._selected_model
        chains = m.chains.chain_ids
        lb = iw = self.iw._sim_basic_mobile_chains_list_box
        lb.clear()
        lb.addItems(chains)
        

    def _selection_changed(self, *_):
        if self._simulation_running:
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
        iw._sim_basic_mobile_whole_model_frame.hide()
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
            self._change_selected_model()
            self._change_selected_chains()
        elif i == 2:
            self._sim_selection_mode = self._sim_selection_modes.whole_model
            iw._sim_basic_mobile_whole_model_frame.show()
            self._change_selected_model()
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
        elif len(seltext):
            self._add_new_map = False
            current_map = self.master_map_list[seltext]
            name, vol, cutoff, coupling = current_map.get_map_parameters()
            iw._em_map_name_field.setText(name)
            iw._em_map_model_combo_box.setCurrentText(vol.id_string())
            iw._em_map_cutoff_spin_box.setValue(cutoff)
            iw._em_map_coupling_spin_box.setValue(coupling)
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
    
    def _change_selected_model(self):
        if len(self._available_models) == 0:
            return
        sm = self._sim_selection_mode.name
        iw = self.iw
        if sm == 'whole_model':
            choice = iw._sim_basic_whole_model_combo_box.currentText()
            if choice == '':
                return
            self._selected_model = self._available_models[choice]
            self.session.selection.clear()
            self._selected_model.selected = True
            self._selected_atoms = self._selected_model.atoms
        elif sm == 'chain':
            choice = iw._sim_basic_by_chain_model_combo_box.currentText()
            if choice == '':
                return
            self._selected_model = self._available_models[choice]
            self.session.selection.clear()
            self._selected_model.selected = True
            self._update_chain_list()
    
    def _change_selected_chains(self,*_):
        if len(self._available_models) == 0:
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
        self._b_and_a_padding = self.iw._sim_basic_mobile_b_and_a_spinbox.value()
        
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
        if self._add_new_map:
            self.add_map(name, model, cutoff, coupling_constant)
        else:
            m = self.master_map_list[name]
            m.change_map_parameters(model, cutoff, coupling_constant)

    def _remove_em_map_from_gui(self, *_):
        name = self.iw._em_map_name_field.text()
        self.remove_map(name)


    def add_map(self, name, vol, cutoff, coupling_constant):
        if name in self.master_map_list:
            raise Exception('Each map must have a unique name!')
        # Check if this model is a unique volumetric map
        if len(vol.models()) !=1 or not hasattr(vol, 'grid_data'):
            raise Exception('vol must be a single volumetric map object')
        
        from .volumetric import IsoldeMap
        new_map = IsoldeMap(self.session, name, vol, cutoff, coupling_constant)
        self.master_map_list[name] = new_map
        self._update_master_map_list_combo_box()
        
    def remove_map(self, name):
        result = self.master_map_list.pop(name, 'Not present')
        if result == 'Not present':
            print(name + ' is an unrecognised key.')
        self._update_master_map_list_combo_box()
    

    ##############################################################
    # Simulation prep
    ##############################################################
    
    
    def start_sim(self):
        if self._logging:
            self._log('Initialising simulation')
            
        print('Simulation should start now')
        if self._simulation_running:
            print('You already have a simulation running!')
            return
        
        sm = self._sim_modes
        if self.sim_mode in [sm.xtal, sm.em]:
            if not len(self.master_map_list):
                errstring = 'You have selected ' + \
                self._human_readable_sim_modes[self.sim_mode] + \
                ' but have not selected any maps. Please choose at least one map.'
                raise Exception(errstring)
                
        self._simulation_running = True
        self._update_sim_control_button_states()
        
        
        # Define final mobile selection
        self._get_final_sim_selection()
        
        # Define "soft shell" of mobile atoms surrounding main selection
        self._soft_shell_atoms = self.get_shell_of_residues(
            self._selected_atoms,
            self.soft_shell_cutoff
            )
        
        
        # Define fixed selection (whole residues with atoms coming within
        # a cutoff of the mobile selection
        total_mobile = self._selected_atoms.merge(self._soft_shell_atoms)
        self._hard_shell_atoms = self.get_shell_of_residues(
            total_mobile,
            self.hard_shell_cutoff
            )
        
        sc = total_mobile.merge(self._hard_shell_atoms)
        sc.selected = True
        from chimerax.core.atomic import selected_atoms, selected_bonds
        sc = self._total_sim_construct = selected_atoms(self.session)
        sb = self._total_sim_bonds = selected_bonds(self.session)

        from . import sim_interface as si
                
        self._sim_map_potential_functions = []
        # Crop down maps and convert to potential fields
        if self.sim_mode in [sm.xtal, sm.em]:
            for m in self.master_map_list:
                vol, mincoor, maxcoor = m.mask_volume_to_selection(
                    total_mobile, invert = True)
                c3d = si.continuous3D_from_volume(vol, mincoor, maxcoor)
                self._sim_map_potential_functions.append(c3d)
                    
        
        
        
        
        if self._logging:
            self._log('Generating topology')
        
        # Generate topology
        self._topology, self._particle_positions = si.openmm_topology_and_coordinates(sc, sb, logging = self._logging, log = self._log)

        if self._logging:
            self._log('Generating forcefield')
        
        forcefield_list = [self._sim_main_ff,
                            self._sim_implicit_solvent_ff,
                            self._sim_water_ff]
        
        self._ff = si.define_forcefield(forcefield_list)
        
        if self._logging:
            self._log('Preparing system')
        
        # Define simulation System
        self._system = si.create_openmm_system(self._topology, self._ff)
        
        # Apply fixed atoms to System
        if self._logging:
            self._log('Applying fixed atoms')

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
        
        
        fixed = len(sc)*[False]
        for i, atom in enumerate(sc):
            if atom in self._hard_shell_atoms:
                fixed[i] = True
                continue
            if self.fix_soft_shell_backbone:
                if atom in self._soft_shell_atoms:
                    if atom.name in ['N','C','O','H','H1','H2','H3']:
                        fixed[i] = True
        for i in range(len(fixed)):
            if fixed[i]:
                self.sim.system.setParticleMass(i, 0)
        if True in fixed:
            self.sim.context.reinitialize()
            
 

        
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
        
        self._event_handler.add_event_handler('do_sim_steps_on_gui_update',
                                              'new frame',
                                              self.do_sim_steps)
        
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
            pad = self._b_and_a_padding
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
                self.iw._sim_pause_button.setText('Resume')
            else:
                self._event_handler.add_event_handler('do_sim_steps_on_gui_update',
                                      'new frame',
                                      self.do_sim_steps)
                self._sim_paused = False
                self.iw._sim_pause_button.setText('Pause')    
    
    def discard_sim(self):
        print("""This function should stop the simulation and revert to
                 the original coordinates""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        if self._saved_positions is not None:
            self._total_sim_construct.coords = self._saved_positions
        self._cleanup_after_sim()
    
    def commit_sim(self):
        print("""This function should stop the simulation and write the
                 coordinates to the target molecule""")
        if not self._simulation_running:
            print('No simulation running!')
            return
        # Write coords back to target here
        self._cleanup_after_sim()
        
    def minimize(self):
        print('Minimisation mode')
        self.simulation_type = 'min'
        self._update_sim_control_button_states()
    
    def equilibrate(self):
        print('Equilibration mode')
        self.simulation_type = 'equil'
        self._update_sim_control_button_states()
    
    def _cleanup_after_sim(self):            
        self._simulation_running = False
        if 'do_sim_steps_on_gui_update' in self._event_handler.list_event_handlers():
            self._event_handler.remove_event_handler('do_sim_steps_on_gui_update')
        self.sim = None
        self._system = None
        self._update_menu_after_sim()
        
    
    #############################################
    # Main simulation functions to be run once per GUI update
    #############################################
    
    def do_sim_steps(self,*_):
        if self._logging:
            self._log('Running ' + str(self.sim_steps_per_update) + ' steps')

        v = self.session.main_view
#        if v.frame_number == self._last_frame_number:
#            return # Make sure we draw a frame before doing another MD calculation
            
        s = self.sim
        c = s.context
        integrator = c.getIntegrator()
        steps = self.sim_steps_per_update
        minsteps = self.min_steps_per_update
        mode = self.simulation_type
        startup = self._sim_startup
        s_count = self._sim_startup_counter
        s_max_count = self._sim_startup_rounds
        pos = self._particle_positions
        sc = self._total_sim_construct


        if self._temperature_changed:
            integrator.setTemperature(self.simulation_temperature)
            c.setVelocitiesToTemperature(self.simulation_temperature)
            self._temperature_changed = False

        
        if startup and s_count:
            start_step = s.currentStep
            s.minimizeEnergy(maxIterations = steps)
            end_step = s.currentStep
            if end_step - start_step < steps:
                # minimisation has converged. We can continue on
                startup = False
            else:
                s_count += 1
        elif self._sim_is_unstable:
            s.minimizeEnergy(maxIterations = steps)
        elif mode == 'min':
            s.minimizeEnergy(maxIterations = minsteps)
        elif mode == 'equil':
            s.step(steps)
        else:
            raise Exception('Unrecognised simulation mode!')
        
        newpos, max_force = self._get_positions_and_max_force()
        if max_force > self._max_allowable_force:
            self._sim_is_unstable = True
            self._oldmode = mode
            if mode == 'equil':
                # revert to the coordinates before this simulation step
                c.setPositions(pos/10)
            self.simulation_type = 'min'
            return
        elif self._sim_is_unstable:
            if max_force < self._max_allowable_force / 2:
                # We're back to stability. We can go back to equilibrating
                self._sim_is_unstable = False
                self.simulation_type = self._oldmode
        if self._logging:
            self._log('Sim 4')
        
        from simtk import unit
        self._particle_positions = newpos
        sc.coords = self._particle_positions
        self._last_frame_number = v.frame_number
        if self._logging:
            self._log('Ran ' + str(self.sim_steps_per_update) + ' steps')
        
            
        
    def _get_positions_and_max_force (self):
        import numpy
        c = self.sim.context
        from simtk.unit import kilojoule_per_mole, nanometer, angstrom
        state = c.getState(getForces = True, getPositions = True)
        forces = state.getForces(asNumpy = True)/(kilojoule_per_mole/nanometer)
        forcesx = forces[:,0]
        forcesy = forces[:,1]
        forcesz = forces[:,2]
        magnitudes =numpy.sqrt(forcesx*forcesx + forcesy*forcesy + forcesz*forcesz)
        pos = state.getPositions(asNumpy = True)/angstrom
        return pos, max(magnitudes)

_openmm_initialized = False
def initialize_openmm():
    # On linux need to set environment variable to find plugins.
    # Without this it gives an error saying there is no "CPU" platform.
    global _openmm_initialized
    if not _openmm_initialized:
        _openmm_initialized = True
        from sys import platform
        if platform == 'linux':
            from os import environ, path
            from chimerax import app_lib_dir
            environ['OPENMM_PLUGIN_DIR'] = path.join(app_lib_dir, 'plugins')

   
class Logger:
    def __init__(self, filename = None):
        self.filename = filename
        self._log_file = None
    def __call__(self, message, close = False):
        if self.filename is None:
            return	# No logging
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
    
    

