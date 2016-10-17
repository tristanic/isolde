
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
        self.session = session
        
        from . import sim_interface
        
        # Dict containing list of all currently loaded atomic models
        self._available_models = {}
        # Selected model on which we are actually going to run a simulation
        self._selected_model = None
        # If we're running by chain(s)
        self._selected_atoms = set()
        self._soft_shell_atoms = set()
        
        # List of forcefields available to the MD package
        self._available_ffs = sim_interface.available_forcefields()
        # Variables holding current forcefield choices
        self._sim_main_ff = None
        self._sim_implicit_solvent_ff = None
        self._sim_water_ff = None 
        # Currently chosen mode for selecting the mobile simulation
        self._sim_selection_mode = None
        
        
        
        
        # Holds the current simulation mode, to be chosen from the GUI
        # drop-down menu or loaded from saved settings.
        self.sim_mode = None
        # Do we have a simulation running right now?
        self._simulation_running = False
        
        # Simulation temperature in Kelvin
        self.simulation_temperature = 100.0
        # Flag to update the temperature of a running simulation
        self._temperature_changed = False
        
        # If a simulation is running, is it paused?
        self._sim_paused = False
        
        # Are we equilibrating or minimising?
        self.simulation_type = 'equil'
        
        
        ####
        # Connect and initialise ISOLDE widget
        ####
        
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
            self.iw._intermediate_frame_example
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
        
        
        from . import eventhandler
        self._event_handler = eventhandler.EventHandler(self.session)
        self._selection_handler = self._event_handler.add_event_handler(
            'update_menu_on_selection', 'selection changed',
            self._selection_changed
            )
        self._selection_changed
        
        self._model_add_handler = self._event_handler.add_event_handler(
            'update_menu_on_model_add', 'add models', 
            self._update_model_list
            )
        
        self._model_add_handler = self._event_handler.add_event_handler(
            'update_menu_on_model_remove', 'remove models', 
            self._update_model_list
            )
        self._update_model_list
        
        
        # Work out menu state based on current ChimeraX session
        self._update_sim_control_button_states()
        self._selection_changed()
        
        
        
    def _populate_menus_and_update_params(self):
        iw = self.iw
        # Clear experience mode placeholders from QT Designer and repopulate
        cb = iw._experience_level_combo_box
        cb.clear()
        for lvl in self._experience_levels:
            cb.addItem(lvl.name, int(lvl))
        
        # Clear simulation mode placeholders from QT Designer and repopulate
        cb = iw._sim_basic_mode_combo_box
        cb.clear()
        for mode in self._sim_modes:
            text = self._human_readable_sim_modes[mode]
            cb.addItem(text, int(mode))
        
        iw._sim_temp_spin_box.setProperty('value', self.simulation_temperature)
    
        # Populate force field combo box with available forcefields
        cb = iw._sim_force_field_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.main_file_descriptions)
        
        # Populate water model combo box with available models
        cb = iw._sim_water_model_combo_box
        cb.clear()
        cb.addItems(self._available_ffs.explicit_water_descriptions)
    
    def _connect_functions(self):
        iw = self.iw
        ####
        # Master switches
        ####
        
        iw._experience_level_combo_box.currentIndexChanged.connect(
            self._change_experience_level_or_sim_mode
            )
        ## Initialise to current level
        #self._change_experience_level_or_sim_mode()
        
        iw._sim_basic_mode_combo_box.currentIndexChanged.connect(
            self._change_experience_level_or_sim_mode
            )            
        # Initialise to selected mode. 
        self._change_experience_level_or_sim_mode()
        
        for button in self._selection_mode_buttons:
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
            
        # Run all connected functions once to initialise
        self._change_force_field()
        self._change_water_model()    
        self._change_selected_model()
        self._change_selected_chains()
        
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
        
    def _update_sim_temperature(self):
        t = self.iw._sim_temp_spin_box.value
        self.simulation_temperature = t
        # So we know to update the temperature in any running simulation
        self._temperature_changed = True
        
    
    ##############################################################
    # Menu control functions to run on key events
    ##############################################################

    def _update_model_list(self, *_):
        self.iw._sim_basic_whole_model_combo_box.clear()
        self.iw._sim_basic_by_chain_model_combo_box.clear()
        models = self.session.models.list()
        model_list = []
        sorted_models = sorted(models, key=lambda m: m.id)
        # Remove non-atomic models from the list
        for i, m in enumerate(sorted_models):
            if not m.atomspec_has_atoms():
                sorted_models.pop(i)
        for model in sorted_models:
            id_str = model.id_string()
            self._available_models[id_str] = model
            model_list.append(id_str)
        self.iw._sim_basic_whole_model_combo_box.addItems(model_list)
        self.iw._sim_basic_by_chain_model_combo_box.addItems(model_list)            

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
        
        i = 0
        for b in self._selection_mode_buttons:
            if b.isChecked():
                break
            i += 1
        
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
        
    
    
    
    
    # Update button states after a simulation has finished
    def _update_menu_after_sim(self):
        self._update_sim_control_button_states()
        
    
    
    ##############################################################
    # Simulation global settings functions
    ##############################################################
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
        lb_sels = lb.selectedItems()
        sel_chain_list = []
        self.session.selection.clear()
        self._selected_atoms = []
        for s in lb_sels:
            sel_chain_list.append(s.text())
        for a in m.atoms:
            if a.chain_id in sel_chain_list:
                a.selected = True
                self._selected_atoms.append(a)

        
        #self.session.selection.clear()
        #self._selected_chains.selected = True
            
    ##############################################################
    # Simulation prep
    ##############################################################
    
    
    def start_sim(self):
        print('Simulation should start now')
        if self._simulation_running:
            print('You already have a simulation running!')
            return
        self._simulation_running = True
        self._update_sim_control_button_states()
        
        # Define final mobile selection
        self._get_final_sim_selection()
        
        # Define "soft shell" of mobile atoms surrounding main selection
        self._get_sim_soft_shell()
        
        
        # Define fixed selection (whole residues with atoms coming within
        # a cutoff of the mobile selection
        
        # Generate topology
        
        # Define simulation system
        
        # Go
    
    def _get_final_sim_selection(self):
        # Get the mobile selection. The method will vary depending on
        # the selection mode chosen in the GUI
        mode = self._sim_selection_mode
        modes = self._sim_selection_modes
        
        if mode == modes.chain or mode == modes.whole_model:
            # Then everything is easy. The selection is already defined
            sel = self._selected_atoms
        elif mode == modes.from_picked_atoms:
            # A bit more complex. Have to work through the model to find
            # the picked atoms (making sure only one model is selected!),
            # then work back and forward from each picked atom to expand
            # the selection by the specified number of residues.            
            m_list = self.session.selection.models()
            for i, m in enumerate(m_list):
                if not hasattr(m, 'num_atoms'):
                    m_list.pop(i)
            if len(m_list) > 1:
                print(len(m_list))
                raise Exception('Selected atoms must all be in the same model!')
            m = m_list[0]
            self._selected_model = m
            pad = self.iw._sim_basic_mobile_b_and_a_spinbox.value()
            for chain in m.chains:
                selected_indices = []
                for r, resid in enumerate(chain.residues):
                    if resid is None:
                        continue
                    for atom in resid.atoms:
                        if atom.selected:
                            for n in range(r-pad, r+pad+1):
                                selected_indices.append(n)
                            # No need to continue with this residue    
                            break
                for s in selected_indices:
                    if s < 0 or s > len(chain.residues):
                        continue
                    sr = chain.residues[s]
                    for atom in sr.atoms:
                        atom.selected = True
            self._selected_atoms = set()
            for atom in m.atoms:
                if atom.selected:
                    self._selected_atoms.add(atom)
                        
        elif mode == modes.custom:
            # relatively simple. Just need to apply the in-built selection
            # text parser
            pass
    
    def _get_sim_soft_shell(self):
        from scipy.spatial.distance import pdist, squareform
        import numpy as np
        allatoms = self._selected_model.atoms
        allcoords = allatoms.coords
        dists = squareform(pdist(allcoords))
        sm = self._selected_model
        
        dist_cutoff = self.iw._sim_basic_mobile_sel_within_spinbox.value()
        close_atoms = set()
        
        for i, atom in enumerate(sm.atoms):
            if atom.selected:
                for j, d in enumerate(dists[i]):
                    if d < dist_cutoff:
                        close_atoms.add(j)
        
        self._soft_shell_atoms = set()
        shell_residues = set()
        for index in close_atoms:
            a = allatoms[index]
            r = a.residue
            if r in shell_residues:
                # This residue's already been done. No need to waste time
                continue
            shell_residues.add(r)
            self._soft_shell_atoms.update(r.atoms)
            r.atoms.selected = True
            
        
        
    ##############################################################
    # Simulation on-the-fly control functions
    ##############################################################
        
        
    
    def pause_sim_toggle(self):
        print('This function should toggle pause/resume of the sim')
        if self._simulation_running:
            if not self._sim_paused:
                self._sim_paused = True
                self.iw._sim_pause_button.setText('Resume')
            else:
                self._sim_paused = False
                self.iw._sim_pause_button.setText('Pause')    
    
    def discard_sim(self):
        print("""This function should stop the simulation and revert to
                 the original coordinates""")
        if not self._simulation_running:
            print('No simulation running!')
            return
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
        self._update_menu_after_sim()
        
    

        


            
    
    
    

