# @Author: Tristan Croll <tic20>
# @Date:   10-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 09-Dec-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

import os, sys
import numpy
from math import inf, degrees, radians, pi
from enum import IntEnum

from Qt.QtWidgets import QFileDialog

from openmm import unit


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

from Qt.QtWidgets import QMessageBox

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
        'center_on_sel_when_masking':           (defaults.CENTER_ON_SEL_WHEN_MASKING, None),
        'rounds_per_map_remask':                (defaults.ROUNDS_PER_MAP_REMASK, None),
            # Radius of mask when isolating a selection
        'map_mask_radius':                      (defaults.STANDARD_MAP_MASK_RADIUS, None),
        'spotlight_radius':                     (defaults.SPOTLIGHT_RADIUS, None),
        'map_shannon_rate':                     (defaults.MAP_SHANNON_RATE, None),
            # Width of the green outline surrounding a selection
        'selection_outline_width':              (defaults.SELECTION_OUTLINE_WIDTH, None),

    }
_root_dir = os.path.dirname(os.path.abspath(__file__))
class Isolde():
    '''
    The master ISOLDE class, primarily designed to be used by the ISOLDE GUI as
    the front-end interface for starting, controlling and analysing interactive
    simulations. Should only be used as a session-level singleton. Once initialised
    will be accessible as :py:attr:`session.isolde`.
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
        default        = 0
        advanced       = 1
        developer      = 2

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
    SELECTED_MODEL_CHANGED = 'selected model changed'
    SIMULATION_STARTED = 'simulation started'
    SIMULATION_TERMINATED = 'simulation terminated'
    SIMULATION_PAUSED = 'simulation paused'
    SIMULATION_RESUMED = 'simulation resumed'
    UNPARAMETERISED_RESIDUE = 'unparameterised residue'
    ISOLDE_CLOSED = 'isolde closed'
    GUI_STARTED = 'gui started'
    trigger_names = (
        SELECTED_MODEL_CHANGED, # Changed the master model selection
        SIMULATION_STARTED,
        SIMULATION_TERMINATED,
        SIMULATION_PAUSED,
        SIMULATION_RESUMED,
        UNPARAMETERISED_RESIDUE,
        ISOLDE_CLOSED,
        GUI_STARTED
        )

    def __init__(self, session, gui=None):
        '''
        Initialises the ISOLDE object and adds it to the ChimeraX session as
        ``session.isolde``.

        Args:
            * gui (default=None):
                - Used by the ISOLDE BundleAPI to prepare the ISOLDE GUI
                  interface (i.e. when running from Tools/General/ISOLDE in the
                  ChimeraX GUI menu).
        '''
        self.session = session
        self.triggers = triggerset.TriggerSet()
        for t in self.trigger_names:
            self.triggers.add_trigger(t)
        self._gui = gui

        self._isolde_events = EventHandler(self)
        self._logging = False
        self._log = Logger('isolde.log')
        self.checkpoint_disabled_reasons = {}

        self.params = IsoldeParams()
        '''
        A :class:`IsoldeParams` instance containing basic parameters controlling
        how the interactive simulation selection is defined and displayed.
        '''

        sp = self.sim_params = SimParams()
        '''
        Parameters controlling the initialisation and behaviour of the
        simulation.
        '''

        from .visualisation import VisualisationStateMgr
        self._vis_mgr = VisualisationStateMgr(self.session, self)
        sp.triggers.add_handler(sp.PARAMETER_CHANGED, self._sim_param_changed_cb)

        from .openmm.forcefields import ForcefieldMgr
        ffmgr = self._ff_mgr = ForcefieldMgr(self.session)

        self._status = self.session.logger.status

        eh = self._event_handler = EventHandler(self.session)
        from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
        eh.add_event_handler('Reinitialize maps on model add', ADD_MODELS, self._model_add_cb)
        eh.add_event_handler('Update selected model if current model deleted', REMOVE_MODELS, self._model_remove_cb)
        # Available pre-defined colors
        from chimerax.core import colors
        import copy

        self._available_models = {}

        self._available_colors = copy.copy(colors.BuiltinColors)
        # Remove duplicates
        for key in list(self._available_colors):
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

        # Placeholder for the openmm_interface.SimHandler object used to perform
        # simulation setup and management
        self._sim_manager = None

        # self.initialize_haptics()



        ####
        # Internal trigger handlers
        ####

        self._menu_prepared = False

        self._ui_panels = []

        if session.ui.is_gui:
            self._prepare_environment()
            from .menu import prepare_isolde_menu
            prepare_isolde_menu(self.session)

        self.selected_model = self.find_next_available_model()

        session.isolde = self
        if session.ui.is_gui:
            self._register_tug_modes()
            session._isolde_tb.isolde_started()
        ffmgr.background_load_ff(sp.forcefield)

        self._event_handler.add_event_handler('ramaplot', 'add tool instance', self._register_ramaplot_callbacks)
        # In case RamaPlot is already open
        self._register_ramaplot_callbacks('', session.tools.list())

    def _prepare_environment(self):
        session = self.session
        from chimerax.std_commands import cofr, camera
        cofr.cofr(session, 'centerOfView', show_pivot=True)
        camera.camera(session, 'ortho')
        self._mouse_modes.register_all_isolde_modes()

        # Increase the visibility of the selection outline
        from chimerax.core.commands import run
        run(session, 'set selectionWidth {}'.format(self.params.selection_outline_width))

    @property
    def gui(self):
        return self._gui
    
    @gui.setter
    def gui(self, gui):
        self._gui = gui
    
    @property
    def gui_mode(self):
        return self.session.ui.is_gui

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
        from chimerax.atomic import Residues
        return Residues([r for r in m.residues if getattr(r, 'isolde_ignore', False)])

    @ignored_residues.setter
    def ignored_residues(self, residues):
        m = self.selected_model
        if m is None:
            return
        for r in m.residues:
            r.isolde_ignore=False
        if residues is None:
            return
        for r in residues:
            r.isolde_ignore=True

    @property
    def sim_handler(self):
        '''
        Returns the :class:`isolde.openmm.SimHandler` instance controlling
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
        Returns the :class:`isolde.openmm.SimManager` instance providing
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
        m = self._selected_model
        if m is not None and not hasattr(m, 'session'):
            # model was deleted
            self.change_selected_model(None)
        return self._selected_model

    @selected_model.setter
    def selected_model(self, model):
        from chimerax.atomic import AtomicStructure
        if model is not None and not isinstance(model, AtomicStructure):
            raise TypeError('Selection must be a single AtomicStructure model!')
        self.change_selected_model(model)

    @property
    def selected_atoms(self):
        '''
        Returns the set of selected atoms in ISOLDE's currently selected model.
        If ISOLDE has no model currently selected, returns an empty :class:`Atoms` object.
        '''
        m = self.selected_model
        if m is None:
            from chimerax.atomic import Atoms
            return Atoms()
        return m.atoms[m.atoms.selecteds]

    ###################################################################
    # GUI related functions
    ###################################################################

    def _register_tug_modes(self):
        from .tugging import TugSingleAtomMode, TugResidueMode, TugSelectionMode
        for mode in (TugSingleAtomMode, TugResidueMode, TugSelectionMode):
            if mode not in self.session.ui.mouse_modes.modes:
                self.session.ui.mouse_modes.add_mode(mode(self.session))


    def _sim_param_changed_cb(self, _, data):
        update_methods = {
            'mouse_tug_spring_constant': self._update_mouse_tug_spring_constant,
            'temperature': self._update_sim_temperature,
            'trajectory_smoothing': self._update_smoothing_state,
            'smoothing_alpha': self._update_smoothing_amount,
        }
        param, value = data
        f = update_methods.get(param, None)
        if f is not None:
            f(value)
    
    def _update_smoothing_state(self, flag):
        if self.simulation_running:
            self.sim_handler.smoothing = flag

    def _update_smoothing_amount(self, alpha):
        if self.simulation_running:
            self.sim_handler.smoothing_alpha = alpha

    def _update_mouse_tug_spring_constant(self, k):
        from .tugging import TugAtomsMode
        for mode in self.session.ui.mouse_modes.modes:
            if isinstance(mode, TugAtomsMode):
                mode.spring_constant = k


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

    def _update_sim_temperature(self, t):
        if self.simulation_running:
            self.sim_handler.temperature = t

    ##############################################################
    # Menu control functions to run on key events
    ##############################################################

    @property
    def available_models(self):
        from chimerax.atomic import AtomicStructure
        models = {}
        for m in self.session.models.list():
            if type(m)==AtomicStructure:
                models[f'{m.id_string}. {m.name}'] = m
        return models

    def _add_cif_template_files_gui(self, *_):
        files = self._choose_cif_template_files(self)
        if files is not None and len(files):
            try:
                from .atomic.template_utils import load_cif_templates
                load_cif_templates(self.session, files)
                self.session.logger.info('You will now be able to add these residues to '
                    'your model with "isolde add ligand {ID}". To be able to simulate '
                    'them, you will need to provide matching ffXML MD parameterisations.')
            except:
                raise

    def _choose_cif_template_files(self, *_):
        options = QFileDialog.Options()
        caption = 'Choose one or more CIF files'
        filetypes = 'CIF template dictionaries (*.cif)'
        dlg = QFileDialog(caption=caption)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilter(filetypes)
        dlg.setFileMode(QFileDialog.ExistingFiles)
        import os
        dlg.setDirectory(os.getcwd())
        if dlg.exec():
            return dlg.selectedFiles()


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
    
    @property
    def selected_model_has_mdff_enabled(self):
        '''
        Returns True if the selected model has at least one MDFF manager present 
        and enabled for MDFF fitting, otherwise False. Read only
        '''
        if not self.selected_model_has_maps:
            return False
        from chimerax.clipper import get_map_mgr
        m = self.selected_model
        mgr = get_map_mgr(m)
        from . import session_extensions as sx
        for v in mgr.all_maps:
            mdff = sx.get_mdff_mgr(m, v)
            if mdff is not None and mdff.enabled:
                return True
        return False

    def _model_add_cb(self, *_):
        m = self.selected_model
        if m is not None:
            self._initialize_maps(m)            

    def _model_remove_cb(self, *_):
        m = self._selected_model
        if m is not None and m.was_deleted:
            self.selected_model = None

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


        from chimerax.clipper.maps import XmapHandler_Live, XmapHandler_Static
        from .session_extensions import get_mdff_mgr
        mdff_mgrs = []
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
            if isinstance(v, XmapHandler_Live):
                # Only XmapHandler_Live objects specifically created to exclude the
                # free set are allowed to act as MDFF potentials
                continue
            new_mgr = False
            mgr = get_mdff_mgr(model, v)
            if mgr is None:
                new_mgr = True
                mgr = get_mdff_mgr(model, v, create=True)
            if isinstance(v, XmapHandler_Static) and new_mgr:
                # Since we don't know the provenance of static maps, it would be
                # a bad idea to simply enable them. But we want people to be
                # able to use them if they really want to. So, we create the
                # manager but leave it disabled by default so the user has to
                # explicitly enable it.
                mgr.enabled = False
            if v.is_difference_map:
                # Similar reasoning for difference maps: there may be occasion to 
                # use them, but this decision should be left up to the user.
                mgr.enabled = False
            mdff_mgrs.append(mgr)
        # recalc_cb.setVisible(live_maps)
        if len(mdff_mgrs):
            return True
        return False

    def _set_live_xmap_recalc(self, state):
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(self.selected_model)
        if mgr is not None:
            for xs in mgr.xmapsets:
                xs.live_update = state


    ####
    # Right mouse button modes
    ####

    def _set_right_mouse_mode_tug_atom(self, *_):
        from chimerax.core.commands import run
        run(self.session, 'ui mousemode right "isolde tug atom"', log=False)


    def _set_right_mouse_mode_tug_residue(self, *_):
        from chimerax.core.commands import run
        run(self.session, 'ui mousemode right "isolde tug residue"', log=False)

    def _set_right_mouse_mode_tug_selection(self, *_):
        from chimerax.core.commands import run
        run(self.session, 'ui mousemode right "isolde tug selection"', log=False)


    ####
    # Validation tab
    ####

    def _register_ramaplot_callbacks(self, trigger_name, tools):
        # Need to delay this to the next frame
        def register_next_frame(*_, tools=tools):
            from chimerax.isolde.validation.ramaplot.tool import Rama_ToolUI
            for tool in tools:
                if isinstance(tool, Rama_ToolUI):
                    break
            else:
                return
            rplot = tool.tool_window
            if rplot.current_model is None:
                rplot.current_model = self.selected_model
            def model_changed_cb(*_, rp=rplot):
                rp.current_model = self.selected_model
            rplot.add_callback(self.triggers, self.SELECTED_MODEL_CHANGED, model_changed_cb)
            # sim start behaviour is managed by the Ramachandran Plot panel on the Validate tab
            # def sim_start_cb(*_, rp=rplot):
            #     if rp.current_model == self.selected_model:
            #         rp.restrict_to_selection(self.sim_manager.sim_construct.mobile_residues, display_text='Mobile residues')
            # rplot.add_callback(self.triggers, self.SIMULATION_STARTED, sim_start_cb)
            def sim_end_cb(*_, rp=rplot):
                if rp.current_model == self.selected_model:
                    rp.display_all_residues()
            rplot.add_callback(self.triggers, self.SIMULATION_TERMINATED, sim_end_cb)
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('new frame', register_next_frame)


    ##############################################################
    # Simulation global settings functions
    ##############################################################

    def find_next_available_model(self):
        from chimerax.clipper import get_symmetry_handler
        from chimerax.atomic import AtomicStructure
        for m in sorted(self.session.models.list(type=AtomicStructure), key=lambda mm: mm.id):
            if get_symmetry_handler(m) is not None:
                return m
        return None

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
        if self.session.ui.is_gui:
            self._change_selected_model(self, model=model, force=True)
        else:
            self._selected_model = model
            self.session.selection.clear()
            if model is not None:
                from .citation import add_isolde_citation
                add_isolde_citation(model)
                # self._selected_model.selected = True
                self._initialize_maps(model)
            self.triggers.activate_trigger(self.SELECTED_MODEL_CHANGED, model)

    def _change_selected_model(self, *_, model=None, force=False):
        m = model
        if self.simulation_running:
            return
        if not hasattr(self, '_model_changes_handler'):
            self._model_changes_handler = None
        session = self.session

        sm = self._selected_model

        if model is None:
            self._selected_model = None
            if sm is not None:
                self.triggers.activate_trigger('selected model changed', data=None)
            return

        with session.triggers.block_trigger('remove models'), session.triggers.block_trigger('add models'):
            if not getattr(m, 'isolde_initialized', False):
                atoms_with_alt_locs = m.atoms[m.atoms.num_alt_locs>0]
                if len(atoms_with_alt_locs):
                    from .dialog import choice_warning
                    result = choice_warning(f'This model contains {len(atoms_with_alt_locs)} atoms with alternate '
                        'conformers. ISOLDE cannot currently see these, but they will be carried through to the '
                        'output model. In most cases it is best to remove them. Would you like to do so now?')
                    if result:
                        m.delete_alt_locs()
                        atoms_with_alt_locs.occupancies = 1
                        self.session.logger.info(f'Removed all altlocs in #{m.id_string} and reset associated occupancies to 1.')
                from .atomic.util import correct_pseudosymmetric_sidechain_atoms
                correct_pseudosymmetric_sidechain_atoms(session, m.residues)
                m.isolde_initialized = True
            from chimerax.clipper.symmetry import get_symmetry_handler
            sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
            from .citation import add_isolde_citation
            add_isolde_citation(m)
            self._selected_model = m
            self.session.selection.clear()
            m.selected = True

            # Load/create validation managers
            from . import session_extensions as sx
            sx.get_rota_annotator(m)
            sx.get_rama_annotator(m)
            self._initialize_maps(m)
            self.triggers.activate_trigger('selected model changed', data=m)


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



    def start_sim(self):
        '''
        Start an interactive simulation based around the currently-selected
        atoms.
        '''
        if self.simulation_running:
            from chimerax.core.errors import UserError
            raise UserError('Simulation already running!')
        self.session.logger.status(
            'Initialising simulation. Please be patient...',
            color='red')
        from .openmm.openmm_interface import SimManager
        sm = self.selected_model
        main_sel = self._last_main_sel = sm.atoms[sm.atoms.selected]
        if not self.selected_model_has_mdff_enabled and getattr(self, '_warn_on_no_maps', True):
            from .dialog import choice_warning
            warn_str = ('You are about to start a simulation with no map fitting forces enabled! '
                'While this is supported by ISOLDE, please be aware that the simulation environment '
                'is not generally suitable for equilibrium simulation and the results should be '
                'interpreted with care. You should strongly consider reinforcing core structures '
                'first with the "isolde restrain distances" command.')
            if self.gui_mode:
                warn_str += '\nAre you sure you want to continue?'
                choice, dont_ask_again = choice_warning(warn_str, allow_dont_ask_again=True)
                self._warn_on_no_maps = not dont_ask_again
                if not choice:
                    return
            else:
                self.session.logger.warning(warn_str)
        try:
            sm = self._sim_manager = SimManager(self, self.selected_model, main_sel,
                self.params, self.sim_params, excluded_residues = self.ignored_residues)
        except ValueError as e:
            err_text = str(e)
            if err_text.startswith("No template found") or \
               err_text.startswith("User-supplied template"):
                # These errors are already handled
                self._sim_end_cb()
                return
            raise
        except RuntimeError as e:
            self._sim_end_cb()
            err_text = str(e)
            if err_text.startswith("Unparameterised"):
                self.triggers.activate_trigger(self.UNPARAMETERISED_RESIDUE, None)
                return
            raise
        except:
            self._sim_end_cb()
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
        m = self.selected_model
        m.atoms.selected = False
        m.bonds.selected = False
        residue.atoms.selected = True
        residue.atoms.intra_bonds.selected = True
        from .view import focus_on_selection
        focus_on_selection(self.session, residue.atoms)
        residue.session.logger.warning('No template found for residue {}{} ({})'.format(
            residue.chain_id, residue.number, residue.name
        ))
        from .dialog import failed_template_warning
        choice = failed_template_warning(residue)
        if choice == 'addh':
            print('Adding hydrogens')
            from chimerax.atomic import AtomicStructures
            from chimerax.addh import cmd
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
        self.session.logger.status('')
        self._set_right_mouse_mode_tug_atom()
        self.triggers.activate_trigger(self.SIMULATION_STARTED, None)
        sh = self.sim_handler
        sh.triggers.add_handler('sim terminated', self._sim_end_cb)
        sh.triggers.add_handler('sim paused', self._sim_pause_cb)
        sh.triggers.add_handler('sim resumed', self._sim_resume_cb)
        self.session.logger.info('ISOLDE: started sim')

    def _sim_end_cb(self, *_):
        for d in self._haptic_devices:
            d.cleanup()
        from chimerax.mouse_modes import TranslateMouseMode
        self.session.ui.mouse_modes.bind_mouse_mode('right', [], TranslateMouseMode(self.session))
        self.triggers.activate_trigger('simulation terminated', None)
        from chimerax.core.commands import run
        run(self.session, f'clipper spot #{self.selected_model.id_string}', log=False)
        from chimerax.clipper import get_map_mgr
        m = self.selected_model
        mmgr = get_map_mgr(m)
        if mmgr is not None:
            for xmapset in mmgr.xmapsets:
                if hasattr(xmapset, 'live_xmap_mgr') and xmapset.live_xmap_mgr is not None:
                    self.session.logger.info('Updating bulk solvent parameters...')
                    xmapset.live_xmap_mgr.bulk_solvent_optimization_needed()
                    xmapset.recalc_needed()
        from .atomic.util import correct_pseudosymmetric_sidechain_atoms
        if self.sim_manager is not None:
            correct_pseudosymmetric_sidechain_atoms(self.session, self.sim_manager.sim_construct.mobile_residues)
        self._sim_manager = None
        self.session.logger.info('ISOLDE: stopped sim')


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

    def pause(self):
        '''
        If simulation is currently running, pause it. Otherwise, do nothing.
        '''
        if self.simulation_running and not self.sim_paused:
            self.pause_sim_toggle()
    
    def resume(self):
        '''
        If an active simulation is currently paused, resume it. Otherwise, do nothing.
        '''
        if self.simulation_running and self.sim_paused:
            self.pause_sim_toggle()


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
        self.triggers.activate_trigger('simulation paused', None)

    def _sim_resume_cb(self, *_):
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
        self.sim_manager.stop_sim(revert=revert_to)

    def commit_sim(self):
        '''
        Stop the simulation and keep the current coordinates.
        '''
        if not self.simulation_running:
            print('No simulation running!')
            return
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
        tugm = sx.get_tuggable_atoms_mgr(self.selected_model, allow_hydrogens=self.sim_params.tug_hydrogens)
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
        tugm = sx.get_tuggable_atoms_mgr(self.selected_model, allow_hydrogens=self.sim_params.tug_hydrogens)
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

    def _cant_checkpoint_error(self):
        if not self.simulation_running:
            self.session.logger.warning('Checkpointing is only available when a simulation is running!')
        else:
            err_str = 'Checkpointing is currently disabled by the '\
                +'following scripts and will be re-enabled when they '\
                +'terminate: \n{}'.format(
                    '\n'.join([r for r in self.checkpoint_disabled_reasons.values()]))
            self.session.logger.warning(err_str)


    def checkpoint(self, *_):
        '''
        Save the current state of the simulation (coordinates and all
        restraints).
        '''
        if self.can_checkpoint:
            self.sim_manager.checkpoint()
        else:
            self._cant_checkpoint_error()

    def revert_to_checkpoint(self, *_):
        '''
        Return the simulation to the last saved checkpoint. Note: a checkpoint
        is automatically saved at the beginning of each simulation.
        '''
        if self.can_checkpoint:
            self.sim_manager.revert_to_checkpoint()
        else:
            self._cant_checkpoint_error()


    @property
    def smoothing(self):
        '''
        Turn smoothing of the simulation trajectory visualisation on or off.
        '''
        return self.sim_params.trajectory_smoothing

    @smoothing.setter
    def smoothing(self, flag):
        self.sim_params.trajectory_smoothing = flag

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
        if alpha < defaults.SMOOTHING_ALPHA_MIN:
            alpha = defaults.SMOOTHING_ALPHA_MIN
        elif alpha > defaults.SMOOTHING_ALPHA_MAX:
            alpha = defaults.SMOOTHING_ALPHA_MAX
        self.sim_params.smoothing_alpha = alpha


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
        omega = pdr_m.add_restraint_by_residue_and_name(res, 'omega')
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
        self.session.logger.status('Closing ISOLDE and cleaning up')

        for p in self._ui_panels:
            p.remove_trigger_handlers()

        if self.simulation_running:
            self.sim_manager.stop_sim()

        self._disable_rebuild_residue_frame()

        # Remove all registered event handlers
        self._event_handler.remove_all_handlers()
        self._isolde_events.remove_all_handlers()

        sm = self.selected_model
        if sm is not None:
            if self._model_changes_handler is not None:
                sm.triggers.remove_handler(self._model_changes_handler)
                

        # Revert mouse modes
        # self._set_chimerax_mouse_mode_panel_enabled(True)
        self._mouse_modes.remove_all_modes()
        self.triggers.activate_trigger('isolde closed', None)
        
        # Remove the ISOLDE object from the session
        delattr(self.session, 'isolde')




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
        fname = os.path.join(self._root_dir, 'docs', 'user', 'isolde.html')
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
        from chimerax.open_command.cmd import provider_open
        model = provider_open(self.session, [model_file])[0]
        from chimerax.std_commands import color
        color.color(self.session, model, color='bychain', target='ac')
        color.color(self.session, model, color='byhetero', target='a')
        from chimerax.clipper import symmetry
        sym_handler = symmetry.SymmetryManager(model,
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
        from chimerax.clipper.util import exclude_nonpolar_hydrogens
        model.atoms[exclude_nonpolar_hydrogens(model.atoms)].displays = True
        from . import view
        view.focus_on_selection(self.session, model.atoms)

    def load_crystal_demo(self):
        m = load_crystal_demo(self.session)
        self.selected_model = m

def load_crystal_demo(session):
    '''
    Load a small protein model with crystallographic maps to explore.
    '''
    from chimerax.open_command.cmd import provider_open
    data_dir = os.path.join(_root_dir, 'demo_data', '3io0')
    before_struct = provider_open(session, [os.path.join(data_dir, 'before.pdb')])[0]
    from chimerax.clipper import symmetry
    sym_handler = symmetry.get_symmetry_handler(before_struct, create=True,
        auto_add_to_session=True)
    xmapset = sym_handler.map_mgr.add_xmapset_from_mtz(os.path.join(data_dir, '3io0-sf.mtz'))
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
    volumecommand.volume(session, [sharp_map], **styleargs)
    styleargs = v.map_style_settings[v.map_styles.solid_t40]
    volumecommand.volume(session, [diff_map], **styleargs)
    # standard_map.set_parameters(surface_levels = (2.5*standard_map.sigma,))
    # sharp_map.set_parameters(surface_levels = (3.0*sharp_map.sigma,))
    diff_map.set_parameters(surface_levels = (-4.0*diff_map.sigma, 4.0*diff_map.sigma))
    from chimerax.std_commands import set
    from chimerax.core.colors import Color
    set.set(session, bg_color=Color([255,255,255,255]))

    if hasattr(session, 'isolde'):
        session.isolde.selected_model = before_struct
    from chimerax.clipper.util import exclude_nonpolar_hydrogens
    before_struct.atoms[exclude_nonpolar_hydrogens(before_struct.atoms)].displays = True
    from . import view
    view.focus_on_selection(session, before_struct.atoms)
    return before_struct

def load_cryo_em_demo(session, model_only=True):
    '''
    Load a high-resolution cryo-EM model and map.
    '''
    from chimerax.open_command.cmd import provider_open
    data_dir = os.path.join(_root_dir, 'demo_data', '6out')
    m = provider_open(session, [os.path.join(data_dir, '6out.pdb')])[0]
    if not model_only:
        mmap = provider_open(session, ['20205'], from_database='emdb')[0]
        from chimerax.core.commands import run
        run(session, f'clipper assoc #{mmap.id_string} to #{m.id_string}', log=False)
    if hasattr(session, 'isolde'):
        session.isolde.selected_model = m


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
