# MolProbity plugin for ChimeraX
# Copyright information here

from PyQt5.QtCore import QObject, pyqtSignal
import chimerax

class MolProbity_GUI:
    '''
    MolProbity is at its core an extensive set of statistical data describing
    what "real" proteins look like, based on a curated non-redundant selection
    of the highest resolution structures in the Protein Data Bank. Secondarily,
    it is a suite of software for comparing your structure to this dataset,
    flagging "unusual" conformations for inspection and, where necessary,
    correction. Always keep in mind that not all outliers are wrong - physically,
    they correspond to strained conformations that would require strong
    supporting interactions with their surroundings to remain stable. Where
    these exist they are usually very clear in the experimental data, and are
    often of biological interest. As a rule of thumb, one can say that if it's
    not strongly supported by the map, your outlier is probably an error.

    This plugin is a port of the core MolProbity functions to run natively in
    ChimeraX, with (wherever possible) optimisations for fast repeated measures
    of the same structure allowing for interactive real-time validation during
    the course of model-building and/or interactive simulation.
    '''
    ###########
    # Environment information
    ###########
    import sys, os
    _root_dir = os.path.dirname(os.path.abspath(__file__))
    _platform = sys.platform

    def __init__(self, gui):
        self.session = gui.session
        self._status = self.session.logger.status

        # Simple class to keep track of ChimeraX trigger handlers we've created
        # to connect events such as model changes to MolProbity methods. The
        # primary role of this is simply to make for easy clean-up when
        # the MolProbity widget is closed.
        from .eventhandler import EventHandler
        self._event_handler = EventHandler(self.session)


        # Selected model on which we're currently focusing
        self._selected_model = None

        # Model object to hold annotations (filled in cis/twisted peptide bonds, etc.)
        from chimerax.core.models import Model
        self._annotations = Model('MolProbity annotations', self.session)
        self.session.models.add([self._annotations])


        # Load in Ramachandran maps
        from . import validation
        # object containing all the Ramachandran contours and lookup functions
        self._status('Preparing Ramachandran contours...')
        self.rama_validator = validation.RamaValidator()
        self._status('')
        # object that handles checking and annotation of peptide bond geometry
        self.omega_validator = validation.OmegaValidator(self._annotations)
        # Generic widget object holding the Ramachandran plot. This is created
        # as an empty widget in Qt5 Designer, and will be linked and filled with
        # a validation.RamaPlot() object once the UI is loaded below.
        self._rama_plot_window = None
        # validation.RamaPlot() Object holding Ramachandran plot information and controls
        self._rama_plot = None
        # Is the Ramachandran plot running live?
        self._update_rama_plot = False
        # dihedrals.Backbone_Dihedrals object holding the protein phi, psi and
        # omega dihedrals for the currently selected model. This gets generated
        # when a model is chosen, and remains in memory until a different model is
        # chosen or MolProbity is closed. If the user wishes to analyse only a
        # subset of residues, a temporary Backbone_Dihedrals object will be
        # created from this one.
        self.backbone_dihedrals = None
        # dihedrals.Backbone_Dihedrals object holding only the backbone dihedrals
        # that are currently selected for analysis
        self._selected_backbone_dihedrals = None

        # Automatically update Ramachandran plot when atoms move?
        self.track_rama = True

        # Start the actual MolProbity widget as a standalone tool
        self.start_gui(gui)

    @property
    def selected_model(self):
        return self._selected_model

    @selected_model.setter
    def selected_model(self, model):
        if not isinstance(model, chimerax.core.atomic.AtomicStructure):
            raise TypeError('Selection must be a single AtomicStructure model!')
        self._change_selected_model(model = model)

    def start_gui(self, gui):
        '''
        Starts MolProbity as a standalone GUI tool. To start as a widget for
        inline use in other tools, use start_widget().
        '''
        self.gui = gui
        self.mw = gui.mw
        self.gui_mode = True

        # Function to remove all event handlers etc. when MolProbity is closed,
        # and return ChimeraX to its standard state.
        self.gui.tool_window.ui_area.destroyed.connect(self._on_close)

        # Any values in the Qt Designer .ui file are placeholders only.
        # Any combo boxes and menus need to be repopulated with their final
        # entries.
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


    def _populate_menus_and_update_params(self):
        mw = self.mw
        # Populate the Ramachandran plot case selector with available
        # cases
        cb = mw._validate_rama_case_combo_box
        cb.clear()
        from . import validation
        # First two keys are N- and C-terminal residues, which we don't plot
        keys = validation.RAMA_CASES[2:]
        for key in reversed(keys):
            cb.addItem(validation.RAMA_CASE_DETAILS[key]['name'], key)

    def _connect_functions(self):
        '''
        Connect PyQt events from the GUI widget to functions.
        '''
        mw = self.mw

        mw._validate_rama_show_button.clicked.connect(
            self._show_rama_plot
            )
        mw._validate_rama_hide_button.clicked.connect(
            self._hide_rama_plot
            )
        mw._validate_rama_case_combo_box.currentIndexChanged.connect(
            self._change_rama_case
            )
        mw._validate_rama_go_button.clicked.connect(
            self._rama_static_plot
            )
        mw._validate_pep_show_button.clicked.connect(
            self._show_peptide_validation_frame
            )
        mw._validate_pep_hide_button.clicked.connect(
            self._hide_peptide_validation_frame
            )
        mw._validate_pep_update_button.clicked.connect(
            self._update_iffy_peptide_lists
            )
        mw._validate_pep_cis_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )
        mw._validate_pep_twisted_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )


    def _prepare_ramachandran_plot(self):
        '''
        Prepare an empty MatPlotLib figure to put the Ramachandran plots in.
        '''
        from . import validation
        mw = self.mw
        container = self._rama_plot_window = mw._validate_rama_plot_layout
        self._rama_plot = validation.RamaPlot(self.session, container, self.rama_validator)



    def _on_close(self):
        self.session.logger.status('Closing ISOLDE and cleaning up')

        # Remove all registered event handlers
        eh_keys = list(self._event_handler.list_event_handlers())
        for e in eh_keys:
            self._event_handler.remove_event_handler(e)
