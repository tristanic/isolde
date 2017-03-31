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

    '''
    Named triggers to simplify handling of changes on key MolProbity events,
    using the same code as ChimeraX's session.triggers. Once a trigger 
    is added to the list below, you can do the following anywhere in the
    code. 
    
    # Adding a handler to run a callback whenever the trigger is fired:
    self._model_changed_handler = self.triggers.add_handler(
        'selected model changed', self._model_changed_callback)
    
    # Setting up the callback function:
    def _model_changed_callback(self, trigger_name, model):
        ... whatever you want to do with the model (e.g. refill a combo
        box with the chain names, rebuild the dihedral topology objects,
        ...)
    
    # Removing the handler to end automatic running of the callback:
    self.triggers.remove_handler(self._model_changed_handler)
    
    # Actually firing the trigger:
    self.triggers.activate_trigger('selected model changed', data = model)
    
    '''    
    from chimerax.core import triggerset
    triggers = triggerset.TriggerSet()
    trigger_names = [
        'selected model changed', # Changed the master model selection
        ]
    for t in trigger_names:
        triggers.add_trigger(t)




    def __init__(self, session, tool_gui = None, widget_container = None):
        '''
        Start the MolProbity plugin, either as a standalone GUI or as
        a widget to be placed within other plugins. 
        Args:
            session: the ChimeraX session
            tool_gui: If launched in standalone mode (i.e. via 
                Tools/General/MolProbity), this argument will be used by
                tool.py. When launched in this way, the widget will be
                wrapped with an outer widget providing a drop-down menu
                for molecule selection and any other features you wish
                to add.
            widget_container: To launch MolProbity as a widget, provide
                a suitable QWidget object to contain it. 
        '''
        self.session = session
        if (tool_gui is not None) == (widget_container is not None):
            errMsg = 'Please provide either the tool_gui or the \
                widget_container argument'
            raise TypeError(errMsg)
        
        if tool_gui:
            self._standalone_mode = True
        else:
            self._standalone_mode = False
        
        self._status = self.session.logger.status

        # Simple class to keep track of ChimeraX trigger handlers we've created
        # to connect events such as model changes to MolProbity methods. The
        # primary role of this is simply to make for easy clean-up when
        # the MolProbity widget is closed.
        from .eventhandler import EventHandler
        self._event_handler = EventHandler(self.session)
        
        # Selected model on which we're currently focusing
        self._selected_model = None
        
        # Callback functions to run whenever the selected model changes.
        # Only place generic functions common to both widget and standalone
        # mode here. Format is:
        # 'descriptive unique key': (function, None)
        # ... where the None is a placeholder for the actual handler when
        # it's generated.
        # Functions that go here should be:
        #   - those involved in pre-analysing the model - for example,
        #     finding, categorising and caching all the backbone dihedrals 
        #     necessary for Ramachandran and Omega analysis. Similar
        #     functions will be needed for rotamers, CaBLAM, suiteness
        #     etc.
        #   - those that change the behaviour of the widget based 
        #     on the model. For example, it might be nice to disable/
        #     hide segments irrelevant to the current model (e.g. 
        #     RNA/DNA analysis for a protein-only structure or vice 
        #     versa).
        self._selected_model_changed_callbacks = {
            'generate backbone dihedrals': [self._backbone_dihedral_cb, None],
        }
        for key, t in self._selected_model_changed_callbacks.items():
            t[1] = self.triggers.add_handler('selected model changed', t[0])
        

        # Model object to hold annotations (filled in cis/twisted peptide bonds, etc.)
        # At the moment this is created at the top-level node, so it's 
        # somewhat divorced from the atomic model. In many ways it would
        # be preferable to group it in with the atomic model, so that 
        # each model can carry around its own annotations.
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
        
        if self._standalone_mode:
            # Start the actual MolProbity widget as a standalone tool
            self.start_gui(tool_gui)
        
        else:
            # Start the widget and place it in the given container
            self.start_widget(widget_container)

    @property
    def selected_model(self):
        return self._selected_model

    @selected_model.setter
    def selected_model(self, model):
        from chimerax.core.atomic import AtomicStructure
        if model is not None and not isinstance(model, AtomicStructure):
            raise TypeError('Selection must be a single AtomicStructure model!')
        self.change_selected_model(model)

    def start_gui(self, gui):
        '''
        Starts MolProbity as a standalone GUI tool. To start as a widget for
        inline use in other tools, use start_widget(container).
        '''
        self.gui = gui
        mm = self.mm = gui.mm
        self.start_widget(self.mm.main_widget_container)
        self.gui_mode = True
        # Dict containing list of all currently loaded atomic models.
        self._available_models = {}


        # Function to remove all event handlers etc. when MolProbity is closed,
        # and return ChimeraX to its standard state.
        self.gui.tool_window.ui_area.destroyed.connect(self._on_close)

        self._event_handler.add_event_handler('update_menu_on_model_add',
                                              'add models',
                                              self._update_model_list_standalone)
        self._event_handler.add_event_handler('update_menu_on_model_remove',
                                              'remove models',
                                              self._update_model_list_standalone)
        self._update_model_list_standalone()

        mm.selected_model_combo_box.currentIndexChanged.connect(
            self._change_selected_model_standalone
            )
        self._change_selected_model_standalone()
            
    def start_widget(self, parent):
        '''
        Start the MolProbity widget proper, and place it in the given
        parent container (e.g. a QWidget, QFrame, QVBoxLayout, ...). 
        Returns a handle for the widget.
        '''
        from PyQt5 import QtWidgets, QtGui
        from . import molprobity_widget
        self.mainwin = QtWidgets.QFrame()
        parent.addWidget(self.mainwin)
        mw = self.mw = molprobity_widget.Ui_Frame()
        mw.setupUi(self.mainwin)
        
       # Frames/widgets that should be hidden at the start
        hidden_at_start = [
            mw._validate_rama_main_frame,
            mw._validate_omega_main_frame,
            ]
        
        for f in hidden_at_start:
            f.hide()
        
        
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
        self._selection_changed()

    def _selection_changed(self, *_):
        '''
        Callback if you want MolProbity to respond in any way if atoms 
        are selected in the main window. For example...
        
        To find all selected atoms (from all currently loaded models):
            from chimerax.core.atomic import selected_atoms
            sel = selected_atoms(self.session)
        
        To find only selected atoms from MolProbity's currently selected
        model:
            m = self._selected_model
            sel = m.atoms.filter(m.atoms.selected)
        
        To get the list of residues which have atoms selected:
            residues = sel.unique_residues
        '''
        pass

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
            self._redraw_rama_plot
            )
        mw._validate_omega_show_button.clicked.connect(
            self._show_peptide_validation_frame
            )
        mw._validate_omega_hide_button.clicked.connect(
            self._hide_peptide_validation_frame
            )
        mw._validate_omega_update_button.clicked.connect(
            self._update_iffy_peptide_lists
            )
        mw._validate_omega_cis_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )
        mw._validate_omega_twisted_list.itemClicked.connect(
            self._show_selected_iffy_peptide
            )




    def _change_selected_model_standalone(self, *_, model = None):
        '''
        Callback that runs if a model is chosen from the "Selected model"
        drop-down menu when MolProbity is run as a standalone tool. 
        Activates the "selected model changed" trigger (and thereby any
        callbacks under its control), and selects the atoms in the model
        to highlight it in the main view window. 
        '''
        if len(self._available_models) == 0:
            return
        mm = self.mm
        mw = self.mw
        if model is not None:
            # Find and select the model in the master combo box, which
            # will automatically call this function again with model = None
            index = mm.selected_model_combo_box.findData(model)
            mm.selected_model_combo_box.setCurrentIndex(index)
            return
        m = mm.selected_model_combo_box.currentData()
        if self._selected_model != m and m is not None:
            self._selected_model = m
            self.session.selection.clear()
            self._selected_model.selected = True
            self.triggers.activate_trigger('selected model changed', data=m)
            self._redraw_rama_plot()
    
    def change_selected_model(self, model):
        '''
        Set the atomic model for MolProbity to work on. Upon selection 
        there may be up to a few seconds delay as MolProbity prepares
        the structure for analysis.
        '''
        if model is None:
            self.clear_model_selection()
            return
            
        if self._standalone_mode:
            self._change_selected_model_standalone(model = model)
        else:
            self._selected_model = model
            self.triggers.activate_trigger('selected model changed', data=model)
            self._redraw_rama_plot()
    
    def clear_model_selection(self):
        self._selected_model = None
        self.triggers.activate_trigger('selected model changed', data=None)
        self._redraw_rama_plot()

    def _update_model_list_standalone(self, *_):
        mm = self.mm
        mw = self.mw
        mm.selected_model_combo_box.clear()
        models = self.session.models.list()
        self._available_models = {}
        atomic_model_list = []
        atomic_model_name_list = []
        sorted_models = sorted(models, key=lambda m: m.id)
        if len(sorted_models) != 0:
            # Find atomic models and sort them into the list
            for i, m in enumerate(sorted_models):
                if m.atomspec_has_atoms():
                    id_str = m.id_string() + ' ' + m.name
                    self._available_models[id_str] = m
                    atomic_model_name_list.append(id_str)
                    atomic_model_list.append(m)
            for l, m in zip(atomic_model_name_list, atomic_model_list):
                mm.selected_model_combo_box.addItem(l, m)
        if not len(self._available_models):
            self.selected_model = None




            
    def _backbone_dihedral_cb(self, trigger_name, model):
        self.generate_backbone_dihedrals(model)
    
    def generate_backbone_dihedrals(self, model):
        '''
        Finds and organises all protein phi, psi and omega dihedrals into
        an object optimised for rapid measurement, scoring and plotting
        of the current conformation.
        '''
        if model is None:
            self.backbone_dihedrals = None
            return
        
        self._status('Finding backbone dihedrals. Please be patient.')
        from . import dihedrals
        self.backbone_dihedrals = dihedrals.Backbone_Dihedrals(self.session, model)
        self._status('')
        return self.backbone_dihedrals            

    def _prepare_ramachandran_plot(self):
        '''
        Prepare an empty MatPlotLib figure to put the Ramachandran plots in.
        '''
        from . import validation
        mw = self.mw
        container = self._rama_plot_window = mw._validate_rama_plot_layout
        self._rama_plot = validation.RamaPlot(self.session, container, self.rama_validator)

    def _show_rama_plot(self, *_):
        mw = self.mw
        mw._validate_rama_stub_frame.hide()
        mw._validate_rama_main_frame.show()
        if self._rama_plot is None:
            # Create the basic MatPlotLib canvas for the Ramachandran plot
            self._prepare_ramachandran_plot()
        self._redraw_rama_plot()

    def _hide_rama_plot(self, *_):
        self.mw._validate_rama_main_frame.hide()
        self.mw._validate_rama_stub_frame.show()

    def _redraw_rama_plot(self, *_):
        '''
        Updates the Ramachandran scatter plot if it's visible.
        '''
        if not self.mw._validate_rama_main_frame.isVisible():
            return
        model = self._selected_model
        if model is not None:
            whole_model = bool(self.mw._validate_rama_sel_combo_box.currentIndex())
            if whole_model:
                self._rama_plot.update_scatter(self.backbone_dihedrals, 
                                                force_update = True)
            else:
                sel = model.atoms.filter(model.atoms.selected)
                residues = sel.unique_residues
                if len(residues):
                    phi, psi, omega = self.backbone_dihedrals.by_residues(residues)
                    from . import dihedrals
                    bd = dihedrals.Backbone_Dihedrals(self.session, 
                                        phi=phi, psi=psi, omega=omega)
                    self._rama_plot.update_scatter(bd, force_update = True)
                else:
                    self._rama_plot.update_scatter(force_update = True)
        else:
            self._rama_plot.update_scatter(force_update = True)
    
    def _change_rama_case(self, *_):
        case_key = self.mw._validate_rama_case_combo_box.currentData()
        self._rama_plot.change_case(case_key)
    
    def _show_peptide_validation_frame(self, *_):
        self.mw._validate_omega_stub_frame.hide()
        self.mw._validate_omega_main_frame.show()
    
    def _hide_peptide_validation_frame(self, *_):
        self.mw._validate_omega_main_frame.hide()
        self.mw._validate_omega_stub_frame.show()    

    def _update_iffy_peptide_lists(self, *_):
        '''
        Finds all cis and twisted peptide bonds in the current model, 
        provides lists in the GUI (with callbacks to display the 
        offending residues when clicked), and draws a pseudo-planar
        trapezoid filling the problematic omega dihedrals.
        '''
        ov = self.omega_validator
        model = self._selected_model
        mw = self.mw
        clist = mw._validate_omega_cis_list
        tlist = mw._validate_omega_twisted_list
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

    

    def _on_close(self):
        self.session.logger.status('Closing ISOLDE and cleaning up')

        # Remove all registered event handlers
        eh_keys = list(self._event_handler.list_event_handlers())
        for e in eh_keys:
            self._event_handler.remove_event_handler(e)
