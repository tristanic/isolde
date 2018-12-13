# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll


import numpy

class RamaPlot:
    WHOLE_MODEL = 0
    SELECTED_ONLY = 1
    MOBILE_ONLY = 2

    mode_dict = {
        WHOLE_MODEL: 'All residues',
        SELECTED_ONLY: 'Selection',
        MOBILE_ONLY: 'Mobile atoms only'
    }
    def __init__(self, session, isolde, container, mode_menu, case_menu,
            restrict_button):
        import numpy
        self.session = session
        from chimerax.isolde import session_extensions as sx
        mgr = self._rama_mgr = sx.get_ramachandran_mgr(session)
        self.isolde = isolde
        cenum = self._case_enum = mgr.Rama_Case
        self.container = container
        self.current_case = None
        self._selection_mode = self.WHOLE_MODEL

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qt5agg import (
            FigureCanvasQTAgg as FigureCanvas,
            NavigationToolbar2QT as NavigationToolbar
            )

        self._current_model = None
        self._current_residues = None
        self._all_current_ramas = None
        self._case_ramas = None

        fig = self.figure = Figure()

        axes = self.axes = fig.add_subplot(111, aspect='equal')


        canvas = self.canvas = FigureCanvas(fig)
        self.resize_cid = self.canvas.mpl_connect('resize_event', self.on_resize)
        self.on_pick_cid = self.canvas.mpl_connect('pick_event', self.on_pick)
        container.addWidget(canvas)

        self.contours = {}

        # Scatter plot needs to always have at least one point, so we'll
        # define a point off the scale for when there's nothing to plot
        self.default_coords = numpy.array([200,200])
        base_size = self.base_scatter_size = 10
        self.default_logscores = numpy.ones(1)
        scatter = self.scatter = self.axes.scatter(
            (200),(200), picker = 2.0, s=base_size,
            edgecolors='black', linewidths = 0.5)
        scatter.set_cmap('RdYlGn_r')

        self.change_case(cenum.GENERAL)

        self._isolde_switch_model_handler = isolde.triggers.add_handler(
            'selected model changed',
            self._isolde_switch_model_cb
        )

        self._model_changes_handler = None
        self.current_model = isolde.selected_model
        #self.on_resize()

        self._sel_restrict_button = restrict_button
        restrict_button.clicked.connect(
            self._restrict_sel_cb
        )

        self._mode_menu = mode_menu
        self._populate_mode_menu(mode_menu)
        mode_menu.currentIndexChanged.connect(
            self._mode_change_cb
        )

        self._case_menu = case_menu
        self._populate_case_menu(case_menu)
        case_menu.currentIndexChanged.connect(self._case_change_cb)

        self.isolde.triggers.add_handler(
            'simulation started', self._sim_start_cb
        )

        self.isolde.triggers.add_handler(
            'simulation terminated', self._sim_end_cb
        )


    def _prepare_tooltip(self):
        ax = self.axes
        annot = self._hover_annotation = ax.annotate(
            '', xy=(0,0),
            xytext=(20,20),
            textcoords='offset points',
            bbox=dict(boxstyle='round', fc='w'),
            arrowprops=dict(arrowstyle='->'))
        annot.set_visible(False)

        def _hover(event):
            sp = self.scatter
            annot = self._hover_annotation
            if event.inaxes == ax:
                cont, ind = sp.contains(event)
                if cont:
                    indices = ind['ind']
                    text = []
                    ramas = self._case_ramas[indices]
                    phipsi = numpy.degrees(ramas.phipsis)
                    x = phipsi[:,0].mean()
                    y = phipsi[:,1].mean()
                    for rama in ramas:
                        res = rama.residue
                        text.append('{} {}{}'.format(res.name, res.chain_id, res.number))
                    annot.set_text('\n'.join(text))
                    annot.xy = x, y
                    if y > 0:
                        annot.xyann=(20,-20)
                    else:
                        annot.xyann=(20,20)
                    annot.set_visible(True)
                else:
                    annot.set_visible(False)
                self.canvas.draw()

        self.canvas.mpl_connect('motion_notify_event', _hover)




    def _sim_start_cb(self, *_):
        sc = self.isolde.sim_manager.sim_construct
        self.set_target_residues(sc.mobile_residues)
        self.selection_mode = self.MOBILE_ONLY

    def _sim_end_cb(self, *_):
        self._mode_change_cb()

    def _restrict_sel_cb(self, *_):
        m = self.current_model
        if m is None:
            self.clear()
            return
        selres = m.atoms[m.atoms.selected].unique_residues
        self.set_target_residues(selres)
        self.update_scatter()

    def _populate_mode_menu(self, menu):
        menu.clear()
        shown_modes = [
            self.WHOLE_MODEL,
            self.SELECTED_ONLY
        ]
        for mode in shown_modes:
            menu.addItem(self.mode_dict[mode])

    def _populate_case_menu(self, menu):
        menu.clear()
        rm = self._rama_mgr
        keys = list(rm.Rama_Case)[1:]
        for key in reversed(keys):
            menu.addItem(rm.RAMA_CASE_DETAILS[key]['name'], key)

    def _mode_change_cb(self, *_):
        mode = self._mode_menu.currentIndex()
        self.selection_mode = mode

    def _case_change_cb(self, *_):
        case_key = self._case_menu.currentData()
        self.change_case(case_key)

    @property
    def current_model(self):
        return self._current_model

    @current_model.setter
    def current_model(self, model):
        if self._model_changes_handler is not None:
            self.current_model.triggers.remove_handler(
                self._model_changes_handler
            )
            self._model_changes_handler = None
        if model:
            self._model_changes_handler = model.triggers.add_handler(
                'changes', self._model_changed_cb
            )
        self._current_model=model
        if self.selection_mode == self.WHOLE_MODEL:
            residues = model.residues
        else:
            residues = model.atoms[model.atoms.selected].unique_residues
        self.set_target_residues(residues)
        self.update_scatter()

    @property
    def selection_mode(self):
        return self._selection_mode

    @selection_mode.setter
    def selection_mode(self, mode):
        if mode not in self.mode_dict.keys():
            raise TypeError('Unrecognised mode!')
        self._selection_mode = mode
        if mode == self.WHOLE_MODEL:
            self.set_target_residues(self.current_model.residues)

        sflag = (mode == self.SELECTED_ONLY)
        self._sel_restrict_button.setEnabled(sflag)

        mflag = (mode != self.MOBILE_ONLY)
        self._mode_menu.setEnabled(mflag)

        self.update_scatter()


    def _isolde_switch_model_cb(self, trigger_name, model):
        self.current_model = model

    def _model_changed_cb(self, trigger_name, changes):
        changes = changes[1]
        update_needed = False
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        added = len(changes.created_atoms())
        deleted = changes.num_deleted_atoms()
        if added or deleted:
            update_needed = True
        if update_needed:
            self.update_scatter()

    @property
    def parent(self):
        return self.container.parentWidget()

    @property
    def visible(self):
        return self.parent.isVisible()

    def _format_axes(self):
        axes = self.axes
        fig = self.figure
        axes.set_xticks([-120,-60,0,60,120])
        axes.set_yticks([-120,-60,0,60,120])
        axes.set_xticklabels([])
        axes.set_yticklabels([])
        axes.set_xlim(-180,180)
        axes.set_ylim(-180,180)
        axes.minorticks_on()
        axes.tick_params(direction='in')
        axes.autoscale(enable=False)
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    def cache_contour_plots(self, key):
        import numpy
        mgr = self._rama_mgr
        grid = numpy.degrees(mgr.interpolator_axes(int(key)))
        values = numpy.rot90(numpy.fliplr(mgr.interpolator_values(key))).astype(float)
        contours = mgr.RAMA_CASE_DETAILS[key]['cutoffs']
        logvalues = numpy.log(values)
        contour_plot = self.axes.contour(*grid, values, contours)
        pcolor_plot = self.axes.pcolormesh(*grid, logvalues, cmap = 'Greys')
        for coll in contour_plot.collections:
           coll.remove()
        pcolor_plot.remove()
        self.contours[key] = (contour_plot, pcolor_plot)


    def on_resize(self, *_):
        axes = self.axes
        f = self.figure
        c = self.canvas
        self._format_axes()
        self.scatter.set_offsets(self.default_coords)
        f.patch.set_facecolor('0.5')
        c.draw()
        self.background = c.copy_from_bbox(axes.bbox)
        self.update_scatter()

    def on_pick(self, event):
        ind = event.ind[0]
        picked_rama = self._case_ramas[ind]
        from .. import view
        view.focus_on_selection(self.session, picked_rama.residue.atoms)

    def change_case(self, case_key):
        import numpy
        mgr = self._rama_mgr
        self.current_case = case_key
        ax = self.axes
        ax.clear()
        try:
            contourplots = self.contours[case_key]
        except KeyError:
            self.cache_contour_plots(case_key)
            contourplots = self.contours[case_key]
        for coll in contourplots[0].collections:
            ax.add_artist(coll)
        ax.add_artist(contourplots[1])
        ax.add_artist(self.scatter)
        from math import log
        contours = mgr.RAMA_CASE_DETAILS[case_key]['cutoffs']
        self.P_limits = [0, -log(contours[0])]
        self.set_target_residues(self._current_residues)
        self.on_resize()
        self._prepare_tooltip()

    def clear(self):
        self.set_target_residues(None)
        self.update_scatter()

    def set_target_residues(self, residues):
        self._current_residues = residues
        case = self.current_case
        mgr = self._rama_mgr
        cenum = self._case_enum
        if residues is not None and len(residues):
            ramas = mgr.get_ramas(residues)
            ramas = ramas[ramas.valids]
            cases = ramas.cases
            self._all_current_ramas = ramas
            self._case_ramas = ramas[cases == case]
        else:
            self._all_current_ramas = None
            self._case_ramas = None

    def update_scatter(self, *_, residues=None):
        if not self.visible:
            return
        if residues is not None:
            self.set_target_residues(residues)
        import numpy
        case = self.current_case
        mgr = self._rama_mgr
        cenum = self._case_enum
        r = self._case_ramas
        if r is None or len(r) == 0:
            phipsi = self.default_coords
            logscores = self.default_logscores
        else:
            phipsi = numpy.degrees(r.phipsis)
            logscores = numpy.log(r.scores)

        c = self.canvas
        s = self.scatter
        axes = self.axes
        c.restore_region(self.background)
        s.set_offsets(phipsi)
        s.set_clim(self.P_limits)

        scales = (-logscores+1)*self.base_scatter_size
        s.set_sizes(scales)
        s.set_array(-logscores)
        axes.draw_artist(s)
        c.blit(axes.bbox)
