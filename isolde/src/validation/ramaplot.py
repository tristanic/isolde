# @Author: Tristan Croll
# @Date:   03-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy

class RamaPlot:
    def __init__(self, session, rama_mgr, container):
        import numpy
        self.session = session
        mgr = self._rama_mgr = rama_mgr
        cenum = self._case_enum = rama_mgr.Rama_Case
        self.container = container
        self.current_case = None

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qt5agg import (
            FigureCanvasQTAgg as FigureCanvas,
            NavigationToolbar2QT as NavigationToolbar
            )

        self._current_residues = None
        self._all_current_ramas = None
        self._case_ramas = None

        fig = self.figure = Figure()

        axes = self.axes = fig.add_subplot(111, aspect='equal')

        # Scatter plot needs to always have at least one point, so we'll
        # define a point off the scale for when there's nothing to plot
        self.default_coords = numpy.array([200,200])
        self.default_logscores = numpy.ones(1)
        self.scatter = None

        self.base_scatter_size = 10
        canvas = self.canvas = FigureCanvas(fig)
        self.resize_cid = self.canvas.mpl_connect('resize_event', self.on_resize)
        self.on_pick_cid = self.canvas.mpl_connect('pick_event', self.on_pick)
        container.addWidget(canvas)

        self.contours = {}
        self.change_case(cenum.GENERAL)
        self.on_resize()

    def _format_axes(self):
        axes = self.axes
        fig = self.figure
        axes.set_xticks([-120,-60,0,60,120])
        axes.set_yticks([-120,-60,0,60,120])
        #axes.set_xticklabels(axes.get_xticklabels(), rotation=60)
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
        pcolor_plot = self.axes.pcolor(*grid, logvalues, cmap = 'Greys')
        for coll in contour_plot.collections:
           coll.remove()
        pcolor_plot.remove()
        self.contours[key] = (contour_plot, pcolor_plot)


    def on_resize(self, *_):
        axes = self.axes
        f = self.figure
        c = self.canvas
        self._format_axes()
        if self.scatter.figure is not None:
            self.scatter.remove()
        f.patch.set_facecolor('0.5')
        c.draw()
        self.background = c.copy_from_bbox(axes.bbox)
        axes.add_artist(self.scatter)
        c.draw()
        self.update_scatter()

    def on_pick(self, event):
        ind = event.ind[0]
        picked_rama = self._case_ramas[ind]
        from .. import view
        view.focus_on_selection(self.session, self.session.main_view, picked_rama.residue.atoms)

    def change_case(self, case_key):
        import numpy
        mgr = self._rama_mgr
        self.current_case = case_key
        self.axes.clear()
        try:
            contourplots = self.contours[case_key]
        except KeyError:
            self.cache_contour_plots(case_key)
            contourplots = self.contours[case_key]
        for coll in contourplots[0].collections:
            self.axes.add_artist(coll)
        self.axes.add_artist(contourplots[1])
        from math import log
        contours = mgr.RAMA_CASE_DETAILS[case_key]['cutoffs']
        self.P_limits = [0, -log(contours[0])]
        scatter = self.scatter = self.axes.scatter(
            (200),(200), picker = 2.0, s=self.base_scatter_size,
            edgecolors='black', linewidths = 0.5)
        scatter.set_cmap('RdYlGn_r')
        self.set_target_residues(self._current_residues)
        self.on_resize()

    def set_target_residues(self, residues):
        self._current_residues = residues
        case = self.current_case
        mgr = self._rama_mgr
        cenum = self._case_enum
        if residues is not None:
            ramas = mgr.get_ramas(residues)
            ramas = ramas[ramas.valids]
            cases = ramas.cases
            self._all_current_ramas = ramas
            self._case_ramas = ramas[cases == case]
        else:
            self._all_current_ramas = None
            self._case_ramas = None

    def update_scatter(self, *_, residues=None):
        if not self.container.parent().isVisible():
            return
        if residues is not None:
            self.set_target_residues(residues)
        import numpy
        case = self.current_case
        #~ rv = self.validator
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
