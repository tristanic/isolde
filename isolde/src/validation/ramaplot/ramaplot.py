# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 08-May-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll


import numpy

from Qt.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSizePolicy
from Qt.QtCore import Qt

class RamaPlot(QWidget):
    def __init__(self, session, manager, rama_case):
        super().__init__(parent=manager.ui_area)
        self._debug=False
        self.session = session
        self.manager = manager
        from chimerax.isolde import session_extensions as sx
        rmgr = self._rama_mgr = sx.get_ramachandran_mgr(session)
        cenum = self._case_enum = rmgr.RamaCase
        self.container = manager.ui_area
        self.case = rama_case
        self.name=rmgr.RAMA_CASE_DETAILS[rama_case]['name']

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qt5agg import (
            FigureCanvasQTAgg as FigureCanvas,
            NavigationToolbar2QT as NavigationToolbar
            )

        self._current_residues = None
        self._all_current_ramas = None
        self._case_ramas = None
        self._displayed_ramas = None

        self.setStyleSheet('color: white; background-color: black;')
        self.setAutoFillBackground(True)

        layout = self.main_layout = QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)

        t = self.title = QLabel(self.name, parent=self)
        layout.addWidget(t)
        t.setAlignment(Qt.AlignHCenter)
        

        fig = self.figure = Figure()

        axes = self.axes = fig.add_subplot(111, aspect='equal')
        axes.autoscale(enable=False)

        canvas = self.canvas = FigureCanvas(fig)

        self.resize_cid = canvas.mpl_connect('resize_event', self.on_resize)
        self.on_pick_cid = canvas.mpl_connect('pick_event', self.on_pick)
        layout.addWidget(canvas)


        # Scatter plot needs to always have at least one point, so we'll
        # define a point off the scale for when there's nothing to plot
        self.default_coords = numpy.array([200,200])
        base_size = self.base_scatter_size
        self.default_logscores = numpy.ones(1)
        scatter = self.scatter = self.axes.scatter(
            (200),(200), picker = 2.0, s=base_size,
            edgecolors='black', linewidths = 0.5)
        scatter.set_cmap('RdYlBu_r')


        self._prepare_contours()

        self._prepare_tooltip()

        self.on_resize()
        self._start_tooltip()


        self._model_changes_handler = None

    @property
    def base_scatter_size(self):
        return self.manager.base_scatter_size

    @property
    def restricted(self):
        from .tool import RamaMainWin
        R = RamaMainWin.Restrictions
        return bool(self.display_mode&R.RESTRICTED)
    
    @property
    def hide_favored(self):
        from .tool import RamaMainWin
        R = RamaMainWin.Restrictions
        return bool(self.display_mode&(R.DISFAVORED_ONLY|R.OUTLIERS_ONLY))
    
    @property
    def outliers_only(self):
        from .tool import RamaMainWin
        R = RamaMainWin.Restrictions
        return bool(self.display_mode&R.OUTLIERS_ONLY)


    @property
    def display_mode(self):
        return self.manager.display_mode
    
    @property
    def current_model(self):
        return self.manager.current_model

    def _selection_changed_cb(self, *_):
        if not self.visible or not self.current_model:
            return
        residues = self.current_model.residues
        sel_res = residues[residues.selected]
        ramas = self._rama_mgr.get_ramas(sel_res)
        cenum = self._case_enum
        if not len(ramas):
            self.update_scatter()
            return
        import numpy
        unique_cases = numpy.unique(ramas.cases)
        unique_cases = unique_cases[unique_cases!=0]
        case = cenum.GENERAL
        if len(unique_cases) == 1:
            case = cenum(unique_cases[0])
        if case != self.case:
            cm = self._case_menu
            cm.setCurrentIndex(cm.findData(case))
        else:
            self.update_scatter()

    def _display_mode_changed_cb(self, _, mode):
        from .tool import RamaMainWin


    def block_tooltip(self):
        self._tooltip_blocked = True
    
    def unblock_tooltip(self):
        self._tooltip_blocked = False
    



    def _prepare_tooltip(self):
        ax = self.axes
        annot = self._hover_annotation = ax.annotate(
            '', xy=(0,0),
            xytext=(20,20),
            textcoords='offset points',
            bbox=dict(boxstyle='round', fc='w'),
            arrowprops=dict(arrowstyle='->'))
        annot.set_visible(False)
        self._annot_cid = None

    def _hover(self, event):
        if self.current_model is None or self.current_model.was_deleted:
            return
        if getattr(self, '_tooltip_blocked', False):
            return
        sp = self.scatter
        annot = self._hover_annotation
        ax = self.axes
        if event.inaxes == ax:
            if self._debug:
                print('In axes')
            cont, ind = sp.contains(event)
            if cont:
                if self._debug:
                    print('In scatter')
                indices = ind['ind']
                text = []
                ramas = self._displayed_ramas[indices]
                phipsi = numpy.degrees(ramas.phipsis)
                x = phipsi[:,0].mean()
                y = phipsi[:,1].mean()
                for rama in ramas:
                    res = rama.residue
                    text.append('{} {}{}'.format(res.name, res.chain_id, res.number))
                if self._debug:
                    print('Annot text: {}'.format('\n'.join(text)))
                annot.set_text('\n'.join(text))
                annot.xy = x, y
                if y > 0:
                    annot.xyann=(20,-20)
                else:
                    annot.xyann=(20,20)
                annot.set_visible(True)
            else:
                annot.set_visible(False)
            self.canvas.restore_region(self.background)
            ax.draw_artist(self.scatter)
            ax.draw_artist(annot)
            self.canvas.blit(ax.bbox)

    def _start_tooltip(self, *_):
        if self._annot_cid is not None:
            self._stop_tooltip()
        self._annot_cid = self.canvas.mpl_connect('motion_notify_event', self._hover)

    def _stop_tooltip(self, *_):
        self.canvas.mpl_disconnect(self._annot_cid)
        self._annot_cid = None


    @property
    def current_model(self):
        return self.manager.current_model

    # @property
    # def selection_mode(self):
    #     return self._selection_mode

    # @selection_mode.setter
    # def selection_mode(self, mode):
    #     if mode not in self.mode_dict.keys():
    #         raise TypeError(f'Unrecognised mode: {mode}!')
    #     self._selection_mode = mode
    #     if self.current_model is None:
    #         self.set_target_residues(None)
    #     elif mode == self.WHOLE_MODEL:
    #         self.set_target_residues(self.current_model.residues)

    #     sflag = (mode == self.SELECTED_ONLY)
    #     self._sel_restrict_button.setEnabled(sflag)

    #     mflag = (mode != self.MOBILE_ONLY)
    #     self._mode_menu.setEnabled(mflag)

    #     self.update_scatter()

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
        return self.isVisible()

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
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    def _prepare_contours(self):
        import numpy
        from math import log
        key = self.case
        mgr = self._rama_mgr
        grid = numpy.degrees(mgr.interpolator_axes(int(key)))
        values = numpy.rot90(numpy.fliplr(mgr.interpolator_values(key))).astype(float)
        contours = mgr.RAMA_CASE_DETAILS[key]['cutoffs']
        self.P_limits = [0, -log(contours[0])]
        logvalues = numpy.log(values)
        self.contour_plot = self.axes.contour(*grid, values, contours)
        self.pcolor_plot = self.axes.pcolormesh(*grid, logvalues, cmap = 'Greys', shading='auto')

    def on_resize(self, *_):
        #background = self._background_cache.get(self.size(), None)
        axes = self.axes
        f = self.figure
        c = self.canvas
        self._hover_annotation.set_visible(False)
        self._format_axes()
        self.scatter.set_offsets(self.default_coords)
        f.patch.set_facecolor('black')
        c.draw()
        self.background = c.copy_from_bbox(axes.bbox)
        self.session.triggers.add_handler('new frame', self._resize_next_frame_cb)

    def _resize_next_frame_cb(self, *_):
        self.update_scatter()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def on_pick(self, event):
        ind = event.ind[0]
        picked_rama = self._displayed_ramas[ind]
        self.session.selection.clear()
        picked_rama.residue.atoms.selected=True
        from chimerax.isolde.navigate import get_stepper
        get_stepper(self.current_model).step_to(picked_rama.residue)
        # from .. import view
        # view.focus_on_selection(self.session, picked_rama.residue.atoms, pad=1.0)

    def clear(self):
        self.set_target_residues(None)
        self.update_scatter()

    def set_target_residues(self, residues):
        self._current_residues = residues
        case = self.case
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
        case = self.case
        mgr = self._rama_mgr
        cenum = self._case_enum
        r = self._case_ramas
        if r is None or len(r) == 0:
            self._displayed_ramas = None
            phipsi = self.default_coords
            logscores = self.default_logscores
            edge_colors='black'
            line_widths=0.5
        else:
            scores = r.scores
            if self.outliers_only or self.hide_favored:
                bins = mgr.bin_scores(r.scores, r.cases)
                if self.outliers_only:
                    mask = (bins == mgr.RamaBin.OUTLIER)
                else:
                    mask = (bins > mgr.RamaBin.FAVORED)
                r = r[mask]
                scores = scores[mask]
            self._displayed_ramas = r
            selecteds = r.ca_atoms.selecteds
            if numpy.any(selecteds):
                # Put the selected residues last so they show on top
                sort_order = numpy.lexsort((r.residues.numbers, r.residues.chain_ids, selecteds))
                r = self._displayed_ramas = r[sort_order]
                scores = scores[sort_order]
                selecteds = selecteds[sort_order]
            phipsi = numpy.degrees(r.phipsis)
            logscores = numpy.log(scores)
            edge_colors = numpy.zeros((len(phipsi),3))
            edge_colors[selecteds] = [0,1,0]
            line_widths = numpy.ones(len(phipsi))*0.5
            line_widths[selecteds] = 1

        c = self.canvas
        s = self.scatter
        axes = self.axes
        c.restore_region(self.background)
        s.set_offsets(phipsi)
        s.set_clim(self.P_limits)

        scales = (-logscores+1)*self.base_scatter_size
        s.set_sizes(scales)
        s.set_array(-logscores)
        s.set_edgecolors(edge_colors)
        s.set_linewidths(line_widths)
        axes.draw_artist(s)
        c.blit(axes.bbox)
