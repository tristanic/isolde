from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, HorizontalLine, ExpertModeSelector
from Qt.QtWidgets import QLabel, QSlider, QGridLayout, QWidget, QDoubleSpinBox, QPushButton, QCheckBox
from Qt.QtCore import Qt
from ..util import slot_disconnected

class NonBondedPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Nonbonded Potentials", expert_level=ExpertModeSelector.DEFAULT, **kwargs)
        nd = self.content = NonbondedDialog(session, isolde, gui, self)
        self.setContentLayout(nd.main_layout)

class NonbondedDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        cbl = DefaultHLayout()
        cb = self._use_softcore_potentials_checkbox = (QCheckBox('Use softcore nonbonded potentials'))
        cb.setChecked(isolde.sim_params.use_softcore_nonbonded_potential)
        cb.stateChanged.connect(self._use_softcore_potentials_checked_cb)
        cb.setToolTip('<span>Choose whether to use soft-core nonbonded potentials (will take effect on next sim start).'
            ' These can be very helpful in escaping from severe clashes, but invoke a ~2-fold penalty in simulation rate.</span>')
        cbl.addWidget(cb)
        ml.addLayout(cbl)
        sliders = self._sliders = []
        slider_layout = DefaultVLayout()
        for i, slider_class in enumerate((
            NbLambdaMinSlider, NbLambdaEquilSlider, NbASlider, NbBSlider, NbCSlider, NbAlphaSlider
        )):
            slider = slider_class(isolde.sim_params, gui.expert_mode_combo_box)
            sliders.append(slider)
            slider_layout.addWidget(slider)
        ml.addLayout(slider_layout)
        hl = DefaultVLayout()
        hl.addWidget(HorizontalLine())
        ml.addLayout(hl)
        
        pl = QGridLayout()
        pl.addWidget(QLabel('Minimization'), 1,0)
        mpi = self._min_potential_indicator = PotentialIndicator(isolde.sim_params, sliders[0],*sliders[2:])
        pl.addWidget(mpi, 2, 0)
        pl.addWidget(QLabel('Equilibration'), 1, 1)
        epi = self._equil_potential_indicator = PotentialIndicator(isolde.sim_params, *sliders[1:])
        pl.addWidget(epi, 2, 1)
        ml.addLayout(pl)
    
    def _use_softcore_potentials_checked_cb(self, checked):
        self.isolde.sim_params.use_softcore_nonbonded_potential = checked
    
    def cleanup(self):
        for slider in self._sliders.values():
            slider.cleanup()
        self._min_potential_indicator.cleanup()
        self._equil_potential_indicator.cleanup()



class ParamSlider(QWidget):
    NAME=None
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.DEFAULT
    PARAM_NAME = None
    def __init__(self, param_mgr, expert_mode_combo_box, *args, **kwargs):
        super().__init__(*args, **kwargs)
        ml = DefaultHLayout()
        self.setLayout(ml)
        l = self._label = QLabel(self.NAME)
        ml.addWidget(l)
        l.setMaximumWidth(60)
        l.setMinimumWidth(60)
        sl = self._slider = QSlider(Qt.Orientation.Horizontal, *args, **kwargs)
        ml.addWidget(sl)
        v = self._value_spinbox = QDoubleSpinBox()
        v.setMinimum(self.MIN_VAL*self.MULTIPLIER)
        v.setDecimals(2)
        ml.addWidget(v)
        v.setMaximumWidth(60)
        v.setMinimumWidth(60)
        v.valueChanged.connect(self._spin_box_changed_cb)
        rb = self._reset_button = QPushButton('Reset')
        rb.clicked.connect(self._reset_value_cb)
        ml.addWidget(rb)


        sl.setMinimum(self.MIN_VAL)
        sl.setMaximum(self.MAX_VAL)
        self._expert_mode_combo_box = expert_mode_combo_box
        self._set_expert_level(expert_mode_combo_box)
        self._param_mgr = param_mgr
        h = self._param_changed_handler = param_mgr.triggers.add_handler(param_mgr.PARAMETER_CHANGED, self._param_changed_cb)
        sl.valueChanged.connect(self._value_changed_cb)
        self._param_changed_cb(None, (self.PARAM_NAME, getattr(param_mgr, self.PARAM_NAME)))

    @property
    def scaled_value(self):
        return self._slider.value()*self.MULTIPLIER
    
    @scaled_value.setter
    def scaled_value(self, value):
        value = max(self.MIN_VAL*self.MULTIPLIER, min(value, self.MAX_VAL*self.MULTIPLIER))
        with self._param_changed_handler.blocked(), slot_disconnected(self._value_spinbox.valueChanged, self._spin_box_changed_cb):
            self._slider.setValue(value/self.MULTIPLIER)
            self._value_spinbox.setValue(value)
            
    def _spin_box_changed_cb(self, value):
        setattr(self._param_mgr, self.PARAM_NAME, value)

    def _set_expert_level(self, emcb):
        el = self.EXPERT_LEVEL
        if el > ExpertModeSelector.DEFAULT:
            stylesheet = ExpertModeSelector.stylesheets[el]
            self.setStyleSheet(stylesheet)
            emcb.currentIndexChanged.connect(self._expert_level_changed_cb)
            self._expert_level_changed_cb()
    
    def _expert_level_changed_cb(self, *_):
        emcb = self._expert_mode_combo_box
        el = emcb.currentData()
        display = (el >= self.EXPERT_LEVEL)
        self.setVisible(display)
    
    def _param_changed_cb(self, _, data):
        key, val = data
        if key == self.PARAM_NAME:
            with slot_disconnected(self._slider.valueChanged, self._value_changed_cb):
                v = val/self.MULTIPLIER
                self._slider.setValue(max(self.MIN_VAL,v, min(self.MAX_VAL,v)))
            with slot_disconnected(self._value_spinbox.valueChanged, self._spin_box_changed_cb):
                self._value_spinbox.setValue(val)
    
    def _value_changed_cb(self, _):
        with self._param_changed_handler.blocked():
            setattr(self._param_mgr, self.PARAM_NAME, self.scaled_value)
        with slot_disconnected(self._value_spinbox.valueChanged, self._spin_box_changed_cb):
            self._value_spinbox.setValue(self.scaled_value)
    
    def _reset_value_cb(self, *_):
        self._param_mgr.set_to_default(self.PARAM_NAME)

    def cleanup(self):
        self._param_changed_handler.remove()
            

class NbLambdaMinSlider(ParamSlider):
    NAME='λ (min)'
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    PARAM_NAME='nonbonded_softcore_lambda_minimize'

class NbLambdaEquilSlider(ParamSlider):
    NAME='λ (equil)'
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    PARAM_NAME='nonbonded_softcore_lambda_equil'

class NbASlider(ParamSlider):
    NAME='a'
    MIN_VAL=1
    MAX_VAL=20
    MULTIPLIER=1/10
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    PARAM_NAME='nonbonded_softcore_a'
    
class NbBSlider(ParamSlider):
    NAME='b'
    MIN_VAL=10
    MAX_VAL=100
    MULTIPLIER=1/25
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    PARAM_NAME='nonbonded_softcore_b'

class NbCSlider(ParamSlider):
    NAME='c'
    MIN_VAL=10
    MAX_VAL=100
    MULTIPLIER=1/10
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    PARAM_NAME='nonbonded_softcore_c'

class NbAlphaSlider(ParamSlider):
    NAME='Alpha'
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    PARAM_NAME='nonbonded_softcore_alpha'



class PotentialIndicator(QWidget):
    FONTSIZE=8
    def __init__(self, param_mgr, lambda_slider, a_slider, b_slider, c_slider, alpha_slider, *args, **kwargs):
        super().__init__(*args, **kwargs)
        import numpy
        from math import log
        radii = self._radii = numpy.exp(numpy.linspace(log(0.5),log(5),80))
        layout = DefaultHLayout()
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib import ticker


        self.setMaximumHeight(200)
        fig = self.figure = Figure()
        fig.set_tight_layout({'pad': 0.25})
        canvas = self.canvas = FigureCanvas(fig)
        loc = self._top_tick_locator = ticker.MultipleLocator(base=1e3)
        ax1, ax2 = fig.subplots(2,1, sharex=True)
        self._axes = (ax1, ax2)
        for ax in (ax1, ax2):
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(self.FONTSIZE)
        layout.addWidget(canvas)
        self.setLayout(layout)

        #ax2.yaxis.set_major_locator(loc)
        ax1.set_ylim([10,1e7])
        ax1.set_yscale('log')
        yticks = [10, 1e3, 1e5, 1e7]
        from math import log10
        ylabels = [f'$10^{{{int(log10(v)):01d}}}$' for v in yticks]
        ax1.set_yticks([10, 1e3, 1e5, 1e7])
        ax1.set_yticklabels(ylabels, fontsize=self.FONTSIZE)
        
        fig.text(0.02, 0.5, 'Energy (kJ/mol)', va='center', rotation='vertical', fontsize=self.FONTSIZE)
        ax1.set_ylabel(' ')
        #ax1.get_yaxis().set_major_formatter(ticker.LogFormatter())
        ax1.spines['bottom'].set_visible(False)
        ax2.set_ylim([-1,10])
        ax2.spines['top'].set_visible(False)
        ax2.set_xlabel('radius (Å)')
        ax2.set_ylabel(' ')

        # Plot unmodified values (lambda=1) as reference
        from chimerax.isolde.openmm.custom_forces import NonbondedSoftcoreForce
        lj, coul = NonbondedSoftcoreForce.potential_values(radii/10, 1, 1, 1, 1, 1)
        for ax in self._axes:
            line1 = ax.plot(radii, lj, 'm-', linewidth=2.5)[0]
            line1.set_label('L-J')
            line2 = ax.plot(radii, coul, 'g-', linewidth=2.5)[0]
            line2.set_label('Coul')
        ax1.legend(fontsize=self.FONTSIZE, loc='upper right')



        # # Diagonal lines to "cut" y-axis

        # d = .0015  # how big to make the diagonal lines in axes coordinates
        # # arguments to pass to plot, just so we don't keep repeating them
        # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        # ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        # ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        # kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        # ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        # ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        self._param_mgr = param_mgr
        self._changes_handler = param_mgr.triggers.add_handler(param_mgr.PARAMETER_CHANGED, self._update_plots_if_necessary)

        self._sliders = {
            'lambda': lambda_slider,
            'a':      a_slider,
            'b':      b_slider,
            'c':      c_slider,
            'alpha':  alpha_slider
        }
        self._update_plots()
    
    def _update_plots_if_necessary(self, _, data):
        name, val = data
        for slider in self._sliders.values():
            if name==slider.PARAM_NAME:
                self._update_plots()
                break
        
    def _update_plots(self, *_):
        from chimerax.isolde.openmm.custom_forces import NonbondedSoftcoreForce
        pmgr = self._param_mgr
        lj, coul = NonbondedSoftcoreForce.potential_values(self._radii/10, *(getattr(pmgr, slider.PARAM_NAME) for slider in self._sliders.values()))
        if not hasattr(self, '_lj_plots'):
            self._lj_plots = [ax.plot(self._radii, lj, 'm-', linewidth=1)[0] for ax in self._axes]
            self._coul_plots = [ax.plot(self._radii, coul, 'g-', linewidth=1)[0] for ax in self._axes]
        else:
            for plot in self._lj_plots:
                plot.set_ydata(lj)
            for plot in self._coul_plots:
                plot.set_ydata(coul)
            self.canvas.draw()
            self.canvas.flush_events()


    def cleanup(self):
        self._changes_handler.remove()



        

