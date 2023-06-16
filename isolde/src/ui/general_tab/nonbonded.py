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
            ' These can be very helpful in escaping from severe clashes, but invoke a 10-20% penalty in simulation rate.</span>')
        cbl.addWidget(cb)
        ml.addLayout(cbl)
        sliders = self._sliders = []
        slider_layout = DefaultVLayout()
        for i, slider_class in enumerate((
            NbLambdaMinSlider, NbLambdaEquilSlider, NbAlphaSlider
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
        mpi = self._min_potential_indicator = PotentialIndicator(isolde.sim_params, sliders[0],sliders[2])
        mpi.setToolTip('<span>Nonbonded potential used during energy minimisation (at simulation startup, after non-simulation-derived coordinate changes, '
            'or on detection of overly-fast atom movements). A lower value of λ improves the ability of the minimiser to escape from severe clashes.</span>')
        pl.addWidget(mpi, 2, 0)
        pl.addWidget(QLabel('Equilibration'), 1, 1)
        epi = self._equil_potential_indicator = PotentialIndicator(isolde.sim_params, *sliders[1:])
        epi.setToolTip('<span>Nonbonded potential used during standard dynamics. In the absence of severe pathologies this should be kept close to '
            'the normal Lennard-Jones potential (recommended λ=0.95). λ=1 replicates the L-J potential exactly.</span>')
        pl.addWidget(epi, 2, 1)
        ml.addLayout(pl)
    
    def _use_softcore_potentials_checked_cb(self, checked):
        self.isolde.sim_params.use_softcore_nonbonded_potential = checked
    
    def cleanup(self):
        for slider in self._sliders:
            slider.cleanup()
        self._min_potential_indicator.cleanup()
        self._equil_potential_indicator.cleanup()



class ParamSlider(QWidget):
    NAME=None
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.DEFAULT
    TOOLTIP=''
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
        self.setToolTip(self.TOOLTIP)


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
                v = round(val/self.MULTIPLIER)
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
    TOOLTIP=('<span>Controlling parameter setting the nonbonded potential used during energy minimisation '
             ' (at simulation startup, after non-simulation-derived coordinate changes, or on detection of overly-fast'
            ' atom movements). A lower value of λ improves the ability of the minimiser to escape from severe clashes.</span>')

class NbLambdaEquilSlider(ParamSlider):
    NAME='λ (equil)'
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    PARAM_NAME='nonbonded_softcore_lambda_equil'
    TOOLTIP=('<span>Controlling parameter setting the nonbonded potential used during standard dynamics. '
             'In most cases this should be kept at a value close to 1 (recommended λ=0.95). If severe clashes are '
             'causing large distortions in the model, temporarily reducing λ to around 0.5 or below can be helpful. '
             'λ=1 replicates the Lennard-Jones potential exactly.</span>')

class NbAlphaSlider(ParamSlider):
    NAME='α'
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    PARAM_NAME='nonbonded_softcore_alpha'
    TOOLTIP=('<span>Provides more subtle control over the shape of the potential within the "clashing" region. '
             'In general it should not be necessary to change this. </span>')



class PotentialIndicator(QWidget):
    FONTSIZE=8
    def __init__(self, param_mgr, lambda_slider, alpha_slider, *args, **kwargs):
        super().__init__(*args, **kwargs)
        import numpy
        from math import log
        radii = self._radii = numpy.exp(numpy.linspace(log(0.5),log(5),80))
        layout = DefaultHLayout()
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

        self.setMaximumHeight(200)
        fig = self.figure = Figure()
        fig.set_tight_layout({'pad': 0.25})
        canvas = self.canvas = FigureCanvas(fig)
        ax1, ax2 = fig.subplots(2,1, sharex=True)
        self._axes = (ax1, ax2)
        for ax in (ax1, ax2):
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(self.FONTSIZE)
        layout.addWidget(canvas)
        self.setLayout(layout)

        ax1.set_ylim([10,1e7])
        ax1.set_yscale('log')
        yticks = [10, 1e3, 1e5, 1e7]
        from math import log10
        ylabels = [f'$10^{{{int(log10(v)):01d}}}$' for v in yticks]
        ax1.set_yticks([10, 1e3, 1e5, 1e7])
        ax1.set_yticklabels(ylabels, fontsize=self.FONTSIZE)
        
        fig.text(0.02, 0.5, 'Energy (kJ/mol)', va='center', rotation='vertical', fontsize=self.FONTSIZE)
        ax1.set_ylabel(' ')
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

        self._param_mgr = param_mgr
        self._changes_handler = param_mgr.triggers.add_handler(param_mgr.PARAMETER_CHANGED, self._update_plots_if_necessary)

        self._sliders = {
            'lambda': lambda_slider,
            'alpha':  alpha_slider
        }
        self._update_plots()
    
    def _update_plots_if_necessary(self, _, data):
        name, val = data
        if name in (
            self._sliders['lambda'].PARAM_NAME,
            'nonbonded_softcore_a',
            'nonbonded_softcore_b',
            'nonbonded_softcore_c',
            'nonbonded_softcore_alpha'
        ):
            self._update_plots()
        
    def _update_plots(self, *_):
        from chimerax.isolde.openmm.custom_forces import NonbondedSoftcoreForce
        pmgr = self._param_mgr
        l = getattr(pmgr, self._sliders['lambda'].PARAM_NAME)
        a = pmgr.nonbonded_softcore_a
        b = pmgr.nonbonded_softcore_b
        c = pmgr.nonbonded_softcore_c
        alpha = pmgr.nonbonded_softcore_alpha
        lj, coul = NonbondedSoftcoreForce.potential_values(self._radii/10, l, a, b, c, alpha)
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



        

