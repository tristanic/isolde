from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, HorizontalLine, ExpertModeSelector
from Qt.QtWidgets import QLabel, QPushButton, QMenu, QToolBar, QSlider, QGridLayout, QWidget
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class NonBondedPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Nonbonded Potentials", expert_level=ExpertModeSelector.ADVANCED, **kwargs)
        nd = self.content = NonbondedDialog(session, isolde, gui, self)
        self.setContentLayout(nd.main_layout)

class NonbondedDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

class ScaledSlider(QSlider):
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.ADVANCED
    def __init__(self, expert_mode_combo_box, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(self.MIN_VAL)
        self.setMaximum(self.MAX_VAL)
        self._expert_mode_combo_box = expert_mode_combo_box
        self._set_expert_level(expert_mode_combo_box)

    @property
    def scaled_value(self):
        return self.value()*self.MULTIPLIER
    
    @scaled_value.setter
    def scaled_value(self, value):
        value = max(self.MIN_VAL*self.MULTIPLIER, min(value, self.MAX_VAL*self.MULTIPLIER))
        self.setValue(value/self.MULTIPLIER)

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


class NbLambdaSlider(ScaledSlider):
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100

class NbASlider(ScaledSlider):
    MIN_VAL=1
    MAX_VAL=20
    MULTIPLIER=1/10
    EXPERT_LEVEL=ExpertModeSelector.DEVELOPER
    
class NbBSlider(ScaledSlider):
    MIN_VAL=10
    MAX_VAL=100
    MULTIPLIER=1/25
    EXPERT_LEVEL=ExpertModeSelector.DEVELOPER

class NbCSlider(ScaledSlider):
    MIN_VAL=10
    MAX_VAL=100
    MULTIPLIER=1/10
    EXPERT_LEVEL=ExpertModeSelector.DEVELOPER

class NbAlphaSlider(ScaledSlider):
    MIN_VAL=1
    MAX_VAL=100
    MULTIPLIER=1/100
    EXPERT_LEVEL=ExpertModeSelector.DEVELOPER



class PotentialIndicator(QWidget):
    def __init__(self, lambda_slider, a_slider, b_slider, c_slider, alpha_slider, *args, **kwargs):
        super().__init__(*args, **kwargs)
        import numpy
        self._radii = numpy.linspace(0.1,5,40)
        layout = DefaultHLayout()
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib import ticker


        self.setMaximumHeight(120)
        fig = self.figure = Figure()
        fig.set_tight_layout({'pad': 0.25})
        canvas = self.canvas = FigureCanvas(fig)
        loc = self._top_tick_locator = ticker.MultipleLocator(base=1e3)
        ax1, ax2 = fig.subplots(2,1)
        for ax in (ax1, ax2):
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)
        layout.addWidget(canvas)
        self.setLayout(layout)

        ax1.yaxis.set_major_locator(loc)
        ax1.set_ylim([-5,20])
        ax2.set_ylim([0,1e8])

        self._lambda_slider = lambda_slider
        self._a_slider = a_slider
        self._b_slider = b_slider
        self._c_slider = c_slider
        self._alpha_slider = alpha_slider

    def _update_plots(self, *_):
        from chimerax.isolde.openmm.custom_forces import NonbondedSoftcoreForce
        lj, coul = NonbondedSoftcoreForce.potential_values(self._radii, self._lambda_slider.scaled_value, 
            self._a_slider.scaled_value, self._b_slider.scaled_value, self._c_slider.scaled_value,
            self._alpha_slider.scaled_value)







        

