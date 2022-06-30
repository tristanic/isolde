from Qt.QtWidgets import QToolButton, QWidget, QFrame, QGridLayout
from Qt.QtCore import Qt, QParallelAnimationGroup, QPropertyAnimation, Signal, QTimer
from Qt.QtWidgets import QSizePolicy

from .ui_base import ExpertModeSelector

class CollapsibleArea(QWidget):
    collapsed = Signal()
    expanded = Signal()
    def __init__(self, gui, parent=None, title='', start_collapsed=True, expert_level=ExpertModeSelector.DEFAULT):
        super().__init__(parent=parent)
        self.gui=gui
        self._is_collapsed = start_collapsed
        self.expert_level = expert_level
        ca = self.content_area = QFrame()
        hl = self.header_line = QFrame()
        tb = self.toggle_button = QToolButton()
        ml = self.main_layout = QGridLayout()

        tb.setStyleSheet('QToolButton: {border: none; }')
        tb.setAutoRaise(True)
        tb.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        tb.setText(title)
        tb.setCheckable(True)
        if start_collapsed:
            tb.setArrowType(Qt.RightArrow)
            tb.setChecked(False)
        else:
            tb.setArrowType(Qt.DownArrow)
            tb.setChecked(True)
        tb.setVisible(True)
        tb.toggled.connect(self._show_or_hide)

        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        hl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        hl.setMinimumHeight(tb.sizeHint().height())
        hl.setMaximumHeight(tb.sizeHint().height())

        ca.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.MinimumExpanding)
        if start_collapsed:
            ca.setVisible(False)
        
        ml.setSpacing(0)
        ml.setContentsMargins(0,0,0,0)

        ml.addWidget(tb, 0,0,1,1, Qt.AlignLeft)
        ml.addWidget(hl, 0,2,1,1)
        ml.addWidget(ca, 1,0,1,3)

        self.setLayout(ml)
        self._set_expert_level()

    @property
    def is_collapsed(self):
        return not self.toggle_button.isChecked()
    
    def expand(self):
        if not self.is_collapsed:
            return
        self.toggle_button.toggle()
    
    def collapse(self):
        if self.is_collapsed:
            return
        self.toggle_button.toggle()

    def _show_or_hide(self, checked):
        forward = checked
        ca = self.content_area

        if forward:
            ca.setVisible(True)
            # Need to delay the event to the next redraw, otherwise some things 
            # (e.g. MatPlotLib redraws) go wrong
            QTimer.singleShot(5, lambda: self.expanded.emit())
        else:
            ca.setVisible(False)
            QTimer.singleShot(5, lambda: self.collapsed.emit())

        arrow_type = Qt.DownArrow if checked else Qt.RightArrow
        self.toggle_button.setArrowType(arrow_type)            

    def setContentLayout(self, layout):
        ca = self.content_area
        ca.setLayout(layout)


    def _set_expert_level(self):
        el = self.expert_level
        emcb = self.gui.expert_mode_combo_box
        if el > ExpertModeSelector.DEFAULT:
            stylesheet = ExpertModeSelector.stylesheets[el]
            self.setStyleSheet(stylesheet)
            emcb.currentIndexChanged.connect(self._expert_level_changed_cb)
            self._expert_level_changed_cb()
    
    def _expert_level_changed_cb(self, *_):
        emcb = self.gui.expert_mode_combo_box
        el = emcb.currentData()
        display = (el >= self.expert_level)
        self.setVisible(display)
