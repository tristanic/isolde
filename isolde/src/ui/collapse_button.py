from Qt.QtWidgets import QToolButton, QWidget, QFrame, QGridLayout
from Qt.QtCore import Qt, QParallelAnimationGroup, QPropertyAnimation, Signal
from Qt.QtWidgets import QSizePolicy

from .ui_base import ExpertModeSelector

class CollapsibleArea(QWidget):
    collapsed = Signal()
    expanded = Signal()
    def __init__(self, gui, parent=None, title='', duration=100, start_collapsed=True, expert_level=ExpertModeSelector.DEFAULT):
        super().__init__(parent=parent)
        self.gui=gui
        self._is_collapsed = start_collapsed
        self.animation_duration=duration
        self.expert_level = expert_level
        anim = self._animator = QParallelAnimationGroup()
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
        tb.toggled.connect(self._do_animation)

        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        hl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        hl.setMinimumHeight(tb.sizeHint().height())

        ca.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        if start_collapsed:
            ca.setMaximumHeight(0)
            ca.setMinimumHeight(0)
        #ca.setVisible(False)

        anim.addAnimation(QPropertyAnimation(self, b'minimumHeight'))
        anim.addAnimation(QPropertyAnimation(self, b'maximumHeight'))
        anim.addAnimation(QPropertyAnimation(ca, b'maximumHeight'))
        anim.finished.connect(self._animation_finished_cb)
        
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

    def _do_animation(self, checked):
        ca = self.content_area
        collapsed_height = self.sizeHint().height() - ca.height()
        content_height = ca.sizeHint().height()
        anim = self._animator
        for i in range(anim.animationCount()-1):
            a = anim.animationAt(i)
            a.setStartValue(collapsed_height)
            a.setEndValue(collapsed_height+content_height)
        a = anim.animationAt(anim.animationCount()-1)
        a.setStartValue(0)
        a.setEndValue(content_height)        

        from Qt.QtCore import QAbstractAnimation
        arrow_type = Qt.DownArrow if checked else Qt.RightArrow
        direction = QAbstractAnimation.Forward if checked else QAbstractAnimation.Backward
        self.toggle_button.setArrowType(arrow_type)
        self._animator.setDirection(direction)
        self._animator.start()
    
    def _animation_finished_cb(self, *_):
        if self.is_collapsed:
            self.collapsed.emit()
        else:
            self.expanded.emit()

    def setContentLayout(self, layout):
        ca = self.content_area
        ca.setLayout(layout)
        anim = self._animator
        for i in range(anim.animationCount()):
            a = anim.animationAt(i)
            a.setDuration(self.animation_duration)


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
