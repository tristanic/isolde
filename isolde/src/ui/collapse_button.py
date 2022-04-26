from Qt.QtWidgets import QToolButton, QWidget, QFrame, QGridLayout
from Qt.QtCore import Qt, QParallelAnimationGroup, QPropertyAnimation, Signal
from Qt.QtWidgets import QSizePolicy


class CollapsibleArea(QWidget):
    collapsed = Signal()
    expanded = Signal()
    def __init__(self, parent=None, title='', duration=100, start_collapsed=True):
        super().__init__(parent=parent)
        self._is_collapsed = start_collapsed
        self.animation_duration=duration
        anim = self._animator = QParallelAnimationGroup()
        ca = self.content_area = QFrame()
        hl = self.header_line = QFrame()
        tb = self.toggle_button = QToolButton()
        ml = self.main_layout = QGridLayout()

        tb.setStyleSheet('QToolButton: {border: none; }')
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
        hl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)

        ca.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        if start_collapsed:
            ca.setMaximumHeight(0)
            ca.setMinimumHeight(0)
        #ca.setVisible(False)

        anim.addAnimation(QPropertyAnimation(self, b'minimumHeight'))
        anim.addAnimation(QPropertyAnimation(self, b'maximumHeight'))
        anim.addAnimation(QPropertyAnimation(ca, b'maximumHeight'))
        anim.finished.connect(self._animation_finished_cb)
        
        ml.setVerticalSpacing(0)
        ml.setContentsMargins(0,0,0,0)

        ml.addWidget(tb, 0,0,1,1, Qt.AlignLeft)
        ml.addWidget(hl, 0,2,1,1)
        ml.addWidget(ca, 1,0,1,3)

        self.setLayout(ml)

    @property
    def is_collapsed(self):
        return not self.toggle_button.isChecked()

    def _do_animation(self, checked):
        from Qt.QtCore import QAbstractAnimation
        arrow_type = Qt.DownArrow if checked else Qt.RightArrow
        direction = QAbstractAnimation.Forward if checked else QAbstractAnimation.Backward
        #self.content_area.setVisible(checked)
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
        collapsed_height = self.sizeHint().height() - ca.maximumHeight()
        content_height = ca.sizeHint().height()
        anim = self._animator
        for i in range(anim.animationCount()-1):
            a = anim.animationAt(i)
            a.setDuration(self.animation_duration)
            a.setStartValue(collapsed_height)
            a.setEndValue(collapsed_height+content_height)
        a = anim.animationAt(anim.animationCount()-1)
        a.setDuration(self.animation_duration)
        a.setStartValue(0)
        a.setEndValue(content_height)        

