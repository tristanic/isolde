from Qt.QtWidgets import QToolButton, QWidget, QScrollArea, QFrame, QGridLayout
from Qt.QtCore import QSize, Qt, QParallelAnimationGroup, QPropertyAnimation
from Qt.QtWidgets import QSizePolicy


class CollapsibleArea(QWidget):
    def __init__(self, parent=None, title='', duration=100):
        super().__init__(parent=parent)
        self.animation_duration=duration
        anim = self._animator = QParallelAnimationGroup()
        ca = self.content_area = QFrame()
        hl = self.header_line = QFrame()
        tb = self.toggle_button = QToolButton()
        ml = self.main_layout = QGridLayout()

        tb.setStyleSheet('QToolButton: {border: none; }')
        tb.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        tb.setArrowType(Qt.RightArrow)
        tb.setText(title)
        tb.setCheckable(True)
        tb.setChecked(False)
        tb.setVisible(True)
        tb.toggled.connect(self._do_animation)

        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        hl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)

        ca.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        # Start out collapsed
        ca.setMaximumHeight(0)
        ca.setMinimumHeight(0)
        #ca.setVisible(False)

        anim.addAnimation(QPropertyAnimation(self, b'minimumHeight'))
        anim.addAnimation(QPropertyAnimation(self, b'maximumHeight'))
        anim.addAnimation(QPropertyAnimation(ca, b'maximumHeight'))
        
        ml.setVerticalSpacing(0)
        ml.setContentsMargins(0,0,0,0)

        ml.addWidget(tb, 0,0,1,1, Qt.AlignLeft)
        ml.addWidget(hl, 0,2,1,1)
        ml.addWidget(ca, 1,0,1,3)

        self.setLayout(ml)

    def _do_animation(self, checked):
        from Qt.QtCore import QAbstractAnimation
        arrow_type = Qt.DownArrow if checked else Qt.RightArrow
        direction = QAbstractAnimation.Forward if checked else QAbstractAnimation.Backward
        #self.content_area.setVisible(checked)
        self.toggle_button.setArrowType(arrow_type)
        self._animator.setDirection(direction)
        self._animator.start()
    
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

