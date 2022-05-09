from chimerax.ui.gui import MainToolWindow

from chimerax.core.tools import ToolInstance
from Qt.QtCore import Qt

from chimerax.isolde.ui.ui_base import (
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem, QComboBox
)

from Qt.QtWidgets import (
    QGridLayout, QFrame, QLabel, QCheckBox
)

class Rama_ToolUI(ToolInstance):
    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        tw = self.tool_window=RamaMainWin(self)
        tw.manage(placement=None, allowed_areas=Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)


class RamaMainWin(MainToolWindow):
    def __init__(self, tool_instance, **kw):
        super().__init__(tool_instance, **kw)
    
        parent = self.ui_area
        main_layout = self.main_layout = QGridLayout(parent)
        parent.setLayout(main_layout)
