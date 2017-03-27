# vim: set expandtab shiftwidth=4 softtabstop=4:

# ToolUI should inherit from ToolInstance if they will be
# registered with the tool state manager.
#
# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
from chimerax.core.tools import ToolInstance


class MolProbity_ToolUI(ToolInstance):

    SESSION_ENDURING = False
    # if SESSION_ENDURING is True, tool instance not deleted at session closure

    def __init__(self, session, tool_name):
        ToolInstance.__init__(self, session, tool_name)
        self.display_name = "MolProbity"
        from chimerax.core.ui.gui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        self.tool_window.manage(placement=None)
        parent = self.tool_window.ui_area
        pp = parent.parent()
        pp.resize(480,850) 

        from PyQt5 import QtWidgets, QtGui
        from . import molprobity_widget as mw
        self.mainwin = QtWidgets.QFrame(parent=parent)
        self.mw = mw.Ui_molprobity_widget()
        self.mw.setupUi(self.mainwin)
        
        import os
        
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.mainwin)
        layout.setStretchFactor(self.mainwin, 1)
        parent.setLayout(layout)
        self.tool_window.manage(placement=None)
        # Should load saved state here
                
        # Any other frames/widgets that should be hidden at the start
        self._hidden_at_start = [
            self.mw._validate_rama_main_frame,
            self.mw._validate_pep_main_frame,
            ]
        
        for f in self._hidden_at_start:
            f.hide()
        
    
