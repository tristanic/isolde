from ..ui_base import IsoldeTab

from Qt.QtWidgets import QPushButton
from chimerax.core.commands import run

class ValidationTab(IsoldeTab):
    def populate(self):
        session = self.session
        isolde = self.isolde
        parent = self.scroll_area
        gui = self.gui

        rpb = QPushButton('Launch Ramachandran plot', parent)
        rpb.clicked.connect(lambda: run(session, 'ui tool show "Ramachandran Plot"'))
        self.addWidget(rpb)

        from .peptide_bond import PeptideBondPanel
        self.addWidget(PeptideBondPanel(session, isolde, parent, gui))