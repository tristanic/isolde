from ..ui_base import IsoldeTab

from Qt.QtWidgets import QPushButton
from chimerax.core.commands import run

class ValidationTab(IsoldeTab):
    def populate(self):
        session = self.session
        isolde = self.isolde
        parent = self.scroll_area
        gui = self.gui

        from .peptide_bond import PeptideBondPanel
        self.addWidget(PeptideBondPanel(session, isolde, parent, gui))
        from .rama import RamaPanel
        self.addWidget(RamaPanel(session, isolde, parent, gui))
        from .rotamer import RotamerPanel
        self.addWidget(RotamerPanel(session, isolde, parent, gui))
        from .clashes import ClashesPanel
        self.addWidget(ClashesPanel(session, isolde, parent, gui))
        from .unparameterised import UnparameterisedResiduesPanel
        self.addWidget(UnparameterisedResiduesPanel(session, isolde, parent, gui))