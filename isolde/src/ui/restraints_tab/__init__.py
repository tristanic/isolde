from ..ui_base import IsoldeTab

class RestraintsTab(IsoldeTab):
    display_name = 'Restraints'
    def populate(self):
        session = self.session
        isolde = self.isolde
        parent = self.scroll_area
        gui = self.gui

        from .position import PositionRestraintsPanel
        self.addWidget(PositionRestraintsPanel(session, isolde, parent, gui))
        from .distance import DistanceRestraintsPanel
        self.addWidget(DistanceRestraintsPanel(session, isolde, parent, gui))
        from .secondary_structure import SecondaryStructurePanel
        self.addWidget(SecondaryStructurePanel(session, isolde, parent, gui))
        from .base_pairing import BasePairingPanel
        self.addWidget(BasePairingPanel(session, isolde, parent, gui))
        from .register_shift import RegisterShiftPanel
        self.addWidget(RegisterShiftPanel(session, isolde, parent, gui))  
        from .reference_model import ReferenceModelPanel
        self.addWidget(ReferenceModelPanel(session, isolde, parent, gui))
        from .manage_restraints import ManageRestraintsPanel
        self.addWidget(ManageRestraintsPanel(session, isolde, parent, gui))      