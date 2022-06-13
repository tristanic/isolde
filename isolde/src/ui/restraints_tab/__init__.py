from ..ui_base import IsoldeTab

class RestraintsTab(IsoldeTab):
    def populate(self):
        session = self.session
        isolde = self.isolde
        parent = self.scroll_area
        gui = self.gui

        from .position import PositionRestraintsPanel
        self.addWidget(PositionRestraintsPanel(session, isolde, parent, gui))
        from .distance import DistanceRestraintsPanel
        self.addWidget(DistanceRestraintsPanel(session, isolde, parent, gui))
        from .register_shift import RegisterShiftPanel
        self.addWidget(RegisterShiftPanel(session, isolde, parent, gui))        