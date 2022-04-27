from ..ui_base import IsoldeTab


class GeneralTab(IsoldeTab):
    def __init__(self, session, isolde, gui, tab_widget, tab_name):
        super().__init__(gui, tab_widget, tab_name)
        self.session = session
        self.isolde = isolde
        self._populate()    

    def _populate(self):
        session = self.session
        isolde = self.isolde
        parent = self.scroll_area
        gui = self.gui
        from .map_settings import MapSettingsPanel
        self.addWidget(MapSettingsPanel(session, isolde, parent, gui))
        from .sim_fidelity import SimFidelityPanel
        self.addWidget(SimFidelityPanel(session, isolde, parent, gui))
