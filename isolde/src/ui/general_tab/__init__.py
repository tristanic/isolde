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
        from .sim_selection import SimSelectionPanel
        self.addWidget(SimSelectionPanel(session, isolde, parent, gui))
        from .add_maps import MapAddPanel
        self.addWidget(MapAddPanel(session, isolde, parent, gui))
        from .param import ParameterisePanel
        self.addWidget(ParameterisePanel(session, isolde, parent, gui))
        from .map_settings import MapSettingsPanel
        self.addWidget(MapSettingsPanel(session, isolde, parent, gui))
        from .sim_fidelity import SimFidelityPanel
        self.addWidget(SimFidelityPanel(session, isolde, parent, gui))
        from .platform import ComputationalPlatformPanel
        self.addWidget(ComputationalPlatformPanel(session, isolde, parent, gui))
        from .mask_settings import MaskAndSpotlightSettingsPanel
        self.addWidget(MaskAndSpotlightSettingsPanel(session, isolde, parent, gui, start_collapsed=False))