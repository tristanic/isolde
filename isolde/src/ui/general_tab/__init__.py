from ..ui_base import IsoldeTab


class GeneralTab(IsoldeTab):

    def populate(self):
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
        from .map_settings import XmapLiveSettingsPanel, XmapStaticSettingsPanel, NXmapSettingsPanel
        self.addWidget(XmapLiveSettingsPanel(session, isolde, parent, gui))
        self.addWidget(XmapStaticSettingsPanel(session, isolde, parent, gui))
        self.addWidget(NXmapSettingsPanel(session, isolde, parent, gui))
        from .sim_fidelity import SimFidelityPanel
        self.addWidget(SimFidelityPanel(session, isolde, parent, gui))
        from .platform import ComputationalPlatformPanel
        self.addWidget(ComputationalPlatformPanel(session, isolde, parent, gui))
        from .mask_settings import MaskAndSpotlightSettingsPanel
        self.addWidget(MaskAndSpotlightSettingsPanel(session, isolde, parent, gui, start_collapsed=False))
        from .sim_runtime import SimRuntimePanel
        self.addWidget(SimRuntimePanel(session, isolde, parent, gui, start_collapsed=False))