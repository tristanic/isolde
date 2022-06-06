from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QPushButton, QRadioButton
from Qt.QtGui import QColor, QBrush
from Qt.QtCore import Qt

class RamaPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Ramachandran Plot", **kwargs)
        rd = self.content = RamaDialog(session, isolde, gui, self)
        self.setContentLayout(rd.main_layout)

class RamaDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=True):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)

        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        tl = DefaultHLayout()
        sb = self.show_button = QPushButton('Launch Ramachandran plot')
        tl.addWidget(sb)
        tl.addStretch()
        ml.addLayout(tl)

        from chimerax.core.commands import run
        sb.clicked.connect(self._launch_ramaplot)

        bl = DefaultHLayout()
        show_all = self.show_all_residues = QRadioButton('Show all residues')
        show_mob = self.show_mobile_residues = QRadioButton('Limit to mobile residues')
        bl.addWidget(show_all)
        bl.addWidget(show_mob)
        bl.addStretch()
        show_mob.setChecked(True)
        show_all.clicked.connect(self._radio_button_cb)
        show_mob.clicked.connect(self._radio_button_cb)
        ml.addLayout(bl)

    def _radio_button_cb(self, *_):
        show_mob = self.show_mobile_residues.isChecked()
        if show_mob:
            self.restrict_to_mobile_selection()
            return
        rplot = self._get_rplot()
        if rplot is not None:
            rplot.display_all_residues()

    def _get_rplot(self):
        from chimerax.isolde.validation.ramaplot.tool import Rama_ToolUI
        tools = self.session.tools.find_by_class(Rama_ToolUI)
        if not len(tools):
            return None
        return tools[0].tool_window

    def _launch_ramaplot(self, *_):
        from chimerax.core.commands import run
        run(self.session, 'ui tool show "Ramachandran Plot"')
        if self.show_mobile_residues:
            self.restrict_to_mobile_selection()


    def restrict_to_mobile_selection(self, *_):
        if not self.isolde.simulation_running:
            return
        rplot = self._get_rplot()
        if rplot is None:
            return
        if rplot.current_model == self.isolde.selected_model:
            rplot.restrict_to_selection(self.isolde.sim_manager.sim_construct.mobile_residues, display_text='Mobile residues')
  
    def sim_start_cb(self, *_):
        if self.show_mobile_residues.isChecked():
            self.restrict_to_mobile_selection()