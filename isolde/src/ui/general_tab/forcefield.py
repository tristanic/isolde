# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Developer-level force-field selector.

Reinstates the (long-dormant) user-selectable force-field framework: a combo box
under the **Developer** experience level letting the user pick the MD force field
for the next simulation. AMBER (``amber14``) is the default; ``garnet`` selects
the experimental garnet-isolde graph-ML force field. Changing the selection emits
``isolde set forcefield "<name>"`` (so the action is logged and scriptable, per
ISOLDE's GUI-logs-its-command convention).
'''
from ..ui_base import (
    DefaultSpacerItem,
    UI_Panel_Base,
    DefaultHLayout,
    ExpertModeSelector,
)
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QLabel, QComboBox


class ForceFieldPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='MD Force Field',
                         expert_level=ExpertModeSelector.DEVELOPER, **kwargs)
        ffd = self.content = ForceFieldDialog(session, isolde, gui, self)
        self.setContentLayout(ffd.main_layout)


class ForceFieldDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area)
        self.container = collapse_area

        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()
        ml.addWidget(QLabel('Force field: ', mf))
        fcb = self.forcefield_combo_box = QComboBox(mf)
        fcb.setToolTip('Force field used to parameterise the next simulation. '
                       '"amber14" is the standard, well-tested option; "garnet" '
                       'is the experimental garnet-isolde graph-ML force field.')
        self._populate_forcefield_combo_box()
        # Connect AFTER populating so the initial setCurrentIndex does not fire the command.
        fcb.currentIndexChanged.connect(self._choose_forcefield_cb)
        ml.addWidget(fcb)
        ml.addItem(DefaultSpacerItem())

    def _populate_forcefield_combo_box(self):
        sim_params = self.session.isolde.sim_params
        names = list(self.session.isolde.forcefield_mgr.available_forcefields)
        # Stable order with the legacy default first.
        names = sorted(names, key=lambda n: (n != 'amber14', n))
        fcb = self.forcefield_combo_box
        fcb.addItems(names)
        current = getattr(sim_params, 'forcefield', 'amber14')
        if current in names:
            fcb.setCurrentIndex(fcb.findText(current))

    def _choose_forcefield_cb(self, *_):
        name = self.forcefield_combo_box.currentText()
        if not name or name == getattr(self.session.isolde.sim_params, 'forcefield', None):
            return
        from chimerax.core.commands import run
        run(self.session, 'isolde set forcefield "{}"'.format(name))
