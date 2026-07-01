# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Forcefield-selection UI panel for the ISOLDE General tab.

Appears as a collapsible panel (ADVANCED expert level) between the
NonBonded settings and the Computational Platform panels.  The combo box
is populated from
:func:`chimerax.isolde.openmm.param_provider.available_parameterisation_backends`.

When ``'amber14+garnet'`` is listed but GARNET cannot be imported it is
shown greyed-out with install instructions in a tooltip.
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
        super().__init__(gui, parent, title='Forcefield',
                         expert_level=ExpertModeSelector.ADVANCED, **kwargs)
        ffd = self.content = ForceFieldDialog(session, isolde, gui, self)
        self.setContentLayout(ffd.main_layout)


class ForceFieldDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area)
        self.container = collapse_area

        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()
        cb = self.combo_box = QComboBox(mf)
        self._populate_combo_box()
        cb.currentIndexChanged.connect(self._combo_changed_cb)
        ml.addWidget(cb)
        ml.addItem(DefaultSpacerItem())

    def _populate_combo_box(self):
        from chimerax.isolde.openmm.param_provider import _garnet_available

        sim_params = self.session.isolde.sim_params
        cb = self.combo_box
        cb.blockSignals(True)
        cb.clear()

        # All backends in display order; ML items are always shown so users
        # know the option exists even when dependencies are absent.
        all_backends = ['amber14', 'amber14+garnet', 'charmm36']
        garnet_installed = _garnet_available()

        for name in all_backends:
            cb.addItem(name)
            if name == 'amber14+garnet' and not garnet_installed:
                model = cb.model()
                item = model.item(cb.count() - 1)
                item.setEnabled(False)
                item.setToolTip(
                    'GARNET is not installed.\n'
                    'To enable, run inside ChimeraX:\n'
                    '  pip install garnetff torch torch_geometric igraph\n'
                    'then restart ChimeraX.')

        current = getattr(sim_params, 'forcefield', 'amber14')
        idx = cb.findText(current)
        cb.setCurrentIndex(idx if idx >= 0 else 0)
        cb.blockSignals(False)

    def _combo_changed_cb(self, *_):
        chosen = self.combo_box.currentText()
        from chimerax.isolde.openmm.param_provider import _garnet_available
        if chosen == 'amber14+garnet' and not _garnet_available():
            return
        self.session.isolde.sim_params.forcefield = chosen
