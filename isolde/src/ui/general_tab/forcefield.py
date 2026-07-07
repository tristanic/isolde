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
        from chimerax.isolde.openmm.param_provider import (
            _garnet_available, _espaloma_available, _espaloma_charge_available,
            _openff_available, _mmff_available)

        sim_params = self.session.isolde.sim_params
        cb = self.combo_box
        cb.blockSignals(True)
        cb.clear()

        # All backends in display order; unavailable items are shown greyed-out
        # so users know the option exists even when dependencies are absent.
        all_backends = ['amber14', 'amber14+garnet', 'amber14+espaloma',
                        'amber14+espaloma-charge', 'amber14+openff',
                        'amber14+mmff', 'charmm36']

        _disabled = {
            'amber14+garnet': (not _garnet_available(),
                'GARNET is not installed.\n'
                'To enable, run inside ChimeraX:\n'
                '  pip install garnetff torch torch_geometric igraph\n'
                'then restart ChimeraX.'),
            'amber14+espaloma': (not _espaloma_available(),
                'Espaloma is not installed.\n'
                'Requires conda — see the ISOLDE docs.'),
            'amber14+espaloma-charge': (not _espaloma_charge_available(),
                'espaloma-charge is not installed.\n'
                'GNN charges (AM1-BCC quality) + MMFF94 bonded terms.\n'
                'To enable, run in an admin cmd:\n'
                '  python.exe -m pip install espaloma-charge dgl\n'
                'then restart ChimeraX.'),
            'amber14+openff': (not _openff_available(),
                'OpenFF Toolkit is not installed.\n'
                'To enable, run inside ChimeraX:\n'
                '  pip install openff-toolkit rdkit\n'
                'then restart ChimeraX.'),
            'amber14+mmff': (not _mmff_available(),
                'RDKit is not installed.\n'
                'To enable, run in an admin cmd:\n'
                '  python.exe -m pip install rdkit\n'
                'then restart ChimeraX.\n'
                '(approximate — prefer amber14+garnet for production)'),
        }

        for name in all_backends:
            cb.addItem(name)
            disabled_flag, tooltip = _disabled.get(name, (False, ''))
            if disabled_flag:
                model = cb.model()
                item = model.item(cb.count() - 1)
                item.setEnabled(False)
                item.setToolTip(tooltip)

        current = getattr(sim_params, 'forcefield', 'amber14')
        idx = cb.findText(current)
        cb.setCurrentIndex(idx if idx >= 0 else 0)
        cb.blockSignals(False)

    def _combo_changed_cb(self, *_):
        chosen = self.combo_box.currentText()
        from chimerax.isolde.openmm.param_provider import (
            _garnet_available, _espaloma_available, _espaloma_charge_available,
            _openff_available, _mmff_available)
        if chosen == 'amber14+garnet'          and not _garnet_available():          return
        if chosen == 'amber14+espaloma'        and not _espaloma_available():        return
        if chosen == 'amber14+espaloma-charge' and not _espaloma_charge_available(): return
        if chosen == 'amber14+openff'          and not _openff_available():          return
        if chosen == 'amber14+mmff'            and not _mmff_available():            return
        self.session.isolde.sim_params.forcefield = chosen
