from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, QSpinBox
from Qt.QtWidgets import QToolBar
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class BasePairingPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Base Pairing", **kwargs)
        cd = self.content = BasePairingDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class BasePairingDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()
        tb = self.toolbar = QToolBar()
        ml.addWidget(tb)
        tb.setMovable(False)
        tb.setIconSize(DEFAULT_ICON_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        from chimerax.core.commands import run
        bp = self.base_pair_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'basepair.png')),
            'Restrain\nbase H-bonds'
        )
        bp.setToolTip('<span>Restrain base-pairing H-bonds for all nucleic acid residues in the selection. '
            'Restraints will only be applied to residues already forming base pairs.</span>')
        bp.triggered.connect(lambda *_:run(session, 'isolde restrain basepair sel'))

        rel = self.release_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'red-x-icon.png')),
            'Release'
        )
        rel.setToolTip('<span>Release all base-pair restraints on the selection</span>')
        rel.triggered.connect(self._release_base_pair_restraints)

        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(
            self.session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed_cb)
        )
        self._selection_changed_cb()
    
    def _release_base_pair_restraints(self, *_):
        sel = self.isolde.selected_atoms.unique_residues
        from chimerax.isolde.restraints.restraint_utils import release_base_pair_restraints
        release_base_pair_restraints(sel)

    def _selection_changed_cb(self, *_):
        sel = self.isolde.selected_atoms.unique_residues
        import numpy
        from chimerax.atomic import Residue
        # The restrain buttons only do anything when at least 3 residues in the selection are 
        # connected... but that would get a bit expensive to detect. Keep it simple.
        if numpy.any(sel.polymer_types==Residue.PT_NUCLEIC):
            self.toolbar.setEnabled(True)
        else:
            self.toolbar.setEnabled(False)
