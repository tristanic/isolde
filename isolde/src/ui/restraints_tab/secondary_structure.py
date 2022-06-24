from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, QSpinBox
from Qt.QtWidgets import QToolBar
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class SecondaryStructurePanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Secondary Structure", **kwargs)
        cd = self.content = SecondaryStructureDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class SecondaryStructureDialog(UI_Panel_Base):
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
        bs = self.beta_strand_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'antiparallel-beta-icon.png')),
            'Beta\nstrand/sheet'
        )
        bs.setToolTip('<span>Restrain all residues in the selection to beta strand geometry. '
            'Existing inter-strand H-bonds will also be restrained.</span>')
        bs.triggered.connect(lambda *_:self._restrain_single_type('strand'))
        ah = self.alpha_helix_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'helix-icon.png')),
            'Alpha\nhelix'
        )
        ah.setToolTip('<span>Restrain all residues in the selection to alpha helix</span>')
        ah.triggered.connect(lambda *_: self._restrain_single_type('helix'))
        mx = self.mixed_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'mixed-ss-icon.png')),
            'Current'
        )
        mx.setToolTip('<span>Restrain all helix/strand residues in the selection to their '
        'currently defined secondary structure</span>')
        mx.triggered.connect(self._restrain_to_current)
        rel = self.release_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'red-x-icon.png')),
            'Release'
        )
        rel.setToolTip('<span>Release all secondary structure restraints on the selection</span>')
        rel.triggered.connect(self._release_ss_restraints)

        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(
            self.session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed_cb)
        )
        self._selection_changed_cb()

    def _restrain_single_type(self, ss_type):
        sel = self.isolde.selected_atoms.unique_residues
        from chimerax.isolde.restraints.restraint_utils import restrain_secondary_structure
        restrain_secondary_structure(self.session, sel, ss_type)

    def _restrain_to_current(self, *_):
        sel = self.isolde.selected_atoms.unique_residues
        from chimerax.atomic import Residue
        strand = sel[sel.ss_types==Residue.SS_STRAND]
        helix = sel[sel.ss_types==Residue.SS_HELIX]
        from chimerax.isolde.restraints.restraint_utils import restrain_secondary_structure
        restrain_secondary_structure(self.session, strand, 'strand')
        restrain_secondary_structure(self.session, helix, 'helix')
    
    def _release_ss_restraints(self, *_):
        sel = self.isolde.selected_atoms.unique_residues
        from chimerax.isolde.restraints.restraint_utils import release_ss_restraints
        release_ss_restraints(sel)

    def _selection_changed_cb(self, *_):
        sel = self.isolde.selected_atoms.unique_residues
        import numpy
        from chimerax.atomic import Residue
        # The restrain buttons only do anything when at least 3 residues in the selection are 
        # connected... but that would get a bit expensive to detect. Keep it simple.
        if numpy.any(sel.polymer_types==Residue.PT_AMINO):
            self.toolbar.setEnabled(True)
        else:
            self.toolbar.setEnabled(False)
