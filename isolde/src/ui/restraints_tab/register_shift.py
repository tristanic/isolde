from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, QSpinBox
from Qt.QtWidgets import QLabel, QRadioButton, QToolButton
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class RegisterShiftPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Register Shift Protein", **kwargs)
        cd = self.content = RegisterShiftDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class RegisterShiftDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultHLayout()
        ml.addWidget(QLabel('Shift by '))
        nsb = self.num_residues_spin_box = QSpinBox()
        nsb.setMinimum(1)
        nsb.setMaximum(100)
        nsb.setValue(1)
        ml.addWidget(nsb)
        rlabel = QLabel(' residues towards ')
        def update_label(val):
            if val==1:
                label = ' residue towards '
            else:
                label = ' residues towards '
        nsb.valueChanged.connect(update_label)
        ml.addWidget(rlabel)
        blayout = DefaultVLayout()
        ntb = self.towards_n_terminus = QRadioButton('N')
        ctb = self.towards_c_terminus = QRadioButton('C')
        ntb.setChecked(True)
        blayout.addWidget(ntb)
        blayout.addWidget(ctb)
        ml.addLayout(blayout)
        ml.addWidget(QLabel(' terminus '))
        gb = self.go_button = QToolButton()
        gb.setIcon(QIcon(os.path.join(icon_dir, 'registershift-icon.png')))
        gb.setIconSize(QSize(48,36))
        gb.setText('Go')
        gb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        gb.clicked.connect(self._start_register_shift)

        sb = self.stop_button = QToolButton()
        sb.setIcon(QIcon(os.path.join(icon_dir, 'red-x-icon.png')))
        sb.setIconSize(DEFAULT_ICON_SIZE)
        sb.setText('Release')
        sb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        sb.clicked.connect(self._stop_button_cb)
        sb.setEnabled(False)

        ml.addWidget(gb)
        ml.addWidget(sb)
        ml.addStretch()

        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(self.session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed_cb))

    def _start_register_shift(self, *_):
        atoms = self.isolde.selected_atoms
        if not len(atoms):
            return
        if self.towards_c_terminus.isChecked():
            direction = 1
        else:
            direction = -1
        n = self.num_residues_spin_box.value() * direction
        from chimerax.isolde.manipulations.register_shift import ProteinRegisterShifter
        shifter = self._shifter = ProteinRegisterShifter(self.session, self.isolde, atoms)
        from chimerax.isolde.selections import extend_selection_along_chains
        residues = atoms.unique_residues
        for i in range(abs(n)):
            residues = extend_selection_along_chains(residues, direction)
        if not self.isolde.simulation_running:
            from chimerax.core.commands import run
            run(self.session, 'isolde sim start sel')
        else:
            mr = self.isolde.sim_manager.sim_construct.mobile_residues
            import numpy
            from chimerax.atomic import selected_residues
            if numpy.any(mr.indices(selected_residues(self.session))==-1):
                from chimerax.core.errors import UserError
                raise UserError('The mobile component of the simulation must include all residues to be shifted, '
                'plus sufficient residues before and after to accommodate the shift')
        self.stop_button.setEnabled(True)
        self.go_button.setEnabled(False)
        shifter.shift_register(n)
    
    def _stop_button_cb(self, *_):
        if getattr(self, '_shifter', None) is None:
            return
        self._shifter.release_all()
        self._shifter = None
        self._selection_changed_cb()

    def _selection_changed_cb(self, *_):
        if getattr(self, '_shifter', None) is not None:
            return
        sel = self.isolde.selected_atoms
        from chimerax.isolde.util import is_continuous_protein_chain
        if is_continuous_protein_chain(sel) and len(sel.unique_residues) > 3:
            self.go_button.setEnabled(True)
        else:
            self.go_button.setEnabled(False)




