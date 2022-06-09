from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QToolBar, QDoubleSpinBox, QLabel, QWidget
from Qt.QtGui import QIcon
from Qt.QtCore import Qt
from .. import icon_dir, DEFAULT_ICON_SIZE

class PositionRestraintsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Position Restraints", **kwargs)
        cd = self.content = PositionRestraintsDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class PositionRestraintsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        tb = self.toolbar = QToolBar()
        tb.setIconSize(DEFAULT_ICON_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        ml.addWidget(tb)
        pa = self.pin_to_current_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'pin-to-current-icon.png')),
            'Pin\nin place'
        )
        pa.setToolTip('<span>Pin all selected non-hydrogen atoms to their current positions</span>')
        pa.triggered.connect(self._pin_selection_in_place)

        pc = self.pin_to_center_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'pin-to-pivot-icon.png')),
            'Pin\n to pivot'
        )
        pc.setToolTip('<span>Pin the selected atom to the current position of the pivot indicator</span>')
        pc.triggered.connect(self._pin_selected_atom_to_pivot)

        ra = self.release_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'red-x-icon.png')),
            'Release'
        )
        ra.setToolTip('<span>Release all position restraints on the selected atoms</span>')
        ra.triggered.connect(self._release_restraints_on_selection)

        tb.addSeparator()
        sbw = QWidget()
        sl = DefaultVLayout()
        sbw.setLayout(sl)
        sl.addWidget(QLabel('Spring constant'))
        ksb = self.spring_constant_spin_box = QDoubleSpinBox()
        ksb.setDecimals(1)
        ksb.setSingleStep(5)
        ksb.setMinimum(0)
        ksb.setMaximum(1000)
        sl.addWidget(ksb)
        sl.addWidget(QLabel("<html><head/><body><p>kJ mol<span style=\" vertical-align:super;\">-1</span> Å<span style=\" vertical-align:super;\">-2</span></p></body></html>"))
        sp = self.isolde.sim_params
        self._param_changed_cb('', ('position_restraint_spring_constant', sp.position_restraint_spring_constant))
        ksb.valueChanged.connect(self._k_spin_box_changed_cb)
        ksb.setKeyboardTracking(False)
        ksb.setToolTip('<span>Spring constant to be applied to newly added restraints (does not affect existing restraints).</span>')
        self._isolde_trigger_handlers.append(sp.triggers.add_handler(
            sp.PARAMETER_CHANGED, self._param_changed_cb
        ))
        tb.addWidget(sbw)

        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(self.session.triggers.add_handler(
            SELECTION_CHANGED, self._selection_changed_cb
        ))
          
    def _pin_selection_in_place(self, *_):
        sel = self.isolde.selected_atoms
        if not len(sel):
            return
        k = self.isolde.sim_params.position_restraint_spring_constant
        from chimerax.isolde.isolde import OPENMM_SPRING_UNIT
        k = k.value_in_unit(OPENMM_SPRING_UNIT)
        from chimerax.isolde import session_extensions as sx
        prm = sx.get_position_restraint_mgr(self.isolde.selected_model)
        prs = prm.add_restraints(sel)
        prs.targets = prs.atoms.coords
        prs.spring_constants = k
        prs.enableds = True
    
    def _pin_selected_atom_to_pivot(self, *_):
        sel = self.isolde.selected_atoms
        if not len(sel):
            return
        if len(sel) > 1:
            from chimerax.core.errors import UserError
            raise UserError('Pinning to the pivot location can only be done for a single atom at a time!')
        sel = sel[0]
        k = self.isolde.sim_params.position_restraint_spring_constant
        from chimerax.isolde.isolde import OPENMM_SPRING_UNIT
        k = k.value_in_unit(OPENMM_SPRING_UNIT)
        from chimerax.isolde import session_extensions as sx
        prm = sx.get_position_restraint_mgr(self.isolde.selected_model)
        pr = prm.add_restraint(sel)
        pr.target = self.session.view.center_of_rotation
        pr.spring_constant = k
        pr.enabled = True

    def _release_restraints_on_selection(self, *_):
        sel = self.isolde.selected_atoms
        if not len(sel):
            return
        from chimerax.isolde import session_extensions as sx
        prm = sx.get_position_restraint_mgr(self.isolde.selected_model)
        prs = prm.get_restraints(sel)
        prs.enableds=False

    
    def _selection_changed_cb(self, *_):
        sel = self.isolde.selected_atoms
        sel = sel[sel.element_names!='H']
        if len(sel):
            self.pin_to_current_action.setEnabled(True)
            self.release_action.setEnabled(True)
        else:
            self.pin_to_current_action.setEnabled(False)
            self.release_action.setEnabled(False)
        if len(sel)==1:
            self.pin_to_center_action.setEnabled(True)
        else:
            self.pin_to_center_action.setEnabled(False)
        
    def _k_spin_box_changed_cb(self, val):
        from chimerax.isolde.isolde import CHIMERAX_SPRING_UNIT, OPENMM_SPRING_UNIT
        self.isolde.sim_params.position_restraint_spring_constant = \
            (val*CHIMERAX_SPRING_UNIT).value_in_unit(OPENMM_SPRING_UNIT)

    
    def _param_changed_cb(self, trigger_name, data):
        param_name, val = data
        if param_name == 'position_restraint_spring_constant':
            from ..util import slot_disconnected
            with slot_disconnected(self.spring_constant_spin_box.valueChanged, self._k_spin_box_changed_cb):
                from chimerax.isolde.isolde import CHIMERAX_SPRING_UNIT
                val = val.value_in_unit(CHIMERAX_SPRING_UNIT)
                self.spring_constant_spin_box.setValue(val)




        