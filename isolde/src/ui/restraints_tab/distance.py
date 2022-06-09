from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultVLayout, QDoubleSpinBox
from Qt.QtWidgets import QToolBar, QLabel, QWidget, QSizePolicy
from Qt.QtGui import QIcon
from Qt.QtCore import Qt
from .. import icon_dir, DEFAULT_ICON_SIZE

class DistanceRestraintsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Simple Distance Restraints", **kwargs)
        cd = self.content = DistanceRestraintsDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class DistanceRestraintsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        tb = self.toolbar = QToolBar()
        tb.setIconSize(DEFAULT_ICON_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        ml.addWidget(tb)

        ra = self.restrain_distance_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'distance-restraint-icon.png')),
            'Add/adjust\nrestraint'
        )
        ra.setToolTip('<span>Add a distance restraint between two selected atoms, or adjust an existing one to the current distance and spring constant.</span>')
        ra.triggered.connect(self._add_or_adjust_distance_restraint)

        rsb = self.release_specific_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'distance_restraint_remove_internal.png')),
            'Release\n(internal)'
        )
        rsb.setToolTip('<span>Release all restraints where *both* atoms are in the current selection</span>')
        rsb.triggered.connect(self._release_intra)

        rab = self.release_permissive_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'distance_restraint_remove_all.png')),
            'Release\n(any)'
        )
        rab.setToolTip('<span>Release all restraints where *either* atom is in the current selection</span>')
        rab.triggered.connect(self._release_all)

        tb.addSeparator()

        kw = QWidget()
        kl = DefaultVLayout()
        kw.setLayout(kl)
        kl.addWidget(QLabel('Spring constant'))
        ksb = self.spring_constant_spin_box = QDoubleSpinBox()
        ksb.setDecimals(1)
        ksb.setSingleStep(5)
        ksb.setMinimum(0)
        ksb.setMaximum(1000)
        ksb.setKeyboardTracking(False)
        kl.addWidget(ksb)
        kl.addWidget(QLabel("<html><head/><body><p>kJ mol<span style=\" vertical-align:super;\">-1</span> Å<span style=\" vertical-align:super;\">-2</span></p></body></html>"))
        sp = self.isolde.sim_params
        from chimerax.isolde.isolde import CHIMERAX_SPRING_UNIT
        ksb.setValue(sp.distance_restraint_spring_constant.value_in_unit(CHIMERAX_SPRING_UNIT))
        ksb.valueChanged.connect(self._k_spin_box_changed_cb)
        self._isolde_trigger_handlers.append(sp.triggers.add_handler(sp.PARAMETER_CHANGED, self._param_changed_cb))
        tb.addWidget(kw)

        # Just to provide a bit of a buffer between spring constant and distance boxes
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)
        
        dw = QWidget()
        dl = DefaultVLayout()
        dw.setLayout(dl)
        dl.addWidget(QLabel('Target'))
        dsb = self.distance_spin_box = QDoubleSpinBox()
        dsb.setDecimals(2)
        dsb.setSingleStep(0.1)
        dsb.setMinimum(0)
        dsb.setValue(3)
        dl.addWidget(dsb)
        dl.addWidget(QLabel('Å'))
        tb.addWidget(dw)

        sa = self.seed_distance_action = tb.addAction(
            QIcon(os.path.join(icon_dir, 'atom_atom_distance.png')),
            'Use current\ndistance'
        )
        sa.setToolTip('<span>Set the target value to the current interatomic distance</span>')
        sa.triggered.connect(self._seed_distance)
        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(self.session.triggers.add_handler(
            SELECTION_CHANGED, self._selection_changed_cb
        ))
        self._selection_changed_cb()


    def _add_or_adjust_distance_restraint(self, *_):
        sel = self.isolde.selected_atoms
        if len(sel) != 2:
            return
        from chimerax.isolde import session_extensions as sx
        drm = sx.get_distance_restraint_mgr(self.isolde.selected_model)
        r = drm.add_restraint(*sel)
        r.target = self.distance_spin_box.value()
        from chimerax.isolde.isolde import OPENMM_SPRING_UNIT
        r.spring_constant = self.isolde.sim_params.distance_restraint_spring_constant.value_in_unit(OPENMM_SPRING_UNIT)
        r.enabled = True
    
    def _release_intra(self, *_):
        sel = self.isolde.selected_atoms
        if not len(sel):
            return
        from chimerax.isolde import session_extensions as sx
        drm = sx.get_distance_restraint_mgr(self.isolde.selected_model)
        drm.intra_restraints(sel).enableds = False

    def _release_all(self, *_):
        sel = self.isolde.selected_atoms
        if not len(sel):
            return
        from chimerax.isolde import session_extensions as sx
        drm = sx.get_distance_restraint_mgr(self.isolde.selected_model)
        drm.atoms_restraints(sel).enableds = False

    def _seed_distance(self, *_):
        sel = self.isolde.selected_atoms
        sel = sel[sel.element_names !='H']
        if len(sel) != 2:
            return
        from chimerax.geometry import distance
        self.distance_spin_box.setValue(distance(*sel.coords))

    def _k_spin_box_changed_cb(self, val):
        from chimerax.isolde.isolde import CHIMERAX_SPRING_UNIT, OPENMM_SPRING_UNIT
        self.isolde.sim_params.distance_restraint_spring_constant = \
            (val*CHIMERAX_SPRING_UNIT).value_in_unit(OPENMM_SPRING_UNIT)

    def _param_changed_cb(self, trigger_name, data):
        param_name, val = data
        if param_name == 'distance_restraint_spring_constant':
            from ..util import slot_disconnected
            with slot_disconnected(self.spring_constant_spin_box.valueChanged, self._k_spin_box_changed_cb):
                from chimerax.isolde.isolde import CHIMERAX_SPRING_UNIT
                val = val.value_in_unit(CHIMERAX_SPRING_UNIT)
                self.spring_constant_spin_box.setValue(val)

    def _selection_changed_cb(self, *_):
        sel = self.isolde.selected_atoms
        sel = sel[sel.element_names!='H']
        if len(sel):
            self.release_specific_action.setEnabled(True)
            self.release_permissive_action.setEnabled(True)
        else:
            self.release_specific_action.setEnabled(False)
            self.release_permissive_action.setEnabled(False)
        if len(sel)==2:
            self.restrain_distance_action.setEnabled(True)
            self.seed_distance_action.setEnabled(True)
        else:
            self.restrain_distance_action.setEnabled(False)
            self.seed_distance_action.setEnabled(False)

