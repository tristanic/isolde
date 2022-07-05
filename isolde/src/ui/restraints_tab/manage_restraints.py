from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QLabel, QPushButton, QMenu, QToolBar
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class ManageRestraintsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Manage Adaptive Restraints", **kwargs)
        cd = self.content = ManageRestraintsDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class ManageRestraintsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        dl = DefaultVLayout()
        dbl1 = DefaultHLayout()
        dbl1.addWidget(QLabel('Distance restraint group: '))
        db = self.distance_restraint_group_menu_button = QPushButton()
        db.setText('(None)')
        dm = self.distance_restraint_group_menu = QMenu()
        db.setMenu(dm)
        dm.aboutToShow.connect(self._update_distance_group_menu)
        dbl1.addWidget(db)
        dbl1.addWidget(QLabel(' Color scheme: '))

        from ..color_button import ThreeColorButton
        bsize =(DEFAULT_ICON_SIZE.width(), DEFAULT_ICON_SIZE.height())
        dcb = self.distance_restraint_color_button = ThreeColorButton(
            max_size=bsize, min_size=bsize, direction='0+', has_alpha_channel=False,
            titles=('Overly compressed', 'Optimally satisfied', 'Overly stretched')
        )
        dcb.setToolTip('<span>Set the colour map for the current restraint group</span>')
        dcb.color_changed.connect(self._distance_restraint_color_button_cb)
        dbl1.addWidget(dcb)

        dbl1.addStretch()
        dl.addLayout(dbl1)
        dtb = self.distance_restraint_toolbar = QToolBar()
        dtb.setIconSize(DEFAULT_ICON_SIZE)
        dtb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        dl.addWidget(dtb)
        from chimerax.core.commands import run

        raa = dtb.addAction(
            QIcon(os.path.join(icon_dir, 'all_selected_distance_restraints.png')),
            'Release\nall selected'
            )
        raa.triggered.connect(lambda *_: 
            run(self.session, f'isolde release distances #{self.isolde.selected_model.id_string}&sel groupName "{self.active_distance_restraint_group.name}"'))
        raa.setToolTip('<span>Release all restraints involving at least one of the selected atoms</span>')
        
        roa = dtb.addAction(
            QIcon(os.path.join(icon_dir, 'outside_distance_restraints.png')),
            'Release\nsurrounding'
        )
        roa.triggered.connect(lambda *_: 
            run(self.session, f'isolde release distances #{self.isolde.selected_model.id_string}&sel externalOnly true groupName "{self.active_distance_restraint_group.name}"'))
        roa.setToolTip('<span>Release all restraints involving exactly one selected atom</span>')        

        ria = dtb.addAction(
            QIcon(os.path.join(icon_dir, 'inside_distance_restraints.png')),
            'Release\ninternal'
        )
        ria.triggered.connect(lambda *_: 
            run(self.session, f'isolde release distances #{self.isolde.selected_model.id_string}&sel internalOnly true groupName "{self.active_distance_restraint_group.name}"'))
        ria.setToolTip('<span>Release all restraints where both endpoint atoms are selected</span>')        

        remd = dtb.addAction(
            QIcon(os.path.join(icon_dir, 'red-x-icon.png')),
            'Remove\nall'
        )
        remd.triggered.connect(lambda *_: run(self.session, f'close #{self.active_distance_restraint_group.id_string}'))


        ml.addLayout(dl)
        
        self.container.expanded.connect(self._expanded_cb)
        self._isolde_trigger_handlers.append(
            self.isolde.triggers.add_handler(self.isolde.SELECTED_MODEL_CHANGED,
            self._selected_model_changed_cb)
        )
        from chimerax.core.models import REMOVE_MODELS, ADD_MODELS
        self._chimerax_trigger_handlers.append(
            self.session.triggers.add_handler(REMOVE_MODELS, self._model_add_or_remove_cb)
        )
        self._chimerax_trigger_handlers.append(
            self.session.triggers.add_handler(ADD_MODELS, self._model_add_or_remove_cb)
        )

    def _distance_restraint_color_button_cb(self):
        dcb = self.distance_restraint_color_button
        self.active_distance_restraint_group.set_colormap(*dcb.color)

    def _expanded_cb(self):
        ag = self.active_distance_restraint_group
        if ag is not None and not ag.was_deleted and ag.parent==self.isolde.selected_model:
            return
        mgrs = self._existing_distance_restraint_mgrs()
        if len(mgrs):
            self.active_distance_restraint_group=mgrs[0]
        else:
            self.active_distance_restraint_group=None
    
    def _selected_model_changed_cb(self, trigger_name, model):
        self._expanded_cb()
    
    def _model_add_or_remove_cb(self, *_):
        self._expanded_cb()

    def _existing_distance_restraint_mgrs(self):
        m = self.isolde.selected_model
        if m is None:
            return []
        from chimerax.isolde.molobject import AdaptiveDistanceRestraintMgr
        return [dm for dm in m.child_models() if isinstance(dm, AdaptiveDistanceRestraintMgr)]

    def _set_distance_restraint_tb_state(self):
        dtb = self.distance_restraint_toolbar
        if self.active_distance_restraint_group is None:
            dtb.setEnabled(False)
        else:
            dtb.setEnabled(True)
            
    def _update_distance_group_menu(self, *_):
        dm = self.distance_restraint_group_menu
        db = self.distance_restraint_group_menu_button
        dm.clear()
        mgrs = self._existing_distance_restraint_mgrs()
        if not len(mgrs):
            return
        for mgr in mgrs:
            a = dm.addAction(mgr.name)
            a.triggered.connect(lambda *_, mm=mgr: self._set_active_distance_restraint_group(mm))
    
    def _get_active_distance_restraint_group(self):
        return getattr(self, '_active_distance_restraint_group', None)
    
    def _set_active_distance_restraint_group(self, group):
        self._active_distance_restraint_group = group
        db = self.distance_restraint_group_menu_button
        if group is None:
            db.setText('(None)')
            self.distance_restraint_color_button.setEnabled(False)
        else:
            db.setText(group.name)
            self.distance_restraint_color_button.setEnabled(True)
            from ..util import slot_disconnected
            with slot_disconnected(self.distance_restraint_color_button.color_changed, self._distance_restraint_color_button_cb):
                self.distance_restraint_color_button.color=group.get_colormap()
        self._set_distance_restraint_tb_state()
    
    active_distance_restraint_group = property(_get_active_distance_restraint_group, _set_active_distance_restraint_group)
        


