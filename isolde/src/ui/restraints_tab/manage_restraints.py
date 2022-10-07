from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, HorizontalLine
from Qt.QtWidgets import QLabel, QPushButton, QMenu, QToolBar, QSlider, QGridLayout
from Qt.QtGui import QIcon
from Qt.QtCore import Qt, QSize
from .. import icon_dir, DEFAULT_ICON_SIZE

class ManageRestraintsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Manage/Release Adaptive Restraints", **kwargs)
        cd = self.content = ManageRestraintsDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class ManageRestraintsDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        ### DISTANCE RESTRAINTS SECTION
        dl = DefaultVLayout()
        dbl1 = DefaultHLayout()
        dbl1.addWidget(QLabel('Distance restraint group: '))
        db = self.distance_restraint_group_menu_button = QPushButton()
        db.setText('(None)')
        dm = self.distance_restraint_group_menu = QMenu()
        db.setMenu(dm)
        dm.aboutToShow.connect(self._update_distance_group_menu)
        dbl1.addWidget(db)
        dl.addLayout(dbl1)
        dchl = DefaultHLayout()
        dchl.addWidget(QLabel(' Colors: '))

        from ..color_button import ThreeColorButton
        bsize =(DEFAULT_ICON_SIZE.width(), DEFAULT_ICON_SIZE.height())
        dcb = self.distance_restraint_color_button = ThreeColorButton(
            max_size=bsize, direction='0+', has_alpha_channel=False,
            titles=('Overly compressed', 'Optimally satisfied', 'Overly stretched')
        )
        dcb.setToolTip('<span>Set the colour map for the current restraint group</span>')
        dcb.color_changed.connect(self._distance_restraint_color_button_cb)
        dchl.addWidget(dcb)
        ddb = self.distance_default_color_button = QPushButton('Reset')
        ddb.setToolTip('<span>Reset colour scheme to default</span>')
        def set_distance_colors_to_default():
            from ..util import slot_disconnected
            group = self.active_distance_restraint_group
            group.set_default_colormap()
            with slot_disconnected(dcb.color_changed, self._distance_restraint_color_button_cb):
                dcb.color=group.get_colormap()


        ddb.clicked.connect(set_distance_colors_to_default)
        dchl.addWidget(ddb)
        dchl.addStretch()
        dl.addLayout(dchl)

        dbl2 = QGridLayout()
        dbl2.addWidget(QLabel('Display threshold'),0,0,1,1)
        ddsl = self.distance_restraint_display_threshold_slider = QSlider(Qt.Orientation.Horizontal)
        ddsl.setMinimum(0)
        ddsl.setMaximum(100)
        ddsl.valueChanged.connect(self._distance_display_slider_value_changed_cb)
        dbl2.addWidget(ddsl, 1,0,1,3)
        dbl2.addWidget(QLabel('All'), 2,0)
        sl = QLabel('Severely strained')
        sl.setAlignment(Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignRight)
        dbl2.addWidget(sl,2,2)
        dl.addLayout(dbl2)

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
        remd.setToolTip('<span>Completely remove (delete) this restraint group</span>')
        remd.triggered.connect(lambda *_: run(self.session, f'close #{self.active_distance_restraint_group.id_string}'))


        ml.addLayout(dl)
        ml.addWidget(HorizontalLine())
        ### TORSION RESTRAINTS SECTION
        tl = DefaultVLayout()
        tbl1 = DefaultHLayout()
        tl.addLayout(tbl1)
        tbl1.addWidget(QLabel('Adaptive Torsion Restraints  '))
        tbl1.addWidget(QLabel('Colors: '))
        tcb = self.torsion_restraint_color_button = ThreeColorButton(
            max_size=bsize, direction='0-', has_alpha_channel=False,
            titles=('Severely strained','Moderately strained','Satisfied')
        )
        tcb.color_changed.connect(self._torsion_restraint_button_color_cb)
        tbl1.addWidget(tcb)
        tdb = self.torsion_default_color_button = QPushButton('Reset')
        def set_default_torsion_colors():
            self.torsion_restraint_mgr.set_default_colors()
            self._update_torsion_restraint_button()
        tdb.clicked.connect(set_default_torsion_colors)
        tbl1.addWidget(tdb)
        tbl1.addStretch()

        ttb = self.torsion_restraints_toolbar = QToolBar()
        ttb.setIconSize(DEFAULT_ICON_SIZE)
        ttb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)

        rata = ttb.addAction(
            QIcon(os.path.join(icon_dir, 'release_all_torsions.png')),
            'Release\nall selected'
        )
        rata.setToolTip('<span>Release all torsion restraints on the selected residues</span>')
        rata.triggered.connect(lambda *_:
            run(self.session, f'isolde release torsions #{self.isolde.selected_model.id_string}&sel'))
        
        rsta = ttb.addAction(
            QIcon(os.path.join(icon_dir, 'release_sidechain_torsions.png')),
            'Release\nselected sidechains'
        )
        rsta.setToolTip('<span>Release sidechain torsion restraints on the selected residues</span>')
        rsta.triggered.connect(lambda *_:
        run(self.session, f'isolde release torsions #{self.isolde.selected_model.id_string} backbone f'))
        remt = ttb.addAction(
            QIcon(os.path.join(icon_dir, 'red-x-icon.png')),
            'Remove\nall'
        )
        remt.setToolTip('<span>Completely remove all restraints</span>')
        remt.triggered.connect(lambda *_: run(self.session, f'close #{self.torsion_restraint_mgr.id_string}'))
        tl.addWidget(ttb)
        ml.addLayout(tl)
        
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

    def _distance_display_slider_value_changed_cb(self, value):
        from chimerax.core.commands import run
        run(self.session, f'isolde adj dist #{self.isolde.selected_model.id_string} displayThreshold {value/50} groupName "{self.active_distance_restraint_group.name}"', log=False)


    def _distance_restraint_color_button_cb(self, *_):
        dcb = self.distance_restraint_color_button
        self.active_distance_restraint_group.set_colormap(*dcb.color)

    def _expanded_cb(self):
        tm = self.torsion_restraint_mgr
        if tm is None:
            self.torsion_restraints_toolbar.setEnabled(False)            
        else:
            self.torsion_restraints_toolbar.setEnabled(True)
        self._update_torsion_restraint_button()
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

    def _update_torsion_restraint_button(self):
        btn = self.torsion_restraint_color_button
        tmgr = self.torsion_restraint_mgr
        if tmgr is None:
            btn.setEnabled(False)
            self.torsion_default_color_button.setEnabled(False)
            return
        btn.setEnabled(True)
        self.torsion_default_color_button.setEnabled(True)
        from ..util import slot_disconnected
        if tmgr is not None:
            with slot_disconnected(btn.color_changed, self._torsion_restraint_button_color_cb):
                btn.color=tmgr.get_color_scale()
                
    
    def _torsion_restraint_button_color_cb(self, *_):
        tmgr = self.torsion_restraint_mgr
        tmgr.set_color_scale(*self.torsion_restraint_color_button.color)

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
            self.distance_default_color_button.setEnabled(False)
            self.distance_restraint_display_threshold_slider.setEnabled(False)
        else:
            db.setText(group.name)
            drcb = self.distance_restraint_color_button
            drcb.setEnabled(True)            
            from ..util import slot_disconnected
            with slot_disconnected(drcb.color_changed, self._distance_restraint_color_button_cb):
                drcb.color=group.get_colormap()
            ddsl = self.distance_restraint_display_threshold_slider
            ddsl.setEnabled(True)
            with slot_disconnected(ddsl.valueChanged, self._distance_display_slider_value_changed_cb):
                ddsl.setValue(min(max(int(self.active_distance_restraint_group.display_threshold*50),0),100))
            self.distance_default_color_button.setEnabled(True)
        self._set_distance_restraint_tb_state()
    
    active_distance_restraint_group = property(_get_active_distance_restraint_group, _set_active_distance_restraint_group)

    @property
    def torsion_restraint_mgr(self):
        if self.isolde.selected_model is None:
            return None
        from chimerax.isolde import session_extensions as sx
        return sx.get_adaptive_dihedral_restraint_mgr(self.isolde.selected_model)    


