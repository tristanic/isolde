from ..ui_base import (
    UI_Panel_Base, DefaultHLayout, DefaultVLayout, DefaultSpacerItem,
    QDoubleSpinBox
)
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QPushButton, QLabel, QMenu, QFileDialog

class MapAddPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, 'Add map(s) to working model', **kwargs)
        mad = self.content = MapAddDialog(session, isolde, gui, self.content_area)
        self.setContentLayout(mad.main_layout)

class MapAddDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, main_frame):
        super().__init__(session, isolde, gui, main_frame)
        self.selected_model_changed_cb('', self.isolde.selected_model)
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        rsmb = self.real_space_map_menu_button = QPushButton(mf)
        rsmb.setText('From loaded volume')
        rsmb.setToolTip('Assign an existing loaded map to the model. The map must be in an '
            'absolute reference frame (that is: move the model to fit the map, rather than '
            'moving the map to fit the model.')
        rsmm = self.real_space_map_menu = QMenu(rsmb)
        rsmm.aboutToShow.connect(self._populate_available_volumes_menu)
        rsmb.setMenu(rsmm)
        ml.addWidget(rsmb)

        hl = DefaultHLayout()

        cxmb = self.crystal_dataset_button = QPushButton(mf)
        cxmb.setText('From crystallographic dataset')
        cxmb.setToolTip('Load a crystallographic dataset or precalculated map structure factors in '
            '.mtz or .cif format. If your model has a space group and cell dimensions '
            'assigned they must match those in the loaded file.')
        cxmb.clicked.connect(self.load_structure_factors)
        hl.addWidget(cxmb)

        hl.addItem(DefaultSpacerItem())
        hl.addWidget(QLabel('Oversampling: ', parent=mf))
        ossb = self.crystal_oversample_spinbox = QDoubleSpinBox(mf)
        ossb.setMinimum(1.5)
        ossb.setMaximum(5)
        ossb.setDecimals(1)
        ossb.setValue(self.isolde.params.map_shannon_rate)
        ossb.valueChanged.connect(self._sampling_rate_spinbox_cb)
        ossb.setSingleStep(0.1)
        hl.addWidget(ossb)
        ml.addLayout(hl)

    def _populate_available_volumes_menu(self):
        rsmm = self.real_space_map_menu
        from chimerax.map import Volume
        # Only use explicit Volume instances, not subclasses
        volumes = [v for v in self.session.models.list() if type(v)==Volume]
        rsmm.clear()
        for v in volumes:
            action = rsmm.addAction(f'#{v.id_string}: {v.name}')
            def add_map_cb(*_, session=self.session, volume=v):
                from chimerax.core.commands import run
                run(session, f'clipper assoc #{volume.id_string} to #{session.isolde.selected_model.parent.id_string}')
            action.triggered.connect(add_map_cb)
        

    def selected_model_changed_cb(self, _, selected_model):
        if selected_model is None:
            self.main_frame.setEnabled(False)
        else:
            self.main_frame.setEnabled(True)

    def _sampling_rate_spinbox_cb(self, value):
        self.isolde.params.map_shannon_rate = value

    def load_structure_factors(self, *_):
        caption = 'Choose a file containing crystallographic structure factors'
        filetypes = "Structure factor files (*.mtz *.cif)"
        dlg = QFileDialog(caption=caption)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilter(filetypes)
        dlg.setFileMode(QFileDialog.ExistingFile)
        selected_file = None
        if dlg.exec():
            selected_file = dlg.selectedFiles()[0]
        if selected_file is not None:
            from chimerax.core.commands import run
            oversamp = self.crystal_oversample_spinbox.value()
            run(self.session, f'clipper open "{selected_file}" structureModel #{self.isolde.selected_model.id_string} overSampling {oversamp:.1f}')
            
        
