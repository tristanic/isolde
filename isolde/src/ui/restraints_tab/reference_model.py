from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultVLayout, DefaultHLayout, QLedLabel
from Qt.QtWidgets import (QLabel, QWidget, QPushButton, QMenu, 
    QTreeWidget, QTreeWidgetItem,
    QFileDialog,
)
from Qt.QtGui import QIcon
from Qt.QtCore import Qt
from .. import icon_dir, DEFAULT_ICON_SIZE

class ReferenceModelPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Reference Models", **kwargs)
        cd = self.content = ReferenceModelDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class ReferenceModelDialog(UI_Panel_Base):

    WORKING_CHAIN_COLUMN=0
    REFERENCE_CHAIN_COLUMN=1
    IDENTITY_COLUMN=2
    DIST_COLUMN=3
    TORS_COLUMN=4

    _column_labels = {
        WORKING_CHAIN_COLUMN: 'Chain',
        REFERENCE_CHAIN_COLUMN: 'Reference chain',
        IDENTITY_COLUMN: 'Identity (%)',
        DIST_COLUMN: 'Distances',
        TORS_COLUMN: 'Torsions',
    }
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        import os
        from chimerax.core.commands import run
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        ml = self.main_layout = DefaultVLayout()
        self._current_reference = None

        rml = DefaultHLayout()
        rml.addWidget(QLabel('Reference model: '))
        rmb = self.reference_model_menu_button = QPushButton('(None)')
        rmm = self.reference_model_menu = QMenu()
        rmm.aboutToShow.connect(self._populate_reference_model_menu)
        rmb.setMenu(rmm)
        rml.addWidget(rmb)
        rml.addStretch()
        ab = self.align_button = QPushButton('Align (optional)')
        ab.clicked.connect(lambda _:run(session, f'match #{self._current_reference.id_string} to #{self.isolde.selected_model.id_string}'))
        rml.addWidget(ab)
        ml.addLayout(rml)

        dl = DefaultVLayout()
        sl = DefaultHLayout()
        ll = self.has_pae_indicator = QLedLabel(None)
        ll.setColor('red')
        sl.addWidget(ll)
        pa = self.pae_assigned_label = QLabel('No PAE matrix assigned')
        pb = self.load_pae_button = QPushButton('Load PAE matrix')
        pb.clicked.connect(self._load_pae_from_file)
        sl.addStretch()
        sl.addWidget(pb)
        dl.addLayout(sl)

        dl.addWidget(QLabel('Assign restraints'))
        dt = self.detail_tree = QTreeWidget()
        dl.addWidget(dt)
        ml.addLayout(dl)

        self._initialize_tree()

        self.container.expanded.connect(self._on_expand)

    def _on_expand(self):
        cm = self._current_reference
        if cm is not None and cm.was_deleted():
            self.set_reference_model(None)
    
    def _initialize_tree(self):
        tree = self.detail_tree
        tree.setHeaderLabels(self._column_labels.values())
    
    def _populate_reference_model_menu(self, *_):
        rmm = self.reference_model_menu
        rmm.clear()
        sm = self.isolde.selected_model
        from chimerax.atomic import AtomicStructure
        models = [m for m in self.session.models if type(m) == AtomicStructure]
        if sm is not None:
            models.pop(models.index(sm))
            a = rmm.addAction('(self)')
            a.triggered.connect(lambda _, m=sm: self.set_reference_model(m))
        for m in models:
            a = rmm.addAction(f'{m.id_string}: {m.name}')
            a.triggered.connect(lambda _, mm=m: self.set_reference_model(mm))
            
    def set_reference_model(self, m):
        self._current_reference = m
        rmb = self.reference_model_menu_button
        if m is None:
            rmb.setText('(None')
        elif m == self.isolde.selected_model:
            rmb.setText('(self)')
        else:
            rmb.setText(f'{m.id_string}')
        self._pae_assigned_check()
        self._populate_detail_tree()

    def _pae_assigned_check(self):
        assigned=False
        cr = self._current_reference
        if cr is not None:
            if getattr(cr, 'alphafold_pae', None) is not None:
                assigned=True
        if assigned:
            self.pae_assigned_label.setText('PAE matrix assigned')
            self.has_pae_indicator.setColor('green')
        else:
            self.pae_assigned_label.setText('No PAE matrix assigned')
            self.has_pae_indicator.setColor('red')


    def _populate_detail_tree(self, *_):
        self.detail_tree.clear()
        sm = self.isolde.selected_model
        cm = self._current_reference
        if sm is None or cm is None:
            return
        
    def _load_pae_from_file(self, *_):
        caption = 'Load the PAE matrix for a structure prediction'
        filetypes = 'JSON files (*.json)'
        dlg = QFileDialog(caption=caption)
        dlg.setNameFilter(filetypes)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setFileMode(QFileDialog.ExistingFile)
        selected_file = None
        if dlg.exec():
            selected_file = dlg.selectedFiles()[0]
        if selected_file is not None:
            from chimerax.core.commands import run
            run(self.session,f'alphafold pae #{self._current_reference.id_string} file {selected_file} plot false')
        self._pae_assigned_check()
