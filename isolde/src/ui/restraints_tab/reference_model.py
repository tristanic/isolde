from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultVLayout, DefaultHLayout, QLedLabel
from Qt.QtWidgets import (QLabel, QWidget, QPushButton, QMenu, 
    QTreeWidget, QTreeWidgetItem,
    QFileDialog,
    QCheckBox, QButtonGroup
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

    IDENTITY_CUTOFF = 0.5

    WORKING_CHAIN_COLUMN=0
    REFERENCE_CHAIN_COLUMN=1
    IDENTITY_COLUMN=2
    ALIGNMENT_LENGTH_COLUMN=3
    DIST_COLUMN=4
    TORS_COLUMN=5

    _column_labels = {
        WORKING_CHAIN_COLUMN: 'Chain',
        REFERENCE_CHAIN_COLUMN: 'Ref chain',
        IDENTITY_COLUMN: 'Identity (%)',
        ALIGNMENT_LENGTH_COLUMN: 'Coverage (%)',
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
        self._restrain_distance_groups = []
        self._restrain_torsions_groups = []

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
            if getattr(cr, 'alphafold_pae', None) is None:
                self._fetch_pae_from_database_if_possible(cr)
            if getattr(cr, 'alphafold_pae', None) is None:
                assigned=True
        if assigned:
            self.pae_assigned_label.setText('PAE matrix assigned')
            self.has_pae_indicator.setColor('green')
        else:
            self.pae_assigned_label.setText('No PAE matrix assigned')
            self.has_pae_indicator.setColor('red')

    def _fetch_pae_from_database_if_possible(self, model):
        from chimerax.alphafold.database import uniprot_id_from_filename
        structure_path = getattr(model, 'filename', None)
        if structure_path is not None:
            uniprot_id = uniprot_id_from_filename(structure_path)
            if uniprot_id is not None:
                self.session.logger.info(f'ISOLDE: attempting to load PAE matrix for model #{model.id_string} from the AlphaFold-EBI '
                    'database. If this fails or you wish to use a different PAE matrix, use the '
                    '"Load PAE matrix" button in ISOLDE\'s reference model restraints widget.')
                from chimerax.core.commands import run
                run(self.session, f'alphafold pae #{model.id_string} uniprot {uniprot_id}')




    def _populate_detail_tree(self, *_):
        self.detail_tree.clear()
        sm = self.isolde.selected_model
        cm = self._current_reference
        self._restrain_distance_groups.clear()
        self._restrain_torsions_groups.clear()
        if sm is None or cm is None:
            return
        from chimerax.atomic import Residue
        model_chains = sm.chains[sm.chains.polymer_types == Residue.PT_AMINO]
        ref_chains = cm.chains[cm.chains.polymer_types == Residue.PT_AMINO]
        
        from chimerax.match_maker.match import defaults, align
        from collections import defaultdict
        alignments = defaultdict(list)
        dssp_cache = {}
        for mc in model_chains:
            for rc in ref_chains:
                score, s1, s2 = align(self.session, mc, rc, defaults['matrix'],
                'nw', defaults['gap_open'], defaults['gap_extend'], dssp_cache
                )
                num_identical=0
                for i, (c1, c2) in enumerate(zip(s1.characters, s2.characters)):
                    if c1==c2:
                        num_identical += 1
                identity = num_identical/i
                coverage = i/len([r for r in mc.residues if r is not None])
                # ref chain, alignment score, fractional identity, coverage, model seq, ref seq
                if identity > self.IDENTITY_CUTOFF:
                    alignments[mc].append((rc, score, identity, coverage, s1, s2))
        # Sort in descending order of score
        for key, rlist in alignments.items():
            alignments[key] = list(sorted(rlist, key=lambda i:i[1], reverse=True))
        
        tree = self.detail_tree
        root = tree.invisibleRootItem()
        for model_chain, refs in alignments.items():
            parent = QTreeWidgetItem(root)
            parent.setText(self.WORKING_CHAIN_COLUMN, f'{model_chain.chain_id}')
            distance_group = QButtonGroup(tree)
            distance_group._model_chain = model_chain
            self._restrain_distance_groups.append(distance_group)
            torsion_group = QButtonGroup(tree)
            torsion_group._model_chain = model_chain
            self._restrain_torsions_groups.append(torsion_group)
            for ref_data in refs:
                rc, score, identity, alignment_length, model_seq, ref_seq = ref_data
                entry = QTreeWidgetItem(parent)
                entry.setText(self.REFERENCE_CHAIN_COLUMN, f'{rc.chain_id}')
                entry.setText(self.IDENTITY_COLUMN, f'{int(round(identity*100, 0))}')
                entry.setText(self.ALIGNMENT_LENGTH_COLUMN, f'{int(round(coverage*100, 0))}')
                w1 = QWidget()
                l = DefaultHLayout()
                w1.setLayout(l)
                l.addStretch()
                dcb = QCheckBox()
                dcb._ref_chain = rc
                distance_group.addButton(dcb)
                l.addWidget(dcb)
                l.addStretch()
                tree.setItemWidget(entry, self.DIST_COLUMN, w1)

                w2 = QWidget()
                l = DefaultHLayout()
                w2.setLayout(l)
                l.addStretch()
                tcb = QCheckBox()
                tcb._ref_chain = rc
                torsion_group.addButton(tcb)
                l.addWidget(tcb)
                l.addStretch()
                tree.setItemWidget(entry, self.TORS_COLUMN, w2)     






        

        
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
