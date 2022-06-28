from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultVLayout, DefaultHLayout, QLedLabel
from Qt.QtWidgets import (QLabel, QWidget, QFrame, QGridLayout, QPushButton, QMenu, 
    QTreeWidget, QTreeWidgetItem,
    QFileDialog,
    QCheckBox, QButtonGroup,
    QSlider
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
    ALIGN_COLUMN=2
    DIST_COLUMN=3
    TORS_COLUMN=4
    IDENTITY_COLUMN=5
    ALIGNMENT_LENGTH_COLUMN=6
    RMSD_COLUMN=7

    _column_labels = {
        WORKING_CHAIN_COLUMN: 'Chain',
        REFERENCE_CHAIN_COLUMN: 'Ref',
        IDENTITY_COLUMN: 'Identity',
        ALIGNMENT_LENGTH_COLUMN: 'Coverage',
        RMSD_COLUMN: 'CA RMSD',
        ALIGN_COLUMN: '',
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
        # ab = self.align_button = QPushButton('Align')
        # ab.clicked.connect(lambda _:run(session, f'match #{self._current_reference.id_string} to #{self.isolde.selected_model.id_string}'))
        # rml.addWidget(ab)
        ll = self.has_pae_indicator = QLedLabel(None)
        ll.setColor('red')
        rml.addWidget(ll)
        pa = self.pae_assigned_label = QLabel('No PAE matrix assigned')
        rml.addWidget(pa)
        pb = self.load_pae_button = QPushButton('Load PAE matrix')
        pb.clicked.connect(self._load_pae_from_file)
        rml.addStretch()
        rml.addWidget(pb)

        ml.addLayout(rml)

        dl = DefaultVLayout()

        dl.addWidget(QLabel('Assign restraints'))
        dt = self.detail_tree = QTreeWidget()
        dl.addWidget(dt)
        ml.addLayout(dl)

        self._initialize_tree()

        options_frame = QFrame()
        options_frame.setContentsMargins(3,3,3,3)
        options_frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        ol = DefaultHLayout()
        opt = self.options = ReferenceModelOptions(gui, options_frame)
        ol.addWidget(opt)
        options_frame.setLayout(ol)
        ml.addWidget(options_frame)



        self.container.expanded.connect(self._on_expand)

    def _on_expand(self):
        cm = self._current_reference
        if cm is not None and cm.was_deleted:
            self.set_reference_model(None)
    
    def _initialize_tree(self):
        tree = self.detail_tree
        items = sorted(self._column_labels.items(), key=lambda i: i[0])
        tree.setHeaderLabels([item[1] for item in items])
        tree.header().setMinimumSectionSize(20)

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
        sm = self.isolde.selected_model
        if m is not None and sm is not None:
            from chimerax.core.commands import run
            run(self.session, f'match #{m.id_string} to #{sm.id_string}')
            run(self.session, f'color #{m.id_string} bychain')
            run(self.session, f'color #{m.id_string} byhet')
            run(self.session, f'color modify #{m.id_string} hue + 50')


        self._pae_assigned_check()
        self._populate_detail_tree()

    def _pae_assigned_check(self):
        assigned=False
        cr = self._current_reference
        if cr is not None:
            if getattr(cr, 'alphafold_pae', None) is None:
                self._fetch_pae_from_database_if_possible(cr)
            if getattr(cr, 'alphafold_pae', None) is not None:
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
                for c1, c2 in zip(s1.characters, s2.characters):
                    if c1==c2:
                        num_identical += 1
                identity = num_identical/len(s1.ungapped())
                coverage = len(s1.ungapped())/len(mc.ungapped())
                print(f's1: {s1.characters}\ns2: {s2.characters}\nmc: {mc.ungapped()}\nidentity: {identity} coverage: {coverage}')
                # ref chain, alignment score, fractional identity, coverage, model seq, ref seq
                if identity > self.IDENTITY_CUTOFF:
                    rmsd = _ca_rmsd(mc.residues, s2.residues)
                    alignments[mc].append((rc, score, identity, coverage, s1, s2, rmsd))
        # Sort in descending order of score
        for key, rlist in alignments.items():
            # Sort by aligned CA-RMSD
            alignments[key] = list(sorted(rlist, key=lambda i:i[-1]))
        
        tree = self.detail_tree
        root = tree.invisibleRootItem()
        for model_chain, refs in alignments.items():
            parent = QTreeWidgetItem(root)
            parent.setText(self.WORKING_CHAIN_COLUMN, f'{model_chain.chain_id}')
            distance_group = ExclusiveOrNoneQButtonGroup(tree)
            distance_group._model_chain = model_chain
            self._restrain_distance_groups.append(distance_group)
            torsion_group = ExclusiveOrNoneQButtonGroup(tree)
            torsion_group._model_chain = model_chain
            self._restrain_torsions_groups.append(torsion_group)
            
            for ref_data in refs:
                rc, score, identity, coverage, model_seq, ref_seq, rmsd = ref_data
                entry = QTreeWidgetItem(parent)
                entry.setText(self.REFERENCE_CHAIN_COLUMN, f'{rc.chain_id}')
                entry.setText(self.IDENTITY_COLUMN, f'{round(identity*100, 0):.0f}%')
                entry.setText(self.ALIGNMENT_LENGTH_COLUMN, f'{round(coverage*100, 0):.0f}%')
                entry.setText(self.RMSD_COLUMN,f'{rmsd:.1f}')
                align_button = QPushButton('Align')
                def align_on_chain(*_,mc=model_chain, rc=rc):
                    from chimerax.core.commands import run
                    run(self.session, f'match #{cm.id_string}/{rc.chain_id} to #{sm.id_string}/{mc.chain_id}')
                    self._populate_detail_tree()
                align_button.clicked.connect(align_on_chain)
                tree.setItemWidget(entry, self.ALIGN_COLUMN, align_button)

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
                    
        tree.expandAll()
        for column in self._column_labels.keys():
            tree.resizeColumnToContents(column)
         






        

        
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


class ReferenceModelOptions(CollapsibleArea):
    def __init__(self, gui, parent, **kwargs):
        super().__init__(gui, parent, title="Options", **kwargs)

        layout = QGridLayout()
        layout.addWidget(QLabel('Distance restraints'),0,0)
        rif = self.restrain_interfaces = QCheckBox('Inter-chain restraints')
        rif.setChecked(True)
        layout.addWidget(rif, 1,0)
        daw = self.adjust_distances_for_confidence = QCheckBox('Adjust by PAE')
        layout.addWidget(daw, 2,0)
        layout.addWidget(QLabel('Overall strength'),3,0)
        dsl = self.distance_strength_slider = DistanceRestraintWeightSlider(Qt.Orientation.Horizontal)
        dsl.setValue(74) # kappa ~ 10
        layout.addWidget(dsl, 4, 0)

        layout.addWidget(QLabel('Torsion restraints'), 0, 1)
        tsc = self.restrain_sidechains = QCheckBox('Restrain sidechains')
        tsc.setChecked(True)
        layout.addWidget(tsc, 1,1)
        taw = self.adjust_torsions_for_confidence = QCheckBox('Adjust by pLDDT')
        taw.setChecked(True)
        layout.addWidget(taw, 2,1)
        layout.addWidget(QLabel('Overall strength'), 3, 1)
        tsl = self.torsion_strength_slider = TorsionRestraintWeightSlider(Qt.Orientation.Horizontal)
        tsl.setValue(82) # spring constant = 250
        layout.addWidget(tsl, 4, 1)



        self.setContentLayout(layout)


class DistanceRestraintWeightSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(196)
    
    def weight(self):
        '''
        Minimum 2, maximum 100
        '''
        from math import exp, floor
        val = 2*exp(self.value()/50)
        if val < 20:
            return round(val, 1)
        return floor(val)
        
class TorsionRestraintWeightSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(200)
    
    def weight(self):
        '''
        Minimum 50, maximum 100
        '''
        from math import exp, floor
        val = 50*exp(self.value()/51)
        if val < 100:
            return floor(val)
        elif val < 200:
            return 5*round(val/5)
        elif val < 500:
            return round(val, -1)
        elif val < 1000:
            return 25*round(val/25)
        else:
            return 50*round(val/50)



class ExclusiveOrNoneQButtonGroup(QButtonGroup):
    '''
    When in exclusive mode, it is impossible to uncheck a button once it's checked, other than by 
    checking another button in the group. This class acts like an exclusive QButtonGroup, but also 
    allows unchecking of the currently-selected button.
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setExclusive(False)
        self.buttonToggled.connect(self._button_toggled_cb)
    
    def _button_toggled_cb(self, button, checked):
        from ..util import slot_disconnected
        with slot_disconnected(self.buttonToggled, self._button_toggled_cb):
            if checked:
                for b in self.buttons():
                    if b != button:
                        b.setChecked(False)

def _ca_rmsd(model_residues, ref_residues):
    from chimerax.atomic import Residues
    model_residues = Residues(model_residues)
    ref_residues = Residues(ref_residues)
    mcas = model_residues.atoms[model_residues.atoms.names=='CA']
    rcas = ref_residues.atoms[ref_residues.atoms.names=='CA']
    import numpy
    return numpy.sqrt(numpy.mean(numpy.sum((mcas.coords-rcas.scene_coords)**2, axis=1)))
