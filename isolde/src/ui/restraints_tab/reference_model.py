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
        opt = self.options = ReferenceModelOptions(self.session, gui, options_frame)
        opt.expanded.connect(lambda: self.gui.restraints_tab.scroll_area.ensureWidgetVisible(opt))
        ol.addWidget(opt)
        options_frame.setLayout(ol)
        ml.addWidget(options_frame)

        bl = DefaultHLayout()
        rab = self.restrain_all_button = QPushButton('Apply')
        rab.clicked.connect(self.apply_restraints)
        bl.addWidget(rab)
        rsb = self.restrain_selected_button = QPushButton('Apply (selected only)')
        rsb.clicked.connect(lambda: self.apply_restraints(selected_only=True))
        bl.addWidget(rsb)
        bl.addStretch()
        ml.addLayout(bl)

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
        tree.setMinimumHeight(150)

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
            if sm != m:
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

    @property
    def pae_assigned(self):
        cr = self._current_reference
        if cr is None or cr.was_deleted:
            return False
        if hasattr(cr, 'alphafold_pae'):
            return True
        return False

    def _fetch_pae_from_database_if_possible(self, model):
        from chimerax.alphafold.database import uniprot_id_from_filename
        structure_path = getattr(model, 'filename', None)
        if structure_path is not None:
            import os
            structure_path = os.path.basename(structure_path)
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
                # ref chain, alignment score, fractional identity, coverage, model seq, ref seq
                if identity > self.IDENTITY_CUTOFF:
                    res1 = []
                    res2 = []
                    for r1, r2 in zip(s1.residues, s2.residues):
                        if r1 is not None and r2 is not None:
                            res1.append(r1)
                            res2.append(r2)
                    rmsd = _ca_rmsd(res1, res2)
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
                if sm != cm:
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

    def apply_restraints(self, selected_only = False):
        sm = self.isolde.selected_model
        cm = self._current_reference
        distance_pairs = {}
        torsion_pairs = {}
        d_groups = self._restrain_distance_groups
        for group in d_groups:
            mc = group._model_chain
            checked = group.checkedButton()
            if checked is not None:
                distance_pairs[mc] = checked._ref_chain
        t_groups = self._restrain_torsions_groups
        for group in t_groups:
            mc = group._model_chain
            checked = group.checkedButton()
            if checked is not None:
                torsion_pairs[mc] = checked._ref_chain
        from chimerax.core.commands import run
        if len(distance_pairs):
            distance_options = self.options.distance_options
            if not selected_only:
                model_sigs = ','.join([f'"#{sm.id_string}/{mc.chain_id}"' for mc in distance_pairs.keys()])
            else:
                model_sigs = ','.join([f'"#{sm.id_string}/{mc.chain_id}&sel"' for mc in distance_pairs.keys()])
            ref_sigs = ','.join([f'"#{cm.id_string}/{rc.chain_id}"' for rc in distance_pairs.values()])

            cmd = (f'isolde restrain distances {model_sigs} template {ref_sigs} '
                f'perchain {not distance_options["restrain interfaces"]} '
                f'adjustForConfidence {distance_options["adjust for confidence"]} '
                f'useCoordinateAlignment {distance_options["use coordinate alignment"]} '
                f'kappa {distance_options["strength"]} '
                f'fallOff {distance_options["falloff"]}'
            )
            run(self.session, cmd)
        
        torsion_options = self.options.torsion_options
        plddt_reasons = {
            'out of range': 'Some b-factors are over 100, but pLDDT values should be provided on a 0-100 or 0-1 scale.',
            'inverted': 'High pLDDT values correspond to high confidence, but the median value for unstructured residues is higher than the median value for structured ones.',
            'different within residues': 'For models with pLDDT scores in the B-factor column, all atoms in a residue typically have the same score assigned. That is not true for this model.'
        }
        if len(torsion_pairs):
            adjust_for_confidence = torsion_options['adjust for confidence']
            if adjust_for_confidence:
                from chimerax.isolde.reference_model.alphafold import b_factors_look_like_plddt
                looks_like, reason = b_factors_look_like_plddt(mc)
                if not looks_like:
                    from chimerax.isolde.dialog import choice_warning
                    msg = ('You have opted to adjust torsion restraints according to pLDDT, but the values in the B-factor column '
                        f'of the reference model do not look like pLDDT values. {plddt_reasons[reason]}\n'
                        'Would you like to turn off this option before continuing?'
                    )
                    disable = choice_warning(msg)
                    if disable:
                        self.options.adjust_torsions_for_confidence.setChecked(False)
                        torsion_options = self.options.torsion_options

        for mc, rc in torsion_pairs.items():
            if selected_only:
                sel_text = "&sel"
            else:
                sel_text = ""

            cmd = (f'isolde restrain torsions #{sm.id_string}/{mc.chain_id}{sel_text} '
                f'template #{cm.id_string}/{rc.chain_id} '
                f'adjustForConfidence {torsion_options["adjust for confidence"]} '
                f'sidechains {torsion_options["restrain sidechains"]} '
                f'springConstant {torsion_options["spring constant"]} '
                f'alpha {torsion_options["alpha"]}'
            )
            run (self.session, cmd)






class ReferenceModelOptions(CollapsibleArea):
    def __init__(self, session, gui, parent, **kwargs):
        self.session = session
        super().__init__(gui, parent, title="Options", **kwargs)

        layout = QGridLayout()
        layout.addWidget(QLabel('Distance restraints'),0,0)
        rif = self.restrain_interfaces = QCheckBox('Inter-chain restraints')
        rif.setChecked(True)
        layout.addWidget(rif, 1,0)
        daw = self.adjust_distances_for_confidence = QCheckBox('Adjust by PAE')
        daw.setChecked(True)
        layout.addWidget(daw, 2,0)
        drb = self.use_rigid_alignment = QCheckBox('Use rigid-body matching')
        drb.setChecked(False)
        drb.setToolTip('<span>Progressively decompose model and template into matching rigid body regions, '
            'and apply restraints within each rigid body. Recommended if not adjusting for PAE '
            'scores, otherwise <strong>not</strong> recommended.</span>')
        layout.addWidget(drb, 3, 0)
        layout.addWidget(QLabel('Overall strength'),4,0)
        dsl = self.distance_strength_slider = DistanceRestraintWeightSlider(Qt.Orientation.Horizontal)
        dsl.setValue(74) # kappa ~ 10
        layout.addWidget(dsl, 5, 0)
        layout.addWidget(QLabel('Fuzziness'), 6,0)
        dfl = self.distance_fuzziness_slider = DistanceRestraintAlphaSlider(Qt.Orientation.Horizontal)
        dfl.setValue(25)
        layout.addWidget(dfl, 7, 0)
        drfi = self.distance_fuzziness_indicator = DistanceRestraintFuzzinessIndicator(session, dsl, dfl)
        layout.addWidget(drfi, 8,0)
        self.expanded.connect(drfi._update_plot)

        layout.addWidget(QLabel('Torsion restraints'), 0, 1)
        tsc = self.restrain_sidechains = QCheckBox('Restrain sidechains')
        tsc.setChecked(True)
        layout.addWidget(tsc, 1,1)
        taw = self.adjust_torsions_for_confidence = QCheckBox('Adjust by pLDDT')
        taw.setChecked(True)
        layout.addWidget(taw, 2,1)
        layout.addWidget(QLabel('Overall strength'), 4, 1)
        tsl = self.torsion_strength_slider = TorsionRestraintWeightSlider(Qt.Orientation.Horizontal)
        tsl.setValue(82) # spring constant = 250
        layout.addWidget(tsl, 5, 1)
        layout.addWidget(QLabel('Fuzziness'), 6, 1)
        tfl = self.torsion_fuzziness_slider = TorsionRestraintAlphaSlider(Qt.Orientation.Horizontal)
        tfl.setValue(80)
        layout.addWidget(tfl, 7, 1)
        trfi = self.torsion_fuzziness_indicator = DihedralRestraintFuzzinessIndicator(session, tsl, tfl)
        layout.addWidget(trfi, 8,1)
        self.expanded.connect(trfi._update_plot)
        self.setContentLayout(layout)
    
    @property
    def distance_options(self):
        return {
            'restrain interfaces': self.restrain_interfaces.isChecked(),
            'adjust for confidence': self.adjust_distances_for_confidence.isChecked(),
            'use coordinate alignment': self.use_rigid_alignment.isChecked(),
            'strength': self.distance_strength_slider.weight,
            'falloff': self.distance_fuzziness_slider.alpha
        }
    
    @property
    def torsion_options(self):
        return {
            'restrain sidechains': self.restrain_sidechains.isChecked(),
            'adjust for confidence': self.adjust_torsions_for_confidence.isChecked(),
            'spring constant': self.torsion_strength_slider.weight,
            'alpha': self.torsion_fuzziness_slider.alpha
        }



class DistanceRestraintWeightSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(196)
    
    @property
    def weight(self):
        '''
        Minimum 2, maximum 100
        '''
        from math import exp, floor
        val = 2*exp(self.value()/50)
        if val < 20:
            return round(val, 1)
        return floor(val)

class DistanceRestraintAlphaSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(50)
    
    @property
    def alpha(self):
        return self.value()/10
        
class TorsionRestraintWeightSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(200)
    
    @property
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

class TorsionRestraintAlphaSlider(QSlider):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(0)
        self.setMaximum(100)
    
    @property
    def alpha(self):
        return 1-self.value()/100

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

class FuzzinessIndicatorBase(QWidget):
    def __init__(self, *args, tick_step=5, **kwargs):
        super().__init__(*args, **kwargs)
        layout = DefaultHLayout()
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib import ticker

        self.setMaximumHeight(120)
        fig = self.figure = Figure()
        fig.set_tight_layout({'pad':0.25})
        canvas = self.canvas = FigureCanvas(fig)
        loc = self._tick_locator = ticker.MultipleLocator(base=tick_step)
        ax = self.axes = fig.add_subplot(111)
        ax.yaxis.set_major_locator(loc)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
        layout.addWidget(canvas)
        self.setLayout(layout)


class DistanceRestraintFuzzinessIndicator(FuzzinessIndicatorBase):
    def __init__(self, session, weight_slider, alpha_slider, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.weight_slider = weight_slider
        self.falloff_slider = alpha_slider
        self.session = session
        self.well_half_width = 0.1
        self.tolerance = 0.025
        self.target_distance = 3
        self._create_dummy_model()

        weight_slider.valueChanged.connect(self._update_plot)
        alpha_slider.valueChanged.connect(self._update_plot)
        axes = self.axes
        axes.set_ylim([0,10])
        axes.set_xlim([0,2])
        axes.set_xlabel('Deviation (Angstroms)')
        axes.set_ylabel('Force (kJ/mol/A)')

        self.canvas.setToolTip('<span>Example applied force for a target distance of 3 Angstroms and a PAE of 100. '
            'Lower-confidence and/or longer-range restraints will be both fuzzier and weaker than this.</span>')


    
    def _update_plot(self, *_):
        td = self.target_distance
        weight = self.weight_slider.weight
        falloff = self.falloff_slider.alpha
        from math import log, sqrt
        alpha = -falloff*log(self.target_distance)

        import numpy
        offsets = numpy.arange(0,2,0.025)
        forces = numpy.empty(offsets.shape)
        dr = self._dummy_restraint
        dr.tolerance = self.tolerance * td
        dr.c = max(sqrt(td)*self.well_half_width, 0.1)
        dr.alpha = alpha
        dr.target = self.target_distance
        dr.kappa = weight

        da = self._moving_dummy_atom
        for i, os in enumerate(offsets):
            da.coord = (0,0,td+os)
            forces[i] = dr.applied_force
        # Switch to Angstrom units
        forces = forces/10
        if forces.max() > 10:
            self.axes.set_ylim([0,25])
        else:self.axes.set_ylim([0,10])
        
        plot = getattr(self, '_plot', None)
        if plot is None:
            plot = self.axes.plot(offsets, forces, color='blue')[0]
            self._plot = plot
        else:
            #self.axes.clear()
            plot.set_ydata(forces)
            self.axes.draw_artist(plot)
            self.canvas.draw()




    def _create_dummy_model(self):
        from chimerax.atomic import AtomicStructure
        dm = self._dummy_model = AtomicStructure(self.session)
        r = dm.new_residue('DUM', 'A', 0)
        a1 = dm.new_atom('C1', 'C')
        a2 = self._moving_dummy_atom = dm.new_atom('C2', 'C')
        r.add_atom(a1)
        r.add_atom(a2)
        a1.coord = (0,0,0)
        a2.coord = (0,0,4)
        from chimerax.isolde import session_extensions as sx
        drm = sx.get_adaptive_distance_restraint_mgr(dm)
        dr = self._dummy_restraint = drm.add_restraint(a1, a2)
        dr.enabled = True

class DihedralRestraintFuzzinessIndicator(FuzzinessIndicatorBase):
    def __init__(self, session, weight_slider, fuzziness_slider, *args, tick_step=250, **kwargs):
        super().__init__(*args, tick_step=tick_step, **kwargs)
        self.weight_slider = weight_slider
        self.fuzziness_slider = fuzziness_slider
        self.session = session
        self._create_dummy_model()
        self.kappa = 5
        self._angle_step = 5

        weight_slider.valueChanged.connect(self._update_plot)
        fuzziness_slider.valueChanged.connect(self._update_plot)
        axes = self.axes
        axes.set_ylim([0,500])
        axes.set_xlim([0,180])
        axes.set_xlabel('Deviation (degrees)')
        axes.set_ylabel('Torque (kJ/mol)')

        self.canvas.setToolTip('<span>Example applied torque around a dihedral for a residue with a pLDDT of 100. '
            'Lower-confidence restraints will be both fuzzier and weaker than this.</span>')
  
    def _update_plot(self, *_):
        self._dummy_atoms.coords = self._starting_coords
        weight = self.weight_slider.weight
        alpha = self.fuzziness_slider.alpha

        import numpy
        offsets = numpy.arange(0, 180, self._angle_step)
        torques = numpy.empty(offsets.shape)
        dr = self._dummy_restraint
        dr.kappa = self.kappa
        dr.alpha = alpha
        dr.spring_constant = weight

        ma = self._moving_atoms
        axis = self._axis
        from chimerax.geometry import rotation
        r = rotation(axis, self._angle_step, center=self._center)
        torques[0] = dr.applied_moment
        for i in range(len(offsets)-1):
            ma.transform(r)
            torques[i+1]=dr.applied_moment
        if torques.max() > 500:
            self.axes.set_ylim([0,2500])
            self._tick_locator.set_params(500)
        else:
            self.axes.set_ylim([0,500])
            self._tick_locator.set_params(250)
        
        plot = getattr(self, '_plot', None)
        if plot is None:
            plot = self.axes.plot(offsets, torques, color='blue')[0]
            self._plot = plot
        else:
            plot.set_ydata(torques)
            self.axes.draw_artist(plot)
            self.canvas.draw()

    def _create_dummy_model(self):
        from chimerax.atomic import AtomicStructure
        dm = self._dummy_model = AtomicStructure(self.session)
        from chimerax.isolde.atomic.building.build_utils import add_amino_acid_residue
        r = add_amino_acid_residue(dm, 'SER', center=(0,0,0), chain_id='A', number=1)
        self._dummy_atoms = r.atoms
        self._starting_coords = r.atoms.coords
        from chimerax.isolde import session_extensions as sx
        adrm = sx.get_adaptive_dihedral_restraint_mgr(dm)
        dr = self._dummy_restraint = adrm.add_restraint_by_residue_and_name(r, 'chi1')
        dr.target = dr.dihedral.angle
        dr.enabled = True
        ab = dr.dihedral.axial_bond
        self._moving_atoms = ab.side_atoms(ab.smaller_side)
        self._axis = ab.atoms[0].coord - ab.atoms[1].coord
        self._center = ab.atoms[1].coord

