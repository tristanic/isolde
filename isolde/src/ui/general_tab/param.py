from inspect import Parameter
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout, DefaultSpacerItem
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import QPushButton, QFileDialog, QCheckBox

class ParameterisePanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, 'Residue Parameterisation', **kwargs)
        prd = self.content = ParameteriseDialog(session, isolde, gui, self)
        self.setContentLayout(prd.main_layout)

class ParameteriseDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area)
        self._dont_warn_again = False
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        ldl = DefaultHLayout()
        ldb = self.load_defs_button = QPushButton('Load residue parameters', mf)
        ldb.setToolTip('<span>Load molecular dynamics parameterisation(s) for one or more'
            'residues. These must be in the ffXML format used by OpenMM - '
            'standard crystallographic .cif restraints are NOT supported.</span>')
        
        ldb.clicked.connect(self._load_params_cb)
        ldl.addWidget(ldb)
        ldl.addItem(DefaultSpacerItem())
        ml.addLayout(ldl)
        prl = DefaultHLayout()
        prb = self.parameterise_residue_button = QPushButton('Parameterise selected residue', mf)
        from chimerax.isolde.openmm.amberff.parameterise import supported_elements
        prb.setToolTip('<span>Parameterise selected residue using ANTECHAMBER. Currently limited to '
            f'druglike molecules composed of the elements {",".join(supported_elements)} '
            'with no covalent links to other residues). This can take a while for large residues.</span>')
        prb.clicked.connect(self._parameterise_sel_cb)
        # Start disabled; only enabled when single qualifying residue selected
        prb.setEnabled(False)
        prl.addWidget(prb)
        pro = self.parameterise_override_checkbox = QCheckBox('Override existing', parent=mf)
        pro.setChecked(False)
        pro.setToolTip('Force generated parameters to override any others with the same residue name')
        prl.addWidget(pro)
        prl.addItem(DefaultSpacerItem())
        ml.addLayout(prl)

        from chimerax.core.selection import SELECTION_CHANGED
        self._chimerax_trigger_handlers.append(self.session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed_cb))
    
    def _load_params_cb(self, *_):
        files = self._choose_ff_files(self.load_defs_button)
        if files is not None and len(files):
            ff = self.isolde.forcefield_mgr[self.isolde.sim_params.forcefield]
            for f in files:
                try:
                    ff.loadFile(f, resname_prefix='USER_')
                except ValueError as e:
                    self.session.logger.warning(f'Failed to add {f}: {str(e)}')
                except:
                    raise

    def _choose_ff_files(self, parent):
        caption = 'Choose one or more ffXML files'
        filetypes = 'ffXML files (*.xml)'
        dlg = QFileDialog(caption=caption, parent=parent)
        dlg.setAcceptMode(QFileDialog.AcceptOpen)
        dlg.setNameFilter(filetypes)
        dlg.setFileMode(QFileDialog.ExistingFiles)
        if dlg.exec():
            return dlg.selectedFiles()
    
    def _parameterise_sel_cb(self, *_):
        from chimerax.isolde.dialog import choice_warning
        go = True
        if not self._dont_warn_again:
            from chimerax.atomic import selected_residues
            r = selected_residues(self.session)[0]
            from chimerax.isolde.dialog import choice_warning
            go, dont_ask_again = choice_warning('This will run ANTECHAMBER to perform a semi-empirical parameterisation on the '
                'selected residue. It is critically important that you make sure that all atoms are present and '
                'hydrogens are correct before going forward. Note that this can take a long time (many minutes) '
                'for very large ligands. If automatic charge assignment fails, try manually running '
                f'the command "isolde param #{r.structure.id_string}{r.atomspec}". \n'
                'Would you like to continue?', allow_dont_ask_again=True
                )
            self._dont_warn_again = dont_ask_again
            
        if go:
            from chimerax.core.commands import run
            run(self.session, f'isolde param sel override {self.parameterise_override_checkbox.isChecked()}')

    def _selection_changed_cb(self, *_):
        from chimerax.atomic import selected_residues
        sel = selected_residues(self.session)
        prb = self.parameterise_residue_button
        enabled = False
        if len(sel) == 1 and len(sel[0].neighbors) == 0:
            from chimerax.isolde.openmm.amberff.parameterise import residue_qualifies
            if residue_qualifies:
                enabled=True
        prb.setEnabled(enabled)

    