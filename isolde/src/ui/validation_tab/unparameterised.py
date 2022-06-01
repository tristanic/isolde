from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import (
    QLabel, 
    QTableWidget, QTableWidgetItem, 
    QTreeWidget, QTreeWidgetItem, 
    QPushButton, QCheckBox
)
from Qt.QtCore import Qt, QTimer
USER_ROLE = Qt.ItemDataRole.UserRole


class UnparameterisedResiduesPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Unparameterised Residues", **kwargs)
        urd = self.content = UnparameterisedResiduesDialog(session, isolde, gui, self)
        self.setContentLayout(urd.main_layout)

class UnparameterisedResiduesDialog(UI_Panel_Base):
    H_TO_HEAVY_ATOM_THRESHOLD_RATIO = 0.5
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        rt = self.residue_table = QTableWidget(mf)
        rt.setColumnCount(3)
        rt.setHorizontalHeaderLabels(['Chain', 'Residue', ''])
        rt.itemClicked.connect(self._table_item_clicked_cb)
        rt.setSizeAdjustPolicy(rt.SizeAdjustPolicy.AdjustToContents)
        rt.setMinimumHeight(100)
        rt.setMaximumHeight(200)
        ml.addWidget(rt)
        ml.addWidget(QLabel('Possible templates'))
        tt = self.template_tree = QTreeWidget(mf)
        tt.setHeaderLabels(['MD Template', 'Match score', 'CCD Template', 'Description'])
        tt.setSizeAdjustPolicy(tt.SizeAdjustPolicy.AdjustToContents)
        tt.setMinimumHeight(150)
        tt.setMaximumHeight(200)
        tt.itemClicked.connect(self._tree_item_clicked_cb)
        ml.addWidget(tt)



        bl = DefaultHLayout()
        fb = self.fix_button = QPushButton('Rebuild residue to template', mf)
        fb.setEnabled(False)
        fb.clicked.connect(self._fix_button_clicked_cb)
        bl.addWidget(fb)
        acb = self.fix_all_check_box = QCheckBox('For all residues of this name', mf)
        bl.addWidget(acb)
        bl.addStretch()
        ub = self.update_button = QPushButton('Update')
        ub.clicked.connect(self._populate_unparameterised_residue_table)
        bl.addWidget(ub)
        
        ml.addLayout(bl)
        self._isolde_trigger_handlers.append(isolde.triggers.add_handler(isolde.UNPARAMETERISED_RESIDUE, self._unparam_res_cb))
        self.container.expanded.connect(self._populate_unparameterised_residue_table)


    def _populate_unparameterised_residue_table(self, *_, ff=None, ambiguous=None, unmatched=None, residues=None):
        if self.container.is_collapsed:
            return
        table = self.residue_table
        ttree = self.template_tree
        table.setRowCount(0)
        ttree.clear()
        self.fix_button.setEnabled(False)
        if ambiguous is None and unmatched is None:
            m = self.isolde.selected_model
            if m is None:
                return
            residues = m.residues
            self._ask_to_add_hydrogens_if_necessary(residues)
            from chimerax.atomic import Residues
            residues = Residues(sorted(residues, key=lambda r:(r.chain_id, r.number, r.insertion_code)))
            if ff is None:
                ffmgr = self.isolde.forcefield_mgr
                ff = ffmgr[self.isolde.sim_params.forcefield]
                ligand_db = ffmgr.ligand_db(self.isolde.sim_params.forcefield)
            from chimerax.isolde.openmm.openmm_interface import find_residue_templates, create_openmm_topology
            template_dict = find_residue_templates(residues, ff, ligand_db=ligand_db, logger=self.session.logger)
            top, residue_templates=create_openmm_topology(residues.atoms, template_dict)
            _, ambiguous, unmatched = ff.assignTemplates(top,
                ignoreExternalBonds=True, explicit_templates=residue_templates)
        row_count = len(unmatched)+len(ambiguous)
        if row_count == 0:
            table.setRowCount(1)
            table.setItem(0,1, QTableWidgetItem('No issues detected'))
        table.setRowCount(row_count)
        count = 0
        for r in unmatched:
            by_name, by_comp = ff.find_possible_templates(r)
            cx_res = residues[r.index]
            data = (
                cx_res.chain_id,
                cx_res.name + ' ' + str(cx_res.number),
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(USER_ROLE, (cx_res, by_name, by_comp))
                table.setItem(count, j, item)
            count += 1
        for r, template_info in ambiguous.items():
            cx_res = residues[r.index]
            data = (
                cx_res.chain_id,
                cx_res.name + ' ' + str(cx_res.number),
                ', '.join([ti[0].name for ti in template_info])
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(USER_ROLE, (cx_res, [], [[ti[0].name,0] for ti in template_info]))
                table.setItem(count, j, item)
            count += 1
        table.resizeColumnsToContents()

    def _table_item_clicked_cb(self, item):
        residue, by_name, by_comp = item.data(USER_ROLE)
        ttree = self.template_tree
        ttree.clear()
        self.fix_all_check_box.setChecked(False)
        parent = ttree.invisibleRootItem()
        if len(by_name):
            np = QTreeWidgetItem(parent)
            np.setText(0, 'By residue name')
            for (tname, score) in by_name:
                item = QTreeWidgetItem(np)
                ccd_template, description = _get_ccd_template_and_name(self.session, tname)
                item.setText(0, tname)
                item.setText(1, f'{score:.3f}')
                if ccd_template is not None:
                    item.setText(2, ccd_template.name)
                    item.setText(3, description)
                else:
                    item.setText(2, 'Not found')
                item.setData(0,USER_ROLE, (residue, tname, ccd_template))
        if len(by_comp):
            cp = QTreeWidgetItem(parent)
            cp.setText(0, 'By topology similarity')
            for tname, score in by_comp:
                item = QTreeWidgetItem(cp)
                ccd_template, description = _get_ccd_template_and_name(self.session, tname)
                item.setText(0, tname)
                item.setText(1, f'{score:.3f}')
                if ccd_template is not None:
                    item.setText(2, ccd_template.name)
                    item.setText(3, description)
                else:
                    item.setText(2, 'Not found')
                item.setData(0, USER_ROLE, (residue, tname, ccd_template))
        ttree.expandAll()
        for i in range(4):
            ttree.resizeColumnToContents(i)
        from chimerax.isolde.view import focus_on_selection
        # show all atoms in the residue and its bonded neighbors, to help
        # diagnose bonding errors
        residue.atoms.displays = True
        for r in residue.neighbors:
            r.atoms.displays = True
        focus_on_selection(self.session, residue.atoms)
                        

    def _tree_item_clicked_cb(self, item):
        data = item.data(0, USER_ROLE)
        enable = (data is not None)
        self.fix_button.setEnabled(enable)
    
    def _fix_button_clicked_cb(self, *_):
        ttree = self.template_tree
        table = self.residue_table
        item = ttree.currentItem()
        index = table.currentRow()
        data = item.data(0, USER_ROLE)
        if data is None:
            return
        residue, template_name, ccd_template = data
        fix_all = self.fix_all_check_box.isChecked()
        if fix_all:
            indices = [i for i in range(table.rowCount()) if table.item(i, 0).data(USER_ROLE)[0].name == residue.name]
            residues = [table.item(i,0).data(USER_ROLE)[0] for i in indices]
        else:
            residues=[residue]
            indices = [index]
        from chimerax.isolde.atomic.template_utils import fix_residue_to_match_md_template
        template = self.isolde.forcefield_mgr[self.isolde.sim_params.forcefield]._templates[template_name]
        for i,r in zip(reversed(indices), reversed(residues)):
            fix_residue_to_match_md_template(self.session, r, template, cif_template=ccd_template)
            table.removeRow(i)
        ttree.clear()
        self.fix_button.setEnabled(False)
        
    def _unparam_res_cb(self, *_):
        from chimerax.isolde.dialog import generic_warning
        warn_str = ('At least one residue in your model does not match any templates '
            'in the MD forcefield. This may be due to missing or superfluous atoms, '
            'or an unusual residue that has not yet been parameterised. Launching '
            'the unparameterised residue widget to help sort this out.')
        generic_warning(warn_str)
        mtw = self.gui.main_tab_widget
        mtw.setCurrentWidget(self.gui.validate_tab)
        self.container.expand()
        QTimer.singleShot(self.container.animation_duration, lambda: self.gui.validate_tab.scroll_area.ensureWidgetVisible(self.main_frame))


    def selected_model_changed_cb(self, *_):
        self.residue_table.clearContents()
        self.template_tree.clear()
            
    def _ask_to_add_hydrogens_if_necessary(self, residues):
            h = residues.atoms[residues.atoms.element_names=='H']
            addh = False
            from chimerax.isolde.dialog import choice_warning
            if not len(h):
                addh = choice_warning('This model does not appear to have hydrogens. Would you like to add them first?')
            elif self.suspiciously_low_h(residues):
                addh = choice_warning('This model has significantly fewer hydrogens than expected for a natural molecule. Would you like to run AddH first?')
            elif self.waters_without_h(residues):
                addh = choice_warning('Some or all waters are missing hydrogens. Would you like to add them first?')
            if addh:
                from chimerax.core.commands import run
                run(self.session, f'addh #{residues.unique_structures[0].id_string}')
                # Occasionally addh will add only one hydrogen to a water (typically when too close to a metal). Catch 
                # and fix to avoid user confusion.
                waters = residues[residues.names=='HOH']
                bad = [w for w in waters if len(w.atoms) < 3]
                from chimerax.build_structure import modify_atom
                for b in bad:
                    o = b.find_atom('O')
                    if o is None:
                        self.session.logger.warning(f'Water /{b.chain_id}:{b.number} is missing its O atom. Deleting.')
                        b.delete()
                    modify_atom(o, o.element, 2)

            
    def suspiciously_low_h(self, residues):
        hydrogens = residues.atoms[residues.atoms.element_names=='H']
        heavy_atoms = residues.atoms[residues.atoms.element_names!='H']
        if len(hydrogens)/len(heavy_atoms) < self.H_TO_HEAVY_ATOM_THRESHOLD_RATIO:
            return True
    
    def waters_without_h(self, residues):
        waters = residues[residues.names=='HOH']
        for w in waters:
            if len(w.atoms) != 3:
                return True
        


def _get_ccd_template_and_name(session, tname):
    from chimerax.isolde.openmm.amberff.template_utils import template_name_to_ccd_name
    from chimerax import mmcif
    ccd_name, extra_info = template_name_to_ccd_name(tname)
    if ccd_name is None:
        return (None, '')
    try:
        tmpl = mmcif.find_template_residue(session, ccd_name)
    except ValueError:
        return (None, "No CCD template found")
    description = tmpl.description
    if extra_info is not None:
        description += ', ' + extra_info
    return (tmpl, description)
