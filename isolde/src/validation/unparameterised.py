# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 28-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

# Make sure core CIF templates are pre-loaded
from chimerax.isolde import atomic

class Unparameterised_Residues_Mgr:

    def __init__(self, isolde):
        self.isolde = isolde
        self.session = isolde.session
        iw = isolde.iw
        mf = self._main_frame = iw._validate_unparameterised_residues_main_frame
        mf.hide()
        self._stub_frame = iw._validate_unparameterised_residues_stub_frame
        self._show_button = iw._validate_unparameterised_residues_show_button
        self._hide_button = iw._validate_unparameterised_residues_hide_button
        self._update_button = iw._validate_unparameterised_residues_update_button
        self._residue_table = iw._validate_unparameterised_residues_table
        self._template_list = iw._validate_possible_templates_list
        self._fix_button = iw._validate_fix_to_template_button
        self._do_for_all = iw._validate_fix_to_template_all_residues_checkbox

        self._unrecognised_residue_handler = None
        self._sim_start_handler = isolde.triggers.add_handler(
            'simulation started', self._sim_start_cb
        )
        self._sim_end_handler = isolde.triggers.add_handler(
            'simulation terminated', self._sim_end_cb
        )
        self._connect_functions()

    def _connect_functions(self):
        self._show_button.clicked.connect(self._show_main_frame)
        self._hide_button.clicked.connect(self._hide_main_frame)
        self._update_button.clicked.connect(self._update_unparameterised_residues_list)
        self._residue_table.itemClicked.connect(self._show_selected_unparameterised_residue)
        self._fix_button.clicked.connect(self._fix_selected_unparameterised_residue)

    def _show_main_frame(self):
        self._stub_frame.hide()
        self._main_frame.show()
        self._update_unparameterised_residues_list()

    def _hide_main_frame(self):
        self._stub_frame.show()
        self._main_frame.hide()

    def _update_unparameterised_residues_list(self, *_, ff=None, ambiguous=None, unmatched=None, residues=None):
        table = self._residue_table
        tlist = self._template_list
        if not table.isVisible():
            return
        table.setRowCount(0)
        tlist.clear()
        if ambiguous is None and unmatched is None:
            if self.isolde.selected_model is None:
                return
            residues = self.isolde.selected_model.residues
            h = residues.atoms[residues.atoms.element_names =='H']
            if not len(h):
                from ..dialog import choice_warning
                addh = choice_warning('This model does not appear to have hydrogens. Would you like to add them first?')
                if addh:
                    from chimerax.core.commands import run
                    run(self.session, 'addh')
            from chimerax.atomic import Residues
            residues = Residues(sorted(residues, key=lambda r:(r.chain_id, r.number, r.insertion_code)))
            if ff is None:
                ffmgr = self.isolde.forcefield_mgr
                ff = ffmgr[self.isolde.sim_params.forcefield]
                ligand_db = ffmgr.ligand_db(self.isolde.sim_params.forcefield)
            from ..openmm.openmm_interface import find_residue_templates, create_openmm_topology
            template_dict = find_residue_templates(residues, ff, ligand_db=ligand_db, logger=self.session.logger)
            top, residue_templates=create_openmm_topology(residues.atoms, template_dict)
            _, ambiguous, unmatched = ff.assignTemplates(top,
                ignoreExternalBonds=True, explicit_templates=residue_templates)
        from PyQt5.QtWidgets import QTableWidgetItem
        table.setRowCount(len(unmatched)+len(ambiguous))
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
                item.data = (cx_res, by_name, by_comp)
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
                item.data = (cx_res, [], [[ti[0].name,0] for ti in template_info])
                table.setItem(count, j, item)
            count += 1
        table.resizeColumnsToContents()

    def _show_selected_unparameterised_residue(self, item):
        residue, by_name, by_comp = item.data
        tlist = self._template_list
        tlist.clear()
        self._do_for_all.setCheckState(False)
        from PyQt5.QtWidgets import QListWidgetItem
        tlist.addItem(QListWidgetItem("Matches by residue name"))
        def get_ccd_template_and_name(session, tname):
            from ..openmm.amberff.template_utils import template_name_to_ccd_name
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


        for (tname, score) in by_name:
            ccd_template, description = get_ccd_template_and_name(self.session, tname)
            entry_text = tname
            entry_text += " (Score: {:.3f})".format(score)
            if len(description):
                entry_text += " " + description
            entry = QListWidgetItem(entry_text)
            entry.data = (residue, tname, ccd_template)
            tlist.addItem(entry)
        tlist.addItem(QListWidgetItem("Matches by similar topology"))
        for (tname, score) in by_comp:
            ccd_template, description = get_ccd_template_and_name(self.session, tname)
            entry_text = tname
            entry_text += " (Score: {:.3f})".format(score)
            if len(description):
                entry_text += " " + description
            entry = QListWidgetItem(entry_text)
            entry.data = (residue, tname, ccd_template)
            tlist.addItem(entry)

        tlist.repaint()
        from ..view import focus_on_selection
        # show all atoms in the residue and its bonded neighbors, to help
        # diagnose bonding errors
        residue.atoms.displays = True
        for r in residue.neighbors:
            r.atoms.displays = True
        focus_on_selection(self.session, residue.atoms)

    def _fix_selected_unparameterised_residue(self, *_):
        tlist = self._template_list
        table = self._residue_table
        item = tlist.currentItem()
        index = table.currentRow()
        if item is None:
            return
        if not hasattr(item, 'data') or item.data is None:
            return
        residue, template_name, ccd_template = item.data
        fix_all = self._do_for_all.isChecked()
        if fix_all:
            indices = [i for i in range(table.rowCount()) if table.item(i, 0).data[0].name == residue.name]
            residues = [table.item(i,0).data[0] for i in indices]
        else:
            residues=[residue]
            indices = [index]
        from chimerax.isolde.atomic.template_utils import fix_residue_to_match_md_template
        template = self.isolde.forcefield_mgr[self.isolde.sim_params.forcefield]._templates[template_name]
        for i,r in zip(reversed(indices), reversed(residues)):
            fix_residue_to_match_md_template(self.session, r, template, cif_template=ccd_template)
            table.removeRow(i)
        tlist.clear()

    def _sim_start_cb(self, *_):
        pass
        # sh = self.isolde.sim_handler
        # sh.triggers.add_handler('unparameterised residue', self._sim_unparam_res_cb)

    def _sim_end_cb(self, *_):
        self._hide_main_frame()

    def _sim_unparam_res_cb(self, *_):
        from ..dialog import generic_warning
        warn_str = ('At least one residue in your model does not match any templates '
            'in the MD forcefield. This may be due to missing or superfluous atoms, '
            'or an unusual residue that has not yet been parameterised. Launching '
            'the unparameterised residue widget to help sort this out.')
        generic_warning(warn_str)
        iw = self.isolde.iw
        st = iw._sim_tab_widget
        st.setCurrentIndex(2)
        self._show_main_frame()
