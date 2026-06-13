from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QTableWidget, QTableWidgetItem, QPushButton
from Qt.QtGui import QColor, QBrush
from Qt.QtCore import Qt

class PeptideBondPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Peptide Bond Validation", **kwargs)
        pbd = self.content = PeptideBondDialog(session, isolde, gui, self)
        self.setContentLayout(pbd.main_layout)

class PeptideBondDialog(UI_Panel_Base):
    NUM_COLUMNS = 4
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=True):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)

        self.container=collapse_area
        mf = self.main_frame
        ml = self.main_layout=DefaultVLayout()
        self._first_rebuild = True
        self._temporary_isolde_handlers = []
        
        table = self.table = QTableWidget(mf)
        ml.addWidget(table)
        hl = DefaultHLayout()
        hl.addStretch()
        b = self.update_button = QPushButton('Update')
        b.clicked.connect(self._populate_table)
        hl.addWidget(b)
        ml.addLayout(hl)
        table.setColumnCount(self.NUM_COLUMNS)
        table.setHorizontalHeaderLabels(['Chain', 'Residues', 'Conformation', ''])
        table.itemClicked.connect(self._item_clicked_cb)

        self.container.expanded.connect(self._populate_table)
        

    def _populate_table(self, *_):
        table = self.table
        table.setRowCount(0)
        m = self.isolde.selected_model
        if m is None or m.deleted:
            return
        from chimerax.isolde.session_extensions import get_proper_dihedral_mgr
        from chimerax.isolde.validation.cmd import classify_peptide_bonds
        pdm = get_proper_dihedral_mgr(self.session)
        if self.isolde.simulation_running:
            residues = self.isolde.sim_manager.sim_construct.mobile_residues
        else:
            residues = m.residues
        iffy = classify_peptide_bonds(pdm, residues)

        table.setRowCount(len(iffy))

        cis_nonpro_color = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        cis_pro_color = QBrush(QColor(100,255,100), Qt.SolidPattern)
        twisted_color = QBrush(QColor(240, 200, 160), Qt.SolidPattern)
        for i, it in enumerate(iffy):
            res1 = it['res1']
            res2 = it['res2']
            angle = it['omega_deg']
            if it['is_cis']:
                conf_text = 'cis'
                color = cis_pro_color if it['is_proline'] else cis_nonpro_color
            else:
                conf_text = 'twisted'
                color = twisted_color
            data = (
                res1.chain_id,
                f'{res1.name}:{res1.number}-{res2.name}:{res2.number}',
                f'{angle:.0f}° ({conf_text})'
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(Qt.ItemDataRole.UserRole, res2)
                item.setBackground(color)
                table.setItem(i, j, item)
            if it['is_cis']:
                from chimerax.core.commands import run
                def cb(_, r=res2):
                    run(self.session, f'isolde cisflip #{m.id_string}/{r.chain_id}:{r.number}')
                b = QPushButton('Flip')
                b.clicked.connect(cb)
                table.setCellWidget(i, 3, b)
        if self._first_rebuild and len(iffy):
            table.resizeColumnsToContents()
            self._first_rebuild=False
        table.resizeRowsToContents()


                
    def sim_start_cb(self, trigger_name, data):
        if not self.container.is_collapsed:
            self._populate_table()

    def sim_end_cb(self, trigger_name, data):
        if not self.container.is_collapsed:
            self._populate_table()
    
    def selected_model_changed_cb(self, trigger_name, m):
        tih = self._temporary_isolde_handlers
        while len(tih):
            tih.pop().remove()
        if m is not None:
            tih.append(m.triggers.add_handler('changes', self._model_changes_cb))
        if not self.container.is_collapsed:
            self._populate_table()

    def _item_clicked_cb(self, item):
        res = item.data(Qt.ItemDataRole.UserRole)
        m = self.isolde.selected_model
        from chimerax.core.commands import run
        run(self.session, f'isolde step #{m.id_string}/{res.chain_id}:{res.number}', log=False)


    def _model_changes_cb(self, trigger_name, changes):
        if self.container.is_collapsed:
            return
        if changes[1].num_deleted_atoms():
            # Rebuild table to purge any deleted residues
            def changes_done_cb(*_):
                if not self.container.is_collapsed:
                    self._populate_table()
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
            from chimerax.atomic import get_triggers
            get_triggers().add_handler('changes done', changes_done_cb)

    def cleanup(self):
        tih = self._temporary_isolde_handlers
        while len(tih):
            tih.pop().remove()
        super().cleanup()

    

