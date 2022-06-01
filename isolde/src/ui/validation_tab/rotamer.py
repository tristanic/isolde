from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QTableWidget, QTableWidgetItem, QPushButton
from Qt.QtGui import QColor, QBrush
from Qt.QtCore import Qt

class RotamerPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Rotamer Validation", **kwargs)
        rd = self.content = RotamerDialog(session, isolde, gui, self)
        self.setContentLayout(rd.main_layout)

class RotamerDialog(UI_Panel_Base):
    COLUMN_LABELS = ('Chain', 'Residue', 'Resname', 'Prior probability (%)')
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
        table.setColumnCount(len(self.COLUMN_LABELS))
        table.setHorizontalHeaderLabels(self.COLUMN_LABELS)
        table.itemClicked.connect(self._item_clicked_cb)

        self.container.expanded.connect(self._populate_table)
        

    def _populate_table(self, *_):
        import numpy
        table = self.table
        table.setRowCount(0)
        if not table.isVisible():
            return
        m = self.isolde.selected_model
        if m is None or m.was_deleted:
            return

        from chimerax.isolde.session_extensions import get_rotamer_mgr
        rota_m = get_rotamer_mgr(self.session)
        if self.isolde.simulation_running:
            residues = self.isolde.sim_manager.sim_construct.mobile_residues
        else:
            residues = m.residues
        rotas = rota_m.get_rotamers(residues)
        iffy, scores = rota_m.non_favored_rotamers(rotas)
        order = numpy.argsort(scores)
        outlier_cutoff = rota_m.cutoffs[1]
        from Qt.QtGui import QColor, QBrush
        from Qt.QtCore import Qt
        badColor = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        table.setRowCount(len(iffy))
        from Qt.QtWidgets import QTableWidgetItem
        for i, index in enumerate(order):
            r = iffy[index]
            score = scores[index]
            res = r.residue
            data = (
                res.chain_id,
                str(res.number),
                res.name,
                '{:.4f}'.format(score*100)
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(Qt.ItemDataRole.UserRole, res)
                if score < outlier_cutoff:
                    item.setBackground(badColor)
                table.setItem(i, j, item)
        if self._first_rebuild and len(order):
            table.resizeColumnsToContents()
            self._first_rebuild = False
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

    

