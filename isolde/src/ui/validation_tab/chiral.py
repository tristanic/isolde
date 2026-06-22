from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QTableWidget, QTableWidgetItem, QPushButton
from Qt.QtGui import QColor, QBrush
from Qt.QtCore import Qt

class ChiralPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Chirality Validation", **kwargs)
        cd = self.content = ChiralDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class ChiralDialog(UI_Panel_Base):
    COLUMN_LABELS = ('Chain', 'Residue', 'Resname', 'Atom', 'State', '')
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
        if m is None or m.deleted:
            return
        if self.isolde.simulation_running:
            atoms = self.isolde.sim_manager.sim_construct.mobile_atoms
        else:
            atoms = m.atoms
        from chimerax.isolde.atomic.chirality import chiral_outliers
        chirals, oriented, severity = chiral_outliers(self.session, atoms)
        order = numpy.argsort(severity)[::-1]  # worst (most inverted) first
        badColor = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        table.setRowCount(len(chirals))
        for row, index in enumerate(order):
            cc = chirals[int(index)]
            o = oriented[index]
            a = cc.chiral_atom
            res = a.residue
            inverted = o < 0
            state = 'Inverted' if inverted else 'Strained'
            data = (res.chain_id, str(res.number), res.name, a.name, state)
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                # Store the chiral *atom* (not the residue) so clicking focuses on
                # the centre itself and the Flip button knows which centre to act on.
                item.setData(Qt.ItemDataRole.UserRole, a)
                if inverted:
                    item.setBackground(badColor)
                table.setItem(row, j, item)
            # Offer a flip only where the handedness is actually wrong: flipping a
            # merely-strained (but correctly-handed) centre would invert it.
            if inverted:
                from chimerax.core.commands import run

                def cb(_, atom=a, mdl=m):
                    r = atom.residue
                    run(
                        self.session, 'isolde chiralflip #{}/{}:{}{}@{}'.format(
                            mdl.id_string, r.chain_id, r.number, r.insertion_code, atom.name
                        )
                    )

                b = QPushButton('Flip')
                b.clicked.connect(cb)
                table.setCellWidget(row, 5, b)
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
        a = item.data(Qt.ItemDataRole.UserRole)
        if a is None or a.deleted:
            return
        from chimerax.atomic import Atoms
        from chimerax.isolde.view import focus_on_selection
        # Frame the chiral centre and its immediate substituents -- not the whole
        # residue, which over-zooms for large ligands -- then highlight the centre.
        disp = Atoms([a] + list(a.neighbors))
        # Force the centre and all four substituents visible. ISOLDE hides
        # non-polar hydrogens by default (display=False), but the pendant H on a
        # chiral centre is exactly what you need to see to interpret a flip, so
        # un-hide it here.
        disp.displays = True
        focus_on_selection(self.session, disp, pad=1.5)
        a.selected = True

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
