from ..collapse_button import CollapsibleArea
from ..ui_base import UI_Panel_Base, DefaultHLayout, DefaultVLayout
from Qt.QtWidgets import QTableWidget, QTableWidgetItem, QPushButton
from Qt.QtGui import QColor, QBrush
from Qt.QtCore import Qt

class ClashesPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title="Clashes", **kwargs)
        cd = self.content = ClashesDialog(session, isolde, gui, self)
        self.setContentLayout(cd.main_layout)

class ClashesDialog(UI_Panel_Base):
    COLUMN_LABELS = ('Atom 1','Atom 2', 'Overlap')
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=True):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        self.container = collapse_area
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        table = self.table = QTableWidget(mf)
        table.setColumnCount(len(self.COLUMN_LABELS))
        table.setHorizontalHeaderLabels(self.COLUMN_LABELS)
        table.itemClicked.connect(self._item_clicked_cb)
        ml.addWidget(table)
        hl = DefaultHLayout()
        hl.addStretch()
        b = self.update_button = QPushButton('Update')
        b.clicked.connect(self._populate_table)
        hl.addWidget(b)
        ml.addLayout(hl)
        self.container.expanded.connect(self._populate_table)

    @property
    def atoms(self):
        sm = self.isolde.selected_model
        if sm is None:
            from chimerax.atomic import Atoms
            return Atoms()
        if self.isolde.simulation_running:
            return self.isolde.sim_manager.sim_construct.mobile_atoms
        return sm.atoms
    
    def sim_start_cb(self, *_):
        sh = self.isolde.sim_handler
        sh.triggers.add_handler('clash detected', self._severe_clash_cb)
    
    def _severe_clash_cb(self, *_):
        self.isolde.sim_manager.pause=True
        from chimerax.isolde.dialog import generic_warning
        msg_string = ('ISOLDE has detected severe clashes in the model that the '
            'minimiser is unable to reconcile on its own. The simulation '
            'is still initialised, but cannot continue until these are '
            'corrected. Clicking OK will open the clash validation panel with '
            'a list of atoms currently experiencing extreme forces. In most '
            'cases these can be dealt with by choosing more appropriate rotamers. '
            'For more extreme issues you may need to stop the simulation and '
            'tell ISOLDE to temporarily ignore some residues using the command '
            '"isolde ignore {selection}". Once you have moved the adjacent '
            'atoms into non-clashing positions, you can stop the simulation and '
            'stop ignoring the residues using "isolde ~ignore". \n'
            'NOTE: the true culprit may not necessarily be at the top of the '
            'list. Look for clashes that have no direct path away from each '
            'other (e.g. bonds threaded through rings). After rearranging atoms '
            'you may check to see if this has solved the problem by pressing '
            'the play button.'
            )
        generic_warning(msg_string)
        self.gui.main_tab_widget.setCurrentWidget(self.gui.validate_tab)
        if self.container.is_collapsed:
            self.container.expand()
        else:
            self._populate_table()

    def _populate_table(self, *_):
        atoms = self.atoms
        t = self.table
        if self.container.is_collapsed:
            return
        t.setRowCount(0)
        if not len(atoms):
            return
        from chimerax.clashes import clashes
        clash_dict = clashes.find_clashes(self.session, atoms, inter_model=False)
        if not len(clash_dict):
            return
        from functools import reduce
        # Clash dict is bi-directional, so number of clashes is half the total count
        clash_count = 1/2 * reduce(lambda x,y: x+y, (len(d.keys()) for d in clash_dict.values()))
        t.setRowCount(clash_count)
        seen = set()
        from chimerax.atomic import Atoms
        # Make list of unique clashing atom pairs and their distances
        clash_list = []
        for a1, clashes in clash_dict.items():
            seen.add(a1)
            for a2, dist in clashes.items():
                if a2 in seen:
                    continue
                clash_list.append((Atoms([a1,a2]), dist))
        # Sort clashes in decreasing order of overlap
        clash_list = sorted(clash_list, key=lambda x: x[1], reverse=True)

        from Qt.QtWidgets import QTableWidgetItem
        for i, (catoms, overlap) in enumerate(clash_list):
            a1, a2 = catoms
            r1, r2 = catoms.residues
            data = (
            "{} {}{}: {}".format(r1.name, r1.chain_id, r1.number, a1.name),
            "{} {}{}: {}".format(r2.name, r2.chain_id, r2.number, a2.name),
            "{:0.2f}".format(overlap)
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(Qt.ItemDataRole.UserRole, catoms)
                t.setItem(i, j, item)
        t.resizeColumnsToContents()


    def _item_clicked_cb(self, item):
        atoms = item.data(Qt.ItemDataRole.UserRole)
        from chimerax.isolde.view import focus_on_selection
        self.session.selection.clear()
        atoms.selecteds=True
        atoms.displays=True
        focus_on_selection(self.session, atoms)

        