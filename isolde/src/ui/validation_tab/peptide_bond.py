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

        self._isolde_trigger_handlers.append(isolde.triggers.add_handler('selected model changed', self.selected_model_changed_cb))
        self.container.expanded.connect(self._populate_table)
        

    def _populate_table(self, *_):
        import numpy
        table = self.table
        table.setRowCount(0)
        m = self.isolde.selected_model
        if m is None or m.was_deleted:
            return
        from chimerax.isolde.session_extensions import get_proper_dihedral_mgr
        pdm = get_proper_dihedral_mgr(self.session)
        if self.isolde.simulation_running:
            residues = self.isolde.sim_manager.sim_construct.mobile_residues
        else:
            residues = m.residues
        omegas = pdm.get_dihedrals(residues, 'omega')
        abs_angles = numpy.abs(omegas.angles)
        from math import pi
        from chimerax.isolde.constants import defaults
        cc = defaults.CIS_PEPTIDE_BOND_CUTOFF
        tc = defaults.TWISTED_PEPTIDE_BOND_DELTA
        cis_mask = abs_angles < cc
        twisted_mask = numpy.logical_and(abs_angles >= cc, abs_angles < pi-tc)
        iffy_mask = numpy.logical_or(cis_mask, twisted_mask)
        iffy = omegas[iffy_mask]
        angles = numpy.degrees(iffy.angles)
        cis_mask = cis_mask[iffy_mask]

        table.setRowCount(len(iffy))

        cis_nonpro_color = QBrush(QColor(255, 100, 100), Qt.SolidPattern)
        cis_pro_color = QBrush(QColor(100,255,100), Qt.SolidPattern)
        twisted_color = QBrush(QColor(240, 200, 160), Qt.SolidPattern)
        for i, (omega, angle, cis) in enumerate(zip(iffy, angles, cis_mask)):
            res1, res2 = omega.atoms.unique_residues
            if cis:
                conf_text = 'cis'
            else:
                conf_text = 'twisted'
            data = (
                res1.chain_id,
                f'{res1.name}:{res1.number}-{res2.name}:{res2.number}',
                f'{angle:.0f}° ({conf_text})'
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.setData(Qt.ItemDataRole.UserRole, res2)
                if cis:
                    if res2.name == 'PRO':
                        color = cis_pro_color
                    else:
                        color = cis_nonpro_color
                else:
                    color = twisted_color
                item.setBackground(color)
                table.setItem(i, j, item)
            if cis:
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
    
    def selected_model_changed_cb(self, trigger_name, data):
        if not self.container.is_collapsed:
            self._populate_table()

    def _item_clicked_cb(self, item):
        res = item.data(Qt.ItemDataRole.UserRole)
        m = self.isolde.selected_model
        from chimerax.core.commands import run
        run(self.session, f'isolde step #{m.id_string}/{res.chain_id}:{res.number}', log=False)





    

