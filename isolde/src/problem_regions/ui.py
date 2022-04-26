
from ..ui.ui_base import UI_Panel_Base

class ProblemAggregatorGUI(UI_Panel_Base):

    def __init__(self, session, isolde, gui, main_frame, category_grid, region_table, bottom_layout, update_button):
        super().__init__(session, isolde, gui, main_frame, sim_sensitive=True)
        from Qt.QtCore import Qt
        self._data_role = Qt.ItemDataRole.UserRole
        cg = self.category_grid = category_grid
        rt = self.region_table = region_table
        rt.itemClicked.connect(self.item_clicked_cb)
        self.update_button = update_button
        update_button.clicked.connect(self.update)

        from .problems import ProblemAggregator
        pa = self.problem_aggregator = ProblemAggregator(session)
        from Qt.QtWidgets import QCheckBox, QLabel, QSpinBox, QDoubleSpinBox, QLabel
        cg.addWidget(QLabel("Unsatisfied restraints"),0,0)
        cg.addWidget(QLabel("Validation issues"),0,1)
        ocb = self.outliers_only_checkbox = QCheckBox("Outliers only")
        cg.addWidget(ocb,0,2)

        self.restraint_checkboxes = []
        for i, name in enumerate(pa.registered_restraint_problem_types):
            cb = QCheckBox(name)
            cb.setChecked(True)
            self.restraint_checkboxes.append(cb)
            cg.addWidget(cb, i+1, 0)
        
        self.validation_checkboxes = []
        for i, name in enumerate(pa.registered_validation_problem_types):
            cb = QCheckBox(name)
            cb.setChecked(True)
            self.validation_checkboxes.append(cb)
            cg.addWidget(cb, i+1, 1)

        csb = self.cutoff_spinbox = QDoubleSpinBox(main_frame)
        csbl = QLabel("Dist cutoff")
        csb.setRange(1.0,10.0)
        csb.setSingleStep(1.0)
        csb.setValue(4.0)
        bottom_layout.insertWidget(0, csb)
        bottom_layout.insertWidget(1, csbl)

        clsb = self.cluster_spinbox = QSpinBox(main_frame)
        clsbl = QLabel("Min cluster size")
        clsb.setRange(2, 20)
        clsb.setValue(5)
        bottom_layout.insertWidget(2, clsb)
        bottom_layout.insertWidget(3, clsbl)


    def sim_end_cb(self, trigger_name, data):
        if self.main_frame.isVisible():
            self.update()
    
    def item_clicked_cb(self, item):
        row = item.row()
        atoms = self.region_table.item(row,0).data(self._data_role)
        if atoms is not None:
            self.session.selection.clear()
            atoms.selected=True
            atoms.intra_bonds.selected=True
            from ..view import focus_on_selection
            focus_on_selection(self.session, atoms)
    
    def update(self, *_):
        t = self.region_table
        t.setRowCount(0)
        m = self.isolde.selected_model
        if m is None:
            return
        from Qt.QtWidgets import QTableWidgetItem
        pa = self.problem_aggregator
        outliers_only = self.outliers_only_checkbox.isChecked()

        restraint_types = [cb.text() for cb in self.restraint_checkboxes if cb.isChecked()]
        validation_types = [cb.text() for cb in self.validation_checkboxes if cb.isChecked()]

        clusters, noise = pa.problem_zones(m, restraint_types=restraint_types, 
            validation_types=validation_types, validation_outliers_only=outliers_only,
            cutoff = self.cutoff_spinbox.value(), min_points=self.cluster_spinbox.value())

        t.clear()
        t.setColumnCount(1+len(restraint_types)+len(validation_types))
        t.setRowCount(len(clusters))
        labels = [t.replace(' ','\n') for t in ['Total']+restraint_types+validation_types]
        t.setHorizontalHeaderLabels(labels)
        for i, cluster in enumerate(clusters):
            atoms = pa.cluster_atoms(cluster)
            issue_count = QTableWidgetItem(str(len(cluster)))
            issue_count.setData(self._data_role,atoms)
            t.setItem(i,0, issue_count)
            for j,it in enumerate(restraint_types+validation_types):
                vtype = pa.registered_type(it)
                vcount = len([v for v in cluster if isinstance(v, vtype)])
                t.setItem(i,j+1,QTableWidgetItem(str(vcount)))

            








