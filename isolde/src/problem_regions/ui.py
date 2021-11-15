
from src.ui.ui_base import UI_Panel_Base
from ..ui.ui_base import UI_Panel_Base

class ProblemAggregatorGUI(UI_Panel_Base):

    def __init__(self, session, isolde, main_frame, category_grid, region_table, update_button):
        super().__init__(session, isolde, main_frame, sim_sensitive=True)
        cg = self.category_grid = category_grid
        self.region_table = region_table
        self.update_button = update_button
        update_button.clicked.connect(self.update)

        from .problems import ProblemAggregator
        pa = self.problem_aggregator = ProblemAggregator(session)
        from Qt.QtWidgets import QCheckBox, QLabel
        cg.addWidget(QLabel("Unsatisfied restraints"),0,0)
        cg.addWidget(QLabel("Validation issues"),0,1)
        ocb = self.outliers_only_checkbox = QCheckBox("Outliers only")
        cg.addWidget(ocb,0,2)

        self.restraint_checkboxes = []
        for i, name in enumerate(pa.registered_restraint_problem_types):
            cb = QCheckBox(name)
            self.restraint_checkboxes.append(cb)
            cg.addWidget(cb, i+1, 0)
        
        self.validation_checkboxes = []
        for i, name in enumerate(pa.registered_validation_problem_types):
            cb = QCheckBox(name)
            self.validation_checkboxes.append(cb)
            cg.addWidget(cb, i+1, 1)


    def sim_end_cb(self, trigger_name, data):
        if self.main_frame.isVisible():
            self.update()
    
    def item_clicked_cb(self, item):
        atoms = item.getData()
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
            validation_types=validation_types, validation_outliers_only=outliers_only)

        t.setRowCount(len(clusters))
        for i, cluster in clusters:
            atoms = pa.cluster_atoms(cluster)
            rank = QTableWidgetItem(str(i+1))
            rank.setData(atoms)
            issue_count = QTableWidgetItem(str(len(cluster)))
            issue_count.setData(atoms)
            t.setItem(i,0, rank)
            t.setItem(i,1, issue_count)

            








