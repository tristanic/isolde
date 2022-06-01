
from ..ui.ui_base import (
    UI_Panel_Base, IsoldeTab,
    DefaultHLayout, DefaultVLayout,
    QSpinBox, QDoubleSpinBox
)
from Qt.QtWidgets import QWidget, QGridLayout, QTableWidget, QPushButton
from Qt.QtCore import Qt
USER_ROLE = Qt.ItemDataRole.UserRole


class ProblemAggregatorTab(IsoldeTab):
    def __init__(self, session, isolde, gui, tab_widget, tab_name='Problem Zones'):
        super().__init__(session, isolde, gui, tab_widget, tab_name)
        aggw = QWidget()
        agg = self.aggregator = ProblemAggregatorGUI(session, isolde, gui, aggw)
        self.addWidget(aggw)

class ProblemAggregatorGUI(UI_Panel_Base):

    def __init__(self, session, isolde, gui, main_frame):
        super().__init__(session, isolde, gui, main_frame, sim_sensitive=True)
        ml = self.main_layout = DefaultVLayout()
        cg = self.category_grid = QGridLayout()
        ml.addLayout(cg)
        main_frame.setLayout(ml)

        tl = DefaultHLayout()
        ml.addLayout(tl)
        rt = self.region_table = QTableWidget()
        rt.setMinimumHeight(250)
        tl.addWidget(rt)


        rt.itemClicked.connect(self.item_clicked_cb)

        from .problems import ProblemAggregator
        pa = self.problem_aggregator = ProblemAggregator(session)
        from ..ui.ui_base import QSpinBox, QDoubleSpinBox
        from Qt.QtWidgets import QCheckBox, QLabel
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

        bottom_layout = DefaultHLayout()
        ml.addLayout(bottom_layout)

        csb = self.cutoff_spinbox = QDoubleSpinBox(main_frame)
        csbl = QLabel("Dist cutoff")
        csb.setRange(1.0,10.0)
        csb.setSingleStep(1.0)
        csb.setValue(4.0)
        bottom_layout.addWidget(csb)
        bottom_layout.addWidget(csbl)

        clsb = self.cluster_spinbox = QSpinBox(main_frame)
        clsbl = QLabel("Min cluster size")
        clsb.setRange(2, 20)
        clsb.setValue(5)
        bottom_layout.addWidget(clsb)
        bottom_layout.addWidget(clsbl)

        bottom_layout.addStretch()
        update_button = self.update_button = QPushButton('Update')
        update_button.clicked.connect(self.update)
        bottom_layout.addWidget(update_button)

        ith = self._isolde_trigger_handlers

    def sim_start_cb(self, *_):
        self.region_table.clear()

    def sim_end_cb(self, trigger_name, data):
        if self.main_frame.isVisible():
            self.update()
    
    def selected_model_changed_cb(self, *_):
        self.region_table.clearContents()

    def item_clicked_cb(self, item):
        row = item.row()
        atoms = self.region_table.item(row,0).data(USER_ROLE)
        if atoms is not None:
            self.session.selection.clear()
            atoms.selected=True
            atoms.intra_bonds.selected=True
            from chimerax.core.commands import run
            run(self.session, 'view sel', log=False)

    
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

        from collections import Counter
        ccounter = Counter([type(v) for cluster in clusters for v in cluster])

        sorted_counts = list(sorted(ccounter.items(), key=lambda c: c[1], reverse=True))
        names = [pa.registered_name(c[0]) for c in sorted_counts]

        t.clear()
        t.setColumnCount(1+len(names))
        t.setRowCount(len(clusters))
        labels = [t.replace(' ','\n') for t in ['Total']+names]
        t.setHorizontalHeaderLabels(labels)
        for i, cluster in enumerate(clusters):
            atoms = pa.cluster_atoms(cluster)
            issue_count = QTableWidgetItem(str(len(cluster)))
            issue_count.setData(USER_ROLE,atoms)
            t.setItem(i,0, issue_count)
            for j,it in enumerate(names):
                vtype = pa.registered_type(it)
                vcount = len([v for v in cluster if isinstance(v, vtype)])
                t.setItem(i,j+1,QTableWidgetItem(str(vcount)))

            








