from ..collapse_button import CollapsibleArea

from ..ui_base import (
    UI_Panel_Base, 
    DefaultHLayout, DefaultVLayout, DefaultSpacerItem,
    ExpertModeSelector,
    QDoubleSpinBox
)

from Qt.QtWidgets import (
    QLabel, QSpinBox
)

class SimSelectionPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Mobilisation behaviour', expert_level=ExpertModeSelector.ADVANCED, **kwargs)
        self.setToolTip('Control how the mobile portion of a simulation is defined based on your selection')
        ssd = self.content = SimSelectionDialog(session, isolde, gui, self)
        self.setContentLayout(ssd.main_layout)

class SimSelectionDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=False)
        self.container = collapse_area
        params = isolde.params
        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()
        scbl = DefaultHLayout()
        scbl.addWidget(QLabel('Mobilise selected residue(s)', parent=mf))
        ml.addLayout(scbl)

        extl = DefaultHLayout()
        extl.addWidget(QLabel('plus up to ', parent=mf))
        con_sb = self.connected_residues_spin_box = QSpinBox(mf)
        con_sb.setMinimum(0)
        con_sb.setMaximum(10)
        con_sb.setValue(params.num_selection_padding_residues)
        con_sb.valueChanged.connect(self._connected_residues_sb_cb)
        extl.addWidget(con_sb)
        extl.addWidget(QLabel(' connected residues.', parent=mf))
        extl.addItem(DefaultSpacerItem())
        ml.addLayout(extl)

        expl = DefaultHLayout()
        expl.addWidget(QLabel('plus residues coming within ', parent=mf))
        buf_sb = self.buffer_residues_spin_box = QDoubleSpinBox(mf)
        buf_sb.setMinimum(2)
        buf_sb.setMaximum(10)
        buf_sb.setDecimals(0)
        buf_sb.setValue(params.soft_shell_cutoff_distance)
        buf_sb.valueChanged.connect(self._buffer_residues_sb_cb)
        expl.addWidget(buf_sb)
        expl.addWidget(QLabel(' Angstroms', parent=mf))
        expl.addItem(DefaultSpacerItem())
        ml.addLayout(expl)
    
    def _connected_residues_sb_cb(self, value):
        self.isolde.params.num_selection_padding_residues = value
    
    def _buffer_residues_sb_cb(self, value):
        self.isolde.params.soft_shell_cutoff_distance = value

        

        
