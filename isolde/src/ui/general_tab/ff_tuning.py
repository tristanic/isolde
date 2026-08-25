# @Author: Tristan Croll
# @Date:   25-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 25-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Developer-gated per-force-field simulation-parameter tuning dashboard.

Exposes the curated tunable subset (openmm/sim_param_meta.py) as editors, scoped
to the currently selected force field: edits, per-parameter reset, a per-field
factory-reset, and an explicit "Save current settings" (persist across sessions).
Every action is issued as an ``isolde simparam ...`` command (GUI-logs-its-command
convention), so the dashboard is a thin, scriptable front-end over that command.
Nothing here persists unless the user hits Save, so tinkering can't jeopardise a
future session.
'''
from ..ui_base import (
    UI_Panel_Base, DefaultVLayout, DefaultHLayout, HorizontalLine, ExpertModeSelector,
)
from ..collapse_button import CollapsibleArea

from Qt.QtWidgets import (
    QLabel, QComboBox, QCheckBox, QDoubleSpinBox, QPushButton, QGridLayout, QWidget,
)


class ForceFieldTuningPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        from .forcefield import FORCEFIELD_UI_EXPERT_LEVEL
        super().__init__(gui, parent, title='Force Field Tuning (advanced)',
                         expert_level=FORCEFIELD_UI_EXPERT_LEVEL, **kwargs)
        d = self.content = ForceFieldTuningDialog(session, isolde, gui, self)
        self.setContentLayout(d.main_layout)


class ForceFieldTuningDialog(UI_Panel_Base):
    def __init__(self, session, isolde, gui, collapse_area):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=False)
        self.container = collapse_area
        from chimerax.isolde.openmm import sim_param_meta as meta
        self._meta = meta

        mf = self.main_frame
        ml = self.main_layout = DefaultVLayout()

        # --- header: active force field + save / factory-reset ---
        hl = DefaultHLayout()
        hl.addWidget(QLabel('Tuning force field:', mf))
        self._ff_label = QLabel('', mf)
        self._ff_label.setStyleSheet('font-weight: bold;')
        hl.addWidget(self._ff_label)
        hl.addStretch()
        save = QPushButton('Save current settings', mf)
        save.setToolTip('Persist this force field\'s current settings so they are '
                        'restored in future sessions. Nothing is saved until you click this.')
        save.clicked.connect(lambda *_: self._run('isolde simparam save'))
        hl.addWidget(save)
        fr = QPushButton('Factory reset', mf)
        fr.setToolTip('Discard all of this force field\'s tuning and restore its '
                      'built-in defaults (this session; not saved until you Save).')
        fr.clicked.connect(lambda *_: self._run('isolde simparam reset'))
        hl.addWidget(fr)
        ml.addLayout(hl)
        ml.addWidget(HorizontalLine())

        # --- per-parameter editor grid ---
        grid = self._grid = QGridLayout()
        self._rows = {}     # name -> dict(widgets)
        for row, name in enumerate(meta.tunable_params()):
            self._build_row(grid, row, name)
        ml.addLayout(grid)

        sp = isolde.sim_params
        self._handler = sp.triggers.add_handler(sp.PARAMETER_CHANGED, self._param_changed_cb)
        self._refresh_all()

    # ---------------------------------------------------------------- build
    def _build_row(self, grid, row, name):
        mf = self.main_frame
        spec = self._meta.spec_for(name)
        mod = QLabel('', mf)
        mod.setStyleSheet('color: #d98014;')           # orange dot when overridden
        mod.setToolTip('Modified from this force field\'s default')
        mod.setFixedWidth(12)
        label = QLabel(spec.label, mf)
        label.setToolTip(spec.tooltip)
        editor = self._make_editor(name, spec)
        reset = QPushButton('↺', mf)              # ↺
        reset.setToolTip('Reset "{}" to this force field\'s default'.format(spec.label))
        reset.setFixedWidth(30)
        reset.clicked.connect(lambda *_, n=name: self._run('isolde simparam {} default'.format(n)))
        grid.addWidget(mod, row, 0)
        grid.addWidget(label, row, 1)
        grid.addWidget(editor, row, 2)
        grid.addWidget(reset, row, 3)
        self._rows[name] = dict(mod=mod, label=label, editor=editor, reset=reset)

    def _make_editor(self, name, spec):
        mf = self.main_frame
        if spec.kind == 'bool':
            w = QCheckBox(mf)
            w.stateChanged.connect(lambda *_, n=name: self._commit(n))
            return w
        if spec.kind == 'enum':
            w = QComboBox(mf)
            w.addItems([tok for tok, _ in spec.choices])
            w.currentIndexChanged.connect(lambda *_, n=name: self._commit(n))
            return w
        w = QDoubleSpinBox(mf)
        if spec.minimum is not None:
            w.setMinimum(float(spec.minimum))
        if spec.maximum is not None:
            w.setMaximum(float(spec.maximum))
        if spec.step is not None:
            w.setSingleStep(float(spec.step))
        w.setDecimals(2 if (spec.step is not None and spec.step < 1) else 1)
        w.setMaximumWidth(90)
        # commit on Enter / focus-out only (not per spin step) to keep the log clean
        w.editingFinished.connect(lambda n=name: self._commit(n))
        return w

    # --------------------------------------------------------------- helpers
    def _run(self, command):
        from chimerax.core.commands import run
        run(self.session, command)

    def _display_value(self, name, spec):
        '''Current SimParams value as a plain editor value (float / bool / token).'''
        from openmm.unit import Quantity
        sp = self.isolde.sim_params
        v = sp[name]
        if spec.kind == 'enum':
            return spec.token_for(v)
        if isinstance(v, Quantity):
            unit = sp._default_params[name][1]
            return float(v.value_in_unit(unit))
        return v

    def _commit(self, name):
        '''An editor changed: log+run the command, unless the value already matches
        (avoids redundant commands from focus-out / programmatic refresh).'''
        spec = self._meta.spec_for(name)
        w = self._rows[name]['editor']
        if spec.kind == 'bool':
            new = bool(w.isChecked())
            if new == bool(self._display_value(name, spec)):
                return
            valstr = 'true' if new else 'false'
        elif spec.kind == 'enum':
            new = w.currentText()
            if new == self._display_value(name, spec):
                return
            valstr = new
        else:
            new = w.value()
            cur = self._display_value(name, spec)
            if abs(new - float(cur)) < 1e-9:
                return
            valstr = repr(new)
        self._run('isolde simparam {} {}'.format(name, valstr))

    # ----------------------------------------------------------- refresh/sync
    def _param_changed_cb(self, _, data):
        key, val = data
        if key == 'forcefield':
            self._refresh_all()
        elif key in self._rows:
            self._refresh_row(key)

    def _refresh_all(self):
        sp = self.isolde.sim_params
        self._ff_label.setText(getattr(sp, 'forcefield', ''))
        for name in self._rows:
            self._refresh_row(name)

    def _refresh_row(self, name):
        from chimerax.isolde.openmm.forcefield_profiles import get_profile
        sp = self.isolde.sim_params
        spec = self._meta.spec_for(name)
        r = self._rows[name]
        w = r['editor']
        # applicability under the active force field
        applies = get_profile(getattr(sp, 'forcefield', None)).applies(name)
        for widget in r.values():
            widget.setVisible(applies)
        if not applies:
            return
        value = self._display_value(name, spec)
        w.blockSignals(True)
        try:
            if spec.kind == 'bool':
                w.setChecked(bool(value))
            elif spec.kind == 'enum':
                idx = w.findText(value)
                if idx >= 0:
                    w.setCurrentIndex(idx)
            else:
                w.setValue(float(value))
        finally:
            w.blockSignals(False)
        r['mod'].setText('●' if sp.is_overridden(name) else '')   # ●

    def cleanup(self):
        self._handler.remove()
