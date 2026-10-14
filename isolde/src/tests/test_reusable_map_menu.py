"""Pure Python map-menu regressions, using Qt and model API doubles.

python3 -m unittest discover -s isolde/src/tests -p test_reusable_map_menu.py
"""
import importlib.util
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


class Volume:
    def __init__(self, data, name='map', id_string='1'):
        self.data, self.name, self.id_string = data, name, id_string


class NXmapHandler(Volume):
    def __init__(self, data, owner, **kw):
        super().__init__(data, **kw)
        self.mapset = SimpleNamespace(structure=owner)


def module(name, **attrs):
    m = ModuleType(name)
    m.__dict__.update(attrs)
    return m


class ReusableMapMenuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ui_base = {name: object for name in ('UI_Panel_Base', 'DefaultHLayout',
            'DefaultVLayout', 'DefaultSpacerItem', 'QDoubleSpinBox')}
        widgets = {name: object for name in ('QPushButton', 'QLabel', 'QMenu', 'QFileDialog')}
        cls.modules = patch.dict(sys.modules, {
            'Qt': module('Qt'), 'Qt.QtWidgets': module('Qt.QtWidgets', **widgets),
            'chimerax': module('chimerax'),
            'chimerax.isolde': module('chimerax.isolde'),
            'chimerax.isolde.ui': module('chimerax.isolde.ui'),
            'chimerax.isolde.ui.general_tab': module('general_tab'),
            'chimerax.isolde.ui.ui_base': module('ui_base', **ui_base),
            'chimerax.isolde.ui.collapse_button': module('collapse_button', CollapsibleArea=object),
            'chimerax.map': module('map', Volume=Volume),
            'chimerax.clipper': module('clipper'),
            'chimerax.clipper.maps': module('maps', NXmapHandler=NXmapHandler),
            'chimerax.clipper.symmetry': module('symmetry',
                get_map_mgr=lambda target, create=False: target.mgr),
        })
        cls.modules.start()
        path = Path(__file__).resolve().parents[1] / 'ui/general_tab/add_maps.py'
        spec = importlib.util.spec_from_file_location('chimerax.isolde.ui.general_tab.add_maps', path)
        cls.ui = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.ui)

    @classmethod
    def tearDownClass(cls):
        cls.modules.stop()

    def setUp(self):
        self.target = SimpleNamespace(id_string='2', name='target', mgr=None)
        self.owner = SimpleNamespace(id_string='3', name='owner', mgr=None)
        self.volumes = []
        self.session = SimpleNamespace(models=SimpleNamespace(list=lambda: self.volumes))

    def test_shared_grid_shown_once_with_owner_and_managed_source(self):
        grid = object()
        ordinary = Volume(grid)
        managed = NXmapHandler(grid, self.owner, id_string='3.1.1')
        self.volumes = [ordinary, managed]
        entries = self.ui._available_real_space_maps(self.session, self.target)
        self.assertEqual(entries, [[managed, ['#3: owner']]])

    def test_switching_model_excludes_only_its_own_grids(self):
        grid = object()
        managed = NXmapHandler(grid, self.owner)
        self.owner.mgr = SimpleNamespace(all_non_xtal_maps=[managed])
        self.volumes = [managed, Volume(grid), Volume(object(), name='unused')]
        self.assertEqual(len(self.ui._available_real_space_maps(self.session, self.target)), 2)
        entries = self.ui._available_real_space_maps(self.session, self.owner)
        self.assertEqual([v.name for v, _ in entries], ['unused'])
        self.assertEqual(self.ui._available_real_space_maps(self.session, None), [])

    def test_crystallographic_subclasses_not_offered(self):
        class CrystalVolume(Volume):
            pass
        self.volumes = [CrystalVolume(object()), Volume(object())]
        self.assertEqual(len(self.ui._available_real_space_maps(self.session, self.target)), 1)

    def test_menu_command_uses_offered_map_and_captured_target(self):
        managed = NXmapHandler(object(), self.owner, id_string='3.1.1')
        self.volumes = [managed]
        actions, commands = [], []
        class Menu:
            def clear(self):
                actions.clear()

            def addAction(self, label):
                action = SimpleNamespace(label=label)
                action.triggered = SimpleNamespace(connect=lambda cb: setattr(action, 'callback', cb))
                actions.append(action)
                return action
        dialog = SimpleNamespace(real_space_map_menu=Menu(), session=self.session,
            isolde=SimpleNamespace(selected_model=self.target))
        self.ui.MapAddDialog._populate_available_volumes_menu(dialog)
        self.assertIn('used by #3: owner', actions[0].label)
        dialog.isolde.selected_model = self.owner
        cx_commands = module('commands', run=lambda session, command: commands.append(command))
        with patch.dict(sys.modules, {'chimerax.core': module('core'),
                'chimerax.core.commands': cx_commands}):
            actions[0].callback()
        self.assertEqual(commands, ['clipper assoc #3.1.1 to #2'])


if __name__ == '__main__':
    unittest.main()
