import os
import numpy

from ..util import slot_disconnected
from ..collapse_button import CollapsibleArea
from ..ui_base import (
    UI_Panel_Base, 
    DefaultHLayout, DefaultVLayout, DefaultHLayout, DefaultSpacerItem,
    QComboBox, QDoubleSpinBox
)

from Qt.QtWidgets import (
    QLabel, QToolButton, QPushButton, QCheckBox, QWidget, 
    QTreeWidget, QTreeWidgetItem, QPushButton, QMenu, QAbstractItemView,
    QSizePolicy, QHeaderView, QItemDelegate
)

from Qt.QtGui import QIcon, QFontMetrics

from Qt.QtCore import Qt, QSize, QRect, QRectF

from .. import icon_dir, DEFAULT_ICON_SIZE

class XmapLiveSettingsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Dynamic Crystallographic Map Settings', **kwargs)
        msd = self.content = XmapLiveSettingsDialog(session, isolde, gui, self)
        self.setContentLayout(msd.main_layout)

class XmapStaticSettingsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Precalculated Crystallographic Map Settings', **kwargs)
        msd = self.content = XmapStaticSettingsDialog(session, isolde, gui, self)
        self.setContentLayout(msd.main_layout)

class NXmapSettingsPanel(CollapsibleArea):
    def __init__(self, session, isolde, parent, gui, **kwargs):
        super().__init__(gui, parent, title='Non-crystallographic Map Settings', **kwargs)
        msd = self.content = NXmapSettingsDialog(session, isolde, gui, self)
        self.setContentLayout(msd.main_layout)


class MapSettingsDialog(UI_Panel_Base):
    '''
    Base class for different map settings dialogues. Not to be instantiated directly.
    '''
    MAP_TYPE=None
    MAP_SET_TYPE=None
    MAP_NAME_PREFIX=None

    MAX_COLUMN_WIDTH=200

    NAME_COLUMN=    0
    ID_COLUMN=      1
    DISPLAY_COLUMN= 2
    STYLE_COLUMN=   3
    COLOR_COLUMN=   4
    DIFF_COLUMN=    5
    MDFF_COLUMN=    6
    WEIGHT_COLUMN=  7

    _labels = {
        NAME_COLUMN:    'Name',
        ID_COLUMN:      'ID',
        DISPLAY_COLUMN: 'Show',
        STYLE_COLUMN:   'Style',
        COLOR_COLUMN:   'Colour',
        DIFF_COLUMN:    'Diff Map?',
        MDFF_COLUMN:    'MDFF?',
        WEIGHT_COLUMN:  'Weight'
    }

    _style_icons = {
        'solid':     os.path.join(icon_dir, 'mapsurf.png'),
        'mesh':    os.path.join(icon_dir, 'mesh.png'),
        'transparent': os.path.join(icon_dir, 'icecube.png')
    }

    ICON_SIZE = (16,16)
    def __init__(self, session, isolde, gui, collapse_area, sim_sensitive=False):
        super().__init__(session, isolde, gui, collapse_area.content_area, sim_sensitive=sim_sensitive)
        import os
        self.container = collapse_area
        mf = self.main_frame
        self._current_maps = []
        ml = self.main_layout = DefaultHLayout()
        self._first_rebuild=True

        self.tree = tree = QTreeWidget(mf)
        ml.addWidget(tree)
        tree.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        tree.setHeader(WordWrapHeader(Qt.Orientation.Horizontal))
        tree.setHeaderLabels(self._labels.values())
        tree.header().setMinimumSectionSize(20)
        tree.setSelectionBehavior(QAbstractItemView.SelectRows)
        tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        tree.setAnimated(True)
        tree.setUniformRowHeights(True)
        tree.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tree.expanded.connect(lambda *_: tree.resizeColumnToContents(self.ID_COLUMN))
        tree.itemClicked.connect(self._row_clicked_cb)

        from chimerax.core.models import (
            ADD_MODELS, REMOVE_MODELS, MODEL_ID_CHANGED,
            MODEL_NAME_CHANGED
        )
        cxh = self._chimerax_trigger_handlers
        for trigger_name in (ADD_MODELS, REMOVE_MODELS, MODEL_ID_CHANGED,
            MODEL_NAME_CHANGED):
            def _model_changes_cb(trigger_name, data):
                # Need to delay the rebuild until the next frame, otherwise the new surfaces may not yet exist
                from chimerax.atomic import get_triggers
                def _new_frame_cb(*_, tname=trigger_name, d=data):
                    self._rebuild_tree_if_necessary(tname, data)
                    from chimerax.core.triggerset import DEREGISTER
                    return DEREGISTER
                self.session.triggers.add_handler('new frame', _new_frame_cb)           
            cxh.append(session.triggers.add_handler(trigger_name, _model_changes_cb))
        
        ith = self._isolde_trigger_handlers
        ith.append(isolde.triggers.add_handler(isolde.SELECTED_MODEL_CHANGED, self._rebuild_tree_if_necessary))
        self.container.expanded.connect(self.rebuild_tree)
        self._temporary_handlers = []
        self.rebuild_tree()


    def _rebuild_tree_if_necessary(self, trigger_name, data):
        from chimerax.core.models import MODEL_NAME_CHANGED, MODEL_ID_CHANGED
        if trigger_name in (MODEL_NAME_CHANGED, MODEL_ID_CHANGED):
            if self.isolde.selected_model is not None:
                if data in self.isolde.selected_model.parent.all_models():
                    return self.rebuild_tree()
        rebuild_needed = self._compare_maps_to_current()
        if rebuild_needed:
            th = self._temporary_handlers
            while len(th):
                th.pop().remove()
        if self.container.is_collapsed:
            return
        if rebuild_needed:
            self.rebuild_tree()

    def init_tree_structure(self, mapsets):
        '''
        Override in derived classes. Should return a dict of {parent TreeWidgetItem: child models}
        '''
        pass


    def rebuild_tree(self):
        th = self._temporary_handlers
        while len(th):
            th.pop().remove()

        tree = self.tree
        tree.clear()
        if not len(self._current_maps):
            return
        from chimerax.clipper import get_map_mgr
        mgr = get_map_mgr(self.isolde.selected_model)
        mapsets = [c for c in mgr.child_models() if isinstance(c, self.MAP_SET_TYPE)]
        if len(mapsets):
            tree_dict = self.init_tree_structure(mapsets)
            for name, (p, maps) in tree_dict.items():
                for map in maps:
                    self.add_tree_entry(p, map, strip_name_prefix=self.MAP_NAME_PREFIX)                    
        tree.expandAll()
        if self._first_rebuild:
            for column in self._labels.keys():
                tree.resizeColumnToContents(column)
                if tree.columnWidth(column) > self.MAX_COLUMN_WIDTH:
                    tree.setColumnWidth(column, self.MAX_COLUMN_WIDTH)
            self._first_rebuild=False

    def _row_clicked_cb(self, item, column):
        v = getattr(item, '_volume', None)
        if v is not None:
            from chimerax.core.commands import run
            run(self.session, f'sel #{v.id_string}', log=False)

    def add_tree_entry(self, parent, map, strip_name_prefix=None):
        ''' Override in derived classes '''
        pass
    
    def add_display_menu_checkbox(self, item, column, map):
        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        cb = QCheckBox()
        cb.setChecked(map.display)
        def set_display(flag):
            map.display = flag
        cb.toggled.connect(set_display)
        th = self._temporary_handlers
        from chimerax.core.models import MODEL_DISPLAY_CHANGED
        from ..util import slot_disconnected
        def display_changed_cb(*_):
            with slot_disconnected(cb.toggled, set_display):
                cb.setChecked(map.display)
        th.append(self.session.triggers.add_handler(MODEL_DISPLAY_CHANGED, display_changed_cb))
        l.addWidget(cb)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)

    def add_style_menu_button(self, item, column, map):
        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        b = QPushButton()
        menu = QMenu(b)
        b.setIconSize(QSize(*self.ICON_SIZE))
        b.setFixedSize(QSize(*self.ICON_SIZE))
        b.setIcon(self._get_icon_for_volume(map))
        b.setStyleSheet('QPushButton::menu-indicator {width:0px;}')
        b.setMenu(menu)
        a = menu.addAction(QIcon(self._style_icons['solid']), 'Opaque surface')
        def set_solid(*_, button=b, volume=map):
            self._solid_map_button_cb(button, volume)
        a.triggered.connect(set_solid)
        a = menu.addAction(QIcon(self._style_icons['transparent']), 'Transparent surface')
        def set_transparent(*_, button=b, volume=map):
            self._transparent_map_button_cb(button, volume)
        a.triggered.connect(set_transparent)
        a = menu.addAction(QIcon(self._style_icons['mesh']), 'Mesh')
        def set_mesh(*_, button=b, volume=map):
            self._mesh_map_button_cb(button, volume)
        a.triggered.connect(set_mesh)
        l.addWidget(b)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)
    
    def add_color_menu_button(self, item, column, map):
        from chimerax.ui.widgets import ColorButton
        from ..color_button import TwoColorButton
        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        # Difference maps have negative contour first, but show positive first in 
        # picker.
        if map.is_difference_map:
            b = TwoColorButton(w, max_size=self.ICON_SIZE, has_alpha_channel=False,
                titles=('Positive contour', 'Negative contour'))
            b.color = [s.color for s in map.surfaces][::-1]
        else:
            b = ColorButton(w, max_size=self.ICON_SIZE, has_alpha_channel=False)
            b.color = map.surfaces[0].color
        def _color_changed_cb(*args, v=map):
            colors = [*args][::-1]
            # Keep existing transparency values
            existing_colors = [s.color for s in v.surfaces]
            for i, e in enumerate(existing_colors):
                colors[i][-1] = e[-1]
            for s, c in zip(v.surfaces, colors):
                s.color = c
        b.color_changed.connect(_color_changed_cb)
        l.addWidget(b)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)

    def add_difference_map_checkbox(self, item, column, map):
        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        cb = QCheckBox(w)
        l.addWidget(cb)
        cb.setChecked(map.is_difference_map)
        l.addStretch()
        def dcb_checked_cb(checked, i=item, v=map):
            v.is_difference_map = checked
            # Need to replace the color picker button
            self.tree.removeItemWidget(i, self.COLOR_COLUMN)
            self.add_color_menu_button(i, self.COLOR_COLUMN, v)
        cb.toggled.connect(dcb_checked_cb)
        self.tree.setItemWidget(item, column, w)

    def add_mdff_checkbox(self, item, column, map, always_allow_mdff=False):
        from chimerax.isolde import session_extensions as sx
        mgr = sx.get_mdff_mgr(self.isolde.selected_model, map, create=always_allow_mdff)
        if mgr is None:
            cb = QLabel('X')
            cb.setStyleSheet('color:red; font-weight:bold;')
            cb.setToolTip('This map is not suitable for MDFF')
        else:
            cb = QCheckBox()
            cb.setChecked(mgr.enabled)
            cb.setToolTip('<span>Allow this map to "pull" on atoms (cannot be set during simulations)</span>')
            def enable_mdff(flag, m=mgr):
                m.enabled=flag
            cb.toggled.connect(enable_mdff)
            self._temporary_handlers.append(self.isolde.triggers.add_handler(self.isolde.SIMULATION_STARTED,
                lambda *_: cb.setEnabled(False)))
            self._temporary_handlers.append(self.isolde.triggers.add_handler(self.isolde.SIMULATION_TERMINATED,
                lambda *_: cb.setEnabled(True)))

        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        l.addWidget(cb)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)

    def add_weight_spinbox(self, item, column, map):
        from chimerax.isolde import session_extensions as sx
        mgr = sx.get_mdff_mgr(self.isolde.selected_model, map)
        if mgr is None:
            return
        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        sb = MapWeightSpinBox()
        sb.setValue(mgr.global_k)
        def set_weight(value, m=mgr):
            mgr.global_k = value
        sb.valueChanged.connect(set_weight)
        l.addWidget(sb)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)

    def _get_icon_for_volume(self, v):
        s = v.surfaces[0]
        style = s.style
        if style=='surface':
            if s.color[-1] < 255:
                style='transparent'
            else:
                style='solid'
        return QIcon(self._style_icons[style])

    
    def _compare_maps_to_current(self):
        m = self.isolde.selected_model
        rebuild_needed = False
        if m is None or m.was_deleted:
            if len(self._current_maps):
                rebuild_needed = True
                self._current_maps = []
            return rebuild_needed
        from chimerax.clipper import get_map_mgr
        mmgr = get_map_mgr(m)
        if (mmgr is None or not len(mmgr.all_maps)):
            if len(self._current_maps):
                rebuild_needed = True
                self._current_maps = []
            return rebuild_needed
        all_maps = list(sorted(mmgr.all_maps, key=lambda m: m.id))
        filtered_maps = [m for m in all_maps if isinstance(m, self.MAP_TYPE)]
        if filtered_maps != self._current_maps:
            rebuild_needed = True
            self._current_maps = filtered_maps      
        return rebuild_needed
            
    def _solid_map_button_cb(self, button, v):
        import os
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [v], **map_style_settings[map_styles.solid_opaque])
        button.setIcon(QIcon(self._style_icons['solid']))

    def _transparent_map_button_cb(self, button, v):
        import os
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        if v.is_difference_map:
            style = map_styles.solid_t40
        else:
            style = map_styles.solid_t20
        volumecommand.volume(self.session, [v], **map_style_settings[style])
        button.setIcon(QIcon(self._style_icons['transparent']))

    def _mesh_map_button_cb(self, button, v):
        import os
        from chimerax.isolde.visualisation import map_styles, map_style_settings
        from chimerax.map import volumecommand
        volumecommand.volume(self.session, [v], **map_style_settings[map_styles.mesh_triangle])
        button.setIcon(QIcon(self._style_icons['mesh']))
    
    def cleanup(self):
        for h in self._temporary_handlers:
            h.remove()
        super().cleanup()

class XmapSettingsDialogBase(MapSettingsDialog):
    def init_tree_structure(self, mapsets):
        ret = {}
        tree = self.tree
        parent = tree.invisibleRootItem()
        for mapset in mapsets:
            maps = [map for map in mapset.all_maps if isinstance(map, self.MAP_TYPE)]
            if not len(maps):
                continue
            ms = QTreeWidgetItem(parent)
            ms.setText(self.NAME_COLUMN, mapset.base_name)
            ms.setText(self.ID_COLUMN, mapset.id_string)
            self.add_display_menu_checkbox(ms, self.DISPLAY_COLUMN, mapset)
            ret[mapset.name] = (ms, maps)
        return ret            

    def add_tree_entry(self, parent, map, strip_name_prefix=None):
        import os
        item = QTreeWidgetItem(parent)
        item._volume = map
        if strip_name_prefix is not None:
            name = map.name.replace(strip_name_prefix, '')
        else:
            name = map.name
        item.setText(self.NAME_COLUMN, name)
        item.setText(self.ID_COLUMN, map.id_string)
        self.add_display_menu_checkbox(item, self.DISPLAY_COLUMN, map)
        self.add_style_menu_button(item, self.STYLE_COLUMN, map)
        self.add_color_menu_button(item, self.COLOR_COLUMN, map)
        self.add_difference_map_checkbox(item, self.DIFF_COLUMN, map)
        self.add_mdff_checkbox(item, self.MDFF_COLUMN, map)
        self.add_weight_spinbox(item, self.WEIGHT_COLUMN, map)


class XmapLiveSettingsDialog(XmapSettingsDialogBase):
    from chimerax.clipper.maps import XmapHandler_Live, XmapSet
    MAP_TYPE=XmapHandler_Live
    MAP_SET_TYPE=XmapSet
    MAP_NAME_PREFIX='(LIVE) '

        

class XmapStaticSettingsDialog(XmapSettingsDialogBase):
    from chimerax.clipper.maps import XmapHandler_Static, XmapSet
    MAP_TYPE=XmapHandler_Static
    MAP_SET_TYPE=XmapSet
    MAP_NAME_PREFIX='(STATIC) '

    dont_ask_again=False

    def add_mdff_checkbox(self, item, column, map):
        from chimerax.isolde import session_extensions as sx
        mgr = sx.get_mdff_mgr(self.isolde.selected_model, map)
        if mgr is None:
            cb = QLabel('X')
            cb.setStyleSheet('color:red; font-weight:bold;')
            cb.setToolTip('This map is not suitable for MDFF')
        else:
            cb = QCheckBox()
            cb.setChecked(mgr.enabled)
            cb.setToolTip('<span>Allow this map to "pull" on atoms (cannot be set during simulations)</span>')
            def enable_mdff(flag, m=mgr):
                if flag and not self.dont_ask_again:
                    from ...dialog import choice_warning
                    go, dag = choice_warning('Since this map was generated from precalculated '
                    'amplitudes and phases, ISOLDE has no way of determining '
                    'whether it is suitable for fitting simulations. If generated '
                    'from crystallographic data, you should ensure that it is a'
                    '2Fo-Fc (or similar) map calculated with the free reflections '
                    'excluded. Are you sure you want to continue?', allow_dont_ask_again=True)
                    if not go:
                        flag = False
                        with slot_disconnected(cb.toggled, enable_mdff):
                            cb.setChecked(False)
                    self.dont_ask_again=dag
                m.enabled=flag
                
            cb.toggled.connect(enable_mdff)
            self._temporary_handlers.append(self.isolde.triggers.add_handler(self.isolde.SIMULATION_STARTED,
                lambda *_: cb.setEnabled(False)))
            self._temporary_handlers.append(self.isolde.triggers.add_handler(self.isolde.SIMULATION_TERMINATED,
                lambda *_: cb.setEnabled(True)))

        w = QWidget()
        l = DefaultHLayout()
        w.setLayout(l)
        l.addStretch()
        l.addWidget(cb)
        l.addStretch()
        self.tree.setItemWidget(item, column, w)



class NXmapSettingsDialog(MapSettingsDialog):
    from chimerax.clipper.maps import NXmapHandler, NXmapSet
    MAP_TYPE=NXmapHandler
    MAP_SET_TYPE=NXmapSet
    MAP_NAME_PREFIX=''
    
    def init_tree_structure(self, mapsets):
        # Should just be one mapset according to current design, but let's future-proof
        tree = self.tree
        ret = {}
        if not any(len(mapset.all_maps) for mapset in mapsets):
            return ret
        parent = tree.invisibleRootItem()
        if len(mapsets)==1:
            return {'null': (parent, mapsets[0].all_maps)}
        for mapset in mapsets:
            if not len(mapset.all_maps):
                continue
            ms = QTreeWidgetItem(parent)
            ms.setText(self.NAME_COLUMN, mapset.name)
            ms.setText(self.ID_COLUMN, mapset.id_string)
            self.add_display_menu_checkbox(ms, self.DISPLAY_COLUMN, mapset)
            ret[mapset.name] = (ms, mapset.all_maps)
        return ret
        
    def add_tree_entry(self, parent, map, strip_name_prefix=None):
        import os
        item = QTreeWidgetItem(parent)
        item._volume = map
        if strip_name_prefix is not None:
            name = map.name.replace(strip_name_prefix, '')
        else:
            name = map.name
        item.setText(self.NAME_COLUMN, name)
        item.setText(self.ID_COLUMN, map.id_string)
        self.add_display_menu_checkbox(item, self.DISPLAY_COLUMN, map)
        self.add_style_menu_button(item, self.STYLE_COLUMN, map)
        self.add_color_menu_button(item, self.COLOR_COLUMN, map)
        self.add_difference_map_checkbox(item, self.DIFF_COLUMN, map)
        self.add_mdff_checkbox(item, self.MDFF_COLUMN, map, always_allow_mdff=True)
        self.add_weight_spinbox(item, self.WEIGHT_COLUMN, map)




class WordWrapHeader(QHeaderView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter|Qt.AlignmentFlag(Qt.TextFlag.TextWordWrap))
    
    def sectionSizeFromContents(self, logicalIndex):
        if self.model():
            headerText = self.model().headerData(logicalIndex, self.orientation(),Qt.ItemDataRole.DisplayRole)
            words = headerText.split()
            if len(words) == 1:
                return super().sectionSizeFromContents(logicalIndex)
            else:
                longest_word = max(words, key=lambda w:len(w))
                metrics = self.fontMetrics()
                max_width = metrics.boundingRect(longest_word).width()
                rect = metrics.boundingRect(QRect(0,0,max_width,5000),
                    self.defaultAlignment(), headerText)
                buffer = QSize(3, 3)
                return rect.size() + buffer
        else:
            return super().sectionSizeFromContents(logicalIndex)
    

class MapWeightSpinBox(QDoubleSpinBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setKeyboardTracking(False)
        self.setMinimumWidth(40)
        self.setMinimumHeight(20)
        self.setMaximum(1e5-1)
        self.setMinimum(1e-6)
        self.valueChanged.connect(self.update_display_and_step)

    def _displayed_decimal_places_and_step(self, number, sig_figs=3):
        from math import log, ceil, floor
        if number <= 0:
            return 2, 0.1
        places = max(ceil(-log(number, 10)+sig_figs-1), 0)
        step = 10**(floor(log(number, 10)-1))
        return places, step

    def update_display_and_step(self, *_):
        v = self.value()
        dps, step = self._displayed_decimal_places_and_step(v)
        self.setDecimals(dps)
        self.setSingleStep(step)
    

