from chimerax.ui.widgets import ColorButton
from numpy import ndarray

from Qt.QtCore import Signal
class TwoColorButton(ColorButton):
    '''
    Like the ChimeraX ColorButton, but shows both colors when two colors are chosen.
    '''
    def __init__(self, *args, titles = None, **kw):
        super().__init__(*args, **kw)
        if titles is None:
            titles = [None, None]
        self._titles = titles

    color_changed = Signal(ndarray, ndarray)
    def set_color(self, colors):
        from chimerax.ui.widgets.color_button import color_to_numpy_rgba8
        colors = [color_to_numpy_rgba8(c) for c in colors]
        color1, color2 = [color_as_hex(c) for c in colors]
        stylesheet = f'background-color: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 {color1}, stop:0.4 {color1}, stop:0.6 {color2}, stop:1.0 {color2});'
        self.setStyleSheet(stylesheet)
        self._color = colors
        self.color_changed.emit(*self.color)
    
    color = property(ColorButton.get_color, set_color)

    def show_color_chooser(self):
        from Qt.QtWidgets import QColorDialog
        from Qt.QtGui import QColor
        colors = []
        cd = QColorDialog(parent=self.window())
        cd.setOption(cd.ColorDialogOption.ShowAlphaChannel, self._has_alpha_channel)
        cd.setOption(cd.ColorDialogOption.NoButtons, True)
        colors.append(cd.getColor(initial=QColor(*tuple(self.color[0])), title=self._titles[0]))
        colors.append(cd.getColor(initial=QColor(*tuple(self.color[1])), title=self._titles[1]))
        if not all(c.isValid() for c in colors):
            return
        self.color = colors


def color_as_hex(rgba8):
    return '#%02x%02x%02x' % tuple(rgba8[:3])
