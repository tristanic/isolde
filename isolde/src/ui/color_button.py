from chimerax.ui.widgets import ColorButton
from numpy import ndarray

from Qt.QtCore import Signal

#TODO: Lots of code duplication here. If wanting to specialise further, need to define a suitable base class.

class TwoColorButton(ColorButton):
    '''
    Like the ChimeraX ColorButton, but shows both colors when two colors are chosen.
    '''
    def __init__(self, *args, titles = None, direction='--', **kw):
        if direction not in '++,+0,+-,-+,-0,--,0+,00,0-'.split(','):
             raise RuntimeError(f'Direction should be one of ++, +0, +-, -+, -=, --. Got {direction}')
        self._direction = direction
        super().__init__(*args, **kw)
        if titles is None:
            titles = [None, None]
        self._titles = titles

    color_changed = Signal(ndarray, ndarray)
    def set_color(self, colors):
        from chimerax.ui.widgets.color_button import color_to_numpy_rgba8
        colors = [color_to_numpy_rgba8(c) for c in colors]
        self._set_stylesheet(colors)
        self._color = colors
        self.color_changed.emit(*self.color)
    
    def _set_stylesheet(self, colors):
        color1, color2 = [color_as_hex(c) for c in colors]
        x1,y1,x2,y2=self._get_gradient_endpoints()
        stylesheet = f'background-color: qlineargradient(x1:{x1}, y1:{y1}, x2:{x2}, y2:{y2}, stop:0 {color1}, stop:0.4 {color1}, stop:0.6 {color2}, stop:1.0 {color2});'
        self.setStyleSheet(stylesheet)


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

    def _get_gradient_endpoints(self):
        xdir, ydir = self._direction
        if xdir=='-':
            x1, x2 = 0,1
        elif xdir=='0':
            x1, x2 = 0,0
        else:
            x1, x2 = 1,0
        if ydir=='-':
            y1, y2 = 0,1
        elif ydir=='0':
            y1, y2 = 0
        else:
            y1, y2 = 1,0
        
        return x1,y1,x2,y2
    
    def changeEvent(self, event):
        if event.type() == event.Type.EnabledChange:
            if self._color is None:
                import numpy
                colors = numpy.array([[127, 127, 127, 127]]*2, numpy.uint8)
            elif self.isEnabled():
                colors = self._color
            else:
                colors = []
                for color in self._color:
                    colors.append([int((c + 218)/2) for c in color])
            self._set_stylesheet(colors)


class ThreeColorButton(ColorButton):
    '''
    Like the ChimeraX ColorButton, but for three colors.
    '''
    def __init__(self, *args, titles = None, direction='--', **kw):
        if direction not in '++,+0,+-,-+,-0,--,0+,00,0-'.split(','):
             raise RuntimeError(f'Direction should be one of ++, +0, +-, -+, -=, --. Got {direction}')
        self._direction = direction
        super().__init__(*args, **kw)
        if titles is None:
            titles = [None, None, None]
        self._titles = titles


    color_changed = Signal(ndarray, ndarray, ndarray)
    def set_color(self, colors):
        from chimerax.ui.widgets.color_button import color_to_numpy_rgba8
        colors = [color_to_numpy_rgba8(c) for c in colors]
        self._set_stylesheet(colors)
        self._color = colors
        self.color_changed.emit(*self.color)
    
    color = property(ColorButton.get_color, set_color)

    def _set_stylesheet(self, colors):
        color1, color2, color3 = [color_as_hex(c) for c in colors]
        x1,y1,x2,y2=self._get_gradient_endpoints()
        stylesheet = f'background-color: qlineargradient(x1:{x1}, y1:{y1}, x2:{x2}, y2:{y2}, stop:0 {color1}, stop:0.25 {color1}, stop: 0.4 {color2}, stop: 0.6 {color2}, stop:0.75 {color3}, stop:1.0 {color3});'
        self.setStyleSheet(stylesheet)

    def _get_gradient_endpoints(self):
        xdir, ydir = self._direction
        if xdir=='-':
            x1, x2 = 0,1
        elif xdir=='0':
            x1, x2 = 0,0
        else:
            x1, x2 = 1,0
        if ydir=='-':
            y1, y2 = 0,1
        elif ydir=='0':
            y1, y2 = 0
        else:
            y1, y2 = 1,0
        
        return x1,y1,x2,y2

    def show_color_chooser(self):
        from Qt.QtWidgets import QColorDialog
        from Qt.QtGui import QColor
        colors = []
        cd = QColorDialog(parent=self.window())
        cd.setOption(cd.ColorDialogOption.ShowAlphaChannel, self._has_alpha_channel)
        cd.setOption(cd.ColorDialogOption.NoButtons, True)
        colors.append(cd.getColor(initial=QColor(*tuple(self.color[0])), title=self._titles[0]))
        colors.append(cd.getColor(initial=QColor(*tuple(self.color[1])), title=self._titles[1]))
        colors.append(cd.getColor(initial=QColor(*tuple(self.color[2])), title=self._titles[2]))
        if not all(c.isValid() for c in colors):
            return
        self.color = colors

    def changeEvent(self, event):
        if event.type() == event.Type.EnabledChange:
            if self._color is None:
                import numpy
                colors = numpy.array([[127, 127, 127, 127]]*3, numpy.uint8)
            elif self.isEnabled():
                colors = self._color
            else:
                colors = []
                for color in self._color:
                    colors.append([int((c + 218)/2) for c in color])
            self._set_stylesheet(colors)



def color_as_hex(rgba8):
    return '#%02x%02x%02x' % tuple(rgba8[:3])
