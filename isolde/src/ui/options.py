from chimerax.ui.options import Option

class RGBA8PaletteOption(Option):
    '''
    Widget with a variable number of color choosers useful for defining custom palettes.

    Supported keyword arguments:

        * num colors: The number of individual color pickers
        * labels: a list of length num_colors or None
        * initial_colors: a list of RGBA8 values of length num_colors or None
    '''
    default_color = [200,200,200,255]
    
    def __init__(self, *args, num_colors=1, **kw):
        self.num_colors = num_colors
        super().__init__(*args, **kw)
    
    def get_value(self):
        colors = [self._color_button[i].color for i in self.num_colors]
    
    def set_value(self, value):
        for i, val in enumerate(value):
            self._color_button[i].color = val
    
    value = property(get_value, set_value)

    def _make_widget(self, **kw):
        nc = self.num_colors
        from chimerax.ui.widgets import MultiColorButton
        from Qt.QtWidgets import QWidget, QHBoxLayout, QLabel
        labels = kw.pop('labels', None)
        if labels is None:
            labels = [None] + ["  "]*(nc-1)
        self.widget = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(0,0,0,0)
        self._color_button = []
        initial_colors = kw.get('initial_colors', [default_color]*nc)
        for i in range(nc):
            label = labels[i]
            if label:
                layout.addWidget(QLabel(label))
            mcb = MultiColorButton(max_size=(16,16), has_alpha_channel=True)
            self._color_button.append(mcb)
            mcb.color = initial_colors[i]
            mcb.color_changed.connect(lambda c, s=self: s.make_callback())
            layout.addWidget(mcb)
        self.widget.setLayout(layout)



            

