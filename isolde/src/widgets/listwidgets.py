
from chimerax.ui.widgets import ModelListWidget, ModelMenuButton
from chimerax.map import Volume

class EligibleVolumeListWidget(ModelListWidget):
    @staticmethod
    def _filter_func(item):
        from chimerax.clipper.maps import MapHandlerBase
        if isinstance(item, MapHandlerBase):
            return False
        return True
    def __init__(self, session, **kw):
        super().__init__(session, class_filter=Volume, filter_func=self._filter_func, **kw)

class EligibleVolumeMenuButton(ModelMenuButton):
    @staticmethod
    def _filter_func(item):
        from chimerax.clipper.maps import MapHandlerBase
        if isinstance(item, MapHandlerBase):
            return False
        return True
    def __init__(self, session, **kw):
        super().__init__(session, class_filter=Volume, filter_func=self._filter_func, **kw)
    

            