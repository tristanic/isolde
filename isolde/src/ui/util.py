
from PyQt5.QtWidgets import QComboBox, QApplication

def win_auto_resize_combo_box_hack(combobox):
    import sys
    if 'win' in sys.platform.lower():
        combobox.__class__ = WinAutoResizeQComboBox

class WinAutoResizeQComboBox(QComboBox):
    '''
    Due to a long-standing (~2011) bug in Qt5 in Windows, the drop-down menu
    from a QComboBox doesn't automatically resize to fit its contents. So, we
    have to do this ugly hack to force it.
    '''
    def _init(self):
        self._width_update_needed=True

    def showPopup(self):
        if self._width_update_needed:
            fix_combo_box_view_width(self)
            self._width_update_needed = False
        super().showPopup()

    def addItem(self, *args, **kwargs):
        super().addItem(*args, **kwargs)
        self._width_update_needed=True

    def addItems(self, texts):
        super().addItems(texts)
        self._width_update_needed=True

    def insertItem(self, *args, **kwargs):
        super().insertItem(*args, **kwargs)
        self._width_update_needed=True

    def insertItems(self, index, ilist):
        super().insertItems(index, ilist)
        self._width_update_needed=True

    def removeItem(self, index):
        super().removeItem(index)
        self._width_update_needed=True

    def setItemText(self, index, text):
        super().setItemText(index, text)
        self._width_update_needed=True

    def clear(self):
        super().clear()
        self._width_update_needed=True



def fix_combo_box_view_width(combobox):
    cb = combobox
    num_entries = cb.count()
    if num_entries <= cb.maxVisibleItems():
        scroll_width = 0
    else:
        style = QApplication.style()
        scroll_width = style.pixelMetric(style.PixelMetric.PM_ScrollBarExtent)

    v = cb.view()
    fm = v.fontMetrics()
    max_width = 0
    for i in range(num_entries):
        width = fm.width(cb.itemText(i))
        if max_width < width:
            max_width = width

    v.setMinimumWidth(max_width+scroll_width)
