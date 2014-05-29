import inspect
from PyQt4.QtGui import *


class ComboBox(QComboBox):

    def __init__(self, parent=None):
        super(ComboBox, self).__init__(parent)
        self.data = []

    def add_classes_from_module(self, module):
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj):
                self.add_item(name, obj)

    def add_item(self, name, item=None):
        super(ComboBox, self).addItem(name)
        self.data.append(item)

    def get_current_item(self):
        if self.currentIndex() == -1:
            return None
        return self.data[self.currentIndex()]

    def get_item(self, i):
        return self.data[i]
