from PyQt4.QtGui import *


class _ListItem(QListWidgetItem):

    def __init__(self, value):
        super(_ListItem, self).__init__()
        self._value = value

    def get_value(self):
        return self._value


class BaseListWidget(QListWidget):

    def __init__(self, parent=None):
        super(BaseListWidget, self).__init__(parent)

    def add_item(self, name, value):
        item = _ListItem(value)
        item.setText(name)
        self.addItem(item)

    def get_current_item(self):
        if self.currentItem() is None:
            return None
        return self.currentItem().get_value()

    def remove_current_item(self):
        item = self.currentItem()
        if item:
            self.takeItem(self.row(item))
