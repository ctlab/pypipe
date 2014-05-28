from PyQt4.QtGui import *
from PyQt4.QtCore import *


class SimpleImmutableItem(QTableWidgetItem):

    def __init__(self, value):
        super(SimpleImmutableItem, self).__init__(value)
        self.setFlags(self.flags() ^ Qt.ItemIsEditable ^ Qt.ItemIsSelectable)


class TypeItem(SimpleImmutableItem):

    def __init__(self, type_):
        super(TypeItem, self).__init__('')
        self.type_ = type_
        self.current_type = type_
        self.change_type('')

    def get_current_type(self):
        return self.current_type

    def get_type(self):
        return self.type_

    def is_mutable(self):
        return type(self.type_) == dict

    def change_type(self, key):
        if type(self.type_) == dict:
            self.current_type = self.type_[key]
        if type(self.current_type) == list:
            value = 'list of ' + self.current_type[0].__name__
        else:
            value = self.current_type.__name__
        self.setText(value)
        return value

