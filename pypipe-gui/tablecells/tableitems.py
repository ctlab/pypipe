from PyQt4.QtGui import *
from PyQt4.QtCore import *


class SimpleImmutableItem(QTableWidgetItem):

    def __init__(self, value):
        super(SimpleImmutableItem, self).__init__(value)
        self.setFlags(self.flags() ^ Qt.ItemIsEditable ^ Qt.ItemIsSelectable)


class TypeItem(SimpleImmutableItem):

    def __init__(self, type_):
        if type(type_) == dict:
            self.current_type = type_['']
        else:
            self.current_type = type_
        if type(self.current_type) == list:
            value = 'list of ' + self.current_type[0].__name__
        else:
            value = self.current_type.__name__
        super(TypeItem, self).__init__(value)
        self.type_ = type_

    def get_current_type(self):
        return self.current_type

    def get_type(self):
        return self.type_

    def is_mutable(self):
        return type(self.type_) == dict

    def change_type(self, key):
        self.current_type = self.type_[key]
        self.setText(self.current_type.__name__)
