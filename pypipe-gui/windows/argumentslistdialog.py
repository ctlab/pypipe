from PyQt4.QtGui import *
from PyQt4.QtCore import *


from widgets.basetablewidget import BaseTableWidget


class ExpandingArgumentsTable(BaseTableWidget):

    value_changed = pyqtSignal()

    def __init__(self, type_, min_len, parent=None):
        super(ExpandingArgumentsTable, self).__init__(parent)
        self.setColumnCount(1)
        headers = ['Value']
        self.setHorizontalHeaderLabels(headers)
        self.type_ = type_
        self.min_len = min_len
        for i in xrange(0, min_len):
            self.expand()

    def expand(self):
        self.insertRow(self.rowCount())
        self.set_widget(self.rowCount() - 1, 0, self.type_)
        self.cellWidget(self.rowCount() - 1, 0).value_changed2.connect(self.value_changed)

    def get_real_values(self):
        values = []
        for i in xrange(0, self.rowCount()):
            w = self.cellWidget(i, 0)
            if w.is_true():
                values.append(w.get_real_value())
        return values

    def get_values(self):
        values = []
        for i in xrange(0, self.rowCount()):
            w = self.cellWidget(i, 0)
            values.append(w.get_value())
        return values


class ArgumentsListDialog(QDialog):

    def __init__(self, type_, parent=None):
        super(ArgumentsListDialog, self).__init__(parent)
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')
        self.plus_button = QPushButton('+')
        self.type_ = type_[0]
        self.min_len = type_[1]
        self.table = ExpandingArgumentsTable(self.type_, self.min_len)
        self.widget = None
        self.old_values = None

        self.setWindowTitle('Add arguments')

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        plus_layout = QHBoxLayout()
        plus_layout.addWidget(self.plus_button)
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        layout.addLayout(plus_layout)
        layout.addLayout(button_layout)
        self.setLayout(layout)

        self.ok_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.reject)
        self.plus_button.clicked.connect(self.table.expand)
        self.rejected.connect(self.set_old_values)
        self.table.value_changed.connect(self.turn_ok_button)

    def turn_ok_button(self):
        l = len(self.get_real_values())
        if l >= self.min_len or l == 0:
            self.ok_button.setEnabled(True)
        else:
            self.ok_button.setEnabled(False)

    def set_old_values(self):
        for i in xrange(0, self.table.rowCount()):
            w = self.table.cellWidget(i, 0)
            w.set_value(self.old_values[i])

    def exec_(self):
        self.turn_ok_button()
        self.old_values = self.table.get_values()
        super(ArgumentsListDialog, self).exec_()

    def get_real_values(self):
        return self.table.get_real_values()