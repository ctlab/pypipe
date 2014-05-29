import inspect
from PyQt4.QtGui import *
from PyQt4.QtCore import *

import pypipe.tools.toolsconfig
import pypipe.utils
from widgets.combobox import ComboBox
from widgets.baselistwidget import BaseListWidget
from widgets.basetablewidget import BaseTableWidget
from tablecells.tableitems import SimpleImmutableItem, TypeItem
from windows.argumentslistdialog import ArgumentsListDialog


class MethodsList(BaseListWidget):

    def __init__(self, parent=None):
        super(MethodsList, self).__init__(parent)

    def generate(self, tool):
        self.clear()
        for name, func in inspect.getmembers(tool):
            if inspect.isfunction(func):
                self.add_item(name, func)


class ArgumentsTable(BaseTableWidget):

    value_changed = pyqtSignal()

    def __init__(self, parent=None):
        super(ArgumentsTable, self).__init__(parent)
        self.setColumnCount(3)
        self.mandatory = []
        self.config = None

    def generate(self, func=None):
        self.clear()
        headers = ['Argument', 'Type', 'Value']
        self.setHorizontalHeaderLabels(headers)
        if func is None:
            return
        self.config = func()
        args = {}
        named = self.config['args']['named']
        for k in named:
            args[k] = named[k]
        unnamed = self.config['args']['unnamed']
        for k, v in unnamed:
            args[k] = v
        args[self.config['log']] = str
        self.setRowCount(len(args))
        i = 0
        for name in args:
            self.setItem(i, 0, SimpleImmutableItem(name))
            type_ = args[name]
            item = TypeItem(type_)
            self.setItem(i, 1, item)
            t = item.get_current_type()
            if type(t) == list:
                self.set_widget(i, 2, t, name, ArgumentsListDialog(t))
            else:
                self.set_widget(i, 2, t, name)
            self.cellWidget(i, 2).value_changed2.connect(self.value_changed)
            i += 1
        self.sortByColumn(0, Qt.AscendingOrder)
        self.link_cells()
        for i in xrange(0, self.rowCount()):
            key = str(self.item(i, 0).text())
            if key[-1] == '*':
                self.mandatory.append(i)

    def link_cells(self):
        mutable_types = [i for i in xrange(0, self.rowCount()) if self.item(i, 1).is_mutable()]
        for i in mutable_types:
            type_item = self.item(i, 1)
            for j in xrange(0, self.rowCount()):
                name = str(self.item(j, 0).text())
                if name in type_item.get_type():
                    widget = self.cellWidget(j, 2)
                    widget.add_receiver(i)
                    widget.value_changed.connect(self.change_type)

    def change_type(self, i, k, flag):
        k = str(k)
        if flag:
            key = k
        else:
            key = ''
        type_cell = self.item(i, 1)
        type_cell.change_type(key)
        t = type_cell.get_current_type()
        if type(t) == list:
            self.set_widget(i, 2, t, self.cellWidget(i, 2).key, ArgumentsListDialog(t))
        else:
            self.set_widget(i, 2, t, self.cellWidget(i, 2).key)

    def all_mandatory_true(self):
        for i in self.mandatory:
            if not self.cellWidget(i, 2).is_true():
                return False
        return True

    def get_parameters(self):
        params = {}
        for i in xrange(0, self.rowCount()):
            w = self.cellWidget(i, 2)
            if w.is_true():
                k = str(self.item(i, 0).text())
                k = k.replace('-', '_').replace('*', '')
                params[k] = w.get_real_value()
        return params


class AddProgramDialog(QDialog):

    def __init__(self, parent=None):
        super(AddProgramDialog, self).__init__(parent)
        self.tools_combo = ComboBox()
        self.methods_list = MethodsList()
        self.arguments_table = ArgumentsTable()
        self.ok_button = QPushButton()
        self.ok_button.setText('OK')
        self.cancel_button = QPushButton()
        self.cancel_button.setText('Cancel')
        self.setWindowTitle('Add program')

        layout = QVBoxLayout()
        top_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel('<b>Tool:</b>'))
        left_layout.addWidget(self.tools_combo)
        left_layout.addWidget(self.methods_list)
        right_layout = QVBoxLayout()
        right_layout.addWidget(QLabel('<b>Parameters:</b>'))
        right_layout.addWidget(self.arguments_table)
        top_layout.addLayout(left_layout)
        top_layout.addLayout(right_layout)
        layout.addLayout(top_layout)
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.ok_button)
        buttons_layout.addWidget(self.cancel_button)
        layout.addLayout(buttons_layout)
        self.setLayout(layout)

        self.tools_combo.add_classes_from_module(pypipe.tools.toolsconfig)
        self.connect_all()

    def connect_all(self):
        self.cancel_button.clicked.connect(self.reject)
        self.ok_button.clicked.connect(self.accept)
        self.tools_combo.currentIndexChanged.connect(
            lambda: self.methods_list.generate(self.tools_combo.get_current_item()))
        self.methods_list.currentItemChanged.connect(
            lambda: self.arguments_table.generate(self.methods_list.get_current_item()))
        self.arguments_table.value_changed.connect(self.turn_ok_button)

    def turn_ok_button(self):
        if self.arguments_table.all_mandatory_true():
            self.ok_button.setEnabled(True)
        else:
            self.ok_button.setEnabled(False)

    def accept(self):
        program_func = pypipe.utils.tool(lambda: self.arguments_table.config)
        program_func(**(self.arguments_table.get_parameters()))
        super(AddProgramDialog, self).accept()

    def exec_(self):
        self.ok_button.setEnabled(False)
        self.tools_combo.currentIndexChanged.emit(0)
        self.arguments_table.generate()
        super(AddProgramDialog, self).exec_()
