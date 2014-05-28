import inspect
from PyQt4.QtGui import *
from PyQt4.QtCore import *

import pypipe.tools.toolsconfig
from widgets.combobox import ComboBox
from widgets.baselistwidget import BaseListWidget
from widgets.basetablewidget import BaseTableWidget
from tablecells.tableitems import SimpleImmutableItem, TypeItem


class MethodsList(BaseListWidget):

    def __init__(self, parent=None):
        super(MethodsList, self).__init__(parent)

    def generate(self, tool):
        self.clear()
        for name, func in inspect.getmembers(tool):
            if inspect.isfunction(func):
                self.add_item(name, func)


class ArgumentsTable(BaseTableWidget):

    def __init__(self, parent=None):
        super(ArgumentsTable, self).__init__(parent)
        self.verticalHeader().setVisible(False)
        self.setColumnCount(3)

    def generate(self, func=None):
        self.clear()
        headers = ['Argument', 'Type', 'Value']
        self.setHorizontalHeaderLabels(headers)
        if func is None:
            return
        config = func()
        args = config['args']['named']
        args.update(config['args']['unnamed'])
        self.setRowCount(len(args))
        i = 0
        for name in args:
            self.setItem(i, 0, SimpleImmutableItem(name))
            type_ = args[name]
            item = TypeItem(type_)
            self.setItem(i, 1, item)
            self.set_widget(i, 2, item.get_current_type(), name)
            i += 1
        self.sortByColumn(0, Qt.AscendingOrder)
        self.link_cells()

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
        self.set_widget(i, 2, type_cell.get_current_type(), self.cellWidget(i, 2).key)


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

    def exec_(self):
        self.tools_combo.currentIndexChanged.emit(0)
        self.arguments_table.generate()
        super(AddProgramDialog, self).exec_()
