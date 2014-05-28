from PyQt4.QtGui import *
from PyQt4.QtCore import *

import pypipe.tools.toolsconfig
from widgets.combobox import ComboBox
from widgets.listwidgets import MethodsList
from widgets.tables import ArgumentsTable


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
