import inspect
from PyQt4.QtGui import *

import formats
from gui.widgets.combobox import ComboBox
from pipeline import pipeline


class AddFileDialog(QDialog):

    def __init__(self, parent=None):
        super(AddFileDialog, self).__init__(parent)
        self.formats_combo = ComboBox()
        self.filename_edit = QLineEdit()
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')

        self.setWindowTitle('Add file')
        layout = QGridLayout()
        layout.addWidget(QLabel('<b>File format:</b>'), 0, 0)
        layout.addWidget(self.formats_combo, 0, 1)
        layout.addWidget(QLabel('<b>File Name:</b>'), 1, 0)
        layout.addWidget(self.filename_edit, 1, 1)
        layout.addWidget(self.ok_button, 2, 0)
        layout.addWidget(self.cancel_button, 2, 1)
        self.setLayout(layout)

        self.formats_combo.add_classes_from_module(formats)
        self.connect_all()

    def connect_all(self):
        self.cancel_button.clicked.connect(self.reject)
        self.filename_edit.textEdited.connect(self.turn_ok_button)
        self.formats_combo.currentIndexChanged.connect(self.turn_ok_button)
        self.ok_button.clicked.connect(self.accept)

    def turn_ok_button(self):
        f = self.get_file()
        if pipeline.can_add_file(f):
            self.ok_button.setEnabled(True)
        else:
            self.ok_button.setEnabled(False)

    def get_file(self):
        init = self.formats_combo.get_current_item()
        path = str(self.filename_edit.text())
        return init(path)

    def exec_(self):
        self.turn_ok_button()
        super(AddFileDialog, self).exec_()

