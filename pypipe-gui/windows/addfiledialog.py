from PyQt4.QtGui import *

import pypipe.formats
import pypipe.basefile
from pypipe.core import pipeline
from widgets.combobox import ComboBox


class AddFileDialog(QDialog):

    def __init__(self, parent=None):
        super(AddFileDialog, self).__init__(parent)
        self.formats_combo = ComboBox()
        self.filename_edit = QLineEdit()
        self.open_button = QPushButton('Open')
        self.ok_button = QPushButton('&OK')
        self.cancel_button = QPushButton('&Cancel')

        self.setWindowTitle('Add file')
        top_layout = QVBoxLayout()
        top_layout.addWidget(QLabel('<b>File format:</b>'))
        top_layout.addWidget(self.formats_combo)
        top_layout.addWidget(QLabel('<b>File Name:</b>'))
        center_layout = QHBoxLayout()
        center_layout.addWidget(self.filename_edit)
        center_layout.addWidget(self.open_button)
        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(self.ok_button)
        bottom_layout.addWidget(self.cancel_button)
        layout = QVBoxLayout()
        layout.addLayout(top_layout)
        layout.addLayout(center_layout)
        layout.addLayout(bottom_layout)
        self.setLayout(layout)

        self.formats_combo.add_classes_from_module(pypipe.formats)
        self.connect_all()

    def connect_all(self):
        self.cancel_button.clicked.connect(self.reject)
        self.filename_edit.textChanged.connect(self.turn_ok_button)
        self.formats_combo.currentIndexChanged.connect(self.turn_ok_button)
        self.ok_button.clicked.connect(self.accept)
        self.open_button.clicked.connect(self.open_file)

    def turn_ok_button(self):
        try:
            f = self.get_file()
            self.ok_button.setEnabled(True)
        except pypipe.basefile.FileNotExistsError:
            self.ok_button.setEnabled(False)
            return
        if pypipe.core.pipeline.can_add_file(f):
            self.ok_button.setEnabled(True)
        else:
            self.ok_button.setEnabled(False)

    def open_file(self):
        file_name = QFileDialog.getOpenFileName(self, 'Open file')
        self.filename_edit.setText(file_name)

    def get_file(self):
        init = self.formats_combo.get_current_item()
        path = str(self.filename_edit.text())
        return init(path)

    def exec_(self):
        self.turn_ok_button()
        super(AddFileDialog, self).exec_()
