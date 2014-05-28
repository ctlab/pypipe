from PyQt4.QtGui import *

class ArgumentsListDialog(QDialog):

    def __init__(self, type_, min_len, parent=None):
        super(ArgumentsListDialog, self).__init__(parent)
        self.type_ = type_
        self.min_len = min_len
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')
        self.table = ExpandingArgumentsTable()

        self.setWindowTitle('Add arguments')

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        layout.addLayout(button_layout)
        self.setLayout(layout)

        self.ok_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.reject)
