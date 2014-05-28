from PyQt4.QtGui import *

from widgets.basetablewidget import BaseTableWidget


class ExpandingArgumentsTable(BaseTableWidget):

    def __init__(self, parent=None):
        super(ExpandingArgumentsTable, self).__init__(parent)


class ArgumentsListDialog(QDialog):

    def __init__(self, parent=None):
        super(ArgumentsListDialog, self).__init__(parent)
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
