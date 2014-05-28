from PyQt4.QtGui import *
from PyQt4.QtCore import *

from widgets.combobox import ComboBox
from pypipe.core import pipeline


class _AbstractCellWidget(QWidget):

    value_changed = pyqtSignal(int, str, bool)

    def __init__(self, key, main_widget, parent=None):
        super(_AbstractCellWidget, self).__init__(parent)
        layout = QHBoxLayout()
        layout.addWidget(main_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignCenter)
        self.setLayout(layout)
        self.indexes = []
        self.key = key

    def is_true(self):
        raise SystemExit('"is_true" method is not implemented')

    def add_receiver(self, i):
        self.indexes.append(i)

    def emit_value_changed(self):
        for i in self.indexes:
            self.value_changed.emit(i, self.key, self.is_true())


class _LineEditCellWidget(_AbstractCellWidget):

    def __init__(self, key, validator=None, parent=None):
        line_edit = QLineEdit()
        if validator:
            validator = QRegExpValidator(QRegExp(validator))
            line_edit.setValidator(validator)
        super(_LineEditCellWidget, self).__init__(key, line_edit, parent)
        line_edit.textEdited.connect(self.emit_value_changed)

    def is_true(self):
        return str(self.children()[1].text()) != ''


class FloatCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(FloatCellWidget, self).__init__(key, '[-+]?\\d*\.?\\d+([eE][-+]?\\d+)?', parent)


class IntCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(IntCellWidget, self).__init__(key, '[-+]?\\d+', parent)


class StrCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(StrCellWidget, self).__init__(key, validator=None, parent=parent)


class FileCellWidget(_AbstractCellWidget):

    def __init__(self, key, type_, parent=None):
        combo = ComboBox()
        valid_files = pipeline.get_files(type_)
        combo.add_item('None')
        for f in valid_files:
            combo.add_item(f.get_name(), f)
        super(FileCellWidget, self).__init__(key, combo, parent)
        combo.currentIndexChanged.connect(self.emit_value_changed)

    def is_true(self):
        return self.children()[1].get_current_item() is not None


class BoolCellWidget(_AbstractCellWidget):

    def __init__(self, key, parent=None):
        checkbox = QCheckBox()
        super(BoolCellWidget, self).__init__(key, checkbox, parent)
        checkbox.stateChanged.connect(self.emit_value_changed)

    def is_true(self):
        return self.children()[1].isChecked()


class ListCellWidget(_AbstractCellWidget):

    def __init__(self, key, type_, parent=None):
        button = QPushButton()
        button.setText('0 values')
        super(ListCellWidget, self).__init__(key, button, parent)

    def is_true(self):
        return True

