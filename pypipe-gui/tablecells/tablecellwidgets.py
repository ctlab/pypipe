import ntpath
from PyQt4.QtGui import *
from PyQt4.QtCore import *

from widgets.combobox import ComboBox
from pypipe.core import pipeline


class _BaseCellWidget(QWidget):

    value_changed = pyqtSignal(int, str, bool)
    value_changed2 = pyqtSignal()

    def __init__(self, key, main_widget, parent=None):
        super(_BaseCellWidget, self).__init__(parent)
        layout = QHBoxLayout()
        layout.addWidget(main_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignCenter)
        self.setLayout(layout)
        self.indexes = []
        self.key = key

    def is_true(self):
        raise SystemExit('"is_true" method is not implemented')

    def get_value(self):
        raise SystemExit('"get_value" method is not implemented')

    def get_real_value(self):
        raise SystemExit('"get_real_value" method is not implemented')

    def set_value(self, v):
        raise SystemExit('"set_value" method is not implemented')

    def add_receiver(self, i):
        self.indexes.append(i)

    def emit_value_changed(self):
        for i in self.indexes:
            self.value_changed.emit(i, self.key, self.is_true())

    def get(self):
        return self.children()[1]


class _LineEditCellWidget(_BaseCellWidget):

    def __init__(self, key, validator=None, parent=None):
        line_edit = QLineEdit()
        if validator:
            validator = QRegExpValidator(QRegExp(validator))
            line_edit.setValidator(validator)
        super(_LineEditCellWidget, self).__init__(key, line_edit, parent)
        line_edit.textEdited.connect(self.emit_value_changed)
        line_edit.textEdited.connect(self.value_changed2)

    def set_value(self, v):
        self.get().setText(v)

    def get_value(self):
        return self.get().text()

    def is_true(self):
        return self.get_value() != ''


class FloatCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(FloatCellWidget, self).__init__(key, '[-+]?\\d*\.?\\d+([eE][-+]?\\d+)?', parent)

    def get_real_value(self):
        return float(self.get_value())


class IntCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(IntCellWidget, self).__init__(key, '[-+]?\\d+', parent)

    def get_real_value(self):
        return int(self.get_value())


class StrCellWidget(_LineEditCellWidget):

    def __init__(self, key, parent=None):
        super(StrCellWidget, self).__init__(key, validator=None, parent=parent)

    def get_real_value(self):
        return str(self.get_value())


class FileCellWidget(_BaseCellWidget):

    def __init__(self, key, type_, parent=None):
        combo = ComboBox()
        valid_files = pipeline.get_files(type_)
        combo.add_item('None')
        for f in valid_files:
            combo.add_item(f.get_name(), f)
        super(FileCellWidget, self).__init__(key, combo, parent)
        combo.currentIndexChanged.connect(self.emit_value_changed)
        combo.currentIndexChanged.connect(self.value_changed2)

    def get_value(self):
        return self.get().currentIndex()

    def get_real_value(self):
        return self.get().get_current_item()

    def set_value(self, v):
        self.get().setCurrentIndex(v)

    def is_true(self):
        return self.get_real_value() is not None


class BoolCellWidget(_BaseCellWidget):

    def __init__(self, key, parent=None):
        checkbox = QCheckBox()
        super(BoolCellWidget, self).__init__(key, checkbox, parent)
        checkbox.stateChanged.connect(self.emit_value_changed)
        checkbox.stateChanged.connect(self.value_changed2)

    def get_value(self):
        return self.is_true()

    def get_real_value(self):
        return self.get_value()

    def set_value(self, v):
        self.get().setEnabled(v)

    def is_true(self):
        return self.get().isChecked()


class ListCellWidget(_BaseCellWidget):

    def __init__(self, key, type_, dialog, parent=None):
        self.button = QPushButton()
        self.button.setText('0 values')
        super(ListCellWidget, self).__init__(key, self.button, parent)
        self.type_ = type_[0]
        self.min_len = type_[1]
        self.dialog = dialog
        self.values_len = 0

        self.button.clicked.connect(self.dialog.exec_)
        self.dialog.accepted.connect(lambda: self.set_text(len(self.dialog.get_real_values())))
        self.dialog.accepted.connect(self.value_changed2)

    def set_text(self, v):
        self.values_len = v
        self.button.setText(str(self.values_len) + ' values')

    def get_value(self):
        return self.dialog.get_real_values()

    def get_real_value(self):
        return self.dialog.get_real_values()

    def is_true(self):
        return self.values_len >= self.min_len


class OutFileCellWidget(_BaseCellWidget):

    def __init__(self, key, parent=None):
        self.button = QPushButton()
        self.button.setText('None')
        super(OutFileCellWidget, self).__init__(key, self.button, parent)
        self.value = None

        self.button.clicked.connect(self.set_file)

    def set_file(self):
        file_name = QFileDialog.getSaveFileName(self, 'Output file')
        if file_name != '':
            self.value = ntpath.basename(str(file_name))
            self.button.setText(self.value)
            self.value_changed2.emit()

    def get_value(self):
        return self.value

    def get_real_value(self):
        return self.value

    def is_true(self):
        return self.value is not None
