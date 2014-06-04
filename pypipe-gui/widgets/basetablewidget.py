import pypipe.basefile
from tablecells.tablecellwidgets import *


class BaseTableWidget(QTableWidget):

    def __init__(self, parent=None):
        super(BaseTableWidget, self).__init__(parent)
        self.verticalHeader().setVisible(False)

    def set_widget(self, i, j, type_, key=None, dialog=None):
        if type_ == 'out':
            self.setCellWidget(i, j, OutFileCellWidget(key))
        elif type(type_) == list:
            self.setCellWidget(i, j, ListCellWidget(key, type_, dialog))
        elif type_ == bool:
            self.setCellWidget(i, j, BoolCellWidget(key))
        elif type_ == float:
            self.setCellWidget(i, j, FloatCellWidget(key))
        elif type_ == int:
            self.setCellWidget(i, j, IntCellWidget(key))
        elif type_ == str:
            self.setCellWidget(i, j, StrCellWidget(key))
        elif issubclass(type_, pypipe.basefile.File):
            self.setCellWidget(i, j, FileCellWidget(key, type_))
        else:
            raise SystemExit('Unknown argument type')
