import pypipe.basefile
from tablecells.tablewidgets import *
from tablecells.tableitems import SimpleImmutableItem, TypeItem


class ArgumentsTable(QTableWidget):

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
            self.set_widget(i, item.get_current_type(), name)
            i += 1
        self.sortByColumn(0, Qt.AscendingOrder)
        self.link_cells()

    def set_widget(self, i, type_, key):
        if type(type_) == list:
            self.setCellWidget(i, 2, ListCellWidget(key, type_))
        elif type_ == bool:
            self.setCellWidget(i, 2, BoolCellWidget(key))
        elif type_ == float:
            self.setCellWidget(i, 2, FloatCellWidget(key))
        elif type_ == int:
            self.setCellWidget(i, 2, IntCellWidget(key))
        elif type_ == str:
            self.setCellWidget(i, 2, StrCellWidget(key))
        elif issubclass(type_, pypipe.basefile.File):
            self.setCellWidget(i, 2, FileCellWidget(key, type_))
        else:
            raise SystemExit('Unknown argument type')

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
        self.set_widget(i, type_cell.get_current_type(), self.cellWidget(i, 2).key)
