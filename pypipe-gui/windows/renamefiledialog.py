from addfiledialog import AddFileDialog


class RenameFileDialog(AddFileDialog):

    def __init__(self, parent=None):
        super(RenameFileDialog, self).__init__(parent)
        self.formats_combo.setEnabled(False)

    def exec_(self, file_):
        type_ = file_.__class__
        i = 0
        for t in self.formats_combo.data:
            if type_ == t:
                self.formats_combo.setCurrentIndex(i)
                break
            i += 1
        super(RenameFileDialog, self).exec_()