from PyQt4.QtGui import *

from windows.addfiledialog import AddFileDialog
from windows.addprogramdialog import AddProgramDialog
from widgets.listwidgets import ListWidget
from widgets.pipelineview import PipelineView

from pypipe.core import pipeline


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.files_list = ListWidget()
        self.add_file_button = QPushButton('Add file')
        self.remove_file_button = QPushButton('remove file')
        self.programs_list = ListWidget()
        self.add_program_button = QPushButton('Add program')
        self.remove_program_button = QPushButton('Remove program')
        self.pipeline_view = PipelineView()
        self.file_menu = FileMenu()
        self.pipeline_menu = PipelineMenu()

        self.menuBar().addMenu(self.file_menu)
        self.menuBar().addMenu(self.pipeline_menu)
        self.setWindowTitle('PyPipe')

        self.add_file_dialog = AddFileDialog(self)
        self.add_program_dialog = AddProgramDialog(self)

        main_layout = QGridLayout()
        panel_layout = QVBoxLayout()
        panel_layout.addWidget(QLabel('<b>Files:</b>'))
        panel_layout.addWidget(self.files_list)
        files_buttons_layout = QHBoxLayout()
        files_buttons_layout.addWidget(self.add_file_button)
        files_buttons_layout.addWidget(self.remove_file_button)
        panel_layout.addLayout(files_buttons_layout)
        panel_layout.addWidget(QLabel('<b>Programs:</b>'))
        panel_layout.addWidget(self.programs_list)
        programs_buttons_layout = QHBoxLayout()
        programs_buttons_layout.addWidget(self.add_program_button)
        programs_buttons_layout.addWidget(self.remove_program_button)
        panel_layout.addLayout(programs_buttons_layout)
        main_layout.addWidget(self.pipeline_view, 0, 0)
        main_layout.addLayout(panel_layout, 0, 1, 0, 6)
        center = QWidget()
        center.setLayout(main_layout)
        self.setCentralWidget(center)

        self.connect_all()

    def connect_all(self):
        self.add_file_button.clicked.connect(self.add_file_dialog.exec_)
        self.remove_file_button.clicked.connect(self.remove_file_from_pipeline)
        self.add_file_dialog.accepted.connect(
            lambda: self.add_new_file(self.add_file_dialog.get_file()))
        self.add_program_button.clicked.connect(self.add_program_dialog.exec_)

    def add_new_file(self, f):
        self.files_list.add_item(f.get_name() + ' (' + f.get_type() + ')', f)
        pipeline.add_file(f)
        self.pipeline_view.draw()

    def remove_file_from_pipeline(self):
        print 'TODO'
        f = self.files_list.get_current_item()
        self.files_list.remove_current_item()


class FileMenu(QMenu):

    def __init__(self):
        super(FileMenu, self).__init__('&File')

        self.save_action = QAction('Save', self)
        self.save_action.setShortcut(QKeySequence.Save)
        self.save_action.setStatusTip('Save pipeline')
        self.addAction(self.save_action)

        self.open_action = QAction('Open', self)
        self.open_action.setShortcut(QKeySequence.Open)
        self.open_action.setShortcut('Open pipeline')
        self.addAction(self.open_action)

        self.close_action = QAction('Quit', self)
        self.close_action.setShortcut(QKeySequence.Quit)
        self.close_action.setStatusTip('Close application')
        self.close_action.triggered.connect(QApplication.quit)
        self.addAction(self.close_action)


class PipelineMenu(QMenu):

    def __init__(self):
        super(PipelineMenu, self).__init__('&Pipeline')

        self.run_action = QAction('Run', self)
        self.addAction(self.run_action)

        self.reset_action = QAction('Reset', self)
        self.addAction(self.reset_action)

        self.reset_all_action = QAction('Reset all', self)
        self.addAction(self.reset_all_action)
