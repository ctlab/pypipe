from PyQt4.QtGui import *

from gui.windows.addfiledialog import AddFileDialog
from gui.windows.addprogramdialog import AddProgramDialog
from gui.widgets.listwidgets import ListWidget
from pipeline import pipeline


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.files_list = ListWidget()
        self.add_file_button = QPushButton('Add file')
        self.remove_file_button = QPushButton('remove file')
        self.programs_list = ListWidget()
        self.add_program_button = QPushButton('Add program')
        self.remove_program_button = QPushButton('Remove program')
        self.pipeline_view = QGraphicsView()
        self.menu = MainMenu()

        self.menuBar().addMenu(self.menu)
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
            lambda: self.add_file_to_pipeline(self.add_file_dialog.get_file()))
        self.add_program_button.clicked.connect(self.add_program_dialog.exec_)

    def add_file_to_pipeline(self, f):
        self.files_list.add_item(f.get_str(), f)
        pipeline.add_file(f)

    def remove_file_from_pipeline(self):
        f = self.files_list.get_current_item()
        self.files_list.remove_current_item()
        pipeline.remove_file(f)


class MainMenu(QMenu):

    def __init__(self):
        super(MainMenu, self).__init__()

        self.save_action = QAction('Save', self)
        self.save_action.setShortcut(QKeySequence.Save)
        self.save_action.setStatusTip('Save pipeline')
        #self.save_action.triggered.connect()
        self.addAction(self.save_action)

        self.open_action = QAction('Open', self)
        self.open_action.setShortcut(QKeySequence.Open)
        self.open_action.setShortcut('Open pipeline')
        #self.open_action.triggered.connect()
        self.addAction(self.open_action)

        self.close_action = QAction('Close', self)
        self.close_action.setShortcut(QKeySequence.Close)
        self.close_action.setStatusTip('Close application')
        #self.close_action.triggered.connect()
        self.addAction(self.close_action)