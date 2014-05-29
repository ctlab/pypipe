import tempfile
import math
import threading
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtSvg import *

from windows.addfiledialog import AddFileDialog
from windows.addprogramdialog import AddProgramDialog
from widgets.baselistwidget import BaseListWidget

from pypipe.core import pipeline


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.files_list = BaseListWidget()
        self.add_file_button = QPushButton('Add file')
        self.rename_file_button = QPushButton('Rename file')
        self.rename_file_button.setEnabled(False)
        self.programs_list = BaseListWidget()
        self.add_program_button = QPushButton('Add program')
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
        files_buttons_layout.addWidget(self.rename_file_button)
        panel_layout.addLayout(files_buttons_layout)
        panel_layout.addWidget(QLabel('<b>Programs:</b>'))
        panel_layout.addWidget(self.programs_list)
        programs_buttons_layout = QHBoxLayout()
        programs_buttons_layout.addWidget(self.add_program_button)
        panel_layout.addLayout(programs_buttons_layout)
        main_layout.addWidget(self.pipeline_view, 0, 0)
        main_layout.addLayout(panel_layout, 0, 1, 0, 6)
        center = QWidget()
        center.setLayout(main_layout)
        self.setCentralWidget(center)

        self.connect_all()

    def connect_all(self):
        self.add_file_button.clicked.connect(self.add_file_dialog.exec_)
        self.rename_file_button.clicked.connect(self.rename_file)
        self.add_file_dialog.accepted.connect(
            lambda: self.add_new_file(self.add_file_dialog.get_file()))
        self.add_program_button.clicked.connect(self.add_program_dialog.exec_)
        self.add_program_dialog.accepted.connect(self.add_new_program)

        self.pipeline_menu.run_action.triggered.connect(self.run_pipeline)
        self.programs_list.currentItemChanged.connect(lambda: self.pipeline_menu.turn_actions(
            self.programs_list.get_current_item()))
        self.pipeline_menu.reset_action.triggered.connect(self.reset_pipeline)
        self.pipeline_menu.reset_all_action.triggered.connect(self.reset_all_pipeline)

        self.file_menu.save_action.triggered.connect(self.save_pipeline)
        self.file_menu.open_action.triggered.connect(self.open_pipeline)

    def save_pipeline(self):
        file_name = QFileDialog.getSaveFileName(self, 'Save pipeline')
        pipeline.save(str(file_name))

    def open_pipeline(self):
        file_name = QFileDialog.getOpenFileName(self, 'Open pipeline')
        try:
            pipeline.load(str(file_name))
        except:
            return
        self.files_list.clear()
        self.programs_list.clear()
        for f in pipeline.files:
            if f.next_programs:
                self.files_list.add_item(f.get_name() + ' (' + f.get_type() + ')' + \
                                 ' (' + str(f.number + 1) + ')', f)
        for p in pipeline.all_programs:
            self.programs_list.add_item(p.name + ' (' + str(p.number + 1) + ')', p)
        self.pipeline_view.draw()

    def run_pipeline(self):
        thread = threading.Thread(target=pipeline.run,
                                  args=(self.programs_list.get_current_item().number,
                                  self.pipeline_view.img_file.name))
        thread.start()

    def reset_pipeline(self):
        pipeline.reset(self.programs_list.get_current_item().number)
        self.pipeline_view.draw()

    def reset_all_pipeline(self):
        pipeline.reset_all()
        self.pipeline_view.draw()

    def add_new_file(self, f):
        pipeline.add_file(f)
        self.files_list.add_item(f.get_name() + ' (' + f.get_type() + ')' + \
                                 ' (' + str(f.number + 1) + ')', f)

    def add_new_program(self):
        p = pipeline.all_programs[-1]
        self.programs_list.add_item(p.name + ' (' + str(p.number + 1) + ')', p)
        self.pipeline_view.draw()

    def rename_file(self):
        print 'TODO'


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
        self.run_action.setEnabled(False)
        self.addAction(self.run_action)

        self.reset_action = QAction('Reset', self)
        self.reset_action.setEnabled(False)
        self.addAction(self.reset_action)

        self.reset_all_action = QAction('Reset all', self)
        self.addAction(self.reset_all_action)

    def turn_actions(self, o):
        if o is not None:
            self.run_action.setEnabled(True)
            self.reset_action.setEnabled(True)
        else:
            self.run_action.setEnabled(False)
            self.reset_action.setEnabled(False)


class PipelineView(QGraphicsView):

    def __init__(self, parent=None):
        super(PipelineView, self).__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.img_file = tempfile.NamedTemporaryFile()
        watcher = QFileSystemWatcher(self)
        watcher.addPath(self.img_file.name)
        watcher.fileChanged.connect(self.update_img)

    def update_img(self):
        self.scene().clear()
        img = QGraphicsSvgItem(self.img_file.name)
        self.scene().addItem(img)

    def draw(self):
        pipeline.draw(self.img_file.name)

    def wheelEvent(self, e):
        factor = math.pow(1.2, e.delta() / 480.0)
        self.scale(factor, factor)
        e.accept()
