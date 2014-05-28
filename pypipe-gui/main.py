import sys
from PyQt4.QtGui import *

from gui.windows.mainwindow import MainWindow


app = QApplication(sys.argv)
w = MainWindow()
w.showMaximized()
app.exec_()