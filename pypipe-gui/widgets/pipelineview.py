import tempfile
import math
from PyQt4.QtGui import *
from PyQt4.QtSvg import *
from PyQt4.QtCore import *

from pypipe.core import pipeline


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
