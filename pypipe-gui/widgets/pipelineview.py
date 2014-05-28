import tempfile
import math
from PyQt4.QtGui import *
from PyQt4.QtSvg import *

from pypipe.core import pipeline


class PipelineView(QGraphicsView):

    def __init__(self, parent=None):
        super(PipelineView, self).__init__(parent)
        self.img_file = tempfile.mktemp()
        self.setScene(QGraphicsScene(self))

    def draw(self):
        pipeline.draw(self.img_file)
        self.scene().clear()
        img = QGraphicsSvgItem(self.img_file)
        self.scene().addItem(img)

    def wheelEvent(self, e):
        factor = math.pow(1.2, e.delta() / 600.0)
        self.scale(factor, factor)
        e.accept()
