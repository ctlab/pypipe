from pypipe.tools.toolsconfig import Freebayes
from pypipe.utils import tool


@tool
def freebayes():
    return Freebayes.freebayes()

