from pypipe.tools.toolsconfig import Bowtie2
from pypipe.utils import tool


@tool
def bowtie2():
    return Bowtie2.bowtie2()

