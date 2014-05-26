from pypipe.tools.toolsconfig import Bwa
from pypipe.utils import tool


@tool
def mem():
    return Bwa.mem()


@tool
def aln():
   return Bwa.aln()


@tool
def samse():
    return Bwa.samse()


@tool
def sampe():
    return Bwa.sampe()


@tool
def bwasw():
    return Bwa.bwasw()

