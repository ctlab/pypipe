from pypipe.tools.toolsconfig import Samtools
from pypipe.utils import tool, check_if_program_exists


check_if_program_exists('samtools')


@tool
def view():
    return Samtools.view()

@tool
def mpileup():
    return Samtools.mpileup()


@tool
def cat():
    return Samtools.cat()


@tool
def sort():
    return Samtools.sort()


@tool
def merge():
    return Samtools.merge()


@tool
def rmdup():
    return Samtools.rmdup()


@tool
def calmd():
    return Samtools.calmd()

