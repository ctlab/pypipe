from pypipe.tools.toolsconfig import Varscan
from pypipe.utils import tool, check_if_program_exists


check_if_program_exists('VarScan.jar')


@tool
def pileup2snp():
    return Varscan.pileup2cns()


@tool
def pileup2indel():
    return Varscan.pileup2indel()


@tool
def pileup2cns():
    return Varscan.pileup2cns()


@tool
def mpileup2snp():
    return Varscan.mpileup2snp()


@tool
def mpileup2indel():
    return Varscan.mpileup2indel()


@tool
def mpileup2cns():
    return Varscan.mpileup2cns()


@tool
def somatic():
    return Varscan.somatic()

