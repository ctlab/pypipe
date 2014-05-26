from pypipe.tools.toolsconfig import Bcftools
from pypipe.utils import install_program, tool


@tool
def view():
    return Bcftools.view()


@tool
def cat():
    return Bcftools.cat()

