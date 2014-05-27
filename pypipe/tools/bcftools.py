from pypipe.tools.toolsconfig import Bcftools
from pypipe.utils import tool, check_if_program_exists


check_if_program_exists('bcftools')


@tool
def view():
    return Bcftools.view()


@tool
def cat():
    return Bcftools.cat()

