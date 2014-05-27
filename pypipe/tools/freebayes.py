from pypipe.tools.toolsconfig import Freebayes
from pypipe.utils import tool, check_if_program_exists


check_if_program_exists('freebayes')


@tool
def freebayes():
    return Freebayes.freebayes()

