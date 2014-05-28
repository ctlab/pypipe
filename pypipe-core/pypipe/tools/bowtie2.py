from pypipe.tools.toolsconfig import Bowtie2
from pypipe.utils import tool, check_if_program_exists


check_if_program_exists('bowtie2')


@tool
def bowtie2():
    return Bowtie2.bowtie2()

