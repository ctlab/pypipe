import sys

import formats
from pipeline import create_program


def sort(align, output, params=None):
    if type(align) != formats.Bam:
        sys.error(align.path + " must be in BAM format");
    cmd = ["samtools", "sort"] + [output]
    if params:
        params = params.split(" ")
        try:
            params.remove("-o")
        except ValueError:
            pass
        cmd += params
    program = create_program(cmd, [align], output)
    return formats.Bam(output, program)

