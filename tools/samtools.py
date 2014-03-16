import sys

import formats
from pipeline import create_program


def index(align):
    if type(align) != formats.Bam:
        msg = "samtools index: " + align.path + " must be in BAM format"
        sys.error(msg)
    cmd = ["samtools", "index", align]
    program = create_program(cmd, [align])


def sort(align, output, params=None):
    if type(align) != formats.Bam:
        msg = "samtools sort: " + align.path + " must be in BAM format"
        sys.error(msg)
    cmd = ["samtools", "sort", align.path, output]
    if params:
        params = params.split(" ")
        params = filter(lambda p: p != "-o", params)
        params = map(lambda p: p.replace("o", ""), params)
        cmd += params
    program = create_program(cmd, [align])
    return formats.Bam(output, program)

