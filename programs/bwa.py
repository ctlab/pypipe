import sys

import formats
from pipeline import create_program


def bwasw(ref, read, output):
    if type(ref) != formats.Fasta:
        sys.exit(ref.path + " must be in FASTA format")
    if type(read) != formats.Fastq:
        sys.exit(read.path + " must be in FASTQ format")
    cmd = ["bwa", "bwasw", ref.path, read.path]
    program = create_program(cmd, [ref, read], output)
    return formats.Sam(output, program)

