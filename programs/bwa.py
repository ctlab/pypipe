import sys

import formats
from pipeline import create_program


def index(ref, params=None):
    if type(ref) != formats.Fasta:
        sys.exit(ref.path + " must be in FASTA format")
    cmd = ["bwa", "index", ref.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref])


def bwasw(ref, read, output, params=None):
    if type(ref) != formats.Fasta:
        sys.exit(ref.path + " must be in FASTA format")
    if type(read) != formats.Fastq:
        sys.exit(read.path + " must be in FASTQ format")
    cmd = ["bwa", "bwasw", ref.path, read.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref, read], output)
    return formats.Sam(output, program)

