import sys

import formats
from pipeline import create_program


def index(ref, params=None):
    if type(ref) != formats.Fasta:
        msg = "bwa index: " + ref.path + " must be in FASTA format"
        sys.exit(msg)
    cmd = ["bwa", "index", ref.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref])


def bwasw(ref, read, output, params=None):
    if type(ref) != formats.Fasta:
        msg = "bwa bwasw: " + ref.path + " must be in FASTA format"
        sys.exit(msg)
    if type(read) != formats.Fastq:
        msg = "bwa bwasw: " + read.path + " must be in FASTQ format"
        sys.exit(msg)
    cmd = ["bwa", "bwasw", ref.path, read.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref, read], output)
    return formats.Sam(output, program)


def samse(ref, index, read, output, params=None):
    if type(ref) != formats.Fasta:
        msg = "bwa samse: " + ref.path + " must be in FASTA format"
        sys.exit(msg)
    if type(index) != formats.Sai:
        msg = "bwa samse: " + index.path + " must be in SAI format"
        sys.exit(msg)
    if type(read) != formats.Fastq:
        msg = "bwa samse: " + read.path + " must be in FASTQ format"
        sys.exit(msg)
    cmd = ["bwa", "samse", ref.path, index.path, read.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref, index, read], output)
    return formats.Sam(output, program)

