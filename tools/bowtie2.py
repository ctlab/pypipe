import sys

import formats
from pipeline import create_program


def unpaired(index, reads, output, params=None):
    for read in reads:
        if type(read) != formats.Fastq
            msg = "bowtie2: " + read.path + " must be in FASTQ format"
    cmd = ["bowtie2", "-x", index, "-U"] + reads + ["-S", output] 
    # params ...
    program = create_program(cmd, reads)
    return formats.Sam(output, program)


def paired(index, reads1, reads2, output, params=None):
    for read in reads1:
        if type(read) != formats.Fastq
            msg = "bowtie2: " + read.path + " must be in FASTQ format"
    for read in reads2:
        if type(read) != formats.Fastq
            msg = "bowtie2: " + read.path + " must be in FASTQ format"
    cmd = ["bowtie2", "-x", index, "-1"] + \
        reads1 + ["-2"] + reads2 + ["-S", output] 
    # params ...
    program = create_program(cmd, reads)
    return formats.Sam(output, program)
