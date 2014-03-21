import sys

import formats
from pipeline import create_program


def bowtie2(x, S, U=None, r1=None, r2=None):
    if (U and r1) or (U and r2) or (not U and not (r1 and r2)):
        sys.exit("bowtie2: Wrong arguments. Use only 'U' or only 'r1' 'r2'")
    if type(x) != formats.Bowtie2Index:
        sys.exit("bowtie2: x argument must be Bowtie2Index")
    if type(S) != str:
        sys.exit("bowtie2: S argument must be string")
    if U:
        if type(U) != formats.Fastq:
            sys.exit("bowtie2: U argument must be in FASTQ format")
    else:
        if type(r1) != formats.Fastq:
            sys.exit("bowtie2: r1 argument must be in FASTQ format")
        if type(r2) != formats.Fastq:
            sys.exit("bowtie2: r2 argument must be in FASTQ format")
    if U:
        cmd = ["bowtie2", "-x", x.path, "-U", U.path, "-S", S]
    else:
        cmd = ["bowtie2", "-x", x.path, "-1", r1.path, "-2", r2.path, "-S", S]
    if U:
        program = create_program(cmd, [U])
    else:
        program = create_program(cmd, [r1, r2])
    return formats.Sam(S, program)
