import sys

import formats
from pipeline import create_program


def bowtie2(x, S, U=None, r1=None, r2=None, q=None, qseq=None, f=None,
            r=None, c=None, s=None, skip=None, u=None, qupto=None,
            trim5=None, trim3=None, phred33=None, phred64=None,
            solexa_quals=None, int_quals=None, very_fast=None,
            fast=None, sensitive=None, very_sensitive=None,
            very_fast_local=None, fast_local=None, sensitive_local=None,
            very_sensitive_local=None):
    # Checking conflicts
    if (U and r1) or (U and r2) or (not U and not (r1 and r2)):
        sys.exit("bowtie2: Wrong options. Use only 'U' or 'r1' and 'r2'")
    if q and (qseq or f or r or c) or qseq and (
            f or r or c) or f and (r or c) or r and c:
        sys.exit("bowtie2: 'q', 'qseq', 'f', 'r', 'c' - this options conflict")
    if s and skip:
        sys.exit("bowtie2: use only 'skip' or 's'")
    if u and qupto:
        sys.exit("bowtie2: use only 'u' or 'qupto'")

    # Checking types
    if type(x) != formats.Bowtie2Index:
        sys.exit("bowtie2: 'x' argument must be Bowtie2Index")
    if type(S) != str:
        sys.exit("bowtie2: 'S' argument must be string")
    if s and type(s) != int:
        sys.exit("bowtie2: 's' argument must be int")
    if skip and type(skip) != int:
        sys.exit("bowtie2: 'skip' argument must be int")
    if u and type(u) != int:
        sys.exit("bowtie2: 'u' argument must be int")
    if qupto and type(qupto) != int:
        sys.exit("bowtie2: 'qupto' argument must be int")
    if trim5 and type(trim5) != int:
        sys.exit("bowtie2: 'trim5' argument must be int")
    if trim3 and type(trim3) != int:
        sys.exit("bowtie2: 'trim3' argument must be int")

    if qseq:
        read_format = formats.Qseq
        str_format = "QSEQ"
    elif f:
        read_format = formats.Fasta
        str_format = "FASTA"
    elif r:
        read_format = formats.Unknown
        str_format = "UNKNOWN"
    elif c:
        read_format = str
        str_format = "string"
    else:
        read_format = formats.Fastq
        str_format = "FASTQ"
    error_msg = "bowtie2: U argument must be array of " + str_format
    if U:
        for read in U:
            if type(read) != read_format:
                sys.exit(error_msg)
    else:
        for read in r1:
            if type(read) != read_format:
                sys.exit(error_msg)
        for read in r2:
            if type(read) != read_format:
                sys.exit(error_msg)

    # Generating basic command
    if U:
        reads = ",".join([read.path for read in U])
        cmd = ["bowtie2", "-x", x.path, "-U", reads, "-S", S]
    else:
        reads1 = ",".join([read.path for read in r1])
        reads2 = ",".join([read.path for read in r2])
        cmd = ["bowtie2", "-x", x.path, "-1", reads1, "-2", reads2, "-S", S]

    # Adding options
    if q:
        cmd.append("-q")
    elif qseq:
        cmd.append("--qseq")
    elif f:
        cmd.append("-f")
    elif r:
        cmd.append("-r")
    elif c:
        cmd.append("-c")
    if s:
        cmd.append("-s")
        cmd.append(str(s))
    elif skip:
        cmd.append("--skip")
        cmd.append(str(skip))
    if u:
        cmd.append("-u")
        cmd.append(str(u))
    elif qupto:
        cmd.append("--qupto")
        cmd.append(str(qupto))
    if trim5:
        cmd.append("--trim5")
        cmd.append(str(trim5))
    if trim3:
        cmd.append("--trim3")
        cmd.append(str(trim3))
    if phred33:
        cmd.append("--phred33")
    if phred64:
        cmd.append("--phred64")
    if solexa_quals:
        cmd.append("--solexa-quals")
    if int_quals:
        cmd.append("--int-quals")
    if very_fast:
        cmd.append("--very-fast")
    if fast:
        cmd.append("--fast")
    if sensitive:
        cmd.append("--sensitive")
    if very_sensitive:
        cmd.append("--very-sensitive")
    if very_fast_local:
        cmd.append("--very-fast-local")
    if fast_local:
        cmd.append("--fast-local")
    if sensitive_local:
        cmd.append("--sensitive-local")
    if very_sensitive_local:
        cmd.append("--very-sensitive-local")

    # Creating program
    print " ".join(cmd)
    if not c:
        files_array = U and U or r1 + r2
    else:
        files_array = []
    files_array.append(x)
    program = create_program(cmd, files_array)
    return formats.Sam(S, program)
