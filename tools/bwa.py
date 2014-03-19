import sys

import formats
from pipeline import create_program


def index(ref, params=None):
    """ bwa index 

    For example: bwa index -a bwtsw ref.fasta
    bwa.index(ref=Fasta("ref.fasta"), params="-a bwtsw")

    This function returns nothing

    """
    if type(ref) != formats.Fasta:
        msg = "bwa index: " + ref.path + " must be in FASTA format"
        sys.exit(msg)
    cmd = ["bwa", "index", ref.path]
    if params:
        params = params.split(" ")
        cmd += params
    program = create_program(cmd, [ref])


def bwasw(ref, read, output, params=None):
    """ bwa bwasw 

    For example: bwa bwasw -t 2 ref.fa read.fq > aln.sam
    aln = bwa.bwasw(Fasta("ref.fa")), Fastq("read.fq"), "aln.sam", "-t 2")

    This function returns SAM file

    """
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
    """ bwa samse 

    For example: bwa samse ref.fa aln_sa.sai read.fq > aln.sam
    ref = Fasta("ref.fa")
    aln_sa = Sai("aln_sa.sai")
    read = Fastq("read.fq")
    aln_sam = bwa.samse(ref=ref, index=aln_sa, read=read, output="aln.sam")

    This function returns SAM file

    """
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

