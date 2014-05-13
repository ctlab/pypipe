import sys

from pypipe import formats
from pypipe.utils import *


install_program("samtools.sh", "samtools")


def view(in_, o, b=None, f=None, F=None, h=None, H=None, l=None, q=None,
         r=None, R=None, S=None, c=None, t=None, u=None, log=None,
         regions=None):
    program = pipeline.add_node("samtools view", log)
    if S:
        program.add_arg(in_, formats.Sam)
    else:
        program.add_arg(in_, formats.Bam)
    program.add_arg(o, str, "-o")
    program.add_arg(b, bool, "-b")
    program.add_arg(f, int, "-f")
    program.add_arg(F, int, "-F")
    program.add_arg(h, bool, "-h")
    program.add_arg(H, bool, "-H")
    program.add_arg(l, str, "-l")
    program.add_arg(q, int, "-q")
    program.add_arg(r, str, "-r")
    program.add_arg(R, formats.Bed, "-R")
    program.add_arg(S, bool, "-S")
    program.add_arg(c, bool, "-c")
    program.add_arg(t, formats.TextFile, "-t")
    program.add_arg(u, bool, "-u")
    program.add_args(regions, str)
    if b:
        return formats.Bam(o, program)
    else:
        return formats.Sam(o, program)


def mpileup(in_, out, _6=None, A=None, B=None, b=None, C=None, d=None, E=None,
            f=None, l=None, q=None, Q=None, r=None, D=None, g=None, S=None,
            u=None, e=None, h=None, I=None, L=None, o=None, P=None, log=None):
    program = pipeline.add_node("samtools mpileup", log, out)
    program.add_args(in_, formats.Bam)
    program.add_arg(_6, bool, "-6")
    program.add_arg(A, bool, "-A")
    program.add_arg(B, bool, "-B")
    program.add_arg(b, formats.TextFile, "-b")
    program.add_arg(C, int, "-C")
    program.add_arg(d, int, "-d")
    program.add_arg(E, bool, "-E")
    program.add_arg(f, formats.Fasta, "-f")
    program.add_arg(l, formats.Bed, "-l")
    program.add_arg(q, int, "-q")
    program.add_arg(Q, int, "-Q")
    program.add_arg(r, str, "-r")
    program.add_arg(D, bool, "-D")
    program.add_arg(g, bool, "-g")
    program.add_arg(S, bool, "-S")
    program.add_arg(u, bool, "-u")
    program.add_arg(e, int, "-e")
    program.add_arg(h, int, "-h")
    program.add_arg(I, bool, "-I")
    program.add_arg(L, int, "-L")
    program.add_arg(o, int, "-o")
    program.add_arg(P, str, "-P")
    if u or g:
        return formats.Bcf(out, program)
    else:
        return formats.Pileup(out, program)


def reheader(in_header, in_, log=None):
    program = pipeline.add_node("samtools reheader", log)
    program.add_arg(in_header, formats.Sam)
    program.add_arg(in_, formats.Bam)
    return formats.Bam(in_.path, program)


def cat(in_, o, h=None, log=None):
    program = pipeline.add_node("samtools cat", log)
    program.add_arg(h, formats.Sam, "-h")
    program.add_arg(o, str, "-o")
    program.add_args(in_, formats.Bam)
    return formats.Bam(o, program)
    

def sort(in_, out, n=None, m=None, log=None):
    program = pipeline.add_node("samtools sort", log)
    program.add_arg(in_, formats.Bam)
    program.add_arg(out, str)
    program.add_arg(n, bool, "-n")
    program.add_arg(m, int, "-m")
    return formats.Bam(out + ".bam", program)


def merge(in_, out, _1=None, f=None, h=None, n=None,
          R=None, r=None, u=None, log=None):
    program = pipeline.add_node("samtools merge", log)
    program.add_arg(out, str)
    program.add_args(in_, formats.Bam, 2)
    program.add_arg(_1, bool, "-1")
    program.add_arg(f, bool, "-f")
    program.add_arg(h, formats.Sam, "-h")
    program.add_arg(n, bool, "-n")
    program.add_arg(R, str, "-R")
    program.add_arg(r, bool, "-r")
    program.add_arg(u, bool, "-u")
    return formats.Bam(out, program)


def index(in_, log=None):
    program = pipeline.add_node("samtools index", log)
    program.add_arg(in_, formats.Bam)
    return formats.Bai(in_.path + ".bai", program)


def faidx(in_, regions=None, log=None):
    program = pipeline.add_node("samtools faidx", log)
    program.add_arg(in_, formats.Fasta)
    program.add_args(regions, str)
    return formats.Fai(in_.path + ".fai", program)


def rmdup(in_, out, s=None, S=None, log=None):
    program = pipeline.add_node("samtools rmdup", log)
    program.add_arg(s, bool, "-s")
    program.add_arg(S, bool, "-S")
    program.add_arg(in_, formats.Bam)
    program.add_arg(out, str)
    return formats.Bam(out, program)


def calmd(in_, ref, out, A=None, e=None, u=None, b=None, S=None, C=None,
          r=None, E=None, log=None):
    program = pipeline.add_node("samtools calmd", log, out)
    program.add_arg(A, bool, "-A")
    program.add_arg(e, bool, "-e")
    program.add_arg(u, bool, "-u")
    program.add_arg(b, bool, "-b")
    program.add_arg(S, bool, "-S")
    program.add_arg(C, int, "-C")
    program.add_arg(r, bool, "-r")
    program.add_arg(E, bool, "-E")
    if S:
        program.add_arg(in_, formats.Sam)
    else:
        program.add_arg(in_, formats.Bam)
    program.add_arg(ref, formats.Fasta)
    if u or b:
        return formats.Bam(out, program)
    else:
        return formats.Sam(out, program)

