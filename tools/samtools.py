import sys

import formats
from pipeline import create_program


# TODO Unknown format
def view(aln, o, b=None, f=None, F=None, h=None, H=None, l=None, q=None,
         r=None, R=None, S=None, c=None, t=None, u=None):
    program = create_program("samtools view")
    if S:
        program.add_argument(aln, formats.Sam)
    else:
        program.add_argument(aln, formats.Bam)
    program.add_argument(o, str, "-o")
    program.add_argument(b, bool, "-b")
    program.add_argument(f, int, "-f")
    program.add_argument(F, int, "-F")
    program.add_argument(h, bool, "-h")
    program.add_argument(H, bool, "-H")
    program.add_argument(l, str, "-l")
    program.add_argument(q, int, "-q")
    program.add_argument(r, str, "-r")
    program.add_argument(R, formats.Unknown, "-R")
    program.add_argument(S, bool, "-S")
    program.add_argument(c, bool, "-c")
    program.add_argument(t, formats.Unknown, "-t")
    program.add_argument(u, bool, "-u")
    print " ".join(program.cmd)  # debug
    if b:
        return formats.Bam(o, program)
    else:
        return formats.Sam(o, program)
    

def sort(aln, out, n=None, m=None):
    program = create_program("samtools sort")
    program.add_argument(aln, formats.Bam)
    program.add_argument(out, str)
    program.add_argument(n, bool, "-n")
    program.add_argument(m, int, "-m")
    print " ".join(program.cmd)  # debug
    return formats.Bam(out + ".bam", program)


def mpileup(_in, out, _6=None, A=None, B=None, b=None, C=None, d=None,
            E=None, f=None, l=None, q=None, Q=None, r=None, D=None, g=None,
            S=None, u=None, e=None, h=None, I=None, L=None, o=None, P=None):
    program = create_program("samtools mpileup", out)
    program.add_arguments(_in, formats.Bam)
    program.add_argument(_6, bool, "-6")
    program.add_argument(A, bool, "-A")
    program.add_argument(B, bool, "-B")
    program.add_argument(b, formats.Unknown, "-b")
    program.add_argument(C, int, "-C")
    program.add_argument(d, int, "-d")
    program.add_argument(E, bool, "-E")
    program.add_argument(f, formats.Fasta, "-f")
    program.add_argument(l, formats.Bed, "-l")
    program.add_argument(q, int, "-q")
    program.add_argument(Q, int, "-Q")
    program.add_argument(r, str, "-r")
    program.add_argument(D, bool, "-D")
    program.add_argument(g, bool, "-g")
    program.add_argument(S, bool, "-S")
    program.add_argument(u, bool, "-u")
    program.add_argument(e, int, "-e")
    program.add_argument(h, int, "-h")
    program.add_argument(I, bool, "-I")
    program.add_argument(L, int, "-L")
    program.add_argument(o, int, "-o")
    program.add_argument(P, str, "-P")
    print " ".join(program.cmd)  # debug
    return formats.Bcf(out, program)


def merge(_in, out, _1=None, f=None, h=None, n=None, R=None, r=None, u=None):
    program = create_program("samtools merge")
    program.add_argument(out, str)
    program.add_arguments(_in, formats.Bam, 2)
    program.add_argument(_1, bool, "-1")
    program.add_argument(f, bool, "-f")
    program.add_argument(h, formats.Sam, "-h")
    program.add_argument(n, bool, "-n")
    program.add_argument(R, str, "-R")
    program.add_argument(r, bool, "-r")
    program.add_argument(u, bool, "-u")
    print " ".join(program.cmd)  # debug
    return formats.Bam(out, program)

