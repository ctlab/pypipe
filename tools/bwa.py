import sys

import formats
from utils import create_program


def aln(ref, read, out, n=None, o=None, e=None, d=None, i=None, l=None, k=None,
        t=None, M=None, O=None, E=None, R=None, c=None, N=None,
        q=None, I=None, B=None, b=None, _0=None, _1=None, _2=None):
    program = create_program("bwa aln", out)
    program.add_arg(ref, formats.Fasta)
    program.add_arg(read, formats.Fastq)
    program.add_arg(n, int, '-n')
    program.add_arg(o, int, '-o')
    program.add_arg(e, int, '-e')
    program.add_arg(d, int, '-d')
    program.add_arg(i, int, '-i')
    program.add_arg(l, int, '-l')
    program.add_arg(k, int, '-k')
    program.add_arg(t, int, '-t')
    program.add_arg(M, int, '-M')
    program.add_arg(O, int, '-O')
    program.add_arg(E, int, '-E')
    program.add_arg(R, int, '-R')
    program.add_arg(c, bool, '-c')
    program.add_arg(N, bool, '-N')
    program.add_arg(q, int, '-q')
    program.add_arg(I, bool, '-I')
    program.add_arg(B, int, '-B')
    program.add_arg(b, bool, '-b')
    program.add_arg(_0, bool, '-0')
    program.add_arg(_1, bool, '-1')
    program.add_arg(_2, bool, '-2')
    return formats.Sai(out, program)


def samse(ref, sai, read, out, n=None, r=None):
    program = create_program("bwa samse", out)
    program.add_arg(ref, formats.Fasta)
    program.add_arg(sai, formats.Sai)
    program.add_arg(read, formats.Fastq)
    program.add_arg(n, int, "-n")
    program.add_arg(r, str, "-r")
    return formats.Sam(out, program)


def sampe(ref, sai1, sai2, in1, in2, out, a=None, o=None,
          P=None, n=None, N=None, r=None):
    program = create_program("bwa sampe", out)
    program.add_arg(ref, formats.Fasta)
    program.add_arg(sai1, formats.Sai)
    program.add_arg(sai2, formats.Sai)
    program.add_arg(in1, formats.Fastq)
    program.add_arg(in2, formats.Fastq)
    program.add_arg(a, int, '-a')
    program.add_arg(o, int, '-o')
    program.add_arg(P, bool, '-P')
    program.add_arg(n, int, '-n')
    program.add_arg(N, int, '-N')
    program.add_arg(r, str, '-r')
    return formats.Sam(out, program)


def bwasw(ref, in1, out, in2=None, a=None, b=None, q=None, r=None, t=None,
          w=None, T=None, c=None, z=None, s=None, N=None):
    program = create_program("bwa bwasw", out)
    program.add_arg(ref, formats.Fasta)
    program.add_arg(in1, formats.Fastq)
    program.add_arg(in2, formats.Fastq)
    program.add_arg(a, int, '-a')
    program.add_arg(b, int, '-b')
    program.add_arg(q, int, '-q')
    program.add_arg(r, int, '-r')
    program.add_arg(t, int, '-t')
    program.add_arg(w, int, '-w')
    program.add_arg(T, int, '-T')
    program.add_arg(c, float, '-c')
    program.add_arg(z, int, '-z')
    program.add_arg(s, int, '-s')
    program.add_arg(N, int, '-N')
    return formats.Sam(out, program)

