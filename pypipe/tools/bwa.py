import sys

from pypipe import formats
from pypipe.utils import create_program


def mem(_ref, _in1, _out, _in2=None, t=None, k=None, w=None, d=None, r=None,
        c=None, P=None, B=None, O=None, E=None, L=None, U=None, p=None,
        R=None, T=None, a=None, C=None, H=None, M=None, v=None, log=None):
    program = create_program("bwa mem", log, _out)
    program.add_arg(t, int, "-t")
    program.add_arg(k, int, "-k")
    program.add_arg(w, int, "-w")
    program.add_arg(d, int, "-d")
    program.add_arg(r, float, "-r")
    program.add_arg(c, int, "-c")
    program.add_arg(P, bool, "-P")
    program.add_arg(a, int, "-a")
    program.add_arg(B, int, "-B")
    program.add_arg(O, int, "-O")
    program.add_arg(E, int, "-E")
    program.add_arg(L, int, "-L")
    program.add_arg(U, int, "-U")
    program.add_arg(p, bool, "-p")
    program.add_arg(R, str, "-R")
    program.add_arg(T, int, "-T")
    program.add_arg(a, bool, "-a")
    program.add_arg(C, bool, "-C")
    program.add_arg(H, bool, "-H")
    program.add_arg(M, bool, "-M")
    program.add_arg(v, int, "-v")
    program.add_arg(_ref, formats.BwaIndex)
    program.add_arg(_in1, formats.Fastq)
    program.add_arg(_in2, formats.Fastq)
    return formats.Sam(_out, program)


def aln(_ref, _in, _out, n=None, o=None, e=None, d=None, i=None, l=None, k=None,
        t=None, M=None, O=None, E=None, R=None, c=None, N=None,
        q=None, I=None, B=None, b=None, _0=None, _1=None, _2=None, log=None):
    program = create_program("bwa aln", log, _out)
    program.add_arg(_ref, formats.BwaIndex)
    program.add_arg(_in, formats.Fastq)
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
    return formats.Sai(_out, program)


def samse(_ref, _sai, _in, _out, n=None, r=None, log=None):
    program = create_program("bwa samse", log, _out)
    program.add_arg(_ref, formats.BwaIndex)
    program.add_arg(_sai, formats.Sai)
    program.add_arg(_in, formats.Fastq)
    program.add_arg(n, int, "-n")
    program.add_arg(r, str, "-r")
    return formats.Sam(_out, program)


def sampe(_ref, _sai1, _sai2, _in1, _in2, _out, a=None, o=None,
          P=None, n=None, N=None, r=None, log=None):
    program = create_program("bwa sampe", log, _out)
    program.add_arg(_ref, formats.BwaIndex)
    program.add_arg(_sai1, formats.Sai)
    program.add_arg(_sai2, formats.Sai)
    program.add_arg(_in1, formats.Fastq)
    program.add_arg(_in2, formats.Fastq)
    program.add_arg(a, int, '-a')
    program.add_arg(o, int, '-o')
    program.add_arg(P, bool, '-P')
    program.add_arg(n, int, '-n')
    program.add_arg(N, int, '-N')
    program.add_arg(r, str, '-r')
    return formats.Sam(_out, program)


def bwasw(_ref, _in1, _out, _in2=None, a=None, b=None, q=None, r=None, t=None,
          w=None, T=None, c=None, z=None, s=None, N=None, log=None):
    program = create_program("bwa bwasw", log, _out)
    program.add_arg(_ref, formats.BwaIndex)
    program.add_arg(_in1, formats.Fastq)
    program.add_arg(_in2, formats.Fastq)
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
    return formats.Sam(_out, program)

