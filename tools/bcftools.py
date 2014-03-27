import sys

import formats
from utils import create_program 


def view(_in, _out, A=None, b=None, D=None, F=None, G=None, l=None, N=None,
         Q=None, s=None, S=None, u=None, c=None, d=None, e=None, g=None,
         i=None, p=None, P=None, t=None, T=None, v=None, _1=None, U=None,
         X=None):
    program = create_program("bcftools view", _out)
    if S:
        program.add_arg(_in, formats.Vcf)
    else:
        program.add_arg(_in, formats.Bcf)
    program.add_arg(A, bool, "-A")
    program.add_arg(b, bool, "-b")
    program.add_arg(D, formats.TextFile, "-D")
    program.add_arg(F, bool, "-F")
    program.add_arg(G, bool, "-G")
    program.add_arg(l, formats.TextFile, "-l")
    program.add_arg(N, bool, "-N")
    program.add_arg(Q, bool, "-Q")
    program.add_arg(s, formats.TextFile, "-s")
    program.add_arg(S, bool, "-S")
    program.add_arg(u, bool, "-u")
    program.add_arg(c, bool, "-c")
    program.add_arg(d, float, "-d")
    program.add_arg(e, bool, "-e")
    program.add_arg(g, bool, "-g")
    program.add_arg(i, float, "-i")
    program.add_arg(p, float, "-p")
    program.add_arg(P, str, "-P")
    program.add_arg(t, float, "-t")
    program.add_arg(T, str, "-T")
    program.add_arg(v, bool, "-v")
    program.add_arg(_1, int, "-1")
    program.add_arg(U, int, "-U")
    program.add_arg(X, float, "-X")
    if b:
        return formats.Bcf(_out, program)
    else:
        return formats.Vcf(_out, program)


def cat(_in, _out):
    program = create_program("bcftools cat", _out)
    program.add_args(_in, formats.Bcf)
    return formats.Bcf(_out, program)

