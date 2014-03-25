import sys

import formats
from pipeline import create_program 


# TODO Unknown files
def view(_in, out, A=None, b=None, D=None, F=None, G=None, l=None, N=None,
         Q=None, s=None, S=None, u=None, c=None, d=None, e=None, g=None,
         i=None, p=None, P=None, t=None, T=None, v=None, _1=None, U=None,
         X=None):
    program = create_program("bcftools view", out)
    if S:
        program.add_arg(_in, formats.Vcf)
    else:
        program.add_arg(_in, formats.Bcf)
    program.add_arg(A, bool, "-A")
    program.add_arg(b, bool, "-b")
    program.add_arg(D, formats.Unknown, "-D")
    program.add_arg(F, bool, "-F")
    program.add_arg(G, bool, "-G")
    program.add_arg(l, formats.Unknown, "-l")
    program.add_arg(N, bool, "-N")
    program.add_arg(Q, bool, "-Q")
    program.add_arg(s, formats.Unknown, "-s")
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
    print " ".join(program.cmd)  # debug
    if b:
        return formats.Bcf(out, program)
    else:
        return formats.Vcf(out, program)


def cat(_in, out):
    program = create_program("bcftools cat", out)
    program.add_args(_in, formats.Bcf)
    print " ".join(program.cmd)  # debug
    return formats.Bcf(out, program)

