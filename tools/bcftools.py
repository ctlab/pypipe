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
        program.add_argument(_in, formats.Vcf)
    else:
        program.add_argument(_in, formats.Bcf)
    program.add_argument(A, bool, "-A")
    program.add_argument(b, bool, "-b")
    program.add_argument(D, formats.Unknown, "-D")
    program.add_argument(F, bool, "-F")
    program.add_argument(G, bool, "-G")
    program.add_argument(l, formats.Unknown, "-l")
    program.add_argument(N, bool, "-N")
    program.add_argument(Q, bool, "-Q")
    program.add_argument(s, formats.Unknown, "-s")
    program.add_argument(S, bool, "-S")
    program.add_argument(u, bool, "-u")
    program.add_argument(c, bool, "-c")
    program.add_argument(d, float, "-d")
    program.add_argument(e, bool, "-e")
    program.add_argument(g, bool, "-g")
    program.add_argument(i, float, "-i")
    program.add_argument(p, float, "-p")
    program.add_argument(P, str, "-P")
    program.add_argument(t, float, "-t")
    program.add_argument(T, str, "-T")
    program.add_argument(v, bool, "-v")
    program.add_argument(_1, int, "-1")
    program.add_argument(U, int, "-U")
    program.add_argument(X, float, "-X")
    print " ".join(program.cmd)  # debug
    if b:
        return formats.Bcf(out, program)
    else:
        return formats.Vcf(out, program)


def cat(_in, out):
    program = create_program("bcftools cat", out)
    program.add_arguments(_in, formats.Bcf)
    print " ".join(program.cmd)  # debug
    return formats.Bcf(out, program)

