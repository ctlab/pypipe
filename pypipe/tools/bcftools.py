import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program("bcftools.sh", "bcftools")


# TODO index
def view(in_, out, A=None, b=None, D=None, F=None, G=None, l=None, N=None,
         Q=None, s=None, S=None, u=None, c=None, d=None, e=None, g=None,
         i=None, p=None, P=None, t=None, T=None, v=None, _1=None, U=None,
         X=None, log=None):
    program = create_program("bcftools view", log, out)
    if S:
        program.add_arg(in_, formats.Vcf)
    else:
        program.add_arg(in_, formats.Bcf)
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
        return formats.Bcf(out, program)
    else:
        return formats.Vcf(out, program)


def cat(in_, out, log=None):
    program = create_program("bcftools cat", log, out)
    program.add_args(in_, formats.Bcf)
    return formats.Bcf(out, program)

