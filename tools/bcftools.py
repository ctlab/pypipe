import sys

import formats
from pipeline import create_program


# TODO Unknown files
def view(_in, out, A=None, b=None, D=None, F=None, G=None, l=None, N=None,
         Q=None, s=None, S=None, u=None, c=None, d=None, e=None, g=None,
         i=None, p=None, P=None, t=None, T=None, v=None, _1=None, U=None,
         X=None):
    if S:
        if type(_in) != formats.Vcf:
            sys.exit("bcftools view: '_in' must be a VCF file")
    else:
        if type(_in) != formats.Bcf:
            sys.exit("bcftools view: '_in' must be a BCF file")
    if type(out) != str:
        sys.exit("bcftools view: 'out' must be a string")
    if D and type(D) != formats.Unknown:
        sys.exit("bcftools view: 'D' must be an Unknown file")
    if l and type(l) != formats.Unknown:
        sys.exit("bcftools view: 'l' must be an Unknown file")
    if s and type(s) != formats.Unknown:
        sys.exit("bcftools view: 's' must be an Unknown file")
    if d and type(d) != float and type(d) != int:
        sys.exit("bcftools view: 'd' must be a float")
    if i and type(i) != float and type(i) != int:
        sys.exit("bcftools view: 'i' must be a float")
    if p and type(p) != float and type(p) != int:
        sys.exit("bcftools view: 'p' must be a float")
    if P and type(P) != str:
        sys.exit("bcftools view: 'P' must be a string")
    if t and type(t) != float and type(t) != int:
        sys.exit("bcftools view: 't' must be a float")
    if T and type(T) != str:
        sys.exit("bcftools view: 'T' must be a string")
    if _1 and type(_1) != int:
        sys.exit("bcftools view: '_1' must be an int")
    if U and type(U) != int:
        sys.exit("bcftools view: 'U' must be an int")
    if X and type(X) != float and type(X) != int:
        sys.exit("bcftools view: 'X' must be a float")
    cmd = ["bcftools", "view", _in.path]
    files_array = [_in]
    if A:
        cmd.append("-A")
    if b:
        cmd.append("-b")
    if D:
        cmd.append("-D")
        cmd.append(D.path)
        files_array.append(D)
    if F:
        cmd.append("-F")
    if G:
        cmd.append("-G")
    if l:
        cmd.append("-l")
        cmd.append(l.path)
        files_array.append(l)
    if N:
        cmd.append("-N")
    if Q:
        cmd.append("-Q")
    if s:
        cmd.append("-s")
        cmd.append(s.path)
        files_array.append(s)
    if S:
        cmd.append("-S")
    if u:
        cmd.append("-u")
    if c:
        cmd.append("-c")
    if d:
        cmd.append("-d")
        cmd.append(str(d))
    if e:
        cmd.append("-e")
    if g:
        cmd.append("-g")
    if i:
        cmd.append("-i")
        cmd.append(str(i))
    if p:
        cmd.append("-p")
        cmd.append(str(p))
    if P:
        cmd.append("-P")
        cmd.append(P)
    if t:
        cmd.append("-t")
        cmd.append(str(t))
    if T:
        cmd.append("-T")
        cmd.append(T)
    if v:
        cmd.append("-v")
    if _1:
        cmd.append("-1")
        cmd.append(str(_1))
    if U:
        cmd.append("-U")
        cmd.append(str(U))
    if X:
        cmd.append("-X")
        cmd.append(str(X))
    print " ".join(cmd)
    program = create_program(cmd, files_array, out)
    if b:
        return formats.Bcf(out, program)
    else:
        return formats.Vcf(out, program)


def cat(_in, out):
    if type(_in) != list and type(_in) != tuple:
        sys.exit("bcftools cat: '_in' must be a list of BCF files")
    for elem in _in:
        if type(elem) != formats.Bcf:
            sys.exit("bcftools cat: '_in' must be a list of BCF files")
    if len(_in) < 1:
        sys.exit("bcftools cat: '_in' length must be >= 1")
    if type(out) != str:
        sys.exit("bcftools cat: 'out' must be a string")
    files_array = list(_in)
    cmd = ["bcftools", "cat"] + [elem.path for elem in files_array]
    program = create_program(cmd, files_array, out)
    return formats.Bcf(out, program)

