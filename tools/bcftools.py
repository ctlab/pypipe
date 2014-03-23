import sys

import formats
from pipeline import create_program


# TODO Unknown files
def view(variant, out, A=None, b=None, D=None, F=None, G=None, l=None, N=None,
         Q=None, s=None, S=None, u=None, c=None, d=None, e=None, g=None,
         i=None, p=None, P=None, t=None, T=None, v=None, _1=None, U=None,
         X=None):
    if S:
        if type(variant) != formats.Vcf:
            sys.exit("bcftools view: 'variant' must be VCF file")
    else:
        if type(variant) != formats.Bcf:
            sys.exit("bcftools view: 'variant' must be BCF file")
    if type(out) != str:
        sys.exit("bcftools view: 'out' must be string")
    if D and type(D) != formats.Unknown:
        sys.exit("bcftools view: 'D' must be Unknown file")
    if l and type(l) != formats.Unknown:
        sys.exit("bcftools view: 'l' must be Unknown file")
    if s and type(s) != formats.Unknown:
        sys.exit("bcftools view: 's' must be Unknown file")
    if d and type(d) != float and type(d) != int:
        sys.exit("bcftools view: 'd' must be float")
    if i and type(i) != float and type(i) != int:
        sys.exit("bcftools view: 'i' must be float")
    if p and type(p) != float and type(p) != int:
        sys.exit("bcftools view: 'p' must be float")
    if P and type(P) != str:
        sys.exit("bcftools view: 'P' must be string")
    if t and type(t) != float and type(t) != int:
        sys.exit("bcftools view: 't' must be float")
    if T and type(T) != str:
        sys.exit("bcftools view: 'T' must be string")
    if _1 and type(_1) != int:
        sys.exit("bcftools view: '_1' must be int")
    if U and type(U) != int:
        sys.exit("bcftools view: 'U' must be int")
    if X and type(X) != float and type(X) != int:
        sys.exit("bcftools view: 'X' must be float")
    cmd = ["bcftools", "view", variant.path]
    files_array = [variant]
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

