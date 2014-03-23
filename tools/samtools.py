import sys

import formats
from pipeline import create_program


# TODO Unknown format
def view(align, o, b=None, f=None, F=None, h=None, H=None, l=None, q=None,
         r=None, R=None, S=None, c=None, t=None, u=None):
    if S and type(align) != formats.Sam:
        sys.exit("samtools view: 'align' must be SAM file")
    if not S and type(align) != formats.Bam:
        sys.exit("samtools view: 'align' must be BAM file")
    if type(o) != str:
        sys.exit("samtools view: 'o' must be string")
    if f and type(f) != int:
        sys.exit("samtools view: 'f' must be int")
    if F and type(F) != int:
        sys.exit("samtools view: 'F' must be int")
    if l and type(l) != str:
        sys.exit("samtools view: 'l' must be string")
    if q and type(q) != int:
        sys.exit("samtools view: 'q' must be int")
    if r and type(r) != str:
        sys.exit("samtools view: 'r' must be string")
    if R and type(R) != formats.Unknown:
        sys.exit("samtools view: 'R' must be Unknown file")
    if t and type(t) != formats.Unknown:
        sys.exit("samtools view: 't' must be Unknown file")
    cmd = ["samtools", "view", align.path, "-o", o]
    files_array = [align]
    if b:
        cmd.append("-b")
    if f:
        cmd.append("-f")
        cmd.append(str(f))
    if F:
        cmd.append("-F")
        cmd.append(str(F))
    if h:
        cmd.append("-h")
    if H:
        cmd.append("-H")
    if l:
        cmd.append("-l")
        cmd.append(l)
    if q:
        cmd.append("-q")
        cmd.append(str(q))
    if r:
        cmd.append("-r")
        cmd.append(r)
    if R:
        cmd.append("-R")
        cmd.append(R.path)
        files_array.append(R)
    if S:
        cmd.append("-S")
    if c:
        cmd.append("-c")
    if t:
        cmd.append("-t")
        cmd.append(t.path)
        files_array.append(t)
    if u:
        cmd.append("-u")
    print " ".join(cmd)
    program = create_program(cmd, files_array)
    if b:
        return formats.Bam(o, program)
    else:
        return formats.Sam(o, program)
    

def sort(align, out, n=None, m=None):
    if type(align) != formats.Bam:
        sys.exit("samtools sort: 'align' must be BAM file")
    if type(out) != str:
        sys.exit("samtools sort: 'out' must be string")
    if m and type(m) != int:
        sys.exit("samtools sort: 'm' must be int")
    cmd = ["samtools", "sort", align.path, out]
    if n:
        cmd.append("-n")
    if m:
        cmd.append("-m")
        cmd.append(str(m))
    print " ".join(cmd)
    program = create_program(cmd, [align])
    return formats.Bam(out + ".bam", program)


def mpileup(align, out, _6=None, A=None, B=None, b=None, C=None, d=None,
            E=None, f=None, l=None, q=None, Q=None, r=None, D=None, g=None,
            S=None, u=None, e=None, h=None, I=None, L=None, o=None, P=None):
    if type(align) != formats.Bam:
        sys.exit("samtools mpileup: 'align' must be BAM file")
    if type(out) != str:
        sys.exit("samtools mpileup: 'out' must be string")
    if b and type(b) != formats.Unknown:
        sys.exit("samtools mpileup: 'b' must be Unknown file")
    if C and type(C) != int:
        sys.exit("samtools mpileup: 'C' must be int")
    if d and type(d) != int:
        sys.exit("samtools mpileup: 'd' must be int")
    if f and type(f) != formats.Fasta:
        sys.exit("samtools mpileup: 'f' must be FASTA file")
    if l and type(l) != formats.Bed:
        sys.exit("samtools mpileup: 'l' must be BED file")
    if q and type(q) != int:
        sys.exit("samtools mpileup: 'q' must be int")
    if Q and type(Q) != int:
        sys.exit("samtools mpileup: 'Q' must be int")
    if r and type(r) != str:
        sys.exit("samtools mpileup: 'w' must be string")
    if e and type(e) != int:
        sys.exit("samtools mpileup: 'e' must be int")
    if h and type(h) != int:
        sys.exit("samtools mpileup: 'h' must be int")
    if L and type(L) != int:
        sys.exit("samtools mpileup: 'L' must be int")
    if o and type(o) != int:
        sys.exit("samtools mpileup: 'o' must be int")
    if P and type(P) != str:
        sys.exit("samtools mpileup: 'P' must be string")
    cmd = ["samtools", "mpileup", align.path]
    files_array = [align]
    if _6:
        cmd.append("-6")
    if B:
        cmd.append("-B")
    if b:
        cmd.append("-b")
        cmd.append(b.path)
        files_array.append(b)
    if C:
        cmd.append("-C")
        cmd.append(str(C))
    if d:
        cmd.append("-d")
        cmd.append(str(d))
    if E:
        cmd.append("-E")
    if f:
        cmd.append("-f")
        cmd.append(f.path)
        files_array.append(f)
    if l:
        cmd.append("-l")
        cmd.append(l.path)
        files_array.append(l)
    if q:
        cmd.append("-q")
        cmd.append(str(q))
    if Q:
        cmd.append("-Q")
        cmd.append(str(Q))
    if r:
        cmd.append("-r")
        cmd.append(Q)
    if D:
        cmd.append("-D")
    if g:
        cmd.append("-g")
    if S:
        cmd.append("-S")
    if u:
        cmd.append("-u")
    if e:
        cmd.append("-e")
        cmd.append(str(e))
    if h:
        cmd.append("-h")
        cmd.append(str(h))
    if I:
        cmd.append("-I")
    if L:
        cmd.append("-L")
        cmd.append(str(L))
    if o:
        cmd.append("-o")
        cmd.append(str(o))
    if P:
        cmd.append("-P")
        cmd.append(P)
    print " ".join(cmd)
    program = create_program(cmd, files_array, out)
    return formats.Bcf(out, program)











