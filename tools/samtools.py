import sys

import formats
from pipeline import create_program


# TODO Unknown format
def view(aln, o, b=None, f=None, F=None, h=None, H=None, l=None, q=None,
         r=None, R=None, S=None, c=None, t=None, u=None):
    if S and type(aln) != formats.Sam:
        sys.exit("samtools view: 'aln' must be a SAM file")
    if not S and type(aln) != formats.Bam:
        sys.exit("samtools view: 'aln' must be a BAM file")
    if type(o) != str:
        sys.exit("samtools view: 'o' must be a string")
    if f and type(f) != int:
        sys.exit("samtools view: 'f' must be an int")
    if F and type(F) != int:
        sys.exit("samtools view: 'F' must be an int")
    if l and type(l) != str:
        sys.exit("samtools view: 'l' must be a string")
    if q and type(q) != int:
        sys.exit("samtools view: 'q' must be an int")
    if r and type(r) != str:
        sys.exit("samtools view: 'r' must be a string")
    if R and type(R) != formats.Unknown:
        sys.exit("samtools view: 'R' must be an Unknown file")
    if t and type(t) != formats.Unknown:
        sys.exit("samtools view: 't' must be an Unknown file")
    cmd = ["samtools", "view", aln.path, "-o", o]
    files_array = [aln]
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
    

def sort(aln, out, n=None, m=None):
    if type(aln) != formats.Bam:
        sys.exit("samtools sort: 'aln' must be a BAM file")
    if type(out) != str:
        sys.exit("samtools sort: 'out' must be a string")
    if m and type(m) != int:
        sys.exit("samtools sort: 'm' must be an int")
    cmd = ["samtools", "sort", aln.path, out]
    if n:
        cmd.append("-n")
    if m:
        cmd.append("-m")
        cmd.append(str(m))
    print " ".join(cmd)
    program = create_program(cmd, [aln])
    return formats.Bam(out + ".bam", program)


def mpileup(aln, out, _6=None, A=None, B=None, b=None, C=None, d=None,
            E=None, f=None, l=None, q=None, Q=None, r=None, D=None, g=None,
            S=None, u=None, e=None, h=None, I=None, L=None, o=None, P=None):
    if type(aln) != tuple and type(aln) != list:
        sys.exit("samtools mpileup: 'aln' must be an array of BAM files")
    for elem in aln:
        if type(elem) != formats.Bam:
            sys.exit("samtools mpileup: 'aln' must be an array of BAM files")
    if len(aln) < 1:
        sys.exit("samtools mpileup: 'aln' length must be >= 1")
    if type(out) != str:
        sys.exit("samtools mpileup: 'out' must be a string")
    if b and type(b) != formats.Unknown:
        sys.exit("samtools mpileup: 'b' must be an Unknown file")
    if C and type(C) != int:
        sys.exit("samtools mpileup: 'C' must be an int")
    if d and type(d) != int:
        sys.exit("samtools mpileup: 'd' must be an int")
    if f and type(f) != formats.Fasta:
        sys.exit("samtools mpileup: 'f' must be a FASTA file")
    if l and type(l) != formats.Bed:
        sys.exit("samtools mpileup: 'l' must be a BED file")
    if q and type(q) != int:
        sys.exit("samtools mpileup: 'q' must be an int")
    if Q and type(Q) != int:
        sys.exit("samtools mpileup: 'Q' must be an int")
    if r and type(r) != str:
        sys.exit("samtools mpileup: 'w' must be a string")
    if e and type(e) != int:
        sys.exit("samtools mpileup: 'e' must be an int")
    if h and type(h) != int:
        sys.exit("samtools mpileup: 'h' must be an int")
    if L and type(L) != int:
        sys.exit("samtools mpileup: 'L' must be an int")
    if o and type(o) != int:
        sys.exit("samtools mpileup: 'o' must be an int")
    if P and type(P) != str:
        sys.exit("samtools mpileup: 'P' must be a string")
    cmd = ["samtools", "mpileup"] + [elem.path for elem in aln]
    files_array = list(aln)
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


def merge(_in, out, _1=None, f=None, h=None, n=None, R=None, r=None, u=None):
    if type(_in) != tuple and type(_in) != list:
        sys.exit("samtools merge: '_in' must be an array of BAM files")
    for elem in _in:
        if type(elem) != formats.Bam:
            sys.exit("samtools merge: '_in' must be an array of BAM files")
    if len(_in) <= 1:
        sys.exit("samtools merge: '_in' length must be >= 2")
    if type(out) != str:
        sys.exit("samtools merge: 'out' must be a string")
    if h and type(h) != formats.Sam:
        sys.exit("samtools merge: 'h' must be a SAM file")
    if R and type(R) != str:
        sys.exit("samtools merge: 'R' must be a string")
    files_array = list(_in)
    cmd = ["samtools", "merge", out] + [elem.path for elem in files_array]
    if _1:
        cmd.append("-1")
    if f:
        cmd.append("-f")
    if h:
        cmd.append("-h")
        cmd.append(h.path)
        files_array.append(h)
    if n:
        cmd.append("-n")
    if R:
        cmd.append("-R")
        cmd.append(R)
    if r:
        cmd.append("-r")
    if u:
        cmd.append("-u")
    program = create_program(cmd, files_array)
    return formats.Bam(out, program)

