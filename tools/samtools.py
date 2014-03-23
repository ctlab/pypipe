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
    return formats.Bam(out, program)
