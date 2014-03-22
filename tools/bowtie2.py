import sys

import formats
from pipeline import create_program


# TODO r?
# TODO -o / --offrate 
# TODO output options
def bowtie2(x, S, U=None, _1=None, _2=None, q=None, qseq=None, f=None,
            r=None, c=None, s=None, skip=None, u=None, qupto=None,
            trim5=None, _5=None, trim3=None, _3=None, phred33=None,
            phred64=None, solexa_quals=None, int_quals=None, very_fast=None,
            fast=None, sensitive=None, very_sensitive=None,
            very_fast_local=None, fast_local=None, sensitive_local=None,
            very_sensitive_local=None, N=None, L=None, i=None, n_ceil=None,
            dpad=None, gbar=None, ignore_quals=None, nofw=None, norc=None,
            no_1mm_upfront=None, end_to_end=None, local=None, ma=None,
            mp=None, np=None, rdg=None, rfg=None, score_min=None, k=None,
            a=None, D=None, R=None, I=None, minins=None, X=None, maxins=None,
            fr=None, rf=None, ff=None, no_mixed=None, no_discordant=None,
            dovetail=None, no_contain=None, no_overlap=None, t=None,
            time=None, un=None, un_gz=None, un_bz2=None, al=None, al_gz=None,
            al_bz2=None, un_conc=None, un_conc_gz=None, un_conc_bz2=None,
            al_conc=None, al_conc_gz=None, al_conc_bz2=None, quiet=None,
            met_file=None, met_stderr=None, met=None, no_unal=None,
            no_hd=None, no_sq=None, rg_id=None, rg=None, omit_sec_seq=None,
            o=None, offrate=None, p=None, threads=None, reorder=None, mm=None,
            qc_filter=None, seed=None, non_deterministic=None):
    # Checking conflicts
    if (U and _1) or (U and _2) or (not U and not (_1 and _2)):
        sys.exit("bowtie2: Wrong options. Use only 'U' or '_1' and '_2'")
    if q and (qseq or f or r or c) or qseq and (
            f or r or c) or f and (r or c) or r and c:
        sys.exit("bowtie2: 'q', 'qseq', 'f', 'r', 'c' - this options conflict")
    if s and skip:
        sys.exit("bowtie2: use only 'skip' or 's'")
    if u and qupto:
        sys.exit("bowtie2: use only 'u' or 'qupto'")
    if _5 and trim5:
        sys.exit("bowtie2: use only '_5' or 'trim5'")
    if _3 and trim3:
        sys.exit("bowtie2: use only '_3' or 'trim3'")
    if I and minins:
        sys.exit("bowtie2: use only 'I' or 'minins'")
    if X and maxins:
        sys.exit("bowtie2: use only 'X' or 'maxins'")
    if o and offrate:
        sys.exit("bowtie2: use only 'o' or 'offrate'")
    if p and threads:
        sys.exit("bowtie2: use only 'p' or 'threads'")


    # Checking types
    if type(x) != formats.Bowtie2Index:
        sys.exit("bowtie2: 'x' argument must be Bowtie2Index")
    if type(S) != str:
        sys.exit("bowtie2: 'S' argument must be string")
    if s and type(s) != int:
        sys.exit("bowtie2: 's' argument must be int")
    if skip and type(skip) != int:
        sys.exit("bowtie2: 'skip' argument must be int")
    if u and type(u) != int:
        sys.exit("bowtie2: 'u' argument must be int")
    if qupto and type(qupto) != int:
        sys.exit("bowtie2: 'qupto' argument must be int")
    if trim5 and type(trim5) != int:
        sys.exit("bowtie2: 'trim5' argument must be int")
    if _5 and type(_5) != int:
        sys.exit("bowtie2: '_5' argument must be int")
    if trim3 and type(trim3) != int:
        sys.exit("bowtie2: 'trim3' argument must be int")
    if _3 and type(_3) != int:
        sys.exit("bowtie2: '_3' argument must be int")
    if N and type(N) != int:
        sys.exit("bowtie2: 'N' argument must be int")
    if L and type(L) != int:
        sys.exit("bowtie2: 'L' argument must be int")
    if i and type(i) != str:
        sys.exit("bowtie2: 'i' argument must be string")
    if n_ceil and type(n_ceil) != str:
        sys.exit("bowtie2: 'n_ceil' argument must be string")
    if dpad and type(dpad) != int:
        sys.exit("bowtie2: 'dpad' argument must be int")
    if gbar and type(gbar) != int:
        sys.exit("bowtie2: 'gbar' argument must be int")
    if k and type(k) != int:
        sys.exit("bowtie2: 'k' argument must be int")
    if D and type(D) != int:
        sys.exit("bowtie2: 'D' argument must be int")
    if R and type(R) != int:
        sys.exit("bowtie2: 'R' argument must be int")
    if ma and type(ma) != int:
        sys.exit("bowtie2: 'ma' argument must be int")
    if mp and (type(mp) != tuple or type(mp[0]) != int or type(mp[1]) != int):
        sys.exit("bowtie2: 'mp' argument must be tuple of two int")
    if np and type(np) != int:
        sys.exit("bowtie2: 'np' argument must be int")
    if rdg and (type(rdg) != tuple or
            type(rdg[0]) != int or type(rdg[1]) != int):
        sys.exit("bowtie2: 'rdg' argument must be tuple of two int")
    if rfg and (type(rfg) != tuple or
            type(rfg[0]) != int or type(rfg[1]) != int):
        sys.exit("bowtie2: 'rfg' argument must be tuple of two int")
    if score_min and type(score_min) != str:
        sys.exit("bowtie2: 'score_min' argument must be string")
    if I and type(I) != int:
        sys.exit("bowtie2: 'I' argument must be int")
    if minins and type(minins) != int:
        sys.exit("bowtie2: 'minins' argument must be int")
    if X and type(X) != int:
        sys.exit("bowtie2: 'X' argument must be int")
    if maxins and type(maxins) != int:
        sys.exit("bowtie2: 'maxins' argument must be int")
    if rg_id and type(rg_id) != str:
        sys.exit("bowtie2: 'rg_id' argument must be string")
    if rg and type(rg) != str:
        sys.exit("bowtie2: 'rg' argument must be string")
    if o and type(o) != int:
        sys.exit("bowtie2: 'o' argument must be int")
    if offrate and type(offrate) != int:
        sys.exit("bowtie2: 'offrate' argument must be int")
    if p and type(p) != int:
        sys.exit("bowtie2: 'p' argument must be int")
    if threads and type(threads) != int:
        sys.exit("bowtie2: 'threads' argument must be int")
    if seed and type(seed) != int:
        sys.exit("bowtie2: 'seed' argument must be int")

    if qseq:
        read_format = formats.Qseq
        str_format = "QSEQ"
    elif f:
        read_format = formats.Fasta
        str_format = "FASTA"
    elif r:
        read_format = formats.Unknown
        str_format = "UNKNOWN"
    elif c:
        read_format = str
        str_format = "string"
    else:
        read_format = formats.Fastq
        str_format = "FASTQ"
    error_msg = "bowtie2: U argument must be array of " + str_format
    if U:
        for read in U:
            if type(read) != read_format:
                sys.exit(error_msg)
    else:
        for read in _1:
            if type(read) != read_format:
                sys.exit(error_msg)
        for read in _2:
            if type(read) != read_format:
                sys.exit(error_msg)

    # Generating basic command
    if U:
        if c:
            reads = ",".join(U)
        else:
            reads = ",".join([read.path for read in U])
        cmd = ["bowtie2", "-x", x.path, "-U", reads, "-S", S]
    else:
        if c:
            reads1 = ",".join(_1)
            reads2 = ",".join(_2)
        else:
            reads1 = ",".join([read.path for read in _1])
            reads2 = ",".join([read.path for read in _2])
        cmd = ["bowtie2", "-x", x.path, "-1", reads1, "-2", reads2, "-S", S]

    # Adding options
    if q:
        cmd.append("-q")
    elif qseq:
        cmd.append("--qseq")
    elif f:
        cmd.append("-f")
    elif r:
        cmd.append("-r")
    elif c:
        cmd.append("-c")
    if s:
        cmd.append("-s")
        cmd.append(str(s))
    elif skip:
        cmd.append("--skip")
        cmd.append(str(skip))
    if u:
        cmd.append("-u")
        cmd.append(str(u))
    elif qupto:
        cmd.append("--qupto")
        cmd.append(str(qupto))
    if trim5:
        cmd.append("--trim5")
        cmd.append(str(trim5))
    if _5:
        cmd.append("-5")
        cmd.append(str(_5))
    if trim3:
        cmd.append("--trim3")
        cmd.append(str(trim3))
    if _3:
        cmd.append("-3")
        cmd.append(str(_3))
    if phred33:
        cmd.append("--phred33")
    if phred64:
        cmd.append("--phred64")
    if solexa_quals:
        cmd.append("--solexa-quals")
    if int_quals:
        cmd.append("--int-quals")
    if very_fast:
        cmd.append("--very-fast")
    if fast:
        cmd.append("--fast")
    if sensitive:
        cmd.append("--sensitive")
    if very_sensitive:
        cmd.append("--very-sensitive")
    if very_fast_local:
        cmd.append("--very-fast-local")
    if fast_local:
        cmd.append("--fast-local")
    if sensitive_local:
        cmd.append("--sensitive-local")
    if very_sensitive_local:
        cmd.append("--very-sensitive-local")
    if N:
        cmd.append("-N")
        cmd.append(str(N))
    if L:
        cmd.append("-L")
        cmd.append(str(L))
    if i:
        cmd.append("-i")
        cmd.append(i)
    if n_ceil:
        cmd.append("--n-ceil")
        cmd.append(n_ceil)
    if dpad:
        cmd.append("--dpad")
        cmd.append(str(dpad))
    if gbar:
        cmd.append("--gbar")
        cmd.append(str(gbar))
    if ignore_quals:
        cmd.append("--ignore-quals")
    if nofw:
        cmd.append("--nofw")
    if norc:
        cmd.append("--norc")
    if no_1mm_upfront:
        cmd.append("--no-1mm-upfront")
    if end_to_end:
        cmd.append("--end-to-end")
    if local:
        cmd.append("--local")
    if k:
        cmd.append("-k")
        cmd.append(str(k))
    if a:
        cmd.append("-a")
    if D:
        cmd.append("-D")
        cmd.append(str(D))
    if R:
        cmd.append("-R")
        cmd.append(str(R))
    if ma:
        cmd.append("--ma")
        cmd.append(str(ma))
    if mp:
        cmd.append("--mp")
        cmd.append(",".join([str(mp[0]), str(mp[1])]))
    if np:
        cmd.append("--np")
        cmd.append(str(np))
    if rdg:
        cmd.append("--rdg")
        cmd.append(",".join([str(rdg[0]), str(rdg[1])]))
    if rfg:
        cmd.append("--rfg")
        cmd.append(",".join([str(rfg[0]), str(rfg[1])]))
    if score_min:
        cmd.append("--score-min")
        cmd.append(score_min)
    if I:
        cmd.append("-I")
        cmd.append(str(I))
    elif minins:
        cmd.append("--minins")
        cmd.append(str(minins))
    if X:
        cmd.append("-X")
        cmd.append(str(X))
    elif maxins:
        cmd.append("--maxins")
        cmd.append(str(maxins))
    if fr:
        cmd.append("--fr")
    if rf:
        cmd.append("--rf")
    if ff:
        cmd.append("--ff")
    if no_mixed:
        cmd.append("--no-mixed")
    if no_discordant:
        cmd.append("--no-discordant")
    if dovetail:
        cmd.append("--dovetail")
    if no_contain:
        cmd.append("--no-contain")
    if no_overlap:
        cmd.append("--no-overlap")
    if no_unal:
        cmd.append("--no-unal")
    if no_hd:
        cmd.append("--no-hd")
    if no_sq:
        cmd.append("--no-sq")
    if rg_id:
        cmd.append("--rg-id")
        cmd.append(rg_id)
    if rg:
        cmd.append("--rg")
        cmd.append(rg)
    if omit_sec_seq:
        cmd.append("--omit-sec-seq")
    if o:
        cmd.append("-o")
        cmd.append(str(o))
    if offrate:
        cmd.append("--offrate")
        cmd.append(str(offrate))
    if p:
        cmd.append("-p")
        cmd.append(str(p))
    if threads:
        cmd.append("--threads")
        cmd.append(str(threads))
    if reorder:
        cmd.append("--reorder")
    if mm:
        cmd.append("--mm")
    if qc_filter:
        cmd.append("--qc-filter")
    if seed:
        cmd.append("--seed")
        cmd.append(str(seed))
    if non_deterministic:
        cmd.append("--non-deterministic")

    # Creating program
    print " ".join(cmd)
    if not c:
        files_array = U and U or _1 + _2
    else:
        files_array = []
    files_array.append(x)
    program = create_program(cmd, files_array)
    return formats.Sam(S, program)
