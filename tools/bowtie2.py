import sys

import formats
from pipeline import create_program


# TODO -r - unknown format
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

    program = create_program("bowtie2")
    program.add_argument(x, formats.Bowtie2Index, "-x")
    if U:
        if qseq:
            program.add_arguments(U, formats.Qseq, 1, ",", "-U")
        elif f:
            program.add_arguments(U, formats.Fasta, 1, ",", "-U")
        elif r:
            program.add_arguments(U, formats.Unknown, 1, ",", "-U")
        elif c:
            program.add_arguments(U, str, 1, ",", "-U")
        else:
            program.add_arguments(U, formats.Fastq, 1, ",", "-U")
    else:
        if qseq:
            program.add_arguments(_1, formats.Qseq, 1, ",", "-1")
            program.add_arguments(_2, formats.Qseq, 1, ",", "-2")
        elif f:
            program.add_arguments(_1, formats.Fasta, 1, ",", "-1")
            program.add_arguments(_2, formats.Fasta, 1, ",", "-2")
        elif r:
            program.add_arguments(_1, formats.Unknown, 1, ",", "-1")
            program.add_arguments(_2, formats.Unknown, 1, ",", "-2")
        elif c:
            program.add_arguments(_1, str, 1, ",", "-1")
            program.add_arguments(_2, str, 1, ",", "-2")
        else:
            program.add_arguments(_1, formats.Fastq, 1, ",", "-1")
            program.add_arguments(_2, formats.Fastq, 1, ",", "-2")
    program.add_argument(S, str, "-S")
    program.add_argument(q, bool, "-q")
    program.add_argument(qseq, bool, "--qseq")
    program.add_argument(f, bool, "-f")
    program.add_argument(r, bool, "-r")
    program.add_argument(c, bool, "-c")
    program.add_argument(s, int, "-s")
    program.add_argument(skip, int, "--skip")
    program.add_argument(u, int, "-u")
    program.add_argument(qupto, int, "--qupto")
    program.add_argument(trim5, int, "--trim5")
    program.add_argument(_5, int, "-5")
    program.add_argument(trim3, int, "--trim3")
    program.add_argument(_3, int, "-3")
    program.add_argument(phred33, bool, "--phred33")
    program.add_argument(phred64, bool, "--phred64")
    program.add_argument(solexa_quals, bool, "--solexa-quals")
    program.add_argument(int_quals, bool, "--int-quals")
    program.add_argument(very_fast, bool, "--very-fast")
    program.add_argument(fast, bool, "--fast")
    program.add_argument(sensitive, bool, "--sensitive")
    program.add_argument(very_sensitive, bool, "--very-sensitive")
    program.add_argument(very_fast_local, bool, "--very-fast-local")
    program.add_argument(fast_local, bool, "--fast-local")
    program.add_argument(sensitive_local, bool, "--sensitive-local")
    program.add_argument(very_sensitive_local, bool, "--very-sensitive-local")
    program.add_argument(N, int, "-N")
    program.add_argument(L, int, "-L")
    program.add_argument(i, str, "-i")
    program.add_argument(n_ceil, str, "--n-ceil")
    program.add_argument(dpad, int, "--dpad")
    program.add_argument(gbar, int, "--gbar")
    program.add_argument(ignore_quals, bool, "--ignore-quals")
    program.add_argument(nofw, bool, "--nofw")
    program.add_argument(norc, bool, "--norc")
    program.add_argument(no_1mm_upfront, bool, "--no-1mm-upfront")
    program.add_argument(end_to_end, bool, "--end-to-end")
    program.add_argument(local, bool, "--local")
    program.add_argument(k, int, "-k")
    program.add_argument(a, bool, "-a")
    program.add_argument(D, int, "-D")
    program.add_argument(R, int, "-R")
    program.add_argument(ma, int, "--ma")
    program.add_arguments(mp, int, 2, ",", "--mp")
    program.add_argument(np, int, "--np")
    program.add_arguments(rdg, int, 2, ",", "--rdg")
    program.add_arguments(rfg, int, 2, ",", "--rfg")
    program.add_argument(score_min, str, "--score-min")
    program.add_argument(I, int, "-I")
    program.add_argument(minins, int, "--minins")
    program.add_argument(X, int, "-X")
    program.add_argument(maxins, int, "--maxins")
    program.add_argument(fr, bool, "--fr")
    program.add_argument(rf, bool, "--rf")
    program.add_argument(ff, bool, "--ff")
    program.add_argument(no_mixed, bool, "--no-mixed")
    program.add_argument(no_discordant, bool, "--no-discordant")
    program.add_argument(dovetail, bool, "--dovetail")
    program.add_argument(no_contain, bool, "--no-contain")
    program.add_argument(no_overlap, bool, "--no-overlap")
    program.add_argument(no_unal, bool, "--no-unal")
    program.add_argument(no_hd, bool, "--no-hd")
    program.add_argument(no_sq, bool, "--no-sq")
    program.add_argument(rg_id, str, "--rg-id")
    program.add_argument(rg, str, "--rg")
    program.add_argument(omit_sec_seq, bool, "--omit-sec-seq")
    program.add_argument(o, int, "-o")
    program.add_argument(offrate, int, "--offrate")
    program.add_argument(p, int, "-p")
    program.add_argument(threads, int, "--threads")
    program.add_argument(reorder, bool, "--reorder")
    program.add_argument(mm, bool, "--mm")
    program.add_argument(qc_filter, bool, "--qc-filter")
    program.add_argument(seed, int, "--seed")
    program.add_argument(non_deterministic, int, "--non-deterministic")
    print " ".join(program.cmd)
    return formats.Sam(S, program)
