import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program(cmd=["git clone https://github.com/BenLangmead/bowtie2",
                     "cd bowtie2",
                     "make",
                     "mv bowtie2 bowtie2* ~/.pypipe"],
                program_name="bowtie2")


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
            qc_filter=None, seed=None, non_deterministic=None, log=None):
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
    if t and time:
        sys.exit("bowtie2: use only 't' or 'time'")

    program = create_program("bowtie2", log)
    program.add_arg(x, formats.Bowtie2Index, "-x")
    if U:
        if qseq:
            program.add_args(U, formats.Qseq, 1, ",", "-U")
        elif f:
            program.add_args(U, formats.Fasta, 1, ",", "-U")
        elif r:
            program.add_args(U, formats.TextFile, 1, ",", "-U")
        elif c:
            program.add_args(U, str, 1, ",", "-U")
        else:
            program.add_args(U, formats.Fastq, 1, ",", "-U")
    else:
        if qseq:
            program.add_args(_1, formats.Qseq, 1, ",", "-1")
            program.add_args(_2, formats.Qseq, 1, ",", "-2")
        elif f:
            program.add_args(_1, formats.Fasta, 1, ",", "-1")
            program.add_args(_2, formats.Fasta, 1, ",", "-2")
        elif r:
            program.add_args(_1, formats.TextFile, 1, ",", "-1")
            program.add_args(_2, formats.TextFile, 1, ",", "-2")
        elif c:
            program.add_args(_1, str, 1, ",", "-1")
            program.add_args(_2, str, 1, ",", "-2")
        else:
            program.add_args(_1, formats.Fastq, 1, ",", "-1")
            program.add_args(_2, formats.Fastq, 1, ",", "-2")
    program.add_arg(S, str, "-S")
    program.add_arg(q, bool, "-q")
    program.add_arg(qseq, bool, "--qseq")
    program.add_arg(f, bool, "-f")
    program.add_arg(r, bool, "-r")
    program.add_arg(c, bool, "-c")
    program.add_arg(s, int, "-s")
    program.add_arg(skip, int, "--skip")
    program.add_arg(u, int, "-u")
    program.add_arg(qupto, int, "--qupto")
    program.add_arg(trim5, int, "--trim5")
    program.add_arg(_5, int, "-5")
    program.add_arg(trim3, int, "--trim3")
    program.add_arg(_3, int, "-3")
    program.add_arg(phred33, bool, "--phred33")
    program.add_arg(phred64, bool, "--phred64")
    program.add_arg(solexa_quals, bool, "--solexa-quals")
    program.add_arg(int_quals, bool, "--int-quals")
    program.add_arg(very_fast, bool, "--very-fast")
    program.add_arg(fast, bool, "--fast")
    program.add_arg(sensitive, bool, "--sensitive")
    program.add_arg(very_sensitive, bool, "--very-sensitive")
    program.add_arg(very_fast_local, bool, "--very-fast-local")
    program.add_arg(fast_local, bool, "--fast-local")
    program.add_arg(sensitive_local, bool, "--sensitive-local")
    program.add_arg(very_sensitive_local, bool, "--very-sensitive-local")
    program.add_arg(N, int, "-N")
    program.add_arg(L, int, "-L")
    program.add_arg(i, str, "-i")
    program.add_arg(n_ceil, str, "--n-ceil")
    program.add_arg(dpad, int, "--dpad")
    program.add_arg(gbar, int, "--gbar")
    program.add_arg(ignore_quals, bool, "--ignore-quals")
    program.add_arg(nofw, bool, "--nofw")
    program.add_arg(norc, bool, "--norc")
    program.add_arg(no_1mm_upfront, bool, "--no-1mm-upfront")
    program.add_arg(end_to_end, bool, "--end-to-end")
    program.add_arg(local, bool, "--local")
    program.add_arg(k, int, "-k")
    program.add_arg(a, bool, "-a")
    program.add_arg(D, int, "-D")
    program.add_arg(R, int, "-R")
    program.add_arg(ma, int, "--ma")
    program.add_args(mp, int, 2, ",", "--mp")
    program.add_arg(np, int, "--np")
    program.add_args(rdg, int, 2, ",", "--rdg")
    program.add_args(rfg, int, 2, ",", "--rfg")
    program.add_arg(score_min, str, "--score-min")
    program.add_arg(I, int, "-I")
    program.add_arg(minins, int, "--minins")
    program.add_arg(X, int, "-X")
    program.add_arg(maxins, int, "--maxins")
    program.add_arg(fr, bool, "--fr")
    program.add_arg(rf, bool, "--rf")
    program.add_arg(ff, bool, "--ff")
    program.add_arg(no_mixed, bool, "--no-mixed")
    program.add_arg(no_discordant, bool, "--no-discordant")
    program.add_arg(dovetail, bool, "--dovetail")
    program.add_arg(no_contain, bool, "--no-contain")
    program.add_arg(no_overlap, bool, "--no-overlap")
    program.add_arg(no_unal, bool, "--no-unal")
    program.add_arg(no_hd, bool, "--no-hd")
    program.add_arg(no_sq, bool, "--no-sq")
    program.add_arg(rg_id, str, "--rg-id")
    program.add_arg(rg, str, "--rg")
    program.add_arg(omit_sec_seq, bool, "--omit-sec-seq")
    program.add_arg(o, int, "-o")
    program.add_arg(offrate, int, "--offrate")
    program.add_arg(p, int, "-p")
    program.add_arg(threads, int, "--threads")
    program.add_arg(reorder, bool, "--reorder")
    program.add_arg(mm, bool, "--mm")
    program.add_arg(qc_filter, bool, "--qc-filter")
    program.add_arg(seed, int, "--seed")
    program.add_arg(non_deterministic, bool, "--non-deterministic")
    program.add_arg(t, bool, "-t")
    program.add_arg(time, bool, "--time")
    program.add_arg(un, str, "--un")
    program.add_arg(un_gz, str, "--un-gz")
    program.add_arg(un_bz2, str, "--un-bz2")
    program.add_arg(al, str, "--al")
    program.add_arg(al_gz, str, "--al-gz")
    program.add_arg(al_bz2, str, "--al-bz2")
    program.add_arg(un_conc, str, "--un-conc")
    program.add_arg(un_conc_gz, str, "--un-conc-gz")
    program.add_arg(un_conc_bz2, str, "--un-conc-bz2")
    program.add_arg(al_conc, str, "--al-conc")
    program.add_arg(al_conc_gz, str, "--al-conc-gz")
    program.add_arg(al_conc_bz2, str, "--al-conc-bz2")
    program.add_arg(quiet, bool, "--quiet")
    program.add_arg(met_file, str, "--met-file")
    program.add_arg(met_stderr, str, "--met-stderr")
    program.add_arg(met, int, "--met")
    return formats.Sam(S, program)

