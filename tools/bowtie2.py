import sys

import formats
from pipeline import create_program


def unpaired(index, reads, output, params=None):
    """ bowtie2 [options]* -x <bt2-idx> -U r -S <res>

    For example: bowtie2 -x "index" -U read1.fq, read2.fq -S out.bam
    index = BowtieIndex("index")
    reads = [Fastq("read1.fq"), Fastq("read2.fq")]
    out_sam = bowtie2.unpaired(index=index, reads=reads, output="out.sam")

    This function returns SAM file

    """
    for read in reads:
        if type(read) != formats.Fastq:
            msg = "bowtie2: " + read.path + " must be in FASTQ format"
            sys.exit(msg)
    reads_arg = ",".join([read.path for read in reads])
    cmd = ["bowtie2", "-x", index.path, "-U", reads_arg, "-S", output] 
    # params ...
    program = create_program(cmd, reads)
    return formats.Sam(output, program)


def paired(index, reads1, reads2, output, params=None):
    """ bowtie2 [options]* -x <bt2-idx> -1 <m1> -2 <m2> -S <res>

    For example: bowtie2 -x "index" -1 "read1.fq" -2 "read2.fq" -S "out.bam"
    index = BowtieIndex("index")
    reads1 = [Fastq("read1.fq")]
    reads2 = [Fastq("read2.fq")]
    bowtie2.paired(index=index, reads1=reads1, reads2=reads2, output="out.sam")

    This function returns SAM file

    """
    reads = reads1 + reads2
    for read in reads:
        if type(read) != formats.Fastq:
            msg = "bowtie2: " + read.path + " must be in FASTQ format"
            sys.exit(msg)
    reads1_arg = ",".join([read.path for read in reads1])
    reads2_arg = ",".join([read.path for read in reads2])
    cmd = ["bowtie2", "-x", index.path, "-1", reads1_arg,
           "-2", reads2_arg, "-S", output]
    # params ...
    program = create_program(cmd, reads)
    return formats.Sam(output, program)

