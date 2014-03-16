from formats import *
from tools import bwa, samtools
from pipeline import run_pipeline

ref = Fasta("GRCh37.fa")
read1 = Fastq("IonXpress_021.fastq")
#read2 = Fastq("IonXpress_022.fastq")
#sam1 = bwa.bwasw(ref=ref, read=read1, output="1.sam")
#sam2 = bwa.bwasw(ref=ref, read=read2, output="2.sam")
#index = Sai("aln_sa.sai")
#bwa.samse(ref=ref, index=index, read=read1, output="out.sam", params="-n 8")
bam = Bam("out.bam")
samtools.sort(align=bam, output="sorted")
run_pipeline()
