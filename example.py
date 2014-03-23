from formats import *
from tools import bwa, bowtie2, samtools
from pipeline import run_pipeline


ref = Fasta("GRCh37.fa")
i = Bowtie2Index("GRCh37")
r1 = [Fastq("IonXpress_021.fastq")]
r2 = [Fastq("IonXpress_022.fastq")]

sam1 = bowtie2.bowtie2(x=i, U=r1, S="out1.sam", p=6)
bam1 = samtools.view(align=sam1, o="out1.bam", b=True, S=True)
s_bam1 = samtools.sort(align=bam1, out="sorted1")
bcf1 = samtools.mpileup(align=[s_bam1], out="out1.bcf", u=True, f=ref)

sam2 = bowtie2.bowtie2(x=i, U=r2, S="out2.sam", p=1)
bam2 = samtools.view(align=sam2, o="out2.bam", b=True, S=True)
s_bam2 = samtools.sort(align=bam2, out="sorted2")
bcf2 = samtools.mpileup(align=[s_bam2], out="out2.bcf", u=True, f=ref)

run_pipeline()
