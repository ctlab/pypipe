from formats import *
from tools import bwa, bowtie2, samtools
from pipeline import run_pipeline


#ref = Fasta("GRCh37.fa")
i = Bowtie2Index("GRCh37")
r1 = [Fastq("IonXpress_021.fastq")]
r2 = [Fastq("IonXpress_022.fastq")]

sam1 = bowtie2.bowtie2(x=i, U=r1, S="OUT1.sam")
bam1 = samtools.view(align=sam1, o="OUT1.bam", b=True, S=True)
sorted_bam1 = samtools.sort(align=bam1, out="SORTED1")

sam2 = bowtie2.bowtie2(x=i, U=r2, S="OUT2.sam")
bam2 = samtools.view(align=sam2, o="OUT2.bam", b=True, S=True)
sorted_bam2 = samtools.sort(align=bam2, out="SORTED2")

run_pipeline()
