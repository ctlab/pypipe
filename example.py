from formats import *
from tools import bwa, bowtie2, samtools, bcftools
from utils import run_pipeline 


ref = Fasta("GRCh37.fa")
i = Bowtie2Index("GRCh37")
r1 = Fastq("IonXpress_021.fastq")
r2 = Fastq("IonXpress_022.fastq")

sam1 = bowtie2.bowtie2(x=i, U=[r1], S="out1.sam", p=1)
bam1 = samtools.view(_in=sam1, o="out1.bam", b=True, S=True)
s_bam1 = samtools.sort(_in=bam1, _out="sorted1")

sam2 = bwa.bwasw(_ref=ref, _in1=r2, _out="out2.sam")
bam2 = samtools.view(_in=sam2, o="out2.bam", b=True, S=True)
s_bam2 = samtools.sort(_in=bam2, _out="sorted2")

bcf = samtools.mpileup(_in=(s_bam1, s_bam2), _out="out.bcf", u=True, f=ref)
vcf = bcftools.view(_in=bcf, _out="out.vcf")

run_pipeline()
