from formats import *
from tools import bowtie2, samtools, bcftools
from utils import run_pipeline 


ref = Fasta("GRCh37.fa")
i = Bowtie2Index("GRCh37")
r1 = [Fastq("IonXpress_021.fastq")]
r2 = [Fastq("IonXpress_022.fastq")]

sam1 = bowtie2.bowtie2(x=i, U=r1, S="out1.sam", p=6)
bam1 = samtools.view(aln=sam1, o="out1.bam", b=True, S=True)
s_bam1 = samtools.sort(aln=bam1, out="sorted1")

sam2 = bowtie2.bowtie2(x=i, U=r2, S="out2.sam", p=1)
bam2 = samtools.view(aln=sam2, o="out2.bam", b=True, S=True)
s_bam2 = samtools.sort(aln=bam2, out="sorted2")

bcf = samtools.mpileup(_in=(s_bam1, s_bam2), out="out.bcf", u=True, f=ref)
vcf = bcftools.view(_in=bcf, out="out.vcf")

run_pipeline()
