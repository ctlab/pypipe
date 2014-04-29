from pypipe.formats import *
from pypipe.tools import bwa, bowtie2, samtools, bcftools
from pypipe.utils import run_pipeline, generate_pipeline_graph


ref = Fasta("GRCh37.fa")
bwa_i = BwaIndex("GRCh37.fa")
i = Bowtie2Index("GRCh37")
r1 = Fastq("IonXpress_021.fastq")
r2 = Fastq("IonXpress_022.fastq")

sam1 = bowtie2.bowtie2(x=i, U=[r1], S="out1.sam", p=1, log="bowtie_log")
bam1 = samtools.view(in_=sam1, o="out1.bam", b=True, S=True)
s_bam1 = samtools.sort(in_=bam1, out="sorted1")

sam2 = bwa.bwasw(ref=bwa_i, in1=r2, out="out2.sam", log="bwa_log")
bam2 = samtools.view(in_=sam2, o="out2.bam", b=True, S=True)
s_bam2 = samtools.sort(in_=bam2, out="sorted2")

bcf = samtools.mpileup(in_=(s_bam1, s_bam2), out="out.bcf", u=True, f=ref)
vcf = bcftools.view(in_=bcf, out="out.vcf")

generate_pipeline_graph("pipeline-graph")
#run_pipeline()


