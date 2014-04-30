from pypipe.formats import *
from pypipe.tools import bwa, bowtie2, samtools, bcftools
from pypipe.utils import run_pipeline, generate_pipeline_graph


# Input files
ref = Fasta("GRCh37.fa")  # reference file
bwa_i = BwaIndex("GRCh37.fa")  # all bwa index files
i = Bowtie2Index("GRCh37")  # all bowtie2 index files
r1 = Fastq("IonXpress_021.fastq")  # read 1
r2 = Fastq("IonXpress_022.fastq")  # read 2

# bowtie2 -x GRCh37 -U IonXpress_021.fastq -S out1.sam -p 1
sam1 = bowtie2.bowtie2(x=i, U=[r1], S="out1.sam", p=1, log="bowtie_log")
# samtools view -bS -o out1.bam out1.sam
bam1 = samtools.view(in_=sam1, o="out1.bam", b=True, S=True)
# samtools sort out1.bam > sorted1
s_bam1 = samtools.sort(in_=bam1, out="sorted1")

# bwa bwasw GRCh37.fa IonXpress_022.fastq > out2.sam
sam2 = bwa.bwasw(ref=bwa_i, in1=r2, out="out2.sam", log="bwa_log")
# samtools view -bS -o out2.bam out2.sam
bam2 = samtools.view(in_=sam2, o="out2.bam", b=True, S=True)
# samtools sort out2.bam > sorted2
s_bam2 = samtools.sort(in_=bam2, out="sorted2")

# samtools mpileup -f GRCh37.fa -u sorted1.bam, sorted2.bam > out.bcf
bcf = samtools.mpileup(in_=(s_bam1, s_bam2), out="out.bcf", u=True, f=ref)
# samtools view out.bcf > out.vcf
vcf = bcftools.view(in_=bcf, out="out.vcf")


#generate_pipeline_graph("pipeline-graph")  # Create visualization picture
run_pipeline(s_bam1)  # Run pipeline

