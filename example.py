from pypipe.tools import bwa, bowtie2, samtools, bcftools


set_input_dir('~/pypipe/input')

# Input files
ref = Fasta("GRCh37.fa")  # reference file
bwa_i = BwaIndex("GRCh37.fa")  # all bwa index files
i = Bowtie2Index("GRCh37")  # all bowtie2 index files
r1 = Fastq("IonXpress_021.fastq")  # read 1
r2 = Fastq("IonXpress_022.fastq")  # read 2

# bowtie2 -x GRCh37 -U IonXpress_021.fastq -S out1.sam -p 1
# samtools view -bS -o out1.bam out1.sam
# samtools sort out1.bam > sorted1
sam1 = bowtie2.bowtie2(_x=i, _U=[r1], _S="out1.sam", _p=4, log="bowtie_log")
bam1 = samtools.view(in_=sam1, _o="out1.bam", _b=True, _S=True)
s_bam1 = samtools.sort(in_=bam1, out="sorted1")

# bwa bwasw GRCh37.fa IonXpress_022.fastq > out2.sam
# samtools view -bS -o out2.bam out2.sam
# samtools sort out2.bam > sorted2
sam2 = bwa.bwasw(ref=bwa_i, in1=r2, out="out2.sam", log="bwa_log")
bam2 = samtools.view(in_=sam2, _o="out2.bam", _b=True, _S=True)
s_bam2 = samtools.sort(in_=bam2, out="sorted2")

# samtools mpileup -f GRCh37.fa -u sorted1.bam sorted2.bam > out.bcf
# samtools view out.bcf > out.vcf
bcf = samtools.mpileup(in_=[s_bam1, s_bam2], out="out.bcf", _u=True, _f=ref)
vcf = bcftools.view(in_=bcf, out="out.vcf")

