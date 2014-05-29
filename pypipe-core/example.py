from pypipe.formats import *
from pypipe.tools import bwa, bowtie2, samtools, bcftools
from pypipe.core import pipeline

#ref = Fasta("GRCh37.fa", check=False)
#pipeline.add_file(ref)
#bwa_i = BwaIndex("GRCh37.fa", check=False) 
#pipeline.add_file(bwa_i)
#i = Bowtie2Index("GRCh37", check=False)
#pipeline.add_file(i)
#r1 = Fastq("IonXpress_021.fastq", check=False)
#pipeline.add_file(r1)
#r2 = Fastq("IonXpress_022.fastq", check=False)
#pipeline.add_file(r2)
#
#sam1 = bowtie2.bowtie2(_x=i, _U=[r1], _S="out1.sam", _p=4, log="bowtie_log")
#pipeline.draw('img1.svg')
#bam1 = samtools.view(in_=sam1, _o="out1.bam", _b=True, _S=True)
#pipeline.draw('img2.svg')
#s_bam1 = samtools.sort(in_=bam1, out="sorted1")
#pipeline.draw('img3.svg')
#
#sam2 = bwa.bwasw(ref=bwa_i, in1=r2, out="out2.sam", log="bwa_log")
#pipeline.draw('img4.svg')
#bam2 = samtools.view(in_=sam2, _o="out2.bam", _b=True, _S=True)
#pipeline.draw('img5.svg')
#s_bam2 = samtools.sort(in_=bam2, out="sorted2")
#pipeline.draw('img6.svg')
#
#bcf = samtools.mpileup(in_=[s_bam1, s_bam2], out="out.bcf", _u=True, _f=ref)
#pipeline.draw('img7.svg')
#vcf = bcftools.view(in_=bcf, out="out.vcf")
#pipeline.draw('img8.svg')

from pypipe.tools import test
r = Fastq('main.py', check=False)
test.one_one(in_=r, out='lolka')
pipeline.draw('img.svg')
