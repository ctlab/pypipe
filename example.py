from formats import Fasta, Fastq
from programs import bwa, samtools, test
from pipeline import run_pipeline

#bwa.bwasw(Fasta("ref"), Fastq("read"), "out.sam")
#sam1 = test.touch("sam1")
#sam2 = test.touch("sam2")
#cp1 = test.cp(sam1)
#cp2 = test.cp(sam2)
seven = test.sleep()
test.cp(seven, Fasta("some_file"))
run_pipeline()
