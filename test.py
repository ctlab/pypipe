from pypipe.tools import test

r1 = Fastq('pypipe.py')
r2 = Fastq('pypipe.py')

fastq = test.one_one(in_=r1, out='ok')
fastq1, fastq2 = test.one_two(in_=fastq, out='ok')

result = test.three_one(in_=[fastq, fastq1, fastq2], out='result')

#from pypipe.core import pipeline
#pipeline.load('db')
#pipeline.run(2)
#pipeline.save('db')
#pipeline.load('db')
#pipeline.run(2)