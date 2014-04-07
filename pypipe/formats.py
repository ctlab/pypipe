import os.path
import sys


class _File(object):

    def __init__(self, path, program):
        self.path = path
        self.program = program
        if not program:
            self.check()

    def check(self):
        if not os.path.isfile(self.path):
            sys.exit("File " + self.path + " doesn't exist")


class TextFile(_File):

    def __init__(self, path, program=None):
        super(TextFile, self).__init__(path, program)


class Fastq(_File):
    
    def __init__(self, path, program=None):
        super(Fastq, self).__init__(path, program)


class Fasta(_File):
    
    def __init__(self, path, program=None):
        super(Fasta, self).__init__(path, program)


class Bam(_File):
    
    def __init__(self, path, program=None):
        super(Bam, self).__init__(path, program)


class Sam(_File):
    
    def __init__(self, path, program=None):
        super(Sam, self).__init__(path, program)


class Bai(_File):
    
    def __init__(self, path, program=None):
        super(Bai, self).__init__(path, program)


class Sai(_File):
    
    def __init__(self, path, program=None):
        super(Sai, self).__init__(path, program)


class Fai(_File):
    
    def __init__(self, path, program=None):
        super(Fai, self).__init__(path, program)


class Bcf(_File):
    
    def __init__(self, path, program=None):
        super(Bcf, self).__init__(path, program)


class Vcf(_File):
    
    def __init__(self, path, program=None):
        super(Vcf, self).__init__(path, program)


class Bed(_File):
    
    def __init__(self, path, program=None):
        super(Bed, self).__init__(path, program)


class Qseq(_File):
    
    def __init__(self, path, program=None):
        super(Qseq, self).__init__(path, program)


class Snp(_File):
    
    def __init__(self, path, program=None):
        super(Snp, self).__init__(path, program)


class Pileup(_File):
    
    def __init__(self, path, program=None):
        super(Pileup, self).__init__(path, program)


class Bowtie2Index(_File):

    def __init__(self, path, program=None):
        super(Bowtie2Index, self).__init__(path, program)

    def check(self):
        files = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"]
        files = [self.path + f + ".bt2" for f in files]
        for f in files:
            if not os.path.isfile(f):
                sys.exit("File " + f + " doesn't exist")


class BwaIndex(_File):

    def __init__(self, path, program=None):
        super(BwaIndex, self).__init__(path, program)

    def check(self):
        files = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        files = [self.path + f for f in files]
        for f in files:
            if not os.path.isfile(f):
                sys.exit("File " + f + " doesn't exist")

