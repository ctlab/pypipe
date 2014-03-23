import os.path
import sys

class Format(object):

    def __init__(self, path, program):
        self.path = path
        self.program = program
        if not program:
            self.check()

    def check(self):
        if not os.path.isfile(self.path):
            sys.exit("File " + self.path + " doesn't exist")


class Unknown(Format):

    def __init__(self, path, program=None):
        super(Unknown, self).__init__(path, program)


class Fastq(Format):
    
    def __init__(self, path, program=None):
        super(Fastq, self).__init__(path, program)


class Fasta(Format):
    
    def __init__(self, path, program=None):
        super(Fasta, self).__init__(path, program)


class Bam(Format):
    
    def __init__(self, path, program=None):
        super(Bam, self).__init__(path, program)


class Sam(Format):
    
    def __init__(self, path, program=None):
        super(Sam, self).__init__(path, program)


class Bai(Format):
    
    def __init__(self, path, program=None):
        super(Bai, self).__init__(path, program)


class Sai(Format):
    
    def __init__(self, path, program=None):
        super(Sai, self).__init__(path, program)


class Bcf(Format):
    
    def __init__(self, path, program=None):
        super(Bcf, self).__init__(path, program)


class Vcf(Format):
    
    def __init__(self, path, program=None):
        super(Vcf, self).__init__(path, program)


class Bed(Format):
    
    def __init__(self, path, program=None):
        super(Bed, self).__init__(path, program)


class Bowtie2Index(Format):

    def __init__(self, path, program=None):
        super(Bowtie2Index, self).__init__(path, program)

    def check(self):
        files = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"]
        files = [self.path + f + ".bt2" for f in files]
        for f in files:
            if not os.path.isfile(f):
                sys.exit("File " + f + " doesn't exist")


