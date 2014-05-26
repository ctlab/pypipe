import abstractfile


class TextFile(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(TextFile, self).__init__(path, program, check)


class Fastq(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Fastq, self).__init__(path, program, check)


class Fasta(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Fasta, self).__init__(path, program, check)


class Bam(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Bam, self).__init__(path, program, check)


class Sam(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Sam, self).__init__(path, program, check)


class Bai(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Bai, self).__init__(path, program, check)


class Sai(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Sai, self).__init__(path, program, check)


class Fai(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Fai, self).__init__(path, program, check)


class Bcf(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Bcf, self).__init__(path, program, check)


class Vcf(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Vcf, self).__init__(path, program, check)


class Bed(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Bed, self).__init__(path, program, check)


class Qseq(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Qseq, self).__init__(path, program, check)


class Snp(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Snp, self).__init__(path, program, check)


class Pileup(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Pileup, self).__init__(path, program, check)


class Indel(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Indel, self).__init__(path, program, check)


class Cns(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Cns, self).__init__(path, program, check)


class Bowtie2Index(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(Bowtie2Index, self).__init__(path, program, check,
            suff=['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'])


class BwaIndex(abstractfile.File):

    def __init__(self, path, program=None, check=True):
        super(BwaIndex, self).__init__(path, program, check,
            suff=['.amb', '.ann', '.bwt', '.pac', '.sa'])
