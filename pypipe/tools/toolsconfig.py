from pypipe import formats


class Bcftools:

    @staticmethod
    def view():
        return {
            'cmd': 'bcftools view',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'b': formats.Bcf, '': formats.Vcf}, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '-A': bool,
                    '-b': bool,
                    '-D': formats.TextFile,
                    '-F': bool,
                    '-G': bool,
                    '-l': formats.TextFile,
                    '-N': bool,
                    '-Q': bool,
                    '-s': formats.TextFile,
                    '-S': bool,
                    '-u': bool,
                    '-c': bool,
                    '-d': float,
                    '-e': bool,
                    '-g': bool,
                    '-i': float,
                    '-p': float,
                    '-P': str,
                    '-t': float,
                    '-T': str,
                    '-v': bool,
                    '-1': int,
                    '-U': int,
                    '-X': float,
                },
                'unnamed': [
                    ('in_*', {'-S': formats.Vcf, '': formats.Bcf}),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def cat():
        return {
            'cmd': 'bcftools cat',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Bcf, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                },
                'unnamed': {
                    'in_*': formats.Bcf,
                    'out*': str,
                },
            }
        }


class Bowtie2:

    @staticmethod
    def bowtie2():
        return {
            'cmd': 'bowtie2',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': '-S', 'type': formats.Sam, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '-x*': formats.Bowtie2Index,
                    '-S*': str,
                    '-U': {
                        '--qseq': [formats.Qseq, 1, ','],
                        '-f':     [formats.Fasta, 1, ','],
                        '-r':     [formats.TextFile, 1, ','],
                        '-c':     str,
                        '':       [formats.Fastq, 1, ',']
                        },
                    '-1': {
                        '--qseq': [formats.Qseq, 1, ','],
                        '-f':     [formats.Fasta, 1, ','],
                        '-r':     [formats.TextFile, 1, ','],
                        '-c':     str,
                        '':       [formats.Fastq, 1, ',']
                        },
                    '-2': {
                        '--qseq': [formats.Qseq, 1, ','],
                        '-f':     [formats.Fasta, 1, ','],
                        '-r':     [formats.TextFile, 1, ','],
                        '-c':     str,
                        '':       [formats.Fastq, 1, ',']
                        },
                    '-q': bool,
                    '--qseq': bool,
                    '-f': bool,
                    '-r': bool,
                    '-c': bool,
                    '-s': int,
                    '-u': int,
                    '-5': int,
                    '-3': int,
                    '--phred33': bool,
                    '--phred64': bool,
                    '--solexa-quals': bool,
                    '--int-quals': bool,
                    '--very-fast': bool,
                    '--fast': bool,
                    '--sensitive': bool,
                    '--very-sensitive': bool,
                    '--very-fast-local': bool,
                    '--fast-local': bool,
                    '--sensitive-local': bool,
                    '--very-sensitive-local': bool,
                    '-N': int,
                    '-L': int,
                    '-i': str,
                    '--n-ceil': str,
                    '--dpad': int,
                    '--gbar': int,
                    '--ignore-quals': bool,
                    '--nofw': bool,
                    '--norc': bool,
                    '--no-1mm-upfront': bool,
                    '--end-to-end': bool,
                    '--local': bool,
                    '-k': int,
                    '-a': bool,
                    '-D': int,
                    '-R': int,
                    '--ma': int,
                    '--mp': [int, 2, ','],
                    '--np': int,
                    '--rdg': [int, 2, ','],
                    '--rfg': [int, 2, ','],
                    '--score-min': str,
                    '-I': int,
                    '-X': int,
                    '--fr': bool,
                    '--rf': bool,
                    '--ff': bool,
                    '--no-mixed': bool,
                    '--no-discordant': bool,
                    '--dovetail': bool,
                    '--no-contain': bool,
                    '--no-overlap': bool,
                    '--no-unal': bool,
                    '--no-hd': bool,
                    '--no-sq': bool,
                    '--rg-id': str,
                    '--rg': str,
                    '--omit-sec-seq': bool,
                    '-o': int,
                    '-p': int,
                    '--reorder': bool,
                    '--mm': bool,
                    '--qc-filter': bool,
                    '--seed': int,
                    '--non-deterministic': bool,
                    '-t': bool,
                    '--un': str,
                    '--un-gz': str,
                    '--un-bz2': str,
                    '--al': str,
                    '--al-gz': str,
                    '--al-bz2': str,
                    '--un-conc': str,
                    '--un-conc-gz': str,
                    '--un-conc-bz2': str,
                    '--al-conc': str,
                    '--al-conc-gz': str,
                    '--al-conc-bz2': str,
                    '--quiet': bool,
                    '--met-file': str,
                    '--met-stderr': str,
                    '--met': int,
                },
                'unnamed': [
                ],
            },
        }


class Bwa:

    @staticmethod
    def mem():
        return {
            'cmd': 'bwa mem',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Sam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-t': int,
                    '-k': int,
                    '-w': int,
                    '-d': int,
                    '-r': float,
                    '-c': int,
                    '-P': bool,
                    '-a': int,
                    '-B': int,
                    '-O': int,
                    '-E': int,
                    '-L': int,
                    '-U': int,
                    '-p': bool,
                    '-R': str,
                    '-T': int,
                    '-C': bool,
                    '-H': bool,
                    '-M': bool,
                    '-v': int,
                },
                'unnamed': [
                    ('ref*', formats.BwaIndex),
                    ('in1*', formats.Fastq),
                    ('in2', formats.Fastq),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def aln():
        return {
            'cmd': 'bwa aln',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Sai, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-n': int,
                    '-o': int,
                    '-e': int,
                    '-d': int,
                    '-i': int,
                    '-l': int,
                    '-k': int,
                    '-t': int,
                    '-M': int,
                    '-O': int,
                    '-E': int,
                    '-R': int,
                    '-c': bool,
                    '-N': bool,
                    '-q': int,
                    '-I': bool,
                    '-B': int,
                    '-b': bool,
                    '-0': bool,
                    '-1': bool,
                    '-2': bool,
                },
                'unnamed': [
                    ('ref*', formats.BwaIndex),
                    ('in_*', formats.Fastq),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def samse():
        return {
            'cmd': 'bwa samse',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Sam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-n': int,
                    '-r': str,
                },
                'unnamed': [
                    ('ref*', formats.BwaIndex),
                    ('sai*', formats.Sai),
                    ('in_*', formats.Fastq),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def sampe():
        return {
            'cmd': 'bwa sampe',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Sam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-a': int,
                    '-o': int,
                    '-P': bool,
                    '-n': int,
                    '-N': int,
                    '-r': str,
                },
                'unnamed': [
                    ('ref*', formats.BwaIndex),
                    ('sai1*', formats.Sai),
                    ('sai2*', formats.Sai),
                    ('in1*', formats.Fastq),
                    ('in2*', formats.Fastq),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def bwasw():
        return {
            'cmd': 'bwa bwasw',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Sam, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '-a': int,
                    '-b': int,
                    '-q': int,
                    '-r': int,
                    '-t': int,
                    '-w': int,
                    '-T': int,
                    '-c': float,
                    '-z': int,
                    '-s': int,
                    '-N': int,
                },
                'unnamed': [
                    ('ref*', formats.BwaIndex),
                    ('in1*', formats.Fastq),
                    ('in2', formats.Fastq),
                    ('out*', str),
                ],
            },
        }

class Freebayes:

    @staticmethod
    def freebayes():
        return {
            'cmd': 'freebayes',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'v', 'type': formats.Vcf, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '-v*': str,
                    '-f*': formats.Fasta,
                    '-b': formats.Bam,
                    '-t': formats.Bed,
                    '-r': str,
                    '-s': formats.TextFile,
                    '--populations': formats.TextFile,
                    '-A': formats.Bed,
                    '--trace': str,
                    '--failed-alleles': formats.Bed,
                    '--variant-input': formats.Vcf,
                    '-l': bool,
                    '--haplotype-basis-alleles': bool,
                    '--report-all-haplotype-alleles': bool,
                    '--report-monorphic': bool,
                    '-P': float,
                    '-T': float,
                    '-p': int,
                    '-J': bool,
                    '-Z': bool,
                    '--reference_quality': [int, 2, ','],
                    '-I': bool,
                    '-i': bool,
                    '-X': bool,
                    '-u': bool,
                    '-n': int,
                    '-E': int,
                    '--max-complex-gap': int,
                    '--haplotype-length': int,
                    '--min-repeat-length': int,
                    '--min-repeat-entropy': int,
                    '--no-partial-observation': bool,
                    '-O': bool,
                    '-4': bool,
                    '-m': int,
                    '-q': int,
                    '-R': int,
                    '-Y': int,
                    '-Q': int,
                    '-U': int,
                    '-z': int,
                    '--read-snp-limit': int,
                    '-e': int,
                    '-0': int,
                    '-F': int,
                    '-C': int,
                    '-3': int,
                    '-G': int,
                    '--min-coverage': int,
                    '-w': bool,
                    '-V': bool,
                    '-a': bool,
                    '--observation-bias': formats.TextFile,
                    '--base-quality-cap': int,
                    '--experimental-gls': bool,
                    '--prob-contamination': float,
                    '--contamination-estimates': formats.TextFile,
                    '--report-genotype-likelihood-max': bool,
                    '-B': int,
                    '--genotype-max-iterations': int,
                    '-W': [int, 2, ','],
                    '-S': int,
                    '-j': bool,
                    '-H': bool,
                    '-D': int,
                    '--genotype-qualities': bool,
                },
                'unnamed': [
                    ('in_*', [formats.Bam, 1, ' ']),
                ]
            }
        }

class Samtools:

    @staticmethod
    def view():
        return {
            'cmd': 'samtools view',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': '-o', 'type': {'-b': formats.Bam, '': formats.Sam}, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '-o*': str,
                    '-b': bool,
                    '-f': int,
                    '-F': int,
                    '-h': bool,
                    '-H': bool,
                    '-l': str,
                    '-q': int,
                    '-r': str,
                    '-R': formats.Bed,
                    '-S': bool,
                    '-c': bool,
                    '-t': formats.TextFile,
                    '-u': bool,
                },
                'unnamed': [
                    ('in_*', {'-S': formats.Sam, '': formats.Bam}),
                ],
            },
        }

    @staticmethod
    def mpileup():
        return {
            'cmd': 'samtools mpileup',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'-u': formats.Bcf, '-g': formats.Bcf, '': formats.Pileup}, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-6': bool,
                    '-A': bool,
                    '-B': bool,
                    '-b': formats.TextFile,
                    '-C': int,
                    '-E': bool,
                    '-f': formats.Fasta,
                    '-l': formats.Bed,
                    '-q': int,
                    '-Q': int,
                    '-r': str,
                    '-D': bool,
                    '-g': bool,
                    '-S': bool,
                    '-u': bool,
                    '-e': int,
                    '-h': int,
                    '-I': bool,
                    '-L': bool,
                    '-o': int,
                    '-P': str,
                },
                'unnamed': [
                    ('in_*', [formats.Bam, 1, ' ']),
                    ('out*', str),
                ]
            }
        }

    @staticmethod
    def cat():
        return {
            'cmd': 'samtools cat',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'o', 'type': formats.Bam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-h': formats.Sam,
                },
                'unnamed': [
                    ('in_*', formats.Bam),
                    ('o*', str)
                ],
            },
        }

    @staticmethod
    def sort():
        return {
            'cmd': 'samtools sort',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': formats.Bam, 'suffix': '.bam'},
                ],
            },
            'args': {
                'named': {
                    '-n': bool,
                    '-m': int,
                },
                'unnamed': [
                    ('in_*', formats.Bam),
                    ('out*', str)
                ],
            },
        }

    @staticmethod
    def merge():
        return {
            'cmd': 'samtools merge',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': formats.Bam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-1': bool,
                    '-f': bool,
                    '-h': formats.Sam,
                    '-n': bool,
                    '-R': str,
                    '-r': bool,
                    '-u': bool,
                },
                'unnamed': [
                    ('out*', str),
                    ('in_*', [formats.Bam, 2, ' ']),
                ],
            },
        }

    @staticmethod
    def rmdup():
        return {
            'cmd': 'samtools rmdup',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': formats.Bam, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-s': bool,
                    '-S': bool,
                },
                'unnamed': [
                    ('in_*', formats.Bam),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def calmd():
        return {
            'cmd': 'samtools calmd',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'-u': formats.Bam, 'b': formats.Bam, '': formats.Sam}, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-A': bool,
                    '-e': bool,
                    '-u': bool,
                    '-b': bool,
                    '-S': bool,
                    '-C': int,
                    '-r': bool,
                    '-E': bool,
                },
                'unnamed': [
                    ('in_*', {'-S': formats.Sam, '': formats.Bam}),
                    ('out*', str),
                ],
            },
        }

class Varscan:

    @staticmethod
    def pileup2snp():
        return {
            'cmd': 'java -jar VarScan.jar pileup2snp',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Snp, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--p-value': float,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def pileup2indel():
        return {
            'cmd': 'java -jar VarScan.jar pileup2indel',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Indel, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--p-value': float,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def pileup2cns():
        return {
            'cmd': 'java -jar VarScan.jar pileup2cns',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': formats.Cns, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--p-value': float,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def mpileup2snp():
        return {
            'cmd': 'java -jar VarScan.jar mpileup2snp',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'--output-vcf': formats.Vcf, '': formats.Snp}, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--min-freq-for-hom': float,
                    '--p-value': float,
                    '--strand-filter': int,
                    '--output-vcf': bool,
                    '--variants': int,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def mpileup2indel():
        return {
            'cmd': 'java -jar VarScan.jar mpileup2indel',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'--output-vcf': formats.Vcf, '': formats.Snp}, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--min-freq-for-hom': float,
                    '--p-value': float,
                    '--strand-filter': int,
                    '--output-vcf': bool,
                    '--variants': int,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def mpileup2cns():
        return {
            'cmd': 'java -jar VarScan.jar mpileup2cns',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': True,
                'return': [
                    {'arg': 'out', 'type': {'--output-vcf': formats.Vcf, '': formats.Snp}, 'suffix': ''},
                ]
            },
            'args': {
                'named': {
                    '--min-coverage': int,
                    '--min-reads2': int,
                    '--min-avg-qual': int,
                    '--min-var-freq': float,
                    '--min-freq-for-hom': float,
                    '--p-value': float,
                    '--strand-filter': int,
                    '--output-vcf': bool,
                    '--variants': int,
                },
                'unnamed': [
                    ('in_*', formats.Pileup),
                    ('out*', str),
                ],
            },
        }

    @staticmethod
    def somatic():
        return {
            'cmd': 'java -jar VarScan.jar somatic',
            'type': 'jar',
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'output', 'type': {'--output-snp': None, '': formats.Snp}, 'suffix': '.indel'},
                    {'arg': 'output', 'type': {'--output-indel': None, '': formats.Indel}, 'suffix': '.snp'},
                ]
            },
            'args': {
                'named': {
                    '--output-snp': bool,
                    '--output-indel': bool,
                    '--min-coverage': int,
                    '--min-coverage-normal': int,
                    '--min-coverage-tumor': int,
                    '--min-var-freq': float,
                    '--min-freq-for-hom': float,
                    '--normal-purity': float,
                    '--tumor-purity': float,
                    '--p-value': float,
                    '--somatic-p-value': float,
                    '--strand-filter': float,
                    '--validation': float,
                },
                'unnamed': [
                    ('normal*', formats.Pileup),
                    ('tumor*', formats.Pileup),
                    ('output*', str),
                ],
            },
        }


class Test:

    @staticmethod
    def one_one():
        return {
            'cmd': 'true',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': {'-1': formats.Fasta, '': formats.Fastq}, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-1': bool,
                },
                'unnamed': [
                    ('in_*', formats.Fastq),
                    ('out*', str),
                ],
            }
        }

    @staticmethod
    def one_two():
        return {
            'cmd': 'true',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': {'-1': formats.Fasta, '': formats.Fastq}, 'suffix': '1'},
                    {'arg': 'out', 'type': {'-1': formats.Fasta, '': formats.Fastq}, 'suffix': '2'},
                ],
            },
            'args': {
                'named': {
                    '-1': bool,
                },
                'unnamed': [
                    ('in_*', formats.Fastq),
                    ('out*', str),
                ],
            }
        }

    @staticmethod
    def three_one():
        return {
            'cmd': 'true',
            'type': None,
            'log': 'log',
            'out': {
                'redirect': False,
                'return': [
                    {'arg': 'out', 'type': {'-1': formats.Fasta, '': formats.Fastq}, 'suffix': ''},
                ],
            },
            'args': {
                'named': {
                    '-1': bool,
                },
                'unnamed': [
                    ('in_*', [formats.Fastq, 3, ',']),
                    ('out*', str),
                ],
            }
        }

