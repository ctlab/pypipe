from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("freebayes.sh", "freebayes")


@tool
def freebayes():
    return {
        'cmd': 'freebayes',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'v', 'type': Vcf, 'suffix': ''},
            ]
        },
        'args': {
            'named': {
                '-v*': str,
                '-f*': Fasta,
                '-b': Bam,
                '-t': Bed,
                '-r': str,
                '-s': TextFile,
                '--populations': TextFile
                '-A': Bed,
                '--trace': str,
                '--failed-alleles': Bed,
                '--variant-input': Vcf,
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
                '--observation-bias': TextFile,
                '--base-quality-cap': int,
                '--experimental-gls': bool,
                '--prob-contamination': float,
                '--contamination-estimates': TextFile,
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
                ('in_*', [Bam, 1, ' ']),
            ]
        }
    }
    
