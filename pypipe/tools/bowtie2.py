from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("bowtie2.sh", "bowtie2")


@tool
def bowtie2():
    return {
        'cmd': 'bowtie2',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': '-S', 'type': Sam, 'suffix': ''},
            ]
        },
        'args': {
            'named': {
                '-x*': Bowtie2Index,
                '-S*': str,
                '-U': {
                    '--qseq': [Qseq, 1, ','],
                    '-f':     [Fasta, 1, ','],
                    '-r':     [TextFile, 1, ','],
                    '-c':     str,
                    '':       [Fastq, 1, ',']
                    },
                '-1': {
                    '--qseq': [Qseq, 1, ','],
                    '-f':     [Fasta, 1, ','],
                    '-r':     [TextFile, 1, ','],
                    '-c':     str,
                    '':       [Fastq, 1, ',']
                    },
                '-2': {
                    '--qseq': [Qseq, 1, ','],
                    '-f':     [Fasta, 1, ','],
                    '-r':     [TextFile, 1, ','],
                    '-c':     str,
                    '':       [Fastq, 1, ',']
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

