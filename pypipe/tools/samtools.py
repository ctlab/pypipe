from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("samtools.sh", "samtools")


@tool
def view():
    return {
        'cmd': 'samtools view',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': '-o', 'type': {'-b': Bam, '': Sam}, 'suffix': ''},
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
                '-R': Bed,
                '-S': bool,
                '-c': bool,
                '-t': TextFile,
                '-u': bool,
            },
            'unnamed': [
                ('in_*', {'-S': Sam, '': Bam}),
            ],
        },
    }


@tool
def mpileup():
    return {
        'cmd': 'samtools mpileup',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'-u': Bcf, '-g': Bcf, '': Pileup}, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
                '-6': bool,
                '-A': bool,
                '-B': bool,
                '-b': TextFile,
                '-C': int,
                '-E': bool,
                '-f': Fasta,
                '-l': Bed,
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
                '-L': int,
                '-o': int,
                '-P': str,
            },
            'unnamed': [
                ('in_*', [Bam, 1, ' ']),
                ('out*', str),
            ]
        }
    }



@tool
def cat():
    return {
        'cmd': 'samtools cat',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'o', 'type': Bam, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
                '-h': Sam,
            },
            'unnamed': [
                ('in_*', Bam),
                ('o*', str)
            ],
        },
    }


@tool
def sort():
    return {
        'cmd': 'samtools sort',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'out', 'type': Bam, 'suffix': '.bam'},
            ],
        },
        'args': {
            'named': {
                '-n': bool,
                '-m': int,
            },
            'unnamed': [
                ('in_*', Bam),
                ('out*', str)
            ],
        },
    }


@tool
def merge():
    return {
        'cmd': 'samtools merge',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'out', 'type': Bam, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
                '-1': bool,
                '-f': bool,
                '-h': Sam,
                '-n': bool,
                '-R': str,
                '-r': bool,
                '-u': bool,
            },
            'unnamed': [
                ('out*', str),
                ('in_*', [Bam, 2, ' ']),
            ],
        },
    }


@tool
def rmdup():
    return {
        'cmd': 'samtools rmdup',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'out', 'type': Bam, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
                '-s': bool,
                '-S': bool,
            },
            'unnamed': [
                ('in_*', Bam),
                ('out*', str),
            ],
        },
    }


@tool
def calmd():
    return {
        'cmd': 'samtools calmd',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'-u': Bam, 'b': Bam, '': Sam}, 'suffix': ''},
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
                ('in_*', {'-S': Sam, '': Bam}),
                ('out*', str),
            ],
        },
    }

