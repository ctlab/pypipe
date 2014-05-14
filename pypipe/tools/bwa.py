from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("bwa.sh", "bwa")


@tool
def mem():
    return {
        'cmd': 'bwa mem',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Sam, 'suffix': ''},
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
                '-a': bool,
                '-C': bool,
                '-H': bool,
                '-M': bool,
                '-v': int,
            },
            'unnamed': [
                ('ref*', BwaIndex),
                ('in1*', Fastq),
                ('in2', Fastq),
                ('out*', str),
            ],
        },
    }



@tool
def aln():
    return {
        'cmd': 'bwa aln',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Sai, 'suffix': ''},
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
                ('ref*', BwaIndex),
                ('in_*', Fastq),
                ('out*', str),
            ],
        },
    }


@tool
def samse():
    return {
        'cmd': 'bwa samse',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Sam, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
                '-n': int,
                '-r': str,
            },
            'unnamed': [
                ('ref*', BwaIndex),
                ('sai*', Sai),
                ('in_*', Fastq),
                ('out*', str),
            ],
        },
    }


@tool
def sampe():
    return {
        'cmd': 'bwa sampe',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Sam, 'suffix': ''},
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
                ('ref*', BwaIndex),
                ('sai1*', Sai),
                ('sai2*', Sai),
                ('in1*', Fastq),
                ('in2*', Fastq),
                ('out*', str),
            ],
        },
    }


@tool
def bwasw():
    return {
        'cmd': 'bwa bwasw',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Sam, 'suffix': ''},
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
                ('ref*', BwaIndex),
                ('in1*', Fastq),
                ('in2', Fastq),
                ('out*', str),
            ],
        },
    }

