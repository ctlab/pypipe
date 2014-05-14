from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("VarScan.sh", "VarScan.jar")


@tool
def pileup2snp():
    return {
        'cmd': 'java -jar VarScan.jar pileup2snp',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Snp, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }


@tool
def pileup2indel():
    return {
        'cmd': 'java -jar VarScan.jar pileup2indel',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Indel, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }


@tool
def pileup2cns():
    return {
        'cmd': 'java -jar VarScan.jar pileup2cns',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Cns, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }


@tool
def mpileup2snp():
    return {
        'cmd': 'java -jar VarScan.jar mpileup2snp',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'--output-vcf': Vcf, '': Snp}, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }


@tool
def mpileup2indel():
    return {
        'cmd': 'java -jar VarScan.jar mpileup2indel',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'--output-vcf': Vcf, '': Snp}, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }


@tool
def mpileup2cns():
    return {
        'cmd': 'java -jar VarScan.jar mpileup2cns',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'--output-vcf': Vcf, '': Snp}, 'suffix': ''},
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
                ('in_*', Pileup),
                ('out*', str),
            ],
        },
    }

@tool
def somatic():
    return {
        'cmd': 'java -jar VarScan.jar somatic',
        'type': 'jar',
        'log': 'log',
        'out': {
            'redirect': False,
            'return': [
                {'arg': 'output', 'type': {'--output-snp': None, '': Snp}, 'suffix': '.indel'},
                {'arg': 'output', 'type': {'--output-indel': None, '': Indel}, 'suffix': '.snp'},
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
                ('normal*', Pileup),
                ('tumor*', Pileup),
                ('output*', str),
            ],
        },
    }

