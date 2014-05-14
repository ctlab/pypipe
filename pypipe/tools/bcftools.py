from pypipe.formats import *
from pypipe.utils import install_program, tool


install_program("bcftools.sh", "bcftools")

@tool
def view():
    return {
        'cmd': 'bcftools view',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': {'b': Bcf, '': Vcf}, 'suffix': ''},
            ]
        },
        'args': {
            'named': {
                '-A': bool,
                '-b': bool,
                '-D': TextFile,
                '-F': bool,
                '-G': bool,
                '-l': TextFile,
                '-N': bool,
                '-Q': bool,
                '-s': TextFile,
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
                ('in_*', {'-S': Vcf, '': Bcf}),
                ('out*', str),
            ],
        },
    }


@tool
def cat(in_, out, log=None):
    return {
        'cmd': 'bcftools cat',
        'type': None,
        'log': 'log',
        'out': {
            'redirect': True,
            'return': [
                {'arg': 'out', 'type': Bcf, 'suffix': ''},
            ],
        },
        'args': {
            'named': {
            },
            'unnamed': {
                'in_*': Bcf,
                'out*': str,
            },
        }
    }
    program = pipeline.add_node("bcftools cat", log, out)
    program.add_args(in_, formats.Bcf)
    return formats.Bcf(out, program)

