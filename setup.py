#!/usr/bin/env python2

import os
from setuptools import setup, find_packages

from pypipe.paths import PYPIPE_DIR, INSTALL_SCRIPTS_DIR
from pypipe.utils import install_program


try:
    os.makedirs(PYPIPE_DIR)
except OSError:
    pass

#  binary file name -> install script name
tools = {
    'bcftools':    'bcftools.sh',
    'bowtie2':     'bowtie2.sh',
    'bwa':         'bwa.sh',
    'freebayes':   'freebayes.sh',
    'samtools':    'samtools.sh',
    'VarScan.jar': 'varscan.sh',
}

for name in tools:
    script = tools[name]
    install_program(script, name)


setup(
    name='pypipe',
    version='0.9',
    description='Bioinformatics pipeline framework',
    author='semkagtn',
    author_email='semkagtn@gmail.com',
    packages=find_packages(),
)
