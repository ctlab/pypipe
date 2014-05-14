import distutils
import os
from setuptools import setup

from pypipe.paths import INSTALL_DIR, INSTALL_DIR_NAME


try:
    os.mkdirs(INSTALL_DIR)
except AttributeError:
    pass
distutils.dir_util.copy_tree(INSTALL_DIR_NAME, INSTALL_DIR)

setup(
    name='pypipe',
    version='0.9',
    description='Bioinformatics pipeline framework',
    author='semkagtn',
    author_email='semkagtn@gmail.com',
    packages=['pypipe', 'pypipe.tools'],
)

