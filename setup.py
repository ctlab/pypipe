import distutils
import os
from setuptools import setup

from pypipe.paths import INSTALL_DIR

try:
    os.mkdirs(INSTALL_DIR)
except AttributeError:
    pass
distutils.dir_util.copy_tree("install-scripts", INSTALL_DIR)

setup(
    name="pypipe",
    version="0.1",
    description="Bioinformatics pipeline framework",
    author="semkagtn",
    author_email="semkagtn@gmail.com",
    packages=["pypipe", "pypipe.tools"],
)

