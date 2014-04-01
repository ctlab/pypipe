import os
import distutils
from setuptools import setup


home_dir = os.environ['HOME']
pypipe_dir = os.path.join(home_dir, ".pypipe")
install_dir = os.path.join(pypipe_dir, "install-scripts")
try:
    os.mkdirs(install_dir)
except AttributeError:
    pass
distutils.dir_util.copy_tree("install-scripts", install_dir)

setup(
    name="pypipe",
    version="0.1",
    description="Bioinformatics pipeline framework",
    author="semkagtn",
    author_email="semkagtn@gmail.com",
    packages=["pypipe", "pypipe.tools"],
)

