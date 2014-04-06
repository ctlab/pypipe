import os
<<<<<<< .merge_file_1iw7eZ
from setuptools import setup


try:
    home_dir = os.environ['HOME']
    pypipe_dir = ".pypipe"
    os.makedirs(os.path.join(home_dir, pypipe_dir))
except:
    pass

=======
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
>>>>>>> .merge_file_VSZChV

setup(
    name="pypipe",
    version="0.1",
    description="Bioinformatics pipeline framework",
    author="semkagtn",
    author_email="semkagtn@gmail.com",
    packages=["pypipe", "pypipe.tools"],
)

