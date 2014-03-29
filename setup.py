import os
from setuptools import setup


try:
    home_dir = os.environ['HOME']
    pypipe_dir = ".pypipe"
    os.makedirs(os.path.join(home_dir, pypipe_dir))
except:
    pass


setup(
    name="pypipe",
    version="0.1",
    description="Bioinformatics pipeline framework",
    author="semkagtn",
    author_email="semkagtn@gmail.com",
    packages=["pypipe", "pypipe.tools"],
)

