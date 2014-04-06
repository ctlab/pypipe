#!/bin/bash

git clone --branch=bcftools+calling git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
mv bcftools "$HOME/.pypipe"
