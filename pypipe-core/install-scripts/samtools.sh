#!/bin/bash

git clone --branch=bcftools+calling git://github.com/samtools/htslib.git
git clone git://github.com/samtools/samtools.git
cd samtools
make
mv samtools "$HOME/.pypipe"
