#!/bin/bash

git clone https://github.com/BenLangmead/bowtie2
cd bowtie2
make
mv bowtie2* "$HOME/.pypipe"
