#!/bin/bash

git clone https://github.com/lh3/bwa.git
cd bwa
make
mv bwa "$HOME/.pypipe"
