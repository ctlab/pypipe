#!/bim/bash

git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
make
mv bin/* "$HOME/.pypipe"
