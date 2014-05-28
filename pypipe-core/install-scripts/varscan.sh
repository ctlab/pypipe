#!/bin/bash

wget "http://sourceforge.net/projects/varscan/files/latest/download?source=files"
mv "download?source=files" "VarScan.jar"
mv "VarScan.jar" "$HOME/.pypipe"
