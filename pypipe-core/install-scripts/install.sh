#!/bin/bash

tmpdir=`mktemp -d`
cd $tmpdir
bash $1
rm -rf $tmpdir
