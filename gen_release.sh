#!/bin/bash

RELEASE=release
if [ -d "$RELEASE" ]; then
    rm -rf $RELEASE
fi

mkdir $RELEASE

cd $RELEASE
cmake ..
make
mkdir kallisto

cp src/kallisto kallisto
cp ../license.txt kallisto
cp ../README.md kallisto

tar -cvvf kallisto.tar kallisto
gzip kallisto.tar
