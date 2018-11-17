#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: ./gen_release.sh OS_NAME"
    exit 1
fi

OS=$1

RELEASE=release
if [ -d "$RELEASE" ]; then
    rm -rf $RELEASE
fi

mkdir $RELEASE

cd $RELEASE
cmake -DLINK=static ..
make VERBOSE=1
mkdir kallisto

cp src/kallisto kallisto
cp -r ../test kallisto

cp ../license.txt kallisto
cp ../README.md kallisto

TAR_NAME="kallisto_${OS}.tar"

tar -cvvf $TAR_NAME kallisto
gzip $TAR_NAME
