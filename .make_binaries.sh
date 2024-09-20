#!/bin/bash

# Just a script to make various kallisto binaries (works on mac; some adjustments may need to be made for other OS's, but this provides a basic framework)

os="mac"



# Standard build

rm -rf kallisto
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake .. -DUSE_HDF5=ON -DUSE_BAM=ON -DBUILD_FUNCTESTING=ON
make
version=$(./src/kallisto version|cut -d' ' -f3|tr -d '\n')
mkdir kallisto
mv ./src/kallisto ./kallisto/
cp -R ../test/ ./kallisto/test/
cp ../README.md ./kallisto/
cp ../license.txt ./kallisto/
tar --no-xattrs --exclude='._*' -czvf kallisto_${os}-v${version}.tar.gz kallisto
cp kallisto_${os}-v${version}.tar.gz ../../
cd ../../

# Long k-mer build

rm -rf kallisto
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake .. -DUSE_HDF5=ON -DUSE_BAM=ON -DBUILD_FUNCTESTING=ON -DMAX_KMER_SIZE=64
make
version=$(./src/kallisto version|cut -d' ' -f3|tr -d '\n')
cp src/kallisto ../../kallisto_${os}-v${version}_kmer64
cd ../../

# Disable opt build

rm -rf kallisto
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake .. -DUSE_HDF5=ON -DUSE_BAM=ON -DBUILD_FUNCTESTING=ON -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF
make
version=$(./src/kallisto version|cut -d' ' -f3|tr -d '\n')
mv src/kallisto ../../kallisto_${os}-v${version}_optoff
cd ../../



# Disable opt build but have long k-mer

rm -rf kallisto
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake .. -DUSE_HDF5=ON -DUSE_BAM=ON -DBUILD_FUNCTESTING=ON -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF -DMAX_KMER_SIZE=64
make
version=$(./src/kallisto version|cut -d' ' -f3|tr -d '\n')
mv src/kallisto ../../kallisto_${os}-v${version}_optoff_k64
cd ../../


