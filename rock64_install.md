---
layout: page
title: "Rock64 Install"
---

{% include JB/setup %}

### Compiling kallisto from source on a Rock64

1. Install automake, autoconf and libhdf5-dev with `sudo apt install automake`, `sudo apt install autoconf` and `sudo apt install hdf5-dev`. You may have to remove a dpkg lock with `sudo rm /var/lib/dpkg/lock`.
2. Download, compile and install `cmake` by following the instructions on the [cmake website](https://cmake.org/install/).
3. Download the kallisto source from https://pachterlab.github.io/kallisto/download.
4. Gunzip and untar the source code file. 
5. Navigate to the source code directory and type `autoreconf ext/htslib/`.
6. Create a new folder called build: `mkdir build`.
7. Navigate into the build directory and type `cmake ..`.
8. Then type `make` followed by `sudo make install` to place kallisto in `/usr/local/bin`.
