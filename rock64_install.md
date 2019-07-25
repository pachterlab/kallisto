---
layout: page
title: "Rock64 Install"
---

{% include JB/setup %}

### Compile kallisto from Source on a Rock64 (ARM architecture)

0. Navigate to home directory: `cd ~`
1. Install automake and autoreconf by the following commands: `sudo apt install automake` and `sudo apt install autoconf`.
2. Download, compile and install `cmake` by following the instructions on the [cmake website](https://cmake.org/install/).
3. From https://pachterlab.github.io/kallisto/download, download most recent version of kallisto_Rock64 file (**kallisto_rock64-v0.46.0.tar.gz**).
4. Uncompress the file: `gunzip kallisto_rock64-v0.46.0.tar.gz`
5. Untar the file: `tar -xvf kallisto_rock64-v0.46.0.tar`
6. Navigate into the Kallisto folder: *cd kallisto*
7. Now, autoreconfigure in the htslib directory: `autoreconf ext/htslib/`
8. Create a new folder called build: `mkdir build`
9. To compile from the parent directory, perform the following command: `cmake ..`
10. To finalize the process, perform the command: `make`
11. To install in user local bin, perform the command: `sudo make install`
