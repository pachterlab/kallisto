---
layout: page
title: "Rock64 Install"
---

{% include JB/setup %}

### How to Compile Kallisto from Source on a Rock64 (Arm64 architecture)

0. Navigate to home directory: `cd ~`
1. Install automake and autoreconf by the following commands: `sudo apt install automake` and `sudo apt install autoreconf`
2. From https://cmake.org/download/, install cmake : `wget cmake-3.15.0-rc4-Linux-x86_64.sh`
Then, perform the command  `bash cmake-3.15.0-rc4-Linux-x86_64.sh` which will run the script that installs cmake, allowing us to compile the kallisto code.
3. From https://pachterlab.github.io/kallisto/download, download most recent version of kallisto_Rock64 file (**kallisto_rock64-v0.46.0.tar.gz**).
4. Uncompress the file: `gunzip kallisto_rock64-v0.46.0.tar.gz`
5. Untar the file: `tar -xvf kallisto_rock64-v0.46.0.tar`
6. Navigate into the Kallisto folder: *cd kallisto*
7. Now, autoreconfigure in the htslib directory: `autoreconf ext/htslib/`
8. Create a new folder called build: `mkdir build`
9. To compile from the parent directory, perform the following command: `cmake ..`
10. To finalize the process, perform the command: `make`
11. To install in user local bin, perform the command: `make install`
