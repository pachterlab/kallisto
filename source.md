---
layout: page
title: "Source"
description: ""
group: navigation
---
{% include JB/setup %}

The __kallisto__ source code is available [here](https://github.com/pachterlab/kallisto). Currently, __kallisto__ can be built on Linux and Mac. If building on Mac, we suggest using a package manager such as [Homebrew](http://brew.sh) to download dependencies. Homebrew is easily installed by copying and pasting the command below at a terminal prompt:

`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`


Other dependencies are either included, or can be installed using package managers on the system.

#### Requirements: 

- __CMake__ version >= 2.8.12
    - Mac: `brew install cmake`
    - Ubuntu: `sudo apt-get install cmake`
    - CentOS: `sudo yum install cmake`
- __zlib__ (should be installed on OSX >= 10.9)
    - Mac: Should be installed by default
    - Ubuntu: `sudo apt-get install zlib1g-dev`
    - CentOS: `sudo yum install zlib-devel`
- __HDF5__ C library
    - Mac: `brew install hdf5`
    - Ubuntu: `sudo apt-get install libhdf5-dev`
    - CentOS: `sudo yum install hdf5-devel`

#### Download

__kallisto__ is hosted on GitHub. The source code can be obtained by cloning the repository as follows:

`git clone git@github.com:pachterlab/kallisto.git`

#### Compile

Begin by moving to the source directory:

`cd kallisto`

Make a build directory and move there:

`mkdir build`

`cd build`

Run cmake:

`cmake ..`

Build the code:

`make`

The __kallisto__ executable is now located in `build/src`. To install __kallisto__  into the cmake install prefix path type:

`make install`

#### Test

The __kallisto__ source code package comes with a small test transcriptome and read files that can be used to test that the package was compiled and installed correctly. The [Getting started](starting.html) page provides a complete description of how to run __kallisto__ on these files.

