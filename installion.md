---
layout: page
title: "Installation"
group: navigation
---

{% include JB/setup %}

The easiest way to install __kallisto__ is by installing a pre-built binary.
Otherwise, you will need to build from source.

## From binary

TODO: include pre-built binaries

## From source

Currently, kallisto can be built on Linux and Mac. We have tested it on Mac OS
X and Linux (Ubuntu and CentOS).

### Requirements

If building on Mac, we suggest using a package manager such as
[Homebrew](http://brew.sh) to deal with dependencies. Other dependencies are
either included, or can be installed using package managers on the system.

- CMake version >= 2.8.8
    - Mac: `brew install cmake`
    - Ubuntu: `sudo apt-get install cmake`
    - CentOS: `sudo yum install cmake`
- zlib (should be installed on OSX >= 10.9)
    - Mac: Should be installed by default
    - Ubuntu: `sudo apt-get install zlib1g-dev`
    - CentOS: `sudo yum install zlib-devel`
- HDF5 C library
    - Mac: `brew install hdf5`
    - Ubuntu: `sudo apt-get install libhdf5-dev`
    - CentOS: `sudo yum install hdf5-devel`
